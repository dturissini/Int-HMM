#!/usr/bin/perl -w
#David Turissini

use DBI;
prepareDbh();


my ($prefix, $run_name, $donor_group, $recip_group, $freq_cutoff, $repeat_regions_table, $gatk_vcf) = @ARGV;

my $time = localtime;
print STDOUT "$time\t$gatk_vcf\n";

my ($group_list, $arms_table) = $dbh->selectrow_array(qq{select group_list, arms_table
                                                         from introgression_snps_params
                                                         where run_name = '$run_name'});



my $donor_query = $dbh->prepare(qq{select ind
                                   from tracking_groups
                                   where group_name = '$donor_group'});
$donor_query->execute();
my %donors = ();
while (my ($ind) = $donor_query->fetchrow_array())
  {$donors{$ind} = 1;}


my $recip_group_query = $dbh->prepare(qq{select ind
                                         from tracking_groups
                                         where group_name = '$recip_group'});
$recip_group_query->execute();
my %recips = ();
while (my ($ind) = $recip_group_query->fetchrow_array())
  {$recips{$ind} = 1;}



my $pop_query = $dbh->prepare(qq{select ind, group_name
                                 from tracking_groups
                                 where group_name in ($group_list)});

$pop_query->execute();

my %pops = ();
my %distinct_pops = ();
while (my ($ind, $group_name) = $pop_query->fetchrow_array())
  {
  $pops{$ind} = $group_name;
  $distinct_pops{$group_name}++;
  }


my ($chrom_name_field, $chrom_len_field) = $dbh->selectrow_array(qq{select distinct chrom_name_field, chrom_len_field
                                                                    from tracking_ref_genomes
                                                                    where ref_db_table = '$arms_table'});

my $chrom_len_query = $dbh->prepare(qq{select $chrom_name_field, $chrom_len_field
                                       from $arms_table
                                       where $chrom_name_field not like '%h'
                                       and $chrom_name_field != '4'});
$chrom_len_query->execute();
my %arm_lengths = ();
while (my ($arm , $arm_len) = $chrom_len_query->fetchrow_array())
  {$arm_lengths{$arm} = $arm_len;}


foreach my $arm (keys %arm_lengths)
  {
  $arm =~ tr/\./_/;
  my $snp_table = 'sites_all_' . $prefix . '_' . $freq_cutoff * 100 . '_' . $arm;
  
  $dbh->do(qq{drop table if exists $snp_table});
  $dbh->do(qq{create table $snp_table
              (vs_id int auto_increment primary key,
               pos int,
               ind varchar(100),
               ref_allele_cov int,
               alt_allele_cov int,
               recip_ref_freq decimal(5,4),
               donor_ref_freq decimal(5,4),
               index(pos))});
  }               




my $cov_quantile_query = $dbh->prepare(qq{select distinct b.ind, cov_99_quantile
                                          from tracking_bams b, tracking_vcf_bams v
                                          where vcf_file = '$gatk_vcf'
                                          and b.bam_file = v.bam_file});
$cov_quantile_query->execute();

my %cov_quantiles = ();
while (my ($ind, $cov_99_quantile) = $cov_quantile_query->fetchrow_array())
  {$cov_quantiles{$ind} = $cov_99_quantile;}


open GATK, "$gatk_vcf";	

my @gt_inds = ();
my %ind_genotypes = ();
my %group_allele_counts = ();
my $counter = 0;
while (my $line = <GATK>)
  {
  if ($counter % 1000000 == 0)
    {
    my $time = localtime;
    print STDOUT "$time\t$counter\n";  
    }
  $counter++;

  chomp $line;
  
  if ($line =~ /^\#/)
    {
    if ($line =~ /\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT/)
      {
      @gt_inds = split /\t/, $line;
      
      for (my $i = 0; $i < scalar @gt_inds; $i++)
        {
        if ($i < 9)
          {push @groups, '';}
        else
          {
          if (exists $donors{$gt_inds[$i]})
            {push @groups, 'donor';}
          elsif (exists $recips{$gt_inds[$i]})
            {push @groups, 'recip';}
          else
            {push @groups, 'ignore';}
          }
        }
      }
    }
  else
    {
    my @values = split /\t/, $line; 

    my @alleles = split /,/, $values[4];
    unshift @alleles, $values[3];

    if (scalar @alleles == 2 and ($values[6] eq 'PASS' or $values[6] eq '.') and exists $arm_lengths{$values[0]})
      {
      my %genotypes = ();
      my %cum_genotypes = (donor => {}, recip => {}, ignore => {});
      my %group_counts = (donor => 0, recip => 0, ignore => 0);
      my %covs = ();
      for (my $i = 9; $i < scalar @values; $i++)	
        {   
        if (exists $pops{$gt_inds[$i]})   
          {  	  
          my @gt = split /\:/, $values[$i];	
          my @gt_alleles = split /\//, $gt[0];	
          my $genotype = $gt_alleles[0] . $gt_alleles[1];     
          unless ($genotype =~ /\./)
            {
            my @gt_covs = split(/,/, $gt[1]); 
                   
            my $cov = 0;
            foreach my $allele_cov (@gt_covs)
              {$cov += $allele_cov;}

            if ($cov <= $cov_quantiles{$gt_inds[$i]})
              {
              $covs{$gt_inds[$i]} = [split(/,/, $gt[1])];
              $genotypes{$groups[$i]}{$gt_inds[$i]} = $genotype;
              $cum_genotypes{$groups[$i]}{$genotype}++; 
              $group_counts{$groups[$i]}++;              
              }    
            }
          }
        }	

      if ($group_counts{donor} == scalar keys %donors and $group_counts{recip} > 0)
        {
        my $donor_freq = 0;
        my $donor_chroms = 0;
        foreach my $donor_gt (keys %{$cum_genotypes{donor}})
          {
          $donor_chroms += 2 * $cum_genotypes{donor}{$donor_gt};
          
          if ($donor_gt eq '01')
            {$donor_freq += $cum_genotypes{donor}{$donor_gt};}
          elsif ($donor_gt eq '00')
            {$donor_freq += 2 * $cum_genotypes{donor}{$donor_gt};}
          }
        $donor_freq = $donor_freq / $donor_chroms;
        
        my $recip_freq = 0;
        my $recip_het_freq = 0;
        my $recip_chroms = 0;
        foreach my $recip_gt (keys %{$cum_genotypes{recip}})
          {
          $recip_chroms += 2 * $cum_genotypes{recip}{$recip_gt};
          
          if ($recip_gt eq '01')
            {
            $recip_freq += $cum_genotypes{recip}{$recip_gt};
            $recip_het_freq += $cum_genotypes{recip}{$recip_gt};
            }
          elsif ($recip_gt eq '00')
            {$recip_freq += 2 * $cum_genotypes{recip}{$recip_gt};}
          }
        $recip_freq = $recip_freq / $recip_chroms;
        $recip_het_freq = $recip_het_freq / ($recip_chroms / 2);
                 
        if (($donor_freq == 1 or $donor_freq == 0) and abs($donor_freq - $recip_freq) >= $freq_cutoff and $recip_het_freq <= .8)
          {
          my $snp_table = 'sites_all_' . $prefix . '_' . $freq_cutoff * 100 . '_' . $values[0];
          foreach my $ind (sort keys %{$genotypes{recip}})
            {
            my $ref_allele_cov = $covs{$ind}[0];
            my $alt_allele_cov = $covs{$ind}[1];
            
            $dbh->do(qq{insert into $snp_table
                        values
                        (null, $values[1], '$ind', $ref_allele_cov, $alt_allele_cov, $recip_freq, $donor_freq)});
            }
          }  
        }
      } #end if (scalar @alleles > 2)
    } #end if ($line =~ /^\#/)
  }


foreach my $arm (keys %arm_lengths)
  {
  my $snp_table = 'sites_all_' . $prefix . '_' . $freq_cutoff * 100 . '_' . $arm;
  my $tmp_repeat_snps_table = 'tmp_' . $prefix . '_' . $freq_cutoff * 100 . '_repeat_snps';

  $dbh->do(qq{drop table if exists $tmp_repeat_snps_table});

#I'm not using a distinct in this query because it changed the mysql explain plan for an unknown reason
  $dbh->do(qq{create table $tmp_repeat_snps_table as
              select straight_join pos
              from $repeat_regions_table r, $snp_table s
              where pos between start and end
              and arm = '$arm';});

  $dbh->do(qq{alter table $tmp_repeat_snps_table add index(pos)});

  $dbh->do(qq{delete from $snp_table
              where pos in (select distinct pos from $tmp_repeat_snps_table)});
              
  $dbh->do(qq{drop table $tmp_repeat_snps_table});
  }

$time = localtime;
print STDOUT "$time\tdone\n";


sub prepareDbh
  {
  my $username = "";
  my $password = "";
  my $database = "";
  my $host = "";
  my $dsn = "DBI:mysql:$database:$host";
  our $dbh = DBI->connect($dsn,$username,$password);
  }


