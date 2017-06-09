#!/usr/bin/perl -w
#David Turissini

use Storable qw(dclone);
use DBI;
prepareDbh();
 

my ($prefix, $introgress_group, $hmm_suffix, $hmm_type, $repeats_table) = @ARGV;


my $combined_tracks_table = $prefix . '_hmm_' . $hmm_type . '_pab_' . $hmm_suffix . '_tracks';
my $combined_tracks_fil_table = $prefix . '_hmm_' . $hmm_type . '_pab_' . $hmm_suffix . '_tracks_fil';
my $combined_unftrk_table = $prefix . '_hmm_' . $hmm_type . '_pab_' . $hmm_suffix . '_unftrk';
my $regions_table = $prefix . '_hmm_' . $hmm_type . '_pab_' . $hmm_suffix . '_regions';
my $regions_inds_table = $prefix . '_hmm_' . $hmm_type . '_pab_' . $hmm_suffix . '_regions_inds';
my $blocks_table = $prefix . '_hmm_' . $hmm_type . '_pab_' . $hmm_suffix . '_blocks';
my $blocks_inds_table = $prefix . '_hmm_' . $hmm_type . '_pab_' . $hmm_suffix . '_blocks_inds';
my $tmp_repeat_table = $prefix . '_hmm_' . $hmm_type . '_pab_' . $hmm_suffix . '_tmprep';



$dbh->do(qq{drop table if exists $combined_tracks_table});
$dbh->do(qq{create table $combined_tracks_table
            (cict_id int auto_increment primary key,
             region_id int,           
             ind varchar(30),
             genotype varchar(20),
             arm varchar(5),
             start int,
             end int,
             total_snps int,
             introgress_snps int,
             track_len int,
             cum_cov int,
             per_repeat decimal(5,2),
             index(ind),
             index(start))});

$dbh->do(qq{drop table if exists $combined_tracks_fil_table});
$dbh->do(qq{create table $combined_tracks_fil_table
            (cict_id int auto_increment primary key,
             region_id int,           
             ind varchar(30),
             genotype varchar(20),
             arm varchar(5),
             start int,
             end int,
             total_snps int,
             introgress_snps int,
             track_len int,
             cum_cov int,
             per_repeat decimal(5,2),
             index(ind),
             index(start))});

$dbh->do(qq{drop table if exists $combined_unftrk_table});
$dbh->do(qq{create table $combined_unftrk_table
            (cuft_id int auto_increment primary key,
             ind varchar(30),
             genotype varchar(20),
             arm varchar(5),
             start int,
             end int,
             total_snps int,
             introgress_snps int,
             track_len int,
             cum_cov int,
             index(ind),
             index(start))});


$dbh->do(qq{drop table if exists $regions_table});
$dbh->do(qq{create table $regions_table
            (region_id int auto_increment primary key,
             arm varchar(5),
             start int,
             end int,
             inds int,
             max_introgress_snps int,
             distinct_haps int,
             index(start))});


$dbh->do(qq{drop table if exists $regions_inds_table});
$dbh->do(qq{create table $regions_inds_table
            (fssfhri_id int auto_increment primary key,
             region_id int,
             ind varchar(100),
             index(region_id),
             index(ind))});


$dbh->do(qq{drop table if exists $blocks_table});
$dbh->do(qq{create table $blocks_table
            (block_id int auto_increment primary key,
             arm varchar(5),
             start int,
             end int,
             inds int,
             block_len int,
             per_repeat decimal(5,2),
             index(start))});


$dbh->do(qq{drop table if exists $blocks_inds_table});
$dbh->do(qq{create table $blocks_inds_table
            (bi_id int auto_increment primary key,
             block_id int,
             ind varchar(100),
             index(block_id),
             index(ind))});



my @arms = ('2L', '2R', '3L', '3R', 'X');

my $ind_query = $dbh->prepare(qq{select distinct ind
                                 from tracking_groups
                                 where group_name in ('$introgress_group')});  
$ind_query->execute();

my $total_inds = 0;
while (my ($ind) = $ind_query->fetchrow_array())
  {
  $total_inds++;
  $ind =~ tr/\./\_/;
  
  my $tracks_table = $prefix . '_' . $ind . '_hmm_' . $hmm_type . '_pab_' . $hmm_suffix . '_tracks';
  $dbh->do(qq{insert into $combined_tracks_table
              select null, null, ind, genotype, arm, start, end, total_snps, introgress_snps, track_len, cum_cov, per_repeat
              from $tracks_table
              order by arm, start});

  my $unftrk_table = $prefix . '_' . $ind . '_hmm_' . $hmm_type . '_pab_' . $hmm_suffix . '_unftrk';
  $dbh->do(qq{insert into $combined_unftrk_table
              select null, ind, genotype, arm, start, end, total_snps, introgress_snps, track_len, cum_cov
              from $unftrk_table
              order by arm, start});



  #tracks_fil table
  foreach my $arm (@arms)
    {
    my $filter_query = $dbh->prepare(qq{select straight_join cict_id cict_id_start, cict_id + 1 cict_id_end, 'homo_1'
                                        from $combined_tracks_table
                                        where cict_id in (select min(cict_id) 
                                                          from $combined_tracks_table
                                                          where arm = '$arm'
                                                          and ind = '$ind')
                                        and genotype in ('het', 'homo_2')
                                        and ind = '$ind'
                                        and (end - start + 1 <= 500
                                             or introgress_snps < 10
                                             or per_repeat >= 30)
                                        union all
                                        select cict_id - 1, cict_id, 'homo_1'
                                        from $combined_tracks_table
                                        where cict_id in (select max(cict_id) 
                                                          from $combined_tracks_table
                                                          where arm = '$arm'
                                                          and ind = '$ind')
                                        and genotype in ('het', 'homo_2')
                                        and ind = '$ind'
                                        and (end - start + 1 <= 500
                                             or introgress_snps < 10
                                             or per_repeat >= 30)
                                        union all
                                        select straight_join t.cict_id, t3.cict_id, 'homo_1'
                                        from $combined_tracks_table t, $combined_tracks_table t3, $combined_tracks_table t2
                                        where t2.cict_id = t.cict_id + 1
                                        and t3.cict_id = t2.cict_id + 1
                                        and t.genotype = 'homo_1'
                                        and t2.genotype in ('het', 'homo_2')
                                        and t3.genotype = 'homo_1'
                                        and t.ind = '$ind'
                                        and t2.ind = t.ind
                                        and t3.ind = t.ind
                                        and (t2.end - t2.start + 1 <= 500
                                             or t2.introgress_snps < 10
                                             or t2.per_repeat >= 30)
                                        order by cict_id_start, cict_id_end});
    
    $filter_query->execute();
    
    my %filtered_cict_ids = ();
    while (my ($cict_id_start, $cict_id_end, $genotype) = $filter_query->fetchrow_array())
      {
      $filtered_cict_ids{$cict_id_start}{end} = $cict_id_end;
      $filtered_cict_ids{$cict_id_start}{genotype} = $genotype;
      }  
      
    my $track_query = $dbh->prepare(qq{select cict_id 
                                       from $combined_tracks_table
                                       where arm = '$arm'
                                       and ind = '$ind'
                                       order by cict_id});
    $track_query->execute();
    
    while (my ($cict_id) = $track_query->fetchrow_array())
      {
      my $filtered = 'N';
      foreach my $cict_id_start (keys %filtered_cict_ids)
        {
        if ($cict_id >= $cict_id_start and $cict_id <= $filtered_cict_ids{$cict_id_start}{end})
          {
          $filtered = 'Y';
          if ($cict_id == $cict_id_start)
            {
            $dbh->do(qq{insert into $combined_tracks_fil_table
                        select null, null, ind, '$filtered_cict_ids{$cict_id_start}{genotype}',
                               arm, min(start), max(end),
                               sum(total_snps), sum(introgress_snps), max(end) - min(start) + 1, sum(cum_cov), null
                        from $combined_tracks_table
                        where arm = '$arm'
                        and ind = '$ind'
                        and cict_id between $cict_id_start and $filtered_cict_ids{$cict_id_start}{end}
                        group by ind, arm, '$filtered_cict_ids{$cict_id_start}{genotype}'}); 
            }                                                                         
          }   
        elsif ($cict_id > $filtered_cict_ids{$cict_id_start}{end})
          {delete $filtered_cict_ids{$cict_id_start};}                                                             
        }
    
    
    
      if ($filtered eq 'N')
        {
        $dbh->do(qq{insert into $combined_tracks_fil_table
                    select null, null, ind, genotype, arm, start, end, total_snps, introgress_snps, track_len, cum_cov, per_repeat
                    from $combined_tracks_table
                    where cict_id = $cict_id});
        }
      }
    }
  }  



foreach my $arm (@arms)
  {
  my $time = localtime;
  print "$time\t$arm\n";
  
  $time = localtime;
  print "$time\nidentify regions\n";
  
  my $region_query = $dbh->prepare(qq{select start, end, count(*), max(introgress_snps)
                                      from $combined_tracks_table 
                                      where arm = '$arm'
                                      and genotype != 'homo_1'
                                      and per_repeat < 30
                                      group by start, end
                                      order by start, end});
  $region_query->execute();
  
  my $region_start = -1;
  my $last_end = -1;
  my $max_introgress_snps = 0;
  my %haps = ();
  while (my ($start, $end, $inds, $introgress_snps) = $region_query->fetchrow_array())
    {
    if ($region_start == -1)
      {$region_start = $start;}
    
    if ($start > $last_end and $last_end != -1)
      {
      my $distinct_haps = scalar keys %haps;
      $dbh->do(qq{insert into $regions_table
                  values
                  (null, '$arm', $region_start, $last_end, null, $max_introgress_snps, $distinct_haps)});
      
      %haps = ();
      $region_start = $start;
      $max_introgress_snps = 0;
      }
    
    if ($end > $last_end)
      {$last_end = $end;}
      
    if ($introgress_snps > $max_introgress_snps)
      {$max_introgress_snps = $introgress_snps} 
    
    my $hap = $start . '-' . $end;
    $haps{$hap}++;  
    } 

  my $distinct_haps = scalar keys %haps;
  $dbh->do(qq{insert into $regions_table
              values
              (null, '$arm', $region_start, $last_end, null, $max_introgress_snps, $distinct_haps)});    



  #make blocks
  my $breakpoint_query = $dbh->prepare(qq{select pos, sum(if(track_end = 'start', 1, 0)) starts, 
                                          sum(if(track_end = 'end', 1, 0)) ends,
                                          sum(if(track_end = 'start' and genotype = 'het', 1, 0)) het_starts, 
                                          sum(if(track_end = 'end' and genotype = 'het', 1, 0)) het_ends
                                          from (select start pos, 'start' track_end, genotype
                                                from $combined_tracks_table
                                                where arm = '$arm' 
                                                and end - start + 1 > 500
                                                and introgress_snps >= 10
                                                and per_repeat < 30
                                                and genotype != 'homo_1'
                                                union all
                                                select end, 'end' track_end, genotype 
                                                from $combined_tracks_table
                                                where arm = '$arm' 
                                                and end - start + 1 > 500
                                                and introgress_snps >= 10
                                                and per_repeat < 30
                                                and genotype != 'homo_1') x
                                          group by pos});
  
  $breakpoint_query -> execute();
  
  my $last_pos = -1;
  my $num_inds = 0;
  my $num_hets = 0;
  while (my ($pos, $num_starts, $num_ends, $num_het_starts, $num_het_ends) = $breakpoint_query->fetchrow_array())
    {
  #print "$pos, $num_starts, $num_ends\n";
  
    unless ($num_inds == 0)
      {
      my $block_start = $last_pos;
      my $block_end = $pos;
      if ($num_starts > 0)
        {$block_end = $pos - 1;}
      
      $dbh->do(qq{insert into $blocks_table
                  values
                  (null, '$arm', $block_start, $block_end, $num_inds, $num_het_inds, $block_end - $block_start + 1, 0)});
      
      
      if ($num_starts > 0 and $num_ends > 0)
        {
        $dbh->do(qq{insert into $blocks_table
                    values
                    (null, '$arm', $pos, $pos, $num_inds + $num_starts, $num_het_inds + $num_het_starts, 1, 0)});
        }
      }
      
    $num_inds += $num_starts - $num_ends;
    $num_het_inds += $num_het_starts - $num_het_ends;
    $last_pos = $pos;
    if ($num_ends > 0)
      {$last_pos++;}
    }
  }


#update block per_repeat
$dbh->do(qq{create table $tmp_repeat_table
            as select t.block_id, sum(rep_len) / (end - start + 1) * 100 per_repeat
            from (select straight_join block_id, r.end - t.start + 1 rep_len
                  from $blocks_table t, $repeats_table r
                  where t.arm = r.arm
                  and r.start < t.start
                  and r.end between t.start and t.end
                  union all
                  select straight_join block_id, t.end - r.start + 1
                  from $blocks_table t, $repeats_table r
                  where t.arm = r.arm
                  and r.start between t.start and t.end
                  and r.end > t.end
                  union all
                  select straight_join block_id, r.end - r.start + 1
                  from $blocks_table t, $repeats_table r
                  where t.arm = r.arm
                  and r.start between t.start and t.end
                  and r.end between t.start and t.end
                  union all
                  select straight_join block_id, t.end - t.start + 1
                  from $repeats_table r, $blocks_table t
                  where t.arm = r.arm
                  and t.start between r.start and r.end
                  and t.end between r.start and r.end
                  and t.start > r.start
                  and t.end < r.end) x, $blocks_table t
            where t.block_id = x.block_id
            group by t.block_id, start, end});

$dbh->do(qq{alter table $tmp_repeat_table add unique index(block_id)});


$dbh->do(qq{update $blocks_table t, $tmp_repeat_table x
            set t.per_repeat = x.per_repeat
            where t.block_id = x.block_id});

$dbh->do(qq{drop table $tmp_repeat_table});



$dbh->do(qq{insert into $blocks_inds_table
            select straight_join null, block_id, ind
            from $combined_tracks_table t, $blocks_table b
            where b.arm = t.arm
            and b.start between t.start and t.end
            and t.end - t.start + 1 > 500
            and introgress_snps >= 10
            and t.per_repeat < 30
            and genotype != 'homo_1'});



#filter het stacks
my $het_block_query = $dbh->prepare(qq{select arm, start, end, inds, het_inds
                                       from $blocks_table 
                                       where het_inds / inds > .9
                                       and inds >= .25 * $total_inds});

$het_block_query->execute();

while (my () = $het_block_query->fetchrow_array())
  {
  
  }





$dbh->do(qq{update $combined_tracks_table t, $regions_table r
            set t.region_id = r.region_id
            where t.arm = r.arm
            and t.start between r.start and r.end
            and t.end between r.start and r.end
            and genotype != 'homo_1'});   
    
  


$time = localtime;
print "$time\tpopulating regions_inds table\n";
$dbh->do(qq{insert into $regions_inds_table
            select null, region_id, ind
            from $combined_tracks_table f
            group by null, region_id, ind});


$dbh->do(qq{update $regions_table r, (select region_id, count(distinct ind) total_inds
                                                 from $regions_inds_table
                                                 group by region_id) x
            set inds = total_inds
            where r.region_id = x.region_id});



$time = localtime;
print "$time\tdone\n";




sub prepareDbh
  {
  my $username = "";
  my $password = "";
  my $database = "";
  my $host = "";
  my $dsn = "DBI:mysql:$database:$host";
  our $dbh = DBI->connect($dsn,$username,$password);
  }



