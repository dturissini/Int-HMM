#!/usr/bin/perl -w
#David Turissini


use DBI;
prepareDbh();

my ($prefix, $recip_group, $pab, $base_a_prob, $err_multiplier, $freq_diff_cutoff) = @ARGV;

my $ind_query = $dbh->prepare(qq{select ind
                                 from tracking_groups
                                 where group_name = '$recip_group'});
$ind_query->execute();

while (my ($ind) = $ind_query->fetchrow_array())
  {
  my $ind_label = $ind;
  $ind_label =~ s/\./\_/g;
  my $hmm_table = $prefix . '_' . $ind_label . '_hmm_dsm_pab_' . '0' * (2 - length(int($pab * 100))) . int($pab * 100) . '_' . int(-1 * log10($base_a_prob) + .5) . '_' . $err_multiplier . '_all_' . int($freq_diff_cutoff * 100);

  if ($pab < .01)
    {$hmm_table = $prefix . '_' . $ind_label . '_hmm_dsm_pab_' . '0' * (3 - length(int($pab * 1000))) . int($pab * 1000) . '_' . int(-1 * log10($base_a_prob) + .5) . '_' . $err_multiplier + '_all_' . int($freq_diff_cutoff * 100);}


  $dbh->do(qq{drop table if exists $hmm_table}); 
  $dbh->do(qq{create table $hmm_table
              (ch_id int auto_increment primary key,
               ind varchar(30),
               arm varchar(5),
               pos int,
               p_homo_1 decimal(5,4),
               p_het decimal(5,4),
               p_homo_2 decimal(5,4),
               p_homo_1_err decimal(5,4),
               p_het_err decimal(5,4),
               p_homo_2_err decimal(5,4),
               a_cov int,
               b_cov int,
               genotype varchar(20),
               index(pos))}); 
  
  
  foreach my $arm ('2L', '2R', '3L', '3R', 'X') 
    { 
    my $err_file = 'logs/err_hmm_introgress_ds_all_' . $freq_diff_cutoff . '_' . $ind . '_' . $arm . '.%J';
    my $out_file = 'logs/out_hmm_introgress_ds_all_' . $freq_diff_cutoff . '_' . $ind . '_' . $arm . '.%J';
    
    system(qq{bsub -e $err_file -o $out_file -M 16 -q week python hmm_introgress_error_states_all_downsample_single_arm.py $ind $prefix $arm $pab $base_a_prob $err_multiplier $freq_diff_cutoff});
    }
  }



sub log10 
  {
  my $n = shift;
  return log($n)/log(10);
  }
  
  
  
sub prepareDbh
  {
  my $username = "";
  my $password = "";
  my $database = "";
  my $host = "";
  my $dsn = "DBI:mysql:$database:$host";
  our $dbh = DBI->connect($dsn,$username,$password);
  }

