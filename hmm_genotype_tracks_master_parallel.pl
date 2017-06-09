#!/usr/bin/perl -w
#David Turissini


use DBI;
prepareDbh();

my ($prefix, $suffix, $hmm_type, $recip_group, $repeats_table) = @ARGV;

my $ind_query = $dbh->prepare(qq{select ind
                                 from tracking_groups
                                 where group_name = '$recip_group'});
$ind_query->execute();


while (my ($ind) = $ind_query->fetchrow_array())
  {
  $ind =~ tr/\./\_/;
  my $err_file = 'logs/err_hmm_genotype_tracks_' . $hmm_type . '_' . $suffix . '_' . $ind . '.%J';
  my $out_file = 'logs/out_hmm_genotype_tracks_' . $hmm_type . '_' . $suffix . '_' . $ind . '.%J';
  
  my $hmm_table = $prefix . '_' . $ind . '_hmm_' . $hmm_type . '_pab_01_9_25000_' . $suffix;
  
  system(qq{bsub -e $err_file -o $out_file python make_hmm_genotype_tracks.py $hmm_table $repeats_table});
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

