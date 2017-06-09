# Int-HMM

Introduction
	This readme describes the HMM used for identifying genomic introgressions described in “Fine scale mapping of genomic introgressions within the Drosophila yakuba clade”, Turissini and Matute (submitted).  The programs are presented as is and are not being actively supported. Scripts will need to be edited to contain connection information (user, password, database, hostname) for MySQL. Also, several sections are hardcoded to use the 5 major chromosome arms found in Drosophila (2L, 2R, 3L, 3R, X). If other chromosome names are required, the programs will need to be edited accordingly.

Prerequisites
Files
o	  Filtered vcf file containing genotype calls from GATK for individuals from both the donor and recipient species

Perl
o	  Modules
•	    DBI 
•	    DBD::mysql 

Python 2
o 	Modules
•	    __future__
•	    datetime
•   	random
•	    math
•	    sys
•	    MySQLdb
•	    os
•	    warnings
•	    Cluster running LSF

MySQL
o	All scripts will need to be edited to contain connection information
o	The following tables need to be created
•	introgression_snps_params
+--------------+--------------+------+-----+---------+-------+
| Field        | Type         | Null | Key | Default | Extra |
+--------------+--------------+------+-----+---------+-------+
| run_name     | varchar(30)  | NO   | PRI | NULL    |       | 
| group_list   | varchar(255) | YES  |     | NULL    |       | 
| arms_table   | varchar(30)  | YES  |     | NULL    |       | 
+--------------+--------------+------+-----+---------+-------+

•	run_name is a descriptive name for the hmm run
•	group_list is a string containing the names of the group_names from tracking_groups that will be part of the run. The list should be in the form 'A', 'B', 'C'
•	arms_table contains the name of the table $arms_table described below



•	tracking_bams
+-----------------+---------------+------+-----+---------+----------------+
| Field           | Type          | Null | Key | Default | Extra          |
+-----------------+---------------+------+-----+---------+----------------+
| tb_id           | int(11)       | NO   | PRI | NULL    | auto_increment | 
| ind             | varchar(100)  | YES  | MUL | NULL    |                | 
| bam_file        | varchar(255)  | YES  |     | NULL    |                | 
| cov_99_quantile | int(11)       | YES  |     | NULL    |                | 
+-----------------+---------------+------+-----+---------+----------------+

•	tb_id is an auto_increment integer
•	ind is the name of the individual
•	bam_file is the path to the bam file containing the mapping results for ind but can be any unique name
•	cov_99_quantile is the 99th quantile of the per site coverage distribution from the bam file



•	tracking_groups
+------------+--------------+------+-----+---------+----------------+
| Field      | Type         | Null | Key | Default | Extra          |
+------------+--------------+------+-----+---------+----------------+
| tg_id      | int(11)      | NO   | PRI | NULL    | auto_increment | 
| ind        | varchar(100) | YES  | MUL | NULL    |                | 
| group_name | varchar(50)  | YES  | MUL | NULL    |                | 
+------------+--------------+------+-----+---------+----------------+

•	tg_id is an auto_increment integer
•	ind is the name of the individual
•	group_name is the name of the group the ind belongs to. The group names should refer to the donor and recipient species



•	tracking_ref_genomes
+----------------------------+--------------+------+-----+---------+-------+
| Field                      | Type         | Null | Key | Default | Extra |
+----------------------------+--------------+------+-----+---------+-------+
| ref_genome                 | varchar(255) | NO   | PRI | NULL    |       | 
| ref_db_table               | varchar(50)  | YES  |     | NULL    |       | 
| chrom_name_field           | varchar(20)  | YES  |     | NULL    |       | 
| chrom_len_field            | varchar(20)  | YES  |     | NULL    |       | 
+----------------------------+--------------+------+-----+---------+-------+

•	ref_genome can be the path to a file containing the reference genome or just a unique name
•	ref_db_table should be the name of the table described below as $arms_table
•	chrom_name_field is the name of the chromosome name field in $arms_table which is arm in the example below
•	chrom_len_field is the name of the field containing the chromosome length in $arms_table



•	tracking_vcf_bams
+----------+--------------+------+-----+---------+----------------+
| Field    | Type         | Null | Key | Default | Extra          |
+----------+--------------+------+-----+---------+----------------+
| tvb_id   | int(11)      | NO   | PRI | NULL    | auto_increment | 
| vcf_file | varchar(255) | YES  |     | NULL    |                | 
| ind      | varchar(30)  | YES  | MUL | NULL    |                | 
| bam_file | varchar(255) | YES  |     | NULL    |                | 
+----------+--------------+------+-----+---------+----------------+

•	tvb is an auto_increment integer
•	vcf_file is the name of the vcf file passed as a parameter into get_site_all_introgress_pop_from_vcf.pl
•	ind is the individual name
•	bam_file should be the same bam_file associated with ind in tracking_bams. A record should be present for every ind in the vcf file.



•	$repeat_regions_table
•	Contains locations of repetitive regions of the genomes which are used for filtering spurious introgressions
+-------+------------+------+-----+---------+----------------+
| Field | Type       | Null | Key | Default | Extra          |
+-------+------------+------+-----+---------+----------------+
| rr_id | int(11)    | NO   | PRI | NULL    | auto_increment | 
| arm   | varchar(5) | YES  |     | NULL    |                | 
| start | int(11)    | YES  | MUL | NULL    |                | 
| end   | int(11)    | YES  | MUL | NULL    |                | 
+-------+------------+------+-----+---------+----------------+

•	rr_id is an auto_increment integer
•	arm is the chromosome name
•	start is the start of the repetitive region
•	end is the end of the repetitive region



•	$arms_table
•	contains the chromsomes present in the genome
+---------+------------+------+-----+---------+-------+
| Field   | Type       | Null | Key | Default | Extra |
+---------+------------+------+-----+---------+-------+
| arm     | varchar(6) | NO   | PRI | NULL    |       | 
| arm_len | int(11)    | YES  |     | NULL    |       | 
+---------+------------+------+-----+---------+-------+
•	arm is the chromosome name and any name can be used as long as the field name is also in tracking_ref_genomes
•	arm_len is the chromosome length and any name can be used as long as it is also in tracking_ref_genomes. 

Running the HMM pipeline

1.	get_site_all_introgress_pop_from_vcf.pl – This program identifies the markers to be used by the HMM. Markers are stored in a separate table for each arm. It takes the following input parameters:
a.	$prefix – the prefix for the run. Note that this should not exceed 7 characters (and maybe less if there are long ind names) such that database tables created by later programs do not exceed the MySQL length limit.
b.	$run_name – the name of the run_name in the introgression_snps_params table that contains the group_list and arms_table for this run.
c.	$donor_group – The group from which the introgression is coming from.
d.	$recip_group – The group into which the introgression is going.
e.	$freq_cutoff – The allele frequency cutoff for markers such that abs(donor_freq – recip_freq) >= $freq_cutoff
f.	$repeat_regions_table – database table containing the unique repetitive regions for the ref genome the reads were mapped to.
g.	$gatk_vcf – vcf file containing genotypes for the groups present in the group_list from introgression_snps_params.

2.	hmm_introgress_all_downsample_single_arm_master.pl – The master program submits individual jobs to run the HMM for each arm/ind combination on a cluster using bsub. The HMM results for each ind are stored in their own database table. The program takes the following input parameters:
a.	$prefix – same prefix as 1a.
b.	$recip_group – same recip_group as 1d.
c.	$pab – sequencing error parameter for hmm: use .01.
d.	$base_a_prob – per base recombination rate used for transition probabilities: use 1e-9
e.	$err_multiplier - multiplier for per base recombination rate used for transition probabilities involving error states: use 25000
f.	$freq_diff_cutoff – same as 1e.

3.	hmm_genotype_tracks_master_parallel.pl – Tracks are contiguous blocks of SNPs with the same most probable genotype as determined by the HMM. This master program submits a job for each ind on a cluster using bsub. The tracks program makes unfiltered and filtered tracks tables as well as a filter log table that is only used for debugging purposes. The program takes the following input parameters:
a.	$prefix - same prefix as 1a.
b.	$suffix – suffix identifying how markers were chosen: use all_30 
c.	$hmm_type - a 3 character name identifying the hmm program that was run. Will be ‘dsm’ most of the time.
d.	$recip_group – same as 1d.
e.	$repeats_table - same table as 1f.

4.	make_hmm_introgress_regions.pl – The regions program should be run after all of the tracks jobs have completed. It creates a series of summary tables that contain information for all of the inds. The unftrk table contains info for all unfiltered tracks. The tracks table contains all of the filtered tracks. The regions table contains all of the regions defined as a set of overlapping tracks. The regions_inds table contains all of the inds with an introgressed track for each region. The blocks table contains blocks which are regions of the genome with the same number of inds with introgressed tracks. The blocks_inds table contains all of the inds with introgressed tracks for each block. The program takes the following input parameters:
a.	$prefix - same prefix as 1a.
b.	$introgress_group – same as 1d.
c.	$hmm_suffix – combines info from pab, err_multiplier, and suffix: use 01_9_25000_all_30 if using the suggested parameters
d.	$hmm_type – same as 1c.
e.	$repeats_table - same table as 1f.


Interpreting the results
	Results can be found in the consolidated tracks table created by make_hmm_introgress_regions.pl with the name [prefix]_hmm_dsm_pab_[hmm_suffix]_tracks. The genotype field will contain 3 values:
•	homo_1 – homozygous for the recipient species allele
•	het – heterozygous for the donor and recipient species alleles indicative of a heterozygous introgression
•	homo _2 – homozygous for the donor species allele indicative of a homozygous introgression
introgress_snps is the count of snps within the track that appear to be introgressed. per_repeat is the percentage of the track that overlaps with repetitive sequence.

