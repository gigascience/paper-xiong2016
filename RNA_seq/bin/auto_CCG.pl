#!/usr/bin/perl
=head1  Name

	auto_CCG.pl -- the pipeline of combine Cufflinks and Glean.

=head1 example

	perl auto_CCG.pl -max_i_len 50000 -min_cds_len 300 -step 123456 genome.fa fq.list para.txt glean.gff

=head1 Information
	
	Author: 	liugeng     		liugeng@genomics.org.cn
	Mender: 	zhouheling          zhouheling@genomics.org.cn
	Mender:     liugeng             liugeng@genomics.org.cn
	Create:	  	2010-11-04
	Updata:		2011-04-02

=head1 Usage

	perl auto_CCG.pl [options] genome_file  fq_list  parameter_file  Glean_annotation_file
	-run <str(qsub/local)>  set running maner of bowtie-build, default local
	-queue <str>            set the queue in qsub, default no
	-pro_code <str>         set the code for project,default no
        -node <str>          	set the compute node,default h=compute-0-190
	-bowtie_vf <str>        set qsub resource in running bowtie-build if -run is qsub, default vf=4G
	-tophat_vf <str>        set qsub resource in running tophat, default vf=5G
    -qual_system <str>      set the quality sysetem type, 33 or 64
	-max_i_len <int>	set max intron length when runing TopHat and Cufflinks, default 500000
	-mismatch <int>		set mismatch number in every mapped reads, mismatch must = 0, 1, 2; default 0
	-cuff_vf <str>          set qsub resource in running Cufflinks, default vf=5G
	-cuff_proc_num <int>	set process number when running Cufflinks, default 1
	-min_cds_len <int>	set min cds len in transcripts ORF
	-orf <str>              set flag of orf, if T, then request complete orf of glean, if F, then not, and try to EXTEND to get a complete orf using RNA-Seq; difault T
	-step <str>		set the start step for running, default 123456
	-species <str>		set the species' name for the prefix of mRNA ID, default CCG
	-help			output help information to screen.

=cut

use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use lib "$Bin/../../common_bin";
use GACP qw(parse_config);

my ($run, $Node,$bowtie_vf, $tophat_vf, $qual_system, $max_i_len, $mismatch, $cuff_vf, $cuff_proc_n, $min_cds_len, $orf, $step, $help,$Queue,$Pro_code);
my $species;
GetOptions
(
	"run:s"=>\$run,
	"queue:s"=>\$Queue,
	"pro_code:s"=>\$Pro_code,
        "node:s"=>\$Node,
	"bowtie_vf:s"=>\$bowtie_vf,
	"tophat_vf:s"=>\$tophat_vf,
    "qual_system:s"=>\$qual_system,
	"max_i_len:i"=>\$max_i_len,
	"mismatch:i"=>\$mismatch,
	"cuff_vf:s"=>\$cuff_vf,
	"cuff_proc_num:i"=>\$cuff_proc_n,
	"min_cds_len:i"=>\$min_cds_len,
	"orf:s"=>\$orf,
	"step:s"=>\$step,
	"species:s"=>\$species,
	"help"=>\$help
);
die `pod2text $0` if($help);
die "Usage: $0  <genome.fa>  <RNAseq-fq.list>  <markov_5.param>  <Glean.gff> [--options]
You can type '--help' to read detailed information about this program.\n" if (@ARGV<2);
my ($fasta, $fq, $para, $gff) = @ARGV;

$run ||= "local";
#$Node ||= "h=compute-0-190";
$bowtie_vf ||= "vf=4G"; $bowtie_vf="vf=$bowtie_vf" unless ($bowtie_vf=~/vf=/);
$tophat_vf ||= "vf=5G"; $tophat_vf="vf=$tophat_vf" unless ($tophat_vf=~/vf=/);
$max_i_len ||= 500000;
$mismatch ||= 0;
$cuff_vf ||= "vf=5G";   $cuff_vf="vf=$cuff_vf" unless ($cuff_vf=~/vf=/);
$cuff_proc_n ||= 1;
$min_cds_len ||= 150;
$orf ||= "T";
$step ||= 123456;
$species ||="CCG";
$qual_system ||="64";
my $config_file = "$Bin/../../config.txt";
my $qsub_sge = parse_config($config_file,"qsub_sge.pl");
my $cds2aa = parse_config($config_file,"cds2aa.pl");
my $getGene = parse_config($config_file,"getGene.pl");
my $auto_tophat_PE = "$Bin/auto_tophat_PE_and_SE.pl";
my $merge_and_sort = "$Bin/merge_and_sort.pl";
my $auto_cufflinks = "$Bin/auto_cufflinks.pl";
my $predict_orf = "$Bin/predict_orf_V2.pl";
my $get_complete_orf = "$Bin/get_complete_orf.pl";
my $get_cds_from_cufflinks = "$Bin/get_cds_from_cufflinks.pl";
my $filter_wrong_cds = "$Bin/filter_wrong_cds.pl";
my $filter_wrong_cds_middle_triple ="$Bin/filter_wrong_cds_middle_triple";
my $get_gtf_UTR="$Bin/get_gtf_UTR.pl";
my $compare = "$Bin/compare.pl";
my $extended_by_cuff = "$Bin/extended_by_cuff.pl";
my $blast_in_tmap = "$Bin/blast_in_tmap.pl";
my $combine_cuff_glean = "$Bin/combine_Cuff_Glean.pl_UTR.change";
my $rebuild="$Bin/rebuild_gff.pl";
my $reorder="$Bin/reorder_id.pl";
my $filter_isoforms = "$Bin/filter_isoforms.pl";

##add by luchx
my $QP_para;
$QP_para.="--queue $Queue " if (defined $Queue);
$QP_para.="--pro_code $Pro_code " if (defined $Pro_code);

my $workdir=`pwd`;
chomp $workdir;

# step 1: run tophat
if( $step =~ /1/ )
{
	system("perl $auto_tophat_PE $fq \"./\" $fasta $max_i_len $mismatch $qual_system >list") == 0 || die $!;
	if (-e "bowtie-build.sh")
	{
		if( $run eq "qsub" )
		{
			system("perl $qsub_sge $QP_para --reqsub --lines 1 --resource $bowtie_vf --convert no  bowtie-build.sh") == 0 || die $!;
		}
		else
		{
			system("sh bowtie-build.sh") == 0 || die $!;
		}
	}
	my @pre_dir=`ls $workdir`;
	for my $i ( 0..$#pre_dir){
		if ( $pre_dir[$i]=~/tophat.sh.\d+.qsub/ || $pre_dir[$i]=~/tophat.sh.\d+.log/){
			system("rm -r $pre_dir[$i]") || die $!;
		}
	}
	system("perl $qsub_sge $QP_para --reqsub --lines 1 --resource $tophat_vf --convert no  tophat.sh") == 0 || die $!;
	my @dir = `ls $workdir`;
	my $qsub_dir;
	for my $i ( 0..$#dir)
	{
		if( $dir[$i] =~ /tophat.sh.\d+.qsub/ )
		{
			$qsub_dir = $dir[$i];
  			chomp $qsub_dir;
  			last;
		}
	}
        open(FH, "list") || die $!;
        open(OUT, ">bam.list") || die $!;
        while(<FH>)
        {
                chomp;
                print OUT "$workdir/$qsub_dir/$_/accepted_hits.bam\n";
        }
        close FH;
        close OUT;
        system("rm list") == 0 || die $!;
	#@dir = `ls $workdir/$qsub_dir`;
	#for my $i ( 0..$#dir)
	#{
	#	if( $dir[$i] =~ /.bam/ )
	#	{
	#		chomp $dir[$i];
 	#		system("cp $workdir/$qsub_dir/$dir[$i] $workdir") == 0 || die $!;
	#		system("rm $workdir/$qsub_dir/$dir[$i]") == 0 || die $!;
	#	}
	#}
}

# step 2: merge and sort sam file (tophat result)
if( $step =~ /2/ )
{
	system("perl $merge_and_sort bam.list merged_sorted.bam") == 0 || die $!;
}

# step 3: run Cufflinks
if( $step =~ /3/ )
{
	if( ($step !~ /2/) and (not -e "merged_sorted.bam") )
	{
                open(FH, "bam.list") || die $!;
                while(<FH>)
                {
                          chomp;
		          system("cat $_ >>merged_sorted.bam") == 0 || die $!;
                }
                close FH;          
	}
	system("perl $auto_cufflinks \"$workdir/cufflinks\" merged_sorted.bam $max_i_len $cuff_proc_n") == 0 || die $!;
	#system("sh cufflinks.sh") == 0 || die $!;
	#system("cat $workdir/cufflinks/part.*.sam/transcripts.gtf >$workdir/transcripts.gtf") == 0 || die $!;
	system("perl $qsub_sge $QP_para --reqsub --lines 1 --resource $cuff_vf --convert no  cufflinks.sh") == 0 || die $!;
	system("cat $workdir/cufflinks/part.*.sam/transcripts.gtf >$workdir/transcripts.gtf") == 0 || die $!;
}

# step 4: predict transcripts orf
if( $step =~ /4/ )
{
	my $gtf = "transcripts.gtf";
	system( "perl $predict_orf $fasta $gtf $para $min_cds_len 1 >info.log" ) == 0 || die $!;
	system( "perl $get_cds_from_cufflinks $fasta 1.orf.gtf 1.orf" ) == 0 ||die $!;
	system( "perl $cds2aa 1.orf.cds >transcripts.orf.pep" ) == 0 || die $!;
	system( "perl $get_complete_orf 1.orf.gtf 1.cds 1.orf 2" ) == 0 || die $!;
	system( "perl $get_cds_from_cufflinks $fasta 2.gtf 3" ) == 0 ||die $!;
	system( "perl $cds2aa -check 3.cds >cds.wrong" ) == 0 || die $!;
	system( "perl $filter_wrong_cds 2.gtf cds.wrong > transcripts.complete.orf.gtf" ) == 0 || die $!;

	rename( "1.orf.gtf", "transcripts.orf.gtf" ) || die $!;
	system( "rm 1.*" ) ==0 || die $!;
	system( "rm 2.*" ) ==0 || die $!;
	system( "rm 3.*" ) ==0 || die $!;
	system( "rm cds.wrong" ) ==0 || die $!;
	system( "rm info.log" ) ==0 || die $!;

	system( "perl $get_gtf_UTR transcripts.gtf transcripts.complete.orf.gtf > transcripts.complete.orf.gtf.UTR") ==0 || die $!;
	system( "perl $get_gtf_UTR transcripts.gtf transcripts.orf.gtf > transcripts.orf.gtf.UTR") ==0 || die $!;
	
}

# step 5: run cuffcompare
if( $step =~ /5/ )
{
	my $glean_gff = $gff;
	my $cuff_gtf = "$workdir/transcripts.orf.gtf.UTR";
	my $cuff_pep = "$workdir/transcripts.orf.pep";
	my $complete_orf = "$workdir/transcripts.complete.orf.gtf";
  
	system( "perl $compare $glean_gff $cuff_gtf" ) == 0 || die $!;
	
	if( $orf eq "F" )
	{
		my $basename = basename($gff);
		system( "perl $getGene $gff $fasta >glean.cds" ) == 0 || die $!;
		system( "perl $cds2aa -check glean.cds >cds.wrong" ) == 0 || die $!;
		system( "perl $extended_by_cuff $gff $complete_orf stdout.loci cds.wrong >$workdir/$basename.extend.gff" ) == 0 || die $!;
		system( "rm cds.wrong" ) == 0 || die $!;
		system( "rm glean.cds" ) == 0 || die $!;
		$glean_gff = "$workdir/$basename.extend.gff";
	}
	my $basename = basename($glean_gff);
	system( "perl $getGene $glean_gff $fasta >$workdir/$basename.cds" ) == 0 || die $!;
	system( "perl $cds2aa $workdir/$basename.cds > $workdir/$basename.cds.pep" ) == 0 || die $!;
	my $glean_pep = "$workdir/$basename.cds.pep";
	
	system( "perl $blast_in_tmap $glean_pep $cuff_pep $workdir/stdout.tmap" ) ==0 ||die $!;
	system( "perl $qsub_sge $QP_para --reqsub --lines 1 --resource vf=1G --convert no --maxjob 50  compare.sh") == 0 || die $!;
	system( "cat $workdir/stdout.tmap.split/*.match.pep >match.pep" ) == 0 || die $!;
}

# step 6: combine Cufflinks and Glean
if( $step =~ /6/ )
{
	my $glean_gff = $gff;
	if( $orf eq "F" )
	{
		my $basename = basename($gff);
		$glean_gff = "$workdir/$basename.extend.gff";
	}
	my $complete_orf = "transcripts.complete.orf.gtf.UTR";
	my $transcripts = "transcripts.gtf";

	my $front_pre;
	if ($species=~/([a-zA-Z])[a-zA-Z]*[\-_\.]([a-zA-Z]{2})\w*/){
		$front_pre=(uc $1).(lc $2);
	}elsif($species=~/(\w)(\w{2})/){
		$front_pre=(uc $1).$2;
	}else{
		$front_pre=$species;
	}

	system( "perl $combine_cuff_glean $glean_gff $complete_orf $transcripts stdout.loci match.pep $front_pre >CCG.gff" ) == 0 || die $!;
	system( "perl $rebuild CCG.gff --out CCG.rebuild.gff --verbose > rebuild.log") ==0 || die $!;
	system( "perl $reorder CCG.rebuild.gff $front_pre > CCG.rebuild.gff.reorder") ==0 || die $!;
	system( "perl $filter_isoforms CCG.rebuild.gff.reorder > CCG.rebuild.gff.reorder.noisoforms.gff" ) == 0 || die $!;
	system( "perl $getGene CCG.rebuild.gff.reorder.noisoforms.gff $fasta >$workdir/CCG.rebuild.gff.reorder.noisoforms.cds" ) == 0 || die $!;
	system( "perl $cds2aa $workdir/CCG.rebuild.gff.reorder.noisoforms.cds > $workdir/CCG.rebuild.gff.reorder.noisoforms.pep" ) == 0 || die $!;
}
