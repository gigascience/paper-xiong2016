#!/usr/bin/perl -w
=head1 Version

	Author:		zhangbo		zhangbo@genomics.cn
	Mender:		zhushilin	zhushilin@genomics.cn
	Mender: 	liugeng		liugeng@genomics.cn
	Mender: 	zhouheling	zhouheling@genomics.cn
	Modified:	2011-10-17
	Modified:	2010-08-25
	Update:		2011-04-02

=head1 Usage

	Generating bash scripts for TopHat mapping using Pair-End RNA-Seq libraries

=head2  liblist format

	FQ1  FQ2  INSERT  SD
	
=cut

use strict;
use FindBin qw($Bin $Script);
use lib "$Bin/../../common_bin";
use GACP qw(parse_config);

die "$0 <FQ_LIST> <OUTDIR> <SCAFF_FASTA> < MAX_INTRON_LENGTH> <MIXMATCH>\n" if (@ARGV != 5);
my ($liblist, $out_dir, $scaff, $max_i_length, $mismatch) = @ARGV;

my $config_file = "$Bin/../../config.txt";
#my $tophat = "/ifs2/PLANT/zhangbo/panfs/tools/tophat-1.0.13/bin/tophat";
#my $samtools = "/ifs2/PLANT/zhangbo/panfs/tools/samtools-0.1.7a/samtools";
my $tophat = parse_config($config_file,"tophat");
my $samtools = parse_config($config_file,"samtools");

#my $tophat_para = "-p 4 --max-intron-length $max_i_length -m $mismatch";

# command line for qsub
my $bowtie = parse_config($config_file,"bowtie");
my $cmd_pe = "export PATH=\"$bowtie/:\$PATH\" ;"
             ."export PATH=\"/opt/blc/genome/biosoft/samtools-0.1.14/:\$PATH\" ;"
             ."$tophat PARAMETERS -o OUTDIR/FQ  SCAFFIDX  RNASEQ;";
#my $cmd_pe = "export PATH=\"$bowtie/:\$PATH\" ;".
#			 "$tophat PARAMETERS -o OUTDIR/FQ  SCAFFIDX  RNASEQ;".
#			 "$samtools view -b -t SCAFFIDX.length -T SCAFFSEQ -o OUTDIR/FQ.bam -S OUTDIR/FQ/accepted_hits.sam;";
#my $cmd_pe = "export PATH=\"/ifs2/PLANT/zhangbo/panfs/tools/bowtie-0.12.5/:\$PATH\" ;".
#             "$tophat PARAMETERS -o OUTDIR/FQ  SCAFFIDX  RNASEQ;".
#             "$samtools view -b -t SCAFFIDX.length -T SCAFFSEQ -o OUTDIR/FQ.bam -S OUTDIR/FQ/accepted_hits.sam;";

$out_dir =~ s/\/$//;
mkdir $out_dir unless ( -e $out_dir);

die "file $scaff not exists!\n" unless ( -f $scaff);   
my ($scaff_idx) = $scaff =~ /(\S+)\.[^\.]+$/;     #bowtie index
my ($pfix) = $out_dir =~ /([^\/]+)$/;    #prefix of bash scripts

# find bowtie index or make one
unless(-f "$scaff_idx.1.ebwt" 
 	&& -f "$scaff_idx.2.ebwt"
  	&& -f "$scaff_idx.rev.1.ebwt"
  	&& -f "$scaff_idx.rev.2.ebwt"){
	open(OUT, ">$out_dir/bowtie-build.sh") || die $!;
	print OUT "$bowtie/bowtie-build -q -f $scaff $scaff_idx";
	close OUT;
	#system ("$bowtie/bowtie-build -q -f $scaff $scaff_idx") == 0 || die $!;
}

# find sequence length list or make one
unless (-f "$scaff_idx.length"){
	system ("perl $Bin/calc_seq_len.pl $scaff > $scaff_idx.length") == 0 || die $!;
}

#make file for qsub
open (my $list, $liblist) || die $!;
open (my $sh , ">tophat.sh") || die $! ;
while (<$list>){
	my $fq;
	my $cmd = $cmd_pe;
	my $id;
	my $read_len;
	my $tophat_para = "-p 4 --max-intron-length $max_i_length -m $mismatch";

	s/\s+$//;
	my @file_te=split;

	# PE (f1 f2 inner_size SD)
	if (@file_te==4){
		my ($fq1, $fq2, $inner, $sd)=@file_te;
		($id) = $fq1 =~ /([^\/]+)_1.fq(?:\.gz)?$/;
		$read_len=length((split("\n", `less $fq1|head -n 2`))[1]);
		$tophat_para .= " -r $inner --mate-std-dev $sd ";
#		$tophat_para .= ' --closure-search ' if $inner<=50;
		$fq="$fq1 $fq2";
	}
	# SE (f)
	elsif(@file_te==1){
		($fq) =@file_te;
		($id) = $fq =~ /([^\/]+).(?:fq|fa|fasta)(?:\.gz)?$/;
		$read_len=length((split("\n", `less $fq|head -n 2`))[1]);
	}else{
		die "Wrong format of list!\nExit!\n";
	}

	$tophat_para .= ' --coverage-search ' if $read_len>=75;
	$tophat_para .= ' --microexon-search ' if $read_len>=50;
	

	$cmd =~ s/PARAMETERS/$tophat_para/;
	$cmd =~ s/OUTDIR/$out_dir/g;
	$cmd =~ s/FQ/$id/g;
	#$cmd =~ s/SCAFFSEQ/$scaff/g;
	$cmd =~ s/SCAFFIDX/$scaff_idx/g;
	$cmd =~ s/RNASEQ/$fq/g;
	
	print $sh "$cmd\n";
	print "$id\n";
	#print "$out_dir/$id/accepted_hits.bam\n"; # print the path of mapping result 	
}
close $list;
close $sh;
