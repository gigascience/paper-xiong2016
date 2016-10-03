#!/usr/bin/perl -w

=head1 Information
	Author: Bob Zhang       zhangbo@genomics.cn
	        liugeng			    liugeng@genomics.org.cn
	Mender: zhouheling      zhouheling@genomics.org.cn
	create: 2010-11-03
	Updata: 2010-12-16
=cut

use strict;
use FindBin qw($Bin $Script);
use lib "$Bin/../../common_bin";
use GACP qw(parse_config);

die "$0 <outdir> <merged_sorted.bam> [max_intron] [process_num]\n" if (@ARGV < 2);

my ($outdir, $bamfile, $max_intron, $proc_n) =  @ARGV;
die "$bamfile not exists!\n" unless (-f $bamfile);

my $config_file = "$Bin/../../config.txt";
my $samtools = parse_config($config_file,"samtools");
my $prefix = "CUFF";
$max_intron ||= 500000; # 100k by default
$proc_n 	||= 1;

#find output diretory
#if exists then use ready-made data, else make new data
$outdir =~ s/\/$//;
my $workdir = `pwd`;
chomp($workdir);
$workdir =~ s/\/$//;
unless (-d $outdir){
	mkdir($outdir);
	my $n = 1;
	my $l = 1;
	
	# split SAM file into $outdir
	mkdir("$outdir/part.$n.sam");
	open (OUT , ">$outdir/part.$n.sam/part.$n.sam") || die $!;
	open (SAM , '-|' ,"$samtools view $bamfile") || die $!;
	while(<SAM>){
		print OUT;
		if ( $l % 1000000 == 0 )
		{
			my $last_chr = (split /\s+/)[2];
			while(<SAM>){
				if ($last_chr eq (split /\s+/)[2]){
					print OUT;
				}else{
					close OUT;
					$n += 1;
					mkdir("$outdir/part.$n.sam");
					open (OUT , ">$outdir/part.$n.sam/part.$n.sam") || die $!;
                                        print OUT;
                                        $l = 1;
					last;
				}
			}
		}
		$l ++;	
	}
	close SAM;
	close OUT;
}

my $cufflinks = parse_config($config_file,"cufflinks");
my $cuff_para = " -I $max_intron -p $proc_n";

#command line for cufflinks qsub
my $cmdx = "$cufflinks PARAMETERS -L PREFIX -o SPLITDIR/SAMDIR/ SPLITDIR/SAMDIR/SAMDIR";

open FH , '>cufflinks.sh';
opendir(DIR, $outdir) || die $!;
for my $sam (grep (/\.sam$/,readdir DIR )){
	die $sam unless (-d "$outdir/$sam");
	my $cmd = $cmdx;
	my ($index) = $sam =~ /\.(\d+)\.sam$/;
	
	$cmd =~ s/PARAMETERS/$cuff_para/;
	$cmd =~ s/SPLITDIR/$outdir/g;
	$cmd =~ s/SAMDIR/$sam/g;
	$cmd =~ s/PREFIX/$prefix$index/;
	#print FH "cd $outdir/$sam\n";
        print FH 'export LD_LIBRARY_PATH=$Bin/../lib ;';
	print FH $cmd."\n";
	#print FH "cd $workdir\n" 
}
close FH;
closedir DIR;
