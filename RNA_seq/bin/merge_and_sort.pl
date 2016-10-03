#!/usr/bin/perl -w

# Author: Bob Zhang
# E-mail: zhangbo@genomics.cn
# create: 2010-09-14

use strict;
use FindBin qw($Bin $Script);
use lib "$Bin/../../common_bin";
use GACP qw(parse_config);

die "$0 <bam.list> [<merged_sorted.bam>] \n" if @ARGV < 1;

my $config_file = "$Bin/../../config.txt";
my $samout = $ARGV[1] || 'merged_sorted.bam';
my $samtools = parse_config($config_file,"samtools");
#my $samtools = '/opt/blc/genome/biosoft/samtools-0.1.14/samtools';

open (FH, $ARGV[0]);
my @bam = ();
while (<FH>){
	chomp;
	s/\s+//g;
	push (@bam , $_) if (/bam$/);
}
close FH;

if( @bam > 1)
{
    system( "$samtools merge MERGED.bam ". join("\t", @bam));
}
else
{
    system( "cat $bam[0] > MERGED.bam");
}
system( "$samtools sort -m 1000000000 MERGED.bam SORTED" );
#system( "$samtools sort -k 3,3 -k 4,4n MERGED.bam SORTED" );
system( "mv SORTED.bam $samout" );
system( "rm -f MERGED.bam");
system( "echo Merge and sort... DONE!");

