#!/usr/bin/perl -w
use strict;
#use Getopt::Long;
use File::Basename qw(basename dirname);
my $file=shift;
my $file_basename = basename($file);
my $gtf = $file_basename.".gtf";
my ($output,$count,$geneid,$gtfid);
open (IN,"< $file") || die("can't open $!");
while(<IN>){
	chomp;
	my ($scaffold,$soft,$type,$start,$end,$strand,$phase,$idline)=(split(/\s+/),$_)[0,1,2,3,4,6,7,8];
	if ($type eq "mRNA"){
		$count=1;
		$geneid=$1 if ($idline=~/ID=([^;]+);/);
		$gtfid=$geneid.".1";
		$output.="$scaffold\t$soft\t$type\t$start\t$end\t0\t$strand\t$phase\tgene_id \"$geneid\"; transcript_id \"$gtfid\";\n";
	}elsif($type eq "CDS"){
		#$geneid=$1 if ($idline=~/Parent=([^;]+);/);#change it if mrna have the different id with cds
		$output.="$scaffold\t$soft\t$type\t$start\t$end\t0\t$strand\t$phase\tgene_id \"$geneid\"; transcript_id \"$gtfid\"; exon_number \"$count\";\n";
		$count++;
	}
}
close IN;
open (OUT,"> $gtf") || die("$!\n");
print OUT $output;
close OUT;

