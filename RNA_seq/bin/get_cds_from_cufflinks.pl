#!/usr/bin/perl -w
use strict;

# Author: Liu geng
# E-mail: liugeng@genomics.cn
# create: 2010-09-09

die "usage: perl $0 <FASTA> <GTF> <PREFIX>\n" if (@ARGV != 3);

my ($fasta,$gtf, $prefix) = @ARGV; 

my %gtf = ();
open(GTF, $gtf)|| die $!;
while (<GTF>){
	next if (/^#/ || /^\s/);
	my @case = split("\t"); 	
	if($case[2] eq 'exon')
	{
		my ($id)  = $case[8] =~ /transcript_id "(\S+)"/;
		$gtf{$case[0]}{$id}{$case[3]} = [ $case[4],$case[6] ];
	}
}
close GTF;

open (FA , $fasta) || die $!;
open (EXON, ">$prefix\.cds") || die $!;
local $/ = "\n>";
while (<FA>){
	chomp;
	s/^>//;
	my ($tag , $seq) = split (/\n/, $_, 2);
	my ($scaff_id)  = $tag =~ /^(\S+)/; 
	$seq =~ s/\W+//g;
	next unless exists $gtf{$scaff_id};
	my $strand;
	for my $transcript_id (keys %{$gtf{$scaff_id}}){
		my $exons = '';
		
		for my $start (sort{$a<=>$b} keys %{$gtf{$scaff_id}{$transcript_id}}){
			$exons .= substr(	$seq,
								$start - 1,
								$gtf{$scaff_id}{$transcript_id}{$start}[0] - $start + 1) ;
		$strand = $gtf{$scaff_id}{$transcript_id}{$start}[1];
		}
		
		print EXON ">$transcript_id  exons\n";
		#$exons =~ s/.{60}/$&\n/g;
		$exons =~ s/\n+$//;
		if( $strand eq "-" )
		{
			$exons = reverse($exons);
			$exons =~ tr/ACGT/TGCA/;
		}
		print EXON "$exons\n";
	}
}
$/ = "\n"; 
close FA;
close EXON;


