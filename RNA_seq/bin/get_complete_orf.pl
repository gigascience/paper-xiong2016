#!/usr/bin/perl
use strict;

die"perl $0 <GTF> <CDS> <ORF_LIST> <PREFIX>\n" if(@ARGV != 4);
my ($gtf, $cds, $orf, $prefix) = @ARGV;

open(ORF, $orf) || die $!;
my %complete_orf_id = ();
while(<ORF>)
{
	chomp;
	my @line = split("\t");
	$complete_orf_id{$line[0]} = 1 if( $line[4] eq "" );
}
close ORF;

open(GTF, $gtf) || die $!;
open(GTF_C, ">$prefix\.gtf") || die $!;
while(<GTF>)
{
	chomp;
	my @line = split("\t");
	my ($id)  = $line[8] =~ /transcript_id "(\S+)"/;
	next unless exists $complete_orf_id{$id};
	print GTF_C $_."\n";
}
close GTF;
close GTF_C;

open(CDS, $cds) || die $!;
open(CDS_C, ">$prefix\.cds") || die $!;
$/="\n>";
while(<CDS>)
{
	chomp;
	s/^>//;
  my ($tag, $seq) = split(/\n/, $_, 2);
  $seq =~ s/\W+//g;
  $seq =~ s/\n//g;
  my ($id) = $tag =~ /^(\S+)/;
  next unless exists $complete_orf_id{$id};
  print CDS_C ">$tag\n";
  print CDS_C $seq."\n";
}
$/="\n";
close CDS;
close CDS_S;