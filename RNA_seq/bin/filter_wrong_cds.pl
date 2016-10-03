#!usr/bin/perl
use strict;

die"Usage:perl $0 <GTF> <WRONG_CDS>\n" if(@ARGV !=2);
my ($gtf, $cds) = @ARGV;

open(CDS, $cds) || die $!;
my %wrong_cds = ();
while(<CDS>)
{
	chomp;
	next if( /ID/ );
	my @line = split("\t");
	$wrong_cds{$line[0]} = 1;
}
close CDS;

open(GTF, $gtf) || die $!;
while(<GTF>)
{
	chomp;
	my @line = split("\t");
	my ($id) = $line[8] =~ /transcript_id "(\S+)"/;
	next if( exists $wrong_cds{$id} );
	print $_."\n";
}
close GTF;
