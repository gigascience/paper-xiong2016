#!/usr/bin/perl
use strict;
use FindBin qw($Bin $Script);

die "Usage:$0 <ANNOTATION> <GENOME>\n" if(@ARGV != 2);
my ($gff, $fa) = @ARGV;

open(GFF, "<$gff") || die $!;

my %gtf = ();
while (<GFF>)
{
	chomp;
	next if( /^#/ || /^\s/ );
	my @line = split("\t");
	
	# Handle annotation from GeneWise
	if( $line[2] eq "CDS" )
	{
		my ($id)  = $line[8] =~ /Parent=(\S+);/;
		$gtf{$line[0]}{$id}{$line[3]} =[$line[4], $line[6]];
	}
	
	# Handel annotation from Cufflinks
	if( $line[2] eq "exon" )
	{
		my ($id)  = $line[8] =~ /transcript_id "(\S+)"/;
		$gtf{$line[0]}{$id}{$line[3]} = [$line[4], $line[6]];
	}
}
close GFF;

open(FA, "<$fa") || die $!;
open (CDS, ">all.cds.tbl") || die $!;
open (INTRON, ">all.intron.tbl") || die $!;

local $/ = "\n>";
while (<FA>)
{
	chomp;
	s/^>//;
	my ($tag , $seq) = split (/\n/, $_, 2);
	my ($scaff_id)  = $tag =~ /^(\S+)/; 
	$seq =~ s/\W+//g;
	next unless exists $gtf{$scaff_id};
	for my $id (keys %{$gtf{$scaff_id}})
	{
		my @start = (sort{$a<=>$b} keys %{$gtf{$scaff_id}{$id}});
		my $exons = "";
		for(my $i=0; $i<=$#start; ++$i)
		{
			$exons .= substr( $seq,
			          $start[$i]- 1,
			          $gtf{$scaff_id}{$id}{$start[$i]}[0] - $start[$i] + 1 );
		}
                if( $gtf{$scaff_id}{$id}{$start[0]}[1] eq "-")
		{		
		        $exons = reverse($exons);
		        $exons =~ tr/ACGT/TGCA/;
		}
		#$exons =~ s/.{60}/$&\n/g;
		$exons =~ s/\n+$//;
		print CDS "model.$scaff_id.$id\_$scaff_id.$id $exons\n";
		
		for(my $i=0; $i<$#start; ++$i)
		{
			my $intron = substr( $seq,
			                     $gtf{$scaff_id}{$id}{$start[$i]}[0],
			                     $start[$i+1] - $gtf{$scaff_id}{$id}{$start[$i]}[0] - 1 );
			if($gtf{$scaff_id}{$id}{$start[$i]}[1] eq "-" )
			{	
				$intron = reverse($intron);
				$intron =~ tr/ACGT/TGCA/;
			}
			#$intron =~ s/.{60}/$&\n/g;
			$intron =~ s/\n+$//;
			my $ind = $i + 1;
			print INTRON "model.$scaff_id.$id.i$ind\_$scaff_id.$id $intron\n";
		}
	}
}
$/ = "\n"; 
close FA;
close CDS;
close INTRON;

system ("sh $Bin/get_markov_order.sh");
system ("sh $Bin/get_markov5_param.sh");

open (PARA, ">markov_5.param") || die $!;
print PARA "# Parameters to Predict Exons: Markov chains(initial - transition)\n";
print PARA "Markov_oligo_logs_file\n";
print PARA "5\n";
print PARA "\n";

print PARA "Markov_Initial_probability_matrix\n";
open (INI,"<set1.cds-intron.5.initial.geneid") || die $!;
while (<INI>)
{
	print PARA $_;
}
close INI;
print PARA "\n";

print PARA "Markov_Transition_probability_matrix\n";
open (TRAN, "<set1.cds-intron.5.transition.geneid") || die $!;
while (<TRAN>)
{
	print PARA $_;
}
close TRAN;

close PARA;

system ("rm set1.*");
print "Done\n";
