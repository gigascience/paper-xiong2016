#!usr/bin/perl
use strict;
use Data::Dumper;
die "Usage:$0 <GFF>\n" if(@ARGV != 1);
my $gff = shift;
open (GFF, $gff) || die $!;
my %gene;my %cac;my $num;my $id;my $cds_len;
while( <GFF> ){
	chomp;
	my @line=split;
	if($line[2] eq "gene"){
		$num=0;
		my @info = split(";", $line[8]);
		($id) = $info[0] =~ /ID=(\S+)/;
	}
	elsif ($line[2] eq "mRNA"){
		$num++;
		$cds_len=0;
#		my $length=$line[4]-$line[3]+1;
		push@{$cac{$id}},$num;
		push @{$gene{$id}{$num}{'txt'}},$_;
#		$gene{$id}{$num}{'len'}=$length;
	}elsif( ($line[2] eq "CDS")){
		push @{$gene{$id}{$num}{'txt'}},$_;
		$cds_len+=$line[4]-$line[3]+1;
		$gene{$id}{$num}{'len'}=$cds_len;
	}elsif(($line[2] eq "UTR_5") || ($line[2] eq "UTR_3") ){
		push @{$gene{$id}{$num}{'txt'}},$_;
	}
}
close GFF;
#print Dumper (\%gene);

foreach my $key(sort keys %cac){
my $t=@{$cac{$key}};
	if ( $t==1){
		for my $i (sort {$a<=>$b} (keys%{$gene{$key}})){
			for my $j ( 0..$#{$gene{$key}{$i}{'txt'}}){
				print  ${$gene{$key}{$i}{'txt'}}[$j]."\n";
			}
		}
	}else{
		my $best_id = "";
		my $len_pre = 0;
		my $len = 0;
			for my $i (sort {$a<=>$b} (keys%{$gene{$key}})){
			$len = $gene{$key}{$i}{'len'};
			if( $len > $len_pre){
				$len_pre = $len;
				$best_id = $i;
			}
		}
		for my $j ( 0..$#{$gene{$key}{$best_id}{'txt'}}){
			print  ${$gene{$key}{$best_id}{'txt'}}[$j]."\n";
		}
	}
}

