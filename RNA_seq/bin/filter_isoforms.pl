#!/usr/bin/perl -w
### zhushilin@genomics.org.cn
use strict;
use Fatal;
use File::Basename;

die "Usage: $0 <gff>\n" if @ARGV<1;

my $gff=shift @ARGV;
my $pre=basename $gff;

my %len;
my %score;
my %infor;
open GFF,"$gff";
while (<GFF>){
	next if (/^#/ || /^\s/);
	s/Cufflinks/Cuff/;
	my @te=split /\t/;
	if ($te[2] eq 'mRNA'){
		my ($Rid)=$te[8]=~/ID=([^;]+)/;
		my $score=($te[1]=~/Cuff/i)? $te[5] : $te[5]*1000;
		$score{$Rid}=$score;
		push @{$infor{$Rid}},$_;
	}elsif($te[2] eq 'CDS'){
		$te[8]=~/Parent=(\w+)\.(\d+)/;
		my $Gid=$1;
		my $Rid="$1.$2";
		my $cds_len=$te[4]-$te[3]+1;
		$len{$Gid}{$Rid}+=$cds_len;
		push @{$infor{$Rid}},$_;
	}elsif($te[2]=~/UTR/i){
		my ($Rid)=$te[8]=~/Parent=([^;]+)/;
		push @{$infor{$Rid}},$_;
	}
}
close GFF;

#open OUT,">$pre.filter";
foreach my $gene (sort keys %len){
	my @sort_id = sort {$len{$gene}{$b} <=> $len{$gene}{$a}} keys %{$len{$gene}};
	my @keep_te;
	my $keep_one=shift @sort_id;
	push @keep_te,$keep_one;
	while(1){
		my $next_one=shift @sort_id;
		last unless ($next_one);
		last if $len{$gene}{$next_one}<$len{$gene}{$keep_one};
		push @keep_te,$next_one;
	}
	my @sort_score=sort {$score{$b} <=> $score{$a}} @keep_te;
	$keep_one=shift @sort_score;
#	print OUT @{$infor{$keep_one}};
        print @{$infor{$keep_one}};
}
#close OUT;
