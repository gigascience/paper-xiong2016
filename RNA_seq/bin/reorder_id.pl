#!/usr/bin/perl -w
### zhushilin@genomics.org.cn
### 2011-05-19
use strict;

die "Usage: $0  <gff>  <pre>\n" if @ARGV<1;

my ($gff,$pre)=@ARGV;
$pre ||="CCG";

my $gene_num=0;
my $mRNA_num=0;
my $gene_id;
my %change_mid;
my %mRNA;
open GFF,"$gff";
while (<GFF>){
	my @te=split /\t/;
	if ($te[2] eq 'gene'){
		$gene_num++;
		$mRNA_num=0;
		$gene_id=$pre."0" x (6-length($gene_num)).$gene_num ;
		$te[8]=~s/ID=[^;\s]+/ID=$gene_id/;
		my $print=join "\t",@te;
		print $print;
	}elsif($te[2] eq 'mRNA'){
		my ($mid,$source)=$te[8]=~/ID=([^;\s]+);source_id=([^;\s]+);/;
		next if (exists $mRNA{$source});
		$mRNA{$source}=1;
		$mRNA_num++;
		my $mRNA_id=$gene_id.".".$mRNA_num;
		$te[8]=~s/ID=$mid/ID=$mRNA_id/;
		$te[8]=~s/Parent=([^;\s]+);//;
		$change_mid{$mid}=$mRNA_id;
		my $print=join "\t",@te;
		print $print;
	}elsif($te[2] eq 'CDS' || $te[2]=~/UTR/i){
		my ($parent)=$te[8]=~/Parent=([^;\s]+)/;
		next if (!exists $change_mid{$parent});
		$te[8]=~s/Parent=$parent/Parent=$change_mid{$parent}/;
		my $print=join "\t",@te;
		print $print;
	}
}
close GFF;
