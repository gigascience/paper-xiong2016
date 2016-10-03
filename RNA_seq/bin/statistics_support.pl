#! /usr/bin/perl 
use strict;
die "Usage:perl $0 <result_gff> \n" if(@ARGV != 1);

my $file=shift;
open IN,"$file"||die "can't open $!";
my $pn;my $in;my $no;my $all;my $source;
while (<IN>){
chomp;
my @line=split;
if($line[2] eq 'mRNA'){
if ($line[1] eq 'GLEAN' ){
my $source=$1 if($line[-1]=~/source_id=(\w+);/);
my $identical=$1 if($line[-1]=~/identical_support_id=(\w+\.\w+\.\w+);/);
if($identical =~/\S+/){
$in++;
}
my $part=$1 if ($line[-1]=~/part_support_id=(\w+\.\w+\.\w+);/);
if($part =~/\S+/){
$pn++;}
if ($identical eq '' && $part eq ''){
$no++;
}
if($identical =~/\S+/ && $part =~/\S+/){
$all++;}
}elsif($line[1] eq 'Cuff' ){
$source++;
}
}
}
close IN;
my $ij=$in-$all;my $pj=$pn-$all;
open OUT1,"> P_I.txt";
print OUT1 "Identical\t$ij\nAll\t$all\nPart\t$pj\nNo\t$no\nCUFF\t$source\n";
close OUT1;
my $R_script=<<RSCRIPT;
pdf(file="CCG.pdf",w=8,h=6)
dat<-read.table("P_I.txt");
ratio<-sprintf("%.2f",100*dat[,2]/sum(dat[,2]));
ratio<-paste(ratio,"%",sep="")
label<-paste(dat[,1],ratio,sep="\n");
pie(dat[,2],col=c("red3","purple3","deepskyblue","royalblue4","yellowgreen"),main="Pie chart for gene sets source",border="gray",labels=label,font=2)
dev.off()
RSCRIPT
open R,">CCG.R";
print R $R_script;
close R;
system("/opt/blc/genome/biosoft/R/bin/R CMD BATCH CCG.R");
