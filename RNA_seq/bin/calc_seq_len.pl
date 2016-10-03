#!perl -w
use strict;
die "$0 <fasta> [<fasta>...]\n" if @ARGV == 0;
my $len = 0;
while(<>){
	if (/^>(\S+)/){
		print "$len\n" if $len;
		$len = 0;
		print "$1\t";
	}else{
		/\w+/;
		$len += length($&);
	}
}
print "$len\n";
