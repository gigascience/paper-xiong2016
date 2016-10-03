#!usr/bin/perl
# Author: liugeng	liugeng@genomics.cn
# Create: 2010-11-09

use strict;
use FindBin qw($Bin $Script);

die "Usage:perl $0 <glean_pep> <cuff_pep> <tmap>\n" if( @ARGV != 3 );
my ($glean_pep, $cuff_pep, $tmap) = @ARGV;

my $_compare_tmap_seqs = "$Bin/_compare_tmap_seqs.pl";

split_file($tmap, 1000);

opendir(DIR, "$tmap.split") || die $!;
open(FH, ">compare.sh") || die $!;
for (readdir DIR)
{
	my ($id) = /(\d+)$/;
	next unless $id;
	print FH "perl $_compare_tmap_seqs $glean_pep $cuff_pep $tmap.split/$_\n";
}

close FH;

#=========================
sub split_file
{
	my ($file, $n) = @_;
	my ($base_name) = $file =~ /([^\/]+)$/;
	my $out_dir = $file.'.split';
	mkdir $out_dir || undef;	
	my $i = 0;
	open (FH, $file) || die $!;
	open (OUT, ">$out_dir/$base_name.1") || die $!;
	while (<FH>)
	{
		$i += 1;		
		if ($i % $n == 0)
		{
			close OUT;
			open (OUT, ">$out_dir/$base_name.".int($i/$n + 1)) || die $!;
		}
		print OUT;
	}	
}
#=========================
