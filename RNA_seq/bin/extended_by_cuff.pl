#!usr/bin/perl
use strict;

# Author  liugeng   
# Email   liugeng@genomics.cn
# create  2010-12-17
# mender  2011-03-21

die "extend glean to get complete orf using transcripts\n\tperl $0 <gff> <gtf> <loci> <cd_check>\n" if( @ARGV != 4 );
my($gff, $gtf, $loci, $cds_check) = @ARGV;

my %loci = read_loci($loci);
my %gff = read_glean_gff($gff);
my %gtf = read_cuff_gtf($gtf);
my %cds_check = read_cds_check($cds_check);

for my $loci_id  (keys %loci) 
{
	my  $glean_id = $loci{$loci_id}{"glean"};
	my @cuff_id_list = @{$loci{$loci_id}{"cuff"}};
	if( exists $cds_check{$glean_id} ){
		if( $cds_check{$glean_id}[0] == 0 ){
			for my $i ( 0..$#cuff_id_list ){
				my $cuff_id=$cuff_id_list[$i];
				next if(!exists $gtf{$cuff_id});
				@{$gff{$glean_id}{CDS}}=sort{$a->[0] <=> $b->[0]} @{$gff{$glean_id}{CDS}};
				@{$gtf{$cuff_id}{exon}}=sort{$a->[0] <=> $b->[0]} @{$gtf{$cuff_id}{exon}};
				if($gff{$glean_id}{mRNA}[3] eq "+" && $gtf{$cuff_id}{transcript}[3] eq "+"){
					my $glean_start = $gff{$glean_id}{CDS}[0][0];
					my $glean_stop = $gff{$glean_id}{CDS}[0][1];
					my $cuff_start = $gtf{$cuff_id}{exon}[0][0];
					my $cuff_stop = $gtf{$cuff_id}{exon}[0][1];
					if( $cuff_start<$glean_start &&  $cuff_stop==$glean_stop){
						$gff{$glean_id}{s_5} = $cuff_start;
						print STDERR "$glean_id\t$cuff_id\t$gff{$glean_id}{s_5}\n";
						last;
					}
				}
				if( $gff{$glean_id}{mRNA}[3] eq "-" && $gtf{$cuff_id}{transcript}[3] eq "-" ){
					my $glean_start = $gff{$glean_id}{CDS}[-1][1];
					my $glean_stop = $gff{$glean_id}{CDS}[-1][0];
					my $cuff_start = $gtf{$cuff_id}{exon}[-1][1];
					my $cuff_stop = $gtf{$cuff_id}{exon}[-1][0];
					if($cuff_start > $glean_start && $cuff_stop == $glean_stop){
						$gff{$glean_id}{s_5} = $cuff_start;
						print STDERR "$glean_id\t$cuff_id\t$gff{$glean_id}{s_5}\n";
						last;
					}
				}
			}
		}
		if( $cds_check{$glean_id}[1] == 0 ){
			for my $i ( 0..$#cuff_id_list ){
				next if( not exists $gtf{$cuff_id_list[$i]} );
				my $cuff_id=$cuff_id_list[$i];
				if($gff{$glean_id}{mRNA}[3] eq "-" && $gtf{$cuff_id}{transcript}[3] eq "-"){
					my $glean_start = $gff{$glean_id}{CDS}[0][1];
					my $glean_stop = $gff{$glean_id}{CDS}[0][0];
					my $cuff_start = $gtf{$cuff_id}{exon}[0][1];
					my $cuff_stop = $gtf{$cuff_id}{exon}[0][0];
					if( $cuff_start == $glean_start && $cuff_stop < $glean_stop ){
						$gff{$glean_id}{e_3} = $cuff_stop;
						print STDERR "$glean_id\t$cuff_id\t$gff{$glean_id}{e_3}\n";
						last;
					}
				}
				if($gff{$glean_id}{mRNA}[3] eq "+" && $gtf{$cuff_id}{transcript}[3] eq "+" ){
					my $glean_start = $gff{$glean_id}{CDS}[-1][0];
					my $glean_stop = $gff{$glean_id}{CDS}[-1][1];
					my $cuff_start = $gtf{$cuff_id}{exon}[-1][0];
					my $cuff_stop = $gtf{$cuff_id}{exon}[-1][1];
					if( $cuff_start==$glean_start  && $cuff_stop>$glean_stop ){
						$gff{$glean_id}{e_3} = $cuff_stop;
						print STDERR "$glean_id\t$cuff_id\t$gff{$glean_id}{e_3}\n";
						last;
					}
				}
			}
		}
	}
}
foreach my $key(keys %gff){
	@{$gff{$key}{CDS}}=sort{$a->[0] <=> $b->[0]} @{$gff{$key}{CDS}};
	if(exists $gff{$key}{s_5} || exists $gff{$key}{e_3}){
		my $count=0;
		foreach(my $i=0;$i<=$#{$gff{$key}{CDS}};$i++){
			$count+=$gff{$key}{CDS}[$i][1]-$gff{$key}{CDS}[$i]+1;
		}
		my $check=$count % 3;
		if(! $check){
			if($gff{$key}{mRNA}[3] eq '+'){
				$gff{$key}{mRNA}[0]=$gff{$key}{s_5} if(exists $gff{$key}{s_5});
				$gff{$key}{CDS}[0][0]=$gff{$key}{s_5} if(exists $gff{$key}{s_5});
				$gff{$key}{mRNA}[1]=$gff{$key}{e_3} if(exists $gff{$key}{e_3});
				$gff{$key}{CDS}[-1][1]=$gff{$key}{e_3} if(exists $gff{$key}{e_3});
			}elsif($gff{$key}{mRNA}[3] eq '-'){
				$gff{$key}{mRNA}[0]=$gff{$key}{e_3} if(exists $gff{$key}{e_3});
				$gff{$key}{CDS}[0][0]=$gff{$key}{e_3} if(exists $gff{$key}{e_3});
				$gff{$key}{mRNA}[1]=$gff{$key}{s_5} if(exists $gff{$key}{s_5});
				$gff{$key}{CDS}[-1][1]=$gff{$key}{s_5} if(exists $gff{$key}{s_5});
			}
		}
	}
	if($gff{$key}{mRNA}[3] eq '+'){
		print "$gff{$key}{mRNA}[5]\t$gff{$key}{mRNA}[6]\tmRNA\t$gff{$key}{mRNA}[0]\t$gff{$key}{mRNA}[1]\t$gff{$key}{mRNA}[2]\t$gff{$key}{mRNA}[3]\t$gff{$key}{mRNA}[4]\tID=$key;\n";
		my $frame=0;
		foreach(my $i=0;$i<=$#{$gff{$key}{CDS}};$i++){
			$gff{$key}{CDS}[$i][4]=$frame;
			my $len = $gff{$key}{CDS}[$i][1] - $gff{$key}{CDS}[$i][0] + 1;
			$frame = (3-(($len-$frame)%3))%3;
		}
		foreach(my $i=0;$i<=$#{$gff{$key}{CDS}};$i++){
			print "$gff{$key}{mRNA}[5]\t$gff{$key}{mRNA}[6]\tCDS\t$gff{$key}{CDS}[$i][0]\t$gff{$key}{CDS}[$i][1]\t.\t$gff{$key}{mRNA}[3]\t$gff{$key}{CDS}[$i][4]\tParent=$key;\n";
		}
	}
	elsif($gff{$key}{mRNA}[3] eq '-'){
		print "$gff{$key}{mRNA}[5]\t$gff{$key}{mRNA}[6]\tmRNA\t$gff{$key}{mRNA}[0]\t$gff{$key}{mRNA}[1]\t$gff{$key}{mRNA}[2]\t$gff{$key}{mRNA}[3]\t$gff{$key}{mRNA}[4]\tID=$key;\n";
		my $frame=0;
		foreach (my $i=$#{$gff{$key}{CDS}};$i>=0;$i--){
			$gff{$key}{CDS}[$i][4]=$frame;
			my $len = $gff{$key}{CDS}[$i][1] - $gff{$key}{CDS}[$i][0] + 1;
			$frame = (3-(($len-$frame)%3))%3;
		}
		foreach (my $i=$#{$gff{$key}{CDS}};$i>=0;$i--){
		print "$gff{$key}{mRNA}[5]\t$gff{$key}{mRNA}[6]\tCDS\t$gff{$key}{CDS}[$i][0]\t$gff{$key}{CDS}[$i][1]\t.\t$gff{$key}{mRNA}[3]\t$gff{$key}{CDS}[$i][4]\tParent=$key;\n";
		}
	}
}
=head
for my $id ( keys %gff ){
	my $frame = 0;
	my @x = @{$gff{$id}{CDS}} ;
	for my $i ( 0..$#x ){
		$gff{$id}{CDS}[$i][4] = $frame;
		my $len = $gff{$id}{CDS}[$i][1] -  $gff{$id}{CDS}[$i][0] + 1;
		$frame = (3-(($len-$frame)%3))%3;
	}
}
=head
my $count;
open(FH, $gff) || die $!;
while (<FHi>){
	chomp;
	my @tmp = split /\t/;
	my $id;
	if( $tmp[2] eq "mRNA" ){
		$id = $1 if($tmp[8] =~ /ID=([^;]+);/);
		$count=-1
	}
	elsif( $tmp[2] eq "CDS" ){
		$id = $1 if($tmp[8]=~ /Parent=([^;]+);/);
		$count++;
	}
	if(! $check{$id}){
		$tmp[3] = $gff{$id}{CDS}[$count][0];
		$tmp[4] = $gff{$id}{CDS}[$count][1];
		$tmp[7] = $gff{$id}{CDS}[$count][4];
	}
	$num++;
	my $tmp = join("\t", @tmp);
	print $tmp."\n";
}
close FH;
=cut
#==================================
sub read_glean_gff{
	open (FH, shift) || die $!;
	my %gff = ();
	my $num;
	while (<FH>){
		chomp;
		next if (/^#/ || /^\s/);
		my @tmp = split(/\s+/,$_);
		if( $tmp[2] eq "mRNA" ){
			my ($id) = $1  if($tmp[8]=~/ID=([^;]+);/);
			$gff{$id}{mRNA} =($tmp[3] < $tmp[4])?[ $tmp[3], $tmp[4], $tmp[5], $tmp[6], $tmp[7], $tmp[0],$tmp[1] ]:[ $tmp[4], $tmp[3], $tmp[5], $tmp[6], $tmp[7], $tmp[0],$tmp[1] ];
		}
		if( $tmp[2] eq "CDS" ){
			my ($id) = $1 if($tmp[8] =~ /Parent=([^;]+);/);
			if($tmp[3]<$tmp[4]){
				push @{$gff{$id}{CDS}},[ $tmp[3], $tmp[4], $tmp[5], $tmp[6], $tmp[7], $tmp[0] ];
			}else{
				push @{$gff{$id}{CDS}},[ $tmp[4], $tmp[3], $tmp[5], $tmp[6], $tmp[7], $tmp[0] ];
			}
		}
	}
	close FH;
	return %gff;
}
#==================================
sub read_cuff_gtf{
	open (FH, shift) || die $!;
	my %gtf = ();
	my $num;
	while (<FH>){
		chomp;
		next if (/^#/ || /^\s/);
		my @tmp = split(/\t/,$_);
		my ($id) = $tmp[8] =~ /transcript_id "(\S+)"/;
		if( $tmp[2] eq "transcript" ){
			$gtf{$id}{transcript} =($tmp[3] < $tmp[4])?[ $tmp[3], $tmp[4], $tmp[5], $tmp[6], $tmp[7] ]:[ $tmp[4], $tmp[3], $tmp[5], $tmp[6], $tmp[7] ];
		}
		if( $tmp[2] eq "exon" ){
			if($tmp[3]<$tmp[4]){
				push @{$gtf{$id}{exon}},[ $tmp[3], $tmp[4], $tmp[5], $tmp[6], $tmp[7] ];
			}else{
				push @{$gtf{$id}{exon}},[ $tmp[4], $tmp[3], $tmp[5], $tmp[6], $tmp[7] ];
			}
		}
	}
	close FH;
	return %gtf;
}
#==================================
sub read_loci{
	open (FH, shift) || die $!;
	my %loci = ();
	my $loci_id = 0;
	while (<FH>){
		my @tmp = split (/\t/,$_);
		next if($tmp[2] eq '-');
		next unless $tmp[3];
		$loci_id=$tmp[0];
		
		my ($scaff_id, $strand, $start, $stop) = $tmp[1] =~ /(\S+)\[(\S+)](\d+)-(\d+)/;
		$loci{$loci_id}{"scaff"} = $scaff_id;
		$loci{$loci_id}{"strand"} = $strand;
		$loci{$loci_id}{"start"} = $start;
		$loci{$loci_id}{"stop"} = $stop;
		
		for my $ref_id ( $tmp[2] =~ /[^,]+\|([^,]+)/g ){
			$loci{$loci_id}{"glean"} = $ref_id;
		}
		for my $cuff_id(split(",", $tmp[3])){
			push @{$loci{$loci_id}{"cuff"}},$cuff_id;
		}
	}
	close FH;
	return %loci;
}
#=============================
sub read_cds_check{
	open(FH, shift) || die $!;
	my %cds_check = ();
        while (<FH>){
		chomp;
		next if( /ID/ );
		my @tmp = split;
		$cds_check{$tmp[0]} = [$tmp[1], $tmp[2]];
	}
	close FH;
	return %cds_check;
}
