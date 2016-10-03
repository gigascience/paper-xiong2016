#!usr/bin/perl 
use strict;

# Author: liugeng	liugeng@genomics.cn
# Create: 2010-11-09

die "Usage:perl $0 <glean_gff> <cuff_gtf>\n" if(@ARGV != 2);
my ($gff, $gtf) = @ARGV;

my %gff = read_glean($gff);
my %gtf = read_cuff($gtf);

open(LOCI, ">stdout.loci") || die $!;
open(TMAP, ">stdout.tmap") || die $!;
my $loci_num = 0;

for my $scaff_id ( keys %gff )
{
	my (@glean_id, @glean_start, @glean_stop);
	if( exists $gff{$scaff_id} )
	{
	  @glean_id = @{$gff{$scaff_id}{"id"}};
	  @glean_start = @{$gff{$scaff_id}{"start"}};
	  @glean_stop = @{$gff{$scaff_id}{"stop"}};
	}
	
	my (@cuff_id, @cuff_start, @cuff_stop, @cuff_tag);
	if( exists $gtf{$scaff_id} )
	{
	  @cuff_id = @{$gtf{$scaff_id}{"id"}};
	  @cuff_start = @{$gtf{$scaff_id}{"start"}};
	  @cuff_stop = @{$gtf{$scaff_id}{"stop"}};
	}
	for my $i ( 0..$#cuff_id )
	{
		$cuff_tag[$i] = 0;
	}
	
	my ($loci, $gene_start, $gene_stop,);	
	
	$loci = "";
	for my $i ( 0..$#glean_id )
	{
		$gene_start = $glean_start[$i];
		$gene_stop = $glean_stop[$i];
		$loci = "XLOCI_"."0" x (6-length($loci_num)).$loci_num."\t"."$scaff_id\[.\]START-STOP\t$glean_id[$i]|$glean_id[$i]\n";
		$loci_num++;
		for my $j ( 0..$#cuff_id )
		{
#			next if( $cuff_tag[$j] == 1 );
			if( (($cuff_start[$j]>=$glean_start[$i])&&($cuff_start[$j]<=$glean_stop[$i])) ||
			    (($cuff_stop[$j]>=$glean_start[$i])&&($cuff_stop[$j]<=$glean_stop[$i])) ||
			    (($cuff_start[$j]>=$glean_start[$i])&&($cuff_stop[$j]<=$glean_stop[$i])) ||
			    (($cuff_start[$j]<=$glean_start[$i])&&($cuff_stop[$j]>=$glean_stop[$i])) )
			{
				$cuff_tag[$j] = 1;
				$loci =~ s/\n/,$cuff_id[$j]\n/;
				$gene_start = $cuff_start[$j] if( $gene_start > $cuff_start[$j] );
				$gene_stop = $cuff_stop[$j] if( $gene_stop < $cuff_stop[$j] );
				print TMAP "$glean_id[$i]\t$cuff_id[$j]\n";
			}
		}
		$loci =~ s/,/\t/;
		$loci =~ s/START/$gene_start/;
		$loci =~ s/STOP/$gene_stop/;
		print LOCI $loci;
		
	}
	
	$loci = "";
	for my $i ( 0..$#cuff_id )
	{
		next if( $cuff_tag[$i] == 1 );
		$gene_start = $cuff_start[$i];
		$gene_stop = $cuff_stop[$i];
		$loci = "XLOCI_"."0" x (6-length($loci_num)).$loci_num."\t"."$scaff_id\[.\]START-STOP\t-\t$cuff_id[$i]\n";
		$loci_num++;
		$cuff_tag[$i] = 1;
		for my $j ( 0..$#cuff_id )
		{
			next if( $cuff_tag[$j] == 1 );
			if( (($cuff_start[$j]>=$cuff_start[$i])&&($cuff_start[$j]<=$cuff_stop[$i])) ||
			    (($cuff_stop[$j]>=$cuff_start[$i])&&($cuff_stop[$j]<=$cuff_stop[$i])) ||
			    (($cuff_start[$j]>=$cuff_start[$i])&&($cuff_stop[$j]<=$cuff_stop[$i])) ||
			    (($cuff_start[$j]<=$cuff_start[$i])&&($cuff_stop[$j]>=$cuff_stop[$i])) )
			{
				$loci =~ s/\n/,$cuff_id[$j]\n/;
				$cuff_tag[$j] = 1;
				$gene_start = $cuff_start[$j] if( $gene_start>$cuff_start[$j] );
				$gene_stop = $cuff_stop[$j] if( $gene_stop<$cuff_stop[$j] );
			}
		}
		$loci =~ s/START/$gene_start/;
		$loci =~ s/STOP/$gene_stop/;
		print LOCI $loci;
	}
}

for my $scaff_id ( keys %gtf )
{
	next if( exists $gff{$scaff_id} );
	
	my @cuff_id = @{$gtf{$scaff_id}{"id"}};
	my @cuff_start = @{$gtf{$scaff_id}{"start"}};
	my @cuff_stop = @{$gtf{$scaff_id}{"stop"}};
	my @cuff_tag;
	for my $i ( 0..$#cuff_id )
	{
		$cuff_tag[$i] = 0;
	}
	my ($loci, $gene_start, $gene_stop,);	
	
	$loci = "";
	for my $i ( 0..$#cuff_id )
	{
		next if( $cuff_tag[$i] == 1 );
		$gene_start = $cuff_start[$i];
		$gene_stop = $cuff_stop[$i];
		$loci = "XLOCI_"."0" x (6-length($loci_num)).$loci_num."\t"."$scaff_id\[.\]START-STOP\t-\t$cuff_id[$i]\n";
		$loci_num++;
		$cuff_tag[$i] = 1;
		for my $j ( 0..$#cuff_id )
		{
			next if( $cuff_tag[$j] == 1 );
			if( (($cuff_start[$j]>=$cuff_start[$i])&&($cuff_start[$j]<=$cuff_stop[$i])) ||
			    (($cuff_stop[$j]>=$cuff_start[$i])&&($cuff_stop[$j]<=$cuff_stop[$i])) ||
			    (($cuff_start[$j]>=$cuff_start[$i])&&($cuff_stop[$j]<=$cuff_stop[$i])) ||
			    (($cuff_start[$j]<=$cuff_start[$i])&&($cuff_stop[$j]>=$cuff_stop[$i])) )
			{
				$loci =~ s/\n/,$cuff_id[$j]\n/;
				$cuff_tag[$j] = 1;
				$gene_start = $cuff_start[$j] if( $gene_start>$cuff_start[$j] );
				$gene_stop = $cuff_stop[$j] if( $gene_stop<$cuff_stop[$j] );
			}
		}
		$loci =~ s/START/$gene_start/;
		$loci =~ s/STOP/$gene_stop/;
		print LOCI $loci;
	}
}

close LOCI;
close TMAP;

#=====================================
sub read_glean
{
	open(FH, shift) || die $!;
	my %gff = ();
	while( <FH> )
	{
		my @tmp = split /\t/;
	  if( $tmp[2] eq "mRNA" )
	  {
	  	my @info = split(";",$tmp[8]);
		  my ($id) = $info[0] =~ /ID=(\S+)/;
		  if( exists $gff{$tmp[0]} )
		  {
		  	push @{$gff{$tmp[0]}{"id"}}, $id;
		  	push @{$gff{$tmp[0]}{"start"}}, $tmp[3];
		  	push @{$gff{$tmp[0]}{"stop"}}, $tmp[4];
		  }
		  else
		  {
		    $gff{$tmp[0]}{"id"} = [ $id ];
		    $gff{$tmp[0]}{"start"} = [ $tmp[3] ];
		    $gff{$tmp[0]}{"stop"} = [ $tmp[4] ];
		  }
	  }
	}
	close FH;
	return %gff
}
#=====================================
sub read_cuff
{
	open(FH, shift) || die $!;
	my %gtf = ();
	while( <FH> )
	{
		my @tmp = split /\t/;
		if( $tmp[2] eq "transcript" )
		{
			my ($id) = $tmp[8] =~ /transcript_id "(\S+)"/;
			if( exists $gtf{$tmp[0]} )
		  {
		  	push(@{$gtf{$tmp[0]}{"id"}}, $id);
		  	push(@{$gtf{$tmp[0]}{"start"}}, $tmp[3]);
		  	push(@{$gtf{$tmp[0]}{"stop"}}, $tmp[4]);
		  }
		  else
		  {
		    $gtf{$tmp[0]}{"id"} = [ $id ];
		    $gtf{$tmp[0]}{"start"} = [ $tmp[3] ];
		    $gtf{$tmp[0]}{"stop"} = [ $tmp[4] ];
		  }
		}
	}
	close FH;
	return %gtf;
}
