#!/usr/bin/perl
use strict;

# Author  liugeng   
# Email   liugeng@genomics.cn
# create  2010-11-04

die"perl $0 <GLEAN.GFF> <CUFF.GTF> <TRANSCRIPT> <LOCI> <MATCH> <PREFIX>\n" if(@ARGV != 6);
my($gff, $gtf, $transcript, $loci, $match, $prefix) = @ARGV;

my %loci = read_loci($loci);
my %match = read_match($match);
my %part = read_match_part($match);
my %gff = read_glean_gff($gff);
my %gtf = read_cuff_gtf($gtf);
my %transcript = read_transcript($transcript);
my %glean_loci = ();

my $gene_num = 0;
for my $loci_id (sort {$a<=>$b} (keys %loci) )
{
	my $scaff_id = $loci{$loci_id}{"scaff"};
	my $strand = $loci{$loci_id}{"strand"};
	my $gene_start = $loci{$loci_id}{"start"};
	my $gene_stop = $loci{$loci_id}{"stop"};
	my @glean_id_list = (keys %{$loci{$loci_id}{"glean"}});
	my @cuff_id_list = (keys %{$loci{$loci_id}{"cuff"}});
  
	my $flag = 0;
	for my $i ( 0..$#glean_id_list )
	{
		if( exists $gff{$glean_id_list[$i]} )
		{
			$flag = 1;
			last;
		}
	}  
	for my $i ( 0..$#cuff_id_list )
	{
		if( exists $gtf{$cuff_id_list[$i]} )
		{
			$flag = 1;
			last;
		}
	}  
  
	if($flag == 1)
	{
		print "$scaff_id\t$prefix\tgene\t$gene_start\t$gene_stop\t.\t$strand\t.\t".
			  "ID=$prefix"."0" x (6-length($gene_num)).$gene_num.";\n";
	}
  
	my $mRNA_num = 0;
	for(my $k=0; $k<=$#glean_id_list; $k++)
	{
		$glean_loci{$glean_id_list[$k]} = 1;
  	
		my $utr_5 = "";
		my $utr_3 = "";
		if( exists $match{"identical"}{$glean_id_list[$k]} )
		{
			my $cuff_id = $match{"identical"}{$glean_id_list[$k]};
			my $glean_cds = ( keys %{$gff{$glean_id_list[$k]}} );
			my ($utr_5_s,$utr_5_e,$utr_3_s,$utr_3_e);
			if( $gff{$glean_id_list[$k]}{0}[3] eq "+" )
			{
				$utr_5_s = $transcript{$cuff_id}[0];
				$utr_5_e = $gff{$glean_id_list[$k]}{1}[0] - 1;
  		
				$utr_3_s = $gff{$glean_id_list[$k]}{$glean_cds - 1}[1] + 1;
				$utr_3_e = $transcript{$cuff_id}[1];
			}
			else
			{
				$utr_3_s = $transcript{$cuff_id}[0];
				$utr_3_e = $gff{$glean_id_list[$k]}{$glean_cds - 1}[0] - 1;
  			
				$utr_5_s = $gff{$glean_id_list[$k]}{1}[1] + 1;
				$utr_5_e = $transcript{$cuff_id}[1];
			}
			$utr_5 = "$scaff_id\tGLEAN\tUTR_5\t$utr_5_s\t$utr_5_e\t.\t.\t.\t".
					 "Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;support_id=$cuff_id;" if $utr_5_e>$utr_5_s;
			$utr_3 = "$scaff_id\tGLEAN\tUTR_3\t$utr_3_s\t$utr_3_e\t.\t.\t.\t".
				     "Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;support_id=$cuff_id;" if $utr_3_e>$utr_3_s;
		}
  	
		if( exists $part{$glean_id_list[$k]} )
		{
			my $glean_cds = ( keys %{$gff{$glean_id_list[$k]}} );
			my ($utr_5_s,$utr_5_e,$utr_3_s,$utr_3_e);
			$utr_5_s = $gff{$glean_id_list[$k]}{0}[1];
			$utr_3_e = 0;
			my @cuff_part = @{$part{$glean_id_list[$k]}};
			for my $i ( 0..$#cuff_part )
			{
				my $cuff_id = $cuff_part[$i];
				next if( not exists $gtf{$cuff_id});
				if( ($gff{$glean_id_list[$k]}{0}[0] == $gtf{$cuff_id}{0}[0]) && ($utr_5_s > $transcript{$cuff_id}[0]) )
				{
					if( $gff{$glean_id_list[$k]}{0}[3] eq "+" )
					{
						$utr_5_s = $transcript{$cuff_id}[0];
						$utr_5_e = $gff{$glean_id_list[$k]}{1}[0] - 1;
						$utr_5 = "$scaff_id\tGLEAN\tUTR_5\t$utr_5_s\t$utr_5_e\t.\t.\t.\t".
						         "Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;support_id=$cuff_id;" if $utr_5_e>$utr_5_s;
					}
					else
					{
						$utr_3_s = $transcript{$cuff_id}[0];
						$utr_3_e = $gff{$glean_id_list[$k]}{$glean_cds - 1}[0] - 1;
						$utr_3 = "$scaff_id\tGLEAN\tUTR_3\t$utr_3_s\t$utr_3_e\t.\t.\t.\t".
						         "Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;support_id=$cuff_id;" if $utr_3_e>$utr_3_s;
					}
					
				}
				if( ($gff{$glean_id_list[$k]}{0}[1] == $gtf{$cuff_id}{0}[1]) && ($utr_3_e < $transcript{$cuff_id}[1]) )
				{
					if( $gff{$glean_id_list[$k]}{0}[3] eq "+" )
					{
						$utr_3_s = $gff{$glean_id_list[$k]}{$glean_cds - 1}[1] + 1;
						$utr_3_e = $transcript{$cuff_id}[1];
						$utr_3 = "$scaff_id\tGLEAN\tUTR_3\t$utr_3_s\t$utr_3_e\t.\t.\t.\t".
						         "Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;support_id=$cuff_id;" if $utr_3_e>$utr_3_s;
					}
					else
					{
						$utr_5_s = $gff{$glean_id_list[$k]}{1}[1] + 1;
						$utr_5_e = $transcript{$cuff_id}[1];
						$utr_5 = "$scaff_id\tGLEAN\tUTR_5\t$utr_5_s\t$utr_5_e\t.\t.\t.\t".
						         "Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;support_id=$cuff_id;" if $utr_5_e>$utr_5_s;
					}
					
				}
			}
		}
  	
  	my $support_info = "";
  	if( exists $match{"identical"}{$glean_id_list[$k]} )
 	  {
  		my $cuff_id = $match{"identical"}{$glean_id_list[$k]};
  		$support_info .= "identical_support_id=$cuff_id;";
  	}
  	if( exists $match{"part"}{$glean_id_list[$k]} )
  	{
  		for my $cuff_id ( $match{"part"}{$glean_id_list[$k]} )
  		{
  			if( exists $gtf{$cuff_id} )
  			{
  				$support_info .= "part_support_id=$cuff_id," if ($support_info eq "" );
  				$support_info .= $cuff_id."," if( $support_info =~ /part_support_id=/ );
  			}
  		}
  		if( $support_info ne "" )
  		{
  		  $support_info = substr($support_info, 0, length($support_info) - 1) . ";";
  		}
  	}
  	if( $support_info =~ /identical/ )
  	{
  		if( $support_info =~ /part/ )
  		{
  			next;
  		}
  		else
  		{
  			$support_info .= "part_support_id=.;"
  		}
  	}
  	else
  	{
  		if( $support_info =~ /part/ )
  		{
  			$support_info = "identical_support_id=.;" . $support_info
  		}
  		else
  		{
  			$support_info .= "identical_support_id=.;part_support_id=.;"
  		}
  	}
  	
  	print "$scaff_id\tGLEAN\tmRNA\t$gff{$glean_id_list[$k]}{0}[0]\t$gff{$glean_id_list[$k]}{0}[1]\t".
  			  "$gff{$glean_id_list[$k]}{0}[2]\t$gff{$glean_id_list[$k]}{0}[3]\t$gff{$glean_id_list[$k]}{0}[4]\t".
  			  "ID=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;".
  			  "source_id=$glean_id_list[$k];$support_info\n";
  	print $utr_5."\n" if($utr_5 ne "");
  	my @array = keys %{$gff{$glean_id_list[$k]}};
  	for my $i ( 1..$#array )
  	{
  		print "$scaff_id\tGLEAN\tCDS\t$gff{$glean_id_list[$k]}{$i}[0]\t$gff{$glean_id_list[$k]}{$i}[1]\t".
  			    "$gff{$glean_id_list[$k]}{$i}[2]\t$gff{$glean_id_list[$k]}{$i}[3]\t$gff{$glean_id_list[$k]}{$i}[4]\t".
  			    "Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;\n";
  	}
  	print $utr_3."\n" if($utr_3 ne "");
  	$mRNA_num++;
  }
  
  for(my $k=0; $k<=$#cuff_id_list; $k++)
  {
  	next if( exists $match{"identical"}{$cuff_id_list[$k]} );
  	next if( exists $match{"part"}{$cuff_id_list[$k]} );
  	next if( not exists $gtf{$cuff_id_list[$k]} );
  	print "$scaff_id\tCufflinks\tmRNA\t$gtf{$cuff_id_list[$k]}{0}[0]\t$gtf{$cuff_id_list[$k]}{0}[1]\t".
  			  "$gtf{$cuff_id_list[$k]}{0}[2]\t$gtf{$cuff_id_list[$k]}{0}[3]\t$gtf{$cuff_id_list[$k]}{0}[4]\t".
  			  "ID=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;".
  			  "source_id=$cuff_id_list[$k];\n";
  	my @array = keys %{$gtf{$cuff_id_list[$k]}};
  	for my $i ( 1..$#array )
  	{
  		print "$scaff_id\tCufflinks\tCDS\t$gtf{$cuff_id_list[$k]}{$i}[0]\t$gtf{$cuff_id_list[$k]}{$i}[1]\t".
  			    "$gtf{$cuff_id_list[$k]}{$i}[2]\t$gtf{$cuff_id_list[$k]}{$i}[3]\t$gtf{$cuff_id_list[$k]}{$i}[4]\t".
  			    "Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;\n";
  	}
  	$mRNA_num++;
  }
  if($flag == 1)
  {
    $gene_num++;
  }
}

#print "#GLEAN Not In Loci\n";
for my $glean_id ( keys %gff )
{
	next if( exists $glean_loci{$glean_id} );
	my $scaff_id = $gff{$glean_id}{0}[5];
	my $mRNA_num = 0;
	for my $i (sort {$a<=>$b} (keys %{$gff{$glean_id}}))
  {
  	if( $i == 0 )
  	{
  		print "$scaff_id\tGLEAN\tmRNA\t$gff{$glean_id}{0}[0]\t$gff{$glean_id}{0}[1]\t".
  			    "$gff{$glean_id}{0}[2]\t$gff{$glean_id}{0}[3]\t$gff{$glean_id}{0}[4]\t".
  			    "ID=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;".
  			    "source_id=$glean_id;identical_support_id=.;part_support_id=.;\n";
  	}
  	else
  	{
  		print "$scaff_id\tGLEAN\tCDS\t$gff{$glean_id}{$i}[0]\t$gff{$glean_id}{$i}[1]\t".
  			    "$gff{$glean_id}{$i}[2]\t$gff{$glean_id}{$i}[3]\t$gff{$glean_id}{$i}[4]\t".
  			    "Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;\n";
  	}
  }
  $gene_num++;
}

#============================
sub read_loci
{
	open (FH, shift) || die $!;
	my %loci = ();
	my $loci_id = 0;
	while (<FH>)
	{
		my @tmp = split /\s+/;
		$loci_id++;
		
		my ($scaff_id, $strand, $start, $stop) = $tmp[1] =~ /(\S+)\[(\S+)](\d+)-(\d+)/;
		$loci{$loci_id}{"scaff"} = $scaff_id;
		$loci{$loci_id}{"strand"} = $strand;
		$loci{$loci_id}{"start"} = $start;
		$loci{$loci_id}{"stop"} = $stop;
		
		for my $ref_id ( $tmp[2] =~ /[^,]+\|([^,]+)/g )
		{
			$loci{$loci_id}{"glean"}{$ref_id} = 1;
		}
		next unless $tmp[3];
		for my $cuff_id(split(",", $tmp[3]))
		{
			$loci{$loci_id}{"cuff"}{$cuff_id} = 1;
		}
	}
	close FH;
	return %loci;
}
#=============================
sub read_match
{
	open (FH, shift) || die $!;
	my %match = ();
	while (<FH>)
	{
		chomp;
		if( /#IDENTICAL/ )
		{
			my @tmp = split("\t");
			$match{"identical"}{$tmp[1]} = $tmp[3];
			$match{"identical"}{$tmp[3]} = $tmp[1];
		}
		if( /#PART/ )
		{
			my @tmp = split("\t");
			$match{"part"}{$tmp[1]} = $tmp[3];
			$match{"part"}{$tmp[3]} = $tmp[1];
		}
	}
	close FH;
	return %match;			
}
#==================================
sub read_match_part
{
	open (FH, shift) || die $!;
	my %part = ();
	while(<FH>)
	{
		chomp;
		if( /#PART/ )
		{
			my @tmp = split("\t");
			if( exists $part{$tmp[1]} )
			{
				push @{$part{$tmp[1]}}, $tmp[3];
			}
			else
			{
				$part{$tmp[1]} = [$tmp[3]];
			}
		}
	}
	close FH;
	return %part;
}
#==================================
sub read_glean_gff
{
	open (FH, shift) || die $!;
	my %gff = ();
	my $num;
	while (<FH>)
	{
		chomp;
		next if (/^#/ || /^\s/);
		my @tmp = split("\t");
		if( $tmp[2] eq "mRNA" )
		{
			$num = 0;
			my @info = split(";",$tmp[8]);
			my ($id) = $info[0] =~ /ID=(\S+)/;
			$gff{$id}{$num++} = [ $tmp[3], $tmp[4], $tmp[5], $tmp[6], $tmp[7], $tmp[0] ];
		}
		if( $tmp[2] eq "CDS" )
		{
			my @info = split(";",$tmp[8]);
			my ($id) = $info[0] =~ /Parent=(\S+)/;
			$gff{$id}{$num++} = [ $tmp[3], $tmp[4], $tmp[5], $tmp[6], $tmp[7], $tmp[0] ];
		}
	}
	close FH;
	return %gff;
}
#==================================
sub read_cuff_gtf
{
	open (FH, shift) || die $!;
	my %gtf = ();
	my $num;
	while (<FH>)
	{
		chomp;
		next if (/^#/ || /^\s/);
		my @tmp = split("\t");
		my ($id) = $tmp[8] =~ /transcript_id "(\S+)"/;
		if( $tmp[2] eq "transcript" )
		{
			$num = 0;
			$gtf{$id}{$num++} = [ $tmp[3], $tmp[4], $tmp[5], $tmp[6], $tmp[7] ];
		}
		if( $tmp[2] eq "exon" )
		{
			$gtf{$id}{$num++} = [ $tmp[3], $tmp[4], $tmp[5], $tmp[6], $tmp[7] ];
		}
	}
	close FH;
	return %gtf;
}
#==================================
sub read_transcript
{
	open(FH, shift) || die $!;
	my %transcript = ();
	while(<FH>)
	{
		chomp;
		my @tmp = split("\t");
		if( $tmp[2] eq "transcript" )
		{
	  	my ($id) = $tmp[8] =~ /transcript_id "(\S+)"/;
	  	$transcript{$id} = [ $tmp[3], $tmp[4] ];
		}
	}
	close FH;
	return %transcript;
}
#==================================
