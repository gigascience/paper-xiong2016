#!/usr/bin/perl
use strict;
use warnings;

# Author  liugeng
# Email   liugeng@genomics.cn
# create  2010-11-04
# mender  2011-03-21

die"perl $0 <GLEAN.GFF> <CUFF.GTF> <TRANSCRIPT> <LOCI> <MATCH> <PREFIX>\n" if(@ARGV != 6);
my($gff, $gtf, $transcript, $loci, $match, $prefix) = @ARGV;

my %loci = read_loci($loci);
my %match = read_match($match);
#my %part = read_match_part($match);
my %gff = read_glean_gff($gff);
my %gtf = read_cuff_gtf($gtf);
my %transcript = read_transcript($transcript);
my %glean_loci = ();
my %record = ();
my %cufflinks = ();
my $gene_num = 0;
for my $loci_id  (keys %loci){
        my $i;
	my $scaff_id = $loci{$loci_id}{"scaff"};
        my $strand = $loci{$loci_id}{"strand"};
        my $gene_start = $loci{$loci_id}{"start"};
        my $gene_stop = $loci{$loci_id}{"stop"};
	next if(!exists $loci{$loci_id}{"glean"});
=head
	if(! exists $loci{$loci_id}{"glean"}){
		print STDERR "$loci_id\n";
		next;
	}
=cut
	my $glean_id=$loci{$loci_id}{"glean"};
	next if(! $glean_id);
	if(!exists $gff{$glean_id}{CDS}){
		print STDERR "$glean_id not has CDS\n";
		next;
	}
=head
	if(! exists $gff{$glean_id}{mRNA} || ! exists $gff{$glean_id}{CDS}){
		print STDERR "$glean_id has no mRNA or CDS\n";
		next;
	}
=cut
	next if(! exists $gff{$glean_id});
	$glean_loci{$glean_id}++;
	my @cuff_id_list=@{$loci{$loci_id}{"cuff"}};
	my $flag = 0;
        for  $i ( 0..$#cuff_id_list ){
                if( exists $gtf{$cuff_id_list[$i]} )
                {
                        $flag = 1;
                        last;
                }
        }
	if($flag == 1){
		print "$scaff_id\t$prefix\tgene\t$gene_start\t$gene_stop\t.\t$strand\t.\t"."ID=$prefix"."0" x (6-length($gene_num)).$gene_num.";\n";
	}
	my $mRNA_num = 0;
	my $utr_5 = "";
	my $utr_3 = "";
	my @tmp;
#	if(!exists $gff{$glean_id}{CDS}){print STDERR "$glean_id\n";}
	my $glean_cds=@{$gff{$glean_id}{CDS}};
	my ($utr_5_s,$utr_5_e,$utr_3_s,$utr_3_e);
	if( exists $match{"identical"}{$glean_id} ){
		@tmp=@{$match{"identical"}{$glean_id}};
		foreach($i=0;$i<=$#tmp;$i++){
			my $cuff_id=$tmp[$i];
			$cufflinks{$cuff_id}=1;
			if( $gff{$glean_id}{mRNA}[3] eq "+"  && $gtf{$cuff_id}{transcript}[3] eq '+' ){
				if(! $utr_5_s){
					$utr_5_s = $transcript{$cuff_id}[0];
					$record{$glean_id}{utr_5_s}=$cuff_id;
				}else{
					($utr_5_s,$record{$glean_id}{utr_5_s}) = ($transcript{$cuff_id}[0]<$utr_5_s)?($transcript{$cuff_id}[0],$cuff_id):($utr_5_s,$record{$glean_id}{utr_5_s});
				}
				if(! $utr_3_e){
					$utr_3_e = $transcript{$cuff_id}[1];
					$record{$glean_id}{utr_3_e} = $cuff_id;
				}else{
					($utr_3_e,$record{$glean_id}{utr_3_e}) = ($transcript{$cuff_id}[1]>$utr_3_e)?($transcript{$cuff_id}[1],$cuff_id):($utr_3_e,$record{$glean_id}{utr_3_e});
				}
				$utr_5_e = $gff{$glean_id}{mRNA}[0] - 1;
				$utr_3_s = $gff{$glean_id}{mRNA}[1] + 1;
			}elsif($gff{$glean_id}{mRNA}[3] eq "-"  && $gtf{$cuff_id}{transcript}[3] eq '-'){
				if(! $utr_5_e){
					$utr_5_s = $transcript{$cuff_id}[1];
					$record{$glean_id}{utr_5_s}=$cuff_id;
				}else{
					($utr_5_s,$record{$glean_id}{utr_5_s}) = ($transcript{$cuff_id}[1]>$utr_5_e)?($transcript{$cuff_id}[1],$cuff_id):($utr_5_e,$record{$glean_id}{utr_5_s});
				}
				if(! $utr_3_e){
					$utr_3_e=$transcript{$cuff_id}[0];
					$record{$glean_id}{utr_3_e} = $cuff_id;
				}else{
					($utr_3_e,$record{$glean_id}{utr_3_e}) = ($transcript{$cuff_id}[0] < $utr_3_e)?($transcript{$cuff_id}[0],$cuff_id):($utr_3_e,$record{$glean_id}{utr_3_e});
				}
				$utr_5_e = $gff{$glean_id}{mRNA}[1] + 1;
				$utr_3_s = $gff{$glean_id}{mRNA}[0] - 1;
			}
		}
	}
	if( exists $match{"part"}{$glean_id} ){
		my @cuff_part = @{$match{"part"}{$glean_id}};
		for $i (0..$#cuff_part){
			my $cuff_id = $cuff_part[$i];
			$cufflinks{$cuff_id}=1;
			next if( not exists $gtf{$cuff_id});
			if($gff{$glean_id}{mRNA}[3] eq "+"  && $gtf{$cuff_id}{transcript}[3] eq '+'){
				if($gff{$glean_id}{mRNA}[0] == $gtf{$cuff_id}{transcript}[0]){
					if(! $utr_5_s){
						$utr_5_s = $transcript{$cuff_id}[0];
						$record{$glean_id}{utr_5_s}=$cuff_id;
					}else{
						($utr_5_s,$record{$glean_id}{utr_5_s}) = ($transcript{$cuff_id}[0]<$utr_5_s)?($transcript{$cuff_id}[0],$cuff_id):($utr_5_s,$record{$glean_id}{utr_5_s});
					}
				}
				if($gff{$glean_id}{mRNA}[1] == $gtf{$cuff_id}{transcript}[1]){
					if(! $utr_3_e){
						$utr_3_e = $transcript{$cuff_id}[1];
						$record{$glean_id}{utr_3_e} = $cuff_id;
					}else{
						($utr_3_e,$record{$glean_id}{utr_3_e}) = ($transcript{$cuff_id}[1]>$utr_3_e)?($transcript{$cuff_id}[1],$cuff_id):($utr_3_e,$record{$glean_id}{utr_3_e});
					}
				}
					$utr_5_e = $gff{$glean_id}{mRNA}[0] - 1;
					$utr_3_s = $gff{$glean_id}{mRNA}[1] + 1;
			}elsif($gff{$glean_id}{mRNA}[3] eq "-"  && $gtf{$cuff_id}{transcript}[3] eq '-'){
				if($gff{$glean_id}{mRNA}[1] == $gtf{$cuff_id}{transcript}[1]){
					if(! $utr_5_s){
						$utr_5_s = $transcript{$cuff_id}[1];
						$record{$glean_id}{utr_5_s}=$cuff_id;
					}else{
						($utr_5_s,$record{$glean_id}{utr_5_s}) = ($transcript{$cuff_id}[1]>$utr_5_e)?($transcript{$cuff_id}[1],$cuff_id):($utr_5_e,$record{$glean_id}{utr_5_s});
					}
				}
				if($gff{$glean_id}{mRNA}[0] == $gtf{$cuff_id}{transcript}[0]){
					if(! $utr_3_e){
						$utr_3_e=$transcript{$cuff_id}[0];
						$record{$glean_id}{utr_3_e} = $cuff_id;
					}else{
						($utr_3_e,$record{$glean_id}{utr_3_e}) = ($transcript{$cuff_id}[0] < $utr_3_e)?($transcript{$cuff_id}[0],$cuff_id):($utr_3_e,$record{$glean_id}{utr_3_e});
					}
				}
				$utr_5_e = $gff{$glean_id}{mRNA}[1] + 1;
				$utr_3_s = $gff{$glean_id}{mRNA}[0] - 1;
			}
		}
	}
	if($gff{$glean_id}{mRNA}[3] eq '+'){
		if($utr_5_s){
			$gff{$glean_id}{mRNA}[0]=$utr_5_s;
			$utr_5 = "$scaff_id\tGLEAN\tUTR_5\t$utr_5_s\t$utr_5_e\t.\t.\t.\t".
				"Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;support_id=$record{$glean_id}{utr_5_s};";
		}
		if($utr_3_e){
			$gff{$glean_id}{mRNA}[1]=$utr_3_e;
			$utr_3 = "$scaff_id\tGLEAN\tUTR_3\t$utr_3_s\t$utr_3_e\t.\t.\t.\t".
				"Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;support_id=$record{$glean_id}{utr_3_e};";
		}
	}elsif($gff{$glean_id}{mRNA}[3] eq '-'){
		if($utr_5_s){
			$gff{$glean_id}{mRNA}[1]=$utr_5_s;
			$utr_5 = "$scaff_id\tGLEAN\tUTR_5\t$utr_5_e\t$utr_5_s\t.\t.\t.\t".
				"Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;support_id=$record{$glean_id}{utr_5_s};";
		}
		if($utr_3_e){
			$gff{$glean_id}{mRNA}[0]=$utr_3_e;
			$utr_3 = "$scaff_id\tGLEAN\tUTR_3\t$utr_3_e\t$utr_3_s\t.\t.\t.\t".
			"Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;support_id=$record{$glean_id}{utr_3_e};";
		}
	}
	my $support_info = "";
	if(exists $match{"identical"}{$glean_id}){
		@tmp=@{$match{"identical"}{$glean_id}};
		for  $i (0..$#tmp){
			$support_info .= "identical_support_id=$tmp[$i];";
		}
	}
	if(exists $match{"part"}{$glean_id}){
		@tmp=@{$match{"part"}{$glean_id}};
		for $i (0..$#tmp){
			if( exists $gtf{$tmp[$i]} ){
				$support_info .= "part_support_id=$tmp[$i]," if ($support_info eq "");
				$support_info .= $tmp[$i]."," if( $support_info =~ /part_support_id=/ );
			}
		}
	}
	if( $support_info ne "" ){
		$support_info = substr($support_info, 0, length($support_info) - 1) . ";";
	}
	if( $support_info =~ /identical/ ){
		if( $support_info =~ /part/ ){
			next;
		}
		else{
			$support_info .= "part_support_id=.;"
		}
	}else{
		if( $support_info =~ /part/ ){
			$support_info = "identical_support_id=.;" . $support_info
		}
		else{
			$support_info .= "identical_support_id=.;part_support_id=.;";
		}
	}
	my @lines;
	if($gff{$glean_id}{mRNA}[3] eq '+'){
		print "$scaff_id\tGLEAN\tmRNA\t$gff{$glean_id}{mRNA}[0]\t$gff{$glean_id}{mRNA}[1]\t"."$gff{$glean_id}{mRNA}[2]\t$gff{$glean_id}{mRNA}[3]\t$gff{$glean_id}{mRNA}[4]\t"."ID=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;"."source_id=$glean_id;$support_info"."\n";
		@lines=sort{$a->[0] <=> $b->[0]} @{$gff{$glean_id}{CDS}};
	}
	elsif($gff{$glean_id}{mRNA}[3] eq '-'){
		print "$scaff_id\tGLEAN\tmRNA\t$gff{$glean_id}{mRNA}[0]\t$gff{$glean_id}{mRNA}[1]\t".
			"$gff{$glean_id}{mRNA}[2]\t$gff{$glean_id}{mRNA}[3]\t$gff{$glean_id}{mRNA}[4]\t".
			"ID=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;".
			"source_id=$glean_id;$support_info"."\n";
		@lines=sort{$b->[0] <=> $a->[0]} @{$gff{$glean_id}{CDS}};
	}
	print $utr_5."\n" if($utr_5 ne "");
	for $i(0..$#lines){
		print "$scaff_id\tGLEAN\tCDS\t$lines[$i][0]\t$lines[$i][1]\t".
		"$lines[$i][2]\t$lines[$i][3]\t$lines[$i][4]\t".
		"Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;\n";
	}
	print $utr_3."\n" if($utr_3 ne "");
	
	$mRNA_num++;
	for(my $k=0; $k<=$#cuff_id_list; $k++){
		next if(exists $cufflinks{$cuff_id_list[$k]});
		next if( not exists $gtf{$cuff_id_list[$k]} );
		if($gtf{$cuff_id_list[$k]}{transcript}[3] eq '+'){
			print "$scaff_id\tCufflinks\tmRNA\t$gtf{$cuff_id_list[$k]}{transcript}[0]\t$gtf{$cuff_id_list[$k]}{transcript}[1]\t".
				"$gtf{$cuff_id_list[$k]}{transcript}[2]\t$gtf{$cuff_id_list[$k]}{transcript}[3]\t$gtf{$cuff_id_list[$k]}{transcript}[4]\t".
				"ID=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;".
				"source_id=$cuff_id_list[$k];\n";
			@{$gtf{$cuff_id_list[$k]}{exon}}=sort{$a->[0] <=> $b->[0]} @{$gtf{$cuff_id_list[$k]}{exon}};
			for  $i(0.. $#{$gtf{$cuff_id_list[$k]}{exon}}){
				print "$scaff_id\tCufflinks\tCDS\t$gtf{$cuff_id_list[$k]}{exon}[$i][0]\t$gtf{$cuff_id_list[$k]}{exon}[$i][1]\t".
					"$gtf{$cuff_id_list[$k]}{exon}[$i][2]\t$gtf{$cuff_id_list[$k]}{exon}[$i][3]\t$gtf{$cuff_id_list[$k]}{exon}[$i][4]\t".
					"Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;\n";
			}
		}
		elsif($gtf{$cuff_id_list[$k]}{transcript}[3] eq '-'){
			print "$scaff_id\tCufflinks\tmRNA\t$gtf{$cuff_id_list[$k]}{transcript}[0]\t$gtf{$cuff_id_list[$k]}{transcript}[1]\t".
				"$gtf{$cuff_id_list[$k]}{transcript}[2]\t$gtf{$cuff_id_list[$k]}{transcript}[3]\t$gtf{$cuff_id_list[$k]}{transcript}[4]\t".
				"ID=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;".
				"source_id=$cuff_id_list[$k];\n";
			@{$gtf{$cuff_id_list[$k]}{exon}}=sort{$b->[0] <=> $a->[0]} @{$gtf{$cuff_id_list[$k]}{exon}};
			for $i(0.. $#{$gtf{$cuff_id_list[$k]}{exon}}){
				print "$scaff_id\tCufflinks\tCDS\t$gtf{$cuff_id_list[$k]}{exon}[$i][0]\t$gtf{$cuff_id_list[$k]}{exon}[$i][1]\t".
					"$gtf{$cuff_id_list[$k]}{exon}[$i][2]\t$gtf{$cuff_id_list[$k]}{exon}[$i][3]\t$gtf{$cuff_id_list[$k]}{exon}[$i][4]\t".
					"Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;\n";
			}
		}
		$mRNA_num++;
		$cufflinks{$cuff_id_list[$k]}=1;
	}
	if($flag == 1){
		$gene_num++;
	}
}
for my $glean_id_remain ( keys %gff ){
	next if( exists $glean_loci{$glean_id_remain} );

	if(! exists $gff{$glean_id_remain}{mRNA} || !exists $gff{$glean_id_remain}{CDS}){
		print STDERR "$glean_id_remain not has mRNA or CDS\n";
		next;
	}
	my $scaff_id = $gff{$glean_id_remain}{mRNA}[5];
	my $mRNA_num = 0;
	print "$scaff_id\tGLEAN\tgene\t$gff{$glean_id_remain}{mRNA}[0]\t$gff{$glean_id_remain}{mRNA}[1]\t".
		"ID=$prefix"."0" x (6-length($gene_num)).$gene_num.";\n";
	print "$scaff_id\tGLEAN\tmRNA\t$gff{$glean_id_remain}{mRNA}[0]\t$gff{$glean_id_remain}{mRNA}[1]\t".
		"$gff{$glean_id_remain}{mRNA}[2]\t$gff{$glean_id_remain}{mRNA}[3]\t$gff{$glean_id_remain}{mRNA}[4]\t".
		"ID=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;".
		"source_id=$glean_id_remain;identical_support_id=.;part_support_id=.;\n";
		my @lines;
		if($gff{$glean_id_remain}{mRNA}[3] eq '+'){
			@lines=sort{$a->[0] <=> $b->[0]} @{$gff{$glean_id_remain}{CDS}};
		}elsif($gff{$glean_id_remain}{mRNA}[3] eq '-'){
			@lines=sort{$b->[0] <=> $a->[0]} @{$gff{$glean_id_remain}{CDS}};
		}
		for my $i(0..$#lines){
			print "$scaff_id\tGLEAN\tCDS\t$lines[$i][0]\t$lines[$i][1]\t".
				"$lines[$i][2]\t$lines[$i][3]\t$lines[$i][4]\t".
				"Parent=$prefix"."0" x (6-length($gene_num)).$gene_num.".$mRNA_num;\n";
		}
	$gene_num++;	
}

#--------------------------------------
sub read_loci
{
        open (FH, shift) || die $!;
        my %loci = ();
#        my $loci_id = 0;
        while (<FH>)
        {
                my @tmp = split (/\t/,$_);
                next unless  $tmp[3];
		next if($tmp[2] eq '-');
#		$loci_id++;
		my $loci_id=$tmp[0];
                my ($scaff_id, $strand, $start, $stop) = $tmp[1] =~ /(\S+)\[(\S+)](\d+)-(\d+)/;
                $loci{$loci_id}{"scaff"} = $scaff_id;
                $loci{$loci_id}{"strand"} = $strand;
                $loci{$loci_id}{"start"} = $start;
                $loci{$loci_id}{"stop"} = $stop;

                for my $ref_id ( $tmp[2] =~ /[^,]+\|([^,]+)/g )
                {
                        $loci{$loci_id}{"glean"}=$ref_id;
                }
                for my $cuff_id(split(",", $tmp[3]))
                {
                        push @{$loci{$loci_id}{"cuff"}},$cuff_id ;
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
                        push @{$match{"identical"}{$tmp[1]}}, $tmp[3];
                        $match{"identical"}{$tmp[3]} = $tmp[1];
                }
                if( /#PART/ )
                {
                        my @tmp = split("\t");
                        push @{$match{"part"}{$tmp[1]}},$tmp[3];
                        $match{"part"}{$tmp[3]} = $tmp[1];
                }
        }
        close FH;
        return %match;
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
                my @tmp = split("\t",$_);
                my ($id) = $tmp[8] =~ /transcript_id "(\S+)"/;
                if( $tmp[2] eq "transcript" )
                {
                        $num = 0;
                        $gtf{$id}{transcript} =($tmp[3]<$tmp[4])?[ $tmp[3], $tmp[4], $tmp[5], $tmp[6], $tmp[7],$tmp[0] ]:[ $tmp[4], $tmp[3], $tmp[5], $tmp[6], $tmp[7] ,$tmp[0]];
                }
                if( $tmp[2] eq "exon" )
                {
			if($tmp[3]<$tmp[4]){
				push @{$gtf{$id}{exon}} ,[ $tmp[3], $tmp[4], $tmp[5], $tmp[6], $tmp[7] ,$tmp[0]];
			}else{
				push @{$gtf{$id}{exon}} ,[ $tmp[4], $tmp[3], $tmp[5], $tmp[6], $tmp[7] ,$tmp[0]];
			}
                }
        }
        close FH;
        return %gtf;
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
		s/^\s+|\s+$//g;
                my @tmp = split(/\t/,$_);
                if( $tmp[2] eq "mRNA" ){
			next if(@tmp <9);
			my $id=$1 if($tmp[8]=~/ID=([^;]+);/);
			$gff{$id}{mRNA} = ($tmp[3]<$tmp[4])?[ $tmp[3], $tmp[4], $tmp[5], $tmp[6], $tmp[7] ,$tmp[0]]:[ $tmp[4], $tmp[3], $tmp[5], $tmp[6], $tmp[7] ,$tmp[0]];
		}elsif( $tmp[2] eq "CDS" ){
			my $id=$1 if($tmp[8]=~/Parent=([^;]+);/);
			next if(! exists $gff{$id}{mRNA});
			if($tmp[3] < $tmp[4]){
                        	push @{$gff{$id}{CDS}} ,[ $tmp[3], $tmp[4], $tmp[5], $tmp[6], $tmp[7], $tmp[0] ,$tmp[0]];
			}else{
				push @{$gff{$id}{CDS}} ,[ $tmp[4], $tmp[3], $tmp[5], $tmp[6], $tmp[7], $tmp[0] ,$tmp[0]];
			}
                }
        }
        close FH;
        return %gff;
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
                $transcript{$id} =($tmp[3]<$tmp[4])?[$tmp[3],$tmp[4]]:[$tmp[4],$tmp[3]];
                }
        }
        close FH;
        return %transcript;
}
#==================================
