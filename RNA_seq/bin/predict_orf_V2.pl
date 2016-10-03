#!/usr/bin/perl 
use strict;

# Author: liugeng
# E-mail: liugeng@genomics.org.cn
# Create: 2010-08-23
# Updata: 2010-11-01

# This program is to predict the phase of exon which was predictd by Cufflinks

die "Usage:$0 <GENOME.FA> <CUFFLINKS.GTF> <PARAMETER> <EXON_LEN_CUTOFF> <PREFIX>\n" if (@ARGV !=5);
my($genome_fasta, $cufflinks_gtf, $parameter, $cutoff, $prefix)=@ARGV;

my $initial_bases=1024; # number of markov initial probability
my $transition_bases=4096; # number of markov transition probability

# Read markov parameter into markov_parameter
my %markov=read_markov_parameter($parameter);

open (CUFFLINKS_GTF,"<$cufflinks_gtf") ||die("can not open cufflinks file!\n");
my %gtf = ();
while (<CUFFLINKS_GTF>){
	next if (/^#/ || /^\s/);
	my @case = split("\t");	
 	
	if($case[2] eq 'exon'){
		my ($id)  = $case[8] =~ /transcript_id "(\S+)"/;

		$gtf{$case[0]}{$id}{$case[3]} = $case[4];
	}
}
close(CUFFLINKS_GTF);

open (GENOME_FASTA,"<$genome_fasta") ||die("can not open genome file!\n");
my %exons=();
$/="\n>";
while(<GENOME_FASTA>)
{
	chomp;
  s/^>//;
  my ($tag, $seq) = split(/\n/, $_, 2);
  $seq =~ s/\W+//g;
  $seq =~ s/\n//g;
  my ($scaff_id) = $tag =~ /^(\S+)/;
  next unless exists $gtf{$scaff_id};
  for my $transcript_id (keys %{$gtf{$scaff_id}})
  {
		$exons{$transcript_id} = '';
		for my $start (sort{$a<=>$b} keys %{$gtf{$scaff_id}{$transcript_id}})
		{
			$exons{$transcript_id} .= substr(	$seq, $start - 1, $gtf{$scaff_id}{$transcript_id}{$start} - $start + 1);
		}
	}
}
$/="\n";
close(GENOME_FASTA);
%gtf=();

open (CUFFLINKS_GTF,"<$cufflinks_gtf") ||die("can not open cufflinks file!\n");
open (CDS, ">$prefix\.cds") || die;
open (ORF, ">$prefix\.orf") ||die;
open (CUFFLINKS_GFF,">$prefix\.orf.gtf") ||die("can not create cufflinks file!\n");
$_ = <CUFFLINKS_GTF>;
while (defined)
{
	chomp;
	# Processing a transcript
	my @cufflinks = ();
	my $num = 0;
	if(/Cufflinks	transcript|TAIR9	mRNA/)
	{
		$cufflinks[$num++] = [split(/\t/, $_)];
		my ($transcript_id)  = $_ =~ /transcript_id "(\S+)"/;
		while (<CUFFLINKS_GTF>)
		{
			chomp;
			if(/Cufflinks	exon/)
			{
				$cufflinks[$num++] = [split(/\t/, $_)];
			}
			else{last;}
		}
		next unless exists $exons{$transcript_id};
		my $exons = $exons{$transcript_id};
		$exons =~ tr/acgt/ACGT/;
		if( $cufflinks[0][6] eq '-')
		{
			$exons = reverse($exons);
		  $exons =~ tr/ACGT/TGCA/;
		}
		
		print ">$transcript_id\n";
		
		my $strand = "+";
		my ($score_tmp, $frame_tmp, $start_tmp, $stop_tmp);
		STRAND:
		my @cds = split("", $exons);
		my @start;
		my @stop;
		my @score;
		for(my $frame=0; $frame<3; ++$frame)
		{
			# Look for stop
			$stop[$frame] = $#cds - (($#cds+1-$frame)%3);
			for(my $i=0+$frame; $i<=$#cds-$frame-3; $i+=3) 
			{
				my $codon = $cds[$i].$cds[$i+1].$cds[$i+2];
				if( (($codon eq "TAA") || ($codon eq "TAG") || ($codon eq "TGA")) && ($i >= int($#cds/2)) )
				{
					$stop[$frame] = $i + 2;
					last;
				}
			}
			# Look for start
			my $utr_stop = $frame;
			for(my $i=$stop[$frame] - 3; $i>=$frame+3; $i-=3)
			{
				my $codon = $cds[$i-2].$cds[$i-1].$cds[$i];
				if( ($codon eq "TAA") || ($codon eq "TAG") || ($codon eq "TGA") )
				{
				  $utr_stop = $i + 1;
				  last;
				}
			}                
		  $start[$frame] = $utr_stop;
			for(my $i=$utr_stop; $i<=int($#cds/2); $i+=3)
			{
				my $codon = $cds[$i].$cds[$i+1].$cds[$i+2];
				if( $codon eq "ATG" )
				{
					$start[$frame] = $i;
					last;
				}
			}
			my @orf = ();
			for(my $i=$start[$frame],my $k=0; $i<=$stop[$frame]; ++$i,++$k)
			{
				$orf[$k] = $cds[$i];
			}
			$score[$frame] = score_markov_probability(@orf);
			#print "$start[$frame]\t$stop[$frame]\t$score[$frame]\n";
			
#			if($utr_stop == $frame)
#			{
#			  @orf = ();
#			  $start[3+$frame]=$frame;
#			  $stop[3+$frame]= $stop[$frame];
#			  for(my $i=$frame,my $k=0; $i<=$stop[$frame]; ++$i,++$k)
#			  {
#				  $orf[$k] = $cds[$i];
#			  }
#			  $score[3+$frame] = score_markov_probability(@orf);
#		  }
#		  else
#		  {
#		  	$start[3+$frame]=0;
#			  $stop[3+$frame]= 0;
#			  $score[3+$frame]=-9999
#		  }
#			#print "$start[3+$frame]\t$stop[3+$frame]\t$score[3+$frame]\n";
		}
		my ($score, $frame) = calculate_phrase(@score);
		my $start = $start[$frame];
		my $stop = $stop[$frame];
		my $frame = $frame % 3;
		
		if($cufflinks[0][6] eq ".")
		{
			if( $strand eq "+")
			{
				 ($score_tmp, $frame_tmp, $start_tmp, $stop_tmp) = ($score, $frame, $start, $stop);
			   $exons = reverse($exons);
		     $exons =~ tr/ACGT/TGCA/;
		     $strand = "-";
			   goto STRAND;
			}
			else
			{
				if( $score > $score_tmp )
				{
					$strand = "-";
				}
				else
				{
					$strand = "+";
					($score, $frame, $start, $stop) = ($score_tmp, $frame_tmp, $start_tmp, $stop_tmp);
					$exons = reverse($exons);
		      $exons =~ tr/ACGT/TGCA/;
				}
				for(my $i=0; $i<$num; ++$i)
				{
					$cufflinks[$i][6] = $strand;
				}
			}
		}
		
		print "$start\t$stop\t$cufflinks[0][6]\t$score\n";
		
		my $info_flag = 0;
		my $start_pos;
		my $stop_pos;
		if($cufflinks[0][6] eq "+")
		{
			$start_pos = 1;
			$stop_pos = $num - 1;
			my @len;
			$len[1] = $cufflinks[1][4]-$cufflinks[1][3]+1;
			for(my $i=2;$i<$num;$i++)
			{
				$len[$i] = $cufflinks[$i][4]-$cufflinks[$i][3]+1;
				$len[$i] +=$len[$i-1];
			}
			for(my $i=1;$i<$num;$i++)
		  {
		  	if($len[$i] == $start+1)
		  	{
		  		$cufflinks[$i][3]=$cufflinks[$i][4];
		  		$start_pos = $i;
		  		last;
		  	}
			  if($len[$i] > $start+1)
			  {
			  	$cufflinks[$i][3]=$cufflinks[$i][4]-($len[$i]-($start+1));
			  	$start_pos=$i;
			  	last;
			  }
		  }
      my $frame = 0;
		  for(my $k=$start_pos;$k<$num;$k++)
		  {
				$cufflinks[$k][7]=$frame;	
			  if($len[$k] >= ($stop+1))
			  {
				 	$cufflinks[$k][4] -= ($len[$k]-($stop+1));
					$stop_pos = $k;
				  last;
			  }
			  my $exon_len=$cufflinks[$k][4]-$cufflinks[$k][3]+1;
			  $frame=(3-(($exon_len-$cufflinks[$k][7])%3))%3;
			}
			$cufflinks[0][3] = $cufflinks[$start_pos][3];
		  $cufflinks[0][4] = $cufflinks[$stop_pos][4];
		}
		else
		{
			$start_pos = $num - 1;
			$stop_pos = 1;
			my @len;
			$len[$num-1] = $cufflinks[$num-1][4]-$cufflinks[$num-1][3]+1;
			for(my $i=$num-2;$i>=1;$i--)
			{
				$len[$i] = $cufflinks[$i][4]-$cufflinks[$i][3]+1;
				$len[$i] += $len[$i+1];
			}
			for(my $i=$num-1;$i>=1;$i--)
		  {
		  	if($len[$i] == $start+1)
		  	{
		  		$cufflinks[$i][4]=$cufflinks[$i][3];
		  		$start_pos = $i;
		  		$info_flag = 1;
		  		last;
		  	}
			  if($len[$i] > $start+1)
			  {
			  	$cufflinks[$i][4]=$cufflinks[$i][3]+($len[$i]-($start+1));
			  	$start_pos=$i;
			  	last;
			  }
		  }
      my $frame = 0;
		  for(my $k=$start_pos;$k>=1;$k--)
		  {
				$cufflinks[$k][7]=$frame;		
			  if($len[$k] > ($stop+1))
			  {
					$cufflinks[$k][3]+=($len[$k]-($stop+1));
					$stop_pos = $k;
				  last;
			  }
			  my $exon_len=$cufflinks[$k][4]-$cufflinks[$k][3]+1;
			  $frame=(3-(($exon_len-$cufflinks[$k][7])%3))%3;
			}
			$cufflinks[0][3] = $cufflinks[$stop_pos][3];
		  $cufflinks[0][4] = $cufflinks[$start_pos][4];
		}
		
		for(my $j=0; $j<=8; ++$j)
		{
			print CUFFLINKS_GFF $cufflinks[0][$j]."\t";
		}
		print CUFFLINKS_GFF "\n";
		if( $cufflinks[0][6] eq '+')
		{
		  for(my $i=$start_pos; $i<=$stop_pos; ++$i)
		  {
			if($cufflinks[$i][4]>$cufflinks[$i][3]){##add by lixin
			  for(my $j=0; $j<=8; ++$j)
			  {
			 	  print CUFFLINKS_GFF $cufflinks[$i][$j]."\t";
			}
			  print CUFFLINKS_GFF "\n";
			}
		}
		}
		else
		{
			for(my $i=$start_pos; $i>=$stop_pos; --$i)
		  {
			if($cufflinks[$i][4]>$cufflinks[$i][3]){##add by lixin
			  for(my $j=0; $j<=8; ++$j)
			  {
			 	  print CUFFLINKS_GFF $cufflinks[$i][$j]."\t";
			  }
			  print CUFFLINKS_GFF "\n";
			}
		}
		}
		
		# ORF_list and CDS sequence
		$exons = substr($exons, $start, $stop - $start + 1);
		my $info = "";
		#$info .= 'START_ERR;' unless $exons =~ /^ATG/;
		#$info .= 'STOP_ERR;' unless $exons =~ /(TAG|TGA|TAA)$/;
		if( ($stop - $start +1) < $cutoff )
		{
			$info .= "SHORT_ERR;"
		}
		if( $score <= -15 )
		{
			$info .= "SCORE_ERR";
		}
		if( ($info eq "") && ($info_flag == 1) )
		{
			$info .= "SCORE_ERR;SHORT_ERR;";
		}
		
		print CDS ">$transcript_id\n";
		print CDS $exons."\n";
		my $offset_start = $start + 1;
		my $offset_stop = $stop + 1;
		print ORF "$transcript_id\t$offset_start\t$offset_stop\t$cufflinks[0][6]\t$info\n";
		
	}
}

close(CUFFLINKS_GFF);
close(CDS);
close(ORF);


# ==========================================
# Read markov parameter into markov_parameter
sub read_markov_parameter{
	my $parameter_file=shift;
	open (PARA_FILE,"<$parameter_file") ||die("can not open parameter file!\n");
	my %markov_parameter;
	my $line;
	my $i;
	my $j;
	
	# Read head of markov initial probability matrix
	$line=<PARA_FILE>;
	while(defined($line))
	{
		chomp($line);
		if($line=~/Markov_Initial_probability_matrix/)
		{
			last;
		}
		$line=<PARA_FILE>;
	}
	
	# Read markov initial probability matrix into markov_parameter
	for($i=0; $i<$initial_bases; ++$i)
	{
		for($j=0; $j<3;++$j)
		{
			chomp($line=<PARA_FILE>);
		  my @line=split / /,$line;
		  if($i==$line[1])
		  {
		  	if($j==$line[2])
		  	{
		  		$markov_parameter{$line[0]}->[$line[2]]=[$line[1],$line[3]];
		  	}
		  	else{die("Reading wrong!\n");}
		  }
		  else{die("Reading wrong!\n");}
		}
	}
	
	# Read head of markov transition probability matrix
	$line=<PARA_FILE>;
	while(defined($line))
	{
		chomp($line);
		if($line=~/Markov_Transition_probability_matrix/)
		{
			last;
		}
		$line=<PARA_FILE>;
	}
	
	# Read markov transition probability matrix into markov_parameter
	for($i=0; $i<$transition_bases; ++$i)
	{
		for($j=0; $j<3;++$j)
		{
			chomp($line=<PARA_FILE>);
		  my @line=split / /,$line;
		  if($i==$line[1])
		  {
		  	if($j==$line[2])
		  	{
		  		$markov_parameter{$line[0]}->[$line[2]]=[$line[1],$line[3]];
		  	}
		  	else{die("Reading wrong!\n");}
		  }
		  else{die("Reading wrong!\n");}
		}
	}
	
	close(PARA_FILE);
	return %markov_parameter;
}

#======================================================
sub score_markov_probability{
	my @exons = @_;
	if( $#exons <= 5 ){return -9999};
	my $score = $markov{$exons[0].$exons[1].$exons[2].$exons[3].$exons[4]}[0][1];
	for(my $i=0; $i<= $#exons - 5; ++$i)
	{
		$score += $markov{$exons[$i].$exons[$i+1].$exons[$i+2].$exons[$i+3].
			             $exons[$i+4].$exons[$i+5]}[(0+$i)%3][1];
	}
	return $score;
}

#======================================================

sub calculate_phrase{
	my @y=sort{$a<=>$b}@_;
	my $i;
	for($i=0; $i<=$#_; ++$i)
	{
		if($y[$#_] == $_[$i]){last;}
	}
	return ($_[$i],$i);
	
}
