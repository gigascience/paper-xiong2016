#!/usr/bin/perl -w
### zhushilin@genomics.org.cn
### 2011-05-18
use strict;

die "Usage: $0  <transcripts.gtf>  <transcripts.complete.orf.gtf>  [<min_len>|1]\n" if @ARGV<2;

my ($gtf,$orf,$min)=@ARGV;
$min=1 if (!defined $min || $min<1);

my %infor;
read_gtf($gtf,\%infor,"g");
read_gtf($orf,\%infor,"o");

my %UTR;
foreach my $id (keys %infor){
	next unless (exists $infor{$id}{'o'} && exists $infor{$id}{'g'});
	my $utr1= ($infor{$id}{'o'}{'str'} eq '+')? 5 : 3;
	my $utr2= ($infor{$id}{'o'}{'str'} eq '+')? 3 : 5;
	if ($infor{$id}{'o'}{'sta'}-$infor{$id}{'g'}{'sta'}>=$min){
		$UTR{$id}{$utr1}{'s'}=$infor{$id}{'g'}{'sta'};
		$UTR{$id}{$utr1}{'e'}=$infor{$id}{'o'}{'sta'}-1;
	}
	if ($infor{$id}{'g'}{'end'}-$infor{$id}{'o'}{'end'}>=$min){
		$UTR{$id}{$utr2}{'s'}=$infor{$id}{'o'}{'end'}+1;
		$UTR{$id}{$utr2}{'e'}=$infor{$id}{'g'}{'end'};
	}
}

open ORF,$orf || die "Fail to open file: $orf\n";
#open OUT,">transcripts.complete.orf.gtf.new";
my %temp;
while (<ORF>){
	my @tem=split /\t/;
	my ($Tid)=$tem[8]=~/transcript_id "([^";]+)";/;
	if ($tem[2] eq 'transcript'){
		if (exists $temp{'sca'} && exists $UTR{$temp{'tid'}}{'3'}){
			print "$temp{'sca'}\tCufflinks\tUTR_3\t$UTR{$temp{'tid'}}{'3'}{'s'}\t$UTR{$temp{'tid'}}{'3'}{'e'}\t$temp{'sco'}\t$temp{'str'}\t.\ttranscript_id \"$temp{'tid'}\"; \n";
		}

		$tem[3]=$infor{$Tid}{'g'}{'sta'};
		$tem[4]=$infor{$Tid}{'g'}{'end'};
		my $print_t=join "\t",@tem;
		print $print_t;

		if (exists $UTR{$Tid}{'5'}){
			print "$tem[0]\tCufflinks\tUTR_5\t$UTR{$Tid}{'5'}{'s'}\t$UTR{$Tid}{'5'}{'e'}\t$tem[5]\t$tem[6]\t.\ttranscript_id \"$Tid\"; \n";
		}

		$temp{'tid'}=$Tid;
		$temp{'sca'}=$tem[0];
		$temp{'sco'}=$tem[5];
		$temp{'str'}=$tem[6];
	}elsif ($tem[2] eq 'exon'){
		print $_;
	}
}
if (exists $UTR{$temp{'tid'}}{'3'}){
	print "$temp{'sca'}\tCufflinks\tUTR_3\t$UTR{$temp{'tid'}}{'3'}{'s'}\t$UTR{$temp{'tid'}}{'3'}{'e'}\t$temp{'sco'}\t$temp{'str'}\t.\ttranscript_id \"$temp{'tid'}\"; \n";
}
close ORF;




#================================================================================
sub read_gtf{
	my ($this_file,$this_hash,$this_tag)=@_;
	
	open FILE,$this_file || die "Fail to open file: $this_file\n";
	while (<FILE>){
		my @te=split /\t/;
		if ($te[2] eq 'transcript'){
			my ($tid)=$te[8]=~/transcript_id "([^";]+)";/;
			unless ($tid){
				warn "Check format of file: $this_file \n$_\n";
				next;
			}
			$this_hash->{$tid}{$this_tag}{'sta'}=$te[3];
			$this_hash->{$tid}{$this_tag}{'end'}=$te[4];
			$this_hash->{$tid}{$this_tag}{'str'}=$te[6];
		}
	}
	close FILE;
}
