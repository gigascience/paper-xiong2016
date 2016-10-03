#!/usr/bin/perl -w

# zhangbo@genomics.cn
# create 2010-04-06
# modify liugeng@genomics.cn 2010-10-09 


use strict;
die "$0 <predict.pep> <rnaseq.pep> <tmap>\n" if @ARGV !=3;

my ($fasta1, $fasta2 , $tmap) = @ARGV;

open (MATCH, ">$tmap.match.pep")|| die $!;
open (BLT, ">$tmap.blast.out") || die $!;

my %seqs = read_fasta_to_hash($fasta1);
%seqs = (%seqs, read_fasta_to_hash($fasta2));

compare_tmap($tmap);

#======================================
sub compare_tmap{
	open (FH, shift) || die $!;
	while(<FH>){
		chomp;
		my ($ref_id, $cuff_id) = (split /\t/)[0,1];
		my ($i, $j, $bl_rslt) = ('','','');
		if (exists $seqs{$ref_id} && exists $seqs{$cuff_id}){
			$i = ">$seqs{$ref_id}[0]\n$seqs{$ref_id}[1]\n";
			open (FI, ">$tmap.i");
			print FI $i;
			close FI;
			
			$j = ">$seqs{$cuff_id}[0]\n$seqs{$cuff_id}[1]\n";
			open (FJ, ">$tmap.j");
			print FJ $j;
			close FJ;
			$bl_rslt =  `/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_annotation/Annotation_v2.0/bin/Annotation_pipeline1_1.0/05.RNAseq_annlysis/bin/bl2seq -F F -p blastp -m T -W 7 -i $tmap.i -j $tmap.j`;
#			$bl_rslt =  `/ifs2/BC_GAG/Bin/Annotation/bin/Annotation_pipeline1_1.0/05.RNAseq_annlysis/bin/bl2seq -F F -p blastp -m T -W 7 -i $tmap.i -j $tmap.j`;

			unless ($bl_rslt =~ /No hits found/){
				print BLT $bl_rslt;
				my ($q_len) = $bl_rslt =~ /\(([\d,]+) letters/;
				my ($s_len) = $bl_rslt =~ /Length = ([\d,]+)/;
				my ($total, $percent) = $bl_rslt =~ /Identities =\s+\d+\/(\d+)\s+\((\S+)\%\),\s+/;
				$q_len =~ s/\D//g;
				$s_len =~ s/\D//g;
				
				if ($percent == 100 && $q_len == $s_len){
					# identical match
					print MATCH "#IDENTICAL\t$ref_id\t.\t$cuff_id\n$i$j\n";	
				}elsif($percent >= 95 && ($total / $s_len >= 0.9 || $total / $q_len >= 0.9) ){
					print MATCH "#PART\t$ref_id\t.\t$cuff_id\n$i$j\n";
				}
			}
		}
	}
	
	system("/bin/rm $tmap.i $tmap.j");
	system("date");
	close FH;
}

#=================================
sub read_fasta_to_hash{
	my %seqHash = ();
	open (FH , shift) || die $!;
	local $/ = "\n>";
	while (<FH>){
		chomp;
		s/^>//;
		my ($tag , $seq) = split (/\n/, $_, 2);
		my ($id)  = $tag =~ /^(\S+)/; 
		$seq =~ s/\W//g;
		$seqHash{$id} = [$tag,$seq];
	}
	$/ = "\n"; 
	close FH;
	return %seqHash;
}

