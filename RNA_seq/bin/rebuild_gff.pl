#!/usr/bin/perl -w
=head1

=head1 Find overlap in gff and gtf files or check overlap in one file.

=head1 By Zhu Shilin <zhushilin@genomics.org.cn> 2011-05-16

=head1 Usage: perl rebuild_gff.pl  <gff>  [<gff>] ...  [--options]

   --gtf  <str>			gtf files
   --printgff <str>		output rebuild gff (yes/no, defoult "yes")
     --out <str>		name of rebuild gff (default "rebuild.gff")
   --conside_strand <str>	conside strand while find overlap (yes/no, default "no")
   --verbose			get verbose information
   --usage			simple usage of this program
   --help			print help information to screen
=cut
use strict;
#use Fatal;
#use File::Basename;
use Getopt::Long;

my (@gtf,$output,$printgff,$conside_str,$verbose,$usage,$help);
GetOptions(
	"gtf:s"=>\@gtf,
	"out:s"=>\$output,
	"printgff:s"=>\$printgff,
	"conside_strand:s"=>\$conside_str,
	"verbose!"=>\$verbose,
	"usage!"=>\$usage,
	"help!"=>\$help
);
$printgff ||="yes";
$output ||="rebuild.gff";
$conside_str ||="no";
$printgff=($printgff=~/^y/i)? 1 : 0;

die `pod2text $0` if ($help);
die "Usage: $0  <gff>  [<gff>] ...   [--gtf <gtf>]  [--gtf <gtf>] ... 
	[--printgff (yes)/no]  [--out <out.gff>]  [--conside_strand (no)/yes]  [--verbose]
You can type '--help' to read detailed information about this program.\n" if ( (@ARGV<1 && @gtf==0) || $usage);

my @gff=@ARGV;

if ($verbose){
	print "Start\n";
	system "date";
	print "\n";
}

my %RNA;
my %element;
### Read gff
print "Read gff file...\n" if ($verbose && @gff!=0);
my $old_gene_num=0;
foreach my $gff (@gff){
	my $this_gff_num=0;
	my $this_gff_gnum=0;
	open GFF,$gff || die "Fail to open $gff\n";
	while (<GFF>){
		next if (/^\s/ || /^#/);
		my @te=split /\t/,$_;
		if($te[2] eq 'gene'){
			$this_gff_gnum++;
		}elsif($te[2] eq 'mRNA'){
			my ($R_id)=$te[8]=~/ID=([^;\s]+)/i;
			die "$_\nWrong format! Last row for mRNA should be \"ID=*\".\n" unless ($R_id);
			my $R_str=($conside_str=~/^n/i)? "." : $te[6];
			$RNA{$te[0]}{$R_str}{$R_id}{'s'}=$te[3];
			$RNA{$te[0]}{$R_str}{$R_id}{'e'}=$te[4];
			$RNA{$te[0]}{$R_str}{$R_id}{'infor'}=$_;
			#$RNA{$te[0]}{$R_str}{$R_id}{'str'}=$te[6];
			$this_gff_num++;
		}elsif( $printgff && ($te[2] eq 'CDS' || $te[2]=~/UTR/i || $te[2] eq 'exon' || $te[2] eq 'intron') ){
			my ($this_parent)=$te[8]=~/Parent=([^;\s]+)/i;
			my @parent_list=split /,/,$this_parent;
			foreach my $every_pa (@parent_list){
				my $deal_infor=$_;
				$deal_infor=~s/Parent=([^;\s]+)/Parent=$every_pa/;
				push @{$element{$every_pa}{'ele'}},$deal_infor;
			}
		}
	}
	$this_gff_gnum=$this_gff_num if $this_gff_gnum==0;
	$old_gene_num+=$this_gff_gnum;
	close GFF;
	print "gene number: $this_gff_gnum\t$gff\n";
	print "mRNA number: $this_gff_num\t$gff\n";
}

### Read gtf
print "Read gtf file...\n" if ($verbose && @gtf!=0);
my $total_trans_num=0;
foreach my $gtf (@gtf){
	my $this_gtf_num=0;
	open GTF,$gtf || die "Fail to open $gtf \n";
	while (<GTF>){
		next if (/^\s/ || /^#/);
		my @temp=split /\t/,$_;
		if ($temp[2] eq 'transcript'){
			my ($T_id)=$temp[8]=~/transcript_id\s+"([^";]+)";/;
			die "$_\nWrong format! Last row for transcript should be \"transcript_id *\".\n" unless ($T_id);
			my $T_str=($conside_str=~/^n/i)? "." : $temp[6];
			$RNA{$temp[0]}{$T_str}{$T_id}{'s'}=$temp[3];
			$RNA{$temp[0]}{$T_str}{$T_id}{'e'}=$temp[4];
			$RNA{$temp[0]}{$T_str}{$T_id}{'infor'}=$_;
			#$RNA{$temp[0]}{$T_str}{$T_id}{'str'}=$temp[6];
			$this_gtf_num++;
		}elsif( $printgff && ($temp[2] eq 'exon' || $temp[2] eq 'intron' || $temp[2] eq 'CDS' || $temp[2]=~/UTR/i) ){
			my ($T_id)=$temp[8]=~/transcript_id\s+"([^";]+)";/;
			push @{$element{$T_id}{'ele'}},$_;
		}
	}
	$total_trans_num+=$this_gtf_num;
	close GTF;
	print "transcript number: $this_gtf_num\t$gtf\n";
}


### Deal loci
print "Deal loci...\n" if ($verbose);
if ($printgff){
	open OUT,">$output" || die "Fail to open output file: $output\n";
}
open OVE,">overlap.infor"   || die "Fail to open output file: overlap.infor\n";
my $gene_num=0;
foreach my $scaf (sort keys %RNA){
	foreach my $str (keys %{$RNA{$scaf}}){
		my @sort_id=sort { $RNA{$scaf}{$str}{$a}{'s'} <=> $RNA{$scaf}{$str}{$b}{'s'} } keys %{$RNA{$scaf}{$str}}; 
		my $id_all_num=scalar @sort_id;

		if ($id_all_num==1){
			my $only_one_R=shift @sort_id;
			$gene_num++;
			print OVE ">gene$gene_num\t1\n";
			print OVE "$only_one_R\n";
			print OUT "$scaf\tBGI\tgene\t$RNA{$scaf}{$str}{$only_one_R}{'s'}\t$RNA{$scaf}{$str}{$only_one_R}{'e'}\t.\t.\t.\tID=$gene_num;\n" if ($printgff);
			my $print_only_one_gene=$RNA{$scaf}{$str}{$only_one_R}{'infor'};
			$print_only_one_gene=~s/Parent=([^;\s]+)/Parent=$gene_num/i;
			print OUT "$print_only_one_gene"     if ($printgff);
			print_element($only_one_R,\%element) if ($printgff);		###########
			next;
		}

		my $RNA_num_control=0;
		my $one_sta;
		my $one_end;
		my @print_R;
		foreach my $each_Rid (@sort_id){
			$RNA_num_control++;
			my $this_start=$RNA{$scaf}{$str}{$each_Rid}{'s'};
			my $this_end=$RNA{$scaf}{$str}{$each_Rid}{'e'};
			unless (defined $one_end){
				$one_sta=$RNA{$scaf}{$str}{$each_Rid}{'s'};
				$one_end=$RNA{$scaf}{$str}{$each_Rid}{'e'};
				push @print_R,$each_Rid;
				next;
			}

			if ($this_start>$one_end && $RNA_num_control<$id_all_num){
				$gene_num++;
				#print OUT "#\n" if ($printgff && @print_R>1);
				print OUT "$scaf\tBGI\tgene\t$one_sta\t$one_end\t.\t.\t.\tID=$gene_num;\n" if ($printgff);
				my $mRNA_num_in_gene=scalar @print_R;
				print OVE ">gene$gene_num\t$mRNA_num_in_gene\n";
				print_RNA(\@print_R,$scaf,$str,$gene_num);
				$one_sta=$this_start;
				@print_R=($each_Rid);
			}elsif($this_start<=$one_end && $RNA_num_control<$id_all_num){
				push @print_R,$each_Rid;
			}elsif($this_start<=$one_end && $RNA_num_control==$id_all_num){
				push @print_R,$each_Rid;
				$one_end=$this_end if ($one_end<$this_end);
				$gene_num++;
				#print OUT "#\n" if ($printgff && @print_R>1);
				my $mRNA_num_in_gene=scalar @print_R;
				print OVE ">gene$gene_num\t$mRNA_num_in_gene\n";
				print OUT "$scaf\tBGI\tgene\t$one_sta\t$one_end\t.\t.\t.\tID=$gene_num;\n" if ($printgff);
				print_RNA(\@print_R,$scaf,$str,$gene_num);
			}elsif($this_start>$one_end && $RNA_num_control==$id_all_num){
				$gene_num++;
				#print OUT "#\n" if ($printgff && @print_R>1);
				my $mRNA_num_in_gene=scalar @print_R;
				print OVE ">gene$gene_num\t$mRNA_num_in_gene\n";
				print OUT "$scaf\tBGI\tgene\t$one_sta\t$one_end\t.\t.\t.\tID=$gene_num;\n" if ($printgff);
				print_RNA(\@print_R,$scaf,$str,$gene_num);

				$gene_num++;
				print OVE ">gene$gene_num\t1\n";
				print OUT "$scaf\tBGI\tgene\t$this_start\t$this_end\t.\t.\t.\tID=$gene_num;\n" if ($printgff);
				@print_R=($each_Rid);
				print_RNA(\@print_R,$scaf,$str,$gene_num);
			}
			$one_end=$this_end if ($one_end<$this_end);
		}
	}
}
close OUT if ($printgff);
close OVE;
print "Old gene number: $old_gene_num\n" if $old_gene_num;
print "Total transcripts number: $total_trans_num\n" if $total_trans_num;
print "New gene number: $gene_num\n\n";
if ($verbose){
	print "DONE!\n";
	system "date";
}

#=========================================================================================
sub print_RNA{
	my ($this_all_R,$this_sca,$this_str,$this_gene_num)=@_;
	foreach my $this_every_R (@$this_all_R){
		print OVE "$this_every_R\n";
		if ($printgff){
			my $print_every_R=$RNA{$this_sca}{$this_str}{$this_every_R}{'infor'};
			$print_every_R=~s/Parent=([^;\s]+)/Parent=$this_gene_num/;
			print OUT "$print_every_R";
			print_element($this_every_R,\%element);
		}
	}
}

#==================================
sub print_element{
	my ($R_name,$e_hash)=@_;

	foreach my $print_e ( @{$e_hash->{$R_name}{'ele'}} ){
		print OUT $print_e;
	}
}
