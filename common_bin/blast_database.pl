#!/usr/bin/perl

=head1 Name

blast_database.pl  --  the pipeline to run blast

=head1 Description

this program is the pipeline to run blast, note that the database
must be formated previously, and the version of formatdb and 
blastall must be same. Also note that protein database and DNA database
can't be used at the same time. 

The option "--program" can only be co-used with "--userdb".

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 2.0,  Date: 2008-5-7
  Note:

=head1 Usage
  
  perl blast_database.pl [options] <sequences.fa>
  --evalue       set the E-val for blast, default 1e-5
  --nr           search against Nr database
  --nt           search against Nt database
  --swissprot    search against SwissProt database	
  --trembl       search against TrEMBL database
  --cog          search against Cog database
  --kegg         search against	Kegg database
  --userdb <str> search against user defined database
  --program      set the program of blastall, only for userdb, default blastp
  --cuts <int>   set the sequence number in each cutted query file, default=100
  --cpu <int>	 set the cpu number to use in parallel, default=3   
  --run <str>    set the parallel type, qsub, or multi, default=qsub
  --resource <str>       set the resource needed, default vf=0.9G
  --node <str>   set the compute node, default h=compute-0-150
  --queue <str>  set the queue
  --outdir <str>  set the result directory, default="."
  --f   <F/T>  set the -F for blast, default 'F';
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl ../bin/blast_database.pl --nr --swissprot --trembl  --cog  --kegg  -cpu 10 ../input/rice_prot100.fa
  perl ../bin/blast_database.pl --nt ../input/rice_cds.fa
  perl ../bin/blast_database.pl --program blastp --userdb ../input/protein_database.fa ../input/rice_protein.fa

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

my ($Node,$Program,$Evalue,$Nr,$Nt,$Swissprot,$TrEMBL,$Cog,$Kegg,$Userdb,$Cuts,$Cpu,$Run,$Outdir,$Resource);
my ($Verbose,$Help,$Queue);
my $F;
GetOptions(
	"program:s"=>\$Program,
	"evalue:s"=>\$Evalue,
	"nr"=>\$Nr,
	"nt"=>\$Nt,
	"swissprot"=>\$Swissprot,
	"trembl"=>\$TrEMBL,
	"cog"=>\$Cog,
	"kegg"=>\$Kegg,
	"userdb:s"=>\$Userdb,
	"cuts:i"=>\$Cuts,
	"cpu:i"=>\$Cpu,
	"run:s"=>\$Run,
	"f:s"=>\$F,
        "resource:s"=>\$Resource,
        "node:s"=>\$Node,
	"queue:s"=>\$Queue,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Program ||= "blastp";
$Evalue ||= 1e-5;
$Cuts ||= 100;
$Cpu ||= 3;
$Run ||= "qsub";
$Outdir ||= ".";
$F ||= 'F';
$Resource ||= "vf=0.9G";
#$Node ||= "h=compute-0-189";
my $Node_para=(defined $Node)?"-node $Node":"";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $seq_file = shift;
my $seq_file_name = basename($seq_file);

my %config;
parse_config("$Bin/../../config.txt",\%config);

my $fastaDeal = $config{"fastaDeal.pl"};
my $blast = $config{"blastall"};
my $blast_parser = $config{"blast_parser"};
my $qsub_sge = $config{"qsub_sge.pl"};
my $multi_process = $config{"multi-process.pl"};
print "$qsub_sge\n";

my $nr_path = "$config{nr}/nr";
my $nt_path = "$config{nt}/nt";
my $swissprot_path = $config{swissprot};
my $trembl_path = $config{trembl};
my $cog_path = $config{cog};
my $kegg_path = $config{kegg};
my $userdb;
$userdb = basename($Userdb) if(defined $Userdb);

my $blast_shell_file = "$Outdir/$seq_file_name.blast.sh";
my @subfiles;

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

`perl $fastaDeal -cuts $Cuts $seq_file -outdir $Outdir`;
@subfiles = glob("$Outdir/$seq_file_name.cut/*.*");

##creat shell file
open OUT,">$blast_shell_file" || die "fail $blast_shell_file";
foreach my $subfile (@subfiles) {
	print OUT "$blast -b 10000 -v 10000 -p $Program -e $Evalue -F $F -d $nr_path -i $subfile -o $subfile.nr.blast; \n" if(defined $Nr);
	print OUT "$blast -b 10000 -v 10000 -p $Program -e $Evalue -F $F -d $nt_path -i $subfile -o $subfile.nt.blast; \n" if(defined $Nt);
	print OUT "$blast -b 10000 -v 10000 -p $Program -e $Evalue -F $F -d $swissprot_path -i $subfile -o $subfile.swissprot.blast; \n" if(defined $Swissprot);
	print OUT "$blast -b 10000 -v 10000 -p $Program -e $Evalue -F $F -d $trembl_path -i $subfile -o $subfile.trembl.blast; \n" if(defined $TrEMBL);
	print OUT "$blast -b 10000 -v 10000 -p $Program -e $Evalue -F $F -d $cog_path -i $subfile -o $subfile.cog.blast; \n" if(defined $Cog);
	print OUT "$blast -b 10000 -v 10000 -p $Program -e $Evalue -F $F -d $kegg_path -i $subfile -o $subfile.kegg.blast; \n" if(defined $Kegg);
	print OUT "$blast  -b 10000 -v 10000 -p $Program -e $Evalue -F $F  -d $Userdb -i $subfile -o $subfile.$userdb.blast; \n" if(defined $Userdb);
}
close OUT;

##run the shell file
my $opt;
$opt=" -queue $Queue " if (defined $Queue);
`perl $qsub_sge $opt  --maxjob $Cpu  --resource $Resource $Node_para $blast_shell_file` if ($Run eq "qsub");
`perl $multi_process -cpu $Cpu $blast_shell_file` if ($Run eq "multi");


##cat together the result
`cat $Outdir/$seq_file_name.cut/*.nr.blast > $Outdir/$seq_file_name.nr.blast;` if(defined $Nr);
`cat $Outdir/$seq_file_name.cut/*.nt.blast > $Outdir/$seq_file_name.nt.blast; ` if(defined $Nt);
`cat $Outdir/$seq_file_name.cut/*.swissprot.blast > $Outdir/$seq_file_name.swissprot.blast; ` if(defined $Swissprot);
`cat $Outdir/$seq_file_name.cut/*.trembl.blast > $Outdir/$seq_file_name.trembl.blast; ` if(defined $TrEMBL);
`cat $Outdir/$seq_file_name.cut/*.cog.blast > $Outdir/$seq_file_name.cog.blast; ` if(defined $Cog);
`cat $Outdir/$seq_file_name.cut/*.kegg.blast > $Outdir/$seq_file_name.kegg.blast; ` if(defined $Kegg);
`cat $Outdir/$seq_file_name.cut/*.$userdb.blast > $Outdir/$seq_file_name.$userdb.blast; ` if(defined $Userdb);

######################################  deal with the result  #########################################################################################
my  $blast_passer_shell_file = "$Outdir/$seq_file_name.blast.parser.sh";
##covert to table format
open (OUT2,">$blast_passer_shell_file") || die "fail:$blast_passer_shell_file";
print OUT2 "perl $blast_parser -nohead $Outdir/$seq_file_name.nr.blast > $Outdir/$seq_file_name.nr.blast.tab;\n" if(defined $Nr);
print OUT2 "perl $blast_parser -nohead  $Outdir/$seq_file_name.nt.blast > $Outdir/$seq_file_name.nt.blast.tab;\n" if(defined $Nt);
print OUT2 "perl $blast_parser -nohead  $Outdir/$seq_file_name.swissprot.blast > $Outdir/$seq_file_name.swissprot.blast.tab;\n" if(defined $Swissprot);
print OUT2 "perl $blast_parser -nohead  $Outdir/$seq_file_name.trembl.blast > $Outdir/$seq_file_name.trembl.blast.tab;\n" if(defined $TrEMBL);
print OUT2 "perl $blast_parser -nohead  $Outdir/$seq_file_name.cog.blast > $Outdir/$seq_file_name.cog.blast.tab;\n" if(defined $Cog);
print OUT2 "perl $blast_parser -nohead  $Outdir/$seq_file_name.kegg.blast > $Outdir/$seq_file_name.kegg.blast.tab;\n" if(defined $Kegg);
print OUT2 "perl $blast_parser -nohead  $Outdir/$seq_file_name.$userdb.blast > $Outdir/$seq_file_name.$userdb.blast.tab;\n" if(defined $Userdb);


##get the the best hit
print OUT2 "perl $blast_parser -nohead -tophit 1 -topmatch 1   $Outdir/$seq_file_name.nr.blast > $Outdir/$seq_file_name.nr.blast.tab.best;\n" if(defined $Nr);
print OUT2 "perl $blast_parser -nohead -tophit 1 -topmatch 1   $Outdir/$seq_file_name.nt.blast > $Outdir/$seq_file_name.nt.blast.tab.best;\n" if(defined $Nt);
print OUT2 "perl $blast_parser -nohead -tophit 1 -topmatch 1   $Outdir/$seq_file_name.swissprot.blast > $Outdir/$seq_file_name.swissprot.blast.tab.best;\n" if(defined $Swissprot);
print OUT2 "perl $blast_parser -nohead -tophit 1 -topmatch 1   $Outdir/$seq_file_name.trembl.blast > $Outdir/$seq_file_name.trembl.blast.tab.best;\n" if(defined $TrEMBL);
print OUT2 "perl $blast_parser -nohead -tophit 1 -topmatch 1   $Outdir/$seq_file_name.cog.blast > $Outdir/$seq_file_name.cog.blast.tab.best;\n" if(defined $Cog);
print OUT2 "perl $blast_parser -nohead -tophit 1 -topmatch 1   $Outdir/$seq_file_name.kegg.blast > $Outdir/$seq_file_name.kegg.blast.tab.best;\n" if(defined $Kegg);
print OUT2 "perl $blast_parser -nohead -tophit 1 -topmatch 1   $Outdir/$seq_file_name.$userdb.blast > $Outdir/$seq_file_name.$userdb.blast.tab.best;\n" if(defined $Userdb);

close OUT2;

`perl $qsub_sge   --maxjob 2 $Node_para  $blast_passer_shell_file` if ($Run eq "qsub");
`perl $multi_process  -cpu 2  $blast_passer_shell_file` if ($Run eq "multi");
####################################################
################### Sub Routines ###################
####################################################


##parse the software.config file, and check the existence of each software
####################################################
sub parse_config{
	my $conifg_file = shift;
	my $config_p = shift;
	
	my $error_status = 0;
	open IN,$conifg_file || die "fail open: $conifg_file";
	while (<IN>) {
		next if(/#/);
		if (/(\S+)\s*=\s*(\S+)/) {
			my ($software_name,$software_address) = ($1,$2);
			$config_p->{$software_name} = $software_address;
			if (! -e $software_address){
				warn "Non-exist:  $software_name  $software_address\n"; 
				$error_status = 1;
			}
		}
	}
	close IN;
	die "\nExit due to error of software configuration\n" if($error_status);
}


