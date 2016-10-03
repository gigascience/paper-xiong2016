#!/usr/bin/perl -w
### Zhu Shilin <silenfree@gmail.com>
use strict;
use File::Basename qw(basename);

die "Usage: $0 <transcripts.gtf>  <transcripts.complete.orf.gtf>\n" if @ARGV < 1;

my ( $gtf, $orf ) = @ARGV;

my $UTR5name = "UTR_5";
my $UTR3name = "UTR_3";

( $gtf =~ /\.gz$/ ) ? open GTF, "gzip -dc $gtf|" : open GTF, $gtf;
( $orf =~ /\.gz$/ ) ? open ORF, "gzip -dc $orf|" : open ORF, $orf;

my %info;
my %exon;
my %mrnaloci;
while (<ORF>) {
	next if ( /^#/ || /^\s/ );
	my @te = split /\t/;

	if ( $te[2] eq "exon" ) {
		my ($tid) = $te[8] =~ /transcript_id "([^\";]+)"/;
		$info{ $te[0] }{$tid}{strand}               = $te[6];
		$info{ $te[0] }{$tid}{cds}{ $te[3] }{e}     = $te[4];
		$info{ $te[0] }{$tid}{cds}{ $te[3] }{infor} = $_;
		push @{ $mrnaloci{ $te[0] }{$tid}{locis} }, $te[3];
		push @{ $mrnaloci{ $te[0] }{$tid}{locis} }, $te[4];
	}
}
close ORF;

while (<GTF>) {
	next if ( /^#/ || /^\s/ );
	my @te = split /\t/;

	my ($tid) = $te[8] =~ /transcript_id "([^\";]+)";/;
	if ( $te[2] eq "transcript" ) {
		@{ $mrnaloci{ $te[0] }{$tid}{infor} } = @te;
	}
	elsif ( $te[2] eq "exon" ) {
		$exon{ $te[0] }{$tid}{ $te[3] } = $te[4];
		push @{ $mrnaloci{ $te[0] }{$tid}{locis} }, $te[3];
		push @{ $mrnaloci{ $te[0] }{$tid}{locis} }, $te[4];
	}
}
close GTF;

foreach my $scaf ( keys %info ) {
	foreach my $tid ( keys %{ $info{$scaf} } ) {

		my @mRNA = @{ $mrnaloci{$scaf}{$tid}{infor} };
		( $mRNA[3], $mRNA[4] ) = ( sort { $a <=> $b } @{ $mrnaloci{$scaf}{$tid}{locis} } )[ 0, -1 ];
		$mRNA[6] = $info{$scaf}{$tid}{strand};
		@{ $mrnaloci{$scaf}{$tid}{locis} } = ();

		if ( exists $exon{$scaf}{$tid} ) {
			if ( $info{$scaf}{$tid}{strand} eq "+" || $info{$scaf}{$tid}{strand} eq "." ) {
				my @sort_exon_start = sort { $a <=> $b } keys %{ $exon{$scaf}{$tid} };
				my ( $first_cds_start, $last_cds_start ) = ( sort { $a <=> $b } keys %{ $info{$scaf}{$tid}{cds} } )[ 0, -1 ];
				foreach my $ev_ex_s (@sort_exon_start) {
					if ( $ev_ex_s < $first_cds_start ) {
						my $utr5_tag_end =
						  ( $exon{$scaf}{$tid}{$ev_ex_s} < $first_cds_start )
						  ? $exon{$scaf}{$tid}{$ev_ex_s}
						  : $first_cds_start - 1;
						$info{$scaf}{$tid}{utr5} .=
						  "$scaf\t$mRNA[1]\t$UTR5name\t$ev_ex_s\t$utr5_tag_end\t$mRNA[5]\t$mRNA[6]\t.\ttranscript_id \"$tid\";\n";
					}
					if ( $exon{$scaf}{$tid}{$ev_ex_s} > $info{$scaf}{$tid}{cds}{$last_cds_start}{e} ) {
						my $utr3_tag_start =
						  ( $ev_ex_s > $info{$scaf}{$tid}{cds}{$last_cds_start}{e} )
						  ? $ev_ex_s
						  : $info{$scaf}{$tid}{cds}{$last_cds_start}{e} + 1;
						$info{$scaf}{$tid}{utr3} .= "$scaf\t$mRNA[1]\t$UTR3name\t$utr3_tag_start\t$exon{$scaf}{$tid}{$ev_ex_s}\t"
						  . "$mRNA[5]\t$mRNA[6]\t.\ttranscript_id \"$tid\";\n";
					}
				}
			}
			else {
				my @sort_exon_start = sort { $b <=> $a } keys %{ $exon{$scaf}{$tid} };
				my ( $first_cds_start, $last_cds_start ) = ( sort { $b <=> $a } keys %{ $info{$scaf}{$tid}{cds} } )[ 0, -1 ];
				foreach my $ev_ex_s (@sort_exon_start) {
					if ( $exon{$scaf}{$tid}{$ev_ex_s} > $info{$scaf}{$tid}{cds}{$first_cds_start}{e} ) {
						my $utr5_tag_start =
						  ( $ev_ex_s > $info{$scaf}{$tid}{cds}{$first_cds_start}{e} )
						  ? $ev_ex_s
						  : $info{$scaf}{$tid}{cds}{$first_cds_start}{e} + 1;
						$info{$scaf}{$tid}{utr5} .= "$scaf\t$mRNA[1]\t$UTR5name\t$utr5_tag_start\t$exon{$scaf}{$tid}{$ev_ex_s}\t"
						  . "$mRNA[5]\t$mRNA[6]\t.\ttranscript_id \"$tid\";\n";
					}
					if ( $ev_ex_s < $last_cds_start ) {
						my $utr3_tag_end =
						  ( $exon{$scaf}{$tid}{$ev_ex_s} < $last_cds_start )
						  ? $exon{$scaf}{$tid}{$ev_ex_s}
						  : $last_cds_start - 1;
						$info{$scaf}{$tid}{utr3} .=
						  "$scaf\t$mRNA[1]\t$UTR3name\t$ev_ex_s\t$utr3_tag_end\t$mRNA[5]\t$mRNA[6]\t.\ttranscript_id \"$tid\";\n";
					}
				}
			}
		}
		%{ $exon{$scaf}{$tid} } = ();

		print join( "\t", @mRNA );
		print $info{$scaf}{$tid}{utr5} if ( exists $info{$scaf}{$tid}{utr5} );
		foreach my $cds_s (
			  ( $info{$scaf}{$tid}{strand} eq "-" )
			? ( sort { $b <=> $a } keys %{ $info{$scaf}{$tid}{cds} } )
			: ( sort { $a <=> $b } keys %{ $info{$scaf}{$tid}{cds} } )
		  )
		{
			print $info{$scaf}{$tid}{cds}{$cds_s}{infor};
		}
		print $info{$scaf}{$tid}{utr3} if ( exists $info{$scaf}{$tid}{utr3} );
	}
}

#========================================================================== Step and Cherish ======================================================================
