#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../";
use SNPiR::config;

################################################################################
################################################################################
# $Revision: $
# Authors: Robert Piskol ( piskol@stanford.edu ), Gokul Ramaswami ( gokulr@stanford.edu )
# Last modification $Author: piskol $
#perl script that runs converts reads mapping to splice reference back to genomic coordinates
my ($INPUTFILE,$OUTPUTFILE,$BAMFILE,$ILLUMINAQUALS,$HELP);
my $QUALOFFSET = 33;
my $MINBASEQUAL = 25;
parse_command_line();

sub parse_command_line {
	my $help;
	
	usage() if (scalar @ARGV==0);

	&GetOptions(
	    "infile=s" => \$INPUTFILE,
	    "outfile=s" => \$OUTPUTFILE,
	    "bamfile=s" => \$BAMFILE,
	    "minbasequal=s" => \$MINBASEQUAL,
	    "illumina" => \$ILLUMINAQUALS,
	    "help" => \$HELP
	    );
	
	usage() if ($help);
	   
	if(!$INPUTFILE || !$OUTPUTFILE || !$BAMFILE){
	    usage();
	}
	if($ILLUMINAQUALS){
		$QUALOFFSET = 64;
	}
}


my $minmismatch = 1;


open (my $INPUT , "<", $INPUTFILE) or die "error opening inputfile: $!\n";
open (my $OUTPUT, ">", $OUTPUTFILE) or die "error opening outputfile: $!\n";
open (my $OUTPUT_failed, ">", $OUTPUTFILE."_failed") or die "error opening outputfile: $!\n";

while (<$INPUT>) {
	chomp;
	my @fields = split;
	my $TEMPNAME = join '', $OUTPUTFILE,'_tmp';
	my ($chr, $position) = ($fields[0], $fields[1]);
	my $bamposition = join '', $chr,':',$position,'-',$position;
	system("$SAMTOOLSEXE view $BAMFILE $bamposition > $TEMPNAME");

	my $editnuc = $fields[4];
	my ($newcov, $newmismatch) = (0,0);
	my $basequalFail;
	my $readPosFail;
	
	open(my $TMPFILE, "<", $TEMPNAME);
	while (<$TMPFILE>) {
		#print $_;
		chomp;
		my @bamfields = split;
		my ($alignment, $readstart, $cigar, $sequence, $qualities) = ($bamfields[1], $bamfields[3], $bamfields[5], $bamfields[9], $bamfields[10]);
		my @sequencebases = split(//,$sequence);
		my @qualscores = split(//,$qualities);
		my ($currentpos, $readpos) = ($readstart,1);
		my $base_readpos;
		my @cigarnums = split(/[MIDNSHP]/, $cigar);
		my @cigarletters = split(/[0-9]+/,$cigar);
		shift @cigarletters;

		for (my $i = 0; $i < @cigarnums; $i++) {
			if ($cigarletters[$i] =~ m/[ISH]/) {
				$readpos = $readpos + $cigarnums[$i];
			}
			elsif ($cigarletters[$i] =~ m/[DN]/) {
				$currentpos = $currentpos + $cigarnums[$i];
			}
			elsif ($cigarletters[$i] =~ m/M/) {
				for (my $j = 0; $j < $cigarnums[$i]; $j++) {
					$base_readpos = $readpos if ($currentpos == $position);
					$currentpos++;
					$readpos++;	
				}	
			}
		}
		next unless $base_readpos;
		my $revstrand = 0;
		$revstrand = 1 if ($alignment & 16);
		if (($revstrand == 0 && $base_readpos > 6) || ($revstrand == 1 && $base_readpos < $readpos - 5)) {
			if (ord($qualscores[$base_readpos-1]) >= $MINBASEQUAL+$QUALOFFSET) {
				$newcov++;
				$newmismatch++ if ($sequencebases[$base_readpos-1] eq $editnuc);
				
			}
			else{
				#print "$bamposition failing quality ".ord($qualscores[$base_readpos-1])."\n";
				$basequalFail=1;
			}
		}
		else{
			$readPosFail=1;
			#print "$bamposition failing position\n";
		}
		#print $base_readpos;
	}

	system("rm $TEMPNAME");
	if ($newmismatch >= $minmismatch) {
		my $varfreq = sprintf("%.3f", $newmismatch/$newcov);
		print $OUTPUT "$fields[0]\t$fields[1]\t$newcov,$newmismatch\t$fields[3]\t$fields[4]\t$varfreq";
		for (my $i = 6; $i < @fields; $i++) {
			print $OUTPUT "\t$fields[$i]";
		}
		print $OUTPUT "\n";
	}
	if ($newmismatch < $minmismatch) {
		my $explain;
		$explain .= ' low basequal;' if($basequalFail);
		$explain .= ' mismatch at readend;' if($readPosFail);
		print $OUTPUT_failed "$fields[0]\t$fields[1]\t$newcov,$newmismatch\t$fields[3]\t$fields[4]\tNAN\treason:$explain";
		for (my $i = 6; $i < @fields; $i++) {
			print $OUTPUT_failed "\t$fields[$i]";
		}
		print $OUTPUT_failed "\n";
	}
}
close $INPUT;	
close $OUTPUT;
close $OUTPUT_failed;

sub usage(){
print<<EOF;

Filter for variants that are caused by mismatches between the first 6 bases of reads, by Robert Piskol (piskol\@stanford.edu) 
							   & Gokul Ramaswami (gokulr\@stanford.edu) 07/25/2013

This program takes a variant file and an indexed bam file. It removes variants that are only supported
by variants within the first 6 positions of mapped reads

usage: $0 -infile FILE -outfile FILE -bamfile FILE [-minbasequal N] [-illumina] 
         

Arguments:
-infile FILE	- File containing a list of variants to be filtered
-outfile FILE	- Output file for filtered variants
-bamfile FILE	- File containing mapped short reads that were used for variant calling
-minbasequal N	- Minimum base quality for mismaches to be taken into account (default: 25)
-illumina	- reads in the bam file are in Illumina 1.3+ FASTQ-like format
-help		- Show this help screen                 


EOF

exit 1;
}					

