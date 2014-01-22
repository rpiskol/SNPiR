#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../";
use snpir::config;

################################################################################
################################################################################
# $Revision: $
# Authors: Robert Piskol ( piskol@stanford.edu ), Gokul Ramaswami ( gokulr@stanford.edu )
# Last modification $Author: piskol $
# script used to filter out variants called due to intronic alignment of RNA-SEQ reads
my ($INPUTFILE,$OUTPUTFILE,$GENEFILE,$HELP);
my $SPLICEDIST = 4;

parse_command_line();

sub parse_command_line {
	my $help;
	
	usage() if (scalar @ARGV==0);

	&GetOptions(
	    "infile=s" => \$INPUTFILE,
	    "outfile=s" => \$OUTPUTFILE,
	    "genefile=s" => \$GENEFILE,
	    "splicedist=s" => \$SPLICEDIST,
	    "help" => \$HELP
	    );
	
	usage() if ($help);
	   
	if(!$INPUTFILE || !$OUTPUTFILE || !$GENEFILE){
	    usage();
	}
}


my %genehash;

open (GENEFILE, $GENEFILE) or die "error opening gene annotation file: $!\n";
open (my $OUTPUT, ">", $OUTPUTFILE) or die "error opening output file: $!\n";;
open (my $OUTPUT_failed, ">", $OUTPUTFILE."_failed") or die "error opening output file: $!\n";

print STDERR "reading genefile...\n";
while (<GENEFILE>) {
	chomp;
	my @fields = split;
	my $chr = $fields[2];
	push(@{$genehash{$chr}},$_);
}
close GENEFILE;

print STDERR "removing variants in $SPLICEDIST bp of splicing junctions...\n";
open (INPUTFILE, $INPUTFILE) or die "error opening input variant list: $!\n"; #open file of candidate editing sites
while (<INPUTFILE>) {
	chomp;
	my $line = $_;
	my @fields = split(/\t/);
	my ($chrom, $position, $found, $chromfound, $splice) = ($fields[0], $fields[1], 0, 0, 0);
	my ($gene, $strand, $gene2) = ('','','');
	foreach my $geneline (@{$genehash{$chrom}}) { #for each variant, loop through gene annotation file and check if variant is in exon and if variant is in intronic sequence near an intron/exon boundary
		my @fieldsref = split(/\t/,$geneline);
		my ($chromref, $txstart, $txend) = ($fieldsref[2], $fieldsref[4], $fieldsref[5]);
		$chromfound = 1;
		if (int($txstart) <= int($position) && int($txend) >= int($position)) {
			$found = 1;
			$strand = join '',$strand ,',' ,$fieldsref[3];
			$gene = join '',$gene ,',' ,$fieldsref[1];
			$gene2 = join '',$gene2 ,',' ,$fieldsref[12];
			my $exoncount = int($fieldsref[8]);
			my @exonstarts = split(/,/, $fieldsref[9]);
			my @exonends = split(/,/, $fieldsref[10]);
			for (my $i = 0; $i < $exoncount; $i++) { #for each exon, check if variant lies in it
				$splice = 1 if ((int($exonstarts[$i])-$SPLICEDIST < int($position) and int($exonstarts[$i])+1 > $position) or (int($exonends[$i]) < $position and int($exonends[$i])+$SPLICEDIST >= $position)); #make sure variant is not within 4 bp to intronic side of ANY intron/exon boundary		
			}
		} 
	}
	if ($splice == 0) { #if variant is within any exon and not within any intronic near splice boundary, print it out
		$gene =~ s/^.//;
		$strand =~ s/^.//;
		$gene2 =~ s/^.//;
		print $OUTPUT "$line\n"; 
	} 
	else{
		print $OUTPUT_failed $line."\n";
	}
}
close INPUTFILE;	
close $OUTPUT;
close $OUTPUT_failed;

sub usage(){
print<<EOF;

Filter for intronic variants close to splicing junctions, by Robert Piskol (piskol\@stanford.edu) 
							   & Gokul Ramaswami (gokulr\@stanford.edu) 07/25/2013

This program takes a variant file and a gene annotation file in UCSC text format
and filters all variants that are in intronic regions in a distance closer than a user selected value.

usage: $0 -infile FILE -outfile FILE -genefile FILE [-splicedist N] 
         

Arguments:
-infile FILE	- File containing a list of variants to be filtered
-outfile FILE	- Output file for filtered variants
-genefile FILE	- File in UCSC txt format (sorted by chomosome and position - i.e. 'sort -k3,3 -5,5n')
-splicedist N	- Maximum filter distance from splicing junction for variants (default: 4)
-help		- Show this help screen                 


EOF

exit 1;
}					
