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
# perl script that screens if the variant is located in a homopolymer
#inputs: input file of edits, output file name
my ($INPUTFILE,$OUTPUTFILE,$REFERENCEGENOME,$HELP);

parse_command_line();

sub parse_command_line {
	my $help;
	
	usage() if (scalar @ARGV==0);

	&GetOptions(
	    "infile=s" => \$INPUTFILE,
	    "outfile=s" => \$OUTPUTFILE,
	    "refgenome=s" => \$REFERENCEGENOME,
	    "help" => \$HELP
	    );
	
	usage() if ($help);
	   
	if(!$INPUTFILE || !$OUTPUTFILE ||  !$REFERENCEGENOME){
	    usage();
	}

}


my ($leftbuffer, $rightbuffer) = (4,4); #sequence buffer length on both sides

my $tmpfile = join '', $INPUTFILE, '.nuctmp';
my $tmpfafile = join '', $INPUTFILE, '.tmpnuctmp';


open (SITES, "<", $INPUTFILE ) or die "error opening inputfile: $!\n"; #Read INPUTFILE
open (OUTPUT, ">", $OUTPUTFILE) or die "error opening outputfile: $!\n";
open (OUTPUT_failed, ">", $OUTPUTFILE."_failed") or die "error opening outputfile: $!\n";


while(<SITES>) { #Read INPUTFILE
	chomp;
	my $line = $_;
	my @fields = split;
	my ($chromname, $position, $editnuc) = ($fields[0], $fields[1], $fields[4]); #### CAN CHANGE $editnuc to $fields[X] for your input file
	my $homopol = 0;
	
	open (TMP, ">", $tmpfile);
	my ($startpos, $endpos) = ($position-$leftbuffer,$position+$rightbuffer+1);
	print TMP "$chromname\t$startpos\t$endpos\n";
	#print "$chromname\t$startpos\t$endpos\n";
	close TMP;

	system("$FASTAFROMBED -fi $REFERENCEGENOME -bed $tmpfile -fo $tmpfafile"); #get sequence
	#print ("/srv/gs1/projects/li/shared-data/software/BEDTools-Version-2.12.0/bin/fastaFromBed -fi /srv/gs1/projects/li/shared-data/reference-genomes/$refspec/softmasked/chr/$chromname.fa -bed $tmpfile -fo $tmpfafile"); #get sequence
	my $seq = '';
	open (TMPFAFILE, "<", $tmpfafile);
	while(<TMPFAFILE>) {
		chomp;
		next if ($_ =~ m/\>/);
		$seq = join '', $seq, $_; #get sequence flanking variant
	} 
	close TMPFAFILE;
	system("rm $tmpfafile");
	system("rm $tmpfile");
	my @splitseq = split(//,$seq);
	$editnuc = uc $editnuc;
	for (my $i = 0; $i < @splitseq; $i++) {
		$splitseq[$i] = uc $splitseq[$i];
	}

	$homopol = 1 if ($editnuc =~ $splitseq[0] && $editnuc =~ $splitseq[1] && $editnuc =~ $splitseq[2] && $editnuc =~ $splitseq[3] );
	$homopol = 1 if ($editnuc =~ $splitseq[1] && $editnuc =~ $splitseq[2] && $editnuc =~ $splitseq[3] && $editnuc =~ $splitseq[5] );
	$homopol = 1 if ($editnuc =~ $splitseq[2] && $editnuc =~ $splitseq[3] && $editnuc =~ $splitseq[5] && $editnuc =~ $splitseq[6] );
	$homopol = 1 if ($editnuc =~ $splitseq[3] && $editnuc =~ $splitseq[5] && $editnuc =~ $splitseq[6] && $editnuc =~ $splitseq[7] );
	$homopol = 1 if ($editnuc =~ $splitseq[5] && $editnuc =~ $splitseq[6] && $editnuc =~ $splitseq[7] && $editnuc =~ $splitseq[8] );

	if(!$homopol){
		print OUTPUT "$line\n"; #check if homopolymer
	} 
	else{
		print OUTPUT_failed "$line\n";
	}
	#print "homopol\n" if $homopol;
}

close SITES;
close OUTPUT;
close OUTPUT_failed;

sub usage(){
print<<EOF;

Homopolymer filter, by Robert Piskol (piskol\@stanford.edu) & Gokul Ramaswami (gokulr\@stanford.edu) 07/25/2013

This program takes a variant file, and a reference genome and removes variants that are located
whtin homopolymers. 

usage: $0 -infile FILE -outfile FILE -refgenome FILE 
         

Arguments:
-infile FILE	- File containing a list of variants to be filtered
-outfile FILE	- Output file for filtered variants
-refgenome FILE	- File in FASTA format containing the reference genome of the species
-help		- Show this help screen                 


EOF

exit 1;
}