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
# perl script that runs blat for all aligned reads for each editing sites
my ($INPUTFILE,$OUTPUTFILE,$BAMFILE,$REFERENCEGENOME,$HELP,$ILLUMINAQUALS);
my $MINBASEQUAL = 25;
my $MINMISMATCH = 1;
my $SCORELIMIT = 0.95;
my $QUALOFFSET = 33;

parse_command_line();

sub parse_command_line {
	my $help;
	
	usage() if (scalar @ARGV==0);

	&GetOptions(
	    "infile=s" => \$INPUTFILE,
	    "outfile=s" => \$OUTPUTFILE,
	    "bamfile=s" => \$BAMFILE,
	    "refgenome=s" => \$REFERENCEGENOME,
	    "minbasequal=s" => \$MINBASEQUAL,
	    "minmismatch=s" => \$MINMISMATCH,
	    "scorelimit=s" => \$SCORELIMIT,
	    "illumina" => \$ILLUMINAQUALS,
	    "help" => \$HELP
	    );
	
	usage() if ($help);
	   
	if(!$INPUTFILE || !$OUTPUTFILE || !$BAMFILE || !$REFERENCEGENOME){
	    usage();
	}
	if($ILLUMINAQUALS){
		$QUALOFFSET = 64;
	}
}


open (my $INPUT , "<", $INPUTFILE) or die "error opening inputfile: $!\n";
open (my $OUTPUT, ">", $OUTPUTFILE) or die "error opening outputfile: $!\n";
open (my $OUTPUT_failed, ">", $OUTPUTFILE."_failed") or die "error opening outputfile: $!\n";

my $fafile = join '', $INPUTFILE, '.fatmp';
my $pslfile = join '', $INPUTFILE, '.psltmp';


#if (! -e $pslfile){
	#if(! -e $fafile){
		print STDERR "preparing blat input...\n";
	    open (FAFILE, ">", $fafile);
		while (<$INPUT>) {
			chomp;
			my $inputline = $_;
			my @fields = split;
			my $TEMPNAME = join '', $OUTPUTFILE,'_tmp';
			my ($chr, $position) = ($fields[0], $fields[1]);
			my $bamposition = join '', $chr,':',$position,'-',$position;
			system("$SAMTOOLSEXE view $BAMFILE $bamposition > $TEMPNAME");
			my $editnuc = $fields[4];
			my $newmismatch = 0;
			my $mismatchreadcount = 0;
		
			open(my $TMPFILE, "<", $TEMPNAME);
			while (<$TMPFILE>) {
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
							$base_readpos = 1 if ($currentpos == $position && $sequencebases[$readpos-1] eq $editnuc && ord($qualscores[$readpos-1]) >= $MINBASEQUAL+$QUALOFFSET);
							$currentpos++;
							$readpos++;	
						}	
					}
				}
				if ($base_readpos) {
					print FAFILE ">$chr\-$position\-$mismatchreadcount\n$sequence\n";
					$mismatchreadcount++;
				}
			}
			close $TMPFILE;
			system("rm $TEMPNAME");
		}
	#}
	close $INPUT;

	print STDERR "blatting reads...\n";
	system("$BLATEXE -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 -noHead $REFERENCEGENOME $fafile $pslfile");
	system("rm $fafile");
#}


#print "reading pslfile\n";
print STDERR "processing reads...\n";
open(PSL, "<", $pslfile );
my %pslhash;
while(<PSL>) {
	chomp;
	my @pslfields = split;
	my $name = $pslfields[9];
	my $blatscore = join '@',$pslfields[0],$pslfields[13],$pslfields[17],$pslfields[18],$pslfields[20];
	if ($pslhash{$name}) {
		$pslhash{$name} = join '-', $pslhash{$name}, $blatscore;
	} elsif ( !$pslhash{$name}) {
		$pslhash{$name} = $blatscore; 
	}
}
close PSL;	

my %sitehash;
my %discardhash;
foreach my $pslkey (keys %pslhash) {
	my @splitkey = split(/\-/, $pslkey);
	my $site = join '_', $splitkey[0],$splitkey[1];
	my @psllines = split(/\-/, $pslhash{$pslkey});
	my $largestscore = 0;
	my $largestscoreline = $psllines[0];
	my @scorearray;
	foreach my $scoreline (@psllines) {
		my @scoresarray = split(/\@/,$scoreline);
		my $linescore = $scoresarray[0];
		push(@scorearray,$linescore);
		if ($linescore > $largestscore) {
			$largestscoreline = $scoreline;
			$largestscore = $linescore;
		}
	}
	@scorearray = sort {$b <=> $a} @scorearray;
	$scorearray[1] = 0 unless ($scorearray[1]);
	my @splitlargestline = split(/\@/,$largestscoreline);
	my $overlapfound = 0;
	if ($splitlargestline[1] eq $splitkey[0] && $scorearray[1] < ($scorearray[0]*$SCORELIMIT)) {
		my ($numblocks, $blocksizes, $blockstarts) = ($splitlargestline[2],$splitlargestline[3],$splitlargestline[4]);
		my @blocksizes = split(/\,/,$blocksizes);
		my @blockstarts = split(/\,/,$blockstarts);
		for (my $i = 0; $i < $numblocks; $i++) {
			my $startpos = $blockstarts[$i]+1;
			my $endpos = $blockstarts[$i] + $blocksizes[$i];
			$overlapfound = 1 if ($splitkey[1] >= $startpos && $splitkey[1] <= $endpos);
		}
		if ($sitehash{$site} && $overlapfound) {
			$sitehash{$site}++;
		} elsif (!$sitehash{$site} && $overlapfound) {
			$sitehash{$site} = 1;
		}
	}
	unless ($overlapfound) {
		if ($discardhash{$site}) {
			$discardhash{$site}++;
		} elsif (!$discardhash{$site}) {
			$discardhash{$site} = 1;
		}
	}
}

open (SITES2, "<", $INPUTFILE ) or die "error opening inputfile: $!\n";
while(<SITES2>) {
	chomp;
	my @fields = split;
	my $inputline = $_;
	my ($cov,$oldalter) = split(/\,/,$fields[2]);
	my $name = join '', $fields[0],'_',$fields[1];
	if ($sitehash{$name}) {
		my $newalter = $sitehash{$name};
		my $discardnum;
		if ($discardhash{$name}) {
			$discardnum = $discardhash{$name};
		} else {
			$discardnum = 0;
		}
		my $newcov = $cov - ($oldalter - $newalter);
		my $neweditfreq = sprintf("%.3f", $newalter/$newcov);
		if ($newalter >= $MINMISMATCH && $newalter > $discardnum){
			print $OUTPUT "$fields[0]\t$fields[1]\t$newcov,$newalter\t$fields[3]\t$fields[4]\t$neweditfreq\n" ;
		}
		else{
			print $OUTPUT_failed join("\t",@fields)."\tfailed_freq #mismatches: $newalter #minimumMismatchesNecessary: $MINMISMATCH #discarded mismatch reads: $discardnum\n";
		}
	}
	else{
		print $OUTPUT_failed join("\t",@fields)."\tfailed_totalcover\n";
	}
}
close SITES2;
close $OUTPUT;
close $OUTPUT_failed;
system("rm $pslfile");



sub usage(){
print<<EOF;

BLAT filter, by Robert Piskol (piskol\@stanford.edu) & Gokul Ramaswami (gokulr\@stanford.edu) 07/25/2013

This program takes a variant file, bamfile with mapped reads (that were used for variant calling)
and reference genome and filters all variants that were erroneously called from potentially mismapped reads.

usage: $0 -infile FILE -outfile FILE -bamfile FILE -refgenome FILE 
		 [-minbasequal N] [-minmismatch N] [-scorelimit N] [-illumina]
         

Arguments:
-infile FILE	- File containing a list of variants to be filtered
-outfile FILE	- Output file for filtered variants
-bamfile FILE	- File containing mapped short reads that were used for variant calling
-refgenome FILE	- File in FASTA format containing the reference genome of the species
-minbasequal N	- minimum basequality of a mismatch to be taken into account (default: 25)
-minmismatch N	- minimum number of mismatches that are supported by correctly mapped reads (default: 1)
-scorelimit N	- fraction of a max read score at which other mapping locations of the same read are
                  considered duplicates (default: 0.95)
                  (e.g. if the scond best mapping of a read has a score>0.95*bestMapping,
                   the read mapping is considered unsafe and the read removed)
-illumina	- reads in the bam file are in Illumina 1.3+ FASTQ-like format
-help		- Show this help screen                 


EOF

exit 1;
}
