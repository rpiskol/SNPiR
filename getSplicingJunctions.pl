#!/usr/bin/env perl


$| = 1;
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Spec;
use File::Basename;
use Data::Dumper;
use FindBin;
#use lib qw(../../../ ../../ /srv/gs1/projects/li/shared-data/software/BioPerl-1.6.1);
#use lib "$FindBin::Bin/../";
#use SNPiR_bak::config;
use Bio::DB::Fasta;

################################################################################
################################################################################
# $Revision: $
# Authors: Robert Piskol ( piskol@stanford.edu )
# Last modification $Author: piskol $
# perl script that creates sequences of a given length around known splicing junctions
# inputs: genome file, gene annotation file
my ($REFGENOME,$GENEFILE,$OUTPUTFILE,$HELP);
my $WINDOWSIZE = 95;

parse_command_line();

sub parse_command_line {
	my $help;
	
	usage() if (scalar @ARGV==0);

	&GetOptions(
	    "refgenome=s" => \$REFGENOME,
	    "outfile=s" => \$OUTPUTFILE,
	    "genefile=s" => \$GENEFILE,
	    "windowsize=s" => \$WINDOWSIZE,
	    "help" => \$HELP
	    );
	
	usage() if ($help);
	   
	if(!$REFGENOME || !$GENEFILE || !$OUTPUTFILE){
	    usage();
	}
}

print STDERR "using window size: $WINDOWSIZE\n";
my $JUNCTIONS = getJunctions($GENEFILE);


############
# write splicing junctions to file
print STDERR "writing splicing junctions to file...\n";
open (OUT, ">", $OUTPUTFILE) or die "error opening output file: $!\n";
foreach my $key(sort(keys %$JUNCTIONS)){
	print OUT ">$key\n";
	print OUT $JUNCTIONS->{$key}."\n";
}
close OUT;


sub getJunctions{
	my $annotfile = shift;
	my $junctions;
	
	############
	# read genome file
	my $seqs;
	print STDERR "reading genome file...\n";
	my $seqDb =  Bio::DB::Fasta->new($REFGENOME);
	#print_dumper($seqDb);
	
	############
	# read annotation file
	print STDERR "reading annotation file...\n";
	my $annot = readAnnotFile($annotfile);
	#print_dumper($annot);
	#exit;
	
	############
	# go through all spliceforms
	print STDERR "generating splicing junctions...\n";
	foreach my $key(sort(keys %$annot)){
		
		#print $key."\n";
		my @exonStarts = @{$annot->{$key}->{exonStarts}};
		my @exonEnds = @{$annot->{$key}->{exonEnds}};
		my ($chromosome) = $key =~ /(.*):.*/;
		
		next if($#exonStarts == 0);
		#print $chromosome."\n";
		#print_dumper($annot->{$key});
		for(my $i = 1; $i <= $#exonStarts; $i++){
			
			#print "#######\n#\n";
			my ($lSeq,$lSeqPos) = getLseq(\@exonStarts,\@exonEnds, $i-1,$seqDb,$chromosome);
			my ($rSeq,$rSeqPos) = getRseq(\@exonStarts,\@exonEnds, $i,$seqDb,$chromosome);

			my $junctionId =  $chromosome."-".$lSeqPos."-".$rSeqPos;

			print STDERR "$junctionId\r";
			$junctions->{$junctionId} = $lSeq.$rSeq;
		}
		
	}
	return($junctions);
	#print_dumper($junctions);
	#print "\n";
}

sub getLseq{
	my ($exonStarts, $exonEnds, $i,$seqDb,$chromosome) = @_;
	
	my $remaining = $WINDOWSIZE;
	my $lSeq = "";
	my @lSeqPosArr = ();
	
	while($remaining>0){
		if($i>=0){
			
			my $winStart = $exonEnds->[$i]-$remaining < $exonStarts->[$i] ? $exonStarts->[$i] : $exonEnds->[$i]-$remaining;
			#$lSeq = uc($seqDb->get_Seq_by_id($chromosome)->subseq($winStart+1=>$exonEnds->[$i])).$lSeq;
			$lSeq = uc(getSeqFromGenome($chromosome,$winStart+1,$exonEnds->[$i],$seqDb)).$lSeq;
			unshift @lSeqPosArr, ($winStart+1)."-".$exonEnds->[$i];
			$remaining -= ($exonEnds->[$i] - $winStart);
			$i -= 1;
		}
		else{
			return ($lSeq,join("-",@lSeqPosArr));
		}
	}
	return ($lSeq,join("-",@lSeqPosArr));	
}

sub getRseq{
	my ($exonStarts, $exonEnds, $i,$seqDb,$chromosome) = @_;
	
	my $remaining = $WINDOWSIZE;
	my $rSeq = "";
	my @rSeqPosArr = ();
	
	while($remaining>0){
		if($i<=scalar(@{$exonStarts})-1){
			my $winEnd = $exonStarts->[$i]+$remaining > $exonEnds->[$i] ? $exonEnds->[$i] : $exonStarts->[$i]+$remaining;
			#$rSeq = $rSeq.uc($seqDb->get_Seq_by_id($chromosome)->subseq($exonStarts->[$i]+1=>$winEnd));
			$rSeq = $rSeq.uc(getSeqFromGenome($chromosome,$exonStarts->[$i]+1,$winEnd,$seqDb));
			push @rSeqPosArr, ($exonStarts->[$i]+1)."-".$winEnd;
			$remaining -= ($winEnd - $exonStarts->[$i]);
			$i += 1;
		}
		else{
			return ($rSeq,join("-",@rSeqPosArr));			
		}
	}
	return ($rSeq,join("-",@rSeqPosArr));		
}

sub getSeqFromGenome{
	my ($chromosome, $start,$end, $seqDb) = @_;
	
	return $seqDb->get_Seq_by_id($chromosome)->subseq($start=>$end);

}

##########################
# get seq from genome using fastaFromBed in the bedtools package
#sub getSeqFromGenome2{
#	my ($chromosome, $start,$end, $seqDb) = @_;
#	
#	open (TMP, ">$$.tmpfile");
#	print TMP "$chromosome\t".($start-1)."\t$end\n";
#	close TMP;
#	
#	system("$FASTAFROMBED -fi $REFGENOME -bed $$.tmpfile -fo $$.tmpfafile"); #get sequence
#	#print ("/srv/gs1/projects/li/shared-data/software/BEDTools-Version-2.12.0/bin/fastaFromBed -fi /srv/gs1/projects/li/shared-data/reference-genomes/$refspec/softmasked/chr/$chromname.fa -bed $tmpfile -fo $tmpfafile"); #get sequence
#	my $seq = '';
#	open (TMPFAFILE, "<$$.tmpfafile");
#	while(<TMPFAFILE>) {
#		chomp;
#		next if ($_ =~ m/\>/);
#		$seq = join '', $seq, $_; #get sequence flanking variant
#	} 
#	close TMPFAFILE;
#	return $seq;
#	#system("rm $$.tmpfafile");
#	#system("rm $$.tmpfile");
#}

sub readAnnotFile{
	my $annotfile = shift;
	my $annot;
	print STDERR $annotfile."\n";
	open(IN, "<".$annotfile);
	my @filecnt = <IN>;
	shift @filecnt;
	#print @filecnt; 
	foreach(@filecnt){
		my $line = $_;
		#print $line."\n";
		my @cnt = split("\t", $line);
		#if($cnt[2] !~ /chr.*_.*/){
			$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5].$cnt[1]}->{exonStarts} = [split(",",$cnt[9])];
			$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5].$cnt[1]}->{exonEnds} = [split(",",$cnt[10])];
			$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5].$cnt[1]}->{gene} = $cnt[1];
			#$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5]}->{exonStarts} = [split(",",$cnt[9])];
			#$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5]}->{exonEnds} = [split(",",$cnt[10])];
			#$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5]}->{gene} = $cnt[1];
		#}
	}
	
	#print_dumper($annot);
	return $annot;
}

sub print_dumper
{
	my $struct = shift;
	
	my $dumper = Data::Dumper -> new([$struct]);
	my $dump  = $dumper->Dump();
	print $dump;
}


sub usage(){
print<<EOF;

Create splicing junctions, by Robert Piskol (me\@robertpiskol.de) 11/07/2014

This program takes a reference genome, gene annotation file in UCSC txt format and a window size to create sequences up/downstream of known splicing junctions.
If an up/downstream exon is shorter than the desired sequence window, the sequence is extended into the next exon.

usage: $0 -refgenome FILE -genefile FILE -windowsize INT -outfile FILE 
         

Arguments:
-windowsize INT - size of the exonic region up/downstream of splicing junction that is included (default: 95)
-refgenome FILE	- refernce genome in fasta format
-genefile FILE	- File in UCSC txt format (sorted by chomosome and position - i.e. 'sort -k3,3 -5,5n')
-outfile FILE   - name of the file that will contain the splicing junctions
-help		- Show this help screen                 


EOF

exit 1;
}

