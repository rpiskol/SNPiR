#########################
## this example demonstrates all filtering steps that can be applied to 
## a set of variant calls from RNA-seq experiments to obtain a confident
## set of genomic variants
##
## - all perl scripts can be run without parameters to display their options
## - variants that fail one of the filters are written into files outfile_failed
## - before running the filters make sure that you have changed the location of 
##   bedtools, samtools and BLAT executables in the file config.pm


##########################
# download supplementary files
#
chrom="chr20"
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/$chrom.fa.gz -O - | gunzip > ./example/$chrom.fa
wget http://www.stanford.edu/~gokulr/database/Human_AG_all_hg19.bed -O ./example/Human_AG_all_hg19.bed


##########################
# convert vcf format into custom SNPiR format and filter variants with quality <20

./convertVCF.sh ./example/$chrom.vcf ./example/$chrom.txt 20


##########################
# filter mismatches at read ends
# note: add the -illumina option if your reads are in Illumina 1.3+ quality format

./filter_mismatch_first6bp.pl -infile ./example/$chrom.txt -outfile ./example/$chrom.rmhex.txt -bamfile ./example/$chrom.bam 


##########################
# filter variants in repetitive regions

awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' ./example/$chrom.rmhex.txt | intersectBed -a stdin -b ./example/hg19_rmsk_chr20_overlap_candidates.bed -v | cut -f1,3-7 > ./example/$chrom.rmhex.rmsk.txt


##########################
# filter intronic sites that are within 4bp of splicing junctions
# make sure your gene annotation file is in UCSC text format and sorted by chromosome and transcript start position

./filter_intron_near_splicejuncts.pl -infile ./example/$chrom.rmhex.rmsk.txt -outfile ./example/$chrom.rmhex.rmsk.rmintron.txt -genefile ./example/hg19_genes_chr20_overlap_candidates.txt


##########################
# filter variants in homopolymers

./filter_homopolymer_nucleotides.pl -infile ./example/$chrom.rmhex.rmsk.rmintron.txt -outfile ./example/$chrom.rmhex.rmsk.rmintron.rmhom.txt -refgenome ./example/$chrom.fa


##########################
# filter variants that were caused by mismapped reads
# this may take a while depending on the number of variants to screen and the size of the reference genome
# note: add the -illumina option if your reads are in Illumina 1.3+ quality format

./BLAT_candidates.pl -infile ./example/$chrom.rmhex.rmsk.rmintron.rmhom.txt -outfile ./example/$chrom.rmhex.rmsk.rmintron.rmhom.rmblat.txt -bamfile ./example/$chrom.bam -refgenome ./example/$chrom.fa


##########################
# remove known RNA editing sites

awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' ./example/$chrom.rmhex.rmsk.rmintron.rmhom.rmblat.txt | intersectBed -a stdin -b ./example/Human_AG_all_hg19.bed -v > ./example/$chrom.rmhex.rmsk.rmintron.rmhom.rmblat.rmedit.bed


##########################
# the final list of variants can be displayed as a track in the igv browser or annotated using Annovar.
# note that the ts/tv ratio of the final set is ~3 and thus larger than the genome average. This agrees with previous findings that ts/tv in exons is 3.0-3.3 (DePristo et al. 2011)   
  