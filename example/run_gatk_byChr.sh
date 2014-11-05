. ~/projects/shared-data/software/.source
mkdir -p /srv/gsfs0/scratch/piskol/ENCODE/GM12878/Gm12878CellPapAlnMerged/STAR/variants_gatk_alt/
cd /srv/gsfs0/scratch/piskol/ENCODE/GM12878/Gm12878CellPapAlnMerged/STAR/variants_gatk_alt/

ln -s /srv/gs1/projects/li/shared-data/reference-genomes/hg19/softmasked/hg19_softmasked.dict
ln -s /srv/gs1/projects/li/shared-data/reference-genomes/hg19/softmasked/hg19_softmasked.fa
ln -s /srv/gs1/projects/li/shared-data/reference-indices/hg19/samtools/hg19_softmasked.fa.fai

############################
# create regions for realignment
if [ ! -e Aligned.out.sorted.bam.$1.intervals ]
then
/usr/java/latest/bin/java -jar /srv/gs1/projects/li/shared-data/software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R hg19_softmasked.fa \
    -U ALLOW_N_CIGAR_READS \
    -I /srv/gsfs0/scratch/piskol/ENCODE/GM12878/Gm12878CellPapAlnMerged/STAR/Aligned.out.sorted.bam \
    -o Aligned.out.sorted.bam.$1.intervals \
    -L $1 
    #-known gold_indels.vcf \ 

fi

############################
# realign regions
if [ ! -e /srv/gsfs0/scratch/piskol/ENCODE/GM12878/Gm12878CellPapAlnMerged/STAR/Aligned.out.sorted.realign.$1.bam ]
then
/usr/java/latest/bin/java -jar /srv/gs1/projects/li/shared-data/software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R hg19_softmasked.fa \
    -I /srv/gsfs0/scratch/piskol/ENCODE/GM12878/Gm12878CellPapAlnMerged/STAR/Aligned.out.sorted.bam \
    -targetIntervals Aligned.out.sorted.bam.$1.intervals \
    -U ALLOW_N_CIGAR_READS \
    -o /srv/gsfs0/scratch/piskol/ENCODE/GM12878/Gm12878CellPapAlnMerged/STAR/Aligned.out.sorted.realign.$1.bam \
    -known /srv/gs1/projects/li/shared-data/annotations/hg19_variants/dbsnp_132_hg19.vcf \
    -L $1

fi

###########################
# determine covariates
#######
#PASS 1
if [ ! -e Aligned.out.sorted.bam.$1.table ]
then
/usr/java/latest/bin/java -jar /srv/gs1/projects/li/shared-data/software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R hg19_softmasked.fa \
    -I /srv/gsfs0/scratch/piskol/ENCODE/GM12878/Gm12878CellPapAlnMerged/STAR/Aligned.out.sorted.realign.$1.bam \
    -o Aligned.out.sorted.bam.$1.table \
    -U ALLOW_N_CIGAR_READS \
    -knownSites /srv/gs1/projects/li/shared-data/annotations/hg19_variants/dbsnp_132_hg19.vcf \
    -L $1
    #-knownSites dbsnp.vcf
    #-knownSites gold_indels.vcf \ 

fi
#######
#PASS 2
if [ ! -e Aligned.out.sorted.bam.$1.table.post ]
then
/usr/java/latest/bin/java -jar /srv/gs1/projects/li/shared-data/software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R hg19_softmasked.fa \
    -I /srv/gsfs0/scratch/piskol/ENCODE/GM12878/Gm12878CellPapAlnMerged/STAR/Aligned.out.sorted.realign.$1.bam \
    -BQSR Aligned.out.sorted.bam.$1.table \
    -knownSites /srv/gs1/projects/li/shared-data/annotations/hg19_variants/dbsnp_132_hg19.vcf \
    -U ALLOW_N_CIGAR_READS \
    -o Aligned.out.sorted.bam.$1.table.post \
    -L $1
    #-knownSites dbsnp.vcf 
    #-knownSites gold_indels.vcf \
fi

###########################
# recalibrate sequences
if [ ! -e /srv/gsfs0/scratch/piskol/ENCODE/GM12878/Gm12878CellPapAlnMerged/STAR/Aligned.out.sorted.realign.recal.$1.bam ]
then
/usr/java/latest/bin/java -jar /srv/gs1/projects/li/shared-data/software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R hg19_softmasked.fa \
    -I /srv/gsfs0/scratch/piskol/ENCODE/GM12878/Gm12878CellPapAlnMerged/STAR/Aligned.out.sorted.realign.$1.bam \
    -BQSR Aligned.out.sorted.bam.$1.table \
    -U ALLOW_N_CIGAR_READS \
    -o /srv/gsfs0/scratch/piskol/ENCODE/GM12878/Gm12878CellPapAlnMerged/STAR/Aligned.out.sorted.realign.recal.$1.bam \
    -L $1
fi

###########################
# call variants using the unified genotyper
if [ ! -e Aligned.out.sorted.realign.recal.snps.raw.$1.vcf ]
then
/usr/java/latest/bin/java -jar /srv/gs1/projects/li/shared-data/software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar \
    -l INFO \
    -R hg19_softmasked.fa \
    -T UnifiedGenotyper \
    -I /srv/gsfs0/scratch/piskol/ENCODE/GM12878/Gm12878CellPapAlnMerged/STAR/Aligned.out.sorted.realign.recal.$1.bam  \
    -o Aligned.out.sorted.realign.recal.snps.raw.$1.vcf \
    -stand_call_conf 0.0 \
    -stand_emit_conf 0.0 \
    -U ALLOW_N_CIGAR_READS \
    -out_mode EMIT_VARIANTS_ONLY \
    -L $1
fi

###########################
# call variants using the haplotype caller
if [ ! -e Aligned.out.sorted.realign.recal.snps.haplotype.raw.$1.vcf ]
then
/usr/java/latest/bin/java -jar /srv/gs1/projects/li/shared-data/software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R hg19_softmasked.fa \
    -I /srv/gsfs0/scratch/piskol/ENCODE/GM12878/Gm12878CellPapAlnMerged/STAR/Aligned.out.sorted.realign.recal.$1.bam \
    -recoverDanglingHeads \
    -dontUseSoftClippedBases \
    -stand_call_conf 20.0 \
    -stand_emit_conf 20.0 \
    -U ALLOW_N_CIGAR_READS \
    -o Aligned.out.sorted.realign.recal.snps.haplotype.raw.$1.vcf \
    -L $1
fi