#!/bin/bash
set -e
set -x

BEDTOOLS=bedtools

#NCBI names:
NCBI_NAMES=NCBI_Names.txt
if [ ! -e GCF_003339765.1_Mmul_10_assembly_report.txt ];then
	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_assembly_report.txt
fi
grep -v '#' GCF_003339765.1_Mmul_10_assembly_report.txt | awk '$2 == "assembled-molecule"' | awk -v OFS='\t' ' { print $7, $3 } ' > $NCBI_NAMES
grep -v '#' GCF_003339765.1_Mmul_10_assembly_report.txt | grep -v 'assembled-molecule' | awk -v OFS='\t' ' { print $7, $5 } ' >> $NCBI_NAMES

NCBI=NCBI.Mmul_10.103.gtf.gz
NCBI_TRANSLATED=NCBI.Mmul_10.103.translated.gtf
NCBI19=NCBI.Mmul_10.103.translated.chr19.gtf
# https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Macaca_mulatta/103/
# Version 103
# 21,121 coding, 
if [ ! -e $NCBI ];then
	wget -O $NCBI https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9544/103/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.gtf.gz
fi

if [ ! -e $NCBI_TRANSLATED ];then
	python3 translateNcbiNames.py $NCBI $NCBI_TRANSLATED $NCBI_NAMES
fi

if [ ! -e $NCBI19 ];then
	zcat $NCBI_TRANSLATED | awk ' $1 == "19" ' > $NCBI19
fi

ENSEMBL=Ensembl.Mmul_10.100.gtf.gz
ENSEMBL_UNZIP=Ensembl.Mmul_10.100.gtf
ENSEMBL19=Ensembl.Mmul_10.100.chr19.gtf
# https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Macaca_mulatta/103/
# Version 103
# 21,748 coding, 64,191 transcripts
if [ ! -e $ENSEMBL ];then
	wget -O $ENSEMBL ftp://ftp.ensembl.org/pub/release-100/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.100.gtf.gz
fi

if ! -e $ENSEMBL_UNZIP ];then
	gunzip $ENSEMBL_UNZIP
fi

if [ ! -e $ENSEMBL19 ];then
	zcat $ENSEMBL | awk ' $1 == "19" ' > $ENSEMBL19
fi

#See: https://ftp.ncbi.nlm.nih.gov/gene/DATA/
GENE2ENSEMBL_MM=gene2ensembl.mmul.txt
if [ ! -e $GENE2ENSEMBL_MM ];then
	wget -O gene2ensembl.gz https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz
	zcat gene2ensembl.gz | head -n 1 > $GENE2ENSEMBL_MM
	zcat gene2ensembl.gz | awk '$1 == 9544 ' >> $GENE2ENSEMBL_MM
fi

limitToType(){
	TYPE=$1
	INPUT=$2
	
	OUT=`basename $INPUT ".gtf"`".${TYPE}.gtf"
	
	cat $INPUT | awk -v T=$TYPE ' $3 == T ' > $OUT
	
	echo $OUT
}

#First find reciprocal intersects for genes:
NCBI_GENES=$(limitToType 'gene' $NCBI_TRANSLATED)
ENSEMBL_GENES=$(limitToType 'gene' $$ENSEMBL_UNZIP)
INTERSECT=GeneIntersectBedtools.txt
$BEDTOOLS intersect -f .9 -r -s -sorted -a $NCBI_GENES -b $ENSEMBL_GENES -wa -wb > $INTERSECT

python3 mergeGtf.py $NCBI_TRANSLATED $ENSEMBL_UNZIP $INTERSECT $GENE2ENSEMBL_MM Ensembl.v100.NCBI.v103 ./output/