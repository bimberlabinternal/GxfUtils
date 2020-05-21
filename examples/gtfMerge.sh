#!/bin/bash
set -e
set -x

BEDTOOLS=bedtools
NCBI_VERSION=103
ENSEMBL_VERSION=100

#Download NCBI/Ensembl nmae translation:
NCBI_NAMES=NCBI_Name_Translation.txt
if [ ! -e GCF_003339765.1_Mmul_10_assembly_report.txt ];then
	wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_assembly_report.txt
fi
grep -v '#' GCF_003339765.1_Mmul_10_assembly_report.txt | awk '$2 == "assembled-molecule"' | awk -v OFS='\t' ' { print $7, $3 } ' > $NCBI_NAMES
grep -v '#' GCF_003339765.1_Mmul_10_assembly_report.txt | grep -v 'assembled-molecule' | awk -v OFS='\t' ' { print $7, $5 } ' >> $NCBI_NAMES

NCBI=NCBI.Mmul_10.${NCBI_VERSION}.gtf.gz
NCBI_TRANSLATED=NCBI.Mmul_10.${NCBI_VERSION}.translated.gtf
if [ ! -e $NCBI ];then
	wget -O $NCBI https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9544/${NCBI_VERSION}/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.gtf.gz
fi

# NCBI encodes their names using RefSeq, not Genbank Accession
if [ ! -e $NCBI_TRANSLATED ];then
	python3 translateNcbiNames.py $NCBI $NCBI_TRANSLATED $NCBI_NAMES
fi

ENSEMBL=Ensembl.Mmul_10.${ENSEMBL_VERSION}.gtf.gz
ENSEMBL_UNZIP=Ensembl.Mmul_10.${ENSEMBL_VERSION}.gtf
if [ ! -e $ENSEMBL_UNZIP ];then
	if [ ! -e $ENSEMBL ];then
		wget -O $ENSEMBL ftp://ftp.ensembl.org/pub/release-${ENSEMBL_VERSION}/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.${ENSEMBL_VERSION}.gtf.gz
	fi

	gunzip $ENSEMBL_UNZIP
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
ENSEMBL_GENES=$(limitToType 'gene' $ENSEMBL_UNZIP)
INTERSECT=GeneIntersectBedtools.txt
$BEDTOOLS intersect -f .9 -r -s -a $NCBI_GENES -b $ENSEMBL_GENES -wa -wb > $INTERSECT
rm $NCBI_GENES
rm $ENSEMBL_GENES

python3 mergeGtf.py $NCBI_TRANSLATED $ENSEMBL_UNZIP $INTERSECT $GENE2ENSEMBL_MM Ensembl.v${ENSEMBL_VERSION}.NCBI.v${NCBI_VERSION} ./output/