#!/usr/bin/bash

SAMPLES=$1
LINEAGE=$2
SET=$3
REF=$4
oriANALYTYPE=$5

baseFOLDER="${HOME}/Proyectos/${LINEAGE}"
#fastqcFOLDER="${HOME}/Programs/FastQC"
#fastpFOLDER="${HOME}/anaconda3/bin"
#multiqcFOLDER="${HOME}/anaconda3/bin"
#bowtie2FOLDER="${HOME}/Programs/bowtie2-2.5.4"
#javaFOLDER="${HOME}/Programs/jre1.8.0_331/bin"
#samtoolsFOLDER="${HOME}/Programs/samtools-0.1.16"
#picardFOLDER="${HOME}/Programs/picard-tools-1.140" 
#GATKfolder="${HOME}/Programs/GenomeAnalysisTK-3.4-46"

    mkdir -p ${baseFOLDER}/FastP_TM_fastq

for i in ${baseFOLDER}/Muestras/*_fastq.gz
do 
    echo "Descomprimiendo muestra $i"
    gunzip "$i"
done 

while read SAMPLE_ID
do
    
    # Ejecutamos fastp
    echo "Procesando muestra $SAMPLE_ID"  
    fastp \
    -i "${baseFOLDER}/Muestras/${SAMPLE_ID}_R1_fastq" -I "${baseFOLDER}/Muestras/${SAMPLE_ID}_R2_fastq" \
    -o "${baseFOLDER}/FastP_TM_fastq/${SAMPLE_ID}_L001_R1_001.trimmed.fastq" \
    -O "${baseFOLDER}/FastP_TM_fastq/${SAMPLE_ID}_L001_R2_001.trimmed.fastq" \
    -h "${baseFOLDER}/FastP_TM_fastq/${SAMPLE_ID}.html" \
    -j "${baseFOLDER}/FastP_TM_fastq/${SAMPLE_ID}.json" \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --qualified_quality_phred 25 \
    --cut_front \
    --cut_tail 
    
done<${SAMPLES}

for i in ${baseFOLDER}/Muestras/*_fastq
do 
    echo "Comprimiendo muestra $i"
    gzip "$i"
done 

    mkdir ${baseFOLDER}/multiqc_trimmed

    multiqc --outdir ${baseFOLDER}/multiqc_trimmed ${baseFOLDER}/FastP_TM_fastq/ 

#sh ${HOME}/Scripts/DANIELA/DO_VARIANTCALLING.sh ${SAMPLES} ${LINEAGE} ${SET} ${REF} ${oriANALYTYPE}
