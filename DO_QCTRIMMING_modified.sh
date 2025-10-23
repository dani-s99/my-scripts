#!/usr/bin/bash

SAMPLES=$1 #Identificador de las muestras
LINEAGE=$2 #Carpeta del proyecto
SET=$3 #Según si estan mis muestras combinadas o las de Biel
REF=$4 #El identificador de la referencia,si luego hace falta el .fasta o .gff poner solo esto
oriANALYTYPE=$5 #Si se hace el análisis con las muestras trimadas (TR), sin trimar (UT) o ambas (BOTH)

baseFOLDER="${HOME}/Proyectos/${LINEAGE}"


    mkdir ${baseFOLDER}/FastQC_raw
    mkdir ${baseFOLDER}/multiQC_raw

while read SAMPLE_ID
do

    for i in 1 2
    do
    
        gunzip ${baseFOLDER}/Muestras/${SAMPLE_ID}_R${i}_fastq.gz

        fastqc -o ${baseFOLDER}/FastQC_raw ${baseFOLDER}/Muestras/${SAMPLE_ID}_R${i}_fastq

        gzip ${baseFOLDER}/Muestras/${SAMPLE_ID}_R${i}_fastq
    
    done

done<$SAMPLES

    multiqc ${baseFOLDER}/FastQC_raw/ -o ${baseFOLDER}/multiQC_raw 

#    sh ${HOME}/Documentos/Scripts_modificados/DO_QCTRIMMING1.sh ${SAMPLES} ${LINEAGE} ${SET} ${REF} ${oriANALYTYPE}