#!/usr/bin/bash

set -euo pipefail #El pipeline se para si hay algún error (-e), si hay alguna variable no definida (-u), si alguno de los comandos de la tuberia falla (-o pipefail)

trap 'echo "❌ Error en la línea $LINENO durante la ejecución del pipeline"; exit 1' ERR #trap activa el mensaje si hay un error (ERR)

SAMPLES=$1
LINEAGE=$2
#SET=$3
#REF=$4  #tengo una única referencia, por lo que no hace falta ponerlo con estas muestras

#Previo a todo el pipeline PONER LAS MUESTRAS EN UNA CARPETA LLAMADA MUESTRA.

baseFOLDER="${HOME}/Proyectos/${LINEAGE}"
originalfastaFOLDER="${HOME}/Documentos/Referencias_originales"
indexed_ref_FOLDER="${HOME}/indexed_libs"
bowtie2_samtools_FOLDER="${HOME}/miniconda3/envs/PipelineIdisba2/bin"
javaFOLDER="${HOME}/miniconda3/envs/java17/bin"
picardFOLDER="${HOME}/miniconda3/envs/PipelineIdisba2/share/picard-3.4.0-0" 
GATKfolder="${HOME}/miniconda3/envs/PipelineIdisba2/share/gatk4-4.6.2.0-0"
snpEffFOLDER="${HOME}/miniconda3/envs/snpeff_only/share/snpeff-5.2-1" #Descargado en un env de conda diferente
snpsiftFOLDER="${HOME}/miniconda3/envs/PipelineIdisba2/share/snpsift-5.2-0"

mkdir -p ${baseFOLDER}/sam_files
mkdir -p ${baseFOLDER}/rawpileup_files/
mkdir -p ${baseFOLDER}/totalpileup_files/
mkdir -p ${baseFOLDER}/annotated_files/
mkdir -p ${baseFOLDER}/curated_files/
mkdir -p ${baseFOLDER}/intermediate_files/


while read SAMPLE_ID
do
    SAMPLE_ID=${SAMPLE_ID//$'\r'/}
    echo "Procesando muestra ${SAMPLE_ID}" 
  
    echo "===BOWTIE2 analysis for sample ${SAMPLE_ID}===" #La carpeta del bowtie2 y samtools es la misma

    gunzip ${baseFOLDER}/Muestras/${SAMPLE_ID}_L001_R1_001.fastq.gz
    gunzip ${baseFOLDER}/Muestras/${SAMPLE_ID}_L001_R2_001.fastq.gz

    ${bowtie2_samtools_FOLDER}/bowtie2 \
    --phred33 \
    -x "${indexed_ref_FOLDER}/ATCC19424/ATCC19424" \
    -q \
    -1 "${baseFOLDER}/Muestras/${SAMPLE_ID}_L001_R1_001.fastq" \
    -2 "${baseFOLDER}/Muestras/${SAMPLE_ID}_L001_R2_001.fastq" \
    -X 500 \
    -S "${baseFOLDER}/sam_files/${SAMPLE_ID}.sam"

    # gzip ${baseFOLDER}/Muestras/${SAMPLE_ID}_L001_R1_001.fastq #Dejar comentado esto para no volver a comprimir
    # gzip ${baseFOLDER}/Muestras/${SAMPLE_ID}_L001_R2_001.fastq


    echo "BOWTIE2 Analysis DONE!"

    echo " "
    echo " "

    echo "===SAMtools Analysis for sample ${SAMPLE_ID}==="

    echo " "
    echo "Generating ${SAMPLE_ID} BAM file" 
    

    ${bowtie2_samtools_FOLDER}/samtools \
    view -bS ${baseFOLDER}/sam_files/${SAMPLE_ID}.sam > ${baseFOLDER}/intermediate_files/${SAMPLE_ID}.bam

    echo "BAM file DONE!"
   
    echo " "
    echo "Sorting ${SAMPLE_ID} sample"

    ${bowtie2_samtools_FOLDER}/samtools \
    sort ${baseFOLDER}/intermediate_files/${SAMPLE_ID}.bam \
    -o ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_mapped.sorted.bam
  
    echo " "
    echo "$SAMPLE_ID sorted" 

    echo " "
    echo "===PICARDTOOLS Analysis==="
    echo "Adding ReadGroups for ${SAMPLE_ID} sample" #En el programa actualizado, primero se deben poner los RG y luego marcar los duplicados
       
    fastq_header=`head -n1 ${baseFOLDER}/Muestras/${SAMPLE_ID}_L001_R1_001.fastq`
    echo ${fastq_header}

    FLOWCELL=$(echo "$fastq_header" | cut -d':' -f3)
    SEQLANE=$(echo "$fastq_header" | cut -d':' -f4)
    INDEX=$(echo "$fastq_header" | cut -d':' -f10)

    ${javaFOLDER}/java -jar $picardFOLDER/picard.jar \
    AddOrReplaceReadGroups \
    -I ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_mapped.sorted.bam \
    -O ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_mapped.sorted_RDG.bam \
    -RGID ${FLOWCELL}.${SEQLANE} \
    -RGLB lib1 \
    -RGPL ILLUMINA \
    -RGPU ${FLOWCELL}.${SEQLANE}.${INDEX} \
    -RGSM ${SAMPLE_ID}

    echo "Marking duplicates for ${SAMPLE_ID} sample" 
    
    ${javaFOLDER}/java -jar $picardFOLDER/picard.jar \
    MarkDuplicates -INPUT ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_mapped.sorted_RDG.bam \
    -OUTPUT ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_mapped.sorted_RDG.dedup.bam \
    -METRICS_FILE ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_markDuplicatesMetrics.txt \
    -ASSUME_SORTED True 

    echo " "
    echo "===Creating Index Needed for ${SAMPLE_ID} sample==="

    ${bowtie2_samtools_FOLDER}/samtools \
    index ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_mapped.sorted_RDG.dedup.bam
    
    echo " "
    echo "===VARIANT CALLING for ${SAMPLE_ID}===" 
    echo "1- Creating reference .fai and .dict"  #Para que funcione, se tiene que haber creado previamente el .fai (samtools faidx) y .dict de la referencia fasta (picard CreateSequenceDictionary)

        
    ${bowtie2_samtools_FOLDER}/samtools \
    faidx ${originalfastaFOLDER}/Neisseria_gonorrhoeae_ATCC_19424.fasta

    ${javaFOLDER}/java -jar $picardFOLDER/picard.jar \
    CreateSequenceDictionary \
    -R ${originalfastaFOLDER}/Neisseria_gonorrhoeae_ATCC_19424.fasta \
    -O ${originalfastaFOLDER}/Neisseria_gonorrhoeae_ATCC_19424.dict  #Si ya está hecho el archivo, da error!

    echo "2- Variant Calling" 

    ${javaFOLDER}/java -jar ${GATKfolder}/gatk-package-4.6.2.0-local.jar HaplotypeCaller \
    -R ${originalfastaFOLDER}/Neisseria_gonorrhoeae_ATCC_19424.fasta \
    -I ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_mapped.sorted_RDG.dedup.bam \
    -O ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_mapped.sorted_RDG.dedup.vcf #llamada de variantes

    echo "3- Creating quality labels" 

    ${javaFOLDER}/java -jar ${GATKfolder}/gatk-package-4.6.2.0-local.jar VariantFiltration \
    -R ${originalfastaFOLDER}/Neisseria_gonorrhoeae_ATCC_19424.fasta \
    -V ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_mapped.sorted_RDG.dedup.vcf \
    -O ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_filtered_variants.vcf \
    --filter-name "QD_filter" --filter-expression "QD < 2.0" --filter-name "FS_filter" --filter-expression "FS > 60.0" #Crear las etiquetas para las calidades de las variantes

    echo "4- Eliminating low quality variants"
    echo " "

    ${javaFOLDER}/java -jar ${GATKfolder}/gatk-package-4.6.2.0-local.jar SelectVariants \
    -R ${originalfastaFOLDER}/Neisseria_gonorrhoeae_ATCC_19424.fasta \
    -V ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_filtered_variants.vcf \
    --exclude-filtered -O ${baseFOLDER}/rawpileup_files/${SAMPLE_ID}_passed_variants.vcf

    echo "Generating ${SAMPLE_ID} totalpileup"
    samtools pileup -c -f ${originalfastaFOLDER}/Neisseria_gonorrhoeae_ATCC_19424.fasta \
    ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_mapped.sorted.dedup_RDG_realignedintervals_realignedID.bam \
    > ${baseFOLDER}/totalpileup_files/${SAMPLE_ID}_${TRIM}.realigned.totalpileup
    
    echo "DONE Variant Calling"
    echo " "
    echo "===ANNOTATING sample ${SAMPLE_ID}==="

    ${javaFOLDER}/java -jar ${snpEffFOLDER}/snpEff.jar \
    NG19424 \
    ${baseFOLDER}/rawpileup_files/${SAMPLE_ID}_passed_variants.vcf > ${baseFOLDER}/annotated_files/${SAMPLE_ID}_ann.vcf

    echo "DONE Annotating sample ${SAMPLE_ID}"

    echo " "
    # echo "===Extracting DATA==="

    # cat ${baseFOLDER}/annotated_files/${SAMPLE_ID}_ann.vcf | perl ${snpEffFOLDER}/scripts/vcfEffOnePerLine.pl > ${HOME}/tmpX
        
    # ${javaFOLDER}/java -jar ${snpsiftFOLDER}/SnpSift.jar extractFields ${HOME}/tmpX CHROM POS REF ALT ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].GENEID ANN[*].HGVS_C ANN[*].HGVS_P ANN[*].CDNA_POS ANN[*].CDNA_LEN ANN[*].AA_POS ANN[*].AA_LEN > ${baseFOLDER}/curated_files/${SAMPLE_ID}_pre_curated

    # cat ${baseFOLDER}/curated_files/${SAMPLE_ID}_pre_curated | grep -v 'downstream_gene_variant' | grep -v 'intergenic_region' | grep -v 'synonymous_variant' | grep -v 'upstream_gene_variant' > ${baseFOLDER}/curated_files/${SAMPLE_ID}_curated

done<$SAMPLES