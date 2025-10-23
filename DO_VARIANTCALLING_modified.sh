#!/usr/bin/bash

SAMPLES=$1
LINEAGE=$2
#SET=$3
#REF=$4  #tengo una única referencia, por lo que no hace falta ponerlo con estas muestras

baseFOLDER="${HOME}/Proyectos/${LINEAGE}"
originalfastaFOLDER="${HOME}/Documentos/Referencias_originales"
indexed_ref_FOLDER="${HOME}/indexed_libs"
#fastqcFOLDER="${HOME}/Programs/FastQC" #Todos estos que están comentados los descargué directamente en el env de conda con conda install
#fastpFOLDER="${HOME}/anaconda3/bin"
#multiqcFOLDER="${HOME}/anaconda3/bin"
#bowtie2FOLDER="${HOME}/Programs/bowtie2-2.5.4"
#javaFOLDER="${HOME}/Programs/jre1.8.0_331/bin"
#samtoolsFOLDER="${HOME}/Programas/samtools-0.1.16" #Este programa es muy viejo y se tiene que compilar el instalador
picardFOLDER="${HOME}/Programas/picard-tools-1.140" 
GATKfolder="${HOME}/Programas/GenomeAnalysisTK-3.4-46"
#snpEffFOLDER="${HOME}/miniconda3/envs/PipelineIdisba/share/snpeff-5.2-1"
oriANALYTYPE="BOTH" #UT/TR/BOTH

while read SAMPLE_ID
do

    if [ ${oriANALYTYPE} = "BOTH" ]
    then

        ANALYTYPE="UT TR"

    fi
    
    echo "Procesando muestra ${SAMPLE_ID}" 

    mkdir ${baseFOLDER}/sam_files
    mkdir ${baseFOLDER}/rawpileup_files/
    mkdir ${baseFOLDER}/totalpileup_files/
    mkdir ${baseFOLDER}/annotated_files/
    mkdir ${baseFOLDER}/curated_files/


    for TRIM in ${ANALYTYPE};
    do

        if [ ${TRIM} = "UT" ]; then           

            echo "BOWTIE2 analysis for UNTRIMMED (UT) sample ${SAMPLE_ID}" 

            gunzip ${baseFOLDER}/Muestras/${SAMPLE_ID}_R1_fastq.gz
            gunzip ${baseFOLDER}/Muestras/${SAMPLE_ID}_R2_fastq.gz

            bowtie2 \
            --phred33 \
            -x "${indexed_ref_FOLDER}/ATCC19424/ATCC19424" \
            -q \
            -1 "${baseFOLDER}/Muestras/${SAMPLE_ID}_R1_fastq" \
            -2 "${baseFOLDER}/Muestras/${SAMPLE_ID}_R2_fastq" \
            -X 500 \
            -S "${baseFOLDER}/sam_files/${SAMPLE_ID}_${TRIM}.sam"

            # gzip ${baseFOLDER}/Muestras/${SAMPLE_ID}_R1_fastq.gz #Hay que quitar el .gz porque si no no lo reconoce
            # gzip ${baseFOLDER}/Muestras/${SAMPLE_ID}_R2_fastq.gz

        elif [ ${TRIM} = "TR" ]; then
        
            echo "BOWTIE2 analysis for TRIMMED (TR) sample ${SAMPLE_ID}"

            bowtie2 \
            --phred33 \
            -x "${indexed_ref_FOLDER}/ATCC19424/ATCC19424" \
            -q \
            -1 "${baseFOLDER}/FastP_TM_fastq/${SAMPLE_ID}_L001_R1_001.trimmed.fastq" \
            -2 "${baseFOLDER}/FastP_TM_fastq/${SAMPLE_ID}_L001_R2_001.trimmed.fastq" \
            -X 500 \
            -S "${baseFOLDER}/sam_files/${SAMPLE_ID}_${TRIM}.sam"

            #gzip ${baseFOLDER}/FastP_TM_fastq/${SAMPLE_ID}_L001_R1_001.trimmed.fastq 
            #gzip ${baseFOLDER}/FastP_TM_fastq/${SAMPLE_ID}_L001_R2_001.trimmed.fastq

        fi

        echo "BOWTIE2 Analysis DONE!"

        echo " "
        echo " "

    done

    for TRIM in ${ANALYTYPE}
    do

        if [ ${TRIM} = "UT" ]
        then

            TRIMESSAGE="UNTRIMMED (UT)"

        fi

        if [ ${TRIM} = "TR" ]

        then

            TRIMESSAGE="TRIMMED (TR)"

        fi

        echo " "
        echo "SAMtools Analysis for ${TRIMESSAGE} sample ${SAMPLE_ID}"

        mkdir ${baseFOLDER}/intermediate_files/
       
        echo " "
        echo "Generating ${SAMPLE_ID}_${TRIM} BAM file" 
    
        #gunzip ${baseFOLDER}/sam_files/${SAMPLE_ID}_${TRIM}.sam.gz #Creo que no se crea el archivo .gz?

        samtools \
        view -bS ${baseFOLDER}/sam_files/${SAMPLE_ID}_${TRIM}.sam > ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}.bam

        #gzip ${baseFOLDER}/sam_files/${SAMPLE_ID}_${TRIM}_map${REF}.sam

        echo " "
        echo "Sorting ${TRIMESSAGE}_${SAMPLE_ID} sample"

        samtools \
        sort ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}.bam \
        "${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_mapped.sorted"
  
        echo " "
        echo "$SAMPLE_ID sorted" 

        echo " "
        echo "Marking duplicates for ${TRIMESSAGE}_${SAMPLE_ID} sample"

        java -jar $picardFOLDER/picard.jar \
        MarkDuplicates INPUT=${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_mapped.sorted.bam \
        OUTPUT=${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_mapped.sorted.dedup.bam \
        METRICS_FILE=${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_markDuplicatesMetrics.txt \
        ASSUME_SORTED=True 

        fastq_header=$(zcat ${baseFOLDER}/Muestras/${SAMPLE_ID}_R1_fastq.gz |  head -1)

        echo " "
        echo "Adding ReadGroups for ${TRIMESSAGE}_${SAMPLE_ID} sample"

        FLOWCELL=$(echo "$fastq_header" | cut -d':' -f3)
        SEQLANE=$(echo "$fastq_header" | cut -d':' -f4)
        INDEX=$(echo "$fastq_header" | cut -d':' -f10)

        java -jar $picardFOLDER/picard.jar \
        AddOrReplaceReadGroups \
        I=${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_mapped.sorted.dedup.bam \
        O=${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_mapped.sorted.dedup_RDG.bam \
        RGID=${FLOWCELL}.${SEQLANE} \
        RGLB=lib1 \
        RGPL=ILLUMINA \
        RGPU=${FLOWCELL}.${LANE}.${INDEX} \
        RGSM=${SAMPLE_ID} 

        echo " "
        echo "Creating Index Needed for ${TRIMESSAGE}_${SAMPLE_ID} sample"

        java -jar ${picardFOLDER}/picard.jar \
        BuildBamIndex I=${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_mapped.sorted.dedup_RDG.bam 
    
        echo " "
        echo "VARIANT CALLING for ${TRIMESSAGE}_${SAMPLE_ID}. Creating Intervals" #Para que funcione, se tiene que haber creado previamente el .fai (samtools faidx) y .dict de la referencia fasta (picard CreateSequenceDictionary)

        java -jar $GATKfolder/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator -R ${originalfastaFOLDER}/Neisseria_gonorrhoeae_ATCC_19424.fasta \
        -I ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_mapped.sorted.dedup_RDG.bam \
        -o ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_mapped.sorted.dedup_RDG_realignedintervals.bed

        echo " "
        echo "VARIANT CALLING for ${TRIMESSAGE}_${SAMPLE_ID} sample. InDels Realignment"

        java -jar $GATKfolder/GenomeAnalysisTK.jar \
        -T IndelRealigner -R ${originalfastaFOLDER}/Neisseria_gonorrhoeae_ATCC_19424.fasta \
        -I ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_mapped.sorted.dedup_RDG.bam \
        -targetIntervals ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_mapped.sorted.dedup_RDG_realignedintervals.bed \
        -o ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_mapped.sorted.dedup_RDG_realignedintervals_realignedID.bam

        echo " "
        echo "Generating ${SAMPLE_ID} totalpileup"
        samtools pileup -c -f ${originalfastaFOLDER}/Neisseria_gonorrhoeae_ATCC_19424.fasta \
        ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_mapped.sorted.dedup_RDG_realignedintervals_realignedID.bam \
        > ${baseFOLDER}/totalpileup_files/${SAMPLE_ID}_${TRIM}.realigned.totalpileup
        

        echo " "
        echo "VARIANT CALLING for ${TRIMESSAGE}_${SAMPLE_ID} sample. RAW pileup" 

        java -jar $GATKfolder/GenomeAnalysisTK.jar \
        -T UnifiedGenotyper -R ${originalfastaFOLDER}/Neisseria_gonorrhoeae_ATCC_19424.fasta \
        -I ${baseFOLDER}/intermediate_files/${SAMPLE_ID}_${TRIM}_mapped.sorted.dedup_RDG_realignedintervals_realignedID.bam \
        -o ${baseFOLDER}/rawpileup_files/${SAMPLE_ID}_${TRIM}.rawpileup -ploidy 1 -stand_call_conf 30 -stand_emit_conf 10 -glm BOTH 

        echo "DONE Analysis for ${TRIMESSAGE} sample ${SAMPLE_ID}"

        # echo " "
        # echo "ANNOTATING sample ${SAMPLE_ID}_${ANALYTYPE}"

        # #El siguiente comando tiene un alias --> alias snpeff='java -jar ~/miniconda3/envs/PipelineIdisba/share/snpeff-5.2-1/snpEff.jar'
        # snpeff NG19424 ${baseFOLDER}/rawpileup_files/${SAMPLE_ID}_${TRIM}.rawpileup > ${baseFOLDER}/annotated_files/${SAMPLE_ID}_${TRIM}.ann.vcf

        # echo "DONE Annotating sample ${SAMPLE_ID}_${ANALYTYPE}"

        # echo " "
        # echo "Extracting DATA"

        # cat ${baseFOLDER}/annotated_files/${SAMPLE_ID}_${TRIM}.ann.vcf | perl ${snpEffFOLDER}/scripts/vcfEffOnePerLine.pl > ${HOME}/tmpX
        
        # #El siguiente comando tiene un alias --> alias SnpSift='java -jar ~/miniconda3/envs/PipelineIdisba/share/snpsift-5.2-0/SnpSift.jar'
        # SnpSift extractFields ${HOME}/tmpX CHROM POS REF ALT ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].GENEID ANN[*].HGVS_C ANN[*].HGVS_P ANN[*].CDNA_POS ANN[*].CDNA_LEN ANN[*].AA_POS ANN[*].AA_LEN > ${baseFOLDER}/curated_files/${SAMPLE_ID}_${TRIM}.pre_curated

        # cat ${baseFOLDER}/curated_files/${SAMPLE_ID}_${TRIM}.pre_curated | grep -v 'downstream_gene_variant' | grep -v 'intergenic_region' | grep -v 'synonymous_variant' | grep -v 'upstream_gene_variant' > ${baseFOLDER}/curated_files/${SAMPLE_ID}_${TRIM}.curated

    done

done<$SAMPLES
