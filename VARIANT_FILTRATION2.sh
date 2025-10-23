#!/usr/bin/env bash

set -euo pipefail
trap 'echo "❌ Error en la línea $LINENO durante la ejecución del pipeline"; exit 1' ERR

MAP_FILE=$1 #Lista de las muestras a comparar con el parental    
LINEAGE=$2
#PARENTAL=$3 #Nombre del parental

baseFOLDER="${HOME}/Proyectos/${LINEAGE}" #Cambiar lo de proyectos en el Linux
bcftools="$HOME/miniconda3/envs/PipelineIdisba2/bin/bcftools"
javaFOLDER="${HOME}/miniconda3/envs/java17/bin"
snpsiftFOLDER="${HOME}/miniconda3/envs/PipelineIdisba2/share/snpsift-5.2-0"
snpEffFOLDER="${HOME}/miniconda3/envs/snpeff_only/share/snpeff-5.2-1"

while read -r SAMPLE_ID PARENTAL; do
    SAMPLE_ID=${SAMPLE_ID//$'\r'/}   # eliminar retorno de carro
    PARENTAL=${PARENTAL//$'\r'/} #CONTINUAR AQUI

    mkdir -p ${baseFOLDER}/${PARENTAL}.anotado
    mkdir -p ${baseFOLDER}/curated_files_${PARENTAL}

    dir_anotado="${baseFOLDER}/${PARENTAL}.anotado"

    # Copiar parental
    cp "${HOME}/Proyectos/${LINEAGE}/annotated_files/${PARENTAL}_ann.vcf" "${dir_anotado}/"

    # Comprimir e indexar parental si no está comprimido

    echo "Comprobando si se comprime e indexa el PARENTAL ${PARENTAL}"
    if [ ! -f "${dir_anotado}/${PARENTAL}_ann.vcf.gz" ]; then
        echo " Comprimiendo ${PARENTAL}..."
        bgzip -f "${dir_anotado}/${PARENTAL}_ann.vcf"
    else
        echo "${PARENTAL} YA ESTÁ COMPRIMIDO"
    fi
    if [ ! -f "${dir_anotado}/${PARENTAL}_ann.vcf.gz.tbi" ]; then
        echo " Indexando ${PARENTAL}..."
        tabix -p vcf "${dir_anotado}/${PARENTAL}_ann.vcf.gz"
    else
        echo "${PARENTAL} YA ESTÁ INDEXADO"
    fi

    echo "=== Procesando muestra ${SAMPLE_ID} ==="

    #Copiar muestras a directorio anotado
    cp "${HOME}/Proyectos/${LINEAGE}/annotated_files/${SAMPLE_ID}_ann.vcf" "${dir_anotado}/"

    sample_file="${dir_anotado}/${SAMPLE_ID}_ann.vcf"

    # Comprimir e indexar la muestra
    echo "Comprobando si se comprime e indexa la MUESTRA ${SAMPLE_ID}"
    if [ ! -f "${sample_file}.gz" ]; then
        echo "Comprimiendo ${SAMPLE_ID}..."
        bgzip -f "${sample_file}"
    else 
        echo "${SAMPLE_ID}.gz YA ESTÁ COMPRIMIDO"
    fi

    if [ ! -f "${sample_file}.gz.tbi" ]; then
        echo "Indexando ${SAMPLE_ID}.gz..."
        tabix -p vcf "${sample_file}.gz"
    else
        echo "${SAMPLE_ID}.gz.tbi YA ESTÁ INDEXADO"
    fi

    #Excluir el parental del análisis (para que no se analice con él mismo)

    if [[ "${sample_file}.gz" == "${PARENTAL}_ann.vcf.gz" ]]; then
            continue
    fi

    #SAMPLE_NAME=$(basename "${sample_file}" _ann.vcf.gz)
    outdir="${dir_anotado}/comparacion_${SAMPLE_ID}"
    mkdir -p "$outdir"

    echo "=== Comparando ${SAMPLE_ID} con parental ${PARENTAL} ==="
    "$bcftools" isec -p "$outdir" \
    -n=1 \
    -w1 \
    "${sample_file}.gz" "${dir_anotado}/${PARENTAL}_ann.vcf.gz"

    #Renombrar el archivo de salida

    if [ -f "$outdir"/0000.vcf ]; then
        mv "$outdir"/0000.vcf "$outdir"/${SAMPLE_ID}_unique.vcf
        echo "Archivo renombrado: ${SAMPLE_ID}_unique.vcf"

        echo "===Extracting DATA==="

        cat "$outdir"/${SAMPLE_ID}_unique.vcf | perl ${snpEffFOLDER}/scripts/vcfEffOnePerLine.pl > ${HOME}/tmpX

        ${javaFOLDER}/java -jar ${snpsiftFOLDER}/SnpSift.jar extractFields ${HOME}/tmpX CHROM POS REF ALT ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].GENEID ANN[*].HGVS_C ANN[*].HGVS_P ANN[*].CDNA_POS ANN[*].CDNA_LEN ANN[*].AA_POS ANN[*].AA_LEN > ${baseFOLDER}/curated_files_${PARENTAL}/${SAMPLE_ID}_pre_curated

        cat ${baseFOLDER}/curated_files_${PARENTAL}/${SAMPLE_ID}_pre_curated | grep -v 'downstream_gene_variant' | grep -v 'intergenic_region' | grep -v 'synonymous_variant' | grep -v 'upstream_gene_variant' > ${baseFOLDER}/curated_files_${PARENTAL}/${SAMPLE_ID}_curated

    else 
        echo "No se encontro variantes únicas para $SAMPLE_ID"
    fi

    echo "Comparación completada para ${SAMPLE_ID}"
    
done < "$MAP_FILE"
