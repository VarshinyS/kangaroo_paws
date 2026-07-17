#!/bin/bash
#SBATCH --job-name=mito_filter
#SBATCH --ntasks=8
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=/home/vselvanayagam/varshiny/logs/mito_filter_%j.out
#SBATCH --error=/home/vselvanayagam/varshiny/logs/mito_filter_%j.err

module load bbmap/39.17

MITO_REF="/group/shooshared/kpdata/anigo_filter/raw_mt/KP9002.fa"
NONCP_DIR="/group/shooshared/kpdata/anigo_filter/filtered/noncp_reads"
OUT_DIR="/home/vselvanayagam/varshiny/filtered_noncp"
mkdir -p $OUT_DIR

SAMPLES=('KP0055' 'KP0052' 'KP0069' 'KP0051')

for SAMPLE in "${SAMPLES[@]}"; do
    MITO_R1="${OUT_DIR}/${SAMPLE}.mito_1.fastq.gz"

    if [ -f "$MITO_R1" ]; then
        echo "  [SKIP] $SAMPLE - already filtered"
        continue
    fi

    echo "Processing: $SAMPLE"

    bbmap.sh \
        ref=$MITO_REF \
        in=${NONCP_DIR}/${SAMPLE}.noncp_1.fastq.gz \
        in2=${NONCP_DIR}/${SAMPLE}.noncp_2.fastq.gz \
        outm=${OUT_DIR}/${SAMPLE}.mito_1.fastq.gz \
        outm2=${OUT_DIR}/${SAMPLE}.mito_2.fastq.gz \
        minid=0.95 \
        threads=8 \
        statsfile=${OUT_DIR}/${SAMPLE}.mito_stats.txt \
        -Xmx60g

    echo "Done: $SAMPLE"
done

echo "All samples complete! Results: $OUT_DIR"
