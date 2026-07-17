#!/bin/bash
#SBATCH --job-name=repeat_filter
#SBATCH --ntasks=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=/home/vselvanayagam/varshiny/logs/repeat_filter_%j.out
#SBATCH --error=/home/vselvanayagam/varshiny/logs/repeat_filter_%j.err

module load bbmap/39.17

REPEAT_REF="/group/shooshared/kpdata/anigo_filter/flavidus_repeat/consensi.fa.classified"
NONCP_DIR="/group/shooshared/kpdata/anigo_filter/filtered/noncp_reads"
OUT_DIR="/home/vselvanayagam/varshiny/filtered_repeats"
mkdir -p $OUT_DIR

SAMPLES=('KP0055' 'KP0052' 'KP0069' 'KP0051')

for SAMPLE in "${SAMPLES[@]}"; do
    REPEAT_R1="${OUT_DIR}/${SAMPLE}.repeat_1.fastq.gz"

    if [ -f "$REPEAT_R1" ]; then
        echo "  [SKIP] $SAMPLE - already filtered"
        continue
    fi

    echo "Processing: $SAMPLE"

    bbmap.sh \
        ref=$REPEAT_REF \
        in=${NONCP_DIR}/${SAMPLE}.noncp_1.fastq.gz \
        in2=${NONCP_DIR}/${SAMPLE}.noncp_2.fastq.gz \
        outm=${OUT_DIR}/${SAMPLE}.repeat_1.fastq.gz \
        outm2=${OUT_DIR}/${SAMPLE}.repeat_2.fastq.gz \
        minid=0.95 \
        threads=8 \
        statsfile=${OUT_DIR}/${SAMPLE}.repeat_stats.txt \
        -Xmx60g

    echo "Done: $SAMPLE"
done

echo "All samples complete! Results: $OUT_DIR"
