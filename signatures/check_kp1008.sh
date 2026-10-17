#!/bin/bash
#SBATCH --job-name=check1008
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/sms049/shared/vselvanayagam/sig20m/check_kp1008.log
module use -p /group/shooshared/modules >/dev/null 2>&1
module load kmertools_s >/dev/null 2>&1
D=/group/shooshared/kpdata/anigo_filter
for f in $D/mt_filter/nuc_reads/KP1008.nuc_1.fastq.gz $D/mt_filter/nuc_reads/KP1008.nuc_2.fastq.gz \
         $D/filtered/noncp_reads20mio/KP1008_1.sub.fastq.gz $D/filtered/noncp_reads20mio/KP1008_2.sub.fastq.gz \
         $D/filtered/noncp_reads/KP1008.noncp_1.fastq.gz $D/filtered/noncp_reads/KP1008.noncp_2.fastq.gz \
         $D/raw/KP1008_1.fastq.gz $D/raw/KP1008_2.fastq.gz; do
  [ -e "$f" ] || { echo -e "$(basename $f)\tmissing"; continue; }
  g=$(gzip -t "$f" 2>&1 && echo ok || echo damaged)
  n=$(seqkit --quiet stats -T "$f" 2>/dev/null | awk -F'\t' 'NR==2{print $4}')
  echo -e "$(basename $f)\tsize=$(ls -l $f | awk '{printf "%.1fGB",$5/1e9}')\tgzip=$g\treads=${n:-unreadable}"
done
echo Done
