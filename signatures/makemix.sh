#!/bin/bash
#SBATCH --job-name=makemix
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-8%4
#SBATCH --output=/scratch/sms049/shared/vselvanayagam/sig20m/logs/makemix_%a.log
set -eo pipefail
module use -p /group/shooshared/modules >/dev/null 2>&1
module load kmertools_s >/dev/null 2>&1
OUT=/scratch/sms049/shared/vselvanayagam/mixreads
NL=/scratch/sms049/shared/vselvanayagam/nuc_links
MIXES=(FR10 FR25 FR50 FR75 FP10 FP25 FP50 FP75)
m=${MIXES[$((SLURM_ARRAY_TASK_ID-1))]}
pair=${m:0:2}; share=${m:2}
A=KP0055; if [ "$pair" = FR ]; then B=KP1067; else B=KP1066; fi
TOTAL=12000000
id=MIX_$m
echo "$id: $A + $B at $share% $B, total $TOTAL pairs, started $(date +%H:%M:%S)"
nA=$(seqkit stats -T $NL/$A.nuc_1.fastq.gz | awk -F'\t' 'NR==2{print $4}')
nB=$(seqkit stats -T $NL/$B.nuc_1.fastq.gz | awk -F'\t' 'NR==2{print $4}')
[ -n "$nA" ] && [ -n "$nB" ] || { echo "could not count reads"; exit 1; }
pA=$(awk -v s=$share -v t=$TOTAL -v n=$nA 'BEGIN{printf "%.8f", (100-s)/100*t/n}')
pB=$(awk -v s=$share -v t=$TOTAL -v n=$nB 'BEGIN{printf "%.8f", s/100*t/n}')
awk -v a=$pA -v b=$pB 'BEGIN{exit !(a<=1 && b<=1)}' || { echo "proportion above 1"; exit 1; }
for mate in 1 2; do
  { seqkit sample -s 42 -p $pA $NL/$A.nuc_$mate.fastq.gz 2>/dev/null
    seqkit sample -s 42 -p $pB $NL/$B.nuc_$mate.fastq.gz 2>/dev/null; } | pigz -p 4 > $OUT/$id.nuc_$mate.fastq.gz
done
seqkit stats -T $OUT/$id.nuc_1.fastq.gz $OUT/$id.nuc_2.fastq.gz
echo "$id done $(date +%H:%M:%S)"
