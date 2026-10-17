#!/bin/bash
#SBATCH --job-name=minid_test
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/sms049/shared/vselvanayagam/sig20m/minid_test.log
module use -p /group/shooshared/modules >/dev/null 2>&1
module load kmertools_s >/dev/null 2>&1
BBMAP=/uwahpc/rocky9/bio-apps/binary/bbmap/39.17/bin/bbmap.sh
LIB=/group/shooshared/kpdata/anigo_filter/flavidus/repeat/flavidus_repeat_library.classfied
D=/group/shooshared/kpdata/anigo_filter/mt_filter/nuc_reads
T=/scratch/sms049/shared/vselvanayagam/sig20m/tmp_minid; mkdir -p $T
for id in KP0052 KP1067; do
  seqkit head -n 1000000 $D/$id.nuc_1.fastq.gz -o $T/${id}_1.fq.gz
  seqkit head -n 1000000 $D/$id.nuc_2.fastq.gz -o $T/${id}_2.fq.gz
  for m in 0.95 0.90 0.85 0.80 0.76; do
    $BBMAP in1=$T/${id}_1.fq.gz in2=$T/${id}_2.fq.gz ref=$LIB minid=$m nodisk overwrite=t fastareadlen=500 t=16 -Xmx24g 2> $T/log_${id}_$m.txt
    echo -e "$id\tminid=$m\tmapped=$(grep -m1 -E '^mapped:' $T/log_${id}_$m.txt | awk '{print $2}')"
  done
done
rm -rf $T
