#!/bin/bash
#SBATCH --job-name=nuc_sketch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-53%8
#SBATCH --output=/scratch/sms049/shared/vselvanayagam/rebuild/logs_nuc/sketch_%a.log

module use -p /group/shooshared/modules
module load kmertools_s

MASH=/group/shooshared/micromamba/envs/kmertools/bin/mash
WORK=/scratch/sms049/shared/vselvanayagam/rebuild
DIR=/group/shooshared/kpdata/anigo_filter/mt_filter/nuc_reads
N=20000000
tag=nuc
cd "$WORK" || exit 1
mkdir -p sketches_nuc logs_nuc

fail() { echo "FAILED: $*"; exit 1; }

IDS=($(cat ids53.txt))
id=${IDS[$((SLURM_ARRAY_TASK_ID-1))]}
echo "sample $id started $(date +%H:%M:%S)"
[ -s "sketches_nuc/$id.msh" ] && { echo "$id already done"; exit 0; }

r1="$DIR/$id.nuc_1.fastq.gz"; r2="$DIR/$id.nuc_2.fastq.gz"
[ -s "$r1" ] || [ -s "$r2" ] || fail "no nuclear reads for $id"

TMP=$WORK/tmp_nuc_$id
rm -rf "$TMP"; mkdir -p "$TMP"

subsample() {
  local f=$1 o=$2 n=$3 tot p
  tot=$(seqkit stats -T "$f" | awk -F'\t' 'NR==2{print $4}')
  if [ -z "$tot" ]; then echo "could not count reads in $f"; return 1; fi
  echo "reads in $(basename "$f"): $tot"
  if [ "$tot" -le "$n" ]; then cp "$f" "$o"; return 0; fi
  p=$(awk -v n="$n" -v t="$tot" 'BEGIN{printf "%.8f", n/t}')
  seqkit sample -p "$p" -s 42 "$f" -o "$o" 2>> "logs_nuc/seqkit_$id.err"
}

ok1=0; ok2=0
n1=""; n2=""
[ -s "$r1" ] && n1=$(seqkit stats -T "$r1" 2>/dev/null | awk -F'\t' 'NR==2{print $4}'); [ -n "$n1" ] && ok1=1
[ -s "$r2" ] && n2=$(seqkit stats -T "$r2" 2>/dev/null | awk -F'\t' 'NR==2{print $4}'); [ -n "$n2" ] && ok2=1
if [ $ok1 -eq 1 ] && [ $ok2 -eq 1 ]; then
  subsample "$r1" "$TMP/a.fq.gz" "$N" || fail "subsample $r1"
  subsample "$r2" "$TMP/b.fq.gz" "$N" || fail "subsample $r2"
  cat "$TMP/a.fq.gz" "$TMP/b.fq.gz" > "$TMP/$id.fq.gz"
elif [ $ok1 -eq 1 ] || [ $ok2 -eq 1 ]; then
  if [ $ok1 -eq 1 ]; then one=$r1; else one=$r2; fi
  echo "$id: one mate is damaged, using $(basename "$one") only" | tee -a damaged_nuc.txt
  subsample "$one" "$TMP/$id.fq.gz" "$((2*N))" || fail "subsample $one"
else
  echo "$id: both mates damaged, sample skipped" | tee -a damaged_nuc.txt
  exit 0
fi

"$MASH" sketch -r -s 500000 -k 21 -o "sketches_nuc/$id" "$TMP/$id.fq.gz" > "logs_nuc/mash_$id.out" 2>&1 || { cat "logs_nuc/mash_$id.out"; fail "mash sketch on $id"; }
[ -s "sketches_nuc/$id.msh" ] || { cat "logs_nuc/mash_$id.out"; fail "empty sketch for $id"; }
rm -rf "$TMP"
echo "$id finished $(date +%H:%M:%S)"
