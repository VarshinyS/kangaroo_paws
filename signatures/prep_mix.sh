#!/bin/bash
#SBATCH --job-name=prep_mix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --array=1-8%4
#SBATCH --output=/scratch/sms049/shared/vselvanayagam/sig20m/logs/prepmix_%a.log

module use -p /group/shooshared/modules
module load kmertools_s

WORK=/scratch/sms049/shared/vselvanayagam/sig20m
BBMAP=/uwahpc/rocky9/bio-apps/binary/bbmap/39.17/bin/bbmap.sh
LIB=/group/shooshared/kpdata/anigo_filter/flavidus/repeat/flavidus_repeat_library.classfied
TARGET=20000000
CI=1
cd "$WORK" || exit 1
mkdir -p kmers stats logs

fail() { echo "FAILED: $*"; exit 1; }
[ -s "$LIB" ] || fail "repeat library not found: $LIB"
[ -x "$BBMAP" ] || fail "bbmap not found: $BBMAP"

IDS=($(cat ids16.txt))
IDS=(MIX_FR10 MIX_FR25 MIX_FR50 MIX_FR75 MIX_FP10 MIX_FP25 MIX_FP50 MIX_FP75)
id=${IDS[$((SLURM_ARRAY_TASK_ID-1))]}
echo "sample $id (task $SLURM_ARRAY_TASK_ID) started $(date +%H:%M:%S)"
[ -s "kmers/$id.kmc_pre" ] && { echo "$id already counted"; exit 0; }

r1=""
for d in /scratch/sms049/shared/vselvanayagam/mixreads; do
  [ -s "$d/$id.nuc_1.fastq.gz" ] && { r1="$d/$id.nuc_1.fastq.gz"; break; }
done
[ -n "$r1" ] || fail "no nuclear reads found for $id"
r2=${r1/.nuc_1./.nuc_2.}
[ -s "$r2" ] || fail "second mate missing for $id"
echo "input: $r1"

TMP=$WORK/tmp_$id
rm -rf "$TMP"; mkdir -p "$TMP"

tot=$(seqkit stats -T "$r1" | awk -F'\t' 'NR==2{print $4}')
[ -n "$tot" ] || fail "could not count reads in $r1"
echo "nuclear read pairs available: $tot"
if [ "$tot" -le "$TARGET" ]; then
  cp "$r1" "$TMP/a1.fq.gz"; cp "$r2" "$TMP/a2.fq.gz"
else
  p=$(awk -v n="$TARGET" -v t="$tot" 'BEGIN{printf "%.8f", n/t}')
  seqkit sample -p "$p" -s 42 "$r1" -o "$TMP/a1.fq.gz" || fail "subsample R1"
  seqkit sample -p "$p" -s 42 "$r2" -o "$TMP/a2.fq.gz" || fail "subsample R2"
fi
sub=$(seqkit stats -T "$TMP/a1.fq.gz" | awk -F'\t' 'NR==2{print $4}')
echo "subsampled read pairs: $sub"

"$BBMAP" in1="$TMP/a1.fq.gz" in2="$TMP/a2.fq.gz" ref="$LIB" minid=0.95 \
  outu1="$TMP/u1.fq.gz" outu2="$TMP/u2.fq.gz" nodisk overwrite=t fastareadlen=500 \
  t=16 -Xmx24g 2> "logs/bbmap_$id.txt" || { tail -5 "logs/bbmap_$id.txt"; fail "bbmap"; }
mapped=$(grep -m1 -E "^mapped:" "logs/bbmap_$id.txt" | awk '{print $2}')
kept=$(seqkit stats -T "$TMP/u1.fq.gz" | awk -F'\t' 'NR==2{print $4}')
echo "repeat-mapped (reported by bbmap): $mapped; repeat-filtered read pairs kept: $kept"

printf '%s\n%s\n' "$TMP/u1.fq.gz" "$TMP/u2.fq.gz" > "$TMP/list.txt"
mkdir -p "$TMP/kmc_tmp"
kmc -k21 -ci$CI -cs65535 -t16 -m24 -fq @"$TMP/list.txt" "kmers/$id" "$TMP/kmc_tmp" > "logs/kmc_$id.txt" 2>&1 || { tail -5 "logs/kmc_$id.txt"; fail "kmc"; }
[ -s "kmers/$id.kmc_pre" ] || fail "kmc produced no database"
uniq_k=$(grep -m1 "No. of unique k-mers" "logs/kmc_$id.txt" | awk -F: '{gsub(/ /,"",$2); print $2}')
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$id" "$tot" "$sub" "$mapped" "$kept" "$uniq_k" > "stats/$id.tsv"
rm -rf "$TMP"
echo "$id finished $(date +%H:%M:%S): unique 21-mers $uniq_k"
