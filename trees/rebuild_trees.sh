#!/bin/bash
#SBATCH --job-name=rebuild_trees
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/sms049/shared/vselvanayagam/rebuild/rebuild_trees.log

module use -p /group/shooshared/modules
module load kmertools_s

MASH=/group/shooshared/micromamba/envs/kmertools/bin/mash
QUICKTREE=/group/shooshared/micromamba/envs/mashtree/bin/quicktree
PY=/scratch/sms049/shared/vselvanayagam/ml_venv/bin/python3
META=/group/shooshared/kpdata/anigo_filter/status/kpgenomics_sequence_inventory_2026-08-13.csv
WORK=/scratch/sms049/shared/vselvanayagam/rebuild
TMP=${SLURM_TMPDIR:-/tmp}/rebuild_$$
mkdir -p "$TMP"
cd "$WORK" || exit 1

fail() { echo "FAILED: $*"; exit 1; }

find_r1() {
  local d=$1 t=$2 id=$3 f
  for f in "$d/${id}.${t}_1.fastq.gz" "$d/${id}_1.fastq.gz" "$d/${id}_1.sub.fastq.gz" "$d/${id}".*_1.fastq.gz; do
    if [ -s "$f" ]; then echo "$f"; return; fi
  done
}

# streaming subsample: first count the reads, then keep a fixed proportion (low memory)
subsample() {
  local f=$1 o=$2 n=$3 tot p
  tot=$(seqkit stats -T "$f" | awk -F'\t' 'NR==2{print $4}')
  if [ -z "$tot" ]; then echo "could not count reads in $f"; return 1; fi
  if [ "$tot" -le "$n" ]; then cp "$f" "$o"; return 0; fi
  p=$(awk -v n="$n" -v t="$tot" 'BEGIN{printf "%.8f", n/t}')
  seqkit sample -p "$p" -s 42 "$f" -o "$o" 2>> "seqkit_$tag.err"
}

run_comp() {
  local tag=$1 dir=$2 n=$3
  mkdir -p "sketches_$tag"
  : > "missing_$tag.txt"
  while read -r id; do
    if [ -s "sketches_$tag/$id.msh" ]; then continue; fi
    r1=$(find_r1 "$dir" "$tag" "$id")
    if [ -z "$r1" ]; then echo "$id" >> "missing_$tag.txt"; continue; fi
    r2=${r1/_1./_2.}
    ok1=0; ok2=0
    n1=$(seqkit stats -T "$r1" 2>/dev/null | awk -F'\t' 'NR==2{print $4}'); [ -n "$n1" ] && ok1=1
    n2=""; [ -s "$r2" ] && n2=$(seqkit stats -T "$r2" 2>/dev/null | awk -F'\t' 'NR==2{print $4}'); [ -n "$n2" ] && ok2=1
    if [ $ok1 -eq 1 ] && [ $ok2 -eq 1 ]; then
      subsample "$r1" "$TMP/${id}_a.fq.gz" "$n" || { tail -5 "seqkit_$tag.err"; fail "subsample $r1"; }
      subsample "$r2" "$TMP/${id}_b.fq.gz" "$n" || { tail -5 "seqkit_$tag.err"; fail "subsample $r2"; }
      cat "$TMP/${id}_a.fq.gz" "$TMP/${id}_b.fq.gz" > "$TMP/$id.fq.gz"
    elif [ $ok1 -eq 1 ] || [ $ok2 -eq 1 ]; then
      if [ $ok1 -eq 1 ]; then one=$r1; else one=$r2; fi
      echo "$id: one mate is damaged, using $(basename "$one") only" | tee -a "damaged_$tag.txt"
      subsample "$one" "$TMP/$id.fq.gz" "$((2*n))" || { tail -5 "seqkit_$tag.err"; fail "subsample $one"; }
    else
      echo "$id: both mates damaged, sample skipped" | tee -a "damaged_$tag.txt"
      echo "$id" >> "missing_$tag.txt"; continue
    fi
    "$MASH" sketch -r -s 500000 -k 21 -o "sketches_$tag/$id" "$TMP/$id.fq.gz" > "mash_$tag.out" 2>&1 || { cat "mash_$tag.out"; fail "mash sketch on $id"; }
    [ -s "sketches_$tag/$id.msh" ] || { cat "mash_$tag.out"; fail "empty sketch for $id"; }
    rm -f "$TMP/${id}_a.fq.gz" "$TMP/${id}_b.fq.gz" "$TMP/$id.fq.gz"
    echo "$tag $id done $(date +%H:%M:%S)"
  done < ids53.txt

  nsk=$(find "sketches_$tag" -name "*.msh" -size +0 | wc -l)
  echo "$tag: non-empty sketches: $nsk; missing: $(wc -l < missing_$tag.txt)"
  [ "$nsk" -eq 53 ] || fail "$tag has $nsk sketches, expected 53"
  rm -f "all_$tag.msh"
  "$MASH" paste "all_$tag" sketches_$tag/*.msh > /dev/null || fail "mash paste $tag"
  "$MASH" triangle "all_$tag.msh" > "dist_$tag.raw" || fail "mash triangle $tag"
  awk 'BEGIN{FS=OFS="\t"} NR==1{print;next} {match($1,/KP[0-9]+/); $1=substr($1,RSTART,RLENGTH); print}' "dist_$tag.raw" > "dist_$tag.phylip"
  [ -s "dist_$tag.phylip" ] || fail "empty distance matrix for $tag"
  "$QUICKTREE" -in m -out t "dist_$tag.phylip" > "tree_${tag}_unrooted.nwk" || fail "quicktree $tag"
  [ -s "tree_${tag}_unrooted.nwk" ] || fail "empty tree for $tag"
  "$PY" tree_tools.py root "tree_${tag}_unrooted.nwk" "tree_${tag}_pure53_rooted.nwk" KP0012 || fail "rooting $tag"
  "$PY" tree_tools.py summary "tree_${tag}_pure53_rooted.nwk" "$META" > "summary_$tag.txt" || fail "summary $tag"
  echo "$tag tree finished"
}

run_comp cp /group/shooshared/kpdata/anigo_filter/filtered/cp_reads 500000
run_comp mt /group/shooshared/kpdata/anigo_filter/mt_filter/mt_reads 60000

rm -rf "$TMP"
echo "=== chloroplast ==="; cat summary_cp.txt
echo "=== mitochondrial ==="; cat summary_mt.txt
echo "Done!"
