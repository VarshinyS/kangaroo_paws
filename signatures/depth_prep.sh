#!/bin/bash
#SBATCH --job-name=depth_prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/sms049/shared/vselvanayagam/sig20m/depth_prep.log

module use -p /group/shooshared/modules >/dev/null 2>&1
module load kmertools_s >/dev/null 2>&1

WORK=${WORK:-/scratch/sms049/shared/vselvanayagam/sig20m}
NUC=${NUC:-/group/shooshared/kpdata/anigo_filter/mt_filter/nuc_reads}
BBMAP=${BBMAP:-/uwahpc/rocky9/bio-apps/binary/bbmap/39.17/bin/bbmap.sh}
LIB=${LIB:-/group/shooshared/kpdata/anigo_filter/flavidus/repeat/flavidus_repeat_library.classfied}
LEVELS=${LEVELS:-"2 4 6 8 10 12"}
UNIT=${UNIT:-1000000}
THREADS=${THREADS:-16}
cd "$WORK" || exit 1
mkdir -p kmers depth stats logs
fail() { echo "FAILED: $*"; exit 1; }
[ -s "$LIB" ] || fail "repeat library not found"

nreads() { seqkit stats -T "$1" 2>/dev/null | awk -F'\t' 'NR==2{print $4}'; }

subsample() {   # file out n
  local f=$1 o=$2 n=$3 tot p
  tot=$(nreads "$f"); [ -n "$tot" ] || { echo "could not count reads in $f"; return 1; }
  if [ "$tot" -le "$n" ]; then cp "$f" "$o"; return 0; fi
  p=$(awk -v n="$n" -v t="$tot" 'BEGIN{printf "%.8f", n/t}')
  seqkit sample -p "$p" -s 42 "$f" -o "$o" 2>> logs/depth_seqkit.err
}

count() {   # db files...
  local db=$1; shift
  printf '%s\n' "$@" > "$T/list.txt"; mkdir -p "$T/kmc_tmp"
  kmc -k21 -ci1 -cs65535 -t$THREADS -m${KMCMEM:-24} -fq @"$T/list.txt" "$db" "$T/kmc_tmp" > "$T/kmc.txt" 2>&1 || { tail -5 "$T/kmc.txt"; fail "kmc $db"; }
  [ -s "$db.kmc_pre" ] || fail "no database $db"
}

filter_sample() {   # id pairs -> sets U1 (and U2 if paired) in $T
  local id=$1 n=$2
  local r1="$NUC/$id.nuc_1.fastq.gz" r2="$NUC/$id.nuc_2.fastq.gz" n1 n2
  n1=$(nreads "$r1"); n2=$(nreads "$r2")
  if [ -n "$n1" ] && [ -n "$n2" ]; then
    subsample "$r1" "$T/a1.fq.gz" "$n" || fail "subsample $r1"
    subsample "$r2" "$T/a2.fq.gz" "$n" || fail "subsample $r2"
    "$BBMAP" in1="$T/a1.fq.gz" in2="$T/a2.fq.gz" ref="$LIB" minid=0.95 outu1="$T/u1.fq.gz" outu2="$T/u2.fq.gz" \
      nodisk overwrite=t fastareadlen=500 t=$THREADS -Xmx24g 2> "logs/depth_bbmap_$id.txt" || fail "bbmap $id"
    U1="$T/u1.fq.gz"; U2="$T/u2.fq.gz"
  elif [ -n "$n1" ] || [ -n "$n2" ]; then
    local one=$r1; [ -n "$n1" ] || one=$r2
    echo "$id: one mate is damaged, using $(basename "$one") only"
    subsample "$one" "$T/a1.fq.gz" "$((2*n))" || fail "subsample $one"
    "$BBMAP" in="$T/a1.fq.gz" ref="$LIB" minid=0.95 outu="$T/u1.fq.gz" \
      nodisk overwrite=t fastareadlen=500 t=$THREADS -Xmx24g 2> "logs/depth_bbmap_$id.txt" || fail "bbmap $id"
    U1="$T/u1.fq.gz"; U2=""
  else
    fail "both mates of $id unreadable"
  fi
  echo "$id repeat-mapped: $(grep -m1 -E '^mapped:' "logs/depth_bbmap_$id.txt" | awk '{print $2}')"
}

# Part A: KP1008 as a reference, from its intact mate
if [ ! -s kmers/KP1008.kmc_pre ]; then
  T=$WORK/tmp_dp_KP1008; rm -rf "$T"; mkdir -p "$T"
  filter_sample KP1008 20000000
  count kmers/KP1008 $U1 $U2
  echo -e "KP1008\tNA\tNA\tboth-mates\t$(nreads "$U1")\tNA" > stats/KP1008.tsv
  rm -rf "$T"; echo "KP1008 reference done $(date +%H:%M:%S)"
fi

# Part B: KP0055 at several depths of repeat-filtered read pairs
T=$WORK/tmp_dp_KP0055; rm -rf "$T"; mkdir -p "$T"
filter_sample KP0055 20000000
all=$(nreads "$U1")
: > depth/levels.tsv
for L in $LEVELS; do
  n=$((L*UNIT))
  [ "$n" -lt "$all" ] || { echo "skip ${L}M (only $all pairs)"; continue; }
  subsample "$U1" "$T/l1.fq.gz" "$n" || fail "level $L"
  if [ -n "$U2" ]; then subsample "$U2" "$T/l2.fq.gz" "$n" || fail "level $L"; count depth/KP0055_${L}M "$T/l1.fq.gz" "$T/l2.fq.gz"
  else count depth/KP0055_${L}M "$T/l1.fq.gz"; fi
  echo -e "KP0055_${L}M\t$(nreads "$T/l1.fq.gz")" >> depth/levels.tsv
  rm -f "$T/l1.fq.gz" "$T/l2.fq.gz"
  echo "level ${L}M done $(date +%H:%M:%S)"
done
count depth/KP0055_all $U1 $U2
echo -e "KP0055_all\t$all" >> depth/levels.tsv
rm -rf "$T"
cat depth/levels.tsv
echo "Done!"
