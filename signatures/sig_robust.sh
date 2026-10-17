#!/bin/bash
#SBATCH --job-name=sig_robust
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=sig_robust.log
set -eo pipefail
module use -p /group/shooshared/modules >/dev/null 2>&1
module load kmertools_s >/dev/null 2>&1
cd /scratch/sms049/shared/vselvanayagam/sig20m
PY=/scratch/sms049/shared/vselvanayagam/ml_venv/bin/python3
BASE_COUNT=10; BASE_PAIRS=16000000; TOPN=1000000
FLAV="KP0042 KP0043 KP0052 KP1008"
RUF="KP1067 KP0058 KP0054"
PUL="KP1066 KP1222 KP1223"
TESTS="KP0055 KP0001 KP0003 KP0006 KP0029 KP0041"
K="kmc_tools -t16 -hp"
O=robust; rm -rf $O; mkdir -p $O/thr $O/sp $O/tab $O/tmp

echo "waiting for kmers/KP1008 $(date +%H:%M:%S)"
while [ ! -s kmers/KP1008.kmc_pre ] || [ ! -s kmers/KP1008.kmc_suf ] || \
      [ $(( $(date +%s) - $(stat -c %Y kmers/KP1008.kmc_suf) )) -lt 180 ]; do sleep 60; done
echo "KP1008 ready $(date +%H:%M:%S)"

printf "sample\tpairs_kept\tmin_count\n" > $O/thresholds.tsv
for id in $FLAV $RUF $PUL $TESTS; do
  [ -s kmers/$id.kmc_pre ] || { echo "missing kmers/$id"; exit 1; }
  kept=$(awk 'END{print $5}' stats/$id.tsv)
  [[ "$kept" =~ ^[0-9]+$ ]] || { echo "no pairs-kept number for $id: $kept"; exit 1; }
  n=$(awk -v k=$kept -v b=$BASE_COUNT -v p=$BASE_PAIRS 'BEGIN{x=int(b*k/p+0.5); if(x<2)x=2; print x}')
  printf "%s\t%s\t%s\n" $id $kept $n >> $O/thresholds.tsv
done
cat $O/thresholds.tsv
minc(){ awk -v s=$1 '$1==s{print $3}' $O/thresholds.tsv; }

for id in $FLAV $RUF $PUL; do
  $K transform kmers/$id -ci$(minc $id) reduce $O/thr/$id >/dev/null
done
echo "cut-offs applied $(date +%H:%M:%S)"

signature(){
  local sp=$1; shift
  local cur=$O/thr/$1; shift
  local i=0 T
  for id in "$@"; do
    i=$((i+1))
    $K simple $cur $O/thr/$id intersect $O/tmp/${sp}_rec$i -ocmin >/dev/null
    cur=$O/tmp/${sp}_rec$i
  done
  $K transform $cur histogram $O/sp/${sp}_hist.txt -ci1 -cx100000 >/dev/null
  T=$(sort -k1,1nr $O/sp/${sp}_hist.txt | awk -v n=$TOPN '$2>0{c+=$2; if(c>=n){print $1; exit}}') || true
  T=${T:-1}
  $K transform $cur -ci$T reduce $O/sp/${sp}_top >/dev/null
  $K transform $O/sp/${sp}_top dump $O/sp/${sp}_top.txt >/dev/null
  echo "$sp: recurrent $(kmc_tools info $cur | grep -i 'total k-mers'), top cut-off >= $T, signature $(wc -l < $O/sp/${sp}_top.txt) $(date +%H:%M:%S)"
}
signature flavidus $FLAV
signature rufus $RUF
signature pulcherrimus $PUL

$K simple $O/sp/flavidus_top $O/sp/rufus_top union $O/tmp/c1 -ocmax >/dev/null
$K simple $O/tmp/c1 $O/sp/pulcherrimus_top union $O/cand -ocmax >/dev/null
$K transform $O/cand dump $O/cand.txt >/dev/null
echo "candidates: $(wc -l < $O/cand.txt)"

for id in $FLAV $RUF $PUL $TESTS; do
  $K simple $O/cand kmers/$id intersect $O/tmp/$id -ocright >/dev/null
  $K transform $O/tmp/$id dump $O/tab/$id.txt >/dev/null
  rm -f $O/tmp/$id.kmc_*
done
echo "counts extracted $(date +%H:%M:%S)"

$PY robust_analysis.py $O
rm -rf $O/tmp
echo "Done!"
