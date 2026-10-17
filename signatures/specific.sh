#!/bin/bash
#SBATCH --job-name=specific
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=specific.log
set -eo pipefail
module use -p /group/shooshared/modules >/dev/null 2>&1
module load kmertools_s >/dev/null 2>&1
cd /scratch/sms049/shared/vselvanayagam/sig20m
K="kmc_tools -t16 -hp"
O=specific; rm -rf $O; mkdir -p $O/tmp
FLAV="KP0042 KP0043 KP0052 KP1008"
RUF="KP1067 KP0058 KP0054"
PUL="KP1066 KP1222 KP1223"
TESTS="KP0055 KP0001 KP0003 KP0006 KP0029 KP0041"
minc(){ awk -v s=$1 '$1==s{print $3}' robust/thresholds.tsv; }
nk(){ kmc_tools info $1 | awk -F: '/total k-mers/{gsub(/ /,"",$2); print $2}'; }

chain(){  # chain <op> <out> <first> <rest...>
  local op=$1 out=$2 cur=$3 i=0; shift 3
  local oc=""; if [ "$op" = intersect ]; then oc="-ocleft"; fi
  for x in "$@"; do
    i=$((i+1)); $K simple $cur $x $op $O/tmp/$(basename $out)_$i $oc >/dev/null
    cur=$O/tmp/$(basename $out)_$i
  done
  $K transform $cur reduce $out >/dev/null
  rm -f $O/tmp/$(basename $out)_*
}
paths(){ local d=$1; shift; for x in "$@"; do echo -n "$d/$x "; done; }

printf "set\trecurrent\tspecific_lenient\tspecific_strict\tabundant_strict\n" > $O/counts.tsv
for g in flavidus rufus pulcherrimus rufpul; do
  case $g in
    flavidus)     IN="$FLAV";     OUT="$RUF $PUL" ;;
    rufus)        IN="$RUF";      OUT="$FLAV $PUL" ;;
    pulcherrimus) IN="$PUL";      OUT="$FLAV $RUF" ;;
    rufpul)       IN="$RUF $PUL"; OUT="$FLAV" ;;
  esac
  chain intersect $O/${g}_rec $(paths robust/thr $IN)
  chain kmers_subtract $O/${g}_len $O/${g}_rec $(paths robust/thr $OUT)
  chain kmers_subtract $O/${g}_str $O/${g}_len $(paths kmers $OUT)
  $K transform $O/${g}_str -ci20 reduce $O/${g}_abund >/dev/null
  printf "%s\t%s\t%s\t%s\t%s\n" $g $(nk $O/${g}_rec) $(nk $O/${g}_len) $(nk $O/${g}_str) $(nk $O/${g}_abund) >> $O/counts.tsv
  echo "$g done $(date +%H:%M:%S)"
done
column -t $O/counts.tsv

printf "set\tlevel\tsize" > $O/tests.tsv; for t in $TESTS; do printf "\t%s" $t >> $O/tests.tsv; done; echo >> $O/tests.tsv
for g in flavidus rufus pulcherrimus rufpul; do
  for lev in len str abund; do
    size=$(nk $O/${g}_$lev)
    printf "%s\t%s\t%s" $g $lev $size >> $O/tests.tsv
    for t in $TESTS; do
      if [ "$size" -gt 0 ]; then
        $K simple $O/${g}_$lev kmers/$t -ci$(minc $t) intersect $O/tmp/p -ocleft >/dev/null
        printf "\t%.3f" $(awk -v a=$(nk $O/tmp/p) -v b=$size 'BEGIN{print a/b}') >> $O/tests.tsv
        rm -f $O/tmp/p.kmc_*
      else
        printf "\tNA" >> $O/tests.tsv
      fi
    done
    echo >> $O/tests.tsv
  done
done
echo "Share of each specific set present in each test plant:"
column -t $O/tests.tsv
rm -rf $O/tmp
echo "Done!"
