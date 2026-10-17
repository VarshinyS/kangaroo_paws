#!/bin/bash
#SBATCH --job-name=venn
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=venn.log
set -eo pipefail
module use -p /group/shooshared/modules >/dev/null 2>&1
module load kmertools_s >/dev/null 2>&1
cd /scratch/sms049/shared/vselvanayagam/sig20m
K="kmc_tools -t16 -hp"
O=specific; mkdir -p $O/tmp
FLAV="KP0042 KP0043 KP0052 KP1008"
RUF="KP1067 KP0058 KP0054"
PUL="KP1066 KP1222 KP1223"
nk(){ kmc_tools info $1 | awk -F: '/total k-mers/{gsub(/ /,"",$2); print $2}'; }
chain(){
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
printf "region\tlenient\tstrict\n" > $O/venn_extra.tsv
chain intersect $O/FR_rec $(paths robust/thr $FLAV $RUF)
chain kmers_subtract $O/FR_len $O/FR_rec $(paths robust/thr $PUL)
chain kmers_subtract $O/FR_str $O/FR_len $(paths kmers $PUL)
printf "FR\t%s\t%s\n" $(nk $O/FR_len) $(nk $O/FR_str) >> $O/venn_extra.tsv
chain intersect $O/FP_rec $(paths robust/thr $FLAV $PUL)
chain kmers_subtract $O/FP_len $O/FP_rec $(paths robust/thr $RUF)
chain kmers_subtract $O/FP_str $O/FP_len $(paths kmers $RUF)
printf "FP\t%s\t%s\n" $(nk $O/FP_len) $(nk $O/FP_str) >> $O/venn_extra.tsv
chain intersect $O/FRP_rec $(paths robust/thr $FLAV $RUF $PUL)
n=$(nk $O/FRP_rec); printf "FRP\t%s\t%s\n" $n $n >> $O/venn_extra.tsv
column -t $O/venn_extra.tsv
/scratch/sms049/shared/vselvanayagam/ml_venv/bin/python3 venn_plot.py
rm -rf $O/tmp
echo "Done!"
