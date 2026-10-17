#!/bin/bash
#SBATCH --job-name=score_hyb
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=score_hybrids.log
set -eo pipefail
module use -p /group/shooshared/modules >/dev/null 2>&1
module load kmertools_s >/dev/null 2>&1
cd /scratch/sms049/shared/vselvanayagam/sig20m
K="kmc_tools -t16 -hp"
HYB="KP0030 KP0068 KP0038 KP0065 KP0040"
nk(){ kmc_tools info $1 | awk -F: '/total k-mers/{gsub(/ /,"",$2); print $2}'; }
printf "sample\tpairs_kept\tmin_count\n" > robust/thresholds_hyb.tsv
for id in $HYB; do
  [ -s kmers/$id.kmc_pre ] || { echo "missing kmers/$id"; exit 1; }
  kept=$(awk 'END{print $5}' stats/$id.tsv)
  n=$(awk -v k=$kept 'BEGIN{x=int(10*k/16000000+0.5); if(x<2)x=2; print x}')
  printf "%s\t%s\t%s\n" $id $kept $n >> robust/thresholds_hyb.tsv
  $K simple robust/cand kmers/$id intersect robust/tmp_$id -ocright >/dev/null
  $K transform robust/tmp_$id dump robust/tab/$id.txt >/dev/null
  rm -f robust/tmp_$id.kmc_*
done
cat robust/thresholds_hyb.tsv
minc(){ awk -v s=$1 '$1==s{print $3}' robust/thresholds_hyb.tsv; }
printf "set\tlevel\tsize" > specific/tests_hyb.tsv; for t in $HYB; do printf "\t%s" $t >> specific/tests_hyb.tsv; done; echo >> specific/tests_hyb.tsv
for g in flavidus rufus pulcherrimus rufpul; do
  for lev in len str abund; do
    size=$(nk specific/${g}_$lev); printf "%s\t%s\t%s" $g $lev $size >> specific/tests_hyb.tsv
    for t in $HYB; do
      if [ "$size" -gt 0 ]; then
        $K simple specific/${g}_$lev kmers/$t -ci$(minc $t) intersect specific/p_$t -ocleft >/dev/null
        printf "\t%.3f" $(awk -v a=$(nk specific/p_$t) -v b=$size 'BEGIN{print a/b}') >> specific/tests_hyb.tsv
        rm -f specific/p_$t.kmc_*
      else printf "\tNA" >> specific/tests_hyb.tsv; fi
    done; echo >> specific/tests_hyb.tsv
  done
done
echo "Share of each specific set present in each hybrid:"; column -t specific/tests_hyb.tsv
/scratch/sms049/shared/vselvanayagam/ml_venv/bin/python3 score_hybrids.py
echo "Done!"
