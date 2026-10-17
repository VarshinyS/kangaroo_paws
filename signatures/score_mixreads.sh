#!/bin/bash
#SBATCH --job-name=score_mix
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=score_mixreads.log
set -eo pipefail
module use -p /group/shooshared/modules >/dev/null 2>&1
module load kmertools_s >/dev/null 2>&1
cd /scratch/sms049/shared/vselvanayagam/sig20m
K="kmc_tools -t16 -hp"
HYB="MIX_FR10 MIX_FR25 MIX_FR50 MIX_FR75 MIX_FP10 MIX_FP25 MIX_FP50 MIX_FP75"
nk(){ kmc_tools info $1 | awk -F: '/total k-mers/{gsub(/ /,"",$2); print $2}'; }
printf "sample\tpairs_kept\tmin_count\n" > robust/thresholds_mix.tsv
for id in $HYB; do
  [ -s kmers/$id.kmc_pre ] || { echo "missing kmers/$id"; exit 1; }
  kept=$(awk 'END{print $5}' stats/$id.tsv)
  n=$(awk -v k=$kept 'BEGIN{x=int(10*k/16000000+0.5); if(x<2)x=2; print x}')
  printf "%s\t%s\t%s\n" $id $kept $n >> robust/thresholds_mix.tsv
  $K simple robust/cand kmers/$id intersect robust/tmp_$id -ocright >/dev/null
  $K transform robust/tmp_$id dump robust/tab/$id.txt >/dev/null
  rm -f robust/tmp_$id.kmc_*
done
cat robust/thresholds_mix.tsv
minc(){ awk -v s=$1 '$1==s{print $3}' robust/thresholds_mix.tsv; }
printf "set\tlevel\tsize" > specific/tests_mix.tsv; for t in $HYB; do printf "\t%s" $t >> specific/tests_mix.tsv; done; echo >> specific/tests_mix.tsv
for g in flavidus rufus pulcherrimus rufpul; do
  for lev in len str abund; do
    size=$(nk specific/${g}_$lev); printf "%s\t%s\t%s" $g $lev $size >> specific/tests_mix.tsv
    for t in $HYB; do
      if [ "$size" -gt 0 ]; then
        $K simple specific/${g}_$lev kmers/$t -ci$(minc $t) intersect specific/p_$t -ocleft >/dev/null
        printf "\t%.3f" $(awk -v a=$(nk specific/p_$t) -v b=$size 'BEGIN{print a/b}') >> specific/tests_mix.tsv
        rm -f specific/p_$t.kmc_*
      else printf "\tNA" >> specific/tests_mix.tsv; fi
    done; echo >> specific/tests_mix.tsv
  done
done
echo "Share of each specific set present in each hybrid:"; column -t specific/tests_mix.tsv
/scratch/sms049/shared/vselvanayagam/ml_venv/bin/python3 score_mixreads.py
echo "Done!"
