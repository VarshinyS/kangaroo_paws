#!/bin/bash
#SBATCH --job-name=panelstats
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=panelstats_%j.log

module use /group/shooshared/modules
module load kmertools_s

W=/scratch/sms049/shared/vselvanayagam/sig20m
S=$W/specific; K=$W/kmers; O=$W/panelstats
PY=/scratch/sms049/shared/vselvanayagam/ml_venv/bin/python3
cd $O; mkdir -p dump tmp

declare -A REFS
REFS[flavidus]="KP0042 KP0043 KP0052 KP1008"
REFS[rufus]="KP1067 KP0058 KP0054"
REFS[pulcherrimus]="KP1066 KP1222 KP1223"
REFS[rufpul]="KP1067 KP0058 KP0054 KP1066 KP1222 KP1223"

> totals.tsv
for R in KP0042 KP0043 KP0052 KP1008 KP1067 KP0058 KP0054 KP1066 KP1222 KP1223; do
  kmc_tools transform $K/$R -ci1 -cx65535 histogram tmp/$R.hist -ci1 -cx65535
  awk -v r=$R '{t+=$1*$2} END{print r"\t"t}' tmp/$R.hist >> totals.tsv
done

for SP in flavidus rufus pulcherrimus rufpul; do
  for R in ${REFS[$SP]}; do
    kmc_tools simple $S/${SP}_str -ci1 $K/$R -ci1 intersect tmp/x -ocright -cs65535
    kmc_tools transform tmp/x -ci1 -cx65535 dump -s dump/${SP}_${R}.txt
    rm -f tmp/x.kmc_pre tmp/x.kmc_suf
  done
done

$PY $O/panel_stats.py
echo done
