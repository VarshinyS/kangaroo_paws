#!/bin/bash
#SBATCH --job-name=calib
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=calib.log
set -eo pipefail
module use -p /group/shooshared/modules >/dev/null 2>&1
module load kmertools_s >/dev/null 2>&1
cd /scratch/sms049/shared/vselvanayagam/sig20m
K="kmc_tools -t16 -hp"
REFS="KP0042 KP0043 KP0052 KP1008 KP1067 KP0058 KP0054 KP1066 KP1222 KP1223"
TESTS="KP0055 KP0001 KP0003 KP0006 KP0029 KP0041"
HYB="KP0030 KP0068 KP0038 KP0065 KP0040"
MIX="MIX_FR10 MIX_FR25 MIX_FR50 MIX_FR75 MIX_FP10 MIX_FP25 MIX_FP50 MIX_FP75"
D=specific/dose; rm -rf $D; mkdir -p $D
for g in flavidus rufus pulcherrimus rufpul; do
  $K transform specific/${g}_abund dump $D/set_$g.txt >/dev/null
done
$K simple specific/flavidus_abund specific/rufus_abund union $D/u1 >/dev/null
$K simple $D/u1 specific/pulcherrimus_abund union $D/u2 >/dev/null
$K simple $D/u2 specific/rufpul_abund union $D/u3 >/dev/null
for s in $REFS $TESTS $HYB $MIX; do
  [ -s kmers/$s.kmc_pre ] || { echo "missing kmers/$s, skipped"; continue; }
  $K simple $D/u3 kmers/$s intersect $D/x -ocright >/dev/null
  $K transform $D/x dump $D/$s.txt >/dev/null
  rm -f $D/x.kmc_*
  echo "$s done $(date +%H:%M:%S)"
done
rm -f $D/u1.kmc_* $D/u2.kmc_* $D/u3.kmc_*
/scratch/sms049/shared/vselvanayagam/ml_venv/bin/python3 calib.py
echo "Done!"
