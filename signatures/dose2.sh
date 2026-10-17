#!/bin/bash
#SBATCH --job-name=dose2
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=dose2.log
set -eo pipefail
module use -p /group/shooshared/modules >/dev/null 2>&1
module load kmertools_s >/dev/null 2>&1
cd /scratch/sms049/shared/vselvanayagam/sig20m
K="kmc_tools -t16 -hp"
D=specific/dose
NEW="KP0004 KP0013 KP0014 KP0015 KP0017 KP0018 KP0019 KP0020 KP0021 KP0022 KP0024 KP0039 KP0046 KP0064 KP0069 KP1204 KP0025 KP0026 KP0027 KP1004 KP1013 KP1058"
$K simple specific/flavidus_abund specific/rufus_abund union $D/u1 >/dev/null
$K simple $D/u1 specific/pulcherrimus_abund union $D/u2 >/dev/null
$K simple $D/u2 specific/rufpul_abund union $D/u3 >/dev/null
for s in $NEW; do
  [ -s kmers/$s.kmc_pre ] || { echo "no k-mers for $s, skipped"; continue; }
  $K simple $D/u3 kmers/$s intersect $D/x -ocright >/dev/null
  $K transform $D/x dump $D/$s.txt >/dev/null
  rm -f $D/x.kmc_*
  echo "$s done $(date +%H:%M:%S)"
done
rm -f $D/u1.kmc_* $D/u2.kmc_* $D/u3.kmc_*
/scratch/sms049/shared/vselvanayagam/ml_venv/bin/python3 dose2.py
echo "Done!"
