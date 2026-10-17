#!/bin/bash
#SBATCH --job-name=spectra2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=spectra2_%j.log

module use /group/shooshared/modules
module load kmertools_s
cd /scratch/sms049/shared/vselvanayagam/sig20m/spectra
mkdir -p tmp
for N in KP0052_cp KP0052_mt; do
  kmc -k21 -ci1 -cs65535 -m24 -t8 @$N.list $N tmp
  kmc_tools transform $N histogram hist/$N.txt -ci1 -cx65535
  rm -f $N.kmc_pre $N.kmc_suf
done
rm -rf tmp
/scratch/sms049/shared/vselvanayagam/ml_venv/bin/python3 plot_spectra.py
echo done
