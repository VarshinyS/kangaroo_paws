#!/bin/bash
#SBATCH --job-name=spectra
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=spectra_%j.log

module use -p /group/shooshared/modules
module load kmertools_s

D=/group/shooshared/kpdata
W=/scratch/sms049/shared/vselvanayagam/sig20m
O=$W/spectra
cd $O
mkdir -p tmp hist

run () {  # name, read1, read2
  echo -e "$2\n$3" > $1.list
  kmc -k21 -ci1 -cs65535 -m24 -t8 @$1.list $1 tmp
  kmc_tools transform $1 histogram hist/$1.txt -ci1 -cx1000
  rm -f $1.kmc_pre $1.kmc_suf
}

run KP0052_cp  $D/anigo_filter/filtered/cp_reads/KP0052.cp_1.fastq.gz  $D/anigo_filter/filtered/cp_reads/KP0052.cp_2.fastq.gz
run KP0052_mt  $D/anigo_filter/mt_filter/mt_reads/KP0052.mt_1.fastq.gz $D/anigo_filter/mt_filter/mt_reads/KP0052.mt_2.fastq.gz
run KP0042_nuc $D/anigo_filter/mt_filter/nuc_reads/KP0042.nuc_1.fastq.gz $D/anigo_filter/mt_filter/nuc_reads/KP0042.nuc_2.fastq.gz
run KP0042_rp  $D/rp_filter_95/rp_reads/KP0042.rp_1.fastq.gz $D/rp_filter_95/rp_reads/KP0042.rp_2.fastq.gz
run KP0042_sc  $D/rp_filter_95/sc_reads/KP0042.sc_1.fastq.gz $D/rp_filter_95/sc_reads/KP0042.sc_2.fastq.gz

for L in 2M 4M 6M 8M 10M 12M all; do
  kmc_tools transform $W/depth/KP0055_$L histogram hist/KP0055_$L.txt -ci1 -cx1000
done

rm -rf tmp
echo done
