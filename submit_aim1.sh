#!/bin/bash
#SBATCH --job-name=aim1_kp
#SBATCH --ntasks=8
#SBATCH --mem=180G
#SBATCH --time=12:00:00
#SBATCH --output=/home/vselvanayagam/varshiny/logs/aim1_%j.out
#SBATCH --error=/home/vselvanayagam/varshiny/logs/aim1_%j.err

export PATH=/group/shooshared/micromamba/envs/skmer/bin:$PATH

/group/shooshared/micromamba/envs/skmer/bin/python3 ~/varshiny/scripts/pilot_kmer.py
