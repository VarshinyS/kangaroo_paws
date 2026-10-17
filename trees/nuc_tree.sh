#!/bin/bash
#SBATCH --job-name=nuc_tree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/sms049/shared/vselvanayagam/rebuild/nuc_tree.log

MASH=/group/shooshared/micromamba/envs/kmertools/bin/mash
QUICKTREE=/group/shooshared/micromamba/envs/mashtree/bin/quicktree
PY=/scratch/sms049/shared/vselvanayagam/ml_venv/bin/python3
META=/group/shooshared/kpdata/anigo_filter/status/kpgenomics_sequence_inventory_2026-08-13.csv
cd /scratch/sms049/shared/vselvanayagam/rebuild || exit 1
fail() { echo "FAILED: $*"; exit 1; }

nsk=$(find sketches_nuc -name "*.msh" -size +0 | wc -l)
echo "nuclear sketches: $nsk"
[ "$nsk" -ge 50 ] || fail "only $nsk sketches (need at least 50)"
[ -s sketches_nuc/KP0012.msh ] || fail "the root sample KP0012 has no sketch"

"$MASH" paste all_nuc sketches_nuc/*.msh > /dev/null || fail "mash paste"
"$MASH" triangle all_nuc.msh > dist_nuc.raw || fail "mash triangle"
awk 'BEGIN{FS=OFS="\t"} NR==1{print;next} {match($1,/KP[0-9]+/); $1=substr($1,RSTART,RLENGTH); print}' dist_nuc.raw > dist_nuc.phylip
[ -s dist_nuc.phylip ] || fail "empty distance matrix"
"$QUICKTREE" -in m -out t dist_nuc.phylip > tree_nuc_unrooted.nwk || fail "quicktree"
[ -s tree_nuc_unrooted.nwk ] || fail "empty tree"
"$PY" tree_tools.py root tree_nuc_unrooted.nwk tree_nuc_pure53_rooted.nwk KP0012 || fail "rooting"
"$PY" tree_tools.py summary tree_nuc_pure53_rooted.nwk "$META" > summary_nuc.txt || fail "summary"
echo "=== nuclear ==="; cat summary_nuc.txt
echo "Done!"
