#!/bin/bash
#SBATCH --job-name=cp_f1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=cp_f1.log
set -eo pipefail
module use -p /group/shooshared/modules >/dev/null 2>&1
module load kmertools_s >/dev/null 2>&1
cd /scratch/sms049/shared/vselvanayagam/rebuild
MASH=/group/shooshared/micromamba/envs/kmertools/bin/mash
CP=/group/shooshared/kpdata/anigo_filter/filtered/cp_reads
O=cp_f1; rm -rf $O; mkdir -p $O/tmp
for id in KP0030 KP0029 KP0068 KP0038 KP0065; do
  n=$(seqkit stats -T $CP/$id.cp_1.fastq.gz | awk -F'\t' 'NR==2{print $4}')
  p=$(awk -v n=$n 'BEGIN{x=500000/n; if(x>1)x=1; printf "%.6f", x}')
  seqkit sample -p $p -s 42 $CP/$id.cp_1.fastq.gz -o $O/tmp/${id}_a.fq.gz 2>/dev/null
  seqkit sample -p $p -s 42 $CP/$id.cp_2.fastq.gz -o $O/tmp/${id}_b.fq.gz 2>/dev/null
  cat $O/tmp/${id}_a.fq.gz $O/tmp/${id}_b.fq.gz > $O/tmp/$id.fq.gz
  $MASH sketch -r -s 500000 -k 21 -o $O/$id $O/tmp/$id.fq.gz >/dev/null 2>&1
  $MASH dist all_cp.msh $O/$id.msh > $O/${id}_dist.tsv
  rm -f $O/tmp/${id}*
  echo "$id done $(date +%H:%M:%S)"
done
/scratch/sms049/shared/vselvanayagam/ml_venv/bin/python3 - << 'PYEOF'
import re
lab = {}
for line in open("itol_labels_samples.txt"):
    f = line.strip().split(None, 1)
    if len(f) == 2 and f[0].startswith("KP"):
        lab[f[0]] = f[1]
names = {"KP0030": "Landscape Orange", "KP0029": "Yellow Gem", "KP0068": "Landscape Scarlet",
         "KP0038": "Landscape Gold", "KP0065": "Orange Cross"}
out = open("cp_f1/summary.txt", "w")
def say(m):
    print(m); out.write(m + "\n")
for q, qn in names.items():
    rows = []
    for line in open(f"cp_f1/{q}_dist.tsv"):
        f = line.split("\t")
        m = re.search(r"(KP\d+)", f[0])
        if m:
            rows.append((float(f[2]), m.group(1)))
    rows.sort()
    say(f"\n{q} {qn}: nearest chloroplasts")
    for d, r in rows[:6]:
        say(f"  {d:.5f}  {r}  {lab.get(r, '?')}")
    grp = {}
    for d, r in rows:
        l = lab.get(r, "").lower()
        for sp in ("flavidus", "rufus", "pulcherrimus"):
            if sp in l:
                grp.setdefault(sp, []).append(d)
    means = {sp: sum(v) / len(v) for sp, v in grp.items()}
    say("  mean distance: " + ", ".join(f"{sp} {means[sp]:.5f} (n={len(grp[sp])})" for sp in means))
    say(f"  closest group: {min(means, key=means.get)}")
PYEOF
echo "Done!"
