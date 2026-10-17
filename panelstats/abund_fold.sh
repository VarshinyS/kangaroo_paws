#!/bin/bash
module use /group/shooshared/modules
module load kmertools_s
S=/scratch/sms049/shared/vselvanayagam/sig20m/specific
cd /scratch/sms049/shared/vselvanayagam/sig20m/panelstats
for SP in flavidus rufus pulcherrimus rufpul; do
  kmc_tools transform $S/${SP}_abund -ci1 -cx65535 dump -s dump/${SP}_abundset.txt
done
/scratch/sms049/shared/vselvanayagam/ml_venv/bin/python3 - << 'PY'
import numpy as np
tot = {k: float(v) for k, v in (l.split() for l in open("totals.tsv"))}
mt = np.mean(list(tot.values()))
refs = {"flavidus": ["KP0042","KP0043","KP0052","KP1008"],
        "rufus": ["KP1067","KP0058","KP0054"],
        "pulcherrimus": ["KP1066","KP1222","KP1223"],
        "rufpul": ["KP1067","KP0058","KP0054","KP1066","KP1222","KP1223"]}
out = ["set\tabundant\tabund_fold_lt2\tpct\tmedian_CV_abund"]
for sp, rr in refs.items():
    km = np.loadtxt(f"dump/{sp}_{rr[0]}.txt", usecols=0, dtype=str)
    raw = np.vstack([np.loadtxt(f"dump/{sp}_{r}.txt", usecols=1) for r in rr])
    nrm = raw * (mt / np.array([tot[r] for r in rr]))[:, None]
    ab = np.isin(km, np.loadtxt(f"dump/{sp}_abundset.txt", usecols=0, dtype=str))
    fold = nrm.max(0) / nrm.min(0); cv = nrm.std(0, ddof=1) / nrm.mean(0)
    n = ab.sum(); k = (ab & (fold < 2)).sum()
    out.append(f"{sp}\t{n}\t{k}\t{100*k/n:.1f}\t{np.median(cv[ab]):.3f}")
open("abund_fold.tsv", "w").write("\n".join(out) + "\n")
print("\n".join(out))
PY
