import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

tot = dict(l.split() for l in open("totals.tsv"))
tot = {k: float(v) for k, v in tot.items()}
mt = np.mean(list(tot.values()))
refs = {"flavidus": ["KP0042","KP0043","KP0052","KP1008"],
        "rufus": ["KP1067","KP0058","KP0054"],
        "pulcherrimus": ["KP1066","KP1222","KP1223"],
        "rufpul": ["KP1067","KP0058","KP0054","KP1066","KP1222","KP1223"]}

out = ["set\tstrict\tabundant_min20\tfold_lt2\tabundant_and_fold_lt2\tmedian_CV\tmedian_fold"]
fig, ax = plt.subplots(2, 4, figsize=(16, 7))
for j, sp in enumerate(refs):
    raw = np.vstack([np.loadtxt(f"dump/{sp}_{r}.txt", usecols=1, dtype=np.int64) for r in refs[sp]])
    nrm = raw * (mt / np.array([tot[r] for r in refs[sp]]))[:, None]
    ab = raw.min(0) >= 20
    mean = nrm.mean(0)
    cv = nrm.std(0, ddof=1) / mean
    fold = nrm.max(0) / nrm.min(0)
    low = fold < 2
    out.append(f"{sp}\t{raw.shape[1]}\t{ab.sum()}\t{low.sum()}\t{(ab & low).sum()}\t{np.median(cv):.3f}\t{np.median(fold):.2f}")
    ax[0, j].hist(mean, bins=np.logspace(0, 4, 60), color="grey")
    ax[0, j].set_xscale("log"); ax[0, j].axvline(20, color="red", ls="--")
    ax[0, j].set_title(sp); ax[0, j].set_xlabel("mean depth-scaled count")
    ax[1, j].hist(np.clip(fold, 1, 10), bins=60, color="grey")
    ax[1, j].axvline(2, color="red", ls="--"); ax[1, j].set_xlabel("fold range across references")
ax[0, 0].set_ylabel("specific k-mers"); ax[1, 0].set_ylabel("specific k-mers")
plt.tight_layout(); plt.savefig("fig_panel_stats.png", dpi=300)
open("panel_stats.tsv", "w").write("\n".join(out) + "\n")
print("\n".join(out))
