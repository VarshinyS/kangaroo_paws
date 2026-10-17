import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
D = "specific/dose"
FLAV = ["KP0042", "KP0043", "KP0052", "KP1008"]
OTHER = {"KP1004": "A. bicolor", "KP1013": "A. humilis", "KP1058": "A. manglesii",
         "KP0025": "'Kings Park Royale'", "KP0026": "'Aussie Spirit'", "KP0027": "'Masquerade'"}
def kept(s):
    return float([l for l in open(f"stats/{s}.tsv") if l.strip()][-1].split()[4])
fk = [l.split()[0] for l in open(f"{D}/set_flavidus.txt")]
idx = {k: i for i, k in enumerate(fk)}
N = len(fk)
def counts(s):
    c = np.zeros(N)
    for line in open(f"{D}/{s}.txt"):
        k, v = line.split()
        i = idx.get(k)
        if i is not None:
            c[i] = int(v)
    return c
F = np.mean([counts(s) / kept(s) * 1e6 for s in FLAV], axis=0)
rank = F.argsort().argsort() / (N - 1) * 100
out = open("specific/overlap_summary.txt", "w")
def say(m):
    print(m, flush=True); out.write(m + "\n")
say(f"Flavidus-specific abundant k-mers: {N:,}; median frequency in flavidus refs {np.median(F):.3f} per million pairs")
pres = {}
for s, name in OTHER.items():
    c = counts(s); k_ = kept(s); n = max(2, int(10 * k_ / 16e6 + 0.5))
    p = c >= n; pres[s] = p
    any1 = c >= 1
    say(f"\n{s} {name}: present above cut-off (>= {n}) {p.sum():,} ({100*p.mean():.2f}%); seen at least once {any1.sum():,} ({100*any1.mean():.1f}%)")
    if p.sum():
        say(f"  overlapping k-mers: median flavidus percentile {np.median(rank[p]):.0f} (all = 50); "
            f"{100*np.mean(rank[p] >= 90):.0f}% in the top 10% most abundant flavidus-specific k-mers")
        say(f"  their count in {name}: median {np.median(c[p]):.0f} (cut-off {n}); "
            f"median ratio to flavidus frequency {np.median((c[p]/k_*1e6)/F[p]):.2f}")
sp = ["KP1004", "KP1013", "KP1058"]
m = np.sum([pres[s] for s in sp], axis=0)
say(f"\nAcross the three short species: in all three {np.sum(m==3):,}, in two {np.sum(m==2):,}, in one {np.sum(m==1):,}, in none {np.sum(m==0):,}")
if np.sum(m == 3):
    say(f"  k-mers in all three: median flavidus percentile {np.median(rank[m==3]):.0f}")
fig, ax = plt.subplots(figsize=(6.5, 4))
bins = np.logspace(np.log10(max(F[F > 0].min(), 1e-3)), np.log10(F.max()), 60)
ax.hist(F, bins=bins, histtype="step", color="#2e8b57", lw=1.5, label="all flavidus-specific", density=True)
for s, col in zip(sp, ["#6a6a6a", "#c0504d", "#e0a030"]):
    ax.hist(F[pres[s]], bins=bins, histtype="step", color=col, lw=1.2, label=f"also in {OTHER[s]}", density=True)
ax.set_xscale("log"); ax.set_xlabel("frequency in A. flavidus references (per million read pairs)")
ax.set_ylabel("density"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig("specific/fig_overlap.png", dpi=200)
say("\nOutputs: specific/overlap_summary.txt, specific/fig_overlap.png")
