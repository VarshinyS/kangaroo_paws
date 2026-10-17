import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import f as fdist

O = "robust"
CV_MAX, FC_MIN, PSEUDO, ALPHA = 0.25, 1.0, 0.05, 0.01
SPECIES = {"flavidus": ["KP0042", "KP0043", "KP0052", "KP1008"],
           "rufus": ["KP1067", "KP0058", "KP0054"],
           "pulcherrimus": ["KP1066", "KP1222", "KP1223"]}
TESTS = ["KP0055", "KP0001", "KP0003", "KP0006", "KP0029", "KP0041"]
HYB = ["MIX_FR10", "MIX_FR25", "MIX_FR50", "MIX_FR75", "MIX_FP10", "MIX_FP25", "MIX_FP50", "MIX_FP75"]
NAME = {"KP0055": "'red' (flavidus)", "KP0001": "'Stripie' 1 (flavidus)", "KP0003": "'green' (flavidus)",
        "KP0006": "'Frosty Yellow' (rufus)", "KP0029": "'Yellow Gem'", "KP0041": "'Stripie' x 'Yellow Gem'",
        "KP0030": "'Landscape Orange'", "KP0068": "'Landscape Scarlet'", "KP0038": "'Landscape Gold'",
        "KP0065": "'Orange Cross'", "KP0040": "'Early Spring' x 'Stripie'"}
SP = list(SPECIES)
REFS = [s for sp in SP for s in SPECIES[sp]]
ALL = REFS + TESTS + HYB
col = {s: j for j, s in enumerate(ALL)}

thr = {}
for fn in ["thresholds.tsv", "thresholds_mix.tsv"]:
    for line in open(f"{O}/{fn}"):
        f = line.split()
        if f[0] != "sample":
            thr[f[0]] = float(f[1])
idx = {}
with open(f"{O}/cand.txt") as fh:
    for line in fh:
        idx[line.split()[0]] = len(idx)
M = len(idx)
C = np.zeros((M, len(ALL)), dtype=np.int64)
for j, s in enumerate(ALL):
    with open(f"{O}/tab/{s}.txt") as fh:
        for line in fh:
            k, c = line.split()
            C[idx[k], j] = int(c)
    print(f"loaded {s}", flush=True)

kept = np.array([thr[s] for s in ALL])
F = C / kept * 1e6
L = np.log2(F + PSEUDO)
stable_all = np.ones(M, bool)
for sp in SP:
    cols = [col[s] for s in SPECIES[sp]]
    m = F[:, cols].mean(1); sd = F[:, cols].std(1, ddof=1)
    stable_all &= np.where(m > 0, sd / np.maximum(m, 1e-12), 0.0) <= CV_MAX
gm = np.stack([L[:, [col[s] for s in SPECIES[sp]]].mean(1) for sp in SP], 1)
ng = np.array([len(SPECIES[sp]) for sp in SP]); N = ng.sum(); kg = 3
grand = L[:, [col[s] for s in REFS]].mean(1)
ssb = (ng * (gm - grand[:, None]) ** 2).sum(1)
ssw = sum(((L[:, [col[s] for s in SPECIES[sp]]] - gm[:, [i]]) ** 2).sum(1) for i, sp in enumerate(SP))
Fr = (ssb / (kg - 1)) / np.maximum(ssw / (N - kg), 1e-9)
fc = np.max(np.abs(gm[:, :, None] - gm[:, None, :]), axis=(1, 2))
p = fdist.sf(Fr, kg - 1, N - kg)
q = np.ones(M); s_idx = np.where(stable_all)[0]; ps = p[s_idx]; o = np.argsort(ps); n = len(ps)
qs = np.minimum.accumulate((ps[o] * n / np.arange(1, n + 1))[::-1])[::-1]
qq = np.empty(n); qq[o] = np.minimum(qs, 1); q[s_idx] = qq
distinct = stable_all & (q <= ALPHA) & (fc >= FC_MIN)
cent = gm[distinct]
print(f"distinct k-mers rebuilt: {distinct.sum():,}", flush=True)

shares, rs = {}, {}
for line in open(f"{O}/mixFR_tests.tsv"):
    f = line.rstrip("\n").split("\t")
    if f[0] == "mix":
        shares.setdefault(float(f[1]), []).append(float(f[9]))
xs = sorted(shares); rm = [np.mean(shares[x]) for x in xs]

out = open(f"{O}/mixreads_summary.txt", "w")
def say(m):
    print(m, flush=True); out.write(m + "\n")
say("Frequency correlation over distinct k-mers, and the rufus share of the simulated flavidus + rufus mix with the same flavidus correlation")
names, eqs = [], []
for s in TESTS + HYB:
    x = L[distinct, col[s]]
    r = [np.corrcoef(x, cent[:, i])[0, 1] for i in range(3)]
    eq = float(np.interp(r[0], rm[::-1], xs[::-1]))
    names.append(NAME.get(s, s)); eqs.append(100 * eq)
    say(f"{s}\t{NAME.get(s, s):28s}\tr flavidus {r[0]:.3f}  rufus {r[1]:.3f}  pulcherrimus {r[2]:.3f}\tmatching mix: {100*eq:.0f}% rufus")

fig, ax = plt.subplots(figsize=(7, 5))
y = np.arange(len(names))[::-1]
colors = ["#2e8b57"] * 3 + ["#ff4500"] + ["0.4"] * 2 + ["0.65"] * len(HYB)
ax.barh(y, eqs, color=colors)
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8)
ax.axvline(50, color="0.7", lw=1, ls="--")
ax.set_xlim(0, 100); ax.set_xlabel("rufus share of the matching simulated flavidus + rufus mix (%)")
fig.tight_layout(); fig.savefig(f"{O}/fig_mixreads.png", dpi=200)
say("Outputs: mixreads_summary.txt, fig_mixreads.png, ../specific/tests_hyb.tsv")
