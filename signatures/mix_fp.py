import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

O = sys.argv[1] if len(sys.argv) > 1 else "robust"
CV_MAX, FC_MIN, PSEUDO, ALPHA = 0.25, 1.0, 0.05, 0.01
BASE_COUNT, BASE_PAIRS = 10, 16e6
SPECIES = {"flavidus": ["KP0042", "KP0043", "KP0052", "KP1008"],
           "rufus": ["KP1067", "KP0058", "KP0054"],
           "pulcherrimus": ["KP1066", "KP1222", "KP1223"]}
TESTS = ["KP0055", "KP0001", "KP0003", "KP0006", "KP0029", "KP0041"]
SP = list(SPECIES)
REFS = [s for sp in SP for s in SPECIES[sp]]
SAMPLES = REFS + TESTS
col = {s: j for j, s in enumerate(SAMPLES)}
SRC_A, SRC_B, TARGET = "KP0055", "KP1066", "KP0041"
FRACS = [0.0, 0.10, 0.25, 0.50, 0.75, 1.0]
REPS = 3

thr = {}
for line in open(f"{O}/thresholds.tsv"):
    f = line.split()
    if f[0] != "sample":
        thr[f[0]] = float(f[1])
idx = {}
with open(f"{O}/cand.txt") as fh:
    for line in fh:
        idx[line.split()[0]] = len(idx)
M = len(idx)
C = np.zeros((M, len(SAMPLES)), dtype=np.int64)
for j, s in enumerate(SAMPLES):
    with open(f"{O}/tab/{s}.txt") as fh:
        for line in fh:
            k, c = line.split()
            C[idx[k], j] = int(c)
    print(f"loaded {s}", flush=True)
top = {}
for sp in SP:
    t = np.zeros(M, bool)
    with open(f"{O}/sp/{sp}_top.txt") as fh:
        for line in fh:
            t[idx[line.split()[0]]] = True
    top[sp] = t
excl = {}
for sp in SP:
    others = np.zeros(M, bool)
    for o in SP:
        if o != sp:
            others |= top[o]
    excl[sp] = top[sp] & ~others

kept = np.array([thr[s] for s in SAMPLES])
F = C / kept * 1e6
L = np.log2(F + PSEUDO)
stable_all = np.ones(M, bool)
for sp in SP:
    cols = [col[s] for s in SPECIES[sp]]
    m = F[:, cols].mean(1)
    sd = F[:, cols].std(1, ddof=1)
    cv = np.where(m > 0, sd / np.maximum(m, 1e-12), 0.0)
    stable_all &= cv <= CV_MAX
gm = np.stack([L[:, [col[s] for s in SPECIES[sp]]].mean(1) for sp in SP], 1)
ng = np.array([len(SPECIES[sp]) for sp in SP])
N = ng.sum()
kg = len(SP)
grand = L[:, [col[s] for s in REFS]].mean(1)
ssb = (ng * (gm - grand[:, None]) ** 2).sum(1)
ssw = sum(((L[:, [col[s] for s in SPECIES[sp]]] - gm[:, [i]]) ** 2).sum(1)
          for i, sp in enumerate(SP))
Fr = (ssb / (kg - 1)) / np.maximum(ssw / (N - kg), 1e-9)
fc = np.max(np.abs(gm[:, :, None] - gm[:, None, :]), axis=(1, 2))
from scipy.stats import f as fdist
p = fdist.sf(Fr, kg - 1, N - kg)
q = np.ones(M)
s_idx = np.where(stable_all)[0]
ps = p[s_idx]
o = np.argsort(ps)
n = len(ps)
qs = np.minimum.accumulate((ps[o] * n / np.arange(1, n + 1))[::-1])[::-1]
qq = np.empty(n)
qq[o] = np.minimum(qs, 1)
q[s_idx] = qq
distinct = stable_all & (q <= ALPHA) & (fc >= FC_MIN)
cent = gm[distinct]
print(f"distinct k-mers rebuilt: {distinct.sum():,}", flush=True)


def score(counts, pairs):
    nmin = max(2, int(BASE_COUNT * pairs / BASE_PAIRS + 0.5))
    l = np.log2(counts / pairs * 1e6 + PSEUDO)
    pres = counts >= nmin
    sig = [pres[top[sp]].mean() for sp in SP]
    ex = [pres[excl[sp]].mean() for sp in SP]
    x = l[distinct]
    r = [np.corrcoef(x, cent[:, i])[0, 1] for i in range(3)]
    d = [np.sqrt(np.mean((x - cent[:, i]) ** 2)) for i in range(3)]
    return sig, ex, r, SP[int(np.argmin(d))]


out = open(f"{O}/mixFP_summary.txt", "w")


def say(m):
    print(m, flush=True)
    out.write(m + "\n")


head = ("sample\tyg_share\trep\t" + "\t".join(f"sig_{s}" for s in SP) + "\t" +
        "\t".join(f"excl_{s}" for s in SP) + "\t" + "\t".join(f"r_{s}" for s in SP) + "\tnearest")
rows = [head]


def row(name, fr, rep, res):
    sig, ex, r, near = res
    rows.append(f"{name}\t{fr}\t{rep}\t" + "\t".join(f"{v:.3f}" for v in sig + ex + r) + f"\t{near}")
    say(f"{name:8s} YG {str(fr):5s} rep {rep}: "
        f"sig F/R/P {sig[0]:.3f}/{sig[1]:.3f}/{sig[2]:.3f} | "
        f"species-only F/R/P {ex[0]:.3f}/{ex[1]:.3f}/{ex[2]:.3f} | "
        f"r F/R/P {r[0]:.3f}/{r[1]:.3f}/{r[2]:.3f} | {near}")


Nt = thr[TARGET]
nA = thr[SRC_A]
nB = thr[SRC_B]
say(f"Mixes of {SRC_A} ('red', flavidus) and {SRC_B} (pulcherrimus), "
    f"thinned to {int(Nt):,} pairs (depth of {TARGET})")
for s in [SRC_A, SRC_B, TARGET, "KP0001"]:
    row(s, "real", "-", score(C[:, col[s]], thr[s]))

rng = np.random.default_rng(42)
A = C[:, col[SRC_A]]
B = C[:, col[SRC_B]]
res_r, res_ex = {}, {}
for fr in FRACS:
    pa = (1 - fr) * Nt / nA
    pb = fr * Nt / nB
    assert pa <= 1 and pb <= 1
    for rep in range(1, REPS + 1):
        mix = rng.binomial(A, pa) + rng.binomial(B, pb)
        res = score(mix, Nt)
        row("mix", fr, rep, res)
        res_r.setdefault(fr, []).append(res[2][0])
        res_ex.setdefault(fr, []).append(res[1][1])
open(f"{O}/mixFP_tests.tsv", "w").write("\n".join(rows) + "\n")

cross = score(C[:, col[TARGET]], Nt)
strip = score(C[:, col["KP0001"]], thr["KP0001"])
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
panels = [(axes[0], res_r, "correlation with flavidus profile", cross[2][0], strip[2][0]),
          (axes[1], res_ex, "share of rufus-only k-mers present", cross[1][1], strip[1][1])]
for ax, res, yl, cval, sval in panels:
    xs = [100 * f for f in FRACS for _ in res[f]]
    ys = [v for f in FRACS for v in res[f]]
    ax.plot(xs, ys, "o", color="0.3", label="simulated mixes")
    ax.plot([100 * f for f in FRACS], [np.mean(res[f]) for f in FRACS], "-", color="0.3")
    ax.axhline(cval, color="k", ls="--", label="cross KP0041")
    ax.axhline(sval, color="#2e8b57", ls=":", label="'Stripie' 1 KP0001")
    ax.axvline(50, color="0.7", lw=1)
    ax.set_xlabel("pulcherrimus share in mix (%)")
    ax.set_ylabel(yl)
axes[0].legend(fontsize=8)
fig.suptitle("Simulated flavidus + pulcherrimus mixes compared with the recorded cross")
fig.tight_layout()
fig.savefig(f"{O}/fig_mixFP.png", dpi=200)
say("Outputs: mixFP_summary.txt, mixFP_tests.tsv, fig_mixFP.png")
