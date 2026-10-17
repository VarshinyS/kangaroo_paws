import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

O = sys.argv[1] if len(sys.argv) > 1 else "robust"
CV_MAX, FC_MIN, PSEUDO, ALPHA = 0.25, 1.0, 0.05, 0.01
SPECIES = {"flavidus": ["KP0042", "KP0043", "KP0052", "KP1008"],
           "rufus": ["KP1067", "KP0058", "KP0054"],
           "pulcherrimus": ["KP1066", "KP1222", "KP1223"]}
TESTS = ["KP0055", "KP0001", "KP0003", "KP0006", "KP0029", "KP0041"]
LABEL = {"KP0055": "'red' (flavidus, held out)", "KP0001": "'Stripie' 1 (flavidus, shallow)",
         "KP0003": "'green' (flavidus, shallow)", "KP0006": "'Frosty Yellow' (rufus, shallow)",
         "KP0029": "Yellow Gem (hybrid)", "KP0041": "Stripie x Yellow Gem (cross)"}
SP = list(SPECIES)
REFS = [s for sp in SP for s in SPECIES[sp]]
SAMPLES = REFS + TESTS
col = {s: j for j, s in enumerate(SAMPLES)}

out = open(f"{O}/summary.txt", "w")
def say(*a):
    m = " ".join(str(x) for x in a); print(m, flush=True); out.write(m + "\n")

thr = {}
for line in open(f"{O}/thresholds.tsv"):
    f = line.split()
    if f[0] != "sample":
        thr[f[0]] = (float(f[1]), int(f[2]))

idx = {}
with open(f"{O}/cand.txt") as fh:
    for line in fh:
        idx[line.split()[0]] = len(idx)
kmers = list(idx)
M = len(idx)
C = np.zeros((M, len(SAMPLES)), dtype=np.float32)
for j, s in enumerate(SAMPLES):
    with open(f"{O}/tab/{s}.txt") as fh:
        for line in fh:
            k, c = line.split()
            C[idx[k], j] = float(c)
    print(f"loaded {s}", flush=True)

kept = np.array([thr[s][0] for s in SAMPLES])
nmin = np.array([thr[s][1] for s in SAMPLES])
F = C / kept * 1e6            # frequency per million read pairs
P = C >= nmin                 # present above the depth-scaled cut-off
L = np.log2(F + PSEUDO)

top = {}
for sp in SP:
    t = np.zeros(M, bool)
    with open(f"{O}/sp/{sp}_top.txt") as fh:
        for line in fh:
            t[idx[line.split()[0]]] = True
    top[sp] = t

say(f"Candidate k-mers (union of the three signatures): {M:,}")
say("Depth-scaled cut-offs (count >= min_count):")
for s in SAMPLES:
    say(f"  {s}\t{int(thr[s][0]):,} pairs\tmin count {thr[s][1]}")

# A. signatures and overlap (presence/absence)
say("\n== A. Signatures and overlap ==")
excl = {}
for sp in SP:
    others = np.zeros(M, bool)
    for o in SP:
        if o != sp:
            others |= top[o]
    excl[sp] = top[sp] & ~others
    say(f"{sp}: {top[sp].sum():,} k-mers; {excl[sp].sum():,} ({100*excl[sp].sum()/top[sp].sum():.1f}%) not in another species' signature")
for a in range(3):
    for b in range(a + 1, 3):
        i = (top[SP[a]] & top[SP[b]]).sum(); u = (top[SP[a]] | top[SP[b]]).sum()
        say(f"Jaccard {SP[a]}-{SP[b]}: {i/u:.3f} ({i:,} shared)")
say("Share of each signature present in each species' references (mean over references; own species = 1 by construction):")
with open(f"{O}/overlap.tsv", "w") as fo:
    fo.write("signature\t" + "\t".join(SP) + "\n")
    for sp in SP:
        row = [np.mean([P[top[sp], col[s]].mean() for s in SPECIES[o]]) for o in SP]
        fo.write(sp + "\t" + "\t".join(f"{x:.3f}" for x in row) + "\n")
        say(f"  {sp} signature: " + ", ".join(f"{o} {x:.3f}" for o, x in zip(SP, row)))

# B. stability within species
say("\n== B. Within-species stability (CV of frequency across references) ==")
mu, cv, stable = {}, {}, {}
for sp in SP:
    cols = [col[s] for s in SPECIES[sp]]
    m = F[:, cols].mean(1); sd = F[:, cols].std(1, ddof=1)
    mu[sp] = m
    cv[sp] = np.where(m > 0, sd / np.maximum(m, 1e-12), 0.0)
    stable[sp] = cv[sp] <= CV_MAX
    own = top[sp]
    say(f"{sp}: median CV in own signature {np.median(cv[sp][own]):.3f}; stable (CV <= {CV_MAX}) {stable[sp][own].sum():,} of {own.sum():,} ({100*stable[sp][own].mean():.1f}%)")
stable_all = stable[SP[0]] & stable[SP[1]] & stable[SP[2]]
say(f"Stable in all three species: {stable_all.sum():,} of {M:,}")

# C. distinct frequencies between species
say("\n== C. Distinct frequencies between species (stable k-mers only) ==")
gm = np.stack([L[:, [col[s] for s in SPECIES[sp]]].mean(1) for sp in SP], 1)
ng = np.array([len(SPECIES[sp]) for sp in SP]); N = ng.sum(); kgrp = len(SP)
grand = L[:, [col[s] for s in REFS]].mean(1)
ssb = (ng * (gm - grand[:, None]) ** 2).sum(1)
ssw = sum(((L[:, [col[s] for s in SPECIES[sp]]] - gm[:, [i]]) ** 2).sum(1) for i, sp in enumerate(SP))
Fr = (ssb / (kgrp - 1)) / np.maximum(ssw / (N - kgrp), 1e-9)
fc = np.max(np.abs(gm[:, :, None] - gm[:, None, :]), axis=(1, 2))
try:
    from scipy.stats import f as fdist
    p = fdist.sf(Fr, kgrp - 1, N - kgrp)
    q = np.ones(M)
    s_idx = np.where(stable_all)[0]
    ps = p[s_idx]; o = np.argsort(ps); n = len(ps)
    qs = ps[o] * n / np.arange(1, n + 1)
    qs = np.minimum.accumulate(qs[::-1])[::-1]
    qq = np.empty(n); qq[o] = np.minimum(qs, 1); q[s_idx] = qq
    sig = q <= ALPHA
    say(f"Test: one-way ANOVA on log2 frequency, Benjamini-Hochberg FDR <= {ALPHA}, fold change >= {2**FC_MIN:.0f}x")
except ImportError:
    q = np.full(M, np.nan)
    sig = Fr >= 9.55   # F critical value, df (2, 7), alpha 0.01
    say("Test: F-ratio >= 9.55 (df 2,7; alpha 0.01, no FDR correction: scipy missing), fold change >= 2x")
distinct = stable_all & sig & (fc >= FC_MIN)
high = gm.argmax(1)
say(f"Distinct k-mers: {distinct.sum():,}")
for i, sp in enumerate(SP):
    say(f"  highest in {sp}: {(distinct & (high == i)).sum():,}")
order = np.where(distinct)[0][np.argsort(-Fr[distinct])]
with open(f"{O}/distinct_kmers.tsv", "w") as fo:
    fo.write("kmer\thighest_in\t" + "\t".join(f"freq_{sp}" for sp in SP) + "\t" +
             "\t".join(f"cv_{sp}" for sp in SP) + "\tF_ratio\tq\n")
    for r in order[:5000]:
        fo.write(f"{kmers[r]}\t{SP[high[r]]}\t" + "\t".join(f"{mu[sp][r]:.4f}" for sp in SP) + "\t" +
                 "\t".join(f"{cv[sp][r]:.3f}" for sp in SP) + f"\t{Fr[r]:.2f}\t{q[r]:.2e}\n")

# D. held-out test plants
say("\n== D. Held-out test plants (not used to choose any k-mer) ==")
cent = gm[distinct]
with open(f"{O}/tests.tsv", "w") as fo:
    fo.write("sample\tlabel\t" + "\t".join(f"sig_{sp}" for sp in SP) + "\t" +
             "\t".join(f"excl_{sp}" for sp in SP) + "\t" + "\t".join(f"r_{sp}" for sp in SP) + "\tnearest\n")
    for s in TESTS:
        j = col[s]
        sh = [P[top[sp], j].mean() for sp in SP]
        ex = [P[excl[sp], j].mean() if excl[sp].sum() else np.nan for sp in SP]
        if distinct.sum() >= 10:
            x = L[distinct, j]
            r = [np.corrcoef(x, cent[:, i])[0, 1] for i in range(3)]
            d = [np.sqrt(np.mean((x - cent[:, i]) ** 2)) for i in range(3)]
            near = SP[int(np.argmin(d))]
        else:
            r = [np.nan] * 3; near = "NA"
        fo.write(f"{s}\t{LABEL[s]}\t" + "\t".join(f"{v:.3f}" for v in sh + ex + r) + f"\t{near}\n")
        say(f"{s} {LABEL[s]}")
        say("   signature present: " + ", ".join(f"{sp} {v:.3f}" for sp, v in zip(SP, sh)))
        say("   species-only part present: " + ", ".join(f"{sp} {v:.3f}" for sp, v in zip(SP, ex)))
        say("   frequency correlation: " + ", ".join(f"{sp} {v:.3f}" for sp, v in zip(SP, r)) + f" | nearest profile: {near}")

# Figures
colors = {"flavidus": "#2e8b57", "rufus": "#ff4500", "pulcherrimus": "#20b2aa"}
fig, ax = plt.subplots(figsize=(6, 4))
for sp in SP:
    v = cv[sp][top[sp] & (mu[sp] > 0)]
    ax.hist(np.clip(v, 0, 1.5), bins=75, histtype="step", lw=1.5, color=colors[sp], label=sp, density=True)
ax.axvline(CV_MAX, ls="--", c="k", lw=1)
ax.set_xlabel("CV of normalised frequency across references"); ax.set_ylabel("density")
ax.set_title("Within-species stability of signature k-mers"); ax.legend()
fig.tight_layout(); fig.savefig(f"{O}/fig_cv.png", dpi=200)

fig, ax = plt.subplots(figsize=(8, 4))
x = np.arange(len(TESTS)); w = 0.27
for i, sp in enumerate(SP):
    ax.bar(x + (i - 1) * w, [P[top[sp], col[s]].mean() for s in TESTS], w, color=colors[sp], label=sp)
ax.set_xticks(x); ax.set_xticklabels([LABEL[s] for s in TESTS], rotation=30, ha="right", fontsize=8)
ax.set_ylim(0, 1); ax.set_ylabel("share of signature present")
ax.set_title("Species signatures in held-out plants"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{O}/fig_tests.png", dpi=200)

if distinct.sum() > 0:
    sel = []
    for i in range(3):
        ii = np.where(distinct & (high == i))[0]
        sel += list(ii[np.argsort(-Fr[ii])][:100])
    Z = L[sel][:, :]
    Z = (Z - Z.mean(1, keepdims=True)) / np.maximum(Z.std(1, keepdims=True), 1e-9)
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(Z, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5, interpolation="nearest")
    ax.set_xticks(range(len(SAMPLES))); ax.set_xticklabels(SAMPLES, rotation=90, fontsize=7)
    for b in np.cumsum([len(SPECIES[sp]) for sp in SP]):
        ax.axvline(b - 0.5, c="k", lw=1)
    ax.set_yticks([]); ax.set_ylabel("distinct k-mers (up to 100 per species)")
    ax.set_title("Normalised frequency of distinct k-mers (row z-score)")
    fig.colorbar(im, ax=ax, fraction=0.03); fig.tight_layout(); fig.savefig(f"{O}/fig_distinct_heatmap.png", dpi=200)
say("\nOutputs: summary.txt, overlap.tsv, distinct_kmers.tsv, tests.tsv, fig_cv.png, fig_tests.png, fig_distinct_heatmap.png")
