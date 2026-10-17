import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HIST = sys.argv[1] if len(sys.argv) > 1 else "hist"
CUT = int(sys.argv[2]) if len(sys.argv) > 2 else 2      # counts up to and including CUT are the low-frequency tail
XMAX = 100
GROUP = {}
for ids, g in ((("KP0042", "KP0043", "KP0052", "KP0055", "KP1008", "KP0003", "KP0001"), "flavidus"),
               (("KP1067", "KP0058", "KP0054", "KP0006"), "rufus"),
               (("KP1066", "KP1222", "KP1223"), "pulcherrimus"),
               (("KP0029", "KP0041"), "hybrid")):
    for i in ids:
        GROUP[i] = g
COL = {"flavidus": "#2E8B57", "rufus": "#FF4500", "pulcherrimus": "#20B2AA", "hybrid": "#555555"}
NAME = {"KP0001": "'Stripie' 1", "KP0042": "'Stripie' 2", "KP0043": "'Stripie' 3", "KP0003": "'green'",
        "KP0052": "'Tall Orange'", "KP0055": "'red'", "KP1008": "F8_002", "KP0006": "'Frosty Yellow'",
        "KP0054": "'Scorcher's Flame'", "KP0058": "'Frosty Red'", "KP1067": "RUF3_006", "KP1066": "PUL6_003",
        "KP1222": "Pul6_006", "KP1223": "Pul8_010", "KP0029": "Yellow Gem", "KP0041": "Stripie x Yellow Gem"}

def load(path, cmax=500):
    h = np.zeros(cmax + 1)
    for line in open(path):
        p = line.split()
        if len(p) == 2 and int(p[0]) <= cmax:
            h[int(p[0])] = float(p[1])
    return h

files = sorted(f for f in os.listdir(HIST) if f.endswith(".txt"))
ids = [f[:-4] for f in files]
data = {i: load(os.path.join(HIST, f"{i}.txt")) for i in ids}

rows = []
n = len(ids)
ncol = 4
nrow = int(np.ceil(n / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(15, 3.1 * nrow))
axes = np.atleast_1d(axes).ravel()
for ax, i in zip(axes, ids):
    h = data[i]
    c = np.arange(len(h))
    mode = CUT + 1 + int(np.argmax(h[CUT + 1:XMAX]))
    ax.plot(c[CUT + 1:XMAX], h[CUT + 1:XMAX], color=COL[GROUP[i]], lw=1.5)
    ax.axvline(mode, color="black", ls="--", lw=0.8)
    ax.set_yscale("log")
    flag = " (at the cut: shallow)" if mode == CUT + 1 else ""
    ax.set_title(f"{NAME.get(i, i)} ({i})\nmost common count above the tail: {mode}{flag}", fontsize=8)
    ax.set_xlabel("k-mer count", fontsize=7); ax.tick_params(labelsize=7)
    tot = h[1:].sum()
    tail = h[1:CUT + 1].sum()
    occ = (c * h)[1:]
    occ_total = occ.sum()
    hi = 3 * mode
    rows.append((i, GROUP[i], int(tot), round(100 * tail / tot, 1), int(tot - tail), mode,
                 round(100 * occ[hi - 1:].sum() / occ_total, 1)))
for ax in axes[n:]:
    ax.axis("off")
fig.suptitle(f"K-mer count spectra of the repeat-filtered nuclear reads with counts 1 to {CUT} removed", fontsize=11)
plt.tight_layout(rect=(0, 0, 1, 0.97))
plt.savefig("kmer_spectra_tail_cut.png", dpi=200)

fig, ax = plt.subplots(figsize=(8, 5.5))
for i in ids:
    h = data[i]
    c = np.arange(len(h))
    y = h / h[CUT + 1:XMAX].max()
    ax.plot(c[CUT + 1:XMAX], y[CUT + 1:XMAX], color=COL[GROUP[i]], lw=1, alpha=0.8, label=GROUP[i])
hnd, lab = ax.get_legend_handles_labels()
seen = dict(zip(lab, hnd))
ax.legend(seen.values(), seen.keys(), fontsize=8)
ax.set_yscale("log"); ax.set_xlabel("k-mer count"); ax.set_ylabel("distinct k-mers (relative to the highest point)")
ax.set_title(f"Spectra with counts 1 to {CUT} removed, all samples")
plt.tight_layout(); plt.savefig("kmer_spectra_overlay_tail_cut.png", dpi=200)

with open("kmer_copynumber_table.tsv", "w") as f:
    f.write("sample\tgroup\tdistinct_kmers\tpct_kmers_in_tail\tdistinct_kmers_above_tail\tmost_common_count_above_tail\tpct_kmer_occurrences_at_3x_or_more_of_that_count\n")
    for r in rows:
        f.write("\t".join(str(x) for x in r) + "\n")
print(open("kmer_copynumber_table.tsv").read())
