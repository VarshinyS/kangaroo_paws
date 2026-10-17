import re, csv, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

rx = re.compile(r"^(KP\d{4})\s+(.+?)\s+(-?[\d.]+)%\s+(-?[\d.]+)%\s+.*\)\s+([\d.]+)\s+([\d.]+)\s*$")
D = {}
for line in open("flavidus_share.txt"):
    m = rx.match(line.strip())
    if m:
        k, n, f, g, r, p = m.groups()
        D[k] = dict(name=n.strip("'"), f=float(f), g=float(g), r=float(r), p=float(p))
H = {r["sample_id"]: (float(r["height_min_m"]), float(r["height_max_m"])) for r in csv.DictReader(open("heights.csv"))}

A = ["KP0029","KP0030","KP0038","KP0065","KP0068","KP0040","KP0041"]
bush = ["KP0004","KP0013","KP0014","KP0015","KP0017","KP0018","KP0019","KP0020","KP0021","KP0022","KP0024","KP0039","KP0046","KP0064","KP0069"]
B = sorted(bush, key=lambda k: -D[k]["f"]) + ["KP0025","KP0026","KP0027"]
bg = np.mean([D[k]["g"] for k in ["KP1004","KP1013","KP1058"]]) / 100
print(f"short-species group background: {bg*100:.1f}%")

cF, cR, cP, cU, cH = "#2a9d3c", "#d95f02", "#f2a541", "#cfcfcf", "#e8b27d"
fig, ax = plt.subplots(2, 2, figsize=(14, 7.4), sharex="col",
                       gridspec_kw={"width_ratios": [7, 18], "height_ratios": [3, 1.1], "hspace": 0.08})

def label(k):
    n = D[k]["name"]
    return n.replace("Bush ", "") if k in bush else n

def heights(a, ids):
    for i, k in enumerate(ids):
        if k in H:
            lo, hi = H[k]
            a.plot([i, i], [lo, hi], color="#444", lw=5, solid_capstyle="butt")
            a.plot(i, (lo + hi) / 2, "o", color="black", ms=3)
        else:
            a.text(i, 0.25, "n.d.", ha="center", fontsize=7, color="grey")
    a.set_ylim(0, 2.2); a.set_yticks([0, 1, 2])

a = ax[0, 0]
for i, k in enumerate(A):
    d = D[k]; f = min(max(d["f"], 0), 100); g = min(max(d["g"], 0), 100 - f)
    col = cR if d["r"] > d["p"] else cP
    a.bar(i, f, color=cF, width=0.85); a.bar(i, g, bottom=f, color=col, width=0.85)
    a.bar(i, 100 - f - g, bottom=f + g, color=cU, width=0.85)
a.set_ylim(0, 100); a.set_ylabel("Calibrated ancestry (%)")
a.set_title("A", loc="left", fontweight="bold", fontsize=14)
heights(ax[1, 0], A); ax[1, 0].set_ylabel("Height (m)")
ax[1, 0].set_xticks(range(len(A))); ax[1, 0].set_xticklabels([label(k) for k in A], rotation=45, ha="right", fontsize=8)

b = ax[0, 1]
for i, k in enumerate(B):
    d = D[k]; f = min(max(d["f"], 0), 100); g = min(d["g"], 100 - f)
    b.bar(i, f, color=cF, width=0.85)
    b.bar(i, g, bottom=f, color=cH, width=0.85, hatch="//", edgecolor="white", linewidth=0)
    b.bar(i, 100 - f - g, bottom=f + g, color=cU, width=0.85)
    e = bg * (100 - f); b.plot([i - 0.42, i + 0.42], [f + e, f + e], color="black", lw=1.4)
    print(f"{k} {label(k):20s} flavidus {d['f']:5.1f} group {d['g']:5.1f} expected {e:5.1f} excess {d['g']-e:+5.1f}")
nb = len(bush) - 0.5
for x in (ax[0, 1], ax[1, 1]): x.axvline(nb, color="grey", ls=":", lw=1)
b.text(nb / 2, 103, "Bush series", ha="center", fontsize=9)
b.text(nb + 1.5, 103, "Celebration range", ha="center", fontsize=9)
b.set_ylim(0, 108); b.set_yticks(range(0, 101, 20))
b.set_title("B", loc="left", fontweight="bold", fontsize=14)
heights(ax[1, 1], B)
ax[1, 1].set_xticks(range(len(B))); ax[1, 1].set_xticklabels([label(k) for k in B], rotation=45, ha="right", fontsize=8)

leg = [Patch(color=cF, label="A. flavidus"), Patch(color=cP, label="A. pulcherrimus (group set)"),
       Patch(color=cR, label="A. rufus (group set)"),
       Patch(facecolor=cH, hatch="//", edgecolor="white", label="Group set reading (includes short-species background)"),
       Line2D([0], [0], color="black", lw=1.4, label="Expected group reading from short species alone"),
       Patch(color=cU, label="Unassigned / other species"),
       Line2D([0], [0], color="#444", lw=5, label="Reported mature height range")]
fig.legend(handles=leg, loc="lower center", ncol=4, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.0))
for row in ax:
    for x in row:
        for s in ["top", "right"]: x.spines[s].set_visible(False)
plt.subplots_adjust(bottom=0.2, wspace=0.12)
plt.savefig("fig6_ancestry_height.png", dpi=300, bbox_inches="tight")
print("saved fig6_ancestry_height.png")
