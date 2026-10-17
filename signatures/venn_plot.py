import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
c = {}
for r in csv.DictReader(open("specific/counts.tsv"), delimiter="\t"):
    c[r["set"]] = (int(r["specific_lenient"]), int(r["specific_strict"]))
for r in csv.DictReader(open("specific/venn_extra.tsv"), delimiter="\t"):
    c[r["region"]] = (int(r["lenient"]), int(r["strict"]))
def f(n):
    return f"{n/1e6:.2f}M" if n >= 1e6 else (f"{n/1e3:.0f}k" if n >= 1e3 else str(n))
def lab(k):
    l, s = c[k]
    return f(l) if k == "FRP" else f"{f(l)}\n(strict {f(s)})"
fig, ax = plt.subplots(figsize=(7, 6.5))
for (x, y), col in [((0, 0.6), "#2e8b57"), ((-0.5, -0.25), "#ff4500"), ((0.5, -0.25), "#20b2aa")]:
    ax.add_patch(Circle((x, y), 0.9, alpha=0.25, color=col))
pos = {"flavidus": (0, 1.15), "rufus": (-0.95, -0.6), "pulcherrimus": (0.95, -0.6),
       "FR": (-0.55, 0.45), "FP": (0.55, 0.45), "rufpul": (0, -0.6), "FRP": (0, 0.05)}
for k, (x, y) in pos.items():
    ax.text(x, y, lab(k), ha="center", va="center", fontsize=9)
for (x, y), name in [((0, 1.65), "A. flavidus"), ((-1.25, -1.25), "A. rufus"), ((1.25, -1.25), "A. pulcherrimus")]:
    ax.text(x, y, name, ha="center", va="center", fontsize=11, style="italic", weight="bold")
ax.set_xlim(-1.8, 1.8); ax.set_ylim(-1.5, 1.85); ax.set_aspect("equal"); ax.axis("off")
fig.tight_layout(); fig.savefig("specific/fig_venn.png", dpi=200)
print("saved specific/fig_venn.png")
