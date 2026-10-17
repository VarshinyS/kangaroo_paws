import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch

COL = {"flavidus": "#2e8b57", "rufus": "#ff4500", "humilis": "#e3b505", "bicolor": "#8a6fd1",
       "manglesii": "#2c7fb8", "viridis": "#17becf", "unknown": "#c8c8c8"}
BUSH = [
    ("Bush Coral", 58.9, "unknown", None, ""),
    ("Bush Fury", 58.0, "putative", {"flavidus": 50, "rufus": 50}, ""),
    ("Bush Ballad", 57.7, "unknown", None, ""),
    ("Bush Inferno", 56.1, "partial", {"flavidus": 50, "unknown": 50}, ""),
    ("Bush Pearl", 52.7, "known", {"flavidus": 50, "humilis": 25, "bicolor": 25}, "a"),
    ("Bush Diamond", 52.7, "known", {"flavidus": 50, "humilis": 25, "bicolor": 25}, "b"),
    ("Bush Blitz", 50.2, "unknown", None, ""),
    ("Bush Glow", 37.6, "conflict", None, "c"),
    ("Bush Fire", 33.1, "unknown", None, ""),
    ("Bush Flare", 25.3, "unknown", None, ""),
    ("Bush Bonanza", 24.2, "conflict", None, "c"),
    ("Bush Tenacity", 14.3, "unknown", None, ""),
    ("Bush Zest", 13.4, "unknown", None, ""),
    ("Bush Elegance", 13.3, "conflict", None, "c"),
]
CELEB = [
    ("Kings Park Royale", 4.5, "known", {"manglesii": 50, "viridis": 50}, ""),
    ("Aussie Spirit", 4.4, "known", {"manglesii": 50, "viridis": 50}, ""),
    ("Masquerade", 4.5, "putative", {"manglesii": 1, "viridis": 1, "humilis": 1}, "d"),
]

def draw(ax, name, meas, kind, comp, mark):
    ax.set_aspect("equal"); ax.axis("off")
    if kind in ("known", "putative", "partial"):
        keys = list(comp)
        wedges, _ = ax.pie([comp[k] for k in keys], colors=[COL[k] for k in keys],
                           startangle=90, counterclock=False, radius=1,
                           wedgeprops=dict(edgecolor="white", linewidth=1.2))
        if kind == "putative":
            for w in wedges:
                w.set_hatch("///"); w.set_alpha(0.75)
        if kind == "partial":
            ax.text(-0.45, -0.35, "?", ha="center", va="center", fontsize=18, color="#555555")
    else:
        ax.add_patch(Circle((0, 0), 1, color=COL["unknown"]))
        ax.text(0, -0.05, "?", ha="center", va="center", fontsize=30, color="white", weight="bold")
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.9, 1.5)
    ax.text(0, 1.3, f"'{name}'" + (f" ({mark})" if mark else ""), ha="center", va="center", fontsize=9)
    ax.text(0, -1.35, f"measured {meas:.0f}%", ha="center", va="center", fontsize=9,
            color=COL["flavidus"], weight="bold")
    ax.text(0, -1.7, "A. flavidus", ha="center", va="center", fontsize=8,
            color=COL["flavidus"], style="italic")

fig, axes = plt.subplots(3, 7, figsize=(13, 7.6))
for ax, row in zip(axes.flat[:14], BUSH):
    draw(ax, *row)
for ax in axes[2]:
    ax.axis("off")
for ax, row in zip(axes[2][2:5], CELEB):
    draw(ax, *row)
fig.text(0.02, 0.965, "Bush series", fontsize=10, weight="bold")
fig.text(0.02, 0.385, "Kings Park Celebration range", fontsize=10, weight="bold")
handles = [Patch(color=COL[k], label=f"A. {k}") for k in ["flavidus", "rufus", "humilis", "bicolor", "manglesii", "viridis"]]
handles += [Patch(color=COL["unknown"], label="unknown / undisclosed"),
            Patch(facecolor="white", edgecolor="0.4", hatch="///", label="tentative record")]
fig.legend(handles=handles, loc="lower center", ncol=8, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, 0.0))
fig.text(0.5, 0.05, "(a) mean of two samples (53.2 and 52.2%); (b) sport of 'Bush Pearl'; "
         "(c) records conflicting or cultivar identity uncertain; (d) species reported, proportions unknown",
         ha="center", fontsize=8, color="0.3")
plt.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.10, wspace=0.05, hspace=0.12)
fig.savefig("specific/fig5B_bush_parentage.png", dpi=200)
print("saved specific/fig5B_bush_parentage.png")
