import os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "specific/dose"
SETS = ["flavidus", "rufus", "pulcherrimus", "rufpul"]
SPECIES = {"flavidus": ["KP0042", "KP0043", "KP0052", "KP1008"],
           "rufus": ["KP1067", "KP0058", "KP0054"],
           "pulcherrimus": ["KP1066", "KP1222", "KP1223"]}
MIX = ["MIX_FR10", "MIX_FR25", "MIX_FR50", "MIX_FR75", "MIX_FP10", "MIX_FP25", "MIX_FP50", "MIX_FP75"]
NAME = {"KP0055": "'red'", "KP0001": "'Stripie' 1", "KP0003": "'green'",
        "KP0029": "'Yellow Gem'", "KP0030": "'Landscape Orange'", "KP0068": "'Landscape Scarlet'",
        "KP0038": "'Landscape Gold'", "KP0065": "'Orange Cross'",
        "KP0040": "'Early Spring' x 'Stripie'", "KP0041": "'Stripie' x 'Yellow Gem'",
        "KP0004": "'Bush Tenacity'", "KP0013": "'Bush Ballad'", "KP0014": "'Bush Blitz'",
        "KP0015": "'Bush Bonanza'", "KP0017": "'Bush Diamond'", "KP0018": "'Bush Fire'",
        "KP0019": "'Bush Fury'", "KP0020": "'Bush Glow'", "KP0021": "'Bush Inferno'",
        "KP0022": "'Bush Pearl' 1", "KP0024": "'Bush Zest'", "KP0039": "'Bush Elegance'",
        "KP0046": "'Bush Coral'", "KP0064": "'Bush Flare'", "KP0069": "'Bush Pearl' 2",
        "KP1204": "'Bush Pearl' 3", "KP0025": "'Kings Park Royale'", "KP0026": "'Aussie Spirit'",
        "KP0027": "'Masquerade'", "KP1004": "A. bicolor (BD6_015)", "KP1013": "A. humilis (H7_009)",
        "KP1058": "A. manglesii (M13_003)"}
GROUPS = [("Pure flavidus (test)", ["KP0055", "KP0001", "KP0003"], "#2e8b57"),
          ("F1 tall hybrids", ["KP0029", "KP0030", "KP0068", "KP0038", "KP0065"], "#7a9e5a"),
          ("Recorded crosses", ["KP0040", "KP0041"], "#a3a3a3"),
          ("Bush series", ["KP0004", "KP0013", "KP0014", "KP0015", "KP0017", "KP0018", "KP0019", "KP0020",
                           "KP0021", "KP0022", "KP0024", "KP0039", "KP0046", "KP0064", "KP0069", "KP1204"], "#e0a030"),
          ("Thought to lack flavidus", ["KP0025", "KP0026", "KP0027"], "#c0504d"),
          ("Other species (interference)", ["KP1004", "KP1013", "KP1058"], "#6a6a6a")]

def kept(s):
    lines = [l for l in open(f"stats/{s}.tsv") if l.strip()]
    return float(lines[-1].split()[4])

member = {}
for i, g in enumerate(SETS):
    for line in open(f"{D}/set_{g}.txt"):
        member[line.split()[0]] = i
size = np.bincount(list(member.values()), minlength=4).astype(float)

dose, pres = {}, {}
for fn in glob.glob(f"{D}/*.txt"):
    s = os.path.basename(fn)[:-4]
    if s.startswith("set_") or not os.path.exists(f"stats/{s}.tsv"):
        continue
    k_ = kept(s); n = max(2, int(10 * k_ / 16e6 + 0.5))
    tot = np.zeros(4); p = np.zeros(4)
    for line in open(fn):
        k, c = line.split(); i = member.get(k)
        if i is None:
            continue
        c = int(c); tot[i] += c
        if c >= n:
            p[i] += 1
    dose[s] = tot / k_ * 1e6; pres[s] = p / size

ref = np.array([np.mean([dose[s][0] for s in SPECIES["flavidus"]]),
                np.mean([dose[s][1] for s in SPECIES["rufus"]]),
                np.mean([dose[s][2] for s in SPECIES["pulcherrimus"]]),
                np.mean([dose[s][3] for s in SPECIES["rufus"] + SPECIES["pulcherrimus"]])])
est = {s: dose[s] / ref for s in dose}
mixes = [m for m in MIX if m in est]
tF = np.array([1 - int(m[-2:]) / 100 for m in mixes])
aF, bF = np.polyfit(tF, [est[m][0] for m in mixes], 1)
aG, bG = np.polyfit(1 - tF, [est[m][3] for m in mixes], 1)

out = open("specific/flavidus_share.txt", "w")
def say(m):
    print(m, flush=True); out.write(m + "\n")
say(f"Calibration (from {len(mixes)} read mixes): flavidus est = {aF:.3f} x true + {bF:.3f}; group est = {aG:.3f} x true + {bG:.3f}")
say("sample\tname\tcalibrated % flavidus\tcalibrated % rufus+pulcherrimus\tflavidus-specific k-mers present\trufus est\tpulcherrimus est")
names, vals, cols = [], [], []
for gname, ids, colr in GROUPS:
    say(f"\n{gname}")
    for s in ids:
        if s not in est:
            say(f"{s}\t{NAME[s]}\tnot available"); continue
        e = est[s]
        f = 100 * (e[0] - bF) / aF; g = 100 * (e[3] - bG) / aG
        say(f"{s}\t{NAME[s]}\t{f:.1f}%\t{g:.1f}%\t{100*pres[s][0]:.1f}% ({int(round(pres[s][0]*size[0])):,} of {int(size[0]):,})\t{e[1]:.3f}\t{e[2]:.3f}")
        names.append(NAME[s]); vals.append(f); cols.append(colr)

fig, ax = plt.subplots(figsize=(7.5, 0.28 * len(names) + 1.2))
y = np.arange(len(names))[::-1]
ax.barh(y, vals, color=cols)
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8)
for x in (0, 50, 100):
    ax.axvline(x, color="0.75", lw=0.8, ls="--")
ax.set_xlabel("calibrated A. flavidus share (%)")
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in GROUPS]
ax.legend(handles, [g for g, _, _ in GROUPS], fontsize=7, loc="lower right")
fig.tight_layout(); fig.savefig("specific/fig_flavidus_share.png", dpi=200)
say("\nOutputs: specific/flavidus_share.txt, specific/fig_flavidus_share.png")
