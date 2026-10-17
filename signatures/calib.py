import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "specific/dose"
SETS = ["flavidus", "rufus", "pulcherrimus", "rufpul"]
SPECIES = {"flavidus": ["KP0042", "KP0043", "KP0052", "KP1008"],
           "rufus": ["KP1067", "KP0058", "KP0054"],
           "pulcherrimus": ["KP1066", "KP1222", "KP1223"]}
REFS = [s for v in SPECIES.values() for s in v]
TESTS = ["KP0055", "KP0001", "KP0003", "KP0006", "KP0029", "KP0041"]
HYB = ["KP0030", "KP0068", "KP0038", "KP0065", "KP0040"]
MIX = ["MIX_FR10", "MIX_FR25", "MIX_FR50", "MIX_FR75", "MIX_FP10", "MIX_FP25", "MIX_FP50", "MIX_FP75"]
NAME = {"KP0055": "'red' (flavidus)", "KP0001": "'Stripie' 1 (flavidus)", "KP0003": "'green' (flavidus)",
        "KP0006": "'Frosty Yellow' (rufus-type)", "KP0029": "'Yellow Gem'", "KP0041": "'Stripie' x 'Yellow Gem'",
        "KP0030": "'Landscape Orange'", "KP0068": "'Landscape Scarlet'", "KP0038": "'Landscape Gold'",
        "KP0065": "'Orange Cross'", "KP0040": "'Early Spring' x 'Stripie'"}

pairs = {}
for fn in ["robust/thresholds.tsv", "robust/thresholds_hyb.tsv", "robust/thresholds_mix.tsv"]:
    if os.path.exists(fn):
        for line in open(fn):
            f = line.split()
            if f[0] != "sample":
                pairs[f[0]] = float(f[1])
member = {}
for i, g in enumerate(SETS):
    for line in open(f"{D}/set_{g}.txt"):
        member[line.split()[0]] = i
size = np.bincount(list(member.values()), minlength=4)

dose = {}
for s in REFS + TESTS + HYB + MIX:
    if not (os.path.exists(f"{D}/{s}.txt") and s in pairs):
        continue
    tot = np.zeros(4)
    for line in open(f"{D}/{s}.txt"):
        k, c = line.split()
        i = member.get(k)
        if i is not None:
            tot[i] += int(c)
    dose[s] = tot / pairs[s] * 1e6
ref = np.array([np.mean([dose[s][0] for s in SPECIES["flavidus"]]),
                np.mean([dose[s][1] for s in SPECIES["rufus"]]),
                np.mean([dose[s][2] for s in SPECIES["pulcherrimus"]]),
                np.mean([dose[s][3] for s in SPECIES["rufus"] + SPECIES["pulcherrimus"]])])
est = {s: dose[s] / ref for s in dose}

out = open("specific/calib_summary.txt", "w")
def say(m):
    print(m, flush=True); out.write(m + "\n")
say("Set sizes (abundant specific k-mers): " + ", ".join(f"{g} {n:,}" for g, n in zip(SETS, size)))
say("Estimated share = summed count of the set's k-mers per million pairs / mean of that species' references")
say("sample\tname\tflavidus\trufus\tpulcherrimus\trufus+pulcherrimus")
for s in REFS + TESTS + HYB + MIX:
    if s in est:
        e = est[s]
        say(f"{s}\t{NAME.get(s, s)}\t{e[0]:.3f}\t{e[1]:.3f}\t{e[2]:.3f}\t{e[3]:.3f}")

mixes = [m for m in MIX if m in est]
cal = None
if len(mixes) >= 4:
    tF = np.array([1 - int(m[-2:]) / 100 for m in mixes]); eF = np.array([est[m][0] for m in mixes])
    tG = 1 - tF; eG = np.array([est[m][3] for m in mixes])
    aF, bF = np.polyfit(tF, eF, 1); aG, bG = np.polyfit(tG, eG, 1)
    r2 = lambda t, e, a, b: 1 - np.sum((e - (a * t + b)) ** 2) / np.sum((e - e.mean()) ** 2)
    say(f"\nCalibration from {len(mixes)} read mixes:")
    say(f"  flavidus set: estimated = {aF:.3f} x true + {bF:.3f} (R2 {r2(tF, eF, aF, bF):.3f})")
    say(f"  rufus+pulcherrimus set: estimated = {aG:.3f} x true + {bG:.3f} (R2 {r2(tG, eG, aG, bG):.3f})")
    for m in mixes:
        say(f"  {m}: true flavidus {1-int(m[-2:])/100:.2f}, estimated {est[m][0]:.3f}; "
            f"rufus {est[m][1]:.3f}, pulcherrimus {est[m][2]:.3f}")
    cal = (aF, bF, aG, bG)
else:
    say("\nFewer than 4 read mixes available: no calibration, uncalibrated estimates only")

say("\nCalibrated shares for test plants and nursery lines (flavidus / rufus+pulcherrimus group; which of the two):")
names, cF, cG = [], [], []
for s in TESTS + HYB:
    if s not in est:
        continue
    e = est[s]
    if cal:
        f = (e[0] - cal[1]) / cal[0]; g = (e[3] - cal[3]) / cal[2]
    else:
        f, g = e[0], e[3]
    side = "pulcherrimus" if e[2] > e[1] else "rufus"
    if e[1] < 0.05 and e[2] < 0.05:
        side = "-"
    say(f"{s}\t{NAME[s]:30s}\tflavidus {100*f:5.1f}%\trufus+pulcherrimus {100*g:5.1f}%\t"
        f"(rufus {e[1]:.3f}, pulcherrimus {e[2]:.3f} -> {side})")
    names.append(NAME[s]); cF.append(100 * f); cG.append(100 * g)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
ax = axes[0]
if len(mixes) >= 4:
    for m in mixes:
        mk = "o" if m.startswith("MIX_FR") else "s"
        ax.plot(100 * (1 - int(m[-2:]) / 100), 100 * est[m][0], mk, color="#2e8b57")
        ax.plot(100 * int(m[-2:]) / 100, 100 * est[m][3], mk, color="#ff4500")
    ax.plot([0, 100], [0, 100], "k--", lw=1)
    ax.plot([], [], "o", color="0.4", label="flavidus + rufus mixes")
    ax.plot([], [], "s", color="0.4", label="flavidus + pulcherrimus mixes")
    ax.plot([], [], "-", color="#2e8b57", label="flavidus set")
    ax.plot([], [], "-", color="#ff4500", label="rufus + pulcherrimus set")
    ax.legend(fontsize=8)
ax.set_xlabel("true share in read mix (%)"); ax.set_ylabel("estimated share (%)")
ax.set_title("Calibration with read mixes")
ax = axes[1]
y = np.arange(len(names))[::-1]
ax.barh(y + 0.2, cF, 0.4, color="#2e8b57", label="flavidus")
ax.barh(y - 0.2, cG, 0.4, color="#ff4500", label="rufus + pulcherrimus")
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8)
ax.axvline(50, color="0.7", ls="--", lw=1); ax.set_xlim(0, 120)
ax.set_xlabel("calibrated share (%)"); ax.set_title("Test plants and nursery lines"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig("specific/fig_calib.png", dpi=200)
say("Outputs: specific/calib_summary.txt, specific/fig_calib.png")
