import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load(f):
    d = np.loadtxt("hist/" + f)
    return d[(d[:, 0] >= 5) & (d[:, 1] > 0)]

fig, ax = plt.subplots(2, 2, figsize=(11, 8))
a = ax.ravel()

d = load("KP0052_cp.txt"); a[0].plot(d[:, 0], d[:, 1], color="green")
a[0].set_title("A  Chloroplast (KP0052)", loc="left")
d = load("KP0052_mt.txt"); a[1].plot(d[:, 0], d[:, 1], color="darkorange")
a[1].set_title("B  Mitochondria (KP0052)", loc="left")

for f, lab, c in [("KP0042_nuc.txt", "Full nuclear", "black"),
                  ("KP0042_rp.txt", "Repeat", "firebrick"),
                  ("KP0042_sc.txt", "Repeat-depleted", "steelblue")]:
    d = load(f); a[2].plot(d[:, 0], d[:, 1], color=c, label=lab)
a[2].set_title("C  Nuclear (KP0042)", loc="left"); a[2].legend(frameon=False)

lv = ["2M", "4M", "6M", "8M", "10M", "12M", "all"]
cols = plt.cm.viridis(np.linspace(0, 0.9, len(lv)))
for L, c in zip(lv, cols):
    d = load(f"KP0055_{L}.txt"); a[3].plot(d[:, 0], d[:, 1], color=c, label=L)
a[3].set_title("D  Nuclear depth series (KP0055)", loc="left")
a[3].legend(frameon=False, title="read pairs", fontsize=8)

for x in a:
    x.set_xscale("log"); x.set_yscale("log")
    x.set_xlabel("k-mer count"); x.set_ylabel("number of distinct k-mers")
plt.tight_layout(); plt.savefig("fig_spectra_panel.png", dpi=300)
print("saved fig_spectra_panel.png")
