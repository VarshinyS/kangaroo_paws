import sys, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tree_tools import parse, leaves

COL = {"Anigozanthos_flavidus": "#2E8B57", "Anigozanthos_humilis": "#FFD700", "Anigozanthos_rufus": "#FF4500",
       "Anigozanthos_kalbarriensis": "#00CED1", "Anigozanthos_bicolor": "#9370DB", "Anigozanthos_viridis": "#32CD32",
       "Anigozanthos_manglesii": "#DAA520", "Anigozanthos_onycis": "#FF69B4", "Anigozanthos_pulcherimus": "#20B2AA",
       "Anigozanthos_preisii": "#F4A460", "Anigozanthos_gabrielae": "#87CEEB", "Conostylis_setosa": "#DDA0DD",
       "Conostylis_candicans": "#EE82EE", "Blancoa_canescens": "#CD5C5C", "Macropidia_fuliginosa": "#8B0000",
       "Anigozanthos_hybrid": "#808080"}
HIGHLIGHT = {"KP1060", "KP1063", "KP1000"}

def order_by(node, pos):
    # put the children of every node in the order of their mean position in the other tree
    for c in node.children:
        order_by(c, pos)
    node.children.sort(key=lambda c: sum(pos[l.name] for l in leaves(c)) / len(leaves(c)))

def number(node, ypos, counter):
    # y position for every node; leaves are numbered in order from the top
    if not node.children:
        ypos[id(node)] = counter[0]
        counter[0] += 1
    else:
        for c in node.children:
            number(c, ypos, counter)
        ypos[id(node)] = sum(ypos[id(c)] for c in node.children) / len(node.children)

def depths(node, d, out):
    out[id(node)] = d
    for c in node.children:
        depths(c, d + c.length, out)

def draw(ax, root, ypos, xdepth, x_of, color="#444444"):
    def rec(n):
        for c in n.children:
            ax.plot([x_of(xdepth[id(n)]), x_of(xdepth[id(c)])], [ypos[id(c)]] * 2, color=color, lw=0.7)
            rec(c)
        if n.children:
            ys = [ypos[id(c)] for c in n.children]
            ax.plot([x_of(xdepth[id(n)])] * 2, [min(ys), max(ys)], color=color, lw=0.7)
    rec(root)

def crossings(a, b):
    n = len(a)
    return sum(1 for i in range(n) for j in range(i + 1, n) if (a[i] - a[j]) * (b[i] - b[j]) < 0)

def main(t1_file, t2_file, meta, labels_file, title1, title2, out):
    species = {}
    with open(meta, newline="") as f:
        for r in csv.DictReader(f):
            species[r["metaindex"].strip()] = r["species"].strip()
    names = {}
    if labels_file != "none":
        seen = False
        for line in open(labels_file):
            if line.strip() == "DATA":
                seen = True
                continue
            if seen and "\t" in line:
                k, v = line.rstrip("\n").split("\t", 1)
                names[k] = v

    A = parse(open(t1_file).read())
    B = parse(open(t2_file).read())
    posA0 = {l.name: i for i, l in enumerate(leaves(A))}
    order_by(B, posA0)
    yA, yB = {}, {}
    number(A, yA, [0])
    number(B, yB, [0])
    dA, dB = {}, {}
    depths(A, 0.0, dA)
    depths(B, 0.0, dB)
    maxA, maxB = max(dA.values()), max(dB.values())

    xa = lambda d: d / maxA * 1.0
    xb = lambda d: 3.4 - d / maxB * 1.0

    n = len(leaves(A))
    fig, ax = plt.subplots(figsize=(15, max(8, n * 0.26)))
    draw(ax, A, yA, dA, xa)
    draw(ax, B, yB, dB, xb)

    leafA = {l.name: yA[id(l)] for l in leaves(A)}
    leafB = {l.name: yB[id(l)] for l in leaves(B)}
    common = [k for k in leafA if k in leafB]
    ordA = [leafA[k] for k in common]
    ordB = [leafB[k] for k in common]
    for k in common:
        c = COL.get(species.get(k, ""), "#000000")
        hi = k in HIGHLIGHT
        ax.plot([1.45, 1.95], [leafA[k], leafB[k]], color="black" if hi else c, lw=2.2 if hi else 1.0,
                alpha=1.0 if hi else 0.8, zorder=3 if hi else 2)
        lab = names.get(k, k)
        ax.text(1.02, leafA[k], lab, va="center", ha="left", fontsize=6.2, color=c, fontweight="bold" if hi else "normal")
        ax.text(2.38, leafB[k], lab, va="center", ha="right", fontsize=6.2, color=c, fontweight="bold" if hi else "normal")
    ax.set_ylim(n, -1.5)
    ax.set_xlim(-0.05, 3.45)
    ax.text(0.5, -1.1, title1, ha="center", fontsize=11, fontweight="bold")
    ax.text(2.9, -1.1, title2, ha="center", fontsize=11, fontweight="bold")
    ax.set_title(f"Line crossings between the two trees: {crossings(ordA, ordB)} "
                 f"(black lines mark M13_010, MQ4_001 and B7_15)", fontsize=9)
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(out, dpi=220)
    print("saved", out, "| samples:", len(common), "| crossings:", crossings(ordA, ordB))

if __name__ == "__main__":
    main(*sys.argv[1:8])
