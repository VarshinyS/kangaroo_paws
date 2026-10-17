import csv, sys
sys.setrecursionlimit(10000)

def parse(s):
    s = s.replace("\n", "").replace("\r", "").strip().rstrip(";")
    adj, label, pos = {}, {}, [0]
    def new():
        n = len(adj); adj[n] = []; return n
    def read_label():
        i = pos[0]
        if i < len(s) and s[i] == "'":
            j = s.index("'", i + 1); pos[0] = j + 1; return s[i+1:j]
        j = i
        while j < len(s) and s[j] not in ":,();": j += 1
        pos[0] = j; return s[i:j].strip()
    def read_len():
        if pos[0] < len(s) and s[pos[0]] == ":":
            j = pos[0] + 1
            while j < len(s) and s[j] not in ",();": j += 1
            v = float(s[pos[0]+1:j]); pos[0] = j; return v
        return 0.0
    def node():
        n = new()
        if s[pos[0]] == "(":
            pos[0] += 1
            while True:
                c = node(); l = read_len()
                adj[n].append((c, l)); adj[c].append((n, l))
                ch = s[pos[0]]; pos[0] += 1
                if ch == ")": break
            read_label()
        else:
            label[n] = read_label()
        return n
    root = node()
    return adj, label, root

def run(infile, outfile, outgroup):
    adj, label, root = parse(open(infile).read())
    tips = set(label.values())
    out = tips & outgroup
    print(f"{infile}: {len(tips)} tips, outgroup tips found: {sorted(out)}")
    sub = {}
    def below(n, p):
        kids = [m for m, l in adj[n] if m != p]
        S = {label[n]} if not kids else set().union(*[below(m, n) for m in kids])
        sub[(p, n)] = S
        return S
    below(root, None)
    def w(n, p):
        kids = [(m, l) for m, l in adj[n] if m != p]
        if not kids: return label[n]
        return "(" + ",".join(f"{w(m, n)}:{l:.8f}" for m, l in kids) + ")"
    for (p, n), S in sub.items():
        if p is None: continue
        L = dict(adj[n])[p]
        if S == out:
            first, second = (n, p), (p, n)
        elif tips - S == out:
            first, second = (p, n), (n, p)
        else:
            continue
        nwk = f"({w(*first)}:{L/2:.8f},{w(*second)}:{L/2:.8f});"
        open(outfile, "w").write(nwk + "\n")
        print(f"  rooted between the two clades -> {outfile}")
        return
    print("  NOT DONE: Conostylis + Blancoa do not form one side of any branch in this tree")

og = set()
with open("meta.csv", newline="", encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        sp = (row.get("species") or "").lower()
        if "conostylis" in sp or "blancoa" in sp:
            og.add(row["metaindex"].strip())

for c in ["cp", "mt", "nuc"]:
    run(f"tree_{c}_unrooted.nwk", f"tree_{c}_pure53_2clade.nwk", og)
