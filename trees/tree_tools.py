import sys, re, csv

class Node:
    def __init__(self, name=None, length=0.0):
        self.name, self.length, self.children, self.parent = name, length, [], None

def parse(text):
    s = re.sub(r"\s+", "", text).rstrip(";")
    name_re = re.compile(r"[^,():;]*")
    len_re = re.compile(r":([-+0-9.eE]+)")
    def read(pos):
        node = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                child, pos = read(pos)
                child.parent = node
                node.children.append(child)
                if s[pos] == ",":
                    pos += 1
                    continue
                pos += 1  # the closing bracket
                break
        m = name_re.match(s, pos)
        node.name = m.group(0) or None
        pos = m.end()
        m = len_re.match(s, pos)
        if m:
            node.length = float(m.group(1))
            pos = m.end()
        return node, pos
    return read(0)[0]

def leaves(n):
    if not n.children:
        return [n]
    out = []
    for c in n.children:
        out += leaves(c)
    return out

def write(n):
    s = n.name if not n.children else "(" + ",".join(write(c) for c in n.children) + ")" + (n.name or "")
    if n.parent is not None:
        s += f":{n.length:.8f}"
    return s

def graph(root):
    adj, names, count = {}, {}, [0]
    def add(n):
        count[0] += 1
        i = count[0]
        names[i] = None if n.children else n.name
        adj.setdefault(i, [])
        for c in n.children:
            j = add(c)
            adj[i].append((j, c.length))
            adj[j].append((i, c.length))
        return i
    add(root)
    return adj, names

def reroot(root, outgroup):
    adj, names = graph(root)
    L = [i for i, n in names.items() if n == outgroup][0]
    nb, length = adj[L][0]
    def build(i, parent, plen):
        node = Node(names[i], plen)
        for j, l in adj[i]:
            if j != parent:
                c = build(j, i, l)
                c.parent = node
                node.children.append(c)
        return node
    a = Node(outgroup, length / 2)
    b = build(nb, L, length / 2)
    new = Node()
    for c in (a, b):
        c.parent = new
        new.children.append(c)
    return new

def distances_from(root, src):
    adj, names = graph(root)
    start = [i for i, n in names.items() if n == src][0]
    dist, stack = {start: 0.0}, [start]
    while stack:
        i = stack.pop()
        for j, l in adj[i]:
            if j not in dist:
                dist[j] = dist[i] + l
                stack.append(j)
    return {names[i]: d for i, d in dist.items() if names[i]}

def clades(n, out):
    ls = frozenset(l.name for l in leaves(n))
    out.append((ls, n))
    for c in n.children:
        clades(c, out)

def summary(tree_file, meta):
    species = {}
    with open(meta, newline="") as f:
        for r in csv.DictReader(f):
            species[r["metaindex"].strip()] = r["species"].strip()
    root = parse(open(tree_file).read())
    names = [l.name for l in leaves(root)]
    print("leaves:", len(names))
    allc = []
    clades(root, allc)
    groups = {}
    for n in names:
        groups.setdefault(species.get(n, "unknown"), []).append(n)
    print("\nMonophyly by species")
    for sp, ids in sorted(groups.items()):
        if len(ids) < 2:
            print(f"  {sp}: 1 sample")
            continue
        target = frozenset(ids)
        mrca = min((c for c, _ in allc if target <= c), key=len)
        extra = sorted(mrca - target)
        if not extra:
            print(f"  {sp} (n={len(ids)}): monophyletic")
        else:
            shown = ", ".join(f"{x} ({species.get(x)})" for x in extra[:6])
            more = f" and {len(extra) - 6} more" if len(extra) > 6 else ""
            print(f"  {sp} (n={len(ids)}): NOT monophyletic; {len(extra)} other samples in its clade: {shown}{more}")
    print("\nClosest relatives of the discordant samples (tree distance)")
    for s in ("KP1060", "KP1063", "KP1000"):
        if s in names:
            d = distances_from(root, s)
            near = sorted((v, k) for k, v in d.items() if k != s)[:4]
            print(f"  {s} ({species.get(s)}): " + ", ".join(f"{k} ({species.get(k)}) {v:.3f}" for v, k in near))
    print("\nOrder of leaves from the root (read the tree on the PDF for the splits)")
    print("  " + " ".join(names))

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "root":
        t = reroot(parse(open(sys.argv[2]).read()), sys.argv[4])
        open(sys.argv[3], "w").write(write(t) + ";\n")
    elif cmd == "summary":
        summary(sys.argv[2], sys.argv[3])
