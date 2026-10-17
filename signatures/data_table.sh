#!/bin/bash
#SBATCH --job-name=data_table
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/sms049/shared/vselvanayagam/sig20m/data_table.log

module use -p /group/shooshared/modules >/dev/null 2>&1
module load kmertools_s >/dev/null 2>&1

D=${D:-/group/shooshared/kpdata/anigo_filter}
OUT=${OUT:-/scratch/sms049/shared/vselvanayagam/sig20m}
SIGST=${SIGST:-$OUT/stats}
META=${META:-/group/shooshared/kpdata/anigo_filter/status/kpgenomics_sequence_inventory_2026-08-13.csv}
PY=${PY:-/scratch/sms049/shared/vselvanayagam/ml_venv/bin/python3}
cd "$OUT" || exit 1
mkdir -p datatab

# reads used in each filtering step, from the BBMap logs (last "Reads Used" line in each log)
used() { awk '/Reads Used:/{r=$3; b=$4} END{gsub(/\(/,"",b); print r"\t"b}' "$1" 2>/dev/null; }
echo -e "sample\tcp_step_reads\tcp_step_bases\tmt_step_reads\tmt_step_bases" > datatab/logs.tsv
for f in $D/filtered/stats/*.bbmap.log; do
  id=$(basename "$f" .bbmap.log)
  echo -e "$id\t$(used "$f")\t$(used "$D/mt_filter/stats/$id.bbmap.log")" >> datatab/logs.tsv
done
echo "logs read: $(($(wc -l < datatab/logs.tsv)-1)) samples"

# mitochondrial pairs = reads in each mt_1 file (mates are kept together); fall back to mt_2 if mt_1 is unreadable
seqkit --quiet stats -T -j 8 $D/mt_filter/mt_reads/*.mt_1.fastq.gz 2>/dev/null | awk -F'\t' 'NR>1{n=split($1,a,"/"); id=a[n]; sub(/\.mt_1\.fastq\.gz$/,"",id); print id"\t"$4"\tmt_1"}' > datatab/mt_counts.tsv
for f in $D/mt_filter/mt_reads/*.mt_1.fastq.gz; do
  id=$(basename "$f" .mt_1.fastq.gz)
  if ! grep -q "^$id	" datatab/mt_counts.tsv; then
    n=$(seqkit --quiet stats -T "$D/mt_filter/mt_reads/$id.mt_2.fastq.gz" 2>/dev/null | awk -F'\t' 'NR==2{print $4}')
    echo -e "$id\t${n:-NA}\tmt_2" >> datatab/mt_counts.tsv
  fi
done
echo "mitochondrial files counted: $(wc -l < datatab/mt_counts.tsv)"

cat > datatab/build.py << 'PYEOF'
import csv, math, os, sys, glob
out, sigst, meta = sys.argv[1], sys.argv[2], sys.argv[3]
species = {}
try:
    for r in csv.DictReader(open(meta, newline="")):
        species[r["metaindex"].strip()] = r["species"].strip()
except Exception as e:
    print("metadata not read:", e)
logs = {r["sample"]: r for r in csv.DictReader(open(f"{out}/datatab/logs.tsv"), delimiter="\t")}
mt = {}
for line in open(f"{out}/datatab/mt_counts.tsv"):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 3 and p[1] not in ("", "NA"):
        mt[p[0]] = (int(p[1]), p[2])
rep = {}
for f in glob.glob(f"{sigst}/*.tsv"):
    p = open(f).read().split("\t")
    try:
        sub, kept = int(p[2]), int(p[4])
        rep[p[0].split()[-1]] = (p[3].strip(), round(100 * (sub - kept) / sub, 2), sub)
    except Exception:
        pass
def num(x):
    try: return int(x)
    except Exception: return None
groups = [20, 160, 200, 300]
rows = []
for s, r in sorted(logs.items()):
    raw_r, raw_b, nc_r = num(r["cp_step_reads"]), num(r["cp_step_bases"]), num(r["mt_step_reads"])
    if not raw_r: continue
    raw = raw_r // 2
    rl = round(raw_b / raw_r, 1) if raw_b else None
    noncp = nc_r // 2 if nc_r else None
    cp = raw - noncp if noncp is not None else None
    m = mt.get(s, (None, ""))[0]
    nuc = noncp - m if (noncp is not None and m is not None) else None
    g = min(groups, key=lambda x: abs(math.log(raw / 1e6 / x))) if raw else None
    pct = lambda x: round(100 * x / raw, 2) if x is not None else None
    rp = rep.get(s, ("", "", ""))
    rows.append([s, species.get(s, ""), g, raw, rl, cp, pct(cp), m, pct(m), nuc, pct(nuc),
                 mt.get(s, (None, ""))[1], rp[0], rp[1], rp[2]])
hdr = ["sample", "species", "depth_group_million_pairs", "raw_pairs", "read_length", "chloroplast_pairs", "chloroplast_pct",
       "mitochondrial_pairs", "mitochondrial_pct", "nuclear_pairs", "nuclear_pct", "mt_count_from",
       "repeat_reads_mapped_pct_in_subsample", "repeat_pairs_removed_pct_in_subsample", "repeat_subsample_pairs"]
with open(f"{out}/datatab/data_table.tsv", "w") as f:
    f.write("\t".join(hdr) + "\n")
    for r in rows: f.write("\t".join("" if v is None else str(v) for v in r) + "\n")
print(f"{len(rows)} samples written to datatab/data_table.tsv")
print("\nsummary by depth group (median, min-max)")
def med(v): v = sorted(v); n = len(v); return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2
for g in groups:
    sel = [r for r in rows if r[2] == g]
    if not sel: continue
    line = f"{g}M group: n={len(sel)}"
    for name, i in (("raw pairs (M)", 3), ("chloroplast %", 6), ("mitochondrial %", 8), ("nuclear %", 10)):
        v = [r[i] / 1e6 if i == 3 else r[i] for r in sel if r[i] is not None]
        if v: line += f" | {name} {med(v):.2f} ({min(v):.2f}-{max(v):.2f})"
    print(line)
PYEOF
"$PY" datatab/build.py "$OUT" "$SIGST" "$META"
echo "Done!"
