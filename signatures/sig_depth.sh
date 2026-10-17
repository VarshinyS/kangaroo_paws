#!/bin/bash
#SBATCH --job-name=sig_depth
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/sms049/shared/vselvanayagam/sig20m/sig_depth.log

module use -p /group/shooshared/modules >/dev/null 2>&1
module load kmertools_s >/dev/null 2>&1

WORK=${WORK:-/scratch/sms049/shared/vselvanayagam/sig20m}
KM=${KM:-$WORK/kmers}; DC=${DC:-$WORK/depth}; ST=${ST:-$WORK/stats}; OUT=${OUT:-$WORK/sigdepth}
PY=${PY:-/scratch/sms049/shared/vselvanayagam/ml_venv/bin/python3}
F_REF=${F_REF:-"KP0042 KP0043 KP0052 KP1008"}
R_REF=${R_REF:-"KP1067 KP0058 KP0054"}
P_REF=${P_REF:-"KP1066 KP1222 KP1223"}
TESTS=${TESTS:-"KP0003 KP0006 KP0029 KP0041"}
TOPN=${TOPN:-1000000}
TH=${THREADS:-8}
mkdir -p "$OUT/w"; cd "$OUT" || exit 1
fail() { echo "FAILED: $*"; exit 1; }
size() { kmc_tools info "$1" | awk -F: '/total k-mers/{gsub(/ /,"",$2); print $2}'; }
run() { kmc_tools -t$TH "$@" > w/last_kmc.txt 2>&1 || { tail -3 w/last_kmc.txt; fail "kmc_tools $*"; }; }

binarise() { [ -s "w/bin_$1.kmc_pre" ] || run transform "$KM/$1" set_counts 1 "w/bin_$1"; }

fold() {   # mode(sum|union) out db...
  local mode=$1 out=$2; shift 2
  local acc=$1; shift
  local p="-ocmax"; [ "$mode" = sum ] && p="-ocsum"
  if [ $# -eq 0 ]; then run simple "$acc" "$acc" union "$out" -ocmax; return; fi
  local i=0 n=$# tgt nx
  for nx in "$@"; do
    i=$((i+1)); tgt="w/fold_tmp_$i"; [ $i -eq $n ] && tgt=$out
    run simple "$acc" "$nx" union "$tgt" $p
    acc=$tgt
  done
  rm -f w/fold_tmp_*.kmc_*
}

for id in $F_REF $R_REF $P_REF; do [ -s "$KM/$id.kmc_pre" ] || fail "missing reference database $id"; binarise $id; done
bins() { for i in "$@"; do printf 'w/bin_%s ' "$i"; done; }

echo -e "species\treferences\tin_all_references\tsignature" > sig_sizes.tsv
for sp in F R P; do
  case $sp in F) mine=$F_REF; others="$R_REF $P_REF";; R) mine=$R_REF; others="$F_REF $P_REF";; P) mine=$P_REF; others="$F_REF $R_REF";; esac
  n=$(echo $mine | wc -w)
  fold sum "w/sum_$sp" $(bins $mine)
  run transform "w/sum_$sp" histogram "consistency_$sp.txt" -ci1 -cx100
  run transform "w/sum_$sp" -ci$n reduce "w/all_$sp"
  fold union "w/oth_$sp" $(bins $others)
  run simple "w/all_$sp" "w/oth_$sp" kmers_subtract "w/sig_$sp"
  echo -e "$sp\t$n\t$(size w/all_$sp)\t$(size w/sig_$sp)" >> sig_sizes.tsv
  echo "signature $sp done $(date +%H:%M:%S)"
done
cat sig_sizes.tsv

# presence: fraction of each signature found in each test sample and each depth level of the held-out sample
echo -e "name\tkind\tpairs\tspecies\tshared\tsignature\tpct" > recovery.tsv
measure() {   # name kind pairs db
  for sp in F R P; do
    run simple "$4" "w/sig_$sp" intersect w/x -ocmax
    s=$(size w/x); g=$(size w/sig_$sp)
    echo -e "$1\t$2\t$3\t$sp\t$s\t$g\t$(awk -v a=$s -v b=$g 'BEGIN{printf "%.2f", (b>0?100*a/b:0)}')" >> recovery.tsv
  done
}
for id in $TESTS; do
  [ -s "$KM/$id.kmc_pre" ] || { echo "no database for test $id"; continue; }
  pairs=$(awk -F'\t' '{print $5}' "$ST/$id.tsv" 2>/dev/null); measure "$id" test "${pairs:-NA}" "$KM/$id"
done
while read -r name pairs; do measure "$name" depth "$pairs" "$DC/$name"; done < "$DC/levels.tsv"
echo "presence done $(date +%H:%M:%S)"

# frequency: are the same high-frequency k-mers recovered at lower depth, after normalising by total k-mers?
full="$DC/KP0055_all"
run transform "$full" histogram w/hist_full.txt -ci1 -cx65535
TTOP=$(awk -v N=$TOPN '{c[$1]=$2; if($1>m)m=$1} END{s=0; for(k=m;k>=1;k--){s+=c[k]; if(s>=N){print k; exit}} print 1}' w/hist_full.txt)
echo "high-frequency threshold: count >= $TTOP"
run transform "$full" -ci$TTOP reduce w/top
run transform w/top dump w/top_full.txt
occ() { awk '{s+=$1*$2} END{printf "%.0f", s}' "$1"; }
echo -e "name\tpairs\ttotal_occurrences\tdump" > freq_inputs.tsv
while read -r name pairs; do
  run transform "$DC/$name" histogram "w/hist_$name.txt" -ci1 -cx65535
  if [ "$name" = "KP0055_all" ]; then d=w/top_full.txt
  else run simple w/top "$DC/$name" intersect "w/top_$name" -ocright; run transform "w/top_$name" dump "w/top_$name.txt"; d="w/top_$name.txt"; fi
  echo -e "$name\t$pairs\t$(occ w/hist_$name.txt)\t$d" >> freq_inputs.tsv
done < "$DC/levels.tsv"
echo "frequency done $(date +%H:%M:%S)"

cat > w/summarise.py << 'PYEOF'
import numpy as np, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
def rankdata(a):
    a = np.asarray(a, dtype=float); s = np.argsort(a, kind="mergesort"); inv = np.empty_like(s); inv[s] = np.arange(len(a))
    srt = a[s]; obs = np.r_[True, srt[1:] != srt[:-1]]; dense = obs.cumsum()[inv]
    cnt = np.r_[np.nonzero(obs)[0], len(obs)]; return 0.5 * (cnt[dense] + cnt[dense - 1] + 1)
def load_dump(p):
    d = {}
    for line in open(p):
        k, c = line.split(); d[k] = float(c)
    return d
rows = list(csv.DictReader(open("freq_inputs.tsv"), delimiter="\t"))
full = [r for r in rows if r["name"] == "KP0055_all"][0]
top = load_dump(full["dump"]); tot_full = float(full["total_occurrences"])
out = []
for r in rows:
    lev = load_dump(r["dump"]); tot = float(r["total_occurrences"])
    shared = [k for k in lev if k in top]
    a = np.array([top[k] / tot_full * 1e6 for k in shared]); b = np.array([lev[k] / tot * 1e6 for k in shared])
    rho = float(np.corrcoef(rankdata(a), rankdata(b))[0, 1]) if len(shared) > 2 else float("nan")
    out.append((r["name"], int(float(r["pairs"])), round(100 * len(shared) / len(top), 2), round(rho, 3),
                round(float(np.median(b / a)), 3) if len(shared) else float("nan")))
out.sort(key=lambda x: x[1])
with open("frequency.tsv", "w") as f:
    f.write("name\tpairs\tpct_high_frequency_kmers_recovered\tspearman_normalised_counts\tmedian_ratio_normalised\n")
    for o in out: f.write("\t".join(str(x) for x in o) + "\n")
print(open("frequency.tsv").read())
rec = list(csv.DictReader(open("recovery.tsv"), delimiter="\t"))
col = {"F": "#2E8B57", "R": "#FF4500", "P": "#20B2AA"}; lab = {"F": "flavidus signature", "R": "rufus signature", "P": "pulcherrimus signature"}
names = {"KP0003": "'green'", "KP0006": "'Frosty Yellow'", "KP0029": "Yellow Gem", "KP0041": "cross"}
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for sp in "FRP":
    d = sorted([(float(r["pairs"]), float(r["pct"])) for r in rec if r["kind"] == "depth" and r["species"] == sp])
    if d: ax[0].plot([x / 1e6 for x, _ in d], [y for _, y in d], "-o", color=col[sp], label=lab[sp] + " (held-out 'red' at each depth)")
    for r in rec:
        if r["kind"] == "test" and r["species"] == sp and r["pairs"] not in ("NA", ""):
            ax[0].scatter(float(r["pairs"]) / 1e6, float(r["pct"]), marker="s", color=col[sp], edgecolor="black", zorder=3)
            ax[0].annotate(names.get(r["name"], r["name"]), (float(r["pairs"]) / 1e6, float(r["pct"])), fontsize=7, xytext=(3, 3), textcoords="offset points")
ax[0].set_xlabel("repeat-filtered read pairs (millions)"); ax[0].set_ylabel("% of signature k-mers present")
ax[0].set_title("Presence: signature recovery against depth"); ax[0].legend(fontsize=7)
x = [o[1] / 1e6 for o in out]
ax[1].plot(x, [o[2] for o in out], "-o", color="#333333", label="% of high-frequency k-mers recovered")
ax2 = ax[1].twinx(); ax2.plot(x, [o[3] for o in out], "-s", color="#1f77b4", label="Spearman correlation of normalised counts")
ax[1].set_xlabel("repeat-filtered read pairs (millions)"); ax[1].set_ylabel("% recovered"); ax2.set_ylabel("Spearman correlation")
ax[1].set_title("Frequency: high-frequency k-mers of the held-out 'red'")
h1, l1 = ax[1].get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax[1].legend(h1 + h2, l1 + l2, fontsize=7, loc="lower right")
plt.tight_layout(); plt.savefig("depth_curve.png", dpi=200)
print("saved depth_curve.png")
PYEOF
"$PY" w/summarise.py || fail "summary"
cat recovery.tsv
echo "Done!"
