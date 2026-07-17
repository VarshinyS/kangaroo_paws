import os
import sys
import subprocess
import logging
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import combinations

CP_READS_DIR    = "/group/shooshared/kpdata/anigo_filter/filtered/cp_reads"
NONCP_READS_DIR = "/group/shooshared/kpdata/anigo_filter/filtered/noncp_reads"
PILOT_DIR    = os.path.expanduser("~/varshiny/pilot")
LOG_DIR      = os.path.expanduser("~/varshiny/logs")

PILOT_SAMPLES = ['KP0055', 'KP0052', 'KP0069', 'KP0051']
PILOT_LABELS  = {
    'KP0055': 'A. flavidus red',
    'KP0052': "A. flavidus 'tall orange'",
    'KP0069': 'A. Bush Pearl',
    'KP0051': 'Pink beauty',
}

KMER_SIZES_CP    = [21, 31, 51, 60, 70, 80, 90, 100]
KMER_SIZES_NONCP = [21, 31, 51]

KMER_MEMORY   = "1G"
KMER_THREADS  = "2"

def get_logger(name):
    os.makedirs(LOG_DIR, exist_ok=True)
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    if not logger.handlers:
        fh = logging.FileHandler(os.path.join(LOG_DIR, f"{name}.log"))
        fh.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
        sh = logging.StreamHandler(sys.stdout)
        sh.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
        logger.addHandler(fh)
        logger.addHandler(sh)
    return logger

def setup_dirs():
    log = get_logger("pilot_setup")
    dirs = [
        PILOT_DIR,
        os.path.join(PILOT_DIR, "kmers", "chloroplast"),
        os.path.join(PILOT_DIR, "kmers", "noncp"),
        os.path.join(PILOT_DIR, "histograms", "chloroplast"),
        os.path.join(PILOT_DIR, "histograms", "noncp"),
        os.path.join(PILOT_DIR, "plots"),
        os.path.join(PILOT_DIR, "results"),
    ]
    for d in dirs:
        os.makedirs(d, exist_ok=True)
        log.info(f"Dir ready: {d}")

def get_read_pair(sample, compartment):
    if compartment == 'chloroplast':
        r1 = os.path.join(CP_READS_DIR, f"{sample}.cp_1.fastq.gz")
        r2 = os.path.join(CP_READS_DIR, f"{sample}.cp_2.fastq.gz")
    else:
        r1 = os.path.join(NONCP_READS_DIR, f"{sample}.noncp_1.fastq.gz")
        r2 = os.path.join(NONCP_READS_DIR, f"{sample}.noncp_2.fastq.gz")
    return r1, r2

def run_jellyfish_pilot(sample, k, compartment):
    log    = get_logger("pilot_jellyfish")
    r1, r2 = get_read_pair(sample, compartment)
    jf     = os.path.join(PILOT_DIR, "kmers", compartment, f"{sample}_k{k}.jf")
    histo  = os.path.join(PILOT_DIR, "histograms", compartment, f"{sample}_k{k}.histo")

    if os.path.exists(histo):
        log.info(f" {sample} k={k} [{compartment}] - histo exists")
        return histo

    cmd_count = (
        f"jellyfish count "
        f"-m {k} -s {KMER_MEMORY} -t {KMER_THREADS} -C "
        f"-o {jf} "
        f"<(zcat {r1}) <(zcat {r2})"
    )
    log.info(f"  jellyfish count: {sample} k={k} [{compartment}]")
    result = subprocess.run(["bash", "-c", cmd_count], capture_output=True, text=True)
    if result.returncode != 0:
        log.error(f"  FAILED: {result.stderr}")
        raise RuntimeError(f"jellyfish failed: {sample} k={k} [{compartment}]")

    cmd_histo = f"jellyfish histo {jf} > {histo}"
    subprocess.run(["bash", "-c", cmd_histo], check=True)

    # Delete .jf immediately after histogram extracted to save disk space
    if os.path.exists(jf):
        os.remove(jf)
        log.info(f"  cleaned .jf: {sample} k={k} [{compartment}]")

    log.info(f"  done: {sample} k={k} [{compartment}]")
    return histo


def run_all_jellyfish():
    log = get_logger("pilot_jellyfish")
    log.info("Pilot: Jellyfish counting")
    outputs = {'chloroplast': {}, 'noncp': {}}
    for k in KMER_SIZES_CP:
        outputs['chloroplast'][k] = {}
        for sample in PILOT_SAMPLES:
            histo = run_jellyfish_pilot(sample, k, 'chloroplast')
            outputs['chloroplast'][k][sample] = {'histo': histo}
    for k in KMER_SIZES_NONCP:
        outputs['noncp'][k] = {}
        for sample in PILOT_SAMPLES:
            histo = run_jellyfish_pilot(sample, k, 'noncp')
            outputs['noncp'][k][sample] = {'histo': histo}
    return outputs

def parse_histo(path):
    return pd.read_csv(path, sep=' ', header=None, names=['frequency', 'count'])


def compute_stats(sample, k, compartment, histo_path):
    df = parse_histo(histo_path)
    total_distinct  = df['count'].sum()
    total_kmer_mass = (df['frequency'] * df['count']).sum()
    unique_kmers    = df[df['frequency'] == 1]['count'].sum()
    error_rate      = unique_kmers / total_distinct if total_distinct > 0 else 0
    return {
        'sample':         sample,
        'variety':        PILOT_LABELS.get(sample, sample),
        'compartment':    compartment,
        'k':              k,
        'total_distinct': total_distinct,
        'total_kmer_mass': total_kmer_mass,
        'unique_kmers':   unique_kmers,
        'error_rate':     round(error_rate, 4),
    }


def compute_all_stats(outputs):
    log = get_logger("pilot_stats")
    log.info("Pilot: Computing stats")
    rows = []
    for compartment, k_dict in outputs.items():
        for k, samples in k_dict.items():
            for sample, paths in samples.items():
                if os.path.exists(paths['histo']):
                    row = compute_stats(sample, k, compartment, paths['histo'])
                    rows.append(row)
                    log.info(f"  {sample} k={k} [{compartment}]: {row['total_distinct']} distinct k-mers")
    df  = pd.DataFrame(rows)
    out = os.path.join(PILOT_DIR, "results", "pilot_kmer_stats.csv")
    df.to_csv(out, index=False)
    log.info(f"Stats saved: {out}")
    return df, out


def compute_jaccard(set_a, set_b):
    intersection = len(set_a & set_b)
    union        = len(set_a | set_b)
    return intersection / union if union > 0 else 0


def load_kmer_set_from_histo(histo_path):
    """Load k-mer counts from histogram for Jaccard — using frequency as proxy."""
    df = parse_histo(histo_path)
    # Use frequency counts as a set proxy — samples with same frequency profile
    # Return set of (frequency, count) tuples for comparison
    return set(zip(df['frequency'].tolist(), df['count'].tolist()))


def compute_all_jaccard(outputs):
    log = get_logger("pilot_jaccard")
    log.info("Pilot: Jaccard similarity (from histograms)")
    rows = []
    for compartment, k_dict in outputs.items():
        for k, samples in k_dict.items():
            log.info(f"  Computing Jaccard k={k} [{compartment}]")
            histo_sets = {}
            for sample, paths in samples.items():
                if os.path.exists(paths['histo']):
                    histo_sets[sample] = load_kmer_set_from_histo(paths['histo'])

            for s1, s2 in combinations(PILOT_SAMPLES, 2):
                if s1 in histo_sets and s2 in histo_sets:
                    j = compute_jaccard(histo_sets[s1], histo_sets[s2])
                    rows.append({
                        'compartment': compartment, 'k': k,
                        'sample1': s1, 'sample2': s2, 'jaccard': round(j, 4)
                    })
                    log.info(f"    {s1} vs {s2}: Jaccard={j:.4f}")

    df  = pd.DataFrame(rows)
    out = os.path.join(PILOT_DIR, "results", "pilot_jaccard.csv")
    df.to_csv(out, index=False)
    log.info(f"Jaccard saved: {out}")
    return df, out


def plot_kmer_distributions(outputs):
    log = get_logger("pilot_plots")
    colours = ['blue', 'orange', 'green', 'red']
    for compartment, k_dict in outputs.items():
        for k, samples in k_dict.items():
            fig, ax = plt.subplots(figsize=(12, 6))
            for (sample, paths), colour in zip(samples.items(), colours):
                if os.path.exists(paths['histo']):
                    df = parse_histo(paths['histo'])
                    df = df[df['frequency'] <= 200]
                    ax.plot(df['frequency'], df['count'],
                            colour, linewidth=1.5, label=sample)
            ax.set_xlabel('K-mer frequency')
            ax.set_ylabel('Number of distinct k-mers')
            ax.set_title(f'{compartment.capitalize()} k-mer distribution - k={k}')
            ax.set_yscale('log')
            ax.legend()
            plt.tight_layout()
            out = os.path.join(PILOT_DIR, "plots", f"kmer_dist_{compartment}_k{k}.png")
            plt.savefig(out, dpi=150)
            plt.close()
            log.info(f"  Saved: {out}")


def plot_kmer_size_comparison(stats_df):
    log = get_logger("pilot_plots")
    for compartment in stats_df['compartment'].unique():
        sub_df = stats_df[stats_df['compartment'] == compartment]
        k_sizes = sorted(sub_df['k'].unique())
        fig, ax = plt.subplots(figsize=(12, 6))
        x      = np.arange(len(PILOT_SAMPLES))
        width  = 0.8 / len(k_sizes)
        colours = plt.cm.viridis(np.linspace(0, 1, len(k_sizes)))
        for i, k in enumerate(k_sizes):
            subset = sub_df[sub_df['k'] == k].set_index('sample')
            vals   = [subset.loc[s, 'total_distinct'] if s in subset.index else 0
                      for s in PILOT_SAMPLES]
            ax.bar(x + i * width, vals, width, label=f'k={k}', color=colours[i])
        ax.set_xlabel('Sample')
        ax.set_ylabel('Total distinct k-mers')
        ax.set_title(f'Total distinct k-mers by k size - {compartment}')
        ax.set_xticks(x + width * (len(k_sizes) - 1) / 2)
        ax.set_xticklabels(PILOT_SAMPLES)
        ax.legend()
        plt.tight_layout()
        out = os.path.join(PILOT_DIR, "plots", f"kmer_size_comparison_{compartment}.png")
        plt.savefig(out, dpi=150)
        plt.close()
        log.info(f"  Saved: {out}")


def plot_jaccard_comparison(jaccard_df):
    log = get_logger("pilot_plots")
    for compartment in jaccard_df['compartment'].unique():
        sub_df = jaccard_df[jaccard_df['compartment'] == compartment]
        mean_j = sub_df.groupby('k')['jaccard'].mean().reset_index()
        fig, ax = plt.subplots(figsize=(8, 5))
        bars = ax.bar([str(k) for k in mean_j['k']], mean_j['jaccard'],
                      color=plt.cm.viridis(np.linspace(0, 1, len(mean_j))))
        ax.set_xlabel('K-mer size (k)')
        ax.set_ylabel('Mean Jaccard similarity')
        ax.set_title(f'Mean Jaccard similarity by k size - {compartment}\n(lower = more distinct = better)')
        for bar, val in zip(bars, mean_j['jaccard']):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.002,
                    f'{val:.4f}', ha='center', va='bottom', fontsize=10)
        plt.tight_layout()
        out = os.path.join(PILOT_DIR, "plots", f"jaccard_comparison_{compartment}.png")
        plt.savefig(out, dpi=150)
        plt.close()
        log.info(f"  Saved: {out}")


def main():
    print("=" * 60)
    print("Pilot K-mer Size Experiment")
    print(f"Samples: {PILOT_SAMPLES}")
    print(f"Chloroplast k sizes: {KMER_SIZES_CP}")
    print(f"NonCP k sizes: {KMER_SIZES_NONCP}")
    print("=" * 60)

    setup_dirs()

    outputs = run_all_jellyfish()
    print("\n Jellyfish done.")

    stats_df, stats_out = compute_all_stats(outputs)
    print(f"\n Stats saved: {stats_out}")

    jaccard_df, jaccard_out = compute_all_jaccard(outputs)
    print(f"\n Jaccard saved: {jaccard_out}")

    plot_kmer_distributions(outputs)
    plot_kmer_size_comparison(stats_df)
    plot_jaccard_comparison(jaccard_df)
    print(f" Plots saved to: {os.path.join(PILOT_DIR, 'plots')}")

    print("\n" + "=" * 60)
    print("PILOT RESULTS SUMMARY")
    print("=" * 60)
    for compartment in ['chloroplast', 'noncp']:
        sub_j = jaccard_df[jaccard_df['compartment'] == compartment]
        if len(sub_j) > 0:
            mean_by_k = sub_j.groupby('k')['jaccard'].mean()
            best_k = mean_by_k.idxmin()
            print(f"\n{compartment.upper()}")
            print(mean_by_k.to_string())
            print(f">>> RECOMMENDED K: k={best_k} <<<")
    print("=" * 60)

if __name__ == "__main__":
    main()
