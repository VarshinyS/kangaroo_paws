import os
import sys
import subprocess
import glob
import logging
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
CP_READS_DIR    = "/group/shooshared/kpdata/anigo_filter/filtered/cp_reads"
NONCP_READS_DIR = "/group/shooshared/kpdata/anigo_filter/filtered/noncp_reads"
KMER_DIR        = os.path.expanduser("~/varshiny/kmers")
HISTO_DIR       = os.path.expanduser("~/varshiny/histograms")
SUB_DIR         = os.path.expanduser("~/varshiny/subsampled")
RESULTS_DIR     = os.path.expanduser("~/varshiny/results")
LOG_DIR         = os.path.expanduser("~/varshiny/logs")
PLOTS_DIR       = os.path.expanduser("~/varshiny/plots")

KMER_SIZE           = "21"
KMER_MEMORY         = "1G"
KMER_THREADS        = "8"
MIN_DEPTH_CHLORO    = "10"
MIN_DEPTH_NONCP     = "2"

GENOME_SIZES = {
    'chloroplast': 165000,
    'noncp':       700000000,
}
READ_LENGTH = 151

SUBSAMPLE_FRACTIONS = [0.1, 0.25, 0.5, 0.75, 1.0]
HIGH_SAMPLE = "/group/shooshared/kpdata/anigo_filter/filtered/cp_reads/KP1063.cp_1.fastq.gz"

TEST_MODE    = True
TEST_SAMPLES = ['KP0001', 'KP0002']

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
    log = get_logger("step1_setup")
    dirs = [
        KMER_DIR,
        os.path.join(KMER_DIR, "chloroplast"),
        os.path.join(KMER_DIR, "noncp"),
        HISTO_DIR,
        os.path.join(HISTO_DIR, "chloroplast"),
        os.path.join(HISTO_DIR, "noncp"),
        SUB_DIR, RESULTS_DIR, LOG_DIR, PLOTS_DIR,
    ]
    for d in dirs:
        os.makedirs(d, exist_ok=True)
        log.info(f"Dir ready: {d}")
    log.info("All directories created.")

def get_samples(compartment):
    if compartment == 'chloroplast':
        pattern = os.path.join(CP_READS_DIR, "*.cp_1.fastq.gz")
        suffix  = ".cp_1.fastq.gz"
    else:
        pattern = os.path.join(NONCP_READS_DIR, "*.noncp_1.fastq.gz")
        suffix  = ".noncp_1.fastq.gz"
    files   = sorted(glob.glob(pattern))
    samples = [os.path.basename(f).replace(suffix, "") for f in files]
    if TEST_MODE:
        samples = [s for s in samples if s in TEST_SAMPLES]
    return samples

def get_read_pair(sample, compartment):
    if compartment == 'chloroplast':
        r1 = os.path.join(CP_READS_DIR,    f"{sample}.cp_1.fastq.gz")
        r2 = os.path.join(CP_READS_DIR,    f"{sample}.cp_2.fastq.gz")
    else:
        r1 = os.path.join(NONCP_READS_DIR, f"{sample}.noncp_1.fastq.gz")
        r2 = os.path.join(NONCP_READS_DIR, f"{sample}.noncp_2.fastq.gz")
    return r1, r2

def count_reads_in_file(filepath):
    cmd    = f"seqtk comp {filepath} | wc -l"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    try:
        return int(result.stdout.strip())
    except ValueError:
        return 0

def calculate_coverage(sample, compartment):
    r1, r2      = get_read_pair(sample, compartment)
    genome_size = GENOME_SIZES[compartment]
    n_r1        = count_reads_in_file(r1)
    n_r2        = count_reads_in_file(r2)
    total_bases = (n_r1 + n_r2) * READ_LENGTH
    depth       = total_bases / genome_size if genome_size > 0 else 0
    return {
        'sample': sample, 'compartment': compartment,
        'r1_reads': n_r1, 'r2_reads': n_r2,
        'total_bases': total_bases, 'genome_size': genome_size,
        'depth': round(depth, 2),
    }

def calculate_all_coverage():
    log = get_logger("step2_coverage")
    log.info("=== STEP 2: Coverage calculation ===")
    rows = []
    for comp in ['chloroplast', 'noncp']:
        for sample in get_samples(comp):
            log.info(f"  {sample} [{comp}]")
            row = calculate_coverage(sample, comp)
            rows.append(row)
            log.info(f"    depth={row['depth']}x")
    df  = pd.DataFrame(rows)
    out = os.path.join(RESULTS_DIR, "coverage_summary.csv")
    df.to_csv(out, index=False)
    log.info(f"Coverage summary saved: {out}")
    return out


def run_jellyfish(sample, compartment):
    log      = get_logger("step3_jellyfish")
    r1, r2   = get_read_pair(sample, compartment)
    jf_out   = os.path.join(KMER_DIR,  compartment, f"{sample}.jf")
    histo_out = os.path.join(HISTO_DIR, compartment, f"{sample}.histo")

    if os.path.exists(jf_out):
        log.info(f"  [SKIP] {sample} [{compartment}] — .jf exists")
        return jf_out, histo_out

    cmd_count = (
        f"jellyfish count "
        f"-m {KMER_SIZE} -s {KMER_MEMORY} -t {KMER_THREADS} -C "
        f"-o {jf_out} "
        f"<(zcat {r1}) <(zcat {r2})"
    )
    log.info(f"  jellyfish count: {sample} [{compartment}]")
    result = subprocess.run(["bash", "-c", cmd_count], capture_output=True, text=True)
    if result.returncode != 0:
        log.error(f"  FAILED:\n{result.stderr}")
        raise RuntimeError(f"jellyfish count failed for {sample} [{compartment}]")

    cmd_histo = f"jellyfish histo {jf_out} > {histo_out}"
    log.info(f"  jellyfish histo: {sample} [{compartment}]")
    result2 = subprocess.run(["bash", "-c", cmd_histo], capture_output=True, text=True)
    if result2.returncode != 0:
        log.error(f"  FAILED:\n{result2.stderr}")
        raise RuntimeError(f"jellyfish histo failed for {sample} [{compartment}]")

    log.info(f"  done: {sample} [{compartment}]")
    return jf_out, histo_out

def run_all_jellyfish():
    log = get_logger("step3_jellyfish")
    log.info("=== STEP 3: Jellyfish counting ===")
    outputs = {}
    for comp in ['chloroplast', 'noncp']:
        outputs[comp] = {}
        for sample in get_samples(comp):
            jf, histo = run_jellyfish(sample, comp)
            outputs[comp][sample] = {'jf': jf, 'histo': histo}
    log.info("Jellyfish done.")
    return outputs

def parse_histo(histo_path):
    return pd.read_csv(histo_path, sep=' ', header=None,
                       names=['frequency', 'count'])

def find_kmer_repeats(sample, compartment, histo_path):
    df = parse_histo(histo_path)
    return {
        'sample':          sample,
        'compartment':     compartment,
        'total_distinct':  df['count'].sum(),
        'total_kmer_mass': (df['frequency'] * df['count']).sum(),
        'unique_kmers':    df[df['frequency'] == 1]['count'].sum(),
        'low_rep_2_10':    df[(df['frequency'] >= 2)  & (df['frequency'] <= 10)]['count'].sum(),
        'mid_rep_11_100':  df[(df['frequency'] > 10)  & (df['frequency'] <= 100)]['count'].sum(),
        'high_rep_100+':   df[df['frequency'] > 100]['count'].sum(),
    }

def run_repeat_analysis(jellyfish_outputs):
    log = get_logger("step4_repeats")
    log.info("=== STEP 4: K-mer repeat analysis ===")
    rows = []
    for comp, samples in jellyfish_outputs.items():
        for sample, paths in samples.items():
            if os.path.exists(paths['histo']):
                row = find_kmer_repeats(sample, comp, paths['histo'])
                rows.append(row)
                log.info(f"  {sample} [{comp}]: {row['total_distinct']} distinct k-mers")
    df  = pd.DataFrame(rows)
    out = os.path.join(RESULTS_DIR, "repeat_summary.csv")
    df.to_csv(out, index=False)
    log.info(f"Repeat summary saved: {out}")
    return out

def plot_kmer_distribution(jellyfish_outputs):
    log = get_logger("step5_plots")
    log.info("=== STEP 5: K-mer distribution plots ===")
    for comp, samples in jellyfish_outputs.items():
        for sample, paths in samples.items():
            if not os.path.exists(paths['histo']):
                continue
            df  = parse_histo(paths['histo'])
            df  = df[df['frequency'] <= 500]
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.plot(df['frequency'], df['count'], linewidth=1)
            ax.set_xlabel('K-mer frequency')
            ax.set_ylabel('Number of distinct k-mers')
            ax.set_title(f'{sample} — {comp} k-mer distribution (k={KMER_SIZE})')
            ax.set_yscale('log')
            plt.tight_layout()
            out = os.path.join(PLOTS_DIR, f"{sample}_{comp}_kmer_dist.png")
            plt.savefig(out, dpi=150)
            plt.close()
            log.info(f"  Saved: {out}")

def plot_compartment_comparison(jellyfish_outputs):
    log = get_logger("step5_plots")
    log.info("=== STEP 5: Compartment comparison plots ===")
    cp_samples    = set(jellyfish_outputs.get('chloroplast', {}).keys())
    noncp_samples = set(jellyfish_outputs.get('noncp', {}).keys())
    for sample in cp_samples & noncp_samples:
        fig, ax = plt.subplots(figsize=(12, 5))
        for comp, colour in [('chloroplast', 'green'), ('noncp', 'purple')]:
            histo = jellyfish_outputs[comp][sample]['histo']
            if os.path.exists(histo):
                df = parse_histo(histo)
                df = df[df['frequency'] <= 500]
                ax.plot(df['frequency'], df['count'], colour, linewidth=1, label=comp)
        ax.set_xlabel('K-mer frequency')
        ax.set_ylabel('Number of distinct k-mers')
        ax.set_title(f'{sample} — chloroplast vs nonCP (k={KMER_SIZE})')
        ax.set_yscale('log')
        ax.legend()
        plt.tight_layout()
        out = os.path.join(PLOTS_DIR, f"{sample}_compartment_comparison.png")
        plt.savefig(out, dpi=150)
        plt.close()
        log.info(f"  Saved: {out}")


def run_subsampling():
    log = get_logger("step6_subsampling")
    log.info("=== STEP 6: Subsampling experiment ===")
    histo_paths = {}
    for frac in SUBSAMPLE_FRACTIONS:
        tag    = str(frac).replace('.', 'p')
        sub_r1 = os.path.join(SUB_DIR, f"KP1063_sub_{tag}_R1.fastq.gz")
        jf_out = os.path.join(SUB_DIR, f"KP1063_sub_{tag}.jf")
        histo  = os.path.join(SUB_DIR, f"KP1063_sub_{tag}.histo")

        if not os.path.exists(sub_r1):
            cmd = f"seqtk sample -s42 {HIGH_SAMPLE} {frac} | gzip > {sub_r1}"
            log.info(f"  seqtk sample {frac}")
            subprocess.run(["bash", "-c", cmd], check=True)

        if not os.path.exists(jf_out):
            cmd_jf = (
                f"jellyfish count "
                f"-m {KMER_SIZE} -s {KMER_MEMORY} -t {KMER_THREADS} -C "
                f"-o {jf_out} <(zcat {sub_r1})"
            )
            log.info(f"  jellyfish count {frac}")
            subprocess.run(["bash", "-c", cmd_jf], check=True)

        if not os.path.exists(histo):
            subprocess.run(["bash", "-c", f"jellyfish histo {jf_out} > {histo}"], check=True)

        histo_paths[frac] = histo
        log.info(f"  Done fraction {frac}")

    fig, ax = plt.subplots(figsize=(12, 6))
    for frac, histo in histo_paths.items():
        if os.path.exists(histo):
            df = parse_histo(histo)
            df = df[df['frequency'] <= 200]
            ax.plot(df['frequency'], df['count'], label=f"{int(frac*100)}% reads")
    ax.set_xlabel('K-mer frequency')
    ax.set_ylabel('Number of distinct k-mers')
    ax.set_title(f'KP1063 chloroplast — subsampling experiment (k={KMER_SIZE})')
    ax.set_yscale('log')
    ax.legend()
    plt.tight_layout()
    out = os.path.join(PLOTS_DIR, "subsampling_comparison.png")
    plt.savefig(out, dpi=150)
    plt.close()
    log.info(f"Subsampling plot saved: {out}")
    return histo_paths

def build_kmer_matrix(compartment, jellyfish_outputs):
    log = get_logger("step7_matrix")
    log.info(f"=== STEP 7: Build k-mer matrix [{compartment}] ===")
    sample_kmer_dicts = {}
    for sample, paths in jellyfish_outputs.get(compartment, {}).items():
        jf = paths['jf']
        if not os.path.exists(jf):
            log.warning(f"  Missing .jf: {jf}")
            continue
        dump_path = os.path.join(KMER_DIR, compartment, f"{sample}_dump.txt")
        if not os.path.exists(dump_path):
            cmd = f"jellyfish dump -c {jf} > {dump_path}"
            log.info(f"  dumping {sample}")
            subprocess.run(["bash", "-c", cmd], check=True)
        kd = {}
        with open(dump_path) as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) == 2:
                    kd[parts[0]] = int(parts[1])
        sample_kmer_dicts[sample] = kd
        log.info(f"  {sample}: {len(kd)} k-mers")

    if not sample_kmer_dicts:
        log.error("No k-mer data found.")
        return None

    df  = pd.DataFrame(sample_kmer_dicts).T.fillna(0).astype(int)
    out = os.path.join(RESULTS_DIR, f"kmer_matrix_{compartment}.csv")
    df.to_csv(out)
    log.info(f"  Matrix shape: {df.shape} — saved: {out}")
    return out


def normalise_matrix(matrix_csv_path, compartment):
    log = get_logger("step8_normalise")
    log.info(f"=== STEP 8: Normalise matrix [{compartment}] ===")
    if matrix_csv_path is None or not os.path.exists(matrix_csv_path):
        log.error(f"Matrix not found: {matrix_csv_path}")
        return None
    df     = pd.read_csv(matrix_csv_path, index_col=0)
    normed = normalize(df.values, norm='l1')
    df_n   = pd.DataFrame(normed, index=df.index, columns=df.columns)
    out    = os.path.join(RESULTS_DIR, f"kmer_matrix_{compartment}_normalised.csv")
    df_n.to_csv(out)
    log.info(f"  Saved: {out}")
    return out


def main():
    print("=" * 60)
    print("Aim 1: K-mer Data Preparation — Kangaroo Paw Project")
    print(f"TEST_MODE = {TEST_MODE}")
    if TEST_MODE:
        print(f"Samples:   {TEST_SAMPLES}")
    print("=" * 60)

    setup_dirs()

    coverage_csv = calculate_all_coverage()
    print(f"\n[Step 2] Coverage summary: {coverage_csv}")

    jellyfish_outputs = run_all_jellyfish()
    print(f"\n[Step 3] Jellyfish done.")

    repeat_csv = run_repeat_analysis(jellyfish_outputs)
    print(f"\n[Step 4] Repeat summary: {repeat_csv}")

    plot_kmer_distribution(jellyfish_outputs)
    plot_compartment_comparison(jellyfish_outputs)
    print(f"\n[Step 5] Plots saved to: {PLOTS_DIR}")

    run_subsampling()
    print(f"\n[Step 6] Subsampling done.")

    normalised_outputs = {}
    for comp in ['chloroplast', 'noncp']:
        matrix_csv = build_kmer_matrix(comp, jellyfish_outputs)
        normed_csv = normalise_matrix(matrix_csv, comp)
        normalised_outputs[comp] = normed_csv

    print("\n" + "=" * 60)
    print("Aim 1 complete!")
    print(f"Results: {RESULTS_DIR}")
    for comp, path in normalised_outputs.items():
        print(f"  [{comp}] normalised matrix: {path}")
    print("=" * 60)

if __name__ == "__main__":
    main()
