#!/usr/bin/env python3
"""Generate the synthetic truth set used by the `test` profile.

Three samples are simulated from the SARS-CoV-2 reference (MN908947.3) with
known variants, so that the pipeline's variant calls and consensus sequences
can be checked against the truth (assets/test/truth_variants.tsv):

  sample_A  12 SNVs and one 6 bp deletion
  sample_B  the variants of sample_A plus 4 further SNVs (a descendant lineage)
  sample_C  an unrelated set of 10 SNVs, with no reads over positions 10001-10500,
            so that the consensus must mask that region with N

Reads are 2 x 150 bp, about 60x coverage, with a 0.2% substitution error rate.
Run from the repository root:  python assets/test/make_test_data.py
"""
import gzip
import random

rng = random.Random(20261117)
REF_PATH = "assets/reference/MN908947.3.fasta"
OUT = "assets/test"
COMP = str.maketrans("ACGT", "TGCA")


def read_fasta(path):
    with open(path) as f:
        return "".join(line.strip() for line in f if not line.startswith(">"))


def pick_snvs(ref, n, exclude):
    snvs = {}
    while len(snvs) < n:
        pos = rng.randint(500, len(ref) - 500)          # 1-based position
        if pos in exclude or any(abs(pos - p) < 30 for p in list(snvs) + list(exclude)):
            continue
        base = ref[pos - 1]
        if base not in "ACGT":
            continue
        snvs[pos] = rng.choice([b for b in "ACGT" if b != base])
    return snvs


def apply_variants(ref, snvs, deletion=None):
    seq = list(ref)
    for pos, alt in snvs.items():
        seq[pos - 1] = alt
    if deletion:
        start, length = deletion                          # delete bases start+1 .. start+length (1-based)
        del seq[start:start + length]
    return "".join(seq)


def simulate(genome, sample, coverage=60, read_len=150, insert_mean=320, skip=None):
    n_pairs = int(len(genome) * coverage / (2 * read_len))
    with gzip.open(f"{OUT}/{sample}_R1.fastq.gz", "wt") as r1, gzip.open(f"{OUT}/{sample}_R2.fastq.gz", "wt") as r2:
        written = 0
        while written < n_pairs:
            insert = max(read_len, int(rng.gauss(insert_mean, 30)))
            start = rng.randint(0, len(genome) - insert)
            if skip and start < skip[1] and start + insert > skip[0]:
                continue
            frag = genome[start:start + insert]
            a, b = frag[:read_len], frag.translate(COMP)[::-1][:read_len]
            a = "".join(rng.choice("ACGT".replace(c, "")) if rng.random() < 0.002 else c for c in a)
            b = "".join(rng.choice("ACGT".replace(c, "")) if rng.random() < 0.002 else c for c in b)
            q = "I" * read_len
            r1.write(f"@{sample}_{written}/1\n{a}\n+\n{q}\n")
            r2.write(f"@{sample}_{written}/2\n{b}\n+\n{q}\n")
            written += 1


def main():
    ref = read_fasta(REF_PATH)
    deletion = (21990, 6)                                  # deletes positions 21991-21996 in the spike gene
    del_pos = set(range(21980, 22010))
    set_a = pick_snvs(ref, 12, del_pos)
    extra_b = pick_snvs(ref, 4, set(set_a) | del_pos)
    set_c = pick_snvs(ref, 10, set(set_a) | set(extra_b) | del_pos | set(range(9950, 10550)))

    samples = {
        "sample_A": (set_a, deletion, None),
        "sample_B": ({**set_a, **extra_b}, deletion, None),
        "sample_C": (set_c, None, (10000, 10500)),
    }
    with open(f"{OUT}/truth_variants.tsv", "w") as truth:
        truth.write("sample\tpos\tref\talt\ttype\n")
        for sample, (snvs, dele, skip) in samples.items():
            genome = apply_variants(ref, snvs, dele)
            simulate(genome, sample, skip=skip)
            for pos in sorted(snvs):
                truth.write(f"{sample}\t{pos}\t{ref[pos - 1]}\t{snvs[pos]}\tSNV\n")
            if dele:
                start, length = dele
                truth.write(f"{sample}\t{start}\t{ref[start - 1:start + length]}\t{ref[start - 1]}\tDEL\n")

    with open(f"{OUT}/samplesheet_test.csv", "w") as sheet:
        sheet.write("sample,fastq_1,fastq_2\n")
        for sample in samples:
            sheet.write(f"{sample},{sample}_R1.fastq.gz,{sample}_R2.fastq.gz\n")


if __name__ == "__main__":
    main()
