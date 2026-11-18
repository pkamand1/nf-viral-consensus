#!/usr/bin/env python3
"""Combine coverage, depth, consensus, and variant outputs into one table per run.

Writes sample_summary.tsv (one row per sample, with a pass/fail call on genome
coverage) and all_variants.tsv (every annotated variant from every sample).
"""
import argparse
import csv
import gzip
import os


def stem(path, suffix):
    base = os.path.basename(path)
    return base[: -len(suffix)] if base.endswith(suffix) else base


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--coverage", nargs="+", required=True)
    p.add_argument("--depth", nargs="+", required=True)
    p.add_argument("--consensus", nargs="+", required=True)
    p.add_argument("--variants", nargs="+", required=True)
    p.add_argument("--min-depth", type=int, default=10)
    p.add_argument("--min-genome-fraction", type=float, default=0.90)
    args = p.parse_args()

    rows = {}
    for path in args.coverage:
        s = stem(path, ".coverage.tsv")
        with open(path) as handle:
            rec = next(csv.DictReader(handle, delimiter="\t"))
        rows[s] = {"sample": s, "mapped_reads": int(rec["numreads"]), "mean_depth": float(rec["meandepth"])}

    for path in args.depth:
        s = stem(path, ".depth.tsv.gz")
        total = covered = 0
        with gzip.open(path, "rt") as handle:
            for line in handle:
                total += 1
                if int(line.split("\t")[2]) >= args.min_depth:
                    covered += 1
        rows[s]["genome_length"] = total
        rows[s][f"pct_genome_depth_{args.min_depth}x"] = round(100 * covered / total, 2) if total else 0

    for path in args.consensus:
        s = stem(path, ".consensus.fasta")
        with open(path) as handle:
            seq = "".join(line.strip() for line in handle if not line.startswith(">"))
        rows[s]["consensus_length"] = len(seq)
        rows[s]["consensus_n"] = seq.upper().count("N")

    all_variants = []
    for path in args.variants:
        s = stem(path, ".variants.tsv")
        with open(path) as handle:
            recs = list(csv.DictReader(handle, delimiter="\t"))
        all_variants.extend(recs)
        rows[s]["n_snv"] = sum(r["type"] == "SNV" for r in recs)
        rows[s]["n_indel"] = sum(r["type"] in ("INS", "DEL") for r in recs)

    key = f"pct_genome_depth_{args.min_depth}x"
    fields = ["sample", "mapped_reads", "mean_depth", "genome_length", key,
              "consensus_length", "consensus_n", "n_snv", "n_indel", "qc_status"]
    with open("sample_summary.tsv", "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for s in sorted(rows):
            r = rows[s]
            r["mean_depth"] = round(r["mean_depth"], 1)
            r["qc_status"] = "PASS" if r.get(key, 0) >= 100 * args.min_genome_fraction else "FAIL"
            w.writerow({k: r.get(k, "") for k in fields})

    with open("all_variants.tsv", "w", newline="") as out:
        if all_variants:
            w = csv.DictWriter(out, fieldnames=list(all_variants[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(all_variants)
        else:
            out.write("sample\tchrom\tpos\tref\talt\ttype\tqual\tdepth\talt_reads\talt_freq\tgene\treview_flag\n")


if __name__ == "__main__":
    main()
