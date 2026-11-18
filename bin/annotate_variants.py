#!/usr/bin/env python3
"""Tabulate the variants in a VCF with read support, the gene that contains each one, and review flags.

review_flag marks calls that deserve a second look: primer_site or near_primer (amplicon data),
indel, frameshift_length, and below_consensus_af (minority variants not applied to the consensus).
"""
import argparse
import gzip


def load_genes(path):
    genes = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            genes.append((f[0], int(f[3]), int(f[4]), attrs.get("Name", attrs.get("gene", "."))))
    return genes


def load_primers(path):
    primers = []
    if path:
        with open(path) as handle:
            for line in handle:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3 and not line.startswith(("#", "track")):
                    primers.append((f[0], int(f[1]) + 1, int(f[2])))   # BED to 1-based inclusive
    return primers


def review_flags(chrom, pos, ref, alt, freq, primers, min_af, window=10):
    flags = []
    end = pos + len(ref) - 1
    if any(c == chrom and start <= end and pos <= stop for c, start, stop in primers):
        flags.append("primer_site")
    elif any(c == chrom and start - window <= end and pos <= stop + window for c, start, stop in primers):
        flags.append("near_primer")
    if len(ref) != len(alt):
        flags.append("indel")
        if (len(alt) - len(ref)) % 3 != 0:
            flags.append("frameshift_length")
    if freq != "" and float(freq) < min_af:
        flags.append("below_consensus_af")
    return ";".join(flags) if flags else "."


def gene_at(genes, chrom, pos):
    hits = [name for c, start, end, name in genes if c == chrom and start <= pos <= end]
    return ",".join(hits) if hits else "intergenic"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True)
    p.add_argument("--gff", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--primer-bed", default=None)
    p.add_argument("--consensus-min-af", type=float, default=0.75)
    args = p.parse_args()

    genes = load_genes(args.gff)
    primers = load_primers(args.primer_bed)
    print("\t".join(["sample", "chrom", "pos", "ref", "alt", "type", "qual", "depth", "alt_reads",
                     "alt_freq", "gene", "review_flag"]))
    opener = gzip.open if args.vcf.endswith(".gz") else open
    with opener(args.vcf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alts, qual = f[0], int(f[1]), f[3], f[4], f[5]
            fmt = dict(zip(f[8].split(":"), f[9].split(":"))) if len(f) > 9 else {}
            ad = [int(x) for x in fmt.get("AD", "").split(",") if x not in ("", ".")]
            depth = sum(ad) if ad else ""
            for i, alt in enumerate(alts.split(",")):
                alt_reads = ad[i + 1] if len(ad) > i + 1 else ""
                freq = f"{alt_reads / depth:.3f}" if ad and depth else ""
                if len(ref) == len(alt) == 1:
                    vtype = "SNV"
                elif len(ref) > len(alt):
                    vtype = "DEL"
                elif len(ref) < len(alt):
                    vtype = "INS"
                else:
                    vtype = "MNV"
                flag = review_flags(chrom, pos, ref, alt, freq, primers, args.consensus_min_af)
                print("\t".join(map(str, [args.sample, chrom, pos, ref, alt, vtype, qual, depth,
                                          alt_reads, freq, gene_at(genes, chrom, pos), flag])))


if __name__ == "__main__":
    main()
