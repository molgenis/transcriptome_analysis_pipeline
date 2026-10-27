#!/usr/bin/env python3
"""
Gene annotation table from the GTF that is also used for counting.

The same table annotates the OUTRIDER and the FRASER results, so both use the
same gene models and nothing is fetched from the internet (no biomaRt) and
nothing depends on installed annotation packages or caches.

Output columns (tab separated, one row per `gene` feature):
  gene_id, gene_id_nover, gene_name, gene_type, chr, start, end, strand
  gene_id_nover = gene_id without the version suffix (ENSG...14 -> ENSG...)
  chr is kept as written in the GTF; start/end are 1-based, inclusive.

Works with GENCODE (gene_type) and Ensembl (gene_biotype) GTFs, plain or .gz.
Standard library only (Python >= 3.6).
"""
import argparse
import csv
import gzip
import re
import sys

ATTR = re.compile(r'(\S+)\s+"([^"]*)"')


def open_text(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path, "rt")


def parse_gtf_genes(path):
    """yield one dict per `gene` line of the GTF"""
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(ATTR.findall(f[8]))
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue
            yield {
                "gene_id": gene_id,
                "gene_id_nover": gene_id.split(".")[0],
                "gene_name": attrs.get("gene_name", ""),
                "gene_type": attrs.get("gene_type", attrs.get("gene_biotype", "")),
                "chr": f[0],
                "start": int(f[3]),
                "end": int(f[4]),
                "strand": f[6] if f[6] in ("+", "-") else "*",
            }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("gtf")
    p.add_argument("--out", default="gene_annotation.tsv")
    a = p.parse_args(argv)

    cols = ["gene_id", "gene_id_nover", "gene_name", "gene_type", "chr",
            "start", "end", "strand"]
    n = 0
    with open(a.out, "wt", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for rec in parse_gtf_genes(a.gtf):
            w.writerow(rec)
            n += 1
    if n == 0:
        print("No gene features found in " + a.gtf, file=sys.stderr)
        return 1
    print("Wrote %d genes to %s" % (n, a.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
