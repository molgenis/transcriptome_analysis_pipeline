#!/usr/bin/env python3
"""
STAR SJ.out.tab -> LeafCutter .junc file (BED12, same layout as regtools).

SJ.out.tab columns: 1 chr, 2 first intron base (1-based), 3 last intron base
(1-based), 4 strand (0 undefined, 1 +, 2 -), 5 motif, 6 annotated,
7 unique reads, 8 multi-mapping reads, 9 max overhang.

LeafCutter takes the intron from the BED12 blocks as
  A = chromStart + blockSize1,  B = chromEnd - blockSize2 + 1.
Writing chromStart = first-2, chromEnd = last+1 and block sizes 1,1 gives
A = first-1 and B = last+1, exactly what regtools produces for the same
intron, so STAR, FRASER and regtools junctions can be mixed.

Junctions with undefined strand (0) are written with "?"; LeafCutter's
clustering skips those, as it does for regtools output.
Standard library only (Python >= 3.6).
"""
import argparse
import gzip
import sys

STRAND = {"1": "+", "2": "-", "0": "?"}


def open_text(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path, "rt")


def convert(sj_path, out_path, with_multi=False, min_count=1):
    n_in = n_out = 0
    with open_text(sj_path) as fin, open(out_path, "wt") as fout:
        for line in fin:
            f = line.split()
            if len(f) < 9:
                continue
            n_in += 1
            first, last = int(f[1]), int(f[2])
            count = int(f[6]) + (int(f[7]) if with_multi else 0)
            if count < min_count or first < 2:
                continue
            start, end = first - 2, last + 1
            strand = STRAND.get(f[3], "?")
            fout.write("\t".join([
                f[0], str(start), str(end), "JUNC", str(count), strand,
                str(start), str(end), "255,0,0", "2", "1,1",
                "0,%d" % (end - start - 1)]) + "\n")
            n_out += 1
    return n_in, n_out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("sj", help="SJ.out.tab (plain or .gz)")
    p.add_argument("--out", required=True, help="output .junc file")
    p.add_argument("--with-multimapping", action="store_true",
                   help="also count multi-mapping reads (default: unique only)")
    p.add_argument("--min-count", type=int, default=1)
    a = p.parse_args(argv)
    n_in, n_out = convert(a.sj, a.out, a.with_multimapping, a.min_count)
    if n_in == 0:
        print("No junctions in " + a.sj, file=sys.stderr)
        return 1
    print("%s: %d of %d junctions written" % (a.sj, n_out, n_in))
    return 0


if __name__ == "__main__":
    sys.exit(main())
