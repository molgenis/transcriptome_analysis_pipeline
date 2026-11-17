#!/usr/bin/env python3
"""
Merge the leafcutterMD.R outputs of all chunks and make per-sample result tables.

Input (from every chunk): <prefix>_pVals.txt, <prefix>_effSize.txt,
<prefix>_clusterPvals.txt as written by leafcutterMD.R (rows = introns or
clusters, columns = samples; empty files are fine).

Per sample the cluster p-values are adjusted with Benjamini-Hochberg over all
tested clusters of that sample. Clusters with padjust <= --padj are written
with the intron that is most affected (smallest intron p-value, then largest
absolute effect size), its effect size, and the overlapping gene from the
gene annotation table (protein coding genes preferred, then the smallest).

Output: <sample>_result_table_leafcutterMD.tsv for every sample of the samplesheet,
combined_samples_result_table_leafcutterMD.tsv, and the merged raw tables
leafcutterMD_{pVals,effSize,clusterPvals}.txt.gz
Standard library only (Python >= 3.6).
"""
import argparse
import bisect
import csv
import glob
import gzip
import os
import sys


# ---------------------------------------------------------------------------
# reading
# ---------------------------------------------------------------------------
def read_table(files):
    """-> (samples, {row_name: [values as strings]}), merged over chunk files"""
    samples = None
    rows = {}
    for f in files:
        if os.path.getsize(f) == 0:
            continue
        with open(f, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            if samples is None:
                samples = header
            elif header != samples:
                raise ValueError("Samples differ between chunks: " + f)
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) == len(samples) + 1:
                    rows[parts[0]] = parts[1:]
    return samples or [], rows


def bh_adjust(pvals):
    """Benjamini-Hochberg adjusted p-values (same order as the input)"""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [0.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * n / rank)
        adj[i] = prev
    return adj


# ---------------------------------------------------------------------------
# genes
# ---------------------------------------------------------------------------
def norm_chr(c):
    c = c[3:] if c.startswith("chr") else c
    return "MT" if c == "M" else c


class Genes:
    def __init__(self, path):
        self.by_chr = {}
        if not path:
            return
        with open(path, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                self.by_chr.setdefault(norm_chr(r["chr"]), []).append(
                    (int(r["start"]), int(r["end"]), r["gene_id"],
                     r["gene_name"], r["gene_type"]))
        self.starts = {}
        self.maxend = {}
        for c, lst in self.by_chr.items():
            lst.sort()
            self.starts[c] = [g[0] for g in lst]
            m, run = [], 0
            for g in lst:
                run = max(run, g[1])
                m.append(run)
            self.maxend[c] = m

    def find(self, chrom, start, end):
        c = norm_chr(chrom)
        if c not in self.by_chr:
            return None
        lst = self.by_chr[c]
        i = bisect.bisect_right(self.starts[c], end) - 1
        hits = []
        while i >= 0 and self.maxend[c][i] >= start:
            g = lst[i]
            if g[0] <= end and g[1] >= start:
                hits.append(g)
            i -= 1
        if not hits:
            return None
        hits.sort(key=lambda g: (g[4] != "protein_coding", g[1] - g[0]))
        return hits[0]


# ---------------------------------------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", default=".", help="directory with the chunk outputs")
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--annotation", default=None, help="gene_annotation.tsv")
    p.add_argument("--padj", type=float, default=0.1,
                   help="only clusters with padjust <= this are written")
    a = p.parse_args(argv)

    out_files = {}
    for kind in ("pVals", "effSize", "clusterPvals"):
        out_files[kind] = sorted(glob.glob(os.path.join(a.dir, "*_%s.txt" % kind)))
        if not out_files[kind]:
            print("No *_%s.txt files found" % kind, file=sys.stderr)
            return 1

    samples, clu_p = read_table(out_files["clusterPvals"])
    s2, int_p = read_table(out_files["pVals"])
    s3, int_e = read_table(out_files["effSize"])

    # merged raw tables
    for kind, (smp, rows) in (("pVals", (s2, int_p)), ("effSize", (s3, int_e)),
                              ("clusterPvals", (samples, clu_p))):
        with gzip.open("leafcutterMD_%s.txt.gz" % kind, "wt") as fh:
            fh.write("\t".join(smp) + "\n")
            for name, vals in rows.items():
                fh.write(name + "\t" + "\t".join(vals) + "\n")

    with open(a.samplesheet, "rt") as fh:
        sheet_ids = [r["sampleID"] for r in csv.DictReader(fh, delimiter="\t")]

    # introns per cluster
    introns = {}
    for name in int_p:
        introns.setdefault(name.split(":")[3], []).append(name)

    genes = Genes(a.annotation)
    cols = ["sampleID", "cluster", "chr", "start", "end", "strand", "nIntrons",
            "effectSize", "intronPvalue", "clusterPvalue", "padjust",
            "hgncSymbol", "EnsemblID"]
    combined = []
    per_sample = {s: [] for s in sheet_ids}
    clusters = list(clu_p)
    n_tested = len(clusters)

    for si, sample in enumerate(samples):
        if sample not in per_sample:
            continue
        pv = []
        for c in clusters:
            try:
                pv.append(float(clu_p[c][si]))
            except ValueError:
                pv.append(1.0)
        adj = bh_adjust(pv)
        for c, p_c, q in zip(clusters, pv, adj):
            if q > a.padj:
                continue
            best = None
            for iname in introns.get(c, []):
                try:
                    ip = float(int_p[iname][si])
                    ef = float(int_e[iname][si])
                except (ValueError, KeyError):
                    continue
                key = (ip, -abs(ef))
                if best is None or key < best[0]:
                    best = (key, iname, ip, ef)
            if best is None:
                continue
            _, iname, ip, ef = best
            chrom, start, end, clu = iname.split(":")[:4]
            gene = genes.find(chrom, int(start), int(end))
            row = {
                "sampleID": sample, "cluster": c, "chr": chrom,
                "start": start, "end": end, "strand": clu.rsplit("_", 1)[-1],
                "nIntrons": len(introns.get(c, [])),
                "effectSize": "%.4f" % ef, "intronPvalue": "%.3g" % ip,
                "clusterPvalue": "%.3g" % p_c, "padjust": "%.3g" % q,
                "hgncSymbol": gene[3] if gene and gene[3] else "NA",
                "EnsemblID": gene[2] if gene else "NA",
            }
            per_sample[sample].append(row)
            combined.append(row)

    def write(path, rows):
        with open(path, "wt", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    for sample, rows in per_sample.items():
        write("%s_result_table_leafcutterMD.tsv" % sample, rows)
    write("combined_samples_result_table_leafcutterMD.tsv", combined)
    print("%d clusters tested in %d samples; %d outlier clusters with padjust <= %g"
          % (n_tested, len(samples), len(combined), a.padj))
    missing = [s for s in sheet_ids if s not in samples]
    if missing:
        print("Samples without MD results (not in the junction data): "
              + ", ".join(missing))
    return 0


if __name__ == "__main__":
    sys.exit(main())
