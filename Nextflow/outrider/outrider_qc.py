#!/usr/bin/env python3
"""
Quality control of the OUTRIDER inputs: finds broken samples and suspicious
strandedness settings BEFORE the counts are merged and fitted.

Two sources are used:
  * the featureCounts summaries (<sampleID>_outrider_count_summary.tsv), only
    present for samples that were counted from a bam. Assigned vs. unassigned
    reads show wrong strandedness and failed/odd libraries.
  * the gene count tables (<sampleID>_outrider_counts.tsv or the expressionFile
    of a sample). These show empty samples, tiny libraries and samples whose
    gene IDs do not match the rest (those become all-zero columns after the
    merge, which breaks the OUTRIDER fit with "L-BFGS-B needs finite values").

Output: a table with metrics, flags and a status (PASS / WARN / FAIL) per
sample. With --stop-on-fail the exit status is 1 when any sample FAILs.

Standard library only (Python >= 3.6).
"""
import argparse
import csv
import gzip
import os
import re
import statistics
import sys

SUMMARY_SUFFIX = "_outrider_count_summary.tsv"

# flags that make a sample unusable for the fit
FAIL_FLAGS = {"MISSING_COUNTS", "ALL_ZERO", "LOW_LIBSIZE", "GENEID_MISMATCH"}


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------
def open_text(path):
    """Open plain or gzipped text, detected from the file content."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_summary(path):
    """featureCounts summary -> dict of status -> reads (first sample column)."""
    stats = {}
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t", quotechar='"')
        next(reader)  # header: "Status", "<bam>"
        for row in reader:
            if len(row) >= 2 and row[0]:
                stats[row[0]] = float(row[1])
    return stats


def summary_metrics(stats):
    """Fractions of the mapped reads (unmapped reads are not informative)."""
    total = sum(stats.values())
    mapped = total - stats.get("Unassigned_Unmapped", 0.0)
    assigned = stats.get("Assigned", 0.0)
    if mapped <= 0:
        return {"assigned": assigned, "mapped": 0.0, "assigned_frac": 0.0,
                "nofeatures_frac": 0.0, "ambiguity_frac": 0.0}
    return {
        "assigned": assigned,
        "mapped": mapped,
        "assigned_frac": assigned / mapped,
        "nofeatures_frac": stats.get("Unassigned_NoFeatures", 0.0) / mapped,
        "ambiguity_frac": stats.get("Unassigned_Ambiguity", 0.0) / mapped,
    }


def read_counts(path):
    """Gene count table (GeneID + one or more sample columns).
    Returns {sample: {"genes": set, "lib": float, "nonzero": int}}."""
    out = {}
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t", quotechar='"')
        header = next(reader)
        samples = header[1:]
        genes = []
        lib = [0.0] * len(samples)
        nonzero = [0] * len(samples)
        for row in reader:
            if not row:
                continue
            genes.append(row[0])
            for i in range(len(samples)):
                try:
                    v = float(row[i + 1])
                except (ValueError, IndexError):
                    v = 0.0
                lib[i] += v
                if v > 0:
                    nonzero[i] += 1
    gene_set = set(genes)
    for i, s in enumerate(samples):
        out[s] = {"genes": gene_set, "n_genes": len(genes), "lib": lib[i],
                  "nonzero": nonzero[i], "file": os.path.basename(path)}
    return out


def read_samplesheet(path):
    with open(path, "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows


def is_set(v):
    return v is not None and str(v).strip() not in ("", "NA")


# --------------------------------------------------------------------------
# statistics helpers
# --------------------------------------------------------------------------
def median(xs):
    return statistics.median(xs) if xs else float("nan")


def mad(xs):
    """Median absolute deviation, scaled to be comparable to a SD."""
    if not xs:
        return float("nan")
    m = statistics.median(xs)
    return 1.4826 * statistics.median([abs(x - m) for x in xs])


def outlier_low(value, ref, k, min_shift):
    """True if value is clearly below the reference group."""
    if len(ref) < 3:
        return False
    m, s = median(ref), mad(ref)
    return value < m - max(k * s, min_shift)


def outlier_high(value, ref, k, min_shift):
    if len(ref) < 3:
        return False
    m, s = median(ref), mad(ref)
    return value > m + max(k * s, min_shift)


# --------------------------------------------------------------------------
# the actual check
# --------------------------------------------------------------------------
def qc_samples(samplesheet_rows, count_tables, summaries, min_lib=1e6,
               min_assigned=5e6, k=4.0, min_group=5, id_overlap_min=0.9):
    """
    samplesheet_rows : list of dicts (sampleID, strandSpecific, ...)
    count_tables     : {sample: {"genes", "n_genes", "lib", "nonzero", ...}}
    summaries        : {sample: {status: reads}}
    returns (records, notes); records is one dict per sample.
    """
    sheet = {r["sampleID"]: r for r in samplesheet_rows}
    sample_ids = [r["sampleID"] for r in samplesheet_rows]

    # reference gene set: IDs present in at least half of the count tables
    tables = list(count_tables.values())
    gene_votes = {}
    for t in {id(t["genes"]): t for t in tables}.values():
        for g in t["genes"]:
            gene_votes[g] = gene_votes.get(g, 0) + 1
    n_sets = len({id(t["genes"]) for t in tables}) or 1
    reference = {g for g, v in gene_votes.items() if v >= 0.5 * n_sets}

    libs = [t["lib"] for t in tables if t["lib"] > 0]
    nz_fracs = [t["nonzero"] / t["n_genes"] for t in tables if t["n_genes"]]
    med_lib = median(libs)
    med_nz = median(nz_fracs)

    # summary metrics, grouped by declared strandedness
    smetrics = {s: summary_metrics(st) for s, st in summaries.items()}
    by_strand = {}
    for s, m in smetrics.items():
        strand = str(sheet.get(s, {}).get("strandSpecific", "0")).strip() or "0"
        by_strand.setdefault(strand, []).append(s)

    def reference_group(sample):
        strand = str(sheet.get(sample, {}).get("strandSpecific", "0")).strip() or "0"
        group = by_strand.get(strand, [])
        if len(group) < min_group:      # too small: compare with everyone
            group = list(smetrics)
        return [g for g in group if g != sample]

    records = []
    for sid in sample_ids:
        flags = []
        rec = {"sampleID": sid,
               "strandSpecific": str(sheet[sid].get("strandSpecific", "")),
               "lib_size": "", "nonzero_frac": "", "geneid_overlap": "",
               "assigned_reads": "", "assigned_frac": "",
               "nofeatures_frac": "", "ambiguity_frac": ""}

        # ---- gene count table checks -------------------------------------
        t = count_tables.get(sid)
        if t is None:
            flags.append("MISSING_COUNTS")
        else:
            rec["lib_size"] = int(t["lib"])
            nz = t["nonzero"] / t["n_genes"] if t["n_genes"] else 0.0
            rec["nonzero_frac"] = round(nz, 4)
            overlap = (len(t["genes"] & reference) / len(reference)
                       if reference else 1.0)
            rec["geneid_overlap"] = round(overlap, 4)
            if t["lib"] <= 0:
                flags.append("ALL_ZERO")
            elif t["lib"] < min_lib:
                flags.append("LOW_LIBSIZE")
            elif med_lib and t["lib"] < 0.1 * med_lib:
                flags.append("LOW_LIBSIZE_REL")
            if overlap < id_overlap_min:
                flags.append("GENEID_MISMATCH")
            if t["n_genes"] and nz < 0.5 * med_nz and "ALL_ZERO" not in flags:
                flags.append("LOW_NONZERO_FRAC")

        # ---- featureCounts summary checks --------------------------------
        m = smetrics.get(sid)
        if m is not None:
            rec["assigned_reads"] = int(m["assigned"])
            rec["assigned_frac"] = round(m["assigned_frac"], 4)
            rec["nofeatures_frac"] = round(m["nofeatures_frac"], 4)
            rec["ambiguity_frac"] = round(m["ambiguity_frac"], 4)
            ref = reference_group(sid)
            if m["assigned"] < min_assigned:
                flags.append("LOW_ASSIGNED_READS")
            if outlier_low(m["assigned_frac"],
                           [smetrics[g]["assigned_frac"] for g in ref], k, 0.10):
                flags.append("LOW_ASSIGNED_FRAC")
            if outlier_high(m["ambiguity_frac"],
                            [smetrics[g]["ambiguity_frac"] for g in ref], k, 0.05):
                flags.append("HIGH_AMBIGUITY")
            if outlier_high(m["nofeatures_frac"],
                            [smetrics[g]["nofeatures_frac"] for g in ref], k, 0.10):
                flags.append("HIGH_NOFEATURES")

        if any(f in FAIL_FLAGS for f in flags):
            status = "FAIL"
        elif flags:
            status = "WARN"
        else:
            status = "PASS"
        rec["flags"] = ",".join(flags) if flags else "-"
        rec["status"] = status
        records.append(rec)

    # ---- group level strandedness hint -----------------------------------
    notes = []
    group_med = {}
    for strand, members in by_strand.items():
        if len(members) >= 3:
            group_med[strand] = (median([smetrics[s]["assigned_frac"] for s in members]),
                                 median([smetrics[s]["ambiguity_frac"] for s in members]),
                                 len(members))
    for strand, (af, amb, n) in sorted(group_med.items()):
        notes.append("strandSpecific=%s: %d samples, median assigned %.2f, "
                     "median ambiguous %.2f" % (strand, n, af, amb))
    if "0" in group_med:
        for strand, (af, amb, n) in group_med.items():
            if strand != "0" and af < 0.75 * group_med["0"][0]:
                notes.append(
                    "CHECK: samples declared strandSpecific=%s have a much lower "
                    "assigned fraction (%.2f) than unstranded ones (%.2f). If "
                    "the library is not really stranded in that direction, set "
                    "strandSpecific in the samplesheet correctly (a wrong "
                    "stranded setting loses ~half of the reads)."
                    % (strand, af, group_med["0"][0]))
    extra = sorted(set(count_tables) - set(sample_ids))
    if extra:
        notes.append("Count columns that are not in the samplesheet: "
                     + ", ".join(extra))
    return records, notes


# --------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------
def sample_from_summary_name(path):
    name = os.path.basename(path)
    return name[:-len(SUMMARY_SUFFIX)] if name.endswith(SUMMARY_SUFFIX) else \
        re.sub(r"\.tsv$", "", name)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--out", default="outrider_qc.tsv")
    p.add_argument("--summaries", nargs="*", default=[],
                   help="featureCounts summary files (<sampleID>%s)" % SUMMARY_SUFFIX)
    p.add_argument("--counts", nargs="+", required=True,
                   help="gene count tables (GeneID + sample column(s))")
    p.add_argument("--stop-on-fail", action="store_true",
                   help="exit with status 1 if any sample FAILs")
    p.add_argument("--min-lib", type=float, default=1e6,
                   help="minimum total gene counts per sample (FAIL below)")
    p.add_argument("--min-assigned", type=float, default=5e6,
                   help="minimum assigned reads in the summary (WARN below)")
    p.add_argument("--mad-k", type=float, default=4.0,
                   help="outlier cutoff in MADs from the group median")
    a = p.parse_args(argv)

    rows = read_samplesheet(a.samplesheet)
    count_tables = {}
    for f in a.counts:
        count_tables.update(read_counts(f))
    summaries = {sample_from_summary_name(f): parse_summary(f)
                 for f in a.summaries}

    records, notes = qc_samples(rows, count_tables, summaries,
                                min_lib=a.min_lib, min_assigned=a.min_assigned,
                                k=a.mad_k)

    cols = ["sampleID", "status", "flags", "strandSpecific", "lib_size",
            "nonzero_frac", "geneid_overlap", "assigned_reads", "assigned_frac",
            "nofeatures_frac", "ambiguity_frac"]
    with open(a.out, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(records)

    n = {s: sum(1 for r in records if r["status"] == s)
         for s in ("PASS", "WARN", "FAIL")}
    print("OUTRIDER input QC: %d PASS, %d WARN, %d FAIL" %
          (n["PASS"], n["WARN"], n["FAIL"]))
    for note in notes:
        print("  " + note)
    for r in records:
        if r["status"] != "PASS":
            print("  %-4s %s: %s" % (r["status"], r["sampleID"], r["flags"]))
    print("Written: " + a.out)

    if a.stop_on_fail and n["FAIL"]:
        print("Stopping: at least one sample FAILED the input QC.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
