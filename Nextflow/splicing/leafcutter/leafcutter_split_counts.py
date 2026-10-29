#!/usr/bin/env python3
"""
Split a LeafCutter counts file (<prefix>_perind_numers.counts.gz) into chunks
of whole clusters, so leafcutterMD.R can run as several independent jobs.
(The outlier test of one cluster does not depend on any other cluster.)

Rows are named chr:start:end:clu_N_strand; all rows of a cluster stay together.
Clusters are spread over the chunks by size (largest first) to balance them.
Writes <prefix>_000.counts.gz, <prefix>_001.counts.gz, ... (never an empty one).
Standard library only (Python >= 3.6).
"""
import argparse
import gzip
import sys


def open_text(path, mode="rt"):
    if mode == "rt":
        with open(path, "rb") as fh:
            gz = fh.read(2) == b"\x1f\x8b"
        return gzip.open(path, "rt") if gz else open(path, "rt")
    return gzip.open(path, mode)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("counts")
    p.add_argument("--chunks", type=int, default=20)
    p.add_argument("--prefix", default="chunk")
    a = p.parse_args(argv)

    with open_text(a.counts) as fh:
        header = fh.readline()
        clusters = {}
        order = []
        for line in fh:
            if not line.strip():
                continue
            name = line.split(" ", 1)[0]
            parts = name.split(":")
            if len(parts) < 4:
                print("Unexpected row name: " + name, file=sys.stderr)
                return 1
            clu = parts[3]
            if clu not in clusters:
                clusters[clu] = []
                order.append(clu)
            clusters[clu].append(line)

    if not clusters:
        print("No clusters in " + a.counts, file=sys.stderr)
        return 1

    n = max(1, min(a.chunks, len(clusters)))
    chunks = [[] for _ in range(n)]
    sizes = [0] * n
    for clu in sorted(order, key=lambda c: -len(clusters[c])):
        i = sizes.index(min(sizes))
        chunks[i].append(clu)
        sizes[i] += len(clusters[clu])

    for i, members in enumerate(chunks):
        name = "%s_%03d.counts.gz" % (a.prefix, i)
        with gzip.open(name, "wt") as out:
            out.write(header)
            mset = set(members)
            for clu in order:           # keep the original order inside a chunk
                if clu in mset:
                    out.writelines(clusters[clu])
    print("%d clusters (%d introns) in %d chunks" %
          (len(clusters), sum(sizes), n))
    return 0


if __name__ == "__main__":
    sys.exit(main())
