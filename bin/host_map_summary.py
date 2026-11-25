#!/usr/bin/env python3
"""Summarise host mapping (minimap2 BAMs) per read set.

For each read set (e.g. passed.bam, failed.bam) it counts primary records only
and reports the host-mapped share of reads and bases, the share at MAPQ >= 20,
and the host:non-host log-ratio (the quantity the multiplicative bias model
works with). It also writes a per-read table so later stages can pair
reads across combinations.

Usage: host_map_summary.py --bam passed=passed.bam --bam failed=failed.bam
"""
import argparse
import math
import subprocess

from fyp_lib import open_text, write_results


def log_ratio(a, b):
    # +0.5 pseudo-count so empty sets do not give log(0); report it in methods.
    return round(math.log((a + 0.5) / (b + 0.5)), 5)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", action="append", required=True, help="read_set=path.bam")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--table", default="host_reads.tsv.gz")
    ap.add_argument("--results", default="results.tsv")
    a = ap.parse_args()

    rows = []
    with open_text(a.table, "wt") as tab:
        tab.write("read_set\tread_id\tlength\tmapped\tmapq\tref\n")
        for spec in a.bam:
            read_set, bam = spec.split("=", 1)
            n = n_map = n_map_q = bases = bases_map = 0
            # -F 0x900 drops secondary and supplementary records: one line per read.
            proc = subprocess.Popen(["samtools", "view", "-F", "0x900", bam], stdout=subprocess.PIPE, text=True)
            for line in proc.stdout:
                f = line.split("\t", 10)
                flag, ref, mapq, L = int(f[1]), f[2], int(f[4]), len(f[9])
                mapped = not flag & 4
                n += 1; bases += L
                if mapped:
                    n_map += 1; bases_map += L
                    if mapq >= a.min_mapq:
                        n_map_q += 1
                tab.write(f"{read_set}\t{f[0]}\t{L}\t{int(mapped)}\t{mapq}\t{ref if mapped else ''}\n")
            if proc.wait() != 0:
                raise SystemExit(f"samtools view failed on {bam}")
            rows += [
                (read_set, "n_reads", n),
                (read_set, "n_host_reads", n_map),
                (read_set, "frac_host_reads", round(n_map / n, 5) if n else ""),
                (read_set, f"frac_host_reads_mapq{a.min_mapq}", round(n_map_q / n, 5) if n else ""),
                (read_set, "frac_host_bases", round(bases_map / bases, 5) if bases else ""),
                (read_set, "log_ratio_host_nonhost_reads", log_ratio(n_map, n - n_map)),
            ]
    write_results(a.results, rows)


if __name__ == "__main__":
    main()
