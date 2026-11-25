#!/usr/bin/env python3
"""Summarise host mapping per read set.

Two input modes, one output format:
  --bam read_set=x.bam                 aligners that write SAM/BAM (minimap2, winnowmap, ngmlr, lra)
  --clean read_set=in.fastq.gz,out.fastq.gz
                                       host-removal tools that write the reads they keep
                                       (hostile): host reads = in - out

For each read set it reports the host share of reads and bases, the host
share at every MAPQ cut-off in --mapq (BAM mode only; one alignment serves all
cut-offs, so MAPQ is a variable without re-mapping), and the host:non-host
log-ratio. A per-read table lets later stages pair reads across combinations.
"""
import argparse
import math
import subprocess

from fyp_lib import open_text, read_fastq, read_id, write_results


def log_ratio(a, b):
    # +0.5 pseudo-count so empty sets do not give log(0); report it in methods.
    return round(math.log((a + 0.5) / (b + 0.5)), 5)


def share(a, n):
    return round(a / n, 5) if n else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", action="append", default=[], help="read_set=path.bam")
    ap.add_argument("--clean", action="append", default=[], help="read_set=input.fastq.gz,cleaned.fastq.gz")
    ap.add_argument("--mapq", default="0,20,60", help="comma-separated MAPQ cut-offs to report")
    ap.add_argument("--table", default="host_reads.tsv.gz")
    ap.add_argument("--results", default="results.tsv")
    a = ap.parse_args()
    cuts = [int(x) for x in a.mapq.split(",") if x.strip()]

    rows = []
    with open_text(a.table, "wt") as tab:
        tab.write("read_set\tread_id\tlength\tmapped\tmapq\tref\n")
        for spec in a.bam:
            read_set, bam = spec.split("=", 1)
            n = n_map = bases = bases_map = 0
            at_cut = dict.fromkeys(cuts, 0)
            # -F 0x900 drops secondary and supplementary records: one line per read.
            proc = subprocess.Popen(["samtools", "view", "-F", "0x900", bam], stdout=subprocess.PIPE, text=True)
            for line in proc.stdout:
                f = line.split("\t", 10)
                flag, ref, mapq, L = int(f[1]), f[2], int(f[4]), len(f[9])
                mapped = not flag & 4
                n += 1
                bases += L
                if mapped:
                    n_map += 1
                    bases_map += L
                    for c in cuts:
                        if mapq >= c:
                            at_cut[c] += 1
                tab.write(f"{read_set}\t{f[0]}\t{L}\t{int(mapped)}\t{mapq}\t{ref if mapped else ''}\n")
            if proc.wait() != 0:
                raise SystemExit(f"samtools view failed on {bam}")
            rows += [(read_set, "n_reads", n), (read_set, "n_host_reads", n_map),
                     (read_set, "frac_host_reads", share(n_map, n)),
                     (read_set, "frac_host_bases", share(bases_map, bases)),
                     (read_set, "log_ratio_host_nonhost_reads", log_ratio(n_map, n - n_map))]
            rows += [(read_set, f"frac_host_reads_mapq{c}", share(at_cut[c], n)) for c in cuts]
        for spec in a.clean:
            read_set, pair = spec.split("=", 1)
            src, kept_file = pair.split(",", 1)
            kept = {read_id(h) for h, _, _, _ in read_fastq(kept_file)}
            n = n_map = bases = bases_map = 0
            for h, seq, _, _ in read_fastq(src):
                rid, L = read_id(h), len(seq)
                host = rid not in kept
                n += 1
                bases += L
                if host:
                    n_map += 1
                    bases_map += L
                tab.write(f"{read_set}\t{rid}\t{L}\t{int(host)}\t\t\n")
            rows += [(read_set, "n_reads", n), (read_set, "n_host_reads", n_map),
                     (read_set, "frac_host_reads", share(n_map, n)),
                     (read_set, "frac_host_bases", share(bases_map, bases)),
                     (read_set, "log_ratio_host_nonhost_reads", log_ratio(n_map, n - n_map))]
    if not rows:
        raise SystemExit("give --bam or --clean")
    write_results(a.results, rows)


if __name__ == "__main__":
    main()
