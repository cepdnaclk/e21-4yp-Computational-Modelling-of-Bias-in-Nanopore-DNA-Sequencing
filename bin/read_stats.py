#!/usr/bin/env python3
"""Per-read table and summary metrics for one FASTQ (used after basecalling and trimming).

Outputs (in --outdir):
  reads_table.tsv.gz  read_id, length, mean_q, gc
  results.tsv         read_set=all, with n_reads, n_bases, N50, mean/median Q and GC,
                      plus the share of reads in each GC bin
"""
import argparse
import os

from fyp_lib import gc_bin, gc_percent, mean_q, open_text, read_fastq, read_id, summarise, write_results


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--read-set", default="all")
    a = ap.parse_args()

    lengths, quals, gcs, bins = [], [], [], {}
    with open_text(os.path.join(a.outdir, "reads_table.tsv.gz"), "wt") as out:
        out.write("read_id\tlength\tmean_q\tgc\n")
        for header, seq, _, qual in read_fastq(a.reads):
            L, q, g = len(seq), mean_q(qual), gc_percent(seq)
            lengths.append(L); quals.append(q); gcs.append(g)
            bins[gc_bin(g)] = bins.get(gc_bin(g), 0) + 1
            out.write(f"{read_id(header)}\t{L}\t{q:.3f}\t{g:.2f}\n")

    rows = [(a.read_set, k, v) for k, v in summarise(lengths, quals, gcs).items()]
    n = max(len(lengths), 1)
    rows += [(a.read_set, f"frac_reads_gc_{b}", round(c / n, 5)) for b, c in sorted(bins.items())]
    write_results(os.path.join(a.outdir, "results.tsv"), rows)


if __name__ == "__main__":
    main()
