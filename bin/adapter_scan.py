#!/usr/bin/env python3
"""Count reads that still carry the LSK110 ligation adapter (lab protocol, Method PDF).

Ligation adapter (LA): top strand 5'-CCTGTACTTCGTTCAGTTACGTATTGCT-3' starts each read;
bottom strand 5'-GCAATACGTAACTGAACGAAGTACAGG-3' follows the insert at the read end.
Each read is searched (edlib, edit distance <= --max-error x adapter length) for:
  5p        top strand within the first --window bases
  3p        bottom strand within the last --window bases
  internal  either strand or its reverse complement between the two windows
            (adapter in the middle of a read = chimera / unsplit read)
On the untrimmed basecalls (trim=none) these are the adapter content; on a trimmed
read set they are what the trimmer left behind.

Appends rows to results.tsv (read_set=all):
  frac_reads_lsk110_5p, frac_reads_lsk110_3p, frac_reads_lsk110_any_end, frac_reads_lsk110_internal
"""
import argparse

import edlib

from fyp_lib import read_fastq

COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def rc(s):
    return s.translate(COMP)[::-1]


def found(adapter, seq, k):
    return bool(seq) and edlib.align(adapter, seq, mode="HW", task="distance", k=k)["editDistance"] >= 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--bottom", required=True)
    ap.add_argument("--window", type=int, default=150)
    ap.add_argument("--max-error", type=float, default=0.2)
    ap.add_argument("--results", default="results.tsv", help="appended to (header written if new)")
    a = ap.parse_args()

    top, bottom = a.top.upper(), a.bottom.upper()
    k_top, k_bot = int(a.max_error * len(top)), int(a.max_error * len(bottom))
    internal = [(top, k_top), (rc(top), k_top), (bottom, k_bot), (rc(bottom), k_bot)]
    n = n5 = n3 = nend = nint = 0
    w = a.window
    for _, seq, _, _ in read_fastq(a.reads):
        seq = seq.upper()
        n += 1
        s5 = found(top, seq[:w], k_top)
        s3 = found(bottom, seq[-w:], k_bot)
        n5 += s5
        n3 += s3
        nend += s5 or s3
        mid = seq[w:-w] if len(seq) > 2 * w else ""
        nint += any(found(ad, mid, k) for ad, k in internal)

    d = max(n, 1)
    rows = [("frac_reads_lsk110_5p", n5 / d), ("frac_reads_lsk110_3p", n3 / d),
            ("frac_reads_lsk110_any_end", nend / d), ("frac_reads_lsk110_internal", nint / d)]
    try:
        new = open(a.results).read() == ""
    except FileNotFoundError:
        new = True
    with open(a.results, "a") as fh:
        if new:
            fh.write("read_set\tmetric\tvalue\n")
        for m, v in rows:
            fh.write(f"all\t{m}\t{v:.5f}\n")
    print(f"adapter_scan: {n} reads; 5' {n5}, 3' {n3}, internal {nint}")


if __name__ == "__main__":
    main()
