#!/usr/bin/env python3
"""Split reads into passed / failed / artefact by mean read Q and length.

Unlike NanoFilt, this keeps every read: failed reads are written out, not
discarded, because the core Phase 1 test is what the QC filter throws away.

  artefact : length > --maxlen (the >50 kb, Q~3 artefact reads; excluded from both sets)
  passed   : mean_q >= --minq and length >= --minlen
  failed   : everything else

Outputs (in --outdir): passed.fastq.gz, failed.fastq.gz, artefact.fastq.gz,
reads_table.tsv.gz (read_id, length, mean_q, gc, qc_class) and results.tsv.
results.tsv reports each read set, the overall pass rate, and the pass rate
per GC bin (a direct check of GC-dependent filtering).
"""
import argparse
import os

from fyp_lib import gc_bin, gc_percent, mean_q, open_text, read_fastq, read_id, summarise, write_results

CLASSES = ("passed", "failed", "artefact")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--minq", type=float, required=True)
    ap.add_argument("--minlen", type=int, required=True)
    ap.add_argument("--maxlen", type=int, default=50000)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    outs = {c: open_text(os.path.join(a.outdir, f"{c}.fastq.gz"), "wt") for c in CLASSES}
    stats = {c: ([], [], []) for c in CLASSES}
    gc_counts = {}  # gc bin -> [passed, failed]

    with open_text(os.path.join(a.outdir, "reads_table.tsv.gz"), "wt") as table:
        table.write("read_id\tlength\tmean_q\tgc\tqc_class\n")
        for header, seq, plus, qual in read_fastq(a.reads):
            L, q, g = len(seq), mean_q(qual), gc_percent(seq)
            if L > a.maxlen:
                cls = "artefact"
            elif q >= a.minq and L >= a.minlen:
                cls = "passed"
            else:
                cls = "failed"
            outs[cls].write(f"{header}\n{seq}\n{plus}\n{qual}\n")
            for lst, v in zip(stats[cls], (L, q, g)):
                lst.append(v)
            if cls != "artefact":
                c = gc_counts.setdefault(gc_bin(g), [0, 0])
                c[0 if cls == "passed" else 1] += 1
            table.write(f"{read_id(header)}\t{L}\t{q:.3f}\t{g:.2f}\t{cls}\n")
    for fh in outs.values():
        fh.close()

    rows = []
    for cls in CLASSES:
        rows += [(cls, k, v) for k, v in summarise(*stats[cls]).items()]
    n_pass, n_fail = len(stats["passed"][0]), len(stats["failed"][0])
    b_pass, b_fail = sum(stats["passed"][0]), sum(stats["failed"][0])
    rows.append(("all", "pass_rate_reads", round(n_pass / max(n_pass + n_fail, 1), 5)))
    rows.append(("all", "pass_rate_bases", round(b_pass / max(b_pass + b_fail, 1), 5)))
    for b, (p, f) in sorted(gc_counts.items()):
        rows.append(("all", f"pass_rate_gc_{b}", round(p / max(p + f, 1), 5)))
        rows.append(("all", f"n_reads_gc_{b}", p + f))
    write_results(os.path.join(a.outdir, "results.tsv"), rows)


if __name__ == "__main__":
    main()
