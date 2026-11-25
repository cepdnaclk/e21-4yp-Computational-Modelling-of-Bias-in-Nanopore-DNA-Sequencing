#!/usr/bin/env python3
"""Dataset inventory from the raw FAST5 files (Week 1, Milestone M2).

Reads only HDF5 attributes (no signal decompression), so it is fast and does
not need the vbz plugin. Writes:
  inventory.csv      one row per FAST5 file: MinKNOW folder (pass/fail), reads,
                     size, multi-read format, flow cell, kit, sampling rate, run id
  raw_reads.tsv.gz   one row per read: read_id, source file, MinKNOW class,
                     channel, mux, start time (s), duration (s), end_reason
The per-read table gives channel and time covariates for the bias models and
lets every basecalled read be joined back to its raw-signal metadata.
"""
import argparse
import csv
import glob
import os

import h5py

from fyp_lib import open_text


def attr(group, key, default=""):
    v = group.attrs.get(key, default)
    return v.decode() if isinstance(v, bytes) else v


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fast5", nargs="+", required=True, help="FAST5 files or globs")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    files = sorted({f for g in a.fast5 for f in glob.glob(g)})
    if not files:
        raise SystemExit("No FAST5 files found")
    inv_rows = []
    with open_text(os.path.join(a.outdir, "raw_reads.tsv.gz"), "wt") as rr:
        rr.write("read_id\tfast5_file\tminknow_class\tchannel\tmux\tstart_time_s\tduration_s\tend_reason\n")
        for f in files:
            cls = "pass" if "_pass_" in os.path.basename(f) else "fail" if "_fail_" in os.path.basename(f) else ""
            with h5py.File(f, "r") as h:
                reads = [k for k in h.keys() if k.startswith("read_")]
                multi = bool(reads)
                info = {}
                for i, key in enumerate(reads):
                    r = h[key]
                    ch, raw = r["channel_id"], r["Raw"]
                    rate = float(attr(ch, "sampling_rate", 4000)) or 4000.0
                    if i == 0:
                        ctx, trk = r["context_tags"], r["tracking_id"]
                        info = {
                            "flow_cell_id": attr(trk, "flow_cell_id"),
                            "flow_cell_product_code": attr(trk, "flow_cell_product_code"),
                            "sequencing_kit": attr(ctx, "sequencing_kit"),
                            "experiment_type": attr(ctx, "experiment_type"),
                            "sample_frequency": attr(ctx, "sample_frequency"),
                            "device_id": attr(trk, "device_id"),
                            "run_id": attr(trk, "run_id"),
                            "sample_id": attr(trk, "sample_id"),
                            "exp_start_time": attr(trk, "exp_start_time"),
                            "minknow_version": attr(trk, "version"),
                        }
                    end_reason = raw.attrs.get("end_reason", "")
                    rr.write(
                        f"{attr(raw, 'read_id') or key[5:]}\t{os.path.basename(f)}\t{cls}\t"
                        f"{attr(ch, 'channel_number')}\t{attr(raw, 'start_mux')}\t"
                        f"{int(attr(raw, 'start_time', 0)) / rate:.2f}\t{int(attr(raw, 'duration', 0)) / rate:.3f}\t"
                        f"{int(end_reason) if end_reason != '' else ''}\n"
                    )
            inv_rows.append({"fast5_file": os.path.basename(f), "minknow_class": cls,
                             "n_reads": len(reads), "size_bytes": os.path.getsize(f),
                             "multi_read": multi, **info, "path": os.path.abspath(f)})

    with open(os.path.join(a.outdir, "inventory.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(inv_rows[0].keys()))
        w.writeheader()
        w.writerows(inv_rows)
    print(f"{len(files)} files, {sum(r['n_reads'] for r in inv_rows)} reads")


if __name__ == "__main__":
    main()
