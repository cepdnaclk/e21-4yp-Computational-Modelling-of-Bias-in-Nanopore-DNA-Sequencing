#!/usr/bin/env python3
"""Gather every node of outputs/ into analysis-ready tables.

Walks <outdir> for metadata.json files (skipping _shared, _runs, _summary),
reads the sibling results.tsv, and writes to <outdir>/_summary/:
  nodes.tsv        one row per node: combo_id, stage, path, parent, one column per factor
  results_long.tsv one row per node x read_set x metric, with the factor columns
  results_long.parquet  the same, for pandas/pyarrow
Factor columns are empty where a factor does not apply (e.g. host_ref on a
classify leaf), so the long table can be filtered by stage and pivoted freely.
"""
import argparse
import json
import os

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    nodes, long_rows = [], []
    for root, dirs, files in os.walk(a.outdir):
        dirs[:] = sorted(d for d in dirs if not d.startswith("_"))
        if "metadata.json" not in files:
            continue
        meta = json.load(open(os.path.join(root, "metadata.json")))
        base = {"combo_id": meta["combo_id"], "stage": meta["stage"], "path": meta["path"],
                "parent_path": meta.get("parent_path", ""), **meta["factors"]}
        nodes.append({**base, "n_outputs": len(meta.get("outputs", {})), "created_utc": meta.get("created_utc", "")})
        res = os.path.join(root, "results.tsv")
        if os.path.exists(res):
            df = pd.read_csv(res, sep="\t", dtype={"value": str})
            for r in df.itertuples(index=False):
                long_rows.append({**base, "read_set": r.read_set, "metric": r.metric, "value": r.value})

    if not nodes:
        raise SystemExit(f"No metadata.json found under {a.outdir}")
    out = os.path.join(a.outdir, "_summary")
    os.makedirs(out, exist_ok=True)
    nodes_df = pd.DataFrame(nodes)
    long_df = pd.DataFrame(long_rows)
    long_df["value_num"] = pd.to_numeric(long_df["value"], errors="coerce")
    nodes_df.to_csv(os.path.join(out, "nodes.tsv"), sep="\t", index=False)
    long_df.to_csv(os.path.join(out, "results_long.tsv"), sep="\t", index=False)
    long_df.to_parquet(os.path.join(out, "results_long.parquet"), index=False)
    print(f"{len(nodes_df)} nodes ({nodes_df['stage'].value_counts().to_dict()}), {len(long_df)} result rows -> {out}")


if __name__ == "__main__":
    main()
