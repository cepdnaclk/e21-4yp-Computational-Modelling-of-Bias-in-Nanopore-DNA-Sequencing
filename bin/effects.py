#!/usr/bin/env python3
"""Estimate factor effects on one metric from outputs/_summary/results_long.tsv.

Two views of the same stored grid:

1. OFAT contrasts: for each factor, vary it while every other factor stays at
   the baseline (from config/experiment.yaml). Prints each level's value and
   its difference from baseline. Easy to read; blind to interactions.

2. Factorial model (needs a full design): ordinary least squares on the metric
   (or its logit, for shares) with every factor as a categorical main effect,
   plus optional two-way interactions. The tools are deterministic and there is
   one sequencing run, so this has no replicate error; treat the coefficients
   as effect sizes and get uncertainty from read-level resampling (bootstrap)
   rather than from these p-values.

Example:
  effects.py --summary outputs/_summary/results_long.tsv --config config/experiment.yaml \
             --stage qc --read-set all --metric pass_rate_reads --logit --interactions
"""
import argparse
import itertools

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import yaml

UPSTREAM = ["basecaller", "model", "trim", "minq", "minlen"]
ANALYSIS = ["aligner", "host_ref", "classifier", "db", "blast_tool", "blast_db"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--stage", required=True)
    ap.add_argument("--metric", required=True)
    ap.add_argument("--read-set", default="all")
    ap.add_argument("--logit", action="store_true", help="model logit(metric); use for shares in (0,1)")
    ap.add_argument("--interactions", action="store_true", help="add all two-way interactions")
    a = ap.parse_args()

    df = pd.read_csv(a.summary, sep="\t", dtype=str)
    df = df[(df.stage == a.stage) & (df.read_set == a.read_set) & (df.metric == a.metric)].copy()
    if df.empty:
        raise SystemExit("No rows match --stage/--read-set/--metric")
    df["y"] = pd.to_numeric(df["value"], errors="coerce")
    factors = [f for f in UPSTREAM + ANALYSIS if f in df and df[f].nunique(dropna=True) > 1]

    cfg = yaml.safe_load(open(a.config))
    base_cfg = (cfg.get("design") or {}).get("baseline") or {}
    bc, model = base_cfg.get("basecall", ":").split(":")
    baseline = {"basecaller": bc, "model": model, **{k: str(v) for k, v in base_cfg.items() if k != "basecall"}}
    for f in factors:  # analysis factors default to their first listed level
        baseline.setdefault(f, sorted(df[f].dropna().unique())[0])

    # ---- 1) OFAT contrasts --------------------------------------------------
    print(f"\n== OFAT contrasts: {a.metric} ({a.read_set}), stage={a.stage}")
    base_mask = np.logical_and.reduce([df[f] == baseline[f] for f in factors])
    base_val = df.loc[base_mask, "y"]
    if len(base_val) != 1:
        print(f"  baseline combination not found once (found {len(base_val)}); OFAT skipped")
    else:
        b = base_val.iloc[0]
        print(f"  baseline {baseline} -> {b:.4f}")
        for f in factors:
            others = [g for g in factors if g != f]
            m = np.logical_and.reduce([df[g] == baseline[g] for g in others]) if others else np.ones(len(df), bool)
            sub = df.loc[m, [f, "y"]].sort_values(f)
            for lvl, y in sub.itertuples(index=False):
                print(f"  {f:<11} {lvl:<20} {y:.4f}   diff {y - b:+.4f}")

    # ---- 2) Factorial model -------------------------------------------------
    d = df.dropna(subset=["y"]).copy()
    if a.logit:
        eps = 1e-4
        d["y"] = np.log(d["y"].clip(eps, 1 - eps) / (1 - d["y"].clip(eps, 1 - eps)))
    terms = [f"C({f}, Treatment('{baseline[f]}'))" for f in factors]
    if a.interactions:
        terms += [f"C({x}):C({y})" for x, y in itertools.combinations(factors, 2)]
    formula = "y ~ " + (" + ".join(terms) if terms else "1")
    n_cells = len(d)
    fit = smf.ols(formula, data=d).fit()
    print(f"\n== Factorial model ({'logit ' if a.logit else ''}{a.metric}); {n_cells} combinations, "
          f"{int(fit.df_model)} effect parameters, {int(fit.df_resid)} residual df")
    print(fit.params.round(4).to_string())
    if fit.df_resid <= 0:
        print("  (saturated: no residual df; drop --interactions or add read-level bootstrap)")


if __name__ == "__main__":
    main()
