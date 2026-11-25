#!/usr/bin/env python3
"""Turn config/experiment.yaml into design.csv: one row per combination (leaf).

The design decides WHICH combinations exist; Nextflow only executes them.
Each row carries its full folder path under outputs/ plus the paths of its
parent stages, so every upstream stage is computed once and shared:

  basecaller=<b>/model=<m>                                         basecall
    trim=<t>                                                       trim
      qfilter=<meanerr|meanphred>/minq=<q>/minlen=<l>              QC (threshold on mean read Q)
      qfilter=filtlong/keep=<pct>/minlen=<l>                       QC (filtlong keeps the best pct)
        analysis=host_map/aligner=<a>/host_ref=<r>                 leaf
        analysis=classify/classifier=<c>/refdb=<r>[/confidence=<x>] leaf (confidence: kraken2 only)
        analysis=blast/blast_tool=<t>/refdb=<r>                     leaf

refdb is the reference database (its content), shared by classify and blast:
resources.refdbs.<r> lists which classifiers have a database built from it and
where its search database is. classifier x refdb pairs without a database are
not part of the design (e.g. a Kraken2-only prebuilt database).

Modes
  full : every combination of every factor level (default)
  ofat : the baseline, plus each factor changed one at a time with the others at baseline

Usage: design.py config/experiment.yaml --out design.csv [--mode ofat] [--strict]
                 [--select minq=9,10] [--select basecall=guppy-6.5.7:hac] [--list]

--select keeps only the listed levels of one factor for this run (repeatable).
Factor names: basecall, trim, qfilter, minq, keep, minlen, analysis, and the analysis
factors aligner, host_ref, classifier, refdb, confidence, blast_tool
(--select refdb=... applies to both classify and blast).
--select analysis=host_map runs only that analysis; --select analysis=none runs
basecall, trim and QC only. --list prints every combination path.
"""
import argparse
import csv
import hashlib
import itertools
import os
import sys

import yaml

UPSTREAM = ["basecall", "trim", "qfilter", "minq", "keep", "minlen"]
NUMERIC = ("minq", "minlen", "keep", "confidence")
THRESHOLD_QFILTERS = ("meanerr", "meanphred")
ANALYSIS_FACTORS = {
    "host_map": ["aligner", "host_ref"],
    "classify": ["classifier", "refdb", "confidence"],
    "blast": ["blast_tool", "refdb"],
}
COLUMNS = [
    "combo_id", "path", "qc_path", "trim_path", "bc_path",
    "basecaller", "model", "bc_type", "bc_bin", "model_spec", "import_fastq",
    "trim", "qfilter", "minq", "keep", "minlen", "analysis",
    "aligner", "host_ref", "host_ref_fasta",
    "classifier", "refdb", "db_path", "confidence",
    "blast_tool", "blast_db_prefix", "blast_db_fasta",
]
TRIM_LEVELS = {"none", "porechop", "porechop_abi", "porechop_lsk110", "cutadapt_lsk110", "headcrop"}
ALIGNERS = {"minimap2", "winnowmap", "ngmlr", "lra", "hostile"}
CLASSIFIERS = {"kraken2", "krakenuniq", "centrifuge", "centrifuger", "kaiju", "metabuli",
               "ganon", "taxor", "kmcp", "sylph", "diamond", "megan"}
BLAST_TOOLS = {"megablast", "dc-megablast", "blastn", "minimap2", "diamond", "last", "mmseqs2"}


def warn(msg):
    print(f"[design] {msg}", file=sys.stderr)


def fmt_num(v):
    v = float(v)
    return str(int(v)) if v.is_integer() else str(v)


def norm(key, v):
    return fmt_num(v) if key.split(".")[-1] in NUMERIC else str(v)


def search_db(ref):
    """(prefix, fasta) of a refdb's read-search database; fasta defaults to <prefix>.fa."""
    prefix = os.path.expanduser(ref.get("search", "") or "")
    fasta = os.path.expanduser(ref.get("search_fasta", "") or (prefix + ".fa" if prefix else ""))
    return prefix, fasta


def blast_db_ready(tool, ref):
    prefix, fasta = search_db(ref)
    if not prefix:
        return False
    if tool == "minimap2":
        return os.path.exists(fasta)
    if tool == "diamond":
        return os.path.exists(prefix + ".dmnd")
    if tool == "last":
        return os.path.exists(prefix + ".lastdb.prj")
    if tool == "mmseqs2":
        return os.path.exists(prefix + ".mmseqs.dbtype")
    return any(os.path.exists(prefix + ext) for ext in (".nal", ".nin", ".ndb", ".00.nin"))


class Design:
    def __init__(self, cfg, strict, select=None):
        self.cfg, self.strict, self.select = cfg, strict, select or {}
        self.res = cfg.get("resources", {})
        self.dropped = set()

    def check(self, ok, what):
        if ok:
            return True
        if self.strict:
            raise SystemExit(f"[design] missing resource: {what}")
        if what not in self.dropped:
            warn(f"skipping (resource not found yet): {what}")
            self.dropped.add(what)
        return False

    # ---- factor levels -------------------------------------------------
    def levels(self):
        f = self.cfg["factors"]
        lv = {
            "basecall": [str(x) for x in f["basecall"]],
            "trim": [str(x) for x in f["trim"]],
            "qfilter": [str(x) for x in f.get("qfilter", ["meanerr"])],
            "minq": [fmt_num(q) for q in f.get("minq", [9])],
            "keep": [fmt_num(k) for k in f.get("keep", [90])],
            "minlen": [fmt_num(x) for x in f["minlen"]],
        }
        bad = set(lv["trim"]) - TRIM_LEVELS
        if bad:
            raise SystemExit(f"[design] unknown trim levels {sorted(bad)}; known: {sorted(TRIM_LEVELS)}")
        bad = set(lv["qfilter"]) - set(THRESHOLD_QFILTERS) - {"filtlong"}
        if bad:
            raise SystemExit(f"[design] unknown qfilter levels {sorted(bad)}")
        for analysis, keys in ANALYSIS_FACTORS.items():
            spec = (self.cfg.get("analyses") or {}).get(analysis)
            if spec:
                for k in keys:
                    if k in spec:
                        lv[f"{analysis}.{k}"] = [norm(k, v) for v in spec[k]]
        for name, wanted in self.select.items():
            if name == "analysis":
                continue
            keys = [k for k in lv if k == name or k.split(".")[-1] == name]
            if not keys:
                raise SystemExit(f"[design] --select: unknown factor {name!r}; factors: {sorted({k.split('.')[-1] for k in lv})}")
            for k in keys:
                want = [norm(k, w) for w in wanted]
                miss = [w for w in want if w not in lv[k]]
                if miss:
                    raise SystemExit(f"[design] --select {name}: {miss} not in the config levels {lv[k]}")
                lv[k] = [v for v in lv[k] if v in want]
        return lv

    def baseline(self, lv):
        b = dict((self.cfg.get("design") or {}).get("baseline") or {})
        out = {}
        for k, levels in lv.items():
            v = norm(k, b.get(k, b.get(k.split(".")[-1], levels[0])))
            if v not in levels:
                if self.select:          # --select removed the baseline level: use the first selected
                    v = levels[0]
                else:
                    raise SystemExit(f"[design] baseline {k}={v} is not one of the listed levels {levels}")
            out[k] = v
        return out

    # ---- upstream combinations (QC nodes) -----------------------------------
    @staticmethod
    def qc_settings(qfilter, minq, keep):
        """minq applies to threshold filters, keep to filtlong; the other is blank."""
        return (minq, "") if qfilter in THRESHOLD_QFILTERS else ("", keep)

    def upstream_grid(self, lv):
        for bc, trim, qf, minlen in itertools.product(lv["basecall"], lv["trim"], lv["qfilter"], lv["minlen"]):
            for v in (lv["keep"] if qf == "filtlong" else lv["minq"]):
                minq, keep = self.qc_settings(qf, v, v)
                yield dict(basecall=bc, trim=trim, qfilter=qf, minq=minq, keep=keep, minlen=minlen)

    def ofat_upstream(self, lv, b):
        base = {k: b[k] for k in UPSTREAM}
        base["minq"], base["keep"] = self.qc_settings(b["qfilter"], b["minq"], b["keep"])
        out = [dict(base)]
        for k in ("basecall", "trim", "minlen"):
            out += [dict(base, **{k: v}) for v in lv[k] if v != base[k]]
        for qf in lv["qfilter"]:          # another filter method at its baseline setting
            if qf != base["qfilter"]:
                minq, keep = self.qc_settings(qf, b["minq"], b["keep"])
                out.append(dict(base, qfilter=qf, minq=minq, keep=keep))
        qf_thr = base["qfilter"] if base["qfilter"] in THRESHOLD_QFILTERS else None
        if qf_thr:                         # each other Q threshold under the baseline method
            out += [dict(base, minq=q) for q in lv["minq"] if q != base["minq"]]
        if "filtlong" in lv["qfilter"]:    # each keep percentage under filtlong
            out += [dict(base, qfilter="filtlong", minq="", keep=k) for k in lv["keep"]]
        return out

    def upstream_row(self, up):
        basecaller, model = up["basecall"].split(":", 1)
        bc = self.res["basecallers"].get(basecaller)
        if bc is None or model not in bc["models"]:
            raise SystemExit(f"[design] no resources.basecallers entry for {basecaller}:{model}")
        mspec = bc["models"][model]
        r = dict.fromkeys(COLUMNS, "")
        r.update(basecaller=basecaller, model=model, bc_type=bc["type"], bc_bin=os.path.expanduser(bc.get("bin", "")),
                 trim=up["trim"], qfilter=up["qfilter"], minq=up["minq"], keep=up["keep"], minlen=up["minlen"])
        if bc["type"] == "import":
            dirs = mspec["fastq"]
            if not self.check(all(os.path.isdir(d) for d in dirs), f"import fastq dirs for {basecaller}:{model}"):
                return None
            r["import_fastq"] = ";".join(dirs)
        else:
            if not self.check(bool(r["bc_bin"]) and os.path.exists(r["bc_bin"]), f"{basecaller} program {r['bc_bin']}"):
                return None
            r["model_spec"] = str(mspec)
            if bc["type"] in ("dorado", "bonito"):
                mdir = os.path.join(os.path.expanduser(bc.get("models_dir", "")), str(mspec))
                if not self.check(os.path.isdir(mdir), f"{basecaller} model {mdir}"):
                    return None
                r["model_spec"] = mdir
        r["bc_path"] = f"basecaller={basecaller}/model={model}"
        r["trim_path"] = f"{r['bc_path']}/trim={r['trim']}"
        if up["qfilter"] == "filtlong":
            r["qc_path"] = f"{r['trim_path']}/qfilter=filtlong/keep={up['keep']}/minlen={up['minlen']}"
        else:
            r["qc_path"] = f"{r['trim_path']}/qfilter={up['qfilter']}/minq={up['minq']}/minlen={up['minlen']}"
        return r

    # ---- analysis leaves ------------------------------------------------------
    def analysis_leaves(self, lv, analysis, fixed=None):
        """Yield the analysis settings; fixed = {'<analysis>.<factor>': value} pins factors (OFAT)."""
        def levels(k):
            key = f"{analysis}.{k}"
            return [fixed[key]] if fixed and key in fixed else lv.get(key, [])
        if analysis == "none":
            yield {}
        elif analysis == "host_map":
            for a, ref in itertools.product(levels("aligner"), levels("host_ref")):
                yield {"aligner": a, "host_ref": ref}
        elif analysis == "classify":
            for clf, ref in itertools.product(levels("classifier"), levels("refdb")):
                if clf not in self.refdb(ref).get("classifiers", []):
                    continue               # no database of this content for this classifier
                for c in ((levels("confidence") or ["0"]) if clf == "kraken2" else [""]):
                    yield {"classifier": clf, "refdb": ref, "confidence": c}
        elif analysis == "blast":
            for t, ref in itertools.product(levels("blast_tool"), levels("refdb")):
                if self.refdb(ref).get("search"):
                    yield {"blast_tool": t, "refdb": ref}

    def refdb(self, name):
        ref = self.res.get("refdbs", {}).get(name)
        if ref is None:
            raise SystemExit(f"[design] refdb {name!r} has no resources.refdbs entry")
        return ref

    def leaf_row(self, base, analysis, an):
        r = dict(base)
        r["analysis"] = "" if analysis == "none" else analysis
        leaf = ""
        if analysis == "host_map":
            if an["aligner"] not in ALIGNERS:
                raise SystemExit(f"[design] unknown aligner {an['aligner']}; known: {sorted(ALIGNERS)}")
            ref = self.res.get("host_refs", {}).get(an["host_ref"], {})
            if not self.check(os.path.exists(ref.get("fasta", "")), f"host_ref {an['host_ref']}"):
                return None
            r.update(aligner=an["aligner"], host_ref=an["host_ref"], host_ref_fasta=ref["fasta"])
            leaf = f"analysis=host_map/aligner={an['aligner']}/host_ref={an['host_ref']}"
        elif analysis == "classify":
            clf, ref = an["classifier"], an["refdb"]
            if clf not in CLASSIFIERS:
                raise SystemExit(f"[design] unknown classifier {clf!r}; known: {sorted(CLASSIFIERS)}")
            path = os.path.join(os.path.expanduser(self.refdb(ref)["root"]), clf)
            if not self.check(os.path.exists(path), f"refdb {ref} for {clf} ({path})"):
                return None
            r.update(classifier=clf, refdb=ref, db_path=path, confidence=an["confidence"])
            leaf = f"analysis=classify/classifier={clf}/refdb={ref}"
            if an["confidence"] != "":
                leaf += f"/confidence={an['confidence']}"
        elif analysis == "blast":
            if an["blast_tool"] not in BLAST_TOOLS:
                raise SystemExit(f"[design] unknown blast_tool {an['blast_tool']}; known: {sorted(BLAST_TOOLS)}")
            ref = self.refdb(an["refdb"])
            if not self.check(blast_db_ready(an["blast_tool"], ref), f"refdb {an['refdb']} for {an['blast_tool']}"):
                return None
            prefix, fasta = search_db(ref)
            r.update(blast_tool=an["blast_tool"], refdb=an["refdb"], blast_db_prefix=prefix, blast_db_fasta=fasta)
            leaf = f"analysis=blast/blast_tool={an['blast_tool']}/refdb={an['refdb']}"
        r["path"] = f"{r['qc_path']}/{leaf}" if leaf else r["qc_path"]
        r["combo_id"] = hashlib.sha1(r["path"].encode()).hexdigest()[:12]
        return r

    # ---- grids -------------------------------------------------------------
    def build(self, mode):
        lv = self.levels()
        analyses = [a for a in ANALYSIS_FACTORS if any(k.startswith(a + ".") for k in lv)]
        if "analysis" in self.select:
            unknown = set(self.select["analysis"]) - set(ANALYSIS_FACTORS) - {"none"}
            if unknown:
                raise SystemExit(f"[design] --select analysis: {sorted(unknown)} not in {sorted(ANALYSIS_FACTORS) + ['none']}")
            analyses = [a for a in analyses if a in self.select["analysis"]]
            if "none" in self.select["analysis"]:
                analyses.append("none")
        rows = []
        if mode == "full":
            for up in self.upstream_grid(lv):
                base = self.upstream_row(up)
                if base:
                    for analysis in analyses:
                        rows += [self.leaf_row(base, analysis, an) for an in self.analysis_leaves(lv, analysis)]
        else:
            b = self.baseline(lv)
            fixed = {k: v for k, v in b.items() if "." in k}
            ups = self.ofat_upstream(lv, b)
            for up in ups:                 # upstream changes: every analysis at its baseline
                base = self.upstream_row(up)
                if base:
                    for analysis in analyses:
                        rows += [self.leaf_row(base, analysis, an) for an in self.analysis_leaves(lv, analysis, fixed)]
            base_qc = self.upstream_row(ups[0])
            if base_qc:                    # analysis changes: upstream at baseline
                for analysis in [a for a in analyses if a != "none"]:
                    for key in [k for k in lv if k.startswith(analysis + ".")]:
                        for v in lv[key]:
                            if v != b[key]:
                                rows += [self.leaf_row(base_qc, analysis, an)
                                         for an in self.analysis_leaves(lv, analysis, dict(fixed, **{key: v}))]
        uniq = {}
        for r in rows:
            if r:
                uniq[r["path"]] = r
        return list(uniq.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=["full", "ofat"], help="overrides design.mode in the config")
    ap.add_argument("--strict", action="store_true", help="fail instead of skipping missing resources")
    ap.add_argument("--select", action="append", default=[], help="factor=level[,level...] (repeatable)")
    ap.add_argument("--list", action="store_true", help="print every combination path")
    a = ap.parse_args()
    select = {}
    for item in a.select:
        if "=" not in item:
            raise SystemExit(f"--select needs factor=level[,level], got {item!r}")
        k, v = item.split("=", 1)
        select[k.strip()] = [x.strip() for x in v.split(",") if x.strip()]

    cfg = yaml.safe_load(open(a.config))
    mode = a.mode or (cfg.get("design") or {}).get("mode", "full")
    rows = Design(cfg, a.strict, select).build(mode)
    if not rows:
        raise SystemExit("[design] no runnable combinations (check resources and --select)")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["path"]))
    if a.list:
        for r in sorted(rows, key=lambda r: r["path"]):
            print(r["path"])
    count = lambda col: len({r[col] for r in rows})
    warn(f"mode={mode}: {len(rows)} leaf combinations | basecall nodes {count('bc_path')}, "
         f"trim nodes {count('trim_path')}, qc nodes {count('qc_path')}")


if __name__ == "__main__":
    main()
