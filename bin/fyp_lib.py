"""Shared helpers for the FYP pipeline scripts (standard library only).

Every script in bin/ imports this module, so they all read FASTQ, compute
mean read quality and GC, and write results.tsv in exactly the same way.
That matters: if two stages computed "mean Q" differently, a change in the
pass rate could come from the formula and not from the factor being tested.
"""
import gzip
import math
import statistics
from typing import Dict, Iterator, List, Tuple

# Mean read Q is the Phred value of the mean per-base error probability,
# the same definition NanoFilt and chopper use. It is lower than the mean
# of the per-base Q values, and is applied identically to every basecaller.
MEANQ_DEFINITION = "-10*log10(mean(10^(-Q_base/10))) over all bases"

GC_BINS = [(0, 30), (30, 40), (40, 50), (50, 60), (60, 70), (70, 101)]


def open_text(path: str, mode: str = "rt"):
    """Open plain or gzip-compressed text based on the file name."""
    if path.endswith(".gz"):
        return gzip.open(path, mode, compresslevel=1) if "w" in mode else gzip.open(path, mode)
    return open(path, mode)


def read_fastq(path: str) -> Iterator[Tuple[str, str, str, str]]:
    """Yield (header_line, seq, plus_line, qual) records, without newlines."""
    with open_text(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().rstrip("\n")
            plus = fh.readline().rstrip("\n")
            qual = fh.readline().rstrip("\n")
            if not header.startswith("@") or len(seq) != len(qual):
                raise ValueError(f"Malformed FASTQ record in {path}: {header[:60]!r}")
            yield header.rstrip("\n"), seq, plus, qual


def read_id(header: str) -> str:
    return header[1:].split()[0]


# Pre-computed error probability for each Phred+33 character.
_ERR = [10 ** (-(i - 33) / 10) if i >= 33 else 1.0 for i in range(256)]


def mean_q(qual: str) -> float:
    if not qual:
        return 0.0
    err = sum(_ERR[c] for c in qual.encode("ascii")) / len(qual)
    return -10 * math.log10(err) if err > 0 else 60.0


MEANPHRED_DEFINITION = "arithmetic mean of the per-base Phred Q values"


def mean_phred(qual: str) -> float:
    """Arithmetic mean of per-base Q. Always >= mean_q(); some tools report this one."""
    if not qual:
        return 0.0
    return sum(qual.encode("ascii")) / len(qual) - 33


def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    s = seq.upper()
    return 100.0 * (s.count("G") + s.count("C")) / len(s)


def gc_bin(gc: float) -> str:
    for lo, hi in GC_BINS:
        if lo <= gc < hi:
            return f"{lo}-{min(hi, 100)}"
    return "70-100"


def n50(lengths: List[int]) -> int:
    if not lengths:
        return 0
    total, run = sum(lengths), 0
    for L in sorted(lengths, reverse=True):
        run += L
        if run * 2 >= total:
            return L
    return 0


def summarise(lengths: List[int], quals: List[float], gcs: List[float]) -> Dict[str, float]:
    """Standard read-set metrics written to results.tsv by several stages."""
    n = len(lengths)
    out = {"n_reads": n, "n_bases": sum(lengths)}
    if n:
        out.update({
            "mean_length": round(sum(lengths) / n, 2),
            "median_length": statistics.median(lengths),
            "n50": n50(lengths),
            "mean_q": round(sum(quals) / n, 3),
            "median_q": round(statistics.median(quals), 3),
            "mean_gc": round(sum(gcs) / n, 3),
            "median_gc": round(statistics.median(gcs), 3),
        })
    return out


def write_results(path: str, rows: List[Tuple[str, str, object]]) -> None:
    """Write results.tsv in the long format every stage shares:
    read_set <tab> metric <tab> value."""
    with open(path, "w") as fh:
        fh.write("read_set\tmetric\tvalue\n")
        for read_set, metric, value in rows:
            fh.write(f"{read_set}\t{metric}\t{value}\n")


def parse_combo_path(path: str) -> Dict[str, str]:
    """'basecaller=guppy/model=hac/...' -> {'basecaller': 'guppy', 'model': 'hac', ...}"""
    factors = {}
    for part in path.strip("/").split("/"):
        if "=" in part:
            k, v = part.split("=", 1)
            factors[k] = v
    return factors
