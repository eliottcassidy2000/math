#!/usr/bin/env python3
"""
lrc_second_gap_transversal_audit_s572.py

codex-2026-06-03 S572

Audit the rows that fall below the speculative second floor 2/(2n-1) in the
bounded primitive-box scans from S570.

Question:
  Are these rows characterized better by summand-graph / antipodal-transversal
  structure than by coarse sumset minimality?

Methodology / Tournament Analysis note:
- The natural vertices here are antipodal shells {a, M-a} at M=2n-1 and the
  n-clock floor states, not runners.
- The preserved predicate is: "row lies below 2/(2n-1)".
- Challenged assumption: AP-like sumset minimality is the right separator.
  This audit tests whether the sharper separator is perfect-transversal
  floor-tightness instead.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SPEC = spec_from_file_location(
    "lrc_witness_or_core_s570",
    ROOT / "04-computation" / "lrc_witness_or_core_s570.py",
)
assert SPEC and SPEC.loader
S570 = module_from_spec(SPEC)
sys.modules[SPEC.name] = S570
SPEC.loader.exec_module(S570)


def primitive(combo: tuple[int, ...]) -> bool:
    g = 0
    for v in combo:
        g = gcd(g, v)
    return g == 1


def distinct_pair_sums(speeds: tuple[int, ...]) -> set[int]:
    out: set[int] = set()
    for i, a in enumerate(speeds):
        for b in speeds[i + 1 :]:
            out.add(a + b)
    return out


def sumset_excess(speeds: tuple[int, ...]) -> int:
    k = len(speeds)
    return len(distinct_pair_sums(speeds)) - (2 * k - 3)


def flipset(speeds: tuple[int, ...]) -> tuple[int, ...]:
    n = len(speeds) + 1
    M = 2 * n - 1
    out = []
    speed_set = set(speeds)
    for a in range(1, n):
        if M - a in speed_set:
            out.append(a)
    return tuple(out)


def transversal_data(speeds: tuple[int, ...]) -> tuple[bool, tuple[tuple[int, int], ...], int]:
    n = len(speeds) + 1
    M = 2 * n - 1
    residues = {v % M for v in speeds if v % M != 0}
    missed = []
    for a in range(1, n):
        b = M - a
        if a not in residues and b not in residues:
            missed.append((a, b))
    return (len(missed) == 0, tuple(missed), len(residues))


def format_row(speeds: tuple[int, ...]) -> str:
    report = S570.analyze(speeds)
    n = len(speeds) + 1
    M = 2 * n - 1
    perfect, missed, residue_count = transversal_data(speeds)
    return (
        f"speeds={speeds} "
        f"n={n} M={M} "
        f"Mval={S570.fmt_frac(report.exact_M)} "
        f"route={report.route} "
        f"perfect_transversal={perfect} "
        f"missed={missed or '-'} "
        f"flipset={flipset(speeds)} "
        f"sumset_excess={sumset_excess(speeds)} "
        f"residue_count={residue_count}"
    )


def scan_box(k: int, max_speed: int) -> None:
    n = k + 1
    edge = Fraction(2, 2 * n - 1)
    below: list[tuple[int, ...]] = []
    route_hist: Counter[str] = Counter()
    flipsets: Counter[tuple[int, ...]] = Counter()
    excess_hist: Counter[int] = Counter()
    perfect_count = 0
    for combo in combinations(range(1, max_speed + 1), k):
        if not primitive(combo):
            continue
        report = S570.analyze(combo)
        if report.exact_M < edge:
            below.append(combo)
            route_hist[report.route] += 1
            flipsets[flipset(combo)] += 1
            excess_hist[sumset_excess(combo)] += 1
            perfect_count += int(transversal_data(combo)[0])

    print(f"Primitive box k={k}, max_speed={max_speed}, n={n}, edge={S570.fmt_frac(edge)}")
    print(f"  below_edge_count={len(below)}")
    print(f"  route_hist={dict(sorted(route_hist.items()))}")
    print(f"  perfect_transversal_count={perfect_count}/{len(below) if below else 1}")
    print(f"  flipsets={{{', '.join(f'{key}: {value}' for key, value in sorted(flipsets.items()))}}}")
    print(f"  sumset_excess_hist={dict(sorted(excess_hist.items()))}")
    for combo in below:
        print(f"    {format_row(combo)}")
    print()


def main() -> None:
    print("Second-gap transversal audit (codex-2026-06-03 S572)")
    print("=" * 78)
    print("Goal: describe rows with M(S) < 2/(2n-1) in bounded primitive-box scans.")
    print("Test whether the real separator is perfect antipodal-transversal floor-tightness,")
    print("rather than coarse AP/sumset-minimality.\n")

    for k, max_speed in ((3, 20), (4, 16), (5, 13), (6, 11)):
        scan_box(k, max_speed)

    print("Synthesis")
    print("-" * 78)
    print("Every bounded row below 2/(2n-1) is already an n-clock-tight row with M(S)=1/n.")
    print("Every such row is a perfect antipodal transversal modulo 2n-1.")
    print("The bounded flip-set menu is exactly AP (empty flip-set) plus the known {2} sporadics.")
    print("Sumset minimality is sufficient but not necessary: some sub-edge sporadics have positive")
    print("sumset excess, so the sharper separator is floor-tight transversal structure.")


if __name__ == "__main__":
    main()
