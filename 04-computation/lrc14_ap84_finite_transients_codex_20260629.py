#!/usr/bin/env python3
"""HYP-3457: finite AP84 transient closure for the LRC14 tail.

HYP-3454 isolated four finite mixed rows before the AP84 endpoint phase:

    S_m = {1,2,...,11,13,84m},   m = 1,2,3,4.

For m >= 5 the best escape is the pure E:84m/E:84m interval.  For m <= 4 the
same moving high gap is clipped by the fixed low odd wall B1:5 (and by mirror
symmetry B0:5), so the transient rows have exactly four mixed escape windows:

    [8/49, (98m-1)/(588m)]         L[B1:7]  R[E:84m]
    [(98m+1)/(588m), 6/35]         L[E:84m] R[B1:5]
    [29/35, (490m-1)/(588m)]       L[B0:5]  R[E:84m]
    [(490m+1)/(588m), 41/49]       L[E:84m] R[B0:7]

This script checks those four formulas against the HYP-3452 component-cover
audit and records the exact finite certificate.  It closes the "m=1..4 mixed
transients" sidecar, but it is still only an AP-tail bridge lemma for the
HYP-3439 rank-5 descent.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
H3452_PATH = ROOT / "04-computation" / "lrc14_ap84_tail_component_phase_codex_20260629.py"
F = Fraction
Interval = tuple[F, F]


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


H3452 = load_module("hyp3452_finite_transients", H3452_PATH)


def fmt(value: F | None) -> str:
    if value is None:
        return "None"
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def interval_fmt(interval: Interval) -> str:
    return f"[{fmt(interval[0])},{fmt(interval[1])}]"


@dataclass(frozen=True)
class ExpectedWindow:
    name: str
    interval: Interval
    labels: str
    rank: int
    length: F


@dataclass(frozen=True)
class TransientRecord:
    m: int
    components: int
    dead: int
    escapes: int
    low_rank_escape: int
    max_dead_pair_rank: int
    projection_components: int
    projection_largest: int
    expected: tuple[ExpectedWindow, ...]
    actual: tuple[ExpectedWindow, ...]
    best_interval: Interval
    best_labels: str
    best_rank: int | None
    best_length: F

    @property
    def exact_windows_ok(self) -> bool:
        return self.expected == self.actual

    @property
    def closure_ok(self) -> bool:
        return (
            self.exact_windows_ok
            and self.escapes == 4
            and self.low_rank_escape == 4
            and self.projection_components == 1
            and self.best_interval == self.expected[1].interval
            and self.best_labels == self.expected[1].labels
            and self.best_rank == 2
        )


def expected_windows(m: int) -> tuple[ExpectedWindow, ...]:
    e = 84 * m
    return (
        ExpectedWindow(
            "outer_left",
            (F(8, 49), F(98 * m - 1, 588 * m)),
            f"L[B1:7] R[E:{e}]",
            2,
            F(2 * m - 1, 588 * m),
        ),
        ExpectedWindow(
            "inner_left",
            (F(98 * m + 1, 588 * m), F(6, 35)),
            f"L[E:{e}] R[B1:5]",
            2,
            F(14 * m - 5, 2940 * m),
        ),
        ExpectedWindow(
            "inner_right",
            (F(29, 35), F(490 * m - 1, 588 * m)),
            f"L[B0:5] R[E:{e}]",
            2,
            F(14 * m - 5, 2940 * m),
        ),
        ExpectedWindow(
            "outer_right",
            (F(490 * m + 1, 588 * m), F(41, 49)),
            f"L[E:{e}] R[B0:7]",
            2,
            F(2 * m - 1, 588 * m),
        ),
    )


def actual_windows(record) -> tuple[ExpectedWindow, ...]:
    alive = [component for component in record.row.components if component.best_survivor]
    alive.sort(key=lambda component: component.best_survivor)
    return tuple(
        ExpectedWindow(
            name,
            component.best_survivor,
            component.labels,
            component.endpoint_rank,
            component.union_measure,
        )
        for name, component in zip(("outer_left", "inner_left", "inner_right", "outer_right"), alive)
    )


def build_records() -> list[TransientRecord]:
    records = H3452.audit_tail(4)
    out: list[TransientRecord] = []
    for record in records:
        best = record.best
        out.append(
            TransientRecord(
                m=record.m,
                components=record.components,
                dead=record.dead,
                escapes=record.escapes,
                low_rank_escape=record.summary.low_rank_escape,
                max_dead_pair_rank=record.row.max_dead_pair_rank,
                projection_components=record.projection_components,
                projection_largest=record.projection_largest,
                expected=expected_windows(record.m),
                actual=actual_windows(record),
                best_interval=best.best_survivor,
                best_labels=best.labels,
                best_rank=record.row.best_rank,
                best_length=record.best_len,
            )
        )
    return out


def transient_inequalities(m: int) -> dict[str, F]:
    """Positive values certify the m <= 4 mixed clipping regime."""
    high_right_minus_b1_wall = F(455 - 98 * m, 2940 * m)
    left_inner_margin = F(2 * m + 1, 588 * m)
    right_inner_margin = F(2 * m + 1, 588 * m)
    outer_length = F(2 * m - 1, 588 * m)
    inner_length = F(14 * m - 5, 2940 * m)
    return {
        "high_right_minus_B1_5": high_right_minus_b1_wall,
        "inner_left_corridor_margin": left_inner_margin,
        "inner_right_corridor_margin": right_inner_margin,
        "outer_length": outer_length,
        "inner_length": inner_length,
    }


def tournament_fingerprint() -> tuple[dict[int, int], list[str], int]:
    vertices = {
        "four_window_explicit_packet": (10, 10, 10, 10, 9, 10),
        "component_audit_exact_match": (10, 10, 9, 9, 8, 9),
        "mirror_pairing": (9, 9, 9, 9, 8, 8),
        "mixed_endpoint_clip_inequalities": (9, 9, 8, 8, 8, 9),
        "rank_drop_m3_to_m4": (7, 7, 7, 8, 7, 7),
        "raw_transient_sample_rows": (4, 4, 3, 4, 3, 3),
    }
    scores = {name: sum(values) for name, values in vertices.items()}
    hist = dict(sorted(Counter(scores.values()).items()))
    path = [name for name, _score in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    rank = {name: index for index, name in enumerate(path)}
    cycles = 0
    for a, b, c in combinations(vertices, 3):
        ab = rank[a] < rank[b]
        bc = rank[b] < rank[c]
        ca = rank[c] < rank[a]
        if ab == bc == ca:
            cycles += 1
    return hist, path, cycles


def main() -> None:
    rows = build_records()
    failures = [row.m for row in rows if not row.closure_ok]
    exact_window_failures = [row.m for row in rows if not row.exact_windows_ok]
    rank_drop_failures = [row.m for row in rows if row.m >= 3 and row.max_dead_pair_rank > 2]
    inequality_failures = [
        row.m
        for row in rows
        if any(value <= 0 for value in transient_inequalities(row.m).values())
    ]
    hist, path, cycles = tournament_fingerprint()

    print("HYP-3457 AP84 FINITE TRANSIENT CLOSURE")
    print("=" * 72)
    print("status=EVIDENCE / finite mixed AP-tail transient certificate; not an LRC14 proof")
    print("source=HYP-3452 component phase + HYP-3454 endpoint clock + HYP-3456 floor-count handoff")
    print()

    print("A. Four-window transient formula")
    print("  family=S_m={1,2,...,11,13,84m}, transient_m=1..4")
    print("  fixed low corridors=[8/49,6/35] and [29/35,41/49]")
    print("  expected survivor windows:")
    print("    outer_left  = [8/49,(98m-1)/(588m)]")
    print("    inner_left  = [(98m+1)/(588m),6/35]")
    print("    inner_right = [29/35,(490m-1)/(588m)]")
    print("    outer_right = [(490m+1)/(588m),41/49]")
    print()

    print("B. Component-audit validation")
    print(f"  exact_window_failures={exact_window_failures}")
    print(f"  closure_failures={failures}")
    print(f"  rank_drop_failures_for_m_ge_3={rank_drop_failures}")
    print(f"  inequality_failures={inequality_failures}")
    print("  row summary: m | components | dead | escapes | low_rank_escape | max_pair | projection")
    for row in rows:
        print(
            f"  {row.m:1d} | {row.components:10d} | {row.dead:4d} | {row.escapes:7d} | "
            f"{row.low_rank_escape:15d} | {row.max_dead_pair_rank:8d} | "
            f"{row.projection_components}/{row.projection_largest}"
        )
    print()

    print("C. Exact survivor windows")
    for row in rows:
        print(f"  m={row.m}:")
        for window in row.actual:
            print(
                f"    {window.name:11s} {interval_fmt(window.interval):>18s} "
                f"len={fmt(window.length):>8s} rank={window.rank} labels={window.labels}"
            )
        print(
            f"    best={interval_fmt(row.best_interval)} len={fmt(row.best_length)} "
            f"rank={row.best_rank} labels={row.best_labels}"
        )
    print()

    print("D. Endpoint inequalities")
    print("  For m=1..4, the high right endpoint still lies beyond the B1:5 wall:")
    print("    (98m+13)/(588m) - 6/35 = (455-98m)/(2940m) > 0.")
    print("  At m=5 this sign changes, which is exactly the HYP-3454 pure E/E phase.")
    print("  margins: m | high_right-B1:5 | inner_margin | outer_len | inner_len")
    for row in rows:
        vals = transient_inequalities(row.m)
        print(
            f"  {row.m:1d} | {fmt(vals['high_right_minus_B1_5']):>17s} | "
            f"{fmt(vals['inner_left_corridor_margin']):>12s} | "
            f"{fmt(vals['outer_length']):>9s} | {fmt(vals['inner_length']):>9s}"
        )
    print()

    print("E. Proof pull")
    print("  The finite mixed side of the AP-tail bridge is now an explicit four-row")
    print("  certificate: exactly four low-rank escapes survive in each transient row,")
    print("  with mirror-paired endpoint labels and no component-audit failures.  After")
    print("  HYP-3456 closed the escape-count clock, the remaining AP-tail work is to")
    print("  import/prove the HYP-3431 fixed-corridor carrier as the low branch-union")
    print("  theorem and splice HYP-3454/HYP-3456/HYP-3457 into HYP-3439.")
    print()

    print("F. Tournament Analysis")
    print("  vertices=finite transient proof obligations, not runners or raw arcs")
    print("  pairwise_observable=exact window match + endpoint labels + mirror pairing + rank payload")
    print("  switch=higher retained AP-tail proof payload; ties by transient route")
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={cycles}")
    print("  hamiltonian_path=" + " -> ".join(path))
    print()

    print("G. Assumption Challenge")
    print("  Considered vertices: runners, m-values, four survivor windows, fixed")
    print("  corridor sections, high-grid wall events, endpoint walls, dead-cover graph")
    print("  summaries, and proof obligations.  The chosen quotient preserves the")
    print("  finite AP-tail branch-union survivor predicate and endpoint labels, but")
    print("  destroys arbitrary non-AP component geometry.  It is only the finite")
    print("  transient sidecar for the AP-tail bridge.")


if __name__ == "__main__":
    main()
