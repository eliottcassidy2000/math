#!/usr/bin/env python3
"""HYP-3470: AP84 phase-color grid bridge for LRC14.

This returns to the older HYP-2593/HYP-2594 phase-color reservoir and asks how
it sees the canonical AP-tail family.  It is the exact CRT-placement sidecar
under HYP-3459's broader AP84 color-packet legality audit, complementary to
HYP-3460's phase-branch color pullback, and downstream of HYP-3461's
colored-extension gate carrier.

    S_m = {1,2,...,11,13,84m}.

In HYP-2593's notation this is

    S = P union {V-e : e in E},  P={1,...,11,13}, E={0}, V=84m.

The continuous color reservoir collapses: color 0 is dead and the other
thirteen colors have identical layers.  The nontrivial signal is therefore not
layer mass.  It is the shifted-grid discrepancy and the closed-boundary hits.

This script records the exact fixed live layer, validates the closed color-grid
count against the original CRT predicate on sample rows, and separates three
clocks:

  * HYP-3456 component escapes have period 35.
  * Total phase-color CRT counts have affine period 385.
  * Individual live colors only become uniformly affine after period 5005.

Tournament Analysis declaration.
  Vertex set: proof carriers for the color-grid bridge, not runners.
  Pairwise observable: retained CRT predicate, color-sidecar payload,
    boundary clock, compatibility with HYP-3456, and scalarization loss.
  Gauge/tie path: higher retained proof payload wins; ties follow the bridge
    route from exact CRT predicate to AP-tail splice.

Assumption challenge.
  Considered vertices: runners, phase colors, color residues, fixed live
  intervals, interval endpoints, boundary hits, AP-tail m residues, HYP-3456
  high gaps, Fourier modes, and proof obligations.  The chosen quotient
  preserves the exact q=14V CRT placement predicate and color-grid discrepancy,
  but destroys branch-cover geometry and component adjacency.  The challenged
  assumption is that the phase-color reservoir's mass should see the AP-tail
  mod-35 clock; for E={0}, mass is blind and the clock moves to boundary/grid
  discrepancy.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "04-computation") not in sys.path:
    sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_colored_discrepancy_bound_codex as disc
import lrc14_phase_color_reservoir_codex as pc


P = tuple(list(range(1, 12)) + [13])
E = (0,)
LIVE_COLORS = tuple(range(1, 14))


def fmt(value: F | int) -> str:
    if isinstance(value, int):
        return str(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def interval_fmt(interval: tuple[F, F]) -> str:
    return f"[{fmt(interval[0])},{fmt(interval[1])}]"


def floor_frac(x: F) -> int:
    return x.numerator // x.denominator


def ceil_frac(x: F) -> int:
    return -((-x.numerator) // x.denominator)


def closed_grid_count(intervals: list[tuple[F, F]], b: int, V: int) -> int:
    """Count t with (b+14t)/(14V) in the closed interval union."""
    total = 0
    for lo, hi in intervals:
        lower = V * lo - F(b, 14)
        upper = V * hi - F(b, 14)
        t0 = max(ceil_frac(lower), 0)
        t1 = min(floor_frac(upper), V - 1)
        if t0 <= t1:
            total += t1 - t0 + 1
    return total


def layer_rows() -> list[dict]:
    return disc.layer_data(P, E)


def color_counts(m: int, rows: list[dict] | None = None) -> list[int]:
    rows = layer_rows() if rows is None else rows
    V = 84 * m
    return [closed_grid_count(row["intervals"], row["b"], V) for row in rows]


def open_color_counts(m: int, rows: list[dict] | None = None) -> list[int]:
    rows = layer_rows() if rows is None else rows
    V = 84 * m
    return [disc.open_grid_count(row["intervals"], row["b"], V) for row in rows]


def total_count(m: int, rows: list[dict] | None = None) -> int:
    return sum(color_counts(m, rows))


def open_total_count(m: int, rows: list[dict] | None = None) -> int:
    return sum(open_color_counts(m, rows))


def direct_crt_count(m: int) -> int:
    V = 84 * m
    return pc.actual_crt_count(pc.build_S(P, E, V), V)


def expected_total(m: int) -> F:
    return F(5112 * m, 385)


def expected_live_color(m: int) -> F:
    return F(5112 * m, 5005)


def residual_total(m: int, rows: list[dict] | None = None) -> int:
    return 5112 * m - 385 * total_count(m, rows)


def direct_validation(rows: list[dict]) -> list[tuple[int, int, int, int, int, F]]:
    out = []
    for m in (1, 2, 3, 4, 5, 7, 12, 35, 36, 70):
        closed = total_count(m, rows)
        open_count = open_total_count(m, rows)
        direct = direct_crt_count(m)
        out.append((m, closed, open_count, direct, closed - open_count, expected_total(m) - closed))
    return out


def total_period_failures(rows: list[dict]) -> list[tuple[int, int]]:
    failures = []
    for m in range(1, 771):
        shift = total_count(m + 385, rows) - total_count(m, rows)
        if shift != 5112:
            failures.append((m, shift))
    return failures


def color_period_failures(rows: list[dict]) -> list[tuple[int, tuple[int, ...]]]:
    failures = []
    for m in range(1, 101):
        before = color_counts(m, rows)
        after = color_counts(m + 5005, rows)
        shift = tuple(after[b] - before[b] for b in range(14))
        expected = (0,) + (5112,) * 13
        if shift != expected:
            failures.append((m, shift))
    return failures


def boundary_bonus(m: int, rows: list[dict]) -> int:
    return total_count(m, rows) - open_total_count(m, rows)


def boundary_rule_failures(rows: list[dict]) -> list[tuple[int, int]]:
    failures = []
    for m in range(1, 771):
        predicted = 0 if m % 7 == 0 else 2
        actual = boundary_bonus(m, rows)
        if actual != predicted:
            failures.append((m, actual))
    return failures


def mirror_failures(rows: list[dict]) -> list[tuple[int, int, int, int]]:
    failures = []
    for m in range(1, 386):
        counts = color_counts(m, rows)
        for b in range(1, 7):
            if counts[b] != counts[14 - b]:
                failures.append((m, b, counts[b], counts[14 - b]))
    return failures


def live_color_discrepancy(rows: list[dict]) -> dict[int, tuple[F, F, F]]:
    out = {}
    for b in LIVE_COLORS:
        deviations = [expected_live_color(m) - color_counts(m, rows)[b] for m in range(1, 5006)]
        out[b] = (max(abs(value) for value in deviations), max(deviations), min(deviations))
    return out


def count_directed_3cycles(path: list[str]) -> int:
    rank = {name: index for index, name in enumerate(path)}
    total = 0
    for a, b, c in combinations(path, 3):
        ab = rank[a] < rank[b]
        bc = rank[b] < rank[c]
        ca = rank[c] < rank[a]
        if ab == bc == ca:
            total += 1
    return total


def tournament_fingerprint() -> tuple[dict[int, int], list[str], int]:
    vertices = {
        "exact_H2593_CRT_color_predicate": (10, 10, 10, 10, 10, 10),
        "closed_color_grid_floor_formula": (10, 10, 10, 9, 9, 10),
        "total_period385_discrepancy_clock": (9, 9, 9, 9, 10, 9),
        "color_vector_period5005_sidecar": (8, 8, 10, 8, 10, 8),
        "mod7_boundary_bonus_gate": (8, 8, 8, 8, 8, 8),
        "H3456_period35_component_bridge": (7, 8, 7, 9, 8, 8),
        "raw_live_layer_mass_scalar": (3, 4, 2, 3, 3, 2),
    }
    scores = {name: sum(values) for name, values in vertices.items()}
    hist = dict(sorted(Counter(scores.values()).items()))
    path = [name for name, _score in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    return hist, path, count_directed_3cycles(path)


def main() -> None:
    rows = layer_rows()
    sigma, K, _ = disc.sigma_K(P, E)
    live_intervals = rows[1]["intervals"]
    live_measure = rows[1]["measure"]
    total_residuals = [residual_total(m, rows) for m in range(1, 386)]
    actual_first35 = [total_count(m, rows) for m in range(1, 36)]
    shift35 = [total_count(m + 35, rows) - total_count(m, rows) for m in range(1, 36)]
    residual_first35 = [residual_total(m, rows) for m in range(1, 36)]
    bonus_hist = dict(sorted(Counter(boundary_bonus(m, rows) for m in range(1, 386)).items()))
    color_shift385 = tuple(color_counts(386, rows)[b] - color_counts(1, rows)[b] for b in range(14))
    color_dev = live_color_discrepancy(rows)
    hist, path, cycles = tournament_fingerprint()

    print("HYP-3470 AP84 PHASE-COLOR GRID BRIDGE")
    print("=" * 72)
    print("status=EVIDENCE / exact AP-tail color-grid discrepancy bridge; not an LRC14 proof")
    print("source=HYP-2593/HYP-2594 reservoir + HYP-3459 color packet + HYP-3460 pullback + HYP-3461 extension carrier + AP84 sidecars")
    print()

    print("A. Phase-color reservoir collapse")
    print("  family=S_m={1,2,...,11,13,84m}")
    print("  H2593 form: P={1,...,11,13}, E={0}, V=84m")
    print(f"  sigma={fmt(sigma)}")
    print(f"  K={K}")
    print(f"  live_color_measure={fmt(live_measure)}")
    print("  color0_measure=0; colors 1..13 have the same fixed live layer")
    print("  live_intervals:")
    for interval in live_intervals:
        print(f"    {interval_fmt(interval)} length={fmt(interval[1] - interval[0])}")
    print("  expected_total_hits=V*sigma=(5112/385)*m")
    print("  expected_live_color_hits=(5112/5005)*m")
    print()

    print("B. Closed grid count validation")
    print("  closed_count uses ceil/floor on fixed live intervals and counts boundary hits.")
    print("  sample: m | closed | open | direct_CRT | boundary_bonus | expected-closed")
    failures = []
    for m, closed, open_count, direct, bonus, deficit in direct_validation(rows):
        if closed != direct:
            failures.append(m)
        print(
            f"  {m:3d} | {closed:6d} | {open_count:4d} | {direct:10d} | "
            f"{bonus:14d} | {fmt(deficit)}"
        )
    print(f"  direct_validation_failures={failures}")
    print()

    print("C. Discrepancy clocks")
    print(f"  actual_first35={actual_first35}")
    print(f"  shift_by_35_first35={shift35}")
    print("  total affine period: A(m+385)-A(m)=5112")
    print(f"  total_period385_failures_m_1_to_770={total_period_failures(rows)}")
    print(
        "  residual R(m)=5112m-385A(m): "
        f"min={min(total_residuals)} at m={[i + 1 for i, v in enumerate(total_residuals) if v == min(total_residuals)]}, "
        f"max={max(total_residuals)} at m={[i + 1 for i, v in enumerate(total_residuals) if v == max(total_residuals)]}"
    )
    print(f"  residual_first35={residual_first35}")
    print()

    print("D. Color-sidecar and boundary resonance")
    print(f"  mirror_color_failures_m_1_to_385={mirror_failures(rows)}")
    print(f"  color_shift_m1_to_m386={color_shift385}")
    print("  individual live colors have uniform affine period 5005:")
    print(f"  color_period5005_failures_m_1_to_100={color_period_failures(rows)}")
    print(f"  boundary_bonus_hist_m_1_to_385={bonus_hist}")
    print("  boundary bonus rule: bonus=0 iff 7 divides m, otherwise bonus=2")
    print(f"  boundary_bonus_rule_failures_m_1_to_770={boundary_rule_failures(rows)}")
    print("  live color discrepancy over one color period: b | max_abs | max | min")
    for b in LIVE_COLORS:
        max_abs, max_dev, min_dev = color_dev[b]
        print(f"  {b:2d} | {fmt(max_abs):>10} | {fmt(max_dev):>10} | {fmt(min_dev):>10}")
    print()

    print("E. Proof pull")
    print("  The old phase-color exactness is useful here, but the continuous mass")
    print("  does not see the AP84 mod-35 component clock: E={0} kills color 0 and")
    print("  leaves thirteen identical live layers.  The AP-tail coloring signal is")
    print("  the shifted grid and boundary sidecar.  HYP-3456's period-35 escape")
    print("  clock remains a component/corridor count; HYP-3459 records the legal")
    print("  AP84 color-packet quotient, HYP-3460 records the phase-branch")
    print("  pullback, and HYP-3461 records the colored-extension carrier,")
    print("  with HYP-3458 retaining the")
    print("  35-state coloring-recursion packet.  HYP-3470 supplies the exact q=14V")
    print("  color-grid placement clock that should be carried if the AP-tail splice")
    print("  uses CRT witnesses rather than branch-union components.")
    print("  Next theorem target: prove the fixed four-interval color-grid formula")
    print("  symbolically and attach its period-385/5005 sidecar to the HYP-3439")
    print("  AP-tail descent only when actual q=14V witness placement is needed.")
    print()

    print("F. Tournament Analysis")
    print("  vertices=color-grid proof carriers, not runners or raw arcs")
    print("  pairwise_observable=CRT exactness + grid formula + discrepancy period + color sidecar + HYP-3456 compatibility")
    print("  switch=higher retained proof payload; ties by AP-tail color-grid bridge route")
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={cycles}")
    print("  hamiltonian_path=" + " -> ".join(path))
    print()

    print("G. Assumption Challenge")
    print("  Considered vertices: runners, phase colors, color residues, fixed live")
    print("  intervals, interval endpoints, boundary hits, AP-tail m residues,")
    print("  HYP-3456 high gaps, Fourier modes, and proof obligations.  The chosen")
    print("  quotient preserves the exact q=14V CRT placement predicate and color-grid")
    print("  discrepancy.  It destroys branch-cover geometry and component adjacency.")
    print("  The challenged assumption is that reservoir mass should detect the AP-tail")
    print("  mod-35 clock; in this AP84 specialization, mass is blind and the useful")
    print("  color information is the boundary/grid discrepancy sidecar.")


if __name__ == "__main__":
    main()
