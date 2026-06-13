#!/usr/bin/env python3
"""
lonely_runner_recursive_metrics_s376.py

codex-2026-05-31 S376

Tournament sessions in this repo treat n as a recursive structure parameter:
class counts, endpoint-transfer ranks, residue defects, and phase changes are
tracked from n to n+1.  This script builds the analogous metric atlas for the
Lonely Runner problem.

For each threshold denominator n, it records five channels:

1. cell-arrangement complexity for the micro-staircase ansatz;
2. the initial-segment unit boundary skeleton;
3. the scalar-ramp puncture moat;
4. one-denominator-gate response;
5. quotient-ladder near-counterexample pressure.

The point is not to prove LRC here.  The point is to expose which quantities
change recursively, and where prime/composite or divisor-layer events create
new behavior.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()
S362 = SourceFileLoader(
    "lonely_runner_bohr_descent_s362",
    str(ROOT / "04-computation" / "lonely_runner_bohr_descent_s362.py"),
).load_module()
S372 = SourceFileLoader(
    "lonely_runner_creative_multiroute_s372",
    str(ROOT / "04-computation" / "lonely_runner_creative_multiroute_s372.py"),
).load_module()


N_MIN = 2
N_MAX = 18
GATE_Q_LIMIT = 32


@dataclass(frozen=True)
class ExactSpeedRow:
    label: str
    speeds: tuple[int, ...]
    classification: str
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    unprotected: int
    first_unprotected: Fraction | None


@dataclass(frozen=True)
class RecursiveRow:
    n: int
    factors: str
    divisor_count: int
    patterns: int
    candidates: int
    break_denominators: int
    unit_count: int
    unit_density: Fraction
    endpoint_count: int
    initial_peel_depth: int
    puncture_missed: int
    puncture_density: Fraction
    puncture_count: int
    puncture_positions: tuple[int, ...]
    puncture_deltas: tuple[int, ...]
    puncture_position_gcds: tuple[int, ...]
    puncture_delta_gcds: tuple[int, ...]
    gate: ExactSpeedRow | None
    ladder: ExactSpeedRow | None


def fmt_frac(x: Fraction | None) -> str:
    if x is None:
        return "-"
    return S356.fmt_frac(x)


def fmt_float(x: Fraction) -> str:
    return f"{float(x):.6f}"


def factorization(n: int) -> str:
    x = n
    parts: list[str] = []
    p = 2
    while p * p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            parts.append(f"{p}^{e}" if e > 1 else str(p))
        p += 1 if p == 2 else 2
    if x > 1:
        parts.append(str(x))
    return "*".join(parts)


def proper_divisors(n: int) -> tuple[int, ...]:
    return tuple(d for d in range(2, n) if n % d == 0)


def totient(n: int) -> int:
    return sum(1 for a in range(1, n) if gcd(a, n) == 1)


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def break_denominator_count(system) -> int:
    denoms = set()
    for lo, hi in system.intervals:
        denoms.add(lo.denominator)
        denoms.add(hi.denominator)
    return len(denoms)


def speed_row(label: str, speeds: tuple[int, ...]) -> ExactSpeedRow:
    summary = S360.summarize(list(speeds))
    return ExactSpeedRow(
        label=label,
        speeds=summary.speeds,
        classification=summary.classification,
        forbidden_length=summary.forbidden_length,
        max_gap=summary.max_gap,
        gap_ratio=summary.max_gap / summary.threshold,
        unprotected=summary.unprotected_count,
        first_unprotected=summary.first_unprotected,
    )


def report_rank(raw_speeds: tuple[int, ...]) -> tuple[Fraction, int, Fraction, tuple[int, ...]]:
    report = S356.report("rank", list(raw_speeds))
    return (
        report.max_gap / report.threshold,
        report.boundary_witness_count,
        1 - report.forbidden_length,
        report.speeds,
    )


def best_gate_response(n: int, q_limit: int = GATE_Q_LIMIT) -> ExactSpeedRow | None:
    """Best replacement of one initial-segment speed by a multiple of n."""

    base = set(range(1, n))
    candidates: list[tuple[tuple[Fraction, int, Fraction, tuple[int, ...]], str, tuple[int, ...]]] = []
    for missing in range(1, n):
        for q in range(1, q_limit + 1):
            speeds = tuple(sorted((base - {missing}) | {n * q}))
            if len(speeds) != n - 1 or not primitive(speeds):
                continue
            rank = report_rank(speeds)
            candidates.append((rank, f"drop={missing},gate={n*q}", speeds))

    if not candidates:
        return None
    _rank, label, speeds = min(candidates, key=lambda row: row[0])
    return speed_row(label, speeds)


def best_quotient_ladder(n: int) -> ExactSpeedRow | None:
    """Best {1} plus d-multiple ladder, where d is a proper divisor of n."""

    candidates: list[tuple[tuple[Fraction, int, Fraction, tuple[int, ...]], str, tuple[int, ...]]] = []
    for d in proper_divisors(n):
        for skip in range(1, n):
            speeds = tuple(sorted({1} | {d * q for q in range(1, n) if q != skip}))
            if len(speeds) != n - 1 or not primitive(speeds):
                continue
            rank = report_rank(speeds)
            candidates.append((rank, f"d={d},skip={skip}", speeds))

    if not candidates:
        return None
    _rank, label, speeds = min(candidates, key=lambda row: row[0])
    return speed_row(label, speeds)


def recursive_row(n: int) -> RecursiveRow:
    system = S372.build_pattern_system(n)
    initial = S362.summarize(list(range(1, n)))
    puncture = S372.puncture_summary(system, radius=1)
    return RecursiveRow(
        n=n,
        factors=factorization(n),
        divisor_count=len(proper_divisors(n)),
        patterns=len(system.patterns),
        candidates=system.candidate_count,
        break_denominators=break_denominator_count(system),
        unit_count=totient(n),
        unit_density=Fraction(totient(n), n - 1),
        endpoint_count=initial.endpoint_count,
        initial_peel_depth=len(initial.peel_layers),
        puncture_missed=puncture.best_missed,
        puncture_density=Fraction(puncture.best_missed, system.candidate_count),
        puncture_count=puncture.best_count,
        puncture_positions=puncture.positions,
        puncture_deltas=puncture.deltas,
        puncture_position_gcds=tuple(sorted({gcd(n, p) for p in puncture.positions})),
        puncture_delta_gcds=tuple(sorted({gcd(n, d) for d in puncture.deltas})),
        gate=best_gate_response(n),
        ladder=best_quotient_ladder(n),
    )


def tuple_text(values: tuple[int, ...], limit: int = 5) -> str:
    if not values:
        return "-"
    shown = ",".join(str(v) for v in values[:limit])
    if len(values) > limit:
        shown += ",..."
    return shown


def speed_summary_text(row: ExactSpeedRow | None) -> str:
    if row is None:
        return "-"
    return (
        f"{row.label}; {row.classification}; "
        f"gap/th={fmt_float(row.gap_ratio)}; "
        f"unprot={row.unprotected}; first={fmt_frac(row.first_unprotected)}"
    )


def print_core_table(rows: list[RecursiveRow]) -> None:
    print("Core recursive LRC metrics")
    print(
        "n  factors  C(n)  dC  denoms  phi  phi/(n-1)  "
        "E0  peel  moat  moat/cand  pos_gcd  delta_gcd"
    )
    prev_patterns = None
    for row in rows:
        d_patterns = "-" if prev_patterns is None else str(row.patterns - prev_patterns)
        prev_patterns = row.patterns
        print(
            f"{row.n:2d} {row.factors:>8s} "
            f"{row.patterns:5d} {d_patterns:>4s} "
            f"{row.break_denominators:6d} "
            f"{row.unit_count:4d} {fmt_float(row.unit_density):>10s} "
            f"{row.endpoint_count:4d} {row.initial_peel_depth:5d} "
            f"{row.puncture_missed:5d} {fmt_float(row.puncture_density):>10s} "
            f"{tuple_text(row.puncture_position_gcds):>7s} "
            f"{tuple_text(row.puncture_delta_gcds):>9s}"
        )
    print()


def print_pressure_table(rows: list[RecursiveRow]) -> None:
    print("Counterexample-pressure channels")
    for row in rows:
        gate = speed_summary_text(row.gate)
        ladder = speed_summary_text(row.ladder)
        print(f"n={row.n:2d} gate:   {gate}")
        print(f"     ladder: {ladder}")
    print()


def print_transition_table(rows: list[RecursiveRow]) -> None:
    print("Transition signatures n -> n+1")
    print("edge  dC  dphi  dmoat  density_shift  event")
    for left, right in zip(rows, rows[1:]):
        d_c = right.patterns - left.patterns
        d_phi = right.unit_count - left.unit_count
        d_moat = right.puncture_missed - left.puncture_missed
        density_shift = right.puncture_density - left.puncture_density
        events = []
        if right.divisor_count == 0:
            events.append("prime")
        else:
            events.append(f"divs={right.divisor_count}")
        if right.puncture_delta_gcds != (1,):
            events.append(f"nonunit_delta={tuple_text(right.puncture_delta_gcds)}")
        if right.ladder is not None and right.ladder.gap_ratio < Fraction(1, 50):
            events.append("tiny_ladder_gap")
        print(
            f"{left.n:2d}->{right.n:<2d} "
            f"{d_c:4d} {d_phi:5d} {d_moat:6d} "
            f"{fmt_float(density_shift):>13s} "
            f"{';'.join(events)}"
        )
    print()


def print_conceptual_synthesis(rows: list[RecursiveRow]) -> None:
    phase_hist = Counter("prime" if row.divisor_count == 0 else "composite" for row in rows)
    nonunit_delta = [row.n for row in rows if row.puncture_delta_gcds != (1,)]
    tiny_ladders = [
        row.n
        for row in rows
        if row.ladder is not None and row.ladder.gap_ratio < Fraction(1, 50)
    ]
    print("Higher-order concepts suggested by the atlas")
    print(
        "  1. C(n) is an arrangement-complexity analogue of tournament class counts; "
        "its increments jump at divisor-rich n."
    )
    print(
        "  2. phi(n)/(n-1) is the unit-boundary skeleton density; primes maximize it, "
        "while composites expose quotient layers."
    )
    print(
        "  3. The scalar-puncture moat is a local curvature metric around the "
        "Dirichlet equality spine."
    )
    print(
        "  4. Nonunit puncture deltas mark divisor resonance. In this range they occur at "
        f"n={nonunit_delta}."
    )
    print(
        "  5. Quotient ladders are the LRC analogue of special tournament families: "
        "they have tiny global gaps but usually large endpoint defects."
    )
    print(
        f"  6. Tiny ladder gaps below 0.02 occur at n={tiny_ladders}; these are the "
        "right places to study endpoint-transfer obstructions."
    )
    print(f"  phase_hist={dict(phase_hist)}")


def main() -> None:
    print("Lonely Runner recursive metric atlas (codex-2026-05-31 S376)")
    print(f"Exact arithmetic, n={N_MIN}..{N_MAX}, gate q<= {GATE_Q_LIMIT}.")
    print()

    rows = [recursive_row(n) for n in range(N_MIN, N_MAX + 1)]
    print_core_table(rows)
    print_pressure_table(rows)
    print_transition_table(rows)
    print_conceptual_synthesis(rows)


if __name__ == "__main__":
    main()
