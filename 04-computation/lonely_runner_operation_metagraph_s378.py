#!/usr/bin/env python3
"""Operation-metagraph metrics for the Lonely Runner residue systems.

This pass connects two recent threads:

* natural-number operation modes, where addition collapses to a transitive
  completion and multiplication leaves a sparse divisor skeleton; and
* the Lonely Runner micro-staircase systems, where scalar ramps are exact
  blockers but non-scalar punctures expose structured witness cells.

The goal is not to prove LRC here.  It is to define tournament-style landscape
features for the finite residue systems as n changes, complementing the S376
metric tower already stored in lonely_runner_recursive_metrics_s376.py.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S372 = SourceFileLoader(
    "lonely_runner_creative_multiroute_s372",
    str(ROOT / "04-computation" / "lonely_runner_creative_multiroute_s372.py"),
).load_module()
OPS = SourceFileLoader(
    "natural_operation_modes_s366",
    str(ROOT / "04-computation" / "natural_operation_modes_s366.py"),
).load_module()


@dataclass(frozen=True)
class MoatRow:
    n: int
    patterns: int
    candidates: int
    scalar_gcd_hist: tuple[tuple[int, int], ...]
    best_missed: int
    best_count: int
    best_positions: tuple[int, ...]
    best_deltas: tuple[int, ...]
    delta_orders: tuple[int, ...]
    delta_gcds: tuple[int, ...]
    missed_density: Fraction
    ps_min_product: int
    ps_min_seed: tuple[int, ...]


@dataclass(frozen=True)
class RepairRow:
    n: int
    old_missed: int
    coord: int
    old_residue: int
    new_residue: int
    gain_old_misses: int
    old_misses_remaining: int
    new_exposed: int
    new_missed: int
    exposure_ratio: Fraction | None


@dataclass(frozen=True)
class EndpointRow:
    n: int
    family: str
    classification: str
    gap_ratio: Fraction
    boundary: int
    unprotected: int
    peel_depth: int
    core_endpoints: int
    speeds: tuple[int, ...]


def fmt_fraction(value: Fraction | None) -> str:
    if value is None:
        return "inf"
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def order_in_cyclic_group(delta: int, n: int) -> int:
    delta %= n
    if delta == 0:
        return 1
    return n // gcd(delta, n)


def prime_factors(n: int) -> tuple[int, ...]:
    out = []
    d = 2
    m = n
    while d * d <= m:
        if m % d == 0:
            out.append(d)
            while m % d == 0:
                m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        out.append(m)
    return tuple(out)


def best_product_sum_witnesses(max_n: int) -> dict[int, object]:
    by_k = OPS.enumerate_witnesses(max_n)
    return {k: rows[0] for k, rows in by_k.items()}


def moat_rows(n_min: int, n_max: int) -> tuple[MoatRow, ...]:
    ps_minima = best_product_sum_witnesses(n_max)
    rows = []
    for n in range(n_min, n_max + 1):
        system = S372.build_pattern_system(n)
        _all_scalar, scalar_gcd_hist = S372.scalar_audit(system)
        puncture = S372.puncture_summary(system, radius=1)
        ps = ps_minima[n]
        rows.append(
            MoatRow(
                n=n,
                patterns=len(system.patterns),
                candidates=system.candidate_count,
                scalar_gcd_hist=scalar_gcd_hist,
                best_missed=puncture.best_missed,
                best_count=puncture.best_count,
                best_positions=puncture.positions,
                best_deltas=puncture.deltas,
                delta_orders=tuple(sorted({order_in_cyclic_group(d, n) for d in puncture.deltas})),
                delta_gcds=tuple(sorted({gcd(d, n) for d in puncture.deltas})),
                missed_density=Fraction(puncture.best_missed, system.candidate_count),
                ps_min_product=ps.product,
                ps_min_seed=ps.seed,
            )
        )
    return tuple(rows)


def best_nonreverting_repair(system, vector: tuple[int, ...]) -> RepairRow:
    old_mask = S372.blocked_mask(system, vector)
    old_missed_mask = system.all_mask ^ old_mask
    old_missed = old_missed_mask.bit_count()
    best: RepairRow | None = None

    for coord, old_residue in enumerate(vector):
        for new_residue in range(system.n):
            if new_residue == old_residue:
                continue
            trial = list(vector)
            trial[coord] = new_residue
            trial_tuple = tuple(trial)
            new_mask = S372.blocked_mask(system, trial_tuple)
            new_missed_mask = system.all_mask ^ new_mask
            new_missed = new_missed_mask.bit_count()
            if new_missed == 0:
                continue

            gain = (new_mask & old_missed_mask).bit_count()
            old_remaining = old_missed - gain
            new_exposed = new_missed - old_remaining
            ratio = None if gain == 0 else Fraction(new_exposed, gain)
            row = RepairRow(
                n=system.n,
                old_missed=old_missed,
                coord=coord + 1,
                old_residue=old_residue,
                new_residue=new_residue,
                gain_old_misses=gain,
                old_misses_remaining=old_remaining,
                new_exposed=new_exposed,
                new_missed=new_missed,
                exposure_ratio=ratio,
            )
            key = (
                row.gain_old_misses,
                -row.new_missed,
                -row.new_exposed,
                -row.coord,
                -row.new_residue,
            )
            if best is None:
                best = row
            else:
                best_key = (
                    best.gain_old_misses,
                    -best.new_missed,
                    -best.new_exposed,
                    -best.coord,
                    -best.new_residue,
                )
                if key > best_key:
                    best = row

    if best is None:
        raise RuntimeError(f"no non-reverting repair found for n={system.n}")
    return best


def repair_rows(n_values: tuple[int, ...]) -> tuple[RepairRow, ...]:
    rows = []
    for n in n_values:
        system = S372.build_pattern_system(n)
        puncture = S372.puncture_summary(system, radius=1)
        rows.append(best_nonreverting_repair(system, puncture.examples[0]))
    return tuple(rows)


def endpoint_rows(n_values: tuple[int, ...]) -> tuple[EndpointRow, ...]:
    rows = []
    for n in n_values:
        families = [
            ("initial", tuple(range(1, n))),
            ("single_gate", tuple(list(range(1, n - 1)) + [n])),
        ]
        for family, speeds in families:
            summary = S372.classify_speed_set(f"n{n} {family}", speeds)
            rows.append(
                EndpointRow(
                    n=n,
                    family=family,
                    classification=summary.classification,
                    gap_ratio=summary.gap_ratio,
                    boundary=summary.boundary_witnesses,
                    unprotected=summary.unprotected,
                    peel_depth=summary.peel_depth,
                    core_endpoints=summary.core_endpoints,
                    speeds=summary.speeds,
                )
            )
    return tuple(rows)


def print_moat_table(rows: tuple[MoatRow, ...]) -> None:
    print("RECURSIVE SCALAR-PUNCTURE MOAT ATLAS")
    print("=" * 78)
    print("n is the forbidden denominator, with k=n-1 speeds.")
    print("ps_min uses product-sum arity n as a denominator-aligned comparison.")
    print("orders are cyclic orders of the best puncture deltas in Z/nZ.")
    print()
    header = (
        "n pfactors patterns candidates best_miss density orders gcds "
        "positions deltas ps_min_n(seed)"
    )
    print(header)
    print("-" * len(header))
    for row in rows:
        print(
            f"{row.n:2d} {str(prime_factors(row.n)):>9} "
            f"{row.patterns:8d} {row.candidates:10d} "
            f"{row.best_missed:9d} {fmt_fraction(row.missed_density):>9} "
            f"{str(row.delta_orders):>8} {str(row.delta_gcds):>8} "
            f"{str(row.best_positions):>12} {str(row.best_deltas):>18} "
            f"{row.ps_min_product}:{row.ps_min_seed}"
        )
    print()


def print_growth_table(rows: tuple[MoatRow, ...]) -> None:
    print("CONSECUTIVE-N GROWTH LEDGER")
    print("=" * 78)
    print("This is the first recursive diagnostic: how much the finite system changes")
    print("when the LRC denominator grows by one.")
    print()
    header = "n-1->n pattern_x candidate_x moat_x density_delta selected_orders"
    print(header)
    print("-" * len(header))
    for prev, row in zip(rows, rows[1:]):
        pattern_x = Fraction(row.patterns, prev.patterns)
        candidate_x = Fraction(row.candidates, prev.candidates)
        moat_x = Fraction(row.best_missed, prev.best_missed)
        density_delta = row.missed_density - prev.missed_density
        print(
            f"{prev.n:2d}->{row.n:<2d} "
            f"{fmt_fraction(pattern_x):>9} {fmt_fraction(candidate_x):>11} "
            f"{fmt_fraction(moat_x):>7} {fmt_fraction(density_delta):>13} "
            f"{row.delta_orders}"
        )
    print()


def print_repair_table(rows: tuple[RepairRow, ...]) -> None:
    print("ONE-STEP NONREVERTING REPAIR DEFICIT")
    print("=" * 78)
    print("For a best scalar puncture, choose the single coordinate change that covers")
    print("the most old missed cells, excluding the scalar-ramp reversion.")
    print()
    header = (
        "n old_miss move gain old_left new_exposed new_missed "
        "new_exposed/gain"
    )
    print(header)
    print("-" * len(header))
    for row in rows:
        move = f"c{row.coord}:{row.old_residue}->{row.new_residue}"
        print(
            f"{row.n:2d} {row.old_missed:8d} {move:>10} "
            f"{row.gain_old_misses:4d} {row.old_misses_remaining:8d} "
            f"{row.new_exposed:11d} {row.new_missed:10d} "
            f"{fmt_fraction(row.exposure_ratio):>16}"
        )
    print()


def print_endpoint_table(rows: tuple[EndpointRow, ...]) -> None:
    print("ENDPOINT-PROTECTION SAMPLE")
    print("=" * 78)
    print("Initial segments are tight boundary systems. Replacing n-1 by n inserts")
    print("the mandatory multiplicative gate and measures the first endpoint response.")
    print()
    header = "n family       class          gap/thresh boundary unprotected peel core"
    print(header)
    print("-" * len(header))
    for row in rows:
        print(
            f"{row.n:2d} {row.family:<12} {row.classification:<14} "
            f"{fmt_fraction(row.gap_ratio):>10} {row.boundary:8d} "
            f"{row.unprotected:11d} {row.peel_depth:4d} {row.core_endpoints:4d}"
        )
    print()


def print_synthesis() -> None:
    print("SYNTHESIS")
    print("=" * 78)
    print(
        "1. The additive natural-number shadow behaves like the scalar LRC spine: "
        "both are complete/equality structures that must be quotiented before "
        "the interesting sparse geometry appears."
    )
    print(
        "2. The multiplicative shadow and the best LRC puncture deltas both select "
        "divisor/torsion layers.  Composite n values prefer low-order jumps such "
        "as order 2 at n=14 and order 3 at n=15, while prime n has no proper "
        "torsion and the best delta layer is broad."
    )
    print(
        "3. The product-sum defect law is a useful analogy for LRC repair: a local "
        "move can absorb old defect only by creating a new defect package.  The "
        "repair table measures this as new_exposed/gain."
    )
    print(
        "4. A tournament-style LRC metagraph should have nodes = normalized residue "
        "classes, edges = coordinate repairs, height = missed cells, and a second "
        "endpoint-protection height for speed-set searches.  The recursive object "
        "is the evolution of this two-height landscape as n grows."
    )


def main() -> None:
    print("Lonely Runner recursive metrics and natural-operation bridge (S378)")
    print()
    rows = moat_rows(4, 22)
    print_moat_table(rows)
    print_growth_table(rows)
    print_repair_table(repair_rows((10, 12, 14, 15, 18, 20, 21, 22)))
    print_endpoint_table(endpoint_rows((8, 10, 12, 14, 15, 16)))
    print_synthesis()


if __name__ == "__main__":
    main()
