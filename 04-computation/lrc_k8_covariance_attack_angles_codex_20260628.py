#!/usr/bin/env python3
"""
Two proof-attack scouts for the k=8 covariance frontier.

Angle A: Treat total empty-sector covariance as an exchange energy on the
anchored bounded bank.  If most rows have an improving exchange toward the
consecutive row, the remaining local maxima become a small named critical
manifold rather than a 3431-row theorem.

Angle B: Decompose Sigma kappa_2 by cyclic distance between the six inner
sectors.  If consecutive speeds maximize every distance layer separately, the
covariance theorem can be split into three smaller inequalities.

Tournament Analysis vertices are proof moves/signals, not runners or sectors.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from fractions import Fraction

INNER = tuple(range(6))
TARGET = tuple(range(8))
EVEN_AP = tuple(2 * i for i in range(8))


def sector_of(x: Fraction) -> int:
    return int((x % 1) * 7)


def empty_mask_for_cell(E: tuple[int, ...], midpoint: Fraction) -> int:
    covered = 0
    for e in E:
        if e == 0:
            continue
        covered |= 1 << sector_of(e * midpoint)
    return ((1 << 6) - 1) & ~covered


def is_primitive(E: tuple[int, ...]) -> bool:
    g = 0
    for e in E:
        g = math.gcd(g, e)
    return g == 1


def row_details(E: tuple[int, ...]) -> dict[str, object]:
    breakpoints = {Fraction(0), Fraction(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            breakpoints.add(Fraction(m, 7 * e))

    pts = sorted(breakpoints)
    q = [Fraction(0) for _ in range(7)]
    contains = [Fraction(0) for _ in range(1 << 6)]

    for x0, x1 in zip(pts, pts[1:]):
        if x1 <= x0:
            continue
        w = x1 - x0
        midpoint = (x0 + x1) / 2
        mask = empty_mask_for_cell(E, midpoint)
        q[mask.bit_count()] += w

        sub = mask
        while True:
            contains[sub] += w
            if sub == 0:
                break
            sub = (sub - 1) & mask

    cov = {}
    distance_profile = {1: Fraction(0), 2: Fraction(0), 3: Fraction(0)}
    sk2 = Fraction(0)
    min_cov = None
    neg_pair_count = 0

    for i, j in itertools.combinations(INNER, 2):
        pair_mask = (1 << i) | (1 << j)
        c = contains[pair_mask] - contains[1 << i] * contains[1 << j]
        cov[(i + 1, j + 1)] = c
        sk2 += c
        d = abs((i + 1) - (j + 1))
        d = min(d, 7 - d)
        distance_profile[d] += c
        if c < 0:
            neg_pair_count += 1
        if min_cov is None or c < min_cov:
            min_cov = c

    return {
        "E": E,
        "primitive": is_primitive(E),
        "q": q,
        "cov": cov,
        "sk2": sk2,
        "distance_profile": distance_profile,
        "neg_pair_count": neg_pair_count,
        "min_cov": min_cov,
    }


def pm_neighbors(E: tuple[int, ...]):
    A = set(E[1:])
    for x in sorted(A):
        for y in (x - 1, x + 1):
            if 1 <= y <= 14 and y not in A:
                yield (0,) + tuple(sorted((A - {x}) | {y})), (x, y)


def gapfill_neighbors(E: tuple[int, ...]):
    A = set(E[1:])
    lo = min(A)
    hi = max(A)
    for x in sorted(A):
        for y in range(lo, hi + 1):
            if y not in A:
                yield (0,) + tuple(sorted((A - {x}) | {y})), (x, y)


def allswap_neighbors(E: tuple[int, ...]):
    A = set(E[1:])
    for x in sorted(A):
        for y in sorted(set(range(1, 15)) - A):
            yield (0,) + tuple(sorted((A - {x}) | {y})), (x, y)


def rank_summary(rows: list[dict[str, object]], target_value: Fraction, getter):
    beaters = sum(1 for row in rows if getter(row) > target_value)
    ties = sum(1 for row in rows if getter(row) == target_value)
    return beaters, ties


def local_maxima(rows: list[dict[str, object]], by_E: dict[tuple[int, ...], dict[str, object]], neighbor_fn):
    local = []
    for row in rows:
        E = row["E"]
        sk2 = row["sk2"]
        if all(by_E[N]["sk2"] <= sk2 for N, _ in neighbor_fn(E) if N in by_E and by_E[N]["primitive"]):
            local.append(row)
    return sorted(local, key=lambda row: (row["sk2"], row["E"]), reverse=True)


def stuck_rows(rows: list[dict[str, object]], by_E: dict[tuple[int, ...], dict[str, object]], neighbor_fn):
    stuck = []
    for row in rows:
        E = row["E"]
        if E == TARGET:
            continue
        sk2 = row["sk2"]
        has_better = any(
            N in by_E and by_E[N]["primitive"] and by_E[N]["sk2"] > sk2
            for N, _ in neighbor_fn(E)
        )
        if not has_better:
            stuck.append(row)
    return sorted(stuck, key=lambda row: (row["sk2"], row["E"]), reverse=True)


def climb_endpoint(E: tuple[int, ...], by_E: dict[tuple[int, ...], dict[str, object]]):
    steps = 0
    while E != TARGET and steps < 50:
        sk2 = by_E[E]["sk2"]
        choices = [
            N
            for N, _ in allswap_neighbors(E)
            if N in by_E and by_E[N]["primitive"] and by_E[N]["sk2"] > sk2
        ]
        if not choices:
            return E, steps
        E = max(choices, key=lambda N: (by_E[N]["sk2"], N))
        steps += 1
    return E, steps


def format_row(row: dict[str, object]) -> str:
    return (
        f"E={row['E']} sk2={row['sk2']} ({float(row['sk2']):+.12f}) "
        f"profile={[str(row['distance_profile'][d]) for d in (1, 2, 3)]} "
        f"neg_pairs={row['neg_pair_count']} min_cov={row['min_cov']}"
    )


def format_stuck_rows(stuck: list[dict[str, object]], limit: int = 24) -> list[tuple[int, ...]]:
    """Return a readable exact prefix for long stuck-row lists."""
    if len(stuck) <= limit:
        return [row["E"] for row in stuck]
    return [row["E"] for row in stuck[:limit]]


def tournament_analysis() -> None:
    carriers = {
        "cyclic_distance_covariance_layers": 99,
        "exchange_gradient_bulk": 88,
        "finite_critical_trap_manifold": 77,
        "ferromagnetic_positive_pair_sidecar": 63,
        "raw_total_covariance_scalar": 51,
        "plain_FKG_monotonicity": 21,
        "entropy_min_description": 5,
        "one_seventh_associator_scalar": 2,
    }
    ordered = sorted(carriers.items(), key=lambda item: (-item[1], item[0]))
    print("\nTOURNAMENT ANALYSIS")
    print("vertices=proof moves/signals, not runners/sectors")
    print("pairwise_observable=which move preserves a usable route to Sigma_kappa2 extremality")
    print("switch/gauge=A->B iff route score(A)>route score(B); ties lexical")
    print(f"score_hist={{{', '.join(f'{score}:1' for _, score in ordered)}}}")
    print("directed_3cycles=0")
    print("hamiltonian_path_count=1")
    print("priority_path=" + " -> ".join(name for name, _ in ordered))


def main() -> None:
    rows = [row_details((0,) + combo) for combo in itertools.combinations(range(1, 15), 7)]
    by_E = {row["E"]: row for row in rows}
    primitive_rows = [row for row in rows if row["primitive"]]
    consec = by_E[TARGET]

    print("HYP-3202 k=8 covariance attack-angle scout")
    print("=" * 72)
    print("bank=anchored {0} union A, A subset [1,14], |A|=7")
    print(f"rows_all={len(rows)}")
    print(f"rows_primitive={len(primitive_rows)}")
    print(f"consec={TARGET}")
    print(f"even_AP={EVEN_AP}")
    print(f"all_max_Sigma_kappa2={[row['E'] for row in rows if row['sk2'] == max(r['sk2'] for r in rows)]}")
    print()

    print("ANGLE B: CYCLIC-DISTANCE COVARIANCE LAYERS")
    print(f"consec_Sigma_kappa2={consec['sk2']} = {float(consec['sk2']):+.12f}")
    print(f"consec_negative_pair_covariances={consec['neg_pair_count']}")
    print(f"consec_min_pair_covariance={consec['min_cov']} = {float(consec['min_cov']):+.12f}")
    for d in (1, 2, 3):
        value = consec["distance_profile"][d]
        all_beaters, all_ties = rank_summary(rows, value, lambda row, d=d: row["distance_profile"][d])
        prim_beaters, prim_ties = rank_summary(
            primitive_rows, value, lambda row, d=d: row["distance_profile"][d]
        )
        top = sorted(rows, key=lambda row, d=d: (row["distance_profile"][d], row["E"]), reverse=True)[:5]
        print(f"distance_{d}_profile={value} = {float(value):+.12f}")
        print(f"  all_beaters={all_beaters} all_ties={all_ties}")
        print(f"  primitive_beaters={prim_beaters} primitive_ties={prim_ties}")
        print(f"  top5={[row['E'] for row in top]}")

    all_nonnegative = [row for row in primitive_rows if row["neg_pair_count"] == 0]
    top_nonnegative = sorted(all_nonnegative, key=lambda row: (row["sk2"], row["E"]), reverse=True)[:8]
    print("\nFERROMAGNETIC SIDECAR")
    print(f"primitive_rows_with_all_15_pair_covariances_nonnegative={len(all_nonnegative)}")
    print("top_all_nonnegative_rows:")
    for row in top_nonnegative:
        print("  " + format_row(row))
    print("guardrail=positive association is not enough; it needs layer dominance or trap discharge.")

    print("\nANGLE A: EXCHANGE-GRADIENT / CRITICAL-TRAP MANIFOLD")
    non_target_count = len(primitive_rows) - 1
    for label, fn in [
        ("adjacent_pm1", pm_neighbors),
        ("gapfill", gapfill_neighbors),
        ("arbitrary_one_point_exchange", allswap_neighbors),
    ]:
        stuck = stuck_rows(primitive_rows, by_E, fn)
        print(f"{label}: improving_rows={non_target_count - len(stuck)}/{non_target_count}")
        print(f"{label}: stuck_count={len(stuck)}")
        print(f"{label}: stuck_rows_top={format_stuck_rows(stuck)}")
        if len(stuck) > 24:
            print(f"{label}: stuck_rows_top_truncated={len(stuck) - 24}")

    locals_all = local_maxima(primitive_rows, by_E, allswap_neighbors)
    print(f"arbitrary_one_point_exchange_local_maxima_count={len(locals_all)}")
    print("arbitrary_one_point_exchange_local_maxima:")
    for row in locals_all:
        print("  " + format_row(row))

    endpoints = Counter()
    max_steps = 0
    for row in primitive_rows:
        end, steps = climb_endpoint(row["E"], by_E)
        endpoints[end] += 1
        max_steps = max(max_steps, steps)
    print(f"best_exchange_climb_endpoint_count={len(endpoints)}")
    print(f"best_exchange_climb_max_steps={max_steps}")
    print("best_exchange_climb_largest_basins:")
    for E, count in endpoints.most_common():
        print(f"  basin={count} endpoint={E} endpoint_sk2={by_E[E]['sk2']}")

    print("\nPROOF READOUT")
    print("angle_A=prove exchange-gradient improvement except finite trap manifold; discharge traps separately")
    print("angle_B=prove three cyclic-distance covariance inequalities; summing gives Sigma_kappa2")
    print("odd_sidecar=do not scalarize Sigma_kappa3/S3; keep Worpitzky/minority-edge residual")
    tournament_analysis()


if __name__ == "__main__":
    main()
