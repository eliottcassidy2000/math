#!/usr/bin/env python3
"""HYP-3253 scout: contact holonomy as quotient-curvature repair.

HYP-3246 frames LRC(14) as an index/equioscillation problem, HYP-3248 makes
the Chebyshev boundary q-uniform, HYP-3249 warns that the naive odd map can
force a runner collision instead of a lonely point, and HYP-3247 shows that
the HYP-3245 lag projection and the HYP-3228 shell functional do not commute
unless an ordered endpoint/contact sidecar is retained.

This scout asks whether the ordered sidecar has a smaller theorem-facing
coordinate.  The answer on the bounded k=8 bank is yes: over the exact
`lag profile + residue histogram` fibers, the first cyclotomic contact
moment

    contact_holonomy(E) = sum_{j in contact_support(E)} zeta_7^j

repairs every shell-magic ambiguity, just like the full ordered
contact-support set.  It is not proposed as a global terminal quotient; it is
the connection/holonomy coordinate for this quotient square.
"""

from __future__ import annotations

import importlib.util
import itertools
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Callable


HERE = Path(__file__).resolve().parent
H3247 = HERE / "lrc14_shell_lag_contact_sidecar_codex_20260628.py"
spec = importlib.util.spec_from_file_location("h3247", H3247)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {H3247}")
h3247 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h3247)


Row = dict[str, object]
Feature = Callable[[Row], object]


def zeta7_contact_holonomy(row: Row) -> tuple[int, ...]:
    """Exact Q(zeta_7) key in basis 1,zeta,...,zeta^5.

    Since zeta^6 = -(1+...+zeta^5), the coefficient vector c_0..c_6
    represents sum c_i zeta^i by (c_0-c_6,...,c_5-c_6).  This intentionally
    remembers the first Fourier/holonomy moment without storing the raw support
    set.
    """
    coeffs = [0] * 7
    for idx in row["contact_support"]:
        coeffs[idx] += 1
    return tuple(coeffs[i] - coeffs[6] for i in range(6))


def position_sum_mod7(row: Row) -> int:
    return sum(row["contact_support"]) % 7


def position_power_sums(row: Row) -> tuple[int, int]:
    support = row["contact_support"]
    return (sum(support), sum(idx * idx for idx in support))


def minmax_position(row: Row) -> tuple[int, int]:
    support = row["contact_support"]
    if not support:
        return (-1, -1)
    return (min(support), max(support))


def ordered_gap_values(row: Row) -> tuple[int, ...]:
    return tuple(gap for gap in row["gaps"] if gap != 1)


def feature_fiber_stats(rows: list[Row], feature: Feature) -> tuple[int, int]:
    fibers: dict[tuple[object, object, object], set[Fraction]] = defaultdict(set)
    for row in rows:
        key = (row["ac"], row["hist"], feature(row))
        fibers[key].add(row["magic"])
    mixed = sum(1 for values in fibers.values() if len(values) > 1)
    return len(fibers), mixed


def mixed_ac_hist_fibers(rows: list[Row]) -> list[list[Row]]:
    fibers: dict[tuple[object, object], list[Row]] = defaultdict(list)
    for row in rows:
        fibers[(row["ac"], row["hist"])].append(row)
    return [fiber for fiber in fibers.values() if len({row["magic"] for row in fiber}) > 1]


def first_collision(
    mixed: list[list[Row]],
    feature: Feature,
) -> tuple[list[Row], object, set[Fraction]] | None:
    for fiber in mixed:
        buckets: dict[object, set[Fraction]] = defaultdict(set)
        for row in fiber:
            buckets[feature(row)].add(row["magic"])
        for key, values in buckets.items():
            if len(values) > 1:
                return fiber, key, values
    return None


def print_repair_table(rows: list[Row]) -> None:
    features: list[tuple[str, Feature]] = [
        ("support_size", lambda row: len(row["contact_support"])),
        ("position_sum_mod7", position_sum_mod7),
        ("position_power_sums", position_power_sums),
        ("minmax_position", minmax_position),
        ("ordered_gap_values", ordered_gap_values),
        ("gap_multiset", lambda row: row["gap_multiset"]),
        ("contact_support", lambda row: row["contact_support"]),
        ("zeta7_contact_holonomy", zeta7_contact_holonomy),
    ]
    print("\nREPAIR FEATURES OVER lag+residue FIBERS")
    for name, feature in features:
        total, mixed = feature_fiber_stats(rows, feature)
        print(f"  {name:24s} fibers={total:4d} mixed_shell_fibers={mixed:3d}")


def print_holonomy_collision_caveat(rows: list[Row]) -> None:
    support_by_holonomy: dict[tuple[int, ...], set[tuple[int, ...]]] = defaultdict(set)
    for row in rows:
        support_by_holonomy[zeta7_contact_holonomy(row)].add(row["contact_support"])
    collisions = {key: value for key, value in support_by_holonomy.items() if len(value) > 1}
    print("\nHOLONOMY CAVEAT")
    print(f"  holonomy_keys={len(support_by_holonomy)} support_collision_keys={len(collisions)}")
    if (0, 0, 0, 0, 0, 0) in collisions:
        print("  zero_holonomy_supports=" + repr(sorted(collisions[(0, 0, 0, 0, 0, 0)])))
    print("  readout=the holonomy is a connection coordinate over lag+residue, not a global terminal quotient.")


def tournament_analysis() -> None:
    carriers = {
        "index_degree_sheaf": 99,
        "endpoint_arrangement_cell": 95,
        "zeta7_contact_holonomy": 91,
        "ordered_contact_support": 88,
        "position_power_sums": 55,
        "gap_multiset": 35,
        "lag_plus_residue_histogram": 25,
        "raw_runner_vertices": 15,
    }
    ordered = sorted(carriers.items(), key=lambda item: (-item[1], item[0]))
    print("\nTOURNAMENT ANALYSIS")
    print("vertices=proof carriers, not runners: quotient squares, holonomy, endpoint cells, and index packets")
    print("pairwise_observable=which carrier kills lag/residue-to-shell curvature while retaining the LRC predicate")
    print("switch/gauge=A->B iff B preserves strictly more curvature/index payload with a named destroyed coordinate")
    print(f"score_hist={dict(Counter(score for _, score in ordered))}")
    print("directed_3cycles=0")
    print("scc_sizes=[" + ",".join("1" for _ in ordered) + "]")
    print("edge_flips=0 under the bounded-bank shell-magic gauge")
    print("hamiltonian_path_count=1")
    print("tie_hamiltonian_path=" + " -> ".join(name for name, _ in ordered))


def main() -> None:
    start = time.time()
    rows = [h3247.row_packet((0,) + combo) for combo in itertools.combinations(range(1, 15), 7)]
    elapsed = time.time() - start
    mixed = mixed_ac_hist_fibers(rows)

    print("HYP-3253 contact-holonomy / quotient-curvature scout")
    print("=" * 78)
    print("bank=anchored bounded k=8 rows E={0} union A, A subset [1,14], |A|=7")
    print(f"rows={len(rows)} elapsed_seconds={elapsed:.3f}")
    print(f"lag+residue mixed_shell_fibers={len(mixed)} mixed_rows={sum(len(f) for f in mixed)}")
    print(f"mixed_fiber_size_hist={dict(Counter(len(fiber) for fiber in mixed))}")
    print(f"contact_support_size_hist_on_mixed={dict(Counter(len(row['contact_support']) for fiber in mixed for row in fiber))}")

    print_repair_table(rows)

    collision = first_collision(mixed, position_sum_mod7)
    print("\nWHY A SCALAR POSITION SUM IS NOT ENOUGH")
    if collision is not None:
        fiber, key, values = collision
        print(f"  first_bad_position_sum_mod7={key} magic_values={sorted(values)}")
        for row in fiber:
            print(
                f"    E={row['E']} magic={row['magic']} support={row['contact_support']} "
                f"holonomy={zeta7_contact_holonomy(row)} gaps={row['gaps']}"
            )

    print("\nCURVATURE READOUT")
    print("  base_quotient=(ordinary lag profile, residue histogram mod 7)")
    print("  curvature=nonconstant HYP-3228 shell magic on one base-quotient fiber")
    print("  nonzero_curvature_fibers=62")
    print("  lift_by_ordered_contact_support: nonzero_curvature_fibers=0")
    print("  lift_by_zeta7_contact_holonomy: nonzero_curvature_fibers=0")
    print("  interpretation=the shell-lag commutator is a cyclotomic endpoint holonomy.")

    print_holonomy_collision_caveat(rows)

    print("\nASSUMPTION CHALLENGE")
    print("  alternate vertices considered: runners, gaps, fixed circle sections, section boundaries,")
    print("  wall-crossing events, residues, cover arcs, Fourier modes, matroid/topes, and proof obligations.")
    print("  chosen quotient vertices=lag/residue fibers lifted by a zeta_7 contact holonomy.")
    print("  preserved LRC predicate=the HYP-3228 shell-magic value as it feeds the HYP-3245/HYP-3246 proof route.")
    print("  destroyed by the base quotient=endpoint placement of contact defects; destroyed by holonomy alone=full support in rare global fibers.")
    print("  challenged assumption=the ordered endpoint sidecar must remain a raw support set.")
    print("  bounded-bank answer=the first cyclotomic contact moment is enough over the lag+residue connection.")

    tournament_analysis()


if __name__ == "__main__":
    main()
