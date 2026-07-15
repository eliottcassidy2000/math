#!/usr/bin/env python3
"""Exact scout for the disconnected-return satellite frontier.

This artifact proves/checks three facts about

    closure(R_U) = closure({d : max_u ||u*d|| < 2/143}).

* Every return component occupies one tooth of B=max(U).  In scaled tooth
  coordinate e=B*d-k it is one explicitly computable interval cut out by the
  balanced residues [u*k]_B.  This is an exact O(B*|U|) classifier.
* If N_R is the number of surviving cells, THM-803's endpoint/cusp selector
  needs at most 2*c_E*N_R + 2W-2g points, rather than 2*c_E^2+2W; here c_E is
  the component count of E_U and g=gcd(x,y).  In particular N_R<=B.
* Connected returns do not follow from primitivity, divisor completeness,
  exact signed complement, and full parity-twisted support.  The family

      B_n=506+360360*n,
      U_n=(1,2,3,4,7,B_n-6,B_n-3,B_n-2,B_n-1,B_n)

  has exactly 3+1440*n return components.

All predicates use integers or Fraction.  A separate interval-intersection
implementation cross-checks the cell classifier on deterministic random rows.

Tournament Analysis is telemetry only.  For U_0, use the three return
components as vertices.  Centre order is the pair observable; component width
is the switch/gauge, and signed cell label is the tie path.  Both tournaments
are transitive even though the two satellites have reciprocal endpoint-owner
handoffs 500->7 and 7->500.  Thus a bare tournament destroys precisely the
signed owner incidence needed by the return classifier.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import reduce
from hashlib import sha256
from itertools import combinations
from math import ceil, floor, gcd
from random import Random
from typing import Sequence


GAMMA = F(2, 143)
RAW_COMPLEMENT = frozenset(range(1, 13)) - {5, 8}
TWISTED_CLASSES = frozenset(range(1, 7))

GLOBAL_ROW = (2, 4, 6, 7, 9, 10, 11, 12, 14, 16)
ALL_GRID_ROW = (45, 48, 50, 54, 55, 62, 85, 95, 105, 116)
EVEN_GRID_ROW = (6, 9, 20, 24, 30, 36, 42, 54, 66, 90)
SATELLITE_ROW = (1, 2, 3, 4, 7, 500, 503, 504, 505, 506)


def fmt(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def norm(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def diamond_value(x: int, y: int, t: F) -> F:
    a = (x + y) // 2
    b = abs(x - y) // 2
    return norm(a * t) + norm(b * t)


def folded_class(value: int) -> int:
    residue = value % 13
    return min(residue, (-residue) % 13)


def parity_twisted_class(value: int) -> int:
    return folded_class(value if value % 2 else value // 2)


def balanced_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    if 2 * residue > modulus:
        residue -= modulus
    return residue


@dataclass(frozen=True)
class ReturnCell:
    label: int
    scaled_left: F
    scaled_right: F
    left: F
    right: F
    left_owners: tuple[int, ...]
    right_owners: tuple[int, ...]


def tooth_constraint(speed: int, label_mod_b: int, B: int) -> tuple[F, F]:
    residue = balanced_residue(speed * label_mod_b, B)
    return ((-B * GAMMA - residue) / speed, (B * GAMMA - residue) / speed)


def return_cells(speeds: tuple[int, ...]) -> tuple[ReturnCell, ...]:
    """Exact max-speed-tooth decomposition of ``closure(R_U)``.

    Cells with zero interior are omitted: R_U is open, so an isolated point in
    the intersection of the closed inequalities is not in its closure.
    """

    speeds = tuple(sorted(speeds))
    B = speeds[-1]
    cells: list[ReturnCell] = []
    for k in range(B):
        left, right = -GAMMA, GAMMA
        constraints = []
        for speed in speeds:
            lo, hi = tooth_constraint(speed, k, B)
            constraints.append((speed, lo, hi))
            left = max(left, lo)
            right = min(right, hi)
            if left >= right:
                break
        if left >= right:
            continue

        signed_label = k if 2 * k <= B else k - B
        interval_left = F(signed_label, B) + left / B
        interval_right = F(signed_label, B) + right / B
        left_owners = tuple(speed for speed, lo, _hi in constraints if lo == left)
        right_owners = tuple(speed for speed, _lo, hi in constraints if hi == right)
        cells.append(
            ReturnCell(
                signed_label,
                left,
                right,
                interval_left,
                interval_right,
                left_owners,
                right_owners,
            )
        )
    return tuple(sorted(cells, key=lambda cell: cell.left))


def intersect_open_intervals(
    left: Sequence[tuple[F, F]], right: Sequence[tuple[F, F]]
) -> list[tuple[F, F]]:
    answer: list[tuple[F, F]] = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            answer.append((lo, hi))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return answer


def direct_return_components(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    """Independent open-interval intersection, returned with closed endpoints."""

    current: list[tuple[F, F]] = [(F(-1, 2), F(1, 2))]
    for speed in sorted(speeds):
        allowed = []
        for integer in range(-speed, speed + 1):
            lo = (F(integer) - GAMMA) / speed
            hi = (F(integer) + GAMMA) / speed
            if hi > F(-1, 2) and lo < F(1, 2):
                allowed.append((max(lo, F(-1, 2)), min(hi, F(1, 2))))
        current = intersect_open_intervals(current, allowed)
    return tuple(current)


def cell_intervals(cells: Sequence[ReturnCell]) -> tuple[tuple[F, F], ...]:
    return tuple((cell.left, cell.right) for cell in cells)


def gcd_label_bound(speeds: tuple[int, ...]) -> int:
    """Necessary-label sieve from |[u*k]_B| <= 4B/143."""

    B = max(speeds)
    bounds = []
    for speed in speeds:
        g = gcd(speed, B)
        bounds.append(g * (2 * floor(F(4 * B, 143 * g)) + 1))
    return min(bounds)


def random_crosscheck() -> tuple[int, Counter[int]]:
    rng = Random(803807)
    rows = 0
    count_histogram: Counter[int] = Counter()
    for B in range(10, 81):
        for _ in range(3):
            sample = tuple(sorted((*rng.sample(range(1, B), 9), B)))
            if reduce(gcd, sample) != 1:
                continue
            classified = return_cells(sample)
            direct = direct_return_components(sample)
            assert cell_intervals(classified) == direct
            assert len(classified) <= B
            assert len(classified) <= gcd_label_bound(sample)
            rows += 1
            count_histogram[len(classified)] += 1
    return rows, count_histogram


def family_B(index: int) -> int:
    return 506 + 360_360 * index


def family_row(index: int) -> tuple[int, ...]:
    B = family_B(index)
    return (1, 2, 3, 4, 7, B - 6, B - 3, B - 2, B - 1, B)


def family_satellite_pairs(B: int) -> int:
    # Positive cell k exists exactly when 7(k-gamma) < gamma(B-6), i.e.
    # k < 2(B+1)/1001.  The strict inequality removes a zero-width cell.
    return ceil(F(2 * (B + 1), 1001)) - 1


def family_positive_interval(B: int, k: int) -> tuple[F, F]:
    assert 1 <= k <= family_satellite_pairs(B)
    return (
        (F(k) - GAMMA) / (B - 6),
        min((F(k) + GAMMA) / B, GAMMA / 7),
    )


def audit_branch_constraints(speeds: tuple[int, ...]) -> None:
    assert len(speeds) == len(set(speeds)) == 10
    assert reduce(gcd, speeds) == 1
    assert all(speed % 13 for speed in speeds)
    assert {speed % 13 for speed in speeds} == RAW_COMPLEMENT
    assert {parity_twisted_class(speed) for speed in speeds} == TWISTED_CLASSES
    assert all(any(speed % modulus == 0 for speed in speeds) for modulus in range(2, 13))


def audit_scalar_taxes(speeds: tuple[int, ...]) -> tuple[F, F]:
    """Certify the THM-772/797 scalar gates for ``(x,y)=(13,5)``.

    The consecutive-speed pigeonhole lemma gives
    ``M({1,2,3,4})=1/5``.  At ``t=1/5`` the extra speed 7 is at distance
    ``2/5``, so ``M({1,2,3,4,7})=1/5`` exactly.  Monotonicity under adding
    speeds therefore gives the metric upper bound used below.
    """

    B = max(speeds)
    x, y = 13, 5
    metric_upper = F(1, 5)
    assert min(norm(F(speed, 5)) for speed in (1, 2, 3, 4)) == metric_upper
    assert norm(F(7, 5)) >= metric_upper
    assert x <= 2 * B - 1 and y <= B - 1
    assert x * B + 2 * x * y <= 2 * B * (x + y)
    rho_upper = (metric_upper - F(1, 13)) / B
    scalar_left = F(1, x * y) + 2 * rho_upper
    scalar_right = F(2, x * x) + F(2, x * y)
    assert scalar_left <= scalar_right
    return scalar_left, scalar_right


def family_audit() -> tuple[tuple[int, int, int], ...]:
    table = []
    for index in range(6):
        B = family_B(index)
        speeds = family_row(index)
        audit_branch_constraints(speeds)
        audit_scalar_taxes(speeds)
        pairs = family_satellite_pairs(B)
        assert pairs == 1 + 720 * index
        components = 1 + 2 * pairs
        assert components == 3 + 1440 * index

        # Exact endpoint and midpoint tests for the first and last positive
        # cells.  The symbolic proof in the companion note excludes all other
        # k by the speed-7 inequality.
        for k in sorted({1, pairs}):
            lo, hi = family_positive_interval(B, k)
            assert lo < hi
            midpoint = (lo + hi) / 2
            assert all(norm(speed * midpoint) < GAMMA for speed in speeds)
            assert max(norm(speed * lo) for speed in speeds) == GAMMA
            assert max(norm(speed * hi) for speed in speeds) == GAMMA
        table.append((B, pairs, components))
    return tuple(table)


def benchmark_audit() -> dict[str, object]:
    for row in (GLOBAL_ROW, ALL_GRID_ROW, EVEN_GRID_ROW):
        assert len(return_cells(row)) == 1

    audit_branch_constraints(SATELLITE_ROW)
    scalar_left, scalar_right = audit_scalar_taxes(SATELLITE_ROW)
    cells = return_cells(SATELLITE_ROW)
    direct = direct_return_components(SATELLITE_ROW)
    assert cell_intervals(cells) == direct
    assert direct == (
        (F(-2, 1001), F(-141, 71500)),
        (F(-1, 36179), F(1, 36179)),
        (F(141, 71500), F(2, 1001)),
    )
    assert tuple(
        (cell.label, cell.left_owners, cell.right_owners) for cell in cells
    ) == (
        (-1, (7,), (500,)),
        (0, (506,), (506,)),
        (1, (500,), (7,)),
    )

    # The family is a return-geometry method limit, not a tight packet.  This
    # exact boundary witness uses only the central return component.  Its
    # positive Q-margin also persists after moving d into the open return set.
    deep_t = F(479, 616)
    central_d = F(1, 36179)
    selector_t = deep_t + central_d
    selector_q = diamond_value(13, 5, selector_t)
    selector_margin = F(11, 13) - selector_q
    assert min(norm(speed * deep_t) for speed in SATELLITE_ROW) == F(1, 11)
    assert central_d == cells[1].right
    assert max(norm(speed * central_d) for speed in SATELLITE_ROW) == GAMMA
    assert selector_t == F(1575487, 2026024)
    assert selector_q == F(226661, 2026024)
    assert selector_margin == F(1487667, 2026024) > 0

    B = max(SATELLITE_ROW)
    L = sum(SATELLITE_ROW)
    N = len(cells)
    W = 13
    g = 1
    old_bound = 2 * L * L + 2 * W
    adaptive_bound = 2 * L * N + 2 * W - 2 * g
    scale_bound = 20 * B * B + 22 * B - 2 * g
    assert old_bound == 12_852_476
    assert adaptive_bound == 15_234
    assert scale_bound == 5_131_850

    centre_order = tuple(cell.label for cell in cells)
    width_order = tuple(
        cell.label
        for cell in sorted(cells, key=lambda cell: (cell.right - cell.left, cell.label))
    )
    assert centre_order == (-1, 0, 1)
    assert width_order == (-1, 1, 0)
    centre_edges = total_order_edges(centre_order)
    width_edges = total_order_edges(width_order)
    flips = len(centre_edges ^ width_edges) // 2
    assert flips == 1

    return {
        "cells": cells,
        "gcd_bound": gcd_label_bound(SATELLITE_ROW),
        "old_bound": old_bound,
        "adaptive_bound": adaptive_bound,
        "scale_bound": scale_bound,
        "centre_order": centre_order,
        "width_order": width_order,
        "edge_flips": flips,
        "scalar_left": scalar_left,
        "scalar_right": scalar_right,
        "deep_t": deep_t,
        "central_d": central_d,
        "selector_t": selector_t,
        "selector_q": selector_q,
        "selector_margin": selector_margin,
    }


def total_order_edges(order: tuple[int, ...]) -> set[tuple[int, int]]:
    rank = {vertex: index for index, vertex in enumerate(order)}
    return {
        (left, right) if rank[left] < rank[right] else (right, left)
        for left, right in combinations(order, 2)
    }


def canonicalize(value: object) -> str:
    if isinstance(value, F):
        return fmt(value)
    if isinstance(value, ReturnCell):
        return (
            f"cell({value.label},{fmt(value.scaled_left)},{fmt(value.scaled_right)},"
            f"{fmt(value.left)},{fmt(value.right)},{value.left_owners},{value.right_owners})"
        )
    if isinstance(value, dict):
        return "{" + ",".join(
            f"{key}:{canonicalize(value[key])}" for key in sorted(value)
        ) + "}"
    if isinstance(value, (tuple, list)):
        return "(" + ",".join(canonicalize(item) for item in value) + ")"
    return repr(value)


def main() -> None:
    random_rows, random_histogram = random_crosscheck()
    family_table = family_audit()
    benchmark = benchmark_audit()

    canonical = canonicalize((random_rows, random_histogram, family_table, benchmark))
    digest = sha256(canonical.encode()).hexdigest()

    print("LRC13 RETURN SATELLITES — EXACT MAX-SPEED CELL CLASSIFIER")
    print("arithmetic=integer+fractions.Fraction floating_point=none")
    print()
    print("CELL_THEOREM")
    print("coordinate=d=(k+e)/B, -2/143<e<2/143")
    print("cell_interval=intersection_u [(-2B/143-[uk]_B)/u,(2B/143-[uk]_B)/u]")
    print("component_count=N_R<=B")
    print("gcd_sieve=N_R<=min_u gcd(u,B)*(2*floor(4B/(143*gcd(u,B)))+1)")
    print(f"independent_random_crosschecks={random_rows} component_histogram={dict(sorted(random_histogram.items()))}")
    print()
    print("ADAPTIVE_SELECTOR")
    print("bound=2*c_E*N_R+2*max(x,y)-2*gcd(x,y)")
    print("coarse_bound=2*L*B+2W-2g<=20B^2+22B-2g")
    print(
        f"satellite_row_old_L2_bound={benchmark['old_bound']} "
        f"adaptive_LN_bound={benchmark['adaptive_bound']} "
        f"coarse_scale_bound={benchmark['scale_bound']}"
    )
    print()
    print("SIGNED_COMPLEMENT_METHOD_LIMIT")
    print(f"U_0={SATELLITE_ROW}")
    print(
        "properties=primitive,divisor_complete_2..12,exact_signed_complement_5_8,"
        "full_parity_twisted_support"
    )
    print(
        f"scalar_taxes_x13_y5=lhs_at_M_upper_1/5_{fmt(benchmark['scalar_left'])}"
        f"<=rhs_{fmt(benchmark['scalar_right'])}"
    )
    for cell in benchmark["cells"]:
        print(
            f"cell_label={cell.label} interval=[{fmt(cell.left)},{fmt(cell.right)}] "
            f"scaled=[{fmt(cell.scaled_left)},{fmt(cell.scaled_right)}] "
            f"owners=({cell.left_owners},{cell.right_owners})"
        )
    print(f"gcd_label_bound={benchmark['gcd_bound']} actual_components=3")
    print(
        f"central_selector_witness=s_{fmt(benchmark['deep_t'])} "
        f"d_{fmt(benchmark['central_d'])} t_{fmt(benchmark['selector_t'])} "
        f"Q_{fmt(benchmark['selector_q'])} margin_{fmt(benchmark['selector_margin'])}"
    )
    print(f"family_table_(B,satellite_pairs,components)={family_table}")
    print("exact_family_law=B=506+360360n => N_R=3+1440n")
    print()
    print("TOURNAMENT_TELEMETRY")
    print("vertices=return_components pair_observable=centre_order switch=width_order tie=signed_cell_label")
    print(
        f"centre_order={benchmark['centre_order']} width_order={benchmark['width_order']} "
        f"edge_flips={benchmark['edge_flips']}"
    )
    print("both_fingerprints=scores_(0,1,2),cycles_0,SCCs_(1,1,1),Hamiltonian_paths_1")
    print("destroyed_by_bare_tournament=reciprocal_satellite_owner_handoffs_500_to_7_and_7_to_500")
    print()
    print(f"sha256={digest}")
    print("PASS: exact satellite cells compress the selector; signed-complement arithmetic permits linearly many satellites.")


if __name__ == "__main__":
    main()
