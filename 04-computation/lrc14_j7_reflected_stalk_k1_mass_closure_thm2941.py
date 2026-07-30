#!/usr/bin/env python3
r"""Exact closure of THM-2941's canonical reflected-stalk ``k=1`` family.

Let ``E`` be any six-element subset of ``{1,...,14}``, put

    L = 14*lcm(E),               Z_E = {L-e : e in E},

and let ``P_(E,Z_E)`` be the lossless projected residual of THM-2941(25g).
If one further aligned multiplier ``a`` covered the original residual, then

    P_(E,Z_E) subset D_a,        hence mu(P_(E,Z_E)) <= 1/7.          (1)

On a body-safe ``1/L`` cell ``t=(j+u)/L``, the reflected label ``L-e``
is dangerous exactly when

    || ((L-e)u-ej)/L || < 1/14.

Write ``U_j`` for the union of these six open subsets of ``0<=u<=1``.
The finite-clause identity THM-2941(25i) gives, up to finitely many seams,

    T \ P_(E,Z_E) = intersection_{j in J_E} U_j.

Consequently any one safe cell with ``mu(U_j)<6/7`` proves
``mu(P_(E,Z_E))>1/7``, contradicting (1) for *every* integer ``a``.

This script gives such a cell for all C(14,6)=3003 roots.  The selector is
particularly small: use the leftmost cell of the first positive safe
component, except for six explicitly frozen roots, where the leftmost cell
of the second component works.  Two algebraically independent descriptions
of every reflected interval are compared exactly, and the union mass is
also rebuilt by an endpoint-slab sweep independent of interval merging.

All arithmetic is integer or Fraction exact.  RuntimeError checks remain
active under ``python -O``.  The positive mass gap makes all open/closed seam
choices immaterial to the conclusion; no finite search over the final label
``a`` and no monochromatic-component/7-bin sidecar are used.
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as Q
from itertools import combinations
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_stalk_k1_mass_closure_thm2941.out"
)

THRESHOLD = Q(6, 7)
EXPECTED_FIRST_EXCEPTIONS = (
    (1, 3, 5, 7, 9, 11),
    (1, 3, 5, 7, 9, 12),
    (1, 3, 5, 7, 10, 12),
    (1, 3, 5, 8, 10, 12),
    (1, 3, 6, 8, 10, 12),
    (1, 4, 6, 8, 10, 12),
)
EXPECTED_MAXIMUM = (
    Q(10028643748, 12527514945),
    (1, 3, 4, 6, 8, 12),
    336,
    24,
    0,
)
EXPECTED_SEMANTIC_SHA256 = (
    "8fc2b276e8899622d54964325adfa4c5f23259d0d34367f812b12dc4d72c2d6c"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def safe_cell_ranges(E: tuple[int, ...]) -> tuple[int, tuple[tuple[int, int], ...]]:
    """Exact positive-length body-safe ranges on the integer ``L`` ruler."""
    L = lcm(*(14 * e for e in E))
    require(L == 14 * lcm(*E), ("body period changed", E, L))
    danger = []
    for e in E:
        radius = L // (14 * e)
        period = L // e
        require(14 * e * radius == L, ("nonintegral body tooth", E, e, L))
        for k in range(e + 1):
            center = k * period
            danger.append(
                (max(0, center - radius), min(L, center + radius))
            )
    danger.sort()
    merged = []
    for left, right in danger:
        if merged and left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    safe = []
    cursor = 0
    for left, right in merged:
        if cursor < left:
            safe.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < L:
        safe.append((cursor, L))
    require(safe, ("empty six-body carrier", E, L))
    require(
        all(0 <= left < right <= L for left, right in safe),
        ("bad safe range", E, L, safe),
    )
    return L, tuple(safe)


def clip(left: Q, right: Q) -> tuple[Q, Q] | None:
    left = max(Q(0), left)
    right = min(Q(1), right)
    return (left, right) if left < right else None


def direct_multiplier_arcs(L: int, a: int, j: int) -> tuple[tuple[Q, Q], ...]:
    """Directly pull back ``||a(j+u)/L||<1/14`` to the cell coordinate."""
    require(0 < a < L, ("multiplier outside direct route", L, a, j))
    lo = (a * j) // L - 1
    hi = (a * (j + 1)) // L + 1
    arcs = []
    for n in range(lo, hi + 1):
        left = Q(L * (14 * n - 1), 14 * a) - j
        right = Q(L * (14 * n + 1), 14 * a) - j
        piece = clip(left, right)
        if piece is not None:
            arcs.append(piece)
    return tuple(sorted(set(arcs)))


def reflected_identity_arcs(L: int, e: int, j: int) -> tuple[tuple[Q, Q], ...]:
    """Independently use ``(L-e)(j+u)=-ej+(L-e)u (mod L)``."""
    z = L - e
    r = (e * j) % L
    lo = (-r) // L - 1
    hi = (z - r) // L + 1
    radius = Q(L, 14 * z)
    arcs = []
    for k in range(lo, hi + 1):
        center = Q(r + k * L, z)
        piece = clip(center - radius, center + radius)
        if piece is not None:
            arcs.append(piece)
    return tuple(sorted(set(arcs)))


def merge_intervals(arcs: tuple[tuple[Q, Q], ...]) -> tuple[tuple[Q, Q], ...]:
    merged = []
    for left, right in sorted(arcs):
        if merged and left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    return tuple(merged)


def merged_mass(arcs: tuple[tuple[Q, Q], ...]) -> Q:
    return sum((right - left for left, right in merge_intervals(arcs)), Q(0))


def slab_sweep_mass(arcs: tuple[tuple[Q, Q], ...]) -> Q:
    """Independent union measure from coverage of consecutive endpoint slabs."""
    endpoints = sorted({Q(0), Q(1), *(x for arc in arcs for x in arc)})
    total = Q(0)
    for left, right in zip(endpoints, endpoints[1:]):
        middle = (left + right) / 2
        if any(a < middle < b for a, b in arcs):
            total += right - left
    return total


def body_cell_is_safe_directly(L: int, E: tuple[int, ...], j: int) -> bool:
    """A second check that the selected open cell misses every body danger set."""
    return all(not direct_multiplier_arcs(L, e, j) for e in E)


def reflected_cell_union(
    L: int, E: tuple[int, ...], j: int
) -> tuple[Q, tuple[tuple[Q, Q], ...], int]:
    direct = []
    reflected = []
    clause_count = 0
    for e in E:
        z = L - e
        route_a = direct_multiplier_arcs(L, z, j)
        route_b = reflected_identity_arcs(L, e, j)
        require(route_a == route_b, ("reflected clause mismatch", E, L, e, j))
        direct.extend(route_a)
        reflected.extend(route_b)
        clause_count += 1
    direct_tuple = tuple(direct)
    reflected_tuple = tuple(reflected)
    require(direct_tuple == reflected_tuple, ("reflected aggregate mismatch", E, j))
    mass_a = merged_mass(direct_tuple)
    mass_b = slab_sweep_mass(reflected_tuple)
    require(mass_a == mass_b, ("union mass mismatch", E, L, j, mass_a, mass_b))
    return mass_a, merge_intervals(direct_tuple), clause_count


def qtext(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()

    records = []
    first_exceptions = []
    clause_checks = 0
    safe_checks = 0
    witness_digest = hashlib.sha256()

    for E in combinations(range(1, 15), 6):
        L, safe = safe_cell_ranges(E)
        first_j = safe[0][0]
        require(body_cell_is_safe_directly(L, E, first_j), ("unsafe first cell", E))
        safe_checks += len(E)
        first_mass, first_union, checks = reflected_cell_union(L, E, first_j)
        clause_checks += checks
        component_index = 0
        witness_j = first_j
        witness_mass = first_mass
        witness_union = first_union
        if first_mass >= THRESHOLD:
            first_exceptions.append(E)
            require(len(safe) >= 2, ("first failure has no second component", E))
            component_index = 1
            witness_j = safe[1][0]
            require(
                body_cell_is_safe_directly(L, E, witness_j),
                ("unsafe second cell", E),
            )
            safe_checks += len(E)
            witness_mass, witness_union, checks = reflected_cell_union(
                L, E, witness_j
            )
            clause_checks += checks
        require(
            witness_mass < THRESHOLD,
            ("reflected-stalk mass survivor", E, L, witness_j, witness_mass),
        )
        record = (
            E,
            L,
            len(safe),
            component_index,
            witness_j,
            witness_mass,
            witness_union,
            first_mass,
        )
        records.append(record)
        witness_digest.update((repr(record) + "\n").encode())

    require(len(records) == 3003, ("body universe changed", len(records)))
    require(
        tuple(first_exceptions) == EXPECTED_FIRST_EXCEPTIONS,
        ("first-component exception set changed", first_exceptions),
    )
    maximum_record = max(records, key=lambda row: (row[5], row[0]))
    maximum = (
        maximum_record[5],
        maximum_record[0],
        maximum_record[1],
        maximum_record[4],
        maximum_record[3],
    )
    require(maximum == EXPECTED_MAXIMUM, ("maximum witness mass changed", maximum))
    uniform_gap = THRESHOLD - maximum[0]
    require(uniform_gap == Q(4964583434, 87692604615), ("gap changed", uniform_gap))

    semantic_payload = (
        tuple(records),
        tuple(first_exceptions),
        clause_checks,
        safe_checks,
        maximum,
        uniform_gap,
        witness_digest.hexdigest(),
    )
    semantic_sha256 = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest changed", semantic_sha256),
        )

    maximum_union = maximum_record[6]
    exception_rows = tuple(
        (row[0], row[1], row[4], row[5]) for row in records if row[3] == 1
    )
    lines = [
        "LRC14 THM-2941 canonical reflected-stalk k=1 mass closure",
        (
            "scope=all six-body roots E subset {1,...,14};"
            "Z_E={L-e:e in E};one remaining aligned label a;"
            "not the whole j=7 sector"
        ),
        (
            "clause_on_safe_cell=(L-e)(j+u)=-ej+(L-e)u mod L;"
            "U_j=union_e {u:||((L-e)u-ej)/L||<1/14}"
        ),
        (
            "projection_law=T\\P_(E,Z_E)=intersection_{j in J_E} U_j a.e.;"
            "therefore mu(P)>=1-mu(U_j)"
        ),
        (
            "final_comb_obstruction=P subset D_a implies mu(P)<=1/7 for every positive integer a;"
            "a witness mu(U_j)<6/7 rules out all a without a label horizon"
        ),
        (
            "selector=leftmost cell of first positive safe component;"
            "on failure use leftmost cell of second positive safe component"
        ),
        f"body_roots={len(records)};first_component_witnesses={len(records)-len(first_exceptions)};second_component_witnesses={len(first_exceptions)}",
        f"first_component_exceptions={tuple(first_exceptions)}",
        f"exception_witness_rows={exception_rows}",
        f"reflected_clause_route_comparisons={clause_checks};direct_body_safe_checks={safe_checks}",
        "union_mass_routes=merged_exact_intervals;independent_endpoint_slab_sweep;all_equal",
        (
            f"maximum_selected_union_mass={qtext(maximum[0])};"
            f"body={maximum[1]};L={maximum[2]};j={maximum[3]};"
            f"safe_component_index={maximum[4]}"
        ),
        f"maximum_selected_union_intervals={maximum_union}",
        f"uniform_gap_6/7_minus_union={qtext(uniform_gap)}",
        (
            f"uniform_residual_floor=1/7+{qtext(uniform_gap)};"
            "strictly_greater_than_final_comb_mass=1/7"
        ),
        (
            "endpoint_status=open danger intervals retained algebraically;"
            "mass uses closed representatives only at finitely many seams;"
            "strict uniform gap makes seam convention irrelevant"
        ),
        (
            "component_sidecar=not invoked;the stronger mass obstruction closes the family "
            "before monochromatic-component or exact 7-bin recursion"
        ),
        f"witness_digest_sha256={witness_digest.hexdigest()}",
        f"semantic_sha256={semantic_sha256}",
        "conclusion=canonical reflected-stalk k=1 family is empty for all 3003 body roots",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
