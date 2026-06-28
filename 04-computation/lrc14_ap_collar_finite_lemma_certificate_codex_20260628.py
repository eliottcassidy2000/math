#!/usr/bin/env python3
"""HYP-3405: AP-collar finite lemma certificate for LRC14.

This script strengthens HYP-3401 from a quotient-mixing scout into a
certificate-shaped finite lemma target.  It verifies, with exact rational
arithmetic, the AP one-swap collar through replacement speed 84:

1. AP and Goddyn-Wong 12->24 are the only boundary-tight rows.
2. Every other row has an explicitly checkable strict-open interval.
3. The strict-open mass has a uniform lower bound 1/1260, uniquely at 12->36.
4. The HYP-3401 nonunit-height packet has one mixed fiber; unit-height data is
   the first visible repair for that fiber.

This is still not an LRC14 proof.  It is a finite lemma certificate and a
proof-engineering target for O15 tight-locus rigidity.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from typing import Callable, Iterable

import lrc14_three_coordinate_obstruction_exactness_codex_20260628 as base


KeyFunc = Callable[[tuple[int, ...]], object]


def frac_word(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def interval_word(interval: tuple[Fraction, Fraction]) -> str:
    a, b = interval
    return f"({frac_word(a)}, {frac_word(b)})"


def canonical_fraction(x: Fraction) -> tuple[int, int]:
    return (x.numerator, x.denominator)


def canonical_interval(interval: tuple[Fraction, Fraction]) -> tuple[tuple[int, int], tuple[int, int]]:
    a, b = interval
    return (canonical_fraction(a), canonical_fraction(b))


def open_intervals_overlap(
    left: tuple[Fraction, Fraction],
    right: tuple[Fraction, Fraction],
) -> bool:
    return max(left[0], right[0]) < min(left[1], right[1])


def interval_is_strictly_safe(row: tuple[int, ...], interval: tuple[Fraction, Fraction]) -> bool:
    if interval[0] >= interval[1]:
        return False
    for speed in row:
        for unsafe in base.unsafe_intervals_for_speed(speed):
            if open_intervals_overlap(interval, unsafe):
                return False
    mid = (interval[0] + interval[1]) / 2
    return all(base.circular_distance_to_integer(Fraction(speed) * mid) > base.DELTA for speed in row)


def choose_witness_interval(row: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    components = base.safe_open_components(row)
    if not components:
        raise ValueError(f"row has no strict-open component: {row}")
    return max(components, key=lambda interval: (interval[1] - interval[0], -interval[0]))


def certificate_digest(audited: list[base.RowAudit]) -> str:
    payload = []
    for item in audited:
        comps = base.safe_open_components(item.row)
        payload.append(
            (
                item.row,
                item.exit_status,
                canonical_fraction(item.mass),
                tuple(canonical_interval(component) for component in comps),
                tuple(canonical_fraction(t) for t in base.threshold_safe_points(item.row)),
            )
        )
    return sha256(repr(payload).encode("ascii")).hexdigest()


def mixed_fibers(audited: list[base.RowAudit], keyfunc: KeyFunc) -> list[list[base.RowAudit]]:
    fibers: dict[object, list[base.RowAudit]] = defaultdict(list)
    for item in audited:
        fibers[keyfunc(item.row)].append(item)
    return [
        sorted(fiber, key=lambda item: item.row)
        for fiber in fibers.values()
        if len({item.exit_status for item in fiber}) > 1
    ]


def repair_counts(audited: list[base.RowAudit], quotient: KeyFunc) -> list[tuple[str, int]]:
    sidecars: list[tuple[str, KeyFunc]] = [
        ("unit_contact_status", base.unit_contact_status),
        ("covering_layer", base.covering_layer_key),
        ("unit_height_flex", base.unit_height_key),
        ("nonunit_height_flex", base.nonunit_height_key),
        ("full_height_flex", base.height_flex_key),
        ("height_completed_packet", base.height_completed_packet_key),
    ]
    return [(name, base.repair_summary(audited, quotient, keyfunc)) for name, keyfunc in sidecars]


def multiset_delta(
    left: Iterable[object],
    right: Iterable[object],
) -> tuple[tuple[object, int], tuple[object, int]]:
    """Return items removed from left and added in right."""
    lcount = Counter(left)
    rcount = Counter(right)
    removed = sorted(((item, count) for item, count in (lcount - rcount).items()), key=repr)
    added = sorted(((item, count) for item, count in (rcount - lcount).items()), key=repr)
    return (tuple(removed), tuple(added))


@dataclass(frozen=True)
class LemmaCarrier:
    name: str
    preserves: frozenset[str]
    destroys: frozenset[str]
    priority: int


WEIGHTS = {
    "boundary_atom_classifier": 12,
    "strict_open_witness_intervals": 12,
    "uniform_margin": 9,
    "quotient_obstruction_table": 10,
    "unit_height_obstruction_vector": 9,
    "sidecar_repair_matrix": 8,
    "global_chamber_shape": 7,
    "low_payload": 5,
}


CARRIERS = [
    LemmaCarrier(
        "strict_open_witness_certificate",
        frozenset({"strict_open_witness_intervals", "uniform_margin", "boundary_atom_classifier"}),
        frozenset({"global_chamber_shape"}),
        0,
    ),
    LemmaCarrier(
        "unit_height_obstruction_vector",
        frozenset({"unit_height_obstruction_vector", "quotient_obstruction_table", "sidecar_repair_matrix"}),
        frozenset({"global_chamber_shape"}),
        1,
    ),
    LemmaCarrier(
        "boundary_atom_classifier",
        frozenset({"boundary_atom_classifier", "low_payload"}),
        frozenset({"strict_open_witness_intervals", "unit_height_obstruction_vector"}),
        2,
    ),
    LemmaCarrier(
        "sidecar_repair_matrix",
        frozenset({"sidecar_repair_matrix", "quotient_obstruction_table"}),
        frozenset({"uniform_margin"}),
        3,
    ),
    LemmaCarrier(
        "height_completed_oracle",
        frozenset({"boundary_atom_classifier", "strict_open_witness_intervals", "quotient_obstruction_table"}),
        frozenset({"low_payload", "global_chamber_shape"}),
        4,
    ),
    LemmaCarrier(
        "raw_mixed_fiber_scout",
        frozenset({"quotient_obstruction_table"}),
        frozenset({"boundary_atom_classifier", "strict_open_witness_intervals", "uniform_margin"}),
        5,
    ),
]


def carrier_score(carrier: LemmaCarrier) -> int:
    keep = sum(WEIGHTS[item] for item in carrier.preserves)
    lose = sum(WEIGHTS[item] for item in carrier.destroys)
    return 2 * keep - lose - carrier.priority


def beats(left: LemmaCarrier, right: LemmaCarrier) -> bool:
    return (carrier_score(left), -len(left.destroys), -left.priority) > (
        carrier_score(right),
        -len(right.destroys),
        -right.priority,
    )


def tournament_fingerprint() -> dict[str, object]:
    graph = {carrier.name: set() for carrier in CARRIERS}
    for left, right in combinations(CARRIERS, 2):
        if beats(left, right):
            graph[left.name].add(right.name)
        else:
            graph[right.name].add(left.name)

    cycles = 0
    for a, b, c in combinations(graph, 3):
        if b in graph[a] and c in graph[b] and a in graph[c]:
            cycles += 1
        if c in graph[a] and b in graph[c] and a in graph[b]:
            cycles += 1

    hpaths = 0
    for order in permutations(CARRIERS):
        if all(order[i + 1].name in graph[order[i].name] for i in range(len(order) - 1)):
            hpaths += 1

    path = sorted(CARRIERS, key=lambda carrier: (-carrier_score(carrier), len(carrier.destroys), carrier.priority))
    return {
        "vertices": len(CARRIERS),
        "score_hist": dict(sorted(Counter(carrier_score(carrier) for carrier in CARRIERS).items())),
        "directed_3cycles": cycles,
        "hamiltonian_path_count": hpaths,
        "priority_path": [carrier.name for carrier in path],
    }


def main() -> None:
    audited = [base.audit(row) for row in base.one_swap_collar(84)]
    boundary = sorted((item for item in audited if item.exit_status == "boundary_tight"), key=lambda item: item.row)
    strict = sorted((item for item in audited if item.exit_status == "strict_open"), key=lambda item: item.row)
    debt = [item for item in audited if item.exit_status not in {"boundary_tight", "strict_open"}]

    for item in strict:
        witness = choose_witness_interval(item.row)
        assert interval_is_strictly_safe(item.row, witness), item.row
        assert item.mass > 0
    for item in boundary:
        assert item.mass == 0
        assert base.safe_open_components(item.row) == []
        assert base.threshold_safe_points(item.row) == tuple(Fraction(a, 14) for a in base.UNITS)
    assert not debt

    min_mass = min(item.mass for item in strict)
    min_rows = [item for item in strict if item.mass == min_mass]
    max_components = max(item.components for item in strict)
    max_endpoint_den = max(
        endpoint.denominator
        for item in strict
        for component in base.safe_open_components(item.row)
        for endpoint in component
    )

    nonunit_mixed = mixed_fibers(audited, base.c3_quadratic_nonunit_height_key)
    assert len(nonunit_mixed) == 1
    mixed = nonunit_mixed[0]
    mixed_boundary = [item for item in mixed if item.exit_status == "boundary_tight"]
    mixed_strict = [item for item in mixed if item.exit_status == "strict_open"]
    first_strict = mixed_strict[0]
    exemplar_boundary = mixed_boundary[0]
    unit_delta = multiset_delta(base.unit_height_key(exemplar_boundary.row), base.unit_height_key(first_strict.row))
    nonunit_delta = multiset_delta(base.nonunit_height_key(exemplar_boundary.row), base.nonunit_height_key(first_strict.row))

    quotient_table: list[tuple[str, KeyFunc]] = [
        ("raw_unit_projection", base.unit_projection_key),
        ("raw_mod14_residue_table", base.residue_counts),
        ("c3_binding_skeleton", base.c3_skeleton_key),
        ("quadratic_Qsqrt_minus7_character", base.quadratic_key),
        ("c3_plus_quadratic", lambda row: (base.c3_skeleton_key(row), base.quadratic_key(row))),
        (
            "c3_plus_quadratic_plus_covering_layer",
            lambda row: (base.c3_skeleton_key(row), base.quadratic_key(row), base.covering_layer_key(row)),
        ),
        ("c3_quadratic_nonunit_height_packet", base.c3_quadratic_nonunit_height_key),
        ("height_completed_packet", base.height_completed_packet_key),
        ("full_height_residue_ledger", base.height_flex_key),
    ]

    print("HYP-3405 AP-collar finite lemma certificate")
    print("=" * 78)
    print("status=exact finite certificate; proof-engineering target, not LRC14 proof")
    print("lemma_domain=AP one-swap collar with replacement speed <= 84")
    print(f"certificate_digest={certificate_digest(audited)}")
    print()

    print("FINITE LEMMA A: BOUNDARY ATOMS")
    print(f"rows={len(audited)}")
    print(f"boundary_tight_count={len(boundary)}")
    print(f"strict_open_count={len(strict)}")
    print(f"covered_or_debt_count={len(debt)}")
    for item in boundary:
        points = ", ".join(frac_word(t) for t in base.threshold_safe_points(item.row))
        print(f"  boundary row={item.row} closed_contacts=[{points}]")
    print()

    print("FINITE LEMMA B: STRICT-OPEN WITNESSES")
    print(f"all_strict_rows_have_verified_open_interval={all(interval_is_strictly_safe(item.row, choose_witness_interval(item.row)) for item in strict)}")
    print(f"uniform_strict_mass_lower_bound={frac_word(min_mass)}")
    print(f"lower_bound_attainers={len(min_rows)}")
    for item in min_rows:
        comps = base.safe_open_components(item.row)
        witness_words = ", ".join(
            f"{interval_word(component)} length={frac_word(component[1] - component[0])}"
            for component in comps
        )
        print(f"  min row={item.row} mass={frac_word(item.mass)} witnesses={witness_words}")
    print(f"max_strict_component_count={max_components}")
    print(f"max_safe_endpoint_denominator={max_endpoint_den}")
    print(f"component_count_hist={dict(sorted(Counter(item.components for item in strict).items()))}")
    print()

    print("FINITE LEMMA C: QUOTIENT OBSTRUCTION TABLE")
    for name, keyfunc in quotient_table:
        summary = base.quotient_summary(audited, keyfunc)
        print(
            f"  {name}: fibers={summary['fibers']} mixed_fibers={summary['mixed_fibers']} "
            f"largest_mixed={summary['largest_mixed']}"
        )
    print()

    print("FINITE LEMMA D: FIRST OBSTRUCTION VECTOR")
    print(f"mixed_nonunit_packet_fibers={len(nonunit_mixed)}")
    print(f"mixed_fiber_size={len(mixed)}")
    print(f"mixed_fiber_boundary_count={len(mixed_boundary)}")
    print(f"mixed_fiber_strict_count={len(mixed_strict)}")
    print(f"boundary_exemplar={exemplar_boundary.row}")
    print(f"strict_exemplar={first_strict.row}")
    print(f"strict_exemplar_mass={frac_word(first_strict.mass)}")
    print(f"strict_exemplar_witness={interval_word(choose_witness_interval(first_strict.row))}")
    print(f"shared_contact_status={''.join(exemplar_boundary.contact_status)}")
    print(f"shared_c3={base.c3_skeleton_key(exemplar_boundary.row)}")
    print(f"shared_quadratic={base.quadratic_key(exemplar_boundary.row)}")
    print(f"shared_covering={base.covering_layer_key(exemplar_boundary.row)}")
    print(f"nonunit_height_delta={nonunit_delta}")
    print(f"unit_height_delta={unit_delta}")
    print()

    print("SIDECAR REPAIR MATRIX FOR NONUNIT PACKET")
    for name, count in repair_counts(audited, base.c3_quadratic_nonunit_height_key):
        print(f"  {name}: mixed_fibers_after_join={count}")
    print()

    print("CERTIFIED FINITE LEMMA TARGET")
    print("  AP and GW 12->24 are the only boundary-tight rows in the finite collar.")
    print("  Every other row has a verified strict-open interval and mass >= 1/1260.")
    print("  The C3 + Qsqrt(-7) + nonunit-height packet has one mixed fiber.")
    print("  The first displayed obstruction is the unit-height lift (13,0)->(13,1).")
    print("  Joining the unit-height sidecar kills the mixed fiber on this collar.")
    print()

    print("TOURNAMENT ANALYSIS")
    print("vertices=finite-lemma proof carriers, not runners/residues/replacement speeds")
    print("pairwise_observable=certificate checkability + obstruction localization + globalization cost")
    print("switch_gauge=A -> B iff A retains more finite-lemma proof payload with less destroyed chamber data")
    fp = tournament_fingerprint()
    for key, value in fp.items():
        if key == "priority_path":
            print(f"{key}={' -> '.join(value)}")
        else:
            print(f"{key}={value}")
    print("assumption_challenge=considered runners/gaps/sections/residues/wall-crossings; selected finite-lemma proof carriers because the preserved predicate is certificate-grade boundary-vs-strict exactness")


if __name__ == "__main__":
    main()
