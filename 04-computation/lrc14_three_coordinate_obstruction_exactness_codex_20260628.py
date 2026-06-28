#!/usr/bin/env python3
"""HYP-3401 scout: three-coordinate obstruction exactness for LRC14.

This is a finite exact scout, not an LRC14 proof.  It instantiates the
HYP-3301 first-obstruction sheaf idea on the concrete HYP-3311 packet:

    C3 binding skeleton
    Q(sqrt(-7)) / quadratic character sidecar
    height/flex ledger

The test bank is the AP one-swap collar with replacement speeds up to 84.  For
each quotient we ask whether its fibers mix boundary-tight rows with strict
open rows.  A mixed fiber is a concrete "first obstruction" to using that
quotient as proof currency without a sidecar.  The key surprise is that the
height/flex ledger must include unit-height lifts too, not only nonunit
covering-layer heights.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
from math import gcd
from typing import Callable, Iterable


DELTA = Fraction(1, 14)
MOD = 14
P = 7
AP = tuple(range(1, 14))
UNITS = (1, 3, 5, 9, 11, 13)
NONZERO_RESIDUES = tuple(range(1, 14))
COMPLEMENT_PAIRS = ((1, 13), (3, 11), (5, 9))
PAIR_INDEX = {pair: i for i, pair in enumerate(COMPLEMENT_PAIRS)}


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def circular_distance_to_integer(x: Fraction) -> Fraction:
    y = frac_part(x)
    return min(y, 1 - y)


def split_circle_interval(a: Fraction, b: Fraction) -> list[tuple[Fraction, Fraction]]:
    while a < 0:
        a += 1
        b += 1
    while a >= 1:
        a -= 1
        b -= 1
    if b <= 1:
        return [(a, b)]
    return [(a, Fraction(1)), (Fraction(0), b - 1)]


def unsafe_intervals_for_speed(v: int) -> tuple[tuple[Fraction, Fraction], ...]:
    radius = DELTA / v
    out: list[tuple[Fraction, Fraction]] = []
    for m in range(v):
        center = Fraction(m, v)
        out.extend(split_circle_interval(center - radius, center + radius))
    return tuple(sorted(out))


def endpoint_candidates_for_speed(v: int) -> tuple[Fraction, ...]:
    out: list[Fraction] = []
    for m in range(v):
        out.append(frac_part((Fraction(m) - DELTA) / v))
        out.append(frac_part((Fraction(m) + DELTA) / v))
    return tuple(sorted(set(out)))


def merge_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not intervals:
        return []
    intervals.sort()
    merged: list[list[Fraction]] = [[intervals[0][0], intervals[0][1]]]
    for a, b in intervals[1:]:
        if a <= merged[-1][1]:
            if b > merged[-1][1]:
                merged[-1][1] = b
        else:
            merged.append([a, b])
    return [(a, b) for a, b in merged]


def complement_intervals(covered: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not covered:
        return [(Fraction(0), Fraction(1))]
    out: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for a, b in covered:
        if cursor < a:
            out.append((cursor, a))
        if b > cursor:
            cursor = b
    if cursor < 1:
        out.append((cursor, Fraction(1)))
    return out


def safe_open_components(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    intervals: list[tuple[Fraction, Fraction]] = []
    for v in speeds:
        intervals.extend(unsafe_intervals_for_speed(v))
    return complement_intervals(merge_intervals(intervals))


def interval_measure(intervals: Iterable[tuple[Fraction, Fraction]]) -> Fraction:
    return sum((b - a for a, b in intervals), start=Fraction(0))


def threshold_safe_points(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    candidates = {Fraction(0)}
    for v in speeds:
        candidates.update(endpoint_candidates_for_speed(v))
    return tuple(
        t
        for t in sorted(candidates)
        if all(circular_distance_to_integer(Fraction(v) * t) >= DELTA for v in speeds)
    )


def replace_speed(row: tuple[int, ...], old: int, new: int) -> tuple[int, ...]:
    return tuple(sorted(new if x == old else x for x in row))


def one_swap_collar(max_new_speed: int = 84) -> list[tuple[int, ...]]:
    rows = {AP}
    for old in AP:
        base = tuple(x for x in AP if x != old)
        for new in range(1, max_new_speed + 1):
            if new in base:
                continue
            rows.add(tuple(sorted(base + (new,))))
    return sorted(rows)


def legendre7(a: int) -> int:
    r = a % P
    if r == 0:
        return 0
    return 1 if r in {1, 2, 4} else -1


def c3_slot_for_residue(r: int) -> int | None:
    r %= MOD
    if gcd(r, MOD) != 1:
        return None
    pair = tuple(sorted((r, (-r) % MOD)))
    return PAIR_INDEX[pair]


def unit_contact_status(row: tuple[int, ...]) -> tuple[str, ...]:
    out = []
    for a in UNITS:
        mind = min(circular_distance_to_integer(Fraction(s * a, MOD)) for s in row)
        if mind < DELTA:
            out.append("K")
        elif mind == DELTA:
            out.append("E")
        else:
            out.append("S")
    return tuple(out)


def residue_counts(row: tuple[int, ...]) -> tuple[int, ...]:
    counts = Counter(s % MOD for s in row)
    return tuple(counts[r] for r in range(MOD))


def unit_projection_key(row: tuple[int, ...]) -> tuple[int, ...]:
    counts = Counter(s % MOD for s in row)
    return tuple(counts[r] for r in UNITS)


def c3_skeleton_key(row: tuple[int, ...]) -> tuple[int, ...]:
    counts = Counter()
    for s in row:
        slot = c3_slot_for_residue(s)
        if slot is not None:
            counts[slot] += 1
    return tuple(counts[i] for i in range(3))


def quadratic_key(row: tuple[int, ...]) -> tuple[int, int, int]:
    counts = Counter(legendre7(s) for s in row)
    return (counts[-1], counts[0], counts[1])


def covering_layer_key(row: tuple[int, ...]) -> tuple[int, int]:
    counts = Counter("apex7" if s % 7 == 0 else "even2U" for s in row if gcd(s, MOD) != 1)
    return (counts["even2U"], counts["apex7"])


def height_flex_key(row: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted((s % MOD, s // MOD) for s in row))


def nonunit_height_key(row: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted((s % MOD, s // MOD) for s in row if gcd(s, MOD) != 1))


def unit_height_key(row: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted((s % MOD, s // MOD) for s in row if gcd(s, MOD) == 1))


def c3_quadratic_nonunit_height_key(row: tuple[int, ...]) -> tuple[object, ...]:
    return (
        c3_skeleton_key(row),
        quadratic_key(row),
        covering_layer_key(row),
        nonunit_height_key(row),
    )


def height_completed_packet_key(row: tuple[int, ...]) -> tuple[object, ...]:
    return (
        c3_skeleton_key(row),
        quadratic_key(row),
        covering_layer_key(row),
        height_flex_key(row),
    )


@dataclass(frozen=True)
class RowAudit:
    row: tuple[int, ...]
    mass: Fraction
    components: int
    closed_points: int
    exit_status: str
    contact_status: tuple[str, ...]


def audit(row: tuple[int, ...]) -> RowAudit:
    comps = safe_open_components(row)
    mass = interval_measure(comps)
    closed = threshold_safe_points(row)
    if mass == 0 and closed:
        status = "boundary_tight"
    elif mass > 0:
        status = "strict_open"
    else:
        status = "covered_or_debt"
    return RowAudit(row, mass, len(comps), len(closed), status, unit_contact_status(row))


def frac_word(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


KeyFunc = Callable[[tuple[int, ...]], object]


def quotient_summary(rows: list[RowAudit], keyfunc: KeyFunc) -> dict[str, object]:
    fibers: dict[object, list[RowAudit]] = defaultdict(list)
    for item in rows:
        fibers[keyfunc(item.row)].append(item)

    mixed = [
        fiber
        for fiber in fibers.values()
        if len({item.exit_status for item in fiber}) > 1
    ]
    largest_mixed = max((len(f) for f in mixed), default=0)
    example = None
    for fiber in sorted(mixed, key=lambda f: (-len(f), [x.row for x in f])):
        boundary = [x for x in fiber if x.exit_status == "boundary_tight"]
        strict = [x for x in fiber if x.exit_status == "strict_open"]
        if boundary and strict:
            example = (boundary[0], strict[0])
            break

    return {
        "fibers": len(fibers),
        "mixed_fibers": len(mixed),
        "largest_mixed": largest_mixed,
        "example": example,
    }


def repair_summary(rows: list[RowAudit], quotient: KeyFunc, sidecar: KeyFunc) -> int:
    def joined(row: tuple[int, ...]) -> tuple[object, object]:
        return (quotient(row), sidecar(row))

    return int(quotient_summary(rows, joined)["mixed_fibers"])


@dataclass(frozen=True)
class Carrier:
    name: str
    preserves: frozenset[str]
    destroys: frozenset[str]
    priority: int


OBLIGATION_WEIGHT = {
    "exit_status": 10,
    "unit_contact_status": 8,
    "c3_binding_skeleton": 8,
    "quadratic_sidecar": 7,
    "covering_layer": 8,
    "height_flex": 9,
    "exactness_test": 9,
    "finite_checkability": 7,
}


CARRIERS = [
    Carrier(
        "height_completed_packet",
        frozenset(OBLIGATION_WEIGHT),
        frozenset(),
        0,
    ),
    Carrier(
        "height_flex_ledger",
        frozenset({"height_flex", "covering_layer", "exit_status", "finite_checkability"}),
        frozenset({"c3_binding_skeleton", "quadratic_sidecar"}),
        1,
    ),
    Carrier(
        "c3_plus_quadratic_field_packet",
        frozenset({"c3_binding_skeleton", "quadratic_sidecar", "exactness_test", "finite_checkability"}),
        frozenset({"height_flex", "covering_layer"}),
        2,
    ),
    Carrier(
        "contact_status_sidecar",
        frozenset({"unit_contact_status", "exit_status", "finite_checkability"}),
        frozenset({"height_flex", "quadratic_sidecar"}),
        3,
    ),
    Carrier(
        "c3_binding_skeleton_only",
        frozenset({"c3_binding_skeleton", "finite_checkability"}),
        frozenset({"quadratic_sidecar", "height_flex", "covering_layer"}),
        4,
    ),
    Carrier(
        "quadratic_character_only",
        frozenset({"quadratic_sidecar", "finite_checkability"}),
        frozenset({"c3_binding_skeleton", "height_flex", "covering_layer"}),
        5,
    ),
    Carrier(
        "raw_residue_table",
        frozenset({"finite_checkability"}),
        frozenset({"height_flex", "exactness_test", "unit_contact_status"}),
        6,
    ),
    Carrier(
        "raw_unit_projection",
        frozenset({"unit_contact_status"}),
        frozenset({"height_flex", "covering_layer", "quadratic_sidecar", "exit_status"}),
        7,
    ),
]


def carrier_score(carrier: Carrier) -> int:
    keep = sum(OBLIGATION_WEIGHT[item] for item in carrier.preserves)
    lose = sum(OBLIGATION_WEIGHT[item] for item in carrier.destroys)
    return 2 * keep - lose


def beats(left: Carrier, right: Carrier) -> bool:
    left_key = (carrier_score(left), -len(left.destroys), -left.priority)
    right_key = (carrier_score(right), -len(right.destroys), -right.priority)
    return left_key > right_key


def tournament_fingerprint() -> dict[str, object]:
    graph = {c.name: set() for c in CARRIERS}
    for a, b in combinations(CARRIERS, 2):
        if beats(a, b):
            graph[a.name].add(b.name)
        else:
            graph[b.name].add(a.name)

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

    ordered = sorted(CARRIERS, key=lambda c: (-carrier_score(c), len(c.destroys), c.priority))
    return {
        "vertices": len(CARRIERS),
        "score_hist": dict(sorted(Counter(carrier_score(c) for c in CARRIERS).items())),
        "directed_3cycles": cycles,
        "hamiltonian_path_count": hpaths,
        "priority_path": [c.name for c in ordered],
    }


def named_rows() -> list[tuple[str, tuple[int, ...]]]:
    return [
        ("AP", AP),
        ("GW_12_to_24", replace_speed(AP, 12, 24)),
        ("positive_near_miss_12_to_36", replace_speed(AP, 12, 36)),
        ("petal_10_to_20", replace_speed(AP, 10, 20)),
        ("same_residue_height_2_to_16", replace_speed(AP, 2, 16)),
        ("same_unit_height_13_to_27", replace_speed(AP, 13, 27)),
    ]


def print_example(pair: tuple[RowAudit, RowAudit] | None) -> str:
    if pair is None:
        return "none"
    a, b = pair
    return (
        f"{a.row} {a.exit_status} mass={frac_word(a.mass)} closed={a.closed_points}"
        f"  ||  {b.row} {b.exit_status} mass={frac_word(b.mass)} closed={b.closed_points}"
    )


def main() -> None:
    audited = [audit(row) for row in one_swap_collar(84)]
    status_counts = Counter(item.exit_status for item in audited)
    boundary_rows = [item for item in audited if item.exit_status == "boundary_tight"]
    strict_rows = [item for item in audited if item.exit_status == "strict_open"]

    print("HYP-3401 LRC14 three-coordinate obstruction exactness scout")
    print("=" * 78)
    print("claim=first lost quotient coordinate should be an exact/sidecar obstruction")
    print("status=finite exact AP-collar scout; not an LRC14 proof")
    print("packet=C3 binding skeleton + Qsqrt(-7) character + height/flex ledger")
    print()

    print("BANK")
    print(f"one_swap_collar_max_new_speed=84")
    print(f"rows={len(audited)}")
    print(f"exit_status_counts={dict(sorted(status_counts.items()))}")
    print(f"boundary_rows={len(boundary_rows)}")
    print(f"strict_open_rows={len(strict_rows)}")
    print()

    print("NAMED ROWS")
    for name, row in named_rows():
        item = audit(row)
        print(
            f"  {name}: row={row} exit={item.exit_status} "
            f"mass={frac_word(item.mass)} components={item.components} "
            f"closed_points={item.closed_points} contacts={''.join(item.contact_status)} "
            f"c3={c3_skeleton_key(row)} qsqrt={quadratic_key(row)} "
            f"cover={covering_layer_key(row)} unit_height={unit_height_key(row)} "
            f"nonunit_height={nonunit_height_key(row)}"
        )
    print()

    quotients: list[tuple[str, KeyFunc]] = [
        ("raw_unit_projection", unit_projection_key),
        ("raw_mod14_residue_table", residue_counts),
        ("c3_binding_skeleton", c3_skeleton_key),
        ("quadratic_Qsqrt_minus7_character", quadratic_key),
        ("c3_plus_quadratic", lambda row: (c3_skeleton_key(row), quadratic_key(row))),
        (
            "c3_plus_quadratic_plus_covering_layer",
            lambda row: (c3_skeleton_key(row), quadratic_key(row), covering_layer_key(row)),
        ),
        ("c3_quadratic_nonunit_height_packet", c3_quadratic_nonunit_height_key),
        ("height_completed_packet", height_completed_packet_key),
        ("full_height_residue_ledger", height_flex_key),
    ]

    sidecars: list[tuple[str, KeyFunc]] = [
        ("unit_contact_status", unit_contact_status),
        ("covering_layer", covering_layer_key),
        ("unit_height_flex", unit_height_key),
        ("nonunit_height_flex", nonunit_height_key),
        ("full_height_flex", height_flex_key),
        ("c3_quadratic_nonunit_height_packet", c3_quadratic_nonunit_height_key),
        ("height_completed_packet", height_completed_packet_key),
    ]

    print("QUOTIENT MIXING AUDIT")
    summaries = {}
    for name, keyfunc in quotients:
        summary = quotient_summary(audited, keyfunc)
        summaries[name] = summary
        print(
            f"  {name}: fibers={summary['fibers']} mixed_fibers={summary['mixed_fibers']} "
            f"largest_mixed={summary['largest_mixed']}"
        )
        print(f"    example={print_example(summary['example'])}")
    print()

    print("SIDECAR REPAIR AUDIT")
    for qname, qfunc in quotients[:6]:
        repairs = []
        for sname, sfunc in sidecars:
            repairs.append((sname, repair_summary(audited, qfunc, sfunc)))
        best = min(value for _, value in repairs)
        repair_word = ", ".join(f"{name}:{value}" for name, value in repairs)
        print(f"  {qname}: best_mixed_after_sidecar={best}; {repair_word}")
    print()

    nonunit_packet_mixed = summaries["c3_quadratic_nonunit_height_packet"]["mixed_fibers"]
    height_completed_mixed = summaries["height_completed_packet"]["mixed_fibers"]
    full_height_mixed = summaries["full_height_residue_ledger"]["mixed_fibers"]
    print("READOUT")
    print(
        "  mixed fibers are concrete first-obstruction cocycles: the quotient "
        "cannot be proof currency until a sidecar kills the mixing."
    )
    print(
        f"  c3_quadratic_nonunit_height_packet_mixed_fibers={nonunit_packet_mixed}; "
        f"height_completed_packet_mixed_fibers={height_completed_mixed}; "
        f"full_height_residue_ledger_mixed_fibers={full_height_mixed}"
    )
    if height_completed_mixed == 0:
        print("  AP-collar result: C3 + quadratic + full height/flex separates boundary-tight from strict-open rows.")
    else:
        print("  AP-collar warning: even the height-completed packet has residual mixed debt.")
    if nonunit_packet_mixed:
        print("  new sidecar debt: nonunit height alone misses unit-height lifts such as 13->27.")
    print()

    print("TOURNAMENT ANALYSIS")
    fp = tournament_fingerprint()
    print("pairwise_observable=retained exit/status coordinates minus destroyed first-obstruction payload")
    print("binary_gauge=A->B iff weighted retained proof payload is larger; ties use fewer destroyed sidecars")
    for key, value in fp.items():
        if key == "priority_path":
            print(f"{key}={' -> '.join(value)}")
        else:
            print(f"{key}={value}")

    print()
    print("NEXT PROOF TARGET")
    print("  Formalize this AP-collar computation as the first finite lemma:")
    print("  any quotient fiber mixing boundary-tight and strict-open rows has a nonzero first obstruction.")
    print("  Partial sidecars locate the obstruction, but the AP-collar repair that kills all mixing is")
    print("  the full C3 + Qsqrt(-7) + all-residue height/flex packet.")


if __name__ == "__main__":
    main()
