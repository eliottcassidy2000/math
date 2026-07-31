#!/usr/bin/env python3
"""Finite-exact heterogeneous-level extension of the reflected k=1 stalk.

For a six-body root ``E`` with ``L=14*lcm(E)``, consider

    Z(E,q) = {q_e L-e : e in E}.

The all-level diagonal theorem treats ``q_e=q``.  Here every level lies in
``1..8`` and at most two coordinates differ from one.  There are exactly

    C(14,6) * (1 + 6*7 + C(6,2)*7^2) = 2,336,334

such labelled packets.  The level-one selector closes all but fourteen on
its chosen body-safe cell.  Those fourteen are genuine selector failures:
their union mass is strictly above ``6/7`` there.  Restoring cell location
closes every exception; an exhaustive exact scan of every body-safe cell
finds a strict sub-``6/7`` cell in all fourteen cases.

Thus this is both a scoped closure theorem and a validity-gate witness: the
common-q selector is not monotone under heterogeneous levels, but movable
cell selection repairs the first two-coordinate level box.  All interval
arithmetic is Fraction exact.  Direct multiplier arcs are compared with the
independent reflected formula, and every exceptional union is recomputed by
an endpoint-slab sweep.  Assertions remain active under ``python -O``.
"""

from __future__ import annotations

import argparse
import hashlib
import multiprocessing as mp
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_two_coordinate_q8_mass_closure_thm2941.out"
)
EXPECTED_BASE_SHA256 = (
    "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
)

LEVEL_MAX = 8
DEVIATION_MAX = 2
THRESHOLD = Q(6, 7)
EXPECTED_VECTOR_COUNT = 778
EXPECTED_PACKET_COUNT = 2_336_334
EXPECTED_CROSSING_KEYS = (
    ((1, 2, 4, 6, 9, 11), (2, 1, 1, 1, 1, 1)),
    ((1, 2, 4, 6, 9, 12), (2, 1, 1, 1, 1, 1)),
    ((1, 2, 4, 6, 10, 12), (2, 1, 1, 1, 1, 1)),
    ((1, 2, 4, 9, 11, 12), (2, 1, 1, 1, 1, 2)),
    ((1, 2, 5, 7, 9, 12), (2, 1, 1, 2, 1, 1)),
    ((1, 3, 4, 6, 10, 12), (2, 2, 1, 1, 1, 1)),
    ((1, 3, 6, 8, 9, 10), (1, 1, 1, 1, 2, 1)),
    ((1, 3, 7, 9, 10, 11), (1, 1, 1, 2, 1, 2)),
    ((1, 3, 8, 9, 10, 11), (1, 1, 1, 2, 1, 2)),
    ((1, 3, 8, 9, 10, 12), (1, 1, 1, 2, 1, 2)),
    ((1, 4, 5, 6, 8, 11), (1, 1, 2, 1, 1, 1)),
    ((1, 4, 5, 6, 8, 12), (1, 1, 2, 1, 1, 1)),
    ((1, 4, 5, 8, 11, 12), (1, 1, 2, 1, 1, 2)),
    ((1, 5, 6, 7, 8, 12), (1, 2, 1, 2, 1, 1)),
)
EXPECTED_PROFILE_SHA256 = (
    "e5394ed809cee45624e8676c6b243db2080f375d4a5c2dde6744ca6d5de833c1"
)
EXPECTED_REPAIR_SHA256 = (
    "fec9966fb452c7868123a59ebf9c2d6603c2b7028701590f267cdacf9732bfee"
)
EXPECTED_SEMANTIC_SHA256 = (
    "8c104172b30a3bceaec3fb7f24a48f92a785cf573ab69931f8f1345258409d05"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qtext(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE) == EXPECTED_BASE_SHA256, "all-q reflected engine changed")
SPEC = spec_from_file_location("reflected_all_q_q8_base", BASE)
require(SPEC is not None and SPEC.loader is not None, "cannot load all-q engine")
R = module_from_spec(SPEC)
SPEC.loader.exec_module(R)


def level_vectors() -> tuple[tuple[int, ...], ...]:
    rows = [(1,) * 6]
    for count in range(1, DEVIATION_MAX + 1):
        for positions in combinations(range(6), count):
            for levels in product(range(2, LEVEL_MAX + 1), repeat=count):
                vector = [1] * 6
                for position, level in zip(positions, levels):
                    vector[position] = level
                rows.append(tuple(vector))
    result = tuple(rows)
    require(len(result) == len(set(result)) == EXPECTED_VECTOR_COUNT, "level box changed")
    return result


VECTORS = level_vectors()


def cell_table(E: tuple[int, ...], L: int, j: int):
    table = {}
    comparisons = 0
    digest = hashlib.sha256()
    for index, e in enumerate(E):
        for level in range(1, LEVEL_MAX + 1):
            multiplier = level * L - e
            direct = R.merge_intervals(R.direct_multiplier_arcs(L, multiplier, j))
            reflected = R.reflected_level_arcs(L, e, level, j)
            require(direct == reflected, ("reflected arc mismatch", E, L, j, e, level))
            mass = R.interval_mass(direct)
            require(
                mass == Q(level * L, 7 * multiplier),
                ("singleton mass law failed", E, L, j, e, level, mass),
            )
            table[index, level] = direct
            digest.update(f"{index}|{e}|{level}|{direct}|{mass}\n".encode())
            comparisons += 1
    return table, comparisons, digest.hexdigest()


def vector_mass(table, vector: tuple[int, ...]) -> Q:
    arcs = tuple(
        arc
        for index, level in enumerate(vector)
        for arc in table[index, level]
    )
    return R.interval_mass(R.merge_intervals(arcs))


def body_profile(E: tuple[int, ...]):
    L, safe = R.safe_cell_ranges(E)
    component = R.selected_component(E)
    require(component < len(safe), ("selected component missing", E, safe))
    j = safe[component][0]
    require(R.body_cell_is_safe(L, E, j), ("selected cell is not body-safe", E, L, j))
    table, comparisons, arc_digest = cell_table(E, L, j)
    values = []
    crossings = []
    increases = 0
    base = vector_mass(table, (1,) * 6)
    for vector in VECTORS:
        mass = vector_mass(table, vector)
        values.append((mass, vector))
        if mass > base:
            increases += 1
        if mass >= THRESHOLD:
            crossings.append((mass, vector))
    maximum = max(values)
    return (
        E,
        L,
        j,
        component,
        base,
        maximum,
        increases,
        tuple(crossings),
        comparisons,
        arc_digest,
    )


def repair_profile(case):
    E, vector, selected_mass = case
    L, safe = R.safe_cell_ranges(E)
    cells = tuple(j for left, right in safe for j in range(left, right))
    require(cells, ("body has no safe grid cell", E, L))
    rows = []
    route_comparisons = 0
    for j in cells:
        arcs = tuple(
            arc
            for e, level in zip(E, vector)
            for arc in R.direct_multiplier_arcs(L, level * L - e, j)
        )
        merged_mass = R.interval_mass(R.merge_intervals(arcs))
        slab_mass = R.slab_sweep_mass(arcs)
        require(merged_mass == slab_mass, ("repair routes disagree", E, vector, j))
        rows.append((merged_mass, j))
        route_comparisons += 1
    by_cell = {j: mass for mass, j in rows}
    require(len(by_cell) == len(cells), ("duplicate safe cell", E, L))
    require(
        all(by_cell[j] == by_cell[L - 1 - j] for j in cells),
        ("cell-reflection symmetry failed", E, vector),
    )
    selected_component = R.selected_component(E)
    selected_j = safe[selected_component][0]
    require(by_cell[selected_j] == selected_mass, ("selected crossing changed", E, vector))
    minimum = min(rows)
    maximum = max(rows)
    average = sum((mass for mass, _j in rows), Q(0)) / len(rows)
    below = sum(mass < THRESHOLD for mass, _j in rows)
    require(minimum[0] < THRESHOLD and below > 0, ("crossing has no repair cell", E, vector))
    return (
        E,
        vector,
        L,
        len(cells),
        selected_j,
        selected_mass,
        minimum,
        maximum,
        average,
        below,
        route_comparisons,
    )


def render(profiles, repairs) -> str:
    require(len(profiles) == comb(14, 6) == 3003, "body universe changed")
    require(len(profiles) * len(VECTORS) == EXPECTED_PACKET_COUNT, "packet count changed")
    crossings = tuple(
        sorted(
            (E, vector, mass)
            for E, _L, _j, _component, _base, _maximum, _increases, rows, _checks, _digest in profiles
            for mass, vector in rows
        )
    )
    crossing_keys = tuple((E, vector) for E, vector, _mass in crossings)
    require(crossing_keys == EXPECTED_CROSSING_KEYS, ("selector crossing set changed", crossing_keys))
    require(len(repairs) == len(crossings) == 14, "repair universe changed")
    require(
        tuple((row[0], row[1]) for row in repairs) == EXPECTED_CROSSING_KEYS,
        "repair order/universe changed",
    )
    require(all(row[6][0] < THRESHOLD for row in repairs), "a crossing was not repaired")
    direct_reflected_checks = sum(row[8] for row in profiles)
    require(direct_reflected_checks == 3003 * 6 * LEVEL_MAX, "arc check count changed")
    repair_route_checks = sum(row[-1] for row in repairs)
    global_maximum = max((row[5][0], row[0], row[5][1]) for row in profiles)
    q1_maximum = max((row[4], row[0]) for row in profiles)
    increased_bodies = sum(row[6] > 0 for row in profiles)
    profile_hash = hashlib.sha256(repr(tuple(profiles)).encode()).hexdigest()
    repair_hash = hashlib.sha256(repr(tuple(repairs)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    if EXPECTED_REPAIR_SHA256 is not None:
        require(repair_hash == EXPECTED_REPAIR_SHA256, "repair digest changed")
    semantic_payload = (
        LEVEL_MAX,
        DEVIATION_MAX,
        len(VECTORS),
        EXPECTED_PACKET_COUNT,
        crossing_keys,
        crossings,
        tuple((row[0], row[1], row[6], row[9]) for row in repairs),
        global_maximum,
        q1_maximum,
        increased_bodies,
        direct_reflected_checks,
        repair_route_checks,
        profile_hash,
        repair_hash,
    )
    semantic_hash = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 reflected heterogeneous two-coordinate q<=8 mass closure",
        f"all_q_engine_sha256={file_sha256(BASE)}",
        (
            "universe=six_body_roots:3003;levels:1..8;"
            "at_most_two_coordinates_differ_from_one"
        ),
        (
            f"level_vectors_per_body={len(VECTORS)};packets={EXPECTED_PACKET_COUNT};"
            f"selected_cell_increased_bodies={increased_bodies}"
        ),
        (
            f"selected_cell_crossings={len(crossings)};"
            f"direct_reflected_singleton_checks={direct_reflected_checks};"
            f"repair_merge_slab_checks={repair_route_checks}"
        ),
        f"selected_global_max={global_maximum};q1_global_max={q1_maximum}",
    ]
    repair_map = {(row[0], row[1]): row for row in repairs}
    for E, vector, selected_mass in crossings:
        row = repair_map[E, vector]
        lines.append(
            f"CROSSING;E={','.join(map(str, E))};q={','.join(map(str, vector))};"
            f"L={row[2]};safe_cells={row[3]};selected_j={row[4]};"
            f"selected_mass={qtext(selected_mass)};selected_excess="
            f"{qtext(selected_mass-THRESHOLD)};best_j={row[6][1]};"
            f"best_mass={qtext(row[6][0])};best_gap={qtext(THRESHOLD-row[6][0])};"
            f"good_cells={row[9]};average={qtext(row[8])};conclusion=CELL-REPAIRED"
        )
    lines.extend(
        (
            "selector_failures=14;repaired=14;survivors=0",
            f"profile_sha256={profile_hash}",
            f"repair_sha256={repair_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count() or 1))
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    bodies = tuple(combinations(range(1, 15), 6))
    if args.workers == 1:
        profiles = [body_profile(E) for E in bodies]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            profiles = list(pool.imap(body_profile, bodies, chunksize=1))
    profiles.sort(key=lambda row: row[0])
    crossing_cases = tuple(
        sorted(
            (E, vector, mass)
            for E, _L, _j, _component, _base, _maximum, _increases, rows, _checks, _digest in profiles
            for mass, vector in rows
        )
    )
    if args.workers == 1:
        repairs = [repair_profile(case) for case in crossing_cases]
    else:
        with mp.get_context("spawn").Pool(min(args.workers, len(crossing_cases))) as pool:
            repairs = list(pool.imap(repair_profile, crossing_cases, chunksize=1))
    repairs.sort(key=lambda row: (row[0], row[1]))
    output = render(tuple(profiles), tuple(repairs))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
