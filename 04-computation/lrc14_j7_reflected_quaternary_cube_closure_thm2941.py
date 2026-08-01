#!/usr/bin/env python3
"""Finite-exact closure of the full reflected level cube q_e in {1,2,3,4}.

For every six-body root E subset {1,...,14}, this audits all 4^6 labelled
packets Z={q_e L-e}.  The canonical all-q selector is
tried first.  Any selector crossing is repaired by searching body-safe cells;
the displayed repair is checked both by interval merging and an independent
endpoint-slab sweep.  All arithmetic is Fraction exact and checks survive -O.
"""

from __future__ import annotations

import argparse
import hashlib
import multiprocessing as mp
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation" / "lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
OUTPUT_PATH = ROOT / "05-knowledge" / "results" / "lrc14_j7_reflected_quaternary_cube_closure_thm2941.out"
EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_BODY_COUNT = 3003
EXPECTED_LEVEL_WORD_COUNT = 4096
EXPECTED_PACKET_COUNT = 12_300_288
EXPECTED_CROSSING_COUNT = 30
EXPECTED_SEMANTIC_SHA256 = "58d4dfcf4c8b9c39f5e473ca7094374ceed928fd983ab20843ea8f9b59248326"
THRESHOLD = Q(6, 7)
LEVELS = tuple(product((1, 2, 3, 4), repeat=6))


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


SPEC = spec_from_file_location("reflected_all_q_ternary_base", BASE)
require(SPEC is not None and SPEC.loader is not None, "cannot import reflected engine")
require(hashlib.sha256(BASE.read_bytes()).hexdigest() == EXPECTED_BASE_SHA256, "reflected base changed")
R = module_from_spec(SPEC)
SPEC.loader.exec_module(R)


def qtext(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def interval_table(E: tuple[int, ...], L: int, j: int):
    table = {}
    checks = 0
    for i, e in enumerate(E):
        for q in (1, 2, 3, 4):
            direct = R.direct_multiplier_arcs(L, q * L - e, j)
            reflected = R.reflected_level_arcs(L, e, q, j)
            require(direct == reflected, ("arc mismatch", E, j, e, q))
            table[i, q] = direct
            checks += 1
    return table, checks


def mass_from_table(table, levels: tuple[int, ...]) -> Q:
    arcs = tuple(arc for i, q in enumerate(levels) for arc in table[i, q])
    return R.interval_mass(R.merge_intervals(arcs))


def body_profile(E: tuple[int, ...]):
    L, safe = R.safe_cell_ranges(E)
    component = R.selected_component(E)
    require(component < len(safe), ("missing selected component", E, safe))
    j = safe[component][0]
    require(R.body_cell_is_safe(L, E, j), ("selector not safe", E, L, j))
    table, checks = interval_table(E, L, j)
    crossings = []
    maximum = (Q(-1), ())
    for levels in LEVELS:
        mass = mass_from_table(table, levels)
        maximum = max(maximum, (mass, levels))
        if mass >= THRESHOLD:
            crossings.append((levels, mass))
    return E, L, j, component, maximum, tuple(crossings), checks


def first_repair(case):
    E, levels, selected_mass = case
    L, safe = R.safe_cell_ranges(E)
    checked = 0
    for left, right in safe:
        for j in range(left, right):
            arcs = tuple(
                arc
                for e, q in zip(E, levels)
                for arc in R.direct_multiplier_arcs(L, q * L - e, j)
            )
            mass = R.interval_mass(R.merge_intervals(arcs))
            checked += 1
            if mass < THRESHOLD:
                slab = R.slab_sweep_mass(arcs)
                require(mass == slab, ("repair route mismatch", E, levels, j))
                reflected = tuple(
                    arc
                    for e, q in zip(E, levels)
                    for arc in R.reflected_level_arcs(L, e, q, j)
                )
                require(
                    R.merge_intervals(arcs) == R.merge_intervals(reflected),
                    ("reflected repair mismatch", E, levels, j),
                )
                return E, levels, selected_mass, L, j, mass, checked
    return E, levels, selected_mass, L, None, None, checked


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count() or 1))
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    bodies = tuple(combinations(range(1, 15), 6))
    require(len(bodies) == EXPECTED_BODY_COUNT, "body universe changed")
    require(len(LEVELS) == EXPECTED_LEVEL_WORD_COUNT, "level-word universe changed")
    with mp.get_context("spawn").Pool(args.workers) as pool:
        profiles = list(pool.imap(body_profile, bodies, chunksize=1))
    profiles.sort(key=lambda row: row[0])
    cases = tuple(
        (E, levels, mass)
        for E, _L, _j, _component, _maximum, crossings, _checks in profiles
        for levels, mass in crossings
    )
    if cases:
        with mp.get_context("spawn").Pool(min(args.workers, len(cases))) as pool:
            repairs = list(pool.imap(first_repair, cases, chunksize=1))
    else:
        repairs = []
    repairs.sort(key=lambda row: (row[0], row[1]))
    survivors = tuple(row for row in repairs if row[4] is None)
    require(len(cases) == EXPECTED_CROSSING_COUNT, ("selector crossing count changed", len(cases)))
    require(len(repairs) == len(cases) and not survivors, "a quaternary packet survived")
    maximum = max((row[4][0], row[0], row[4][1], row[2]) for row in profiles)
    payload = (
        len(bodies),
        len(LEVELS),
        len(bodies) * len(LEVELS),
        maximum,
        cases,
        tuple(repairs),
        survivors,
        sum(row[-1] for row in profiles),
    )
    semantic = hashlib.sha256(repr(payload).encode()).hexdigest()
    require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest changed", semantic))
    lines = [
        "LRC14 reflected full quaternary level-cube closure",
        f"all_q_engine_sha256={hashlib.sha256(BASE.read_bytes()).hexdigest()}",
        f"universe=bodies:{len(bodies)};level_words:{len(LEVELS)};packets:{len(bodies)*len(LEVELS)}",
        f"selected_crossings={len(cases)};repaired={len(repairs)-len(survivors)};survivors={len(survivors)}",
        f"selected_global_max={maximum}",
        f"direct_reflected_selector_checks={sum(row[-1] for row in profiles)}",
    ]
    for E, levels, selected_mass, L, j, mass, checked in repairs:
        lines.append(
            "REPAIR;E=" + ",".join(map(str, E))
            + ";q=" + ",".join(map(str, levels))
            + f";L={L};selected_mass={qtext(selected_mass)};"
            + (f"j={j};mass={qtext(mass)};gap={qtext(THRESHOLD-mass)};" if j is not None else "j=NONE;")
            + f"cells_checked={checked}"
        )
    lines.extend((
        "projection_consequence=for every packet some body-safe cell has union mass below 6/7;the projected residual has mass above 1/7 and cannot fit in one aligned danger comb",
        "scope_boundary=finite q_e in 1..4 only;not arbitrary reflected k=1 and not LRC14",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ))
    output = "\n".join(lines) + "\n"
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
