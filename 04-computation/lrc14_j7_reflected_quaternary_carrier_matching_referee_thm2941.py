#!/usr/bin/env python3
"""Carrier-averaged perfect-matching referee for reflected k=1 packets.

For body-safe L-cells J and reflected drifts z_i=q_i L-e_i, let A_i(j)
be the danger set in the normalized coordinate of cell j.  Put

    omega_ij = |J|^{-1} sum_j mu(A_i(j) cap A_j(j))
             = mu(C_E cap D_zi cap D_zj) / mu(C_E).

For every perfect matching M of K6, pointwise Hunter pairing gives

    mu(union_i A_i(j))
      <= sum_i mu(A_i(j))-sum_(ij in M) mu(A_i(j) cap A_j(j)).

Each singleton mass is 1/7+e_i/[7(q_i L-e_i)].  Consequently

    sum_(ij in M) omega_ij > sum_i e_i/[7(q_i L-e_i)]

forces the average union below 6/7, hence supplies a body-safe repair cell.
This exact referee applies the criterion to the thirty fixed-selector
exceptions in the full q_i in {1,2,3,4} cube.  It uses the canonical reduced-
period carrier-intersection engine and directly reconstructs the cell average
for two controls.  It is a sufficient matching certificate, not an iff test
and not an all-level reflected-k=1 theorem.
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CUBE_SOURCE = ROOT / "04-computation" / "lrc14_j7_reflected_quaternary_cube_closure_thm2941.py"
CUBE_OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_reflected_quaternary_cube_closure_thm2941.out"
PAIR_SOURCE = ROOT / "04-computation" / "lrc14_j7_critical_scalar_wall_independent_thm2941.py"
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_quaternary_carrier_matching_referee_thm2941.out"
)
EXPECTED_CASE_COUNT = 30
EXPECTED_MINIMUM = Q(
    1055296376296906628142003007621,
    17950320242912613039384288143938,
)
EXPECTED_MINIMUM_KEY = (
    (1, 3, 7, 9, 10, 11),
    (1, 1, 1, 2, 1, 2),
    ((0, 5), (1, 2), (3, 4)),
)
DIRECT_CONTROL_KEYS = (
    ((1, 2, 3, 9, 10, 12), (4, 1, 4, 1, 2, 2)),
    ((1, 3, 7, 9, 10, 11), (1, 1, 1, 2, 1, 2)),
)
EXPECTED_SEMANTIC_SHA256 = (
    "a65796c9d9cf8cfb29c693c18f7023c140d682385703ae3d51f196c4b8d5b127"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qtext(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def load(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CUBE = load("reflected_quaternary_matching_cube", CUBE_SOURCE)
PAIR = load("reflected_quaternary_matching_pair", PAIR_SOURCE)


def parse_cases() -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    rows = []
    for line in CUBE_OUTPUT.read_text(encoding="utf-8").splitlines():
        if not line.startswith("REPAIR;"):
            continue
        fields = dict(piece.split("=", 1) for piece in line.split(";")[1:] if "=" in piece)
        rows.append(
            (
                tuple(map(int, fields["E"].split(","))),
                tuple(map(int, fields["q"].split(","))),
            )
        )
    result = tuple(rows)
    require(len(result) == len(set(result)) == EXPECTED_CASE_COUNT, "exception universe changed")
    return result


def perfect_matchings(vertices: tuple[int, ...]):
    if not vertices:
        return ((),)
    first = vertices[0]
    rows = []
    for index in range(1, len(vertices)):
        second = vertices[index]
        rest = vertices[1:index] + vertices[index + 1 :]
        for matching in perfect_matchings(rest):
            rows.append(((first, second),) + matching)
    return tuple(rows)


MATCHINGS = perfect_matchings(tuple(range(6)))
require(len(MATCHINGS) == 15, "K6 matching universe changed")


def interval_intersection(first, second):
    rows = []
    i = j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            rows.append((left, right))
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return tuple(rows)


def profile(case):
    E, levels = case
    L = 14 * lcm(*E)
    labels = tuple(q * L - e for e, q in zip(E, levels))
    carrier = PAIR.carrier_for(E)
    h = Q(sum(right - left for left, right in carrier), PAIR.RULER)
    require(h > 0, ("empty carrier", E))
    weights = {
        pair: PAIR.restricted_pair_intersection(
            carrier, labels[pair[0]], labels[pair[1]]
        )
        / h
        for pair in combinations(range(6), 2)
    }
    rows = tuple(
        (
            sum((weights[tuple(sorted(pair))] for pair in matching), Q(0)),
            matching,
        )
        for matching in MATCHINGS
    )
    best_weight, best_matching = min(rows, key=lambda row: (-row[0], row[1]))
    excess = sum(
        (Q(e, 7 * (q * L - e)) for e, q in zip(E, levels)),
        Q(0),
    )
    margin = best_weight - excess
    require(margin > 0, ("matching did not close", E, levels, best_weight, excess))
    return E, levels, L, labels, h, tuple(sorted(weights.items())), best_matching, best_weight, excess, margin


def direct_cell_control(row):
    E, levels, L, labels, h, _weights, matching, best_weight, _excess, _margin = row
    direct_L, safe = CUBE.R.safe_cell_ranges(E)
    require(direct_L == L, ("ruler mismatch", E, direct_L, L))
    cells = tuple(j for left, right in safe for j in range(left, right))
    require(Q(len(cells), L) == h, ("cell/carrier mass mismatch", E, len(cells), h))
    total = Q(0)
    for j in cells:
        arcs = {
            index: CUBE.R.merge_intervals(
                CUBE.R.direct_multiplier_arcs(L, labels[index], j)
            )
            for index in {vertex for pair in matching for vertex in pair}
        }
        total += sum(
            (
                CUBE.R.interval_mass(interval_intersection(arcs[first], arcs[second]))
                for first, second in matching
            ),
            Q(0),
        )
    average = total / len(cells)
    require(average == best_weight, ("cell/global matching mismatch", E, levels, average, best_weight))
    return E, levels, len(cells), matching, average


def render(rows, controls):
    minimum = min((row[9], row[0], row[1], row[6]) for row in rows)
    require(minimum[0] == EXPECTED_MINIMUM, ("minimum margin changed", minimum))
    require(minimum[1:] == EXPECTED_MINIMUM_KEY, ("minimum witness changed", minimum))
    require(tuple((row[0], row[1]) for row in controls) == DIRECT_CONTROL_KEYS, "control order changed")
    payload = (
        digest(CUBE_SOURCE),
        digest(CUBE_OUTPUT),
        digest(PAIR_SOURCE),
        tuple(rows),
        tuple(controls),
        minimum,
    )
    semantic = hashlib.sha256(repr(payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest changed", semantic))
    lines = [
        "LRC14 reflected quaternary carrier-averaged K6 matching referee",
        f"cube_source_sha256={digest(CUBE_SOURCE)}",
        f"cube_output_sha256={digest(CUBE_OUTPUT)}",
        f"pair_engine_sha256={digest(PAIR_SOURCE)}",
        "criterion=for perfect matching M,average_union<=6/7+epsilon-sum_(ij in M)omega_ij",
        "omega_ij=mu(C_E intersect D_zi intersect D_zj)/mu(C_E)=average_body_safe_cell_pair_overlap",
        f"exceptions={len(rows)};perfect_matchings_per_case={len(MATCHINGS)};closed={len(rows)};survivors=0",
        f"minimum_margin={qtext(minimum[0])};minimum_E={minimum[1]};minimum_q={minimum[2]};minimum_matching={minimum[3]}",
    ]
    for row in rows:
        lines.append(
            f"CASE;E={','.join(map(str,row[0]))};q={','.join(map(str,row[1]))};L={row[2]};"
            f"matching={row[6]};weight={qtext(row[7])};epsilon={qtext(row[8])};margin={qtext(row[9])}"
        )
    for control in controls:
        lines.append(
            f"DIRECT_CONTROL;E={','.join(map(str,control[0]))};q={','.join(map(str,control[1]))};"
            f"cells={control[2]};matching={control[3]};average={qtext(control[4])};routes=cell_sweep/reduced_period;agree"
        )
    lines.extend(
        (
            "consequence=all 30 fixed-selector exceptions have a body-safe cell with union mass below 6/7",
            "scope_boundary=sufficient averaged matching criterion;finite quaternary exception set only;not arbitrary k1",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    rows = tuple(profile(case) for case in parse_cases())
    by_key = {(row[0], row[1]): row for row in rows}
    controls = tuple(direct_cell_control(by_key[key]) for key in DIRECT_CONTROL_KEYS)
    text = render(rows, controls)
    args.output.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")


if __name__ == "__main__":
    main()
