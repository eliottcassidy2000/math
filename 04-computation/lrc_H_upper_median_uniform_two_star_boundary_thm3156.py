#!/usr/bin/env python3
r"""Exact companion for THM-3156: a coherent uniform two-star boundary law.

Fix H=(1,2,3,4,6,12), ruler 168, body-safe cell 90, and the nine-edge
two-star on the first two labels.  For each 6<=D<=100 this referee considers
every injection of the six labels into [D,2D].  It minimizes the margin

    (1/9) sum_edge intersection_mass - exact_singleton_debt.

The exact minimization is used through D=12.  Thereafter directed rounding
gives an integer lower bound for 9*SCALE times the exact margin.  A top-four
matching lemma makes the four-leaf injective minimization exhaustive: in an
optimal matching each leaf uses one of its four cheapest available levels,
since at most three cheaper choices can be occupied by the other leaves.

This is a finite exact theorem inside the reflected THM-2941 sufficient
family.  It is not a physical-survivor theorem or a proof of LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import permutations
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
BASE = ROOT / "04-computation/lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
OUTPUT = ROOT / "05-knowledge/results/lrc_H_upper_median_uniform_two_star_boundary_thm3156.out"

EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_GLOBAL_MINIMUM = F(2942257539652583, 320491868062929255)
EXPECTED_GLOBAL_D = 6
EXPECTED_GLOBAL_ASSIGNMENT = (6, 9, 12, 11, 8, 7)
EXPECTED_CONSERVATIVE_MINIMUM = 85295963189446
EXPECTED_CONSERVATIVE_D = 24

H = (1, 2, 3, 4, 6, 12)
L = 168
CELL = 90
PIVOTS = (0, 1)
LEAVES = (2, 3, 4, 5)
EDGES = ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
         (1, 2), (1, 3), (1, 4), (1, 5))
MIN_D = 6
MAX_D = 100
EXACT_THROUGH = 12
SCALE = 10**15


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def floor_scaled(value: F) -> int:
    return value.numerator * SCALE // value.denominator


def ceil_scaled(value: F) -> int:
    return -((-value.numerator * SCALE) // value.denominator)


def intersection_mass(first, second) -> F:
    i = 0
    j = 0
    total = F(0)
    while i < len(first) and j < len(second):
        total += max(
            F(0),
            min(first[i][1], second[j][1])
            - max(first[i][0], second[j][0]),
        )
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return total


def four_leaf_minimum(choices, costs):
    """Exhaust the at most 4^4 top-four choices with injectivity."""
    best = None
    for q2 in choices[0]:
        for q3 in choices[1]:
            if q3 == q2:
                continue
            for q4 in choices[2]:
                if q4 == q2 or q4 == q3:
                    continue
                for q5 in choices[3]:
                    if q5 == q2 or q5 == q3 or q5 == q4:
                        continue
                    assignment = (q2, q3, q4, q5)
                    value = sum(costs[leaf, q] for leaf, q in zip(LEAVES, assignment))
                    row = (value, assignment)
                    if best is None or row < best:
                        best = row
    require(best is not None, ("empty top-four matching", choices))
    return best


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    require(sha256(BASE) == EXPECTED_BASE_SHA256,
            ("proved arc engine changed", sha256(BASE), EXPECTED_BASE_SHA256))
    R = import_module("uniform_two_star_proved_arc_engine", BASE)
    ruler, safe = R.safe_cell_ranges(H)
    require(ruler == L and R.body_cell_is_safe(L, H, CELL),
            ("unsafe fixed cell", ruler, safe, CELL))

    all_levels = tuple(range(MIN_D, 2 * MAX_D + 1))
    arcs = {}
    debts = {}
    rounded_debts = {}
    for i, e in enumerate(H):
        for q in all_levels:
            reflected = R.reflected_level_arcs(L, e, q, CELL)
            direct = R.direct_multiplier_arcs(L, q * L - e, CELL)
            expected_mass = F(q * L, 7 * (q * L - e))
            require(reflected == direct, ("arc routes", i, e, q))
            require(R.interval_mass(reflected) == expected_mass,
                    ("singleton mass", i, e, q))
            arcs[i, q] = reflected
            debt = F(e, 7 * (q * L - e))
            debts[i, q] = debt
            rounded_debts[i, q] = ceil_scaled(debt)

    gain_cache = {}
    floor_cache = {}

    def gain(i: int, j: int, p: int, q: int) -> F:
        key = (i, j, p, q)
        value = gain_cache.get(key)
        if value is None:
            value = intersection_mass(arcs[i, p], arcs[j, q])
            merge_value = (
                R.interval_mass(arcs[i, p])
                + R.interval_mass(arcs[j, q])
                - R.interval_mass(R.merge_intervals(arcs[i, p] + arcs[j, q]))
            )
            require(value == merge_value, ("intersection routes", key, value, merge_value))
            gain_cache[key] = value
        return value

    def floor_gain(i: int, j: int, p: int, q: int) -> int:
        key = (i, j, p, q)
        value = floor_cache.get(key)
        if value is None:
            value = floor_scaled(gain(i, j, p, q))
            floor_cache[key] = value
        return value

    rows = []
    data_digest = hashlib.sha256()
    for D in range(MIN_D, MAX_D + 1):
        levels = tuple(range(D, 2 * D + 1))
        exact = D <= EXACT_THROUGH
        best = None
        for q0, q1 in permutations(levels, 2):
            if exact:
                base = gain(0, 1, q0, q1) / 9 - debts[0, q0] - debts[1, q1]
            else:
                base = (
                    floor_gain(0, 1, q0, q1)
                    - 9 * (rounded_debts[0, q0] + rounded_debts[1, q1])
                )
            costs = {}
            choices = []
            for leaf in LEAVES:
                ranked = []
                for q in levels:
                    if q in (q0, q1):
                        continue
                    if exact:
                        value = (
                            gain(0, leaf, q0, q) + gain(1, leaf, q1, q)
                        ) / 9 - debts[leaf, q]
                    else:
                        value = (
                            floor_gain(0, leaf, q0, q)
                            + floor_gain(1, leaf, q1, q)
                            - 9 * rounded_debts[leaf, q]
                        )
                    costs[leaf, q] = value
                    ranked.append((value, q))
                choices.append(tuple(q for _, q in sorted(ranked)[:4]))
            leaf_value, leaf_assignment = four_leaf_minimum(tuple(choices), costs)
            row = (base + leaf_value, (q0, q1) + leaf_assignment)
            if best is None or row < best:
                best = row
        require(best is not None and best[0] > 0, ("nonpositive or empty", D, best))
        score, assignment = best
        gain_sum = sum(
            (gain(i, j, assignment[i], assignment[j]) for i, j in EDGES),
            F(0),
        )
        debt_sum = sum((debts[i, assignment[i]] for i in range(6)), F(0))
        representative_margin = gain_sum / 9 - debt_sum
        if exact:
            require(score == representative_margin, ("exact score", D, score, representative_margin))
            kind = "exact_minimum"
        else:
            rounded_score = (
                sum(floor_gain(i, j, assignment[i], assignment[j]) for i, j in EDGES)
                - 9 * sum(rounded_debts[i, assignment[i]] for i in range(6))
            )
            require(score == rounded_score <= 9 * SCALE * representative_margin,
                    ("rounding direction", D, score, rounded_score, representative_margin))
            kind = "conservative_floor_for_9S_margin"
        line = (
            f"D={D};m={D};pivot_checks={D * (D + 1)};{kind}={qtext(score) if exact else score};"
            f"representative_margin={qtext(representative_margin)};levels={assignment};"
            f"nine_edge_gain_sum={qtext(gain_sum)};debt={qtext(debt_sum)}"
        )
        rows.append((D, exact, score, assignment, representative_margin, gain_sum, debt_sum, line))
        data_digest.update((line + "\n").encode())

    exact_global = min((row[2], row[0], row[3]) for row in rows if row[1])
    conservative_global = min((row[2], row[0]) for row in rows if not row[1])
    require(exact_global == (EXPECTED_GLOBAL_MINIMUM, EXPECTED_GLOBAL_D,
                             EXPECTED_GLOBAL_ASSIGNMENT), exact_global)
    require(conservative_global == (EXPECTED_CONSERVATIVE_MINIMUM,
                                    EXPECTED_CONSERVATIVE_D), conservative_global)
    conservative_margin_floor = F(EXPECTED_CONSERVATIVE_MINIMUM, 9 * SCALE)
    require(conservative_margin_floor > EXPECTED_GLOBAL_MINIMUM,
            (conservative_margin_floor, EXPECTED_GLOBAL_MINIMUM))

    lines = [
        "LRC H UNIFORM TWO-STAR BOUNDARY -- THM-3156 EXACT COMPANION",
        f"status=FINITE-EXACT;body={H};ruler={L};cell={CELL};edges={EDGES}",
        f"universe=D:{MIN_D}..{MAX_D};m=D;levels=[D,2D];six-label injections;rows={len(rows)}",
        f"law=uniform weight 1/9 on all nine fixed two-star edges;scale={SCALE}",
        "rounding=integer conservative score is a lower bound for 9*SCALE*((sum_edge I_edge)/9-debt)",
        f"exact_through={EXACT_THROUGH};global_exact_minimum={qtext(EXPECTED_GLOBAL_MINIMUM)}@D={EXPECTED_GLOBAL_D},levels={EXPECTED_GLOBAL_ASSIGNMENT}",
        f"smallest_conservative_9S_floor={EXPECTED_CONSERVATIVE_MINIMUM}@D={EXPECTED_CONSERVATIVE_D};converted={qtext(conservative_margin_floor)};strictly_above_global_exact_minimum",
        "consequence=uniform expected overlap exceeds full singleton debt, so some two-star edge does;pairwise Bonferroni puts the six-arc union strictly below 6/7 at cell 90",
        "scope=fixed reflected H packet and finite boundary D=6..100 inside THM-2941;not arbitrary D;not physical LRC14",
        "reproduction=python3 04-computation/lrc_H_upper_median_uniform_two_star_boundary_thm3156.py",
        *[row[-1] for row in rows],
        f"data_semantic_sha256={data_digest.hexdigest()}",
        f"proved_arc_engine_sha256={EXPECTED_BASE_SHA256}",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"source_sha256={sha256(HERE)}",
        "all_exact_controls=PASS",
    ]
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
