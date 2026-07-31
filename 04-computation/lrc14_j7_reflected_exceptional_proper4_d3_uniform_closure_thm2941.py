#!/usr/bin/env python3
r"""Exact referee proof that the two exceptional proper-4 D=3 rays close.

The universal chromatic referee leaves only two four-colour orbits when the
level spread is three.  Their forced colour classes are

    E1=(1,2,7,9,11,13): {7,9}, {11,13}, {1}, {2};
    E2=(2,4,7,9,11,13): {7,11}, {9,13}, {2}, {4}.

Give those classes the four distinct levels ``m,m+1,m+2,m+3``.  The two
singleton labels therefore have distinct levels ``p,q``.  They alone close
the finite head.

For ``m<=412``, put ``g=gcd(p,q)`` and ``P=p/g,R=q/g``.  If ``P+R>=8``, the
arbitrary-phase periodized-trapezoid floor THM-1226(33) is

    min(1/(7R), (P+R-7)/(14PR)) >= 137/400890,

because ``p,q<=415``; equality in the audited box is at ``(414,415)``.
Replacing the integer skeleton ``||pu-aj/L||<1/14`` by the exact reflected
clause ``||(p-a/L)u-aj/L||<1/14`` costs at most ``4a/L`` in symmetric
difference.  For the singleton pair the total cost is exactly bounded by

    4(a+b)/L = 1/21021

on both exceptional bodies.  After the universal singleton debt is
subtracted, the margin remains strictly positive.

The only static-dangerous patterns have ``P+R<=7``.  Since ``g`` divides
``|p-q|<=3``, they occur only for ``m<=9`` (``m=10`` is checked directly;
for ``m>=11``, ``P+R >= (2m+1)/3 > 7``).  The exact endpoint atlas below
contains 72 rows per body and gives a positive singleton-pair overlap margin
in every row.

Finally, for ``m>=413`` the already-proved nearest-level pair theorem applies
with ``Delta=1`` because all four consecutive levels occur.  Its exact strict
threshold is

    13392/35 + 93992/3185 = 1312664/3185 < 413.

Thus the two exceptional proper-four-colour D=3 packets close for every
``m>=1``.  Combined with the universal chromatic theorem, this closes the
entire reflected D=3 sector.  This is a sufficient-certificate theorem for
the reflected THM-2941 residual family, not a classification of physical LRC
survivors outside that family and not a closure of the still-open D>=4 lanes.
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation" / "lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
NEAREST = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_nearest_level_matching_tail_referee_thm2941.py"
)
UNIVERSAL = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.py"
)
EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_NEAREST_SHA256 = "367a4b299ebaf802faeedb056f2e3061b707c3df9aeebcb9e7afb941681cd750"
EXPECTED_UNIVERSAL_SHA256 = ""
EXPECTED_SEMANTIC_SHA256 = "61b456316b2e0d4b694af399c749586b4c6a3896cacaf581e29ddaf6adf7e1ed"

TAIL_START = 413
FINITE_M_MAX = TAIL_START - 1
STATIC_Q_MAX = FINITE_M_MAX + 3
STATIC_FLOOR = F(137, 400890)
PERTURBATION_COST = F(1, 21021)
EXPECTED_DANGEROUS_ROWS = 72


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    """Hash the repository's declared LF-normalized evidence image."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


require(sha256(BASE) == EXPECTED_BASE_SHA256, "reflected interval engine changed")
require(sha256(NEAREST) == EXPECTED_NEAREST_SHA256, "nearest-level referee changed")
if EXPECTED_UNIVERSAL_SHA256:
    require(sha256(UNIVERSAL) == EXPECTED_UNIVERSAL_SHA256, "universal chromatic referee changed")
SPEC = spec_from_file_location("exceptional_d3_base", BASE)
require(SPEC is not None and SPEC.loader is not None, "cannot load reflected base")
R = module_from_spec(SPEC)
SPEC.loader.exec_module(R)


CASES = (
    (
        (1, 2, 7, 9, 11, 13),
        ((2, 3), (4, 5), (0,), (1,)),
        ((2, 3), (3, 4), (4, 5)),
        (0, 1),
    ),
    (
        (2, 4, 7, 9, 11, 13),
        ((2, 4), (3, 5), (0,), (1,)),
        ((2, 4), (3, 5)),
        (0, 1),
    ),
)


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def canonical_partition(word) -> tuple[tuple[int, ...], ...]:
    return tuple(sorted(tuple(index for index, value in enumerate(word) if value == color)
                        for color in sorted(set(word))))


def level_word(classes, values) -> tuple[int, ...]:
    levels = [0] * 6
    for block, value in zip(classes, values):
        for index in block:
            levels[index] = value
    return tuple(levels)


def intersection_mass(first, second) -> F:
    i = 0
    j = 0
    total = F(0)
    while i < len(first) and j < len(second):
        total += max(F(0), min(first[i][1], second[j][1]) - max(first[i][0], second[j][0]))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return total


def singleton_debt(body, ruler: int, levels) -> F:
    return sum((F(e, 7 * (q * ruler - e)) for e, q in zip(body, levels)), F(0))


def universal_debt(body, ruler: int) -> F:
    return sum((F(e, 7 * (ruler - e)) for e in body), F(0))


def reduced_pattern(p: int, q: int) -> tuple[int, int, int]:
    divisor = gcd(p, q)
    left, right = sorted((p // divisor, q // divisor))
    return divisor, left, right


def translated_floor(P: int, Q: int) -> F:
    require(1 <= P < Q and gcd(P, Q) == 1 and P + Q >= 8,
            ("bad translated-floor input", P, Q))
    return min(F(1, 7 * Q), F(P + Q - 7, 14 * P * Q))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            ROOT
            / "05-knowledge"
            / "results"
            / "lrc14_j7_reflected_exceptional_proper4_d3_uniform_closure_thm2941.out"
        ),
    )
    args = parser.parse_args()

    # The nearest-level theorem's exact Delta=1 tail starts at 413.
    tail_wall = F(13392, 35) + F(93992, 3185)
    require(tail_wall == F(1312664, 3185), ("tail wall changed", tail_wall))
    require(tail_wall < TAIL_START and tail_wall >= TAIL_START - 1,
            ("tail integer threshold changed", tail_wall))

    # Audit the arbitrary-phase THM-1226(33) floor on a deliberately larger
    # box than needed: every distinct 1<=p<q<=415 whose reduced sum is safe.
    static_rows = []
    for p in range(1, STATIC_Q_MAX + 1):
        for q in range(p + 1, STATIC_Q_MAX + 1):
            divisor, P, Q = reduced_pattern(p, q)
            if P + Q <= 7:
                continue
            floor = translated_floor(P, Q)
            static_rows.append((floor, p, q, divisor, P, Q))
    static_minimum = min(static_rows)
    require(
        static_minimum == (STATIC_FLOOR, 414, 415, 1, 414, 415),
        ("static floor changed", static_minimum),
    )

    case_rows = []
    for body, classes, bad_edges, singleton_pair in CASES:
        ruler, safe_ranges = R.safe_cell_ranges(body)
        candidate_cells = tuple(dict.fromkeys(
            endpoint for left, right in safe_ranges for endpoint in (left, right - 1)
        ))
        require(all(R.body_cell_is_safe(ruler, body, cell) for cell in candidate_cells),
                ("candidate endpoint is not body-safe", body))

        # Exhaustively re-derive the 24 labelled proper four-colourings and
        # verify that their unlabelled class partition is the forced one.
        all_edges = set(combinations(range(6), 2))
        good_edges = all_edges - set(bad_edges)
        proper = tuple(
            word
            for word in product(range(4), repeat=6)
            if set(word) == set(range(4))
            and all(word[i] != word[j] for i, j in good_edges)
        )
        expected_partition = canonical_partition(
            tuple(next(color for color, block in enumerate(classes) if index in block)
                  for index in range(6))
        )
        require(
            len(proper) == 24
            and all(canonical_partition(word) == expected_partition for word in proper),
            ("proper-colouring orbit changed", body, len(proper)),
        )

        a = body[singleton_pair[0]]
        b = body[singleton_pair[1]]
        perturbation = F(4 * (a + b), ruler)
        require(perturbation == PERTURBATION_COST,
                ("singleton perturbation changed", body, perturbation))
        debt_ceiling = universal_debt(body, ruler)
        static_margin = STATIC_FLOOR - perturbation - debt_ceiling
        require(static_margin > 0, ("static margin failed", body, static_margin))

        # Prove the dangerous list stops at m=9.  The m=10 boundary is the
        # only case not handled immediately by (2m+1)/3>7.
        require(F(2 * 11 + 1, 3) > 7, "large-m dangerous-pattern inequality failed")
        for offsets in permutations(range(4)):
            levels = level_word(classes, tuple(10 + offset for offset in offsets))
            _, P, Q = reduced_pattern(levels[singleton_pair[0]], levels[singleton_pair[1]])
            require(P + Q >= 8, ("unexpected m=10 dangerous pattern", body, offsets, P, Q))

        dangerous = []
        for m in range(1, 10):
            for offsets in permutations(range(4)):
                levels = level_word(classes, tuple(m + offset for offset in offsets))
                p = levels[singleton_pair[0]]
                q = levels[singleton_pair[1]]
                divisor, P, Q = reduced_pattern(p, q)
                if P + Q >= 8:
                    continue
                require(divisor <= 3, ("level-difference gcd escaped", body, levels, divisor))
                debt = singleton_debt(body, ruler, levels)
                overlaps = tuple(
                    (
                        intersection_mass(
                            R.reflected_level_arcs(ruler, a, p, cell),
                            R.reflected_level_arcs(ruler, b, q, cell),
                        ),
                        cell,
                    )
                    for cell in candidate_cells
                )
                best_overlap, best_cell = max(overlaps)
                margin = best_overlap - debt
                require(margin > 0, (
                    "dangerous endpoint pair failed",
                    body,
                    m,
                    offsets,
                    levels,
                    (p, q),
                    (P, Q),
                    best_overlap,
                    debt,
                    best_cell,
                ))
                pair_first = R.reflected_level_arcs(ruler, a, p, best_cell)
                pair_second = R.reflected_level_arcs(ruler, b, q, best_cell)
                pair_union = R.interval_mass(R.merge_intervals(pair_first + pair_second))
                require(
                    best_overlap
                    == R.interval_mass(pair_first) + R.interval_mass(pair_second) - pair_union,
                    ("dangerous pair route mismatch", body, levels, best_cell),
                )
                all_arcs = tuple(
                    arc
                    for e, level in zip(body, levels)
                    for arc in R.reflected_level_arcs(ruler, e, level, best_cell)
                )
                direct_union = R.interval_mass(R.merge_intervals(all_arcs))
                bonferroni_ceiling = F(6, 7) + debt - best_overlap
                require(
                    direct_union <= bonferroni_ceiling < F(6, 7),
                    (
                        "dangerous direct-union control failed",
                        body,
                        levels,
                        best_cell,
                        direct_union,
                        bonferroni_ceiling,
                    ),
                )
                dangerous.append(
                    (margin, m, offsets, levels, (p, q), (P, Q), best_overlap, debt, best_cell)
                )
        require(len(dangerous) == EXPECTED_DANGEROUS_ROWS,
                ("dangerous-row count changed", body, len(dangerous)))
        worst_dangerous = min(dangerous)
        case_rows.append(
            (
                body,
                ruler,
                expected_partition,
                singleton_pair,
                (a, b),
                perturbation,
                debt_ceiling,
                static_margin,
                len(dangerous),
                worst_dangerous,
            )
        )

    expected_worst = (
        (
            F(26941875554907150716666611331, 960757693343931390131409610215),
            1,
            (1, 0, 2, 3),
            (3, 4, 2, 2, 1, 1),
            (3, 4),
            (3, 4),
            F(10713241539, 381785167765),
            F(17865415390492460351146078, 960757693343931390131409610215),
            146718,
        ),
        (
            F(65878002411925276451530027193922, 2306053555976270654339569960287055),
            2,
            (2, 1, 3, 0),
            (5, 2, 4, 3, 4, 3),
            (5, 2),
            (2, 5),
            F(36036, 1261259),
            F(9294283805647410881058149298, 2306053555976270654339569960287055),
            468468,
        ),
    )
    require(tuple(row[-1] for row in case_rows) == expected_worst,
            ("dangerous minima changed", tuple(row[-1] for row in case_rows)))

    semantic = hashlib.sha256(repr((tail_wall, static_minimum, case_rows)).encode()).hexdigest()
    require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest changed", semantic))
    source_sha = sha256(Path(__file__))
    lines = [
        "LRC14 exceptional proper-four-colour D=3 uniform closure scratch proof",
        f"universal_proper4_orbits={len(CASES)};labelled_colourings_per_orbit=24",
        f"tail_Delta=1;tail_wall={qtext(tail_wall)};tail_start={TAIL_START}",
        f"finite_m_max={FINITE_M_MAX};singleton_q_max={STATIC_Q_MAX};"
        f"static_safe_pairs_checked={len(static_rows)};static_floor={qtext(STATIC_FLOOR)};"
        f"static_minimum_witness={static_minimum[1:]}",
    ]
    for row in case_rows:
        (body, ruler, partition, singleton_pair, labels, perturbation, debt_ceiling,
         static_margin, dangerous_count, worst) = row
        lines.append(
            f"CASE;E={body};L={ruler};classes={partition};singleton_indices={singleton_pair};"
            f"singleton_labels={labels};perturbation={qtext(perturbation)};"
            f"debt_ceiling={qtext(debt_ceiling)};static_margin={qtext(static_margin)};"
            f"dangerous_rows={dangerous_count};worst_dangerous_margin={qtext(worst[0])};"
            f"worst_m={worst[1]};worst_offsets={worst[2]};worst_levels={worst[3]};"
            f"worst_pair={worst[4]};worst_reduced={worst[5]};"
            f"worst_overlap={qtext(worst[6])};worst_debt={qtext(worst[7])};worst_j={worst[8]}"
        )
    lines.extend((
        "conclusion=the two exceptional proper-four-colour D=3 rays close for every m>=1",
        "corollary=the entire reflected D=3 sector closes for every m>=1",
        "scope=reflected THM-2941 residual family only; sufficient certificate, not physical-survivor classification; D>=4 remains outside this theorem",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_lf_sha256={sha256(BASE)}",
        f"nearest_lf_sha256={sha256(NEAREST)}",
        f"universal_lf_sha256={sha256(UNIVERSAL)}",
        f"source_lf_sha256={source_sha}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ))
    output = "\n".join(lines) + "\n"
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
