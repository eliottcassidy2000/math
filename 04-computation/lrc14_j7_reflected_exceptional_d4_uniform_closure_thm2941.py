#!/usr/bin/env python3
r"""Exact referee proof that the entire reflected ``D=4`` sector closes.

The universal chromatic theorem attaches a good graph to each six-label body.
For 3001 of the 3003 bodies this graph is ``K6``.  Six labels taking values in
the five consecutive levels ``m,...,m+4`` therefore have a repeated good edge,
which closes immediately.  Only the two exceptional bodies can remain:

    E1=(1,2,7,9,11,13),   bad path 7--9--11--13;
    E2=(2,4,7,9,11,13),   bad matching {7--11,9--13}.

After removing words with a repeated good edge, exact enumeration leaves 432
and 312 proper labelled words over the offsets ``0,...,4`` with both endpoints
used.  They use four or five colours.  In particular the two singleton labels
``(1,2)`` for E1 and ``(2,4)`` for E2 always have distinct levels ``p,q``.

The four-colour rows are a necessary guardrail.  A tempting scout inference
was that spread four forces all five offsets to occur.  It does not: the
proper words ``(0,1,2,2,4,4)`` and ``(0,1,2,4,2,4)`` are explicit witnesses,
and there are exactly 72 four-colour rows on each exceptional body.  Thus a
``Delta=1`` tail is false here; the proof below retains the correct
``Delta<=2`` bound.

Every proper word uses at least four offsets, including zero and four, so its
nearest positive offset ``Delta`` is at most two.  The nearest-level theorem
therefore closes the tail whenever

    m > 2*(13392/35) + 93992/3185 = 2531336/3185,

hence for every ``m>=795``.

For ``m<=794`` the singleton levels satisfy ``p,q<=798``.  Put
``g=gcd(p,q)`` and sort ``P=p/g,Q=q/g``.  When ``P+Q>=8``, the exact
arbitrary-phase periodized-trapezoid floor from THM-1226(33) is

    min(1/(7Q), (P+Q-7)/(14PQ)) >= 397/2226021.

The displayed minimum is attained at ``(p,q)=(797,798)``.  Notice that the
tempting ``1/(7*798)`` is slightly too large; the exact finite audit below
keeps the smaller sharp value.  Replacing the integer skeleton by the two
reflected singleton clauses costs at most

    4(a+b)/L = 1/21021

on both exceptional bodies.  Subtracting this cost and the universal excess
debt still leaves a strict positive margin.

If ``P+Q<=7``, then ``g`` divides ``|p-q|<=4``.  Thus ``m>=14`` is impossible
because ``P+Q>=(2m+1)/4>7``; an exact boundary check eliminates ``m=13``.
The remaining finite endpoint/direct-union atlas has respectively 1272 and
924 rows, all at ``m<=12``, and every row has a positive singleton-overlap
margin after its exact six-clause debt.

Consequently both exceptional proper-word lanes close, and hence every
reflected ``D=4`` packet closes for every ``m>=1``.  This is a sufficient
certificate theorem inside the reflected THM-2941 residual family.  It is not
a classification of physical LRC survivors outside that family, and says
nothing about the still-open ``D>=5`` sectors.
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
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
UNIVERSAL_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.out"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_exceptional_d4_uniform_closure_thm2941.out"
)

EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_NEAREST_SHA256 = "367a4b299ebaf802faeedb056f2e3061b707c3df9aeebcb9e7afb941681cd750"
EXPECTED_UNIVERSAL_SHA256 = "a6f58c1a52dfc1fca61a239068dbe0b216bac41f1622b98748bc4a6d213fb6e8"
EXPECTED_UNIVERSAL_OUTPUT_SHA256 = "7364d5866171405fa90539a9ad76727c0c52f020ac1a104a1ab4f0276aedd115"
EXPECTED_SEMANTIC_SHA256 = ""

TAIL_START = 795
FINITE_M_MAX = TAIL_START - 1
STATIC_Q_MAX = FINITE_M_MAX + 4
STATIC_FLOOR = F(397, 2226021)
PERTURBATION_COST = F(1, 21021)
EXPECTED_CASE_COUNTS = (
    (432, ((4, 72), (5, 360)), ((1, 408), (2, 24)), 1272),
    (312, ((4, 72), (5, 240)), ((1, 288), (2, 24)), 924),
)
EXPECTED_FOUR_COLOUR_WITNESSES = (
    (0, 1, 2, 2, 4, 4),
    (0, 1, 2, 4, 2, 4),
)
EXPECTED_WORST = (
    (
        F(67358106244020066311889523979, 2401914229426762691547047386035),
        1,
        (2, 3, 4, 1, 0, 0),
        (3, 4, 5, 2, 1, 1),
        (3, 4),
        (1, 3, 4),
        F(10713241539, 381785167765),
        F(41807289679466787302696362, 2401914229426762691547047386035),
        146718,
        F(217432167327592319879772826618, 343130604203823241649578198005),
    ),
    (
        F(823316684250880241382046206698, 28824745742245431840091991187315),
        1,
        (4, 1, 2, 0, 2, 0),
        (5, 2, 3, 1, 3, 1),
        (5, 2),
        (1, 2, 5),
        F(36036, 1261259),
        F(248132783175715237137477562, 28824745742245431840091991187315),
        468468,
        F(62908110591836276488867223, 106110258982162392793980435),
    ),
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


require(sha256(BASE) == EXPECTED_BASE_SHA256, "reflected interval engine changed")
require(sha256(NEAREST) == EXPECTED_NEAREST_SHA256, "nearest-level referee changed")
require(sha256(UNIVERSAL) == EXPECTED_UNIVERSAL_SHA256, "universal chromatic referee changed")
require(
    sha256(UNIVERSAL_OUTPUT) == EXPECTED_UNIVERSAL_OUTPUT_SHA256,
    "universal chromatic stored output changed",
)
SPEC = spec_from_file_location("exceptional_d4_base", BASE)
require(SPEC is not None and SPEC.loader is not None, "cannot load reflected base")
R = module_from_spec(SPEC)
SPEC.loader.exec_module(R)


CASES = (
    (
        (1, 2, 7, 9, 11, 13),
        ((2, 3), (3, 4), (4, 5)),
        (0, 1),
    ),
    (
        (2, 4, 7, 9, 11, 13),
        ((2, 4), (3, 5)),
        (0, 1),
    ),
)


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


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


def excess_debt(body, ruler: int, levels) -> F:
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
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    # The nearest-level theorem with Delta<=2 begins at the sharp next integer.
    tail_wall = 2 * F(13392, 35) + F(93992, 3185)
    require(tail_wall == F(2531336, 3185), ("tail wall changed", tail_wall))
    require(tail_wall < TAIL_START and tail_wall >= TAIL_START - 1,
            ("tail integer threshold changed", tail_wall))

    # Audit the THM-1226 floor on the full singleton box needed in the head.
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
        static_minimum == (STATIC_FLOOR, 797, 798, 1, 797, 798),
        ("static floor changed", static_minimum),
    )

    all_edges = set(combinations(range(6), 2))
    case_rows = []
    for case_index, (body, bad_edges, singleton_pair) in enumerate(CASES):
        ruler, safe_ranges = R.safe_cell_ranges(body)
        candidate_cells = tuple(dict.fromkeys(
            endpoint for left, right in safe_ranges for endpoint in (left, right - 1)
        ))
        require(all(R.body_cell_is_safe(ruler, body, cell) for cell in candidate_cells),
                ("candidate endpoint is not body-safe", body))

        # Exact proper-word universe for D=4: offsets 0 and 4 both occur.
        good_edges = all_edges - set(bad_edges)
        proper = tuple(
            word
            for word in product(range(5), repeat=6)
            if min(word) == 0
            and max(word) == 4
            and all(word[i] != word[j] for i, j in good_edges)
        )
        colour_counts = tuple(sorted(Counter(len(set(word)) for word in proper).items()))
        nearest_counts = tuple(sorted(Counter(
            min(value for value in set(word) if value > 0) for word in proper
        ).items()))
        expected_count, expected_colours, expected_nearest, expected_dangerous = (
            EXPECTED_CASE_COUNTS[case_index]
        )
        require(len(proper) == expected_count, ("proper-word count changed", body, len(proper)))
        require(colour_counts == expected_colours,
                ("proper colour support changed", body, colour_counts))
        require(nearest_counts == expected_nearest,
                ("nearest positive offset changed", body, nearest_counts))
        require(all(word[singleton_pair[0]] != word[singleton_pair[1]] for word in proper),
                ("singleton levels ceased to be distinct", body))
        four_colour_witness = EXPECTED_FOUR_COLOUR_WITNESSES[case_index]
        require(
            four_colour_witness in proper
            and len(set(four_colour_witness)) == 4
            and min(four_colour_witness) == 0
            and max(four_colour_witness) == 4,
            ("spread-four four-colour guardrail failed", body, four_colour_witness),
        )

        a = body[singleton_pair[0]]
        b = body[singleton_pair[1]]
        perturbation = F(4 * (a + b), ruler)
        require(perturbation == PERTURBATION_COST,
                ("singleton perturbation changed", body, perturbation))
        debt_ceiling = universal_debt(body, ruler)
        static_margin = STATIC_FLOOR - perturbation - debt_ceiling
        require(static_margin > 0, ("static margin failed", body, static_margin))

        # The divisibility argument handles m>=14; audit the m=13 boundary.
        require(F(2 * 14 + 1, 4) > 7, "large-m dangerous-pattern inequality failed")
        for word in proper:
            p = 13 + word[singleton_pair[0]]
            q = 13 + word[singleton_pair[1]]
            _, P, Q = reduced_pattern(p, q)
            require(P + Q >= 8, ("unexpected m=13 dangerous pattern", body, word, P, Q))

        # Enumerate the entire finite head, not merely the analytically possible
        # first twelve bases, so the printed count inherits the stated universe.
        dangerous_inputs = []
        for m in range(1, FINITE_M_MAX + 1):
            for word in proper:
                levels = tuple(m + value for value in word)
                p = levels[singleton_pair[0]]
                q = levels[singleton_pair[1]]
                divisor, P, Q = reduced_pattern(p, q)
                if P + Q <= 7:
                    dangerous_inputs.append((m, word, levels, p, q, divisor, P, Q))
        require(len(dangerous_inputs) == expected_dangerous,
                ("dangerous-row count changed", body, len(dangerous_inputs)))
        require(max(row[0] for row in dangerous_inputs) == 12,
                ("dangerous m boundary changed", body, max(row[0] for row in dangerous_inputs)))

        # Cache exact pair profiles: many proper words induce the same ordered
        # singleton levels, while their six-clause debts remain distinct.
        pair_profiles = {}
        for _, _, _, p, q, _, _, _ in dangerous_inputs:
            if (p, q) in pair_profiles:
                continue
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
            pair_profiles[(p, q)] = max(overlaps)

        dangerous = []
        for m, word, levels, p, q, divisor, P, Q in dangerous_inputs:
            debt = excess_debt(body, ruler, levels)
            best_overlap, best_cell = pair_profiles[(p, q)]
            margin = best_overlap - debt
            require(margin > 0, (
                "dangerous endpoint pair failed",
                body,
                m,
                word,
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
                (
                    margin,
                    m,
                    word,
                    levels,
                    (p, q),
                    (divisor, P, Q),
                    best_overlap,
                    debt,
                    best_cell,
                    direct_union,
                )
            )

        worst_dangerous = min(dangerous)
        require(worst_dangerous == EXPECTED_WORST[case_index],
                ("dangerous minimum changed", body, worst_dangerous))
        case_rows.append(
            (
                body,
                ruler,
                len(proper),
                colour_counts,
                nearest_counts,
                four_colour_witness,
                singleton_pair,
                (a, b),
                perturbation,
                debt_ceiling,
                static_margin,
                len(dangerous),
                len(pair_profiles),
                len(candidate_cells),
                worst_dangerous,
            )
        )

    semantic = hashlib.sha256(repr((tail_wall, static_minimum, case_rows)).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest changed", semantic))
    source_sha = sha256(Path(__file__))
    lines = [
        "LRC14 exceptional-body reflected D=4 uniform closure exact proof",
        "complete_good_graph_bodies=3001;automatic_pigeonhole_closure=six_labels_into_five_levels",
        f"tail_Delta_max=2;tail_wall={qtext(tail_wall)};tail_start={TAIL_START}",
        f"finite_m_max={FINITE_M_MAX};singleton_q_max={STATIC_Q_MAX};"
        f"static_safe_pairs_checked={len(static_rows)};static_floor={qtext(STATIC_FLOOR)};"
        f"static_minimum_witness={static_minimum[1:]}",
    ]
    for row in case_rows:
        (
            body,
            ruler,
            proper_count,
            colour_counts,
            nearest_counts,
            four_colour_witness,
            singleton_pair,
            labels,
            perturbation,
            debt_ceiling,
            static_margin,
            dangerous_count,
            pair_profile_count,
            endpoint_count,
            worst,
        ) = row
        lines.append(
            f"CASE;E={body};L={ruler};proper_words={proper_count};"
            f"colour_counts={colour_counts};nearest_positive_counts={nearest_counts};"
            f"four_colour_spread4_witness={four_colour_witness};"
            f"singleton_indices={singleton_pair};singleton_labels={labels};"
            f"perturbation={qtext(perturbation)};debt_ceiling={qtext(debt_ceiling)};"
            f"static_margin={qtext(static_margin)};dangerous_rows={dangerous_count};"
            f"pair_profiles={pair_profile_count};endpoint_cells={endpoint_count};"
            f"worst_dangerous_margin={qtext(worst[0])};worst_m={worst[1]};"
            f"worst_word={worst[2]};worst_levels={worst[3]};worst_pair={worst[4]};"
            f"worst_reduced={worst[5]};worst_overlap={qtext(worst[6])};"
            f"worst_debt={qtext(worst[7])};worst_j={worst[8]};"
            f"worst_direct_union={qtext(worst[9])}"
        )
    lines.extend((
        "conclusion=both exceptional proper-word D=4 lanes close for every m>=1",
        "corollary=the entire reflected D=4 sector closes for every m>=1",
        "guardrail=spread four does not force all five offsets:72 proper four-colour rows occur on each exceptional body;Delta can equal 2",
        "scope=reflected THM-2941 residual family only; sufficient certificate, not physical-survivor classification; D>=5 remains outside this theorem",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_sha256={sha256(BASE)}",
        f"nearest_sha256={sha256(NEAREST)}",
        f"universal_sha256={sha256(UNIVERSAL)}",
        f"universal_output_sha256={sha256(UNIVERSAL_OUTPUT)}",
        f"source_sha256={source_sha}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ))
    output = "\n".join(lines) + "\n"
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
