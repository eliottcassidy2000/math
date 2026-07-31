#!/usr/bin/env python3
r"""Universal repeated-level closure from the signed-pair chromatic graph.

For a six-label body ``E`` let ``L=14*lcm(E)``.  The signed same-level pair
theorem in the pinned D=1 referee assigns some pairs ``a<b`` a body-safe cell
and a level-independent overlap floor

    P_Q > F_E(a,b)                       for every common level Q>=1.

Put

    eps_max(E) = sum_(e in E) e/[7(L-e)].

For arbitrary positive levels ``q_e``, the exact singleton excess satisfies
``eps(E,q)<=eps_max(E)`` term by term.  Call a pair *good* when
``F_E(a,b)>eps_max(E)`` and form the undirected graph ``G_E`` on the six
labels.  Whenever a good pair shares a level, one-edge Bonferroni gives

    mu(union A_e) <= 6/7+eps(E,q)-P_Q < 6/7,

so the THM-2941 reflected residual closes.

The exact census below computes all 3003 graphs and all 45045 labelled pairs.
For 3001 bodies, ``G_E=K6``.  The two exceptions are still 4-chromatic:

* ``E=(1,2,7,9,11,13)`` has bad-edge path
  ``7--9--11--13``.  Slots ``{0,1,2,4}`` give a good K4 and
  ``(0,1,2,2,3,3)`` is a proper four-colouring.
* ``E=(2,4,7,9,11,13)`` has bad matching ``{7--11,9--13}``.
  Slots ``{0,1,2,3}`` give a good K4 and ``(0,1,2,3,2,3)`` is a proper
  four-colouring.

The five bad edges are structural, not near-threshold accidents.  Their label
differences are respectively 2 and 4, and an exact 593992-cell incidence
census gives no body-safe cell meeting ``||difference*t||<1/7``.  Hence the
two corresponding same-level clauses are disjoint for every common level.
This explains why the chromatic boundary is genuinely four for this method.

Thus every assignment using at most three distinct levels has a monochromatic
good edge.  In particular the entire reflected ``D=2`` sector closes for all
``m>=1``.  The proved statement is stronger: any reflected k=1 packet with at
most three distinct positive levels closes, whether or not those levels are
consecutive.  Bodywise it is sharper still.  On the 3001 complete graphs, a
packet not closed by this one-pair mechanism must use six pairwise-distinct
levels.  On either exceptional body its level word must be a proper colouring
of the displayed good graph, hence uses at least four levels and can repeat
only along the displayed bad path or matching.  This is a sufficient
certificate reduction, not a physical survivor classification of those
proper-colouring residues.
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from math import lcm
from pathlib import Path


HERE = Path(__file__).resolve()
BASE_RELATIVE = Path("04-computation") / "lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
ROOT = next(parent for parent in HERE.parents if (parent / BASE_RELATIVE).is_file())
BASE = ROOT / BASE_RELATIVE
SIGNED_D1 = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_d1_signed_pair_complete_closure_thm2941.py"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.out"
)
EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
# Filled after the signed-D1 referee receives its final semantic pin.
EXPECTED_SIGNED_D1_SHA256 = "7650fe26c0047aa1c667b5be4fd0345eae556a647393d537cbfaefe56e0a09b2"
EXPECTED_SEMANTIC_SHA256 = "c51287141ba366c1c0ca12cfb15997b3821cf5371029f59ee57ff3ead40280d6"

BODY_COUNT = 3003
PAIR_COUNT = 45045
COMPLETE_GRAPH_COUNT = 3001
EXPECTED_CHI_COUNTS = ((4, 2), (6, 3001))
EXPECTED_EXCEPTIONS = (
    (
        (1, 2, 7, 9, 11, 13),
        ((2, 3), (3, 4), (4, 5)),
        (0, 1, 2, 4),
        (0, 1, 2, 2, 3, 3),
    ),
    (
        (2, 4, 7, 9, 11, 13),
        ((2, 4), (3, 5)),
        (0, 1, 2, 3),
        (0, 1, 2, 3, 2, 3),
    ),
)
EXPECTED_WEAKEST_GOOD = (
    F(321692483216665992817789205, 24976183956869361050998460388),
    (2, 3, 6, 7, 11, 13),
    (3, 5),
    (F(33, 2548), 68249, (-1, F(1089), F(4, 7), 5, 5)),
    F(4896776765901250909162, 68615889991399343546699067),
)
EXPECTED_BLIND_SUPPORT_CHECKS = 593992


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


require(sha256(BASE) == EXPECTED_BASE_SHA256, "reflected interval engine changed")
if EXPECTED_SIGNED_D1_SHA256:
    require(sha256(SIGNED_D1) == EXPECTED_SIGNED_D1_SHA256, "signed-D1 referee changed")
SPEC = spec_from_file_location("universal_pair_signed_d1", SIGNED_D1)
require(SPEC is not None and SPEC.loader is not None, "cannot load signed-D1 referee")
D1 = module_from_spec(SPEC)
SPEC.loader.exec_module(D1)


EDGES = tuple(combinations(range(6), 2))


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def universal_debt(body: tuple[int, ...], ruler: int) -> F:
    return sum((F(e, 7 * (ruler - e)) for e in body), F(0))


def best_q0_profile(
    body: tuple[int, ...],
    ruler: int,
    ranges: tuple[tuple[int, int], ...],
    i: int,
    k: int,
):
    """Best complete signed-chart floor at Q0=1 for one labelled pair."""
    a, b = body[i], body[k]
    cells = D1.candidate_cells(ruler, ranges, a, b)
    breakpoint_checks = D1.breakpoint_piece_audit(ruler, ranges, a, b, cells)
    ranked = []
    chart_count = 0
    for cell in cells:
        chart = D1.signed_chart(ruler, a, b, cell)
        if chart is None:
            continue
        chart_count += 1
        floor = D1.uniform_floor(ruler, 1, chart)
        if floor <= 0 or chart[1] + chart[2] <= 0:
            continue
        ranked.append((floor, cell, chart))
    return max(ranked, default=None), len(cells), chart_count, breakpoint_checks


def chromatic_number(good_edges: tuple[tuple[int, int], ...]):
    edge_set = set(good_edges)
    checks = 0
    for colors in range(1, 7):
        for word in product(range(colors), repeat=6):
            checks += 1
            if all(word[a] != word[b] for a, b in edge_set):
                return colors, word, checks
    raise RuntimeError("six-vertex graph not six-colourable")


def is_clique(vertices: tuple[int, ...], good_edges: tuple[tuple[int, int], ...]) -> bool:
    edge_set = set(good_edges)
    return all(tuple(sorted((a, b))) in edge_set for a, b in combinations(vertices, 2))


def tent_support_cell_hit(ruler: int, difference: int, cell: int) -> bool:
    """Whether ``||difference*t||<1/7`` meets this whole ``1/L`` cell."""
    require(ruler % 7 == 0 and difference > 0, ("bad tent-support input", ruler, difference))
    left = difference * cell
    right = difference * (cell + 1)
    radius = ruler // 7
    quotient = left // ruler
    return any(
        max(left, integer * ruler - radius) < min(right, integer * ruler + radius)
        for integer in range(quotient - 1, quotient + 3)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bodies = tuple(combinations(range(1, 15), 6))
    require(len(bodies) == BODY_COUNT, ("body census changed", len(bodies)))
    chi_counts: Counter[int] = Counter()
    graph_rows = []
    exception_rows = []
    weakest_good = None
    pair_count = 0
    direct_checks = 0
    debt_checks = 0
    candidate_count = 0
    chart_count = 0
    breakpoint_checks = 0
    coloring_checks = 0
    profile_digest = hashlib.sha256()
    blind_support_checks = 0
    blind_support_hits = 0

    for body in bodies:
        ruler, ranges = D1.R.safe_cell_ranges(body)
        require(ruler == 14 * lcm(*body), ("body ruler changed", body, ruler))
        debt = universal_debt(body, ruler)
        for e in body:
            for level in range(1, 6):
                term = F(e, 7 * (level * ruler - e))
                ceiling_term = F(e, 7 * (ruler - e))
                require(term <= ceiling_term, ("universal debt term failed", body, e, level))
                debt_checks += 1

        good = []
        bad = []
        for i, k in EDGES:
            profile, candidates, charts, pieces = best_q0_profile(
                body, ruler, ranges, i, k
            )
            pair_count += 1
            candidate_count += candidates
            chart_count += charts
            breakpoint_checks += pieces
            margin = None if profile is None else profile[0] - debt
            if margin is not None and margin > 0:
                good.append((i, k))
                weakest_row = (margin, body, (i, k), profile, debt)
                if weakest_good is None or weakest_row < weakest_good:
                    weakest_good = weakest_row
            else:
                bad.append((i, k))

            if profile is not None:
                floor, cell, chart = profile
                a, b = body[i], body[k]
                require(D1.R.body_cell_is_safe(ruler, body, cell),
                        ("selected chart cell unsafe", body, (i, k), cell))
                for level in (1, 2, 5):
                    actual = D1.intersection_mass(
                        D1.R.reflected_level_arcs(ruler, a, level, cell),
                        D1.R.reflected_level_arcs(ruler, b, level, cell),
                    )
                    formula = D1.overlap_formula(ruler, a, b, level, chart)
                    require(actual == formula and actual >= floor,
                            ("signed chart control failed", body, (i, k), level,
                             actual, formula, floor))
                    direct_checks += 1
            profile_digest.update(
                f"{body}|{ruler}|{(i,k)}|{debt}|{profile}|{margin}|{candidates}|{charts}|{pieces}\n".encode()
            )

        good_tuple = tuple(good)
        bad_tuple = tuple(bad)
        chi, coloring, checks = chromatic_number(good_tuple)
        coloring_checks += checks
        require(chi >= 4, ("three-colourable good graph", body, good_tuple, coloring))
        chi_counts[chi] += 1
        graph_rows.append((body, ruler, debt, good_tuple, bad_tuple, chi, coloring))
        if bad_tuple:
            expected = next(
                (row for row in EXPECTED_EXCEPTIONS if row[0] == body),
                None,
            )
            require(expected is not None, ("unexpected incomplete graph", body, bad_tuple))
            _, expected_bad, clique, expected_coloring = expected
            require(
                bad_tuple == expected_bad
                and chi == 4
                and is_clique(clique, good_tuple)
                and all(expected_coloring[a] != expected_coloring[b] for a, b in good_tuple),
                ("exception graph changed", body, bad_tuple, chi, coloring),
            )
            # These are not weak positive edges.  They are exact blind
            # corridors: the difference tent has empty support on every
            # body-safe cell.  Consequently the corresponding same-level
            # clauses are disjoint for every common level Q by the circle
            # triangle inequality.
            for i, k in bad_tuple:
                difference = body[k] - body[i]
                for left, right in ranges:
                    for cell in range(left, right):
                        hit = tent_support_cell_hit(ruler, difference, cell)
                        blind_support_hits += int(hit)
                        blind_support_checks += 1
            exception_rows.append((body, bad_tuple, clique, expected_coloring))
        else:
            require(len(good_tuple) == 15 and chi == 6,
                    ("complete graph chromatic audit failed", body, chi))

    require(pair_count == PAIR_COUNT, ("pair census changed", pair_count))
    require(
        sum(1 for row in graph_rows if not row[4]) == COMPLETE_GRAPH_COUNT,
        "complete graph count changed",
    )
    require(tuple(sorted(chi_counts.items())) == EXPECTED_CHI_COUNTS,
            ("chromatic histogram changed", chi_counts))
    require(tuple(exception_rows) == EXPECTED_EXCEPTIONS,
            ("exception ledger changed", exception_rows))
    require(weakest_good == EXPECTED_WEAKEST_GOOD,
            ("weakest good margin changed", weakest_good))
    require(
        blind_support_checks == EXPECTED_BLIND_SUPPORT_CHECKS
        and blind_support_hits == 0,
        ("blind-corridor support census changed", blind_support_checks, blind_support_hits),
    )

    semantic_payload = (
        BODY_COUNT,
        PAIR_COUNT,
        EXPECTED_BASE_SHA256,
        sha256(SIGNED_D1),
        tuple(graph_rows),
        tuple(sorted(chi_counts.items())),
        tuple(exception_rows),
        weakest_good,
        pair_count,
        direct_checks,
        debt_checks,
        candidate_count,
        chart_count,
        breakpoint_checks,
        coloring_checks,
        blind_support_checks,
        blind_support_hits,
        profile_digest.hexdigest(),
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic changed", semantic))

    lines = [
        "LRC14 universal repeated-level signed-pair chromatic closure",
        f"base_engine_sha256={sha256(BASE)}",
        f"signed_D1_referee_sha256={sha256(SIGNED_D1)}",
        f"bodies={BODY_COUNT};labelled_pairs={pair_count};Q0=1",
        f"universal_debt=epsilon(q)<=epsilon_max(E)=sum_e e/[7(L-e)];term_checks={debt_checks}",
        f"signed_profiles=candidates={candidate_count};charts={chart_count};breakpoint_piece_checks={breakpoint_checks};direct_formula_floor_checks={direct_checks}",
        f"graph_census=K6:{COMPLETE_GRAPH_COUNT};chi_histogram={tuple(sorted(chi_counts.items()))};coloring_checks={coloring_checks}",
        f"exception_blind_corridors=tent_support_cell_checks={blind_support_checks};hits={blind_support_hits};bad same-level pairs are exactly disjoint for every Q",
        "chromatic_mechanism=good edge means signed same-level overlap floor>epsilon_max;chi(G_E)>=4 forces a monochromatic good edge in every <=3-level assignment",
    ]
    for body, bad_edges, clique, coloring in exception_rows:
        labels_bad = tuple((body[i], body[k]) for i, k in bad_edges)
        labels_clique = tuple(body[i] for i in clique)
        lines.append(
            f"EXCEPTION;E={body};bad_edges_slots={bad_edges};bad_edges_labels={labels_bad};"
            f"chi=4;K4_slots={clique};K4_labels={labels_clique};proper4={coloring}"
        )
    require(weakest_good is not None, "no good edge found")
    margin, body, edge, profile, debt = weakest_good
    lines.append(
        f"weakest_good_margin={qtext(margin)};E={body};pair_slots={edge};"
        f"pair_labels={(body[edge[0]],body[edge[1]])};floor={qtext(profile[0])};"
        f"cell={profile[1]};chart={profile[2]};epsilon_max={qtext(debt)}"
    )
    lines.extend(
        (
            "consequence=every reflected k=1 packet using at most three distinct positive levels closes",
            "corollary=the full normalized D=2 sector q_e in {m,m+1,m+2} is empty for every m>=1",
            "bodywise_reduction=on 3001 K6 bodies any uncertified word has six pairwise-distinct levels;on the two exceptions it is a proper coloring of the displayed chi4 graph",
            "scope_boundary=proper-coloring residues are certificate failures only,not physical survivors",
            f"profile_digest_sha256={profile_digest.hexdigest()}",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
