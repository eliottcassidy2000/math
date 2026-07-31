#!/usr/bin/env python3
r"""Low-phase ratio graph, sharp periodized floor, and robust-body closure.

This exact referee isolates the projective level graph behind THM-1226.  Two
distinct positive levels ``p,q`` are *low-phase adjacent* when, after writing

    p=gP, q=gQ,  gcd(P,Q)=1,

one has ``P+Q<=7``.  Its clique number is exactly five.  After scaling the
least vertex to one, the only five-cliques are

    {1,3/2,2,3,6},       {1,2,3,4,6}.                 (A)

The normalized clique counts in sizes one through five are ``1,8,13,8,2``.
In particular six distinct levels always contain a pair with ``P+Q>=8``.

The positive-phase lower bound can be made coefficient independent.  Put

    T_s(z)=sum_(n in Z) (s-|z+n|)_+.

The exact periodized THM-1226 fiber is

    F_(P,Q)(z)=[T_((P+Q)/14)(z)-T_((Q-P)/14)(z)]/(PQ). (B)

The Fourier series of the periodic tent is

    T_s(z)=s^2+2 sum_(n>=1) [sin(pi*n*s)/(pi*n)]^2 cos(2*pi*n*z).

If ``r={s}``, the absolute Fourier tail is at most
``r(1-r)<=1/4``.  Hence

    F_(P,Q)(z) >= 1/49-1/(2PQ).                         (C)

For ``PQ>=46`` this is strictly above ``1/105``.  For ``PQ<46``, all
breakpoints lie on the fourteen-grid; the exact 63-pair bank below proves

    min_z F_(P,Q)(z) >= 1/105       whenever P+Q>=8,    (D)

with equality only for channel ``(3,5)``.  Notice that the height
``1/(7Q)`` belongs to one *unperiodized* trapezoid and does not upper-bound
the sum (B): for example ``min F_(414,415)=12272/601335``.  The much smaller
``137/400890`` used in the D=3 proof was a valid one-lift lower bound, not
the exact minimum.

Now let ``E`` be a six-body, ``L=14*lcm(E)``, and

    eps_max(E)=sum_(e in E) e/[7(L-e)].

Call a label pair ``a,b`` *robust* when

    1/105 - 4(a+b)/L > eps_max(E).                      (R)

For levels ``p=gP,q=gQ``, subdivide a body-safe cell by
``u=(r+x)/g``.  The exact phases are

    P*x-a(gj+r+x)/(gL),   Q*x-b(gj+r+x)/(gL).

Dropping the two ``x/(gL)`` terms gives one primitive fiber (B), while the
two clause symmetric differences total at most ``4(a+b)/(gL)``.  Thus (D),
``g>=1``, and (R) make every robust high-phase pair beat the exact singleton
debt and close by one-edge Bonferroni.

Exactly 2,217 of the 3,003 bodies have robust graph ``K6``.  Two of these are
the exceptional same-level bodies from the universal chromatic theorem.
On each of the other 2,215 bodies, every arbitrary positive level assignment
closes: a repeat lies on a same-level-good edge, while six distinct levels
have a high-phase pair by (A), and every label pair is robust.

The boundary is exact for this uncoloured pair mechanism.  The integer level
set

    {2,3,4,6,8,12}                                      (H)

induces ``K6`` minus exactly the edge ``3--8`` in the low-phase graph.  Given
any noncomplete robust graph, put ``3,8`` on one of its nonedges and the
other four values on the remaining labels.  Every robust edge is then
low-phase, so the mechanism has no certified pair.  This is a hostile to the
certificate, not a physical survivor.

For the two exceptional robust-K6 bodies, any uncertified word must be a
proper colouring of the displayed same-level-good graph and its distinct
values must form a low-phase clique.  Hence it has four or five distinct
values.  There are eight normalized four-cliques and the two five-cliques
in (A), giving exactly 1,584 labelled scale rays before further certificates:
384 four-colour rays and 1,200 five-colour rays.  This is an exact next finite
obstruction, not a claim that those rays survive.
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation" / "lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
UNIVERSAL = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.py"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_low_phase_clique_robust_body_closure_thm2941.out"
)
EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_UNIVERSAL_SHA256 = "dc6f23a201e817dd9134e8660d35e83d3053c67d26fc271ce3eae07f0f857689"
EXPECTED_SEMANTIC_SHA256 = "0063d2beb7c627f1836576081689720d41f45a078b4b7b3e06aeb81574039248"

LOW_RATIOS = frozenset(
    (F(4, 3), F(3, 2), F(2), F(5, 2), F(3), F(4), F(5), F(6))
)
NORMALIZED_VERTICES = (F(1),) + tuple(sorted(LOW_RATIOS))
EXPECTED_CLIQUE_COUNTS = ((1, 1), (2, 8), (3, 13), (4, 8), (5, 2))
EXPECTED_MAX_CLIQUES = (
    (F(1), F(3, 2), F(2), F(3), F(6)),
    (F(1), F(2), F(3), F(4), F(6)),
)
HOSTILE_SIX = (F(2), F(3), F(4), F(6), F(8), F(12))
HOSTILE_HIGH_EDGE = (F(3), F(8))
FIBER_FLOOR = F(1, 105)
PRODUCT_TAIL_START = 46
FINITE_PAIR_COUNT = 63
BODY_COUNT = 3003
ROBUST_COMPLETE_COUNT = 2217
ARBITRARY_LEVEL_CLOSED_COUNT = 2215
EXPECTED_ROBUST_CLIQUE_HIST = ((1, 211), (2, 144), (3, 159), (4, 152), (5, 120), (6, 2217))
EXPECTED_ROBUST_EDGE_HIST = (
    (0, 211), (1, 88), (2, 41), (3, 53), (4, 53), (5, 36),
    (6, 46), (7, 33), (8, 21), (9, 35), (10, 32), (11, 15),
    (12, 37), (13, 26), (14, 59), (15, 2217),
)
EXCEPTIONS = (
    ((1, 2, 7, 9, 11, 13), ((2, 3), (3, 4), (4, 5))),
    ((2, 4, 7, 9, 11, 13), ((2, 4), (3, 5))),
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def low_adjacent(x: F, y: F) -> bool:
    require(x > 0 and y > 0 and x != y, (x, y))
    return max(x, y) / min(x, y) in LOW_RATIOS


def is_low_clique(vertices: tuple[F, ...]) -> bool:
    return all(low_adjacent(x, y) for x, y in combinations(vertices, 2))


def trapezoid(p: int, q: int, value: F) -> F:
    require(1 <= p <= q, (p, q))
    value = abs(value)
    plateau = F(q - p, 14)
    support = F(p + q, 14)
    if value <= plateau:
        return F(1, 7 * q)
    if value < support:
        return (support - value) / (p * q)
    return F(0)


def direct_periodization(p: int, q: int, z: F) -> F:
    bound = (p + q) // 14 + 3
    return sum((trapezoid(p, q, z + n) for n in range(-bound, bound + 1)), F(0))


def periodic_tent(radius: F, z: F) -> F:
    bound = radius.numerator // radius.denominator + 3
    return sum((max(F(0), radius - abs(z + n)) for n in range(-bound, bound + 1)), F(0))


def fiber_density(p: int, q: int, z: F) -> F:
    require(1 <= p <= q and gcd(p, q) == 1, (p, q))
    return (
        periodic_tent(F(p + q, 14), z)
        - periodic_tent(F(q - p, 14), z)
    ) / (p * q)


def fiber_minimum(p: int, q: int) -> tuple[F, tuple[int, ...]]:
    rows = tuple((fiber_density(p, q, F(k, 14)), k) for k in range(14))
    minimum = min(value for value, _ in rows)
    return minimum, tuple(k for value, k in rows if value == minimum)


def universal_debt(body: tuple[int, ...], ruler: int) -> F:
    return sum((F(e, 7 * (ruler - e)) for e in body), F(0))


def robust_edges(body: tuple[int, ...]) -> tuple[int, F, tuple[tuple[int, int], ...]]:
    ruler = 14 * lcm(*body)
    debt = universal_debt(body, ruler)
    edges = tuple(
        (i, j)
        for i, j in combinations(range(6), 2)
        if FIBER_FLOOR - F(4 * (body[i] + body[j]), ruler) > debt
    )
    return ruler, debt, edges


def clique_number(edges: tuple[tuple[int, int], ...]) -> int:
    edge_set = set(edges)
    return max(
        len(vertices)
        for size in range(1, 7)
        for vertices in combinations(range(6), size)
        if all(tuple(sorted(edge)) in edge_set for edge in combinations(vertices, 2))
    )


def proper_surjections(bad_edges: tuple[tuple[int, int], ...], colors: int):
    good_edges = set(combinations(range(6), 2)) - set(bad_edges)
    return tuple(
        word
        for word in product(range(colors), repeat=6)
        if len(set(word)) == colors and all(word[i] != word[j] for i, j in good_edges)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    require(sha256(BASE) == EXPECTED_BASE_SHA256, "reflected interval engine changed")
    require(sha256(UNIVERSAL) == EXPECTED_UNIVERSAL_SHA256, "universal chromatic theorem changed")

    normalized_cliques = {
        size: tuple(
            row
            for row in combinations(NORMALIZED_VERTICES, size)
            if row[0] == 1 and is_low_clique(row)
        )
        for size in range(1, 7)
    }
    clique_counts = tuple((size, len(normalized_cliques[size])) for size in range(1, 6))
    require(clique_counts == EXPECTED_CLIQUE_COUNTS, (clique_counts, normalized_cliques))
    require(not normalized_cliques[6], normalized_cliques[6])
    require(normalized_cliques[5] == EXPECTED_MAX_CLIQUES, normalized_cliques[5])

    hostile_edges = tuple(
        pair for pair in combinations(HOSTILE_SIX, 2) if low_adjacent(*pair)
    )
    hostile_missing = tuple(pair for pair in combinations(HOSTILE_SIX, 2) if pair not in hostile_edges)
    require(hostile_missing == (HOSTILE_HIGH_EDGE,), hostile_missing)

    tail_margin = F(1, 49) - F(1, 2 * PRODUCT_TAIL_START) - FIBER_FLOOR
    require(tail_margin == F(1, 67620) > 0, tail_margin)
    for residue in range(14):
        r = F(residue, 14)
        require(r * (1 - r) <= F(1, 4), (residue, r))

    finite_rows = []
    route_checks = 0
    for q in range(1, PRODUCT_TAIL_START):
        for p in range(1, q + 1):
            if gcd(p, q) != 1 or p + q < 8 or p * q >= PRODUCT_TAIL_START:
                continue
            minimum, minimizers = fiber_minimum(p, q)
            require(minimum >= FIBER_FLOOR, (p, q, minimum, minimizers))
            for k in range(14):
                z = F(k, 14)
                require(fiber_density(p, q, z) == direct_periodization(p, q, z),
                        (p, q, z))
                route_checks += 1
            finite_rows.append((minimum, p, q, minimizers))
    require(len(finite_rows) == FINITE_PAIR_COUNT, len(finite_rows))
    require(min(finite_rows) == (FIBER_FLOOR, 3, 5, (6, 7, 8)), min(finite_rows))
    equality_channels = tuple((p, q) for minimum, p, q, _ in finite_rows if minimum == FIBER_FLOOR)
    require(equality_channels == ((3, 5),), equality_channels)
    hostile_large = fiber_minimum(414, 415)
    require(hostile_large == (F(12272, 601335), (3, 4, 5, 6, 7, 8, 9, 10, 11)), hostile_large)
    for k in range(14):
        z = F(k, 14)
        require(fiber_density(414, 415, z) == direct_periodization(414, 415, z),
                (414, 415, z))
        route_checks += 1
    crude_large = min(F(1, 7 * 415), F(414 + 415 - 7, 14 * 414 * 415))
    require(crude_large == F(137, 400890) < FIBER_FLOOR < hostile_large[0], crude_large)

    bodies = tuple(combinations(range(1, 15), 6))
    require(len(bodies) == BODY_COUNT, len(bodies))
    edge_hist: Counter[int] = Counter()
    clique_hist: Counter[int] = Counter()
    complete_rows = []
    incomplete_hostile_checks = 0
    body_digest = hashlib.sha256()
    for body in bodies:
        ruler, debt, edges = robust_edges(body)
        edge_hist[len(edges)] += 1
        omega = clique_number(edges)
        clique_hist[omega] += 1
        if len(edges) == 15:
            complete_rows.append((body, ruler, debt))
        else:
            edge_set = set(edges)
            nonedge = next(edge for edge in combinations(range(6), 2) if edge not in edge_set)
            assignment = [None] * 6
            assignment[nonedge[0]], assignment[nonedge[1]] = HOSTILE_HIGH_EDGE
            remaining_values = iter(value for value in HOSTILE_SIX if value not in HOSTILE_HIGH_EDGE)
            for index in range(6):
                if assignment[index] is None:
                    assignment[index] = next(remaining_values)
            require(
                all(low_adjacent(assignment[i], assignment[j]) for i, j in edges),
                (body, edges, nonedge, assignment),
            )
            incomplete_hostile_checks += 1
        body_digest.update(f"{body}|{ruler}|{debt}|{edges}|{omega}\n".encode())
    require(tuple(sorted(edge_hist.items())) == EXPECTED_ROBUST_EDGE_HIST, edge_hist)
    require(tuple(sorted(clique_hist.items())) == EXPECTED_ROBUST_CLIQUE_HIST, clique_hist)
    require(len(complete_rows) == ROBUST_COMPLETE_COUNT, len(complete_rows))
    require(all(any(body == exception for body, _, _ in complete_rows) for exception, _ in EXCEPTIONS),
            "exception missing from robust-K6 bank")
    require(len(complete_rows) - len(EXCEPTIONS) == ARBITRARY_LEVEL_CLOSED_COUNT,
            len(complete_rows))
    require(incomplete_hostile_checks == BODY_COUNT - ROBUST_COMPLETE_COUNT,
            incomplete_hostile_checks)

    exception_rows = []
    four_types = len(normalized_cliques[4])
    five_types = len(normalized_cliques[5])
    for body, bad_edges in EXCEPTIONS:
        four_words = proper_surjections(bad_edges, 4)
        five_words = proper_surjections(bad_edges, 5)
        require(len(four_words) == 24, (body, len(four_words)))
        if body == EXCEPTIONS[0][0]:
            require(len(five_words) == 360, (body, len(five_words)))
        else:
            require(len(five_words) == 240, (body, len(five_words)))
        exception_rows.append(
            (body, bad_edges, len(four_words), len(five_words),
             four_types * len(four_words), five_types * len(five_words))
        )
    four_rays = sum(row[-2] for row in exception_rows)
    five_rays = sum(row[-1] for row in exception_rows)
    require((four_rays, five_rays, four_rays + five_rays) == (384, 1200, 1584),
            (four_rays, five_rays))

    semantic_payload = (
        clique_counts,
        normalized_cliques[4],
        normalized_cliques[5],
        hostile_missing,
        tail_margin,
        tuple(sorted(finite_rows)),
        hostile_large,
        crude_large,
        tuple(sorted(edge_hist.items())),
        tuple(sorted(clique_hist.items())),
        tuple(complete_rows),
        incomplete_hostile_checks,
        tuple(exception_rows),
        four_rays,
        five_rays,
        body_digest.hexdigest(),
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, (semantic, EXPECTED_SEMANTIC_SHA256))

    source_sha = sha256(Path(__file__))
    lines = [
        "LRC14 reflected low-phase clique and robust-body closure exact referee",
        f"low_ratio_graph=omega_5;normalized_clique_counts={clique_counts};"
        f"maximum_cliques={EXPECTED_MAX_CLIQUES}",
        f"normalized_four_cliques={normalized_cliques[4]}",
        f"periodized_fiber_floor={qtext(FIBER_FLOOR)};tail_product_start={PRODUCT_TAIL_START};"
        f"tail_boundary_margin={qtext(tail_margin)};finite_pairs={len(finite_rows)};"
        f"finite_route_checks={route_checks};equality_channel={(3,5)}",
        f"hostile_414_415_exact_min={qtext(hostile_large[0])};minimizers={hostile_large[1]};"
        f"one_lift_lower={qtext(crude_large)}",
        f"robust_complete_bodies={len(complete_rows)}/{BODY_COUNT};"
        f"arbitrary_level_closed_nonexceptional={ARBITRARY_LEVEL_CLOSED_COUNT};"
        f"robust_edge_hist={tuple(sorted(edge_hist.items()))};"
        f"robust_clique_hist={tuple(sorted(clique_hist.items()))}",
        f"sharp_certificate_hostile_levels={HOSTILE_SIX};unique_high_pair={HOSTILE_HIGH_EDGE};"
        f"incomplete_body_relabellings_checked={incomplete_hostile_checks}",
        f"exception_residuals={tuple(exception_rows)};four_colour_scale_rays={four_rays};"
        f"five_colour_scale_rays={five_rays};total_scale_rays={four_rays+five_rays}",
        "scope=closes arbitrary reflected k=1 levels on 2215 bodies; the 786 noncomplete-robust bodies and 1584 exceptional low-clique rays are certificate residues, not physical survivors",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_sha256={sha256(BASE)}",
        f"universal_sha256={sha256(UNIVERSAL)}",
        f"body_digest={body_digest.hexdigest()}",
        f"source_sha256={source_sha}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    output = "\n".join(lines) + "\n"
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
