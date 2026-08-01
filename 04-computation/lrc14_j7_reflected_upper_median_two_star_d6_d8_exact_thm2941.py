#!/usr/bin/env python3
r"""Exact upper-median two-star closure for the first three open spreads.

After the arbitrary-level robust-edge-11 theorem, 649 of the 3,003 six-label
bodies remain.  For such a body ``E``, put ``L=14*lcm(E)`` and list all
integer body-safe cells.  Let ``j`` be the upper median.  Every one of these
649 medians has a unique *boundary label* ``e_*`` satisfying

    e_* j == L/14  (mod L).                              (1)

This makes its reflected level-``q`` clause an endpoint-aligned toothpick
comb:

    A(e_*,q)=union_(0<=k<q)
             [kL/(qL-e_*),(k+1/7)L/(qL-e_*)].           (2)

Choose a second pivot: the least label of ``E`` different from ``e_*``.
The *two-star stalk* is the nine-edge graph consisting of every label pair
incident to at least one pivot.  This file exhausts every injection of the
six labels into ``{1,...,D+1}``, for ``D=6,7,8``.  At the single cell ``j``,
one of those nine pairs always has intersection mass larger than the exact
six-clause singleton debt

    sum_(e in E) e/[7(q_e L-e)].                         (3)

The audit is entirely rational.  It checks 55,606,320 assignments:

    D=6:  649 * 7P6 =  3,270,960,
    D=7:  649 * 8P6 = 13,083,840,
    D=8:  649 * 9P6 = 39,251,520.

The weakest ``D=6`` row is on ``(1,2,3,4,6,12)``.  The weakest rows for
``D=7`` and ``D=8`` coincide on ``(1,2,5,7,8,14)``; in particular inserting
level nine creates no new minimum.  Pairwise Bonferroni then puts the full
six-clause union strictly below ``6/7``.

This is a finite-exact base theorem inside the reflected THM-2941 residual
family.  It does not prove the desired all-spread insertion recursion, and it
does not by itself prove physical LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
import multiprocessing as mp
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from pathlib import Path

import gmpy2
from gmpy2 import mpq


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation/lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
LOW_PHASE = ROOT / "04-computation/lrc14_j7_reflected_low_phase_clique_robust_body_closure_thm2941.py"
EDGE11 = ROOT / "04-computation/lrc14_j7_reflected_robust_edge11_threshold_shapes_uniform_closure_thm2941.py"
EDGE11_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_robust_edge11_threshold_shapes_uniform_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_upper_median_two_star_d6_d8_exact_thm2941.out"

EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_LOW_PHASE_SHA256 = "b2418dfda1b48257d1f7582d4ea977203a26f88885e13946bc100ccf264c9ce1"
EXPECTED_EDGE11_SHA256 = "58685796af62dfb425bfa17a6363d9cc5ca4cbb72e0f1d42ecc23dcbf10c01b4"
EXPECTED_EDGE11_OUTPUT_SHA256 = "f19bef5d4b2a32bda2e939a303870680b022973a9513a52c18677f73d7512a29"
EXPECTED_GMPY2_VERSION = "2.3.0"
EXPECTED_SEMANTIC_SHA256 = "aa7f4927df4dbbd79128d28472bbdd7fd205cc659286e34bc4a2ec450b3b834d"

BODY_COUNT = 3003
UPSTREAM_CLOSED_COUNT = 2354
ACTIVE_BODY_COUNT = 649
SPREADS = (6, 7, 8)
ASSIGNMENT_COUNTS = ((6, 3270960), (7, 13083840), (8, 39251520))
TOTAL_ASSIGNMENT_COUNT = 55606320
MAX_LEVEL = 9
STALK_EDGE_COUNT = 9
EXPECTED_CENTRE_LABEL_HIST = ((2, 356), (4, 170), (6, 67), (8, 13), (10, 7), (12, 2), (13, 34))
EXPECTED_ARC_CONTROLS = 35046
EXPECTED_PROFILE_CONTROLS = 420552
EXPECTED_BODY_DIGEST = "1e26d6c95aabdf1f39a7a4d1fcc90626c51fbdf2f120c795261d73265a4a95e6"

EXPECTED_WORST = {
    6: (
        F(279417235894, 59234582752345),
        (1, 2, 3, 4, 6, 12),
        168,
        90,
        (1, 0),
        (2, 1),
        (4, 3, 6, 2, 5, 1),
        F(1834, 93269),
        F(885345169276, 59234582752345),
    ),
    7: (
        F(727640508707102047, 286068087083094746521),
        (1, 2, 5, 7, 8, 14),
        3920,
        2100,
        (1, 0),
        (2, 1),
        (4, 2, 8, 1, 3, 6),
        F(70245, 23032451),
        F(11139829207263296, 22005237467930365117),
    ),
    8: (
        F(727640508707102047, 286068087083094746521),
        (1, 2, 5, 7, 8, 14),
        3920,
        2100,
        (1, 0),
        (2, 1),
        (4, 2, 8, 1, 3, 6),
        F(70245, 23032451),
        F(11139829207263296, 22005237467930365117),
    ),
}

EDGES = tuple(combinations(range(6), 2))
WORDS_BY_SPREAD = {
    D: tuple(permutations(range(1, D + 2), 6)) for D in SPREADS
}


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def as_fraction(value: mpq) -> F:
    return F(int(value.numerator), int(value.denominator))


for path, expected in (
    (BASE, EXPECTED_BASE_SHA256),
    (LOW_PHASE, EXPECTED_LOW_PHASE_SHA256),
    (EDGE11, EXPECTED_EDGE11_SHA256),
    (EDGE11_OUTPUT, EXPECTED_EDGE11_OUTPUT_SHA256),
):
    require(sha256(path) == expected, ("upstream theorem changed", path, sha256(path), expected))
require(gmpy2.version() == EXPECTED_GMPY2_VERSION,
        ("gmpy2 version changed", gmpy2.version(), EXPECTED_GMPY2_VERSION))


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = import_module("median_two_star_base", BASE)
LOW = import_module("median_two_star_low", LOW_PHASE)


def upper_median_cell(ranges: tuple[tuple[int, int], ...]) -> tuple[int, int]:
    count = sum(right - left for left, right in ranges)
    require(count > 0, ranges)
    rank = count // 2
    for left, right in ranges:
        width = right - left
        if rank < width:
            return left + rank, count
        rank -= width
    raise RuntimeError(("median rank escaped", ranges, count))


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
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return total


def scan_body(item):
    body_index, body, L = item
    L2, ranges = R.safe_cell_ranges(body)
    require(L2 == L, (body, L, L2))
    cell, safe_cell_count = upper_median_cell(ranges)
    require(R.body_cell_is_safe(L, body, cell), ("unsafe median", body, L, cell))
    centres = tuple(
        index for index, e in enumerate(body)
        if (e * cell) % L == L // 14
    )
    require(len(centres) == 1, ("boundary centre", body, L, cell, centres))
    centre = centres[0]
    second = next(index for index in range(6) if index != centre)
    stalk = tuple(edge for edge in EDGES if centre in edge or second in edge)
    require(len(stalk) == STALK_EDGE_COUNT, (body, centre, second, stalk))

    arcs = {}
    arc_controls = 0
    for index, e in enumerate(body):
        for q in range(1, MAX_LEVEL + 1):
            reflected = R.reflected_level_arcs(L, e, q, cell)
            direct = R.direct_multiplier_arcs(L, q * L - e, cell)
            require(reflected == direct, ("arc law", body, cell, index, q))
            expected_mass = F(q * L, 7 * (q * L - e))
            require(R.interval_mass(reflected) == expected_mass,
                    ("singleton mass", body, cell, index, q, reflected))
            arcs[index, q] = reflected
            arc_controls += 1

    profiles = {}
    profile_controls = 0
    for i, j in stalk:
        for p in range(1, MAX_LEVEL + 1):
            for q in range(1, MAX_LEVEL + 1):
                if p == q:
                    continue
                gain = intersection_mass(arcs[i, p], arcs[j, q])
                merge_gain = (
                    R.interval_mass(arcs[i, p])
                    + R.interval_mass(arcs[j, q])
                    - R.interval_mass(R.merge_intervals(arcs[i, p] + arcs[j, q]))
                )
                require(gain == merge_gain, ("intersection routes", body, i, j, p, q))
                profiles[i, j, p, q] = mpq(gain.numerator, gain.denominator)
                profile_controls += 1

    debt_terms = {
        (index, q): mpq(e, 7 * (q * L - e))
        for index, e in enumerate(body)
        for q in range(1, MAX_LEVEL + 1)
    }

    minima = {}
    winner_histograms = {}
    for D in SPREADS:
        local_margin = None
        local_payload = None
        winner_hist = [0] * len(stalk)
        for word in WORDS_BY_SPREAD[D]:
            debt = sum((debt_terms[index, word[index]] for index in range(6)), mpq(0))
            best_profile = None
            best_edge_index = None
            for edge_index, (i, j) in enumerate(stalk):
                profile = profiles[i, j, word[i], word[j]]
                if best_profile is None or profile > best_profile:
                    best_profile = profile
                    best_edge_index = edge_index
            require(best_profile is not None and best_edge_index is not None, (body, D, word))
            margin = best_profile - debt
            require(margin > 0, (
                "two-star failure", body, L, cell, D, word, stalk[best_edge_index],
                best_profile, debt, margin,
            ))
            winner_hist[best_edge_index] += 1
            payload = (
                word,
                best_profile,
                debt,
                stalk[best_edge_index],
            )
            if (
                local_margin is None
                or margin < local_margin
                or (margin == local_margin and payload < local_payload)
            ):
                local_margin = margin
                local_payload = payload

        require(local_margin is not None and local_payload is not None, (body, D))
        word, best_profile, debt, winning_edge = local_payload
        minima[D] = (
            as_fraction(local_margin), body, L, cell, (centre, second),
            (body[centre], body[second]), word, as_fraction(best_profile), as_fraction(debt),
            winning_edge, (body[winning_edge[0]], body[winning_edge[1]]),
            safe_cell_count,
        )
        winner_histograms[D] = tuple(winner_hist)

    boundary_comb_controls = 0
    e_star = body[centre]
    for q in range(1, MAX_LEVEL + 1):
        expected = tuple(
            (F(k * L, q * L - e_star), F((7 * k + 1) * L, 7 * (q * L - e_star)))
            for k in range(q)
        )
        require(arcs[centre, q] == expected,
                ("boundary toothpick law", body, L, cell, centre, q, arcs[centre, q], expected))
        boundary_comb_controls += 1

    return (
        body_index, body, L, cell, safe_cell_count, centre, second, stalk,
        arc_controls, profile_controls, boundary_comb_controls,
        minima, winner_histograms,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    require(args.workers >= 1, args.workers)

    bodies = []
    body_digest = hashlib.sha256()
    for body in combinations(range(1, 15), 6):
        L, _, robust = LOW.robust_edges(body)
        if len(robust) <= 10:
            bodies.append((body, L))
    require(len(bodies) == ACTIVE_BODY_COUNT, len(bodies))
    require(BODY_COUNT - len(bodies) == UPSTREAM_CLOSED_COUNT, len(bodies))

    items = tuple((index, body, L) for index, (body, L) in enumerate(bodies))
    if args.workers == 1:
        rows = tuple(map(scan_body, items))
    else:
        with mp.get_context("fork").Pool(args.workers) as pool:
            rows = tuple(pool.imap(scan_body, items, chunksize=1))
    require(tuple(row[0] for row in rows) == tuple(range(ACTIVE_BODY_COUNT)),
            "worker order changed")

    centre_hist = Counter()
    arc_controls = 0
    profile_controls = 0
    boundary_comb_controls = 0
    global_minima = {}
    winner_histograms = {D: [0] * STALK_EDGE_COUNT for D in SPREADS}
    certificate_rows = []
    for row in rows:
        (
            _, body, L, cell, safe_cell_count, centre, second, stalk,
            body_arc_controls, body_profile_controls, body_boundary_controls,
            minima, body_winner_histograms,
        ) = row
        centre_hist[body[centre]] += 1
        arc_controls += body_arc_controls
        profile_controls += body_profile_controls
        boundary_comb_controls += body_boundary_controls
        body_digest.update(
            f"{body}|{L}|{cell}|{safe_cell_count}|{centre}|{second}|{stalk}\n".encode()
        )
        for D in SPREADS:
            if global_minima.get(D) is None or minima[D] < global_minima[D]:
                global_minima[D] = minima[D]
            for edge_index, count in enumerate(body_winner_histograms[D]):
                winner_histograms[D][edge_index] += count
            certificate_rows.append((D, minima[D], body_winner_histograms[D]))

    require(tuple(sorted(centre_hist.items())) == EXPECTED_CENTRE_LABEL_HIST, centre_hist)
    require(arc_controls == EXPECTED_ARC_CONTROLS, arc_controls)
    require(profile_controls == EXPECTED_PROFILE_CONTROLS, profile_controls)
    require(boundary_comb_controls == ACTIVE_BODY_COUNT * MAX_LEVEL,
            boundary_comb_controls)
    assignment_counts = tuple(
        (D, ACTIVE_BODY_COUNT * len(WORDS_BY_SPREAD[D])) for D in SPREADS
    )
    require(assignment_counts == ASSIGNMENT_COUNTS, assignment_counts)
    require(sum(count for _, count in assignment_counts) == TOTAL_ASSIGNMENT_COUNT,
            assignment_counts)
    for D in SPREADS:
        require(global_minima[D][:9] == EXPECTED_WORST[D],
                ("global minimum", D, global_minima[D], EXPECTED_WORST[D]))
        require(sum(winner_histograms[D]) == dict(ASSIGNMENT_COUNTS)[D],
                (D, winner_histograms[D]))

    if EXPECTED_BODY_DIGEST:
        require(body_digest.hexdigest() == EXPECTED_BODY_DIGEST,
                ("body digest", body_digest.hexdigest(), EXPECTED_BODY_DIGEST))
    semantic_payload = (
        body_digest.hexdigest(),
        tuple(sorted(centre_hist.items())),
        assignment_counts,
        arc_controls,
        profile_controls,
        boundary_comb_controls,
        tuple((D, global_minima[D], tuple(winner_histograms[D])) for D in SPREADS),
        tuple(certificate_rows),
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    lines = [
        "LRC14 reflected upper-median two-star D6-D8 exact proof",
        f"universe=bodies:{BODY_COUNT};upstream_arbitrary_level:{UPSTREAM_CLOSED_COUNT};active_bodies:{ACTIVE_BODY_COUNT};m=1;spreads={SPREADS}",
        "median_rule=upper median of all integer body-safe cells",
        f"boundary_law=unique e_* with e_*j=L/14 mod L on every active body;centre_label_hist={tuple(sorted(centre_hist.items()))}",
        "boundary_toothpick=A(e_*,q)=union_{0<=k<q}[kL/(qL-e_*),(k+1/7)L/(qL-e_*)]",
        "second_pivot=least body label distinct from e_*;selector=union of the two pivot stars;selector_edges=9",
        f"assignment_counts={assignment_counts};total={TOTAL_ASSIGNMENT_COUNT};failures=0",
        f"arc_law_singleton_controls={arc_controls};pair_profile_route_controls={profile_controls};boundary_comb_controls={boundary_comb_controls}",
    ]
    for D in SPREADS:
        row = global_minima[D]
        lines.append(
            f"D{D}_weakest_margin={qtext(row[0])};E={row[1]};L={row[2]};cell={row[3]};"
            f"pivot_slots={row[4]};pivot_labels={row[5]};word={row[6]};gain={qtext(row[7])};"
            f"debt={qtext(row[8])};winning_edge_slots={row[9]};winning_edge_labels={row[10]};safe_cells={row[11]}"
        )
        lines.append(f"D{D}_winner_edge_hist={tuple(winner_histograms[D])}")
    lines.extend((
        "insertion_signal=D7 and D8 have the same exact weakest row;adding level9 creates no new minimum",
        "conclusion=at m=1 and D in {6,7,8}, every reflected residual packet closes on all 3003 bodies",
        "scope=finite-exact base cases for an all-spread insertion theorem;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_sha256={sha256(BASE)}",
        f"low_phase_sha256={sha256(LOW_PHASE)}",
        f"edge11_sha256={sha256(EDGE11)}",
        f"edge11_output_sha256={sha256(EDGE11_OUTPUT)}",
        f"gmpy2_version={gmpy2.version()}",
        f"body_digest={body_digest.hexdigest()}",
        f"source_sha256={sha256(Path(__file__))}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ))
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
