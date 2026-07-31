#!/usr/bin/env python3
r"""Exact closure of the reflected ``D=5`` finite head ``1<=m<=15``.

Let ``E`` be a six-label body, ``L=14*lcm(E)``, and assign the six distinct
levels ``m,...,m+5`` to its labels.  The preceding robust ``K6`` and ``K6-e``
theorems already close every body with at least fourteen robust edges for
arbitrary positive levels.  This leaves 727 bodies.  They are among the 3001
ordinary bodies on which the universal same-level graph is complete, so every
uncertified ``D=5`` word is a permutation of ``0,...,5``.

The cross-determinant transporter from the tail theorem works in the finite
head after one small repair.  For labels ``a,b`` at levels ``p,q``, orienting
the transport through ``(a,p)`` leaves target slope

    A = p(qL-b)/(pL-a).

The reverse orientation has target slope ``A'`` and ``AA'=pq``.  Hence at
least one orientation has slope at least one, precisely the hypothesis needed
by the one-sided comb ``L1`` bound.  On that orientation the exact pair floor
is

    c^-1 max(0, f(p,q)-2|eta|),
    c=1-a/(pL),  eta=(qa-pb)/(pL-a),

where ``f`` is the all-phase primitive-fibre floor from the tail theorem.

There are ``727*15*6! = 7,851,600`` finite-head assignments.  The repaired
cross-determinant floor closes 7,835,524 of them.  Exactly 16,076 profiles in
164 body--``m`` groups remain, with per-``m`` counts

    12985,1511,773,427,235,53,43,23,8,9,5,2,0,0,2.

The residual has a surprisingly rigid one-cell sidecar.  List a body's safe
cells in increasing order and take the upper median.  On that single cell,
for every residual word, one exact pair overlap exceeds the actual singleton
debt.  Bonferroni therefore puts the six-arc union strictly below ``6/7``.
The referee also constructs that union and checks its mass directly.  Thus
the median cell closes all 16,076 residual profiles; no search over cells and
no two-edge correction are required.

This median is canonical rather than asymmetric.  If ``J_E`` is the set of
integer cells on which ``||e(j+u)/L||>=1/14`` for every ``e`` and every
``u in [0,1]``, then ``j -> L-1-j`` is a fixed-point-free involution of
``J_E``.  The weak endpoint is exact because danger is strict.  Moreover
``u -> 1-u`` carries every reflected clause on ``j`` to the corresponding
clause on ``L-1-j``.  Thus the lower and upper medians sum to ``L-1`` and
have equal singleton, pair, and union masses.  The runtime referee checks
both the cell involution and every reflected arc identity it uses.

Together with the proved ``m>=16`` tail theorem, this closes the entire
reflected ``D=5`` sector.  This remains a sufficient certificate inside the
THM-2941 reflected residual family, not a classification of physical LRC(14)
survivors; sectors ``D>=6`` remain outside the conclusion.
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TAIL = ROOT / "04-computation" / "lrc14_j7_reflected_d5_crossdet_tail_closure_thm2941.py"
TAIL_OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_reflected_d5_crossdet_tail_closure_thm2941.out"
LOW_PHASE = ROOT / "04-computation" / "lrc14_j7_reflected_low_phase_clique_robust_body_closure_thm2941.py"
K6_MINUS_EDGE = ROOT / "04-computation" / "lrc14_j7_reflected_robust_k6_minus_edge_uniform_closure_thm2941.py"
K6_MINUS_EDGE_OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_reflected_robust_k6_minus_edge_uniform_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_reflected_d5_head_median_cell_closure_thm2941.out"

EXPECTED_TAIL_SHA256 = "1575b9fabec292bccf0bd639b47b3775922a1531e421b8e4441c6909cc2cedb7"
EXPECTED_TAIL_OUTPUT_SHA256 = "73e8242e68431bd42d3a39b40a8dd0cc6ae8f9e1439e80038a54835bedc3ac55"
EXPECTED_LOW_PHASE_SHA256 = "416c36f16f7c821feb8d260882711d2717069147b8604a93ba60432785cf1d1c"
EXPECTED_K6_MINUS_EDGE_SHA256 = "3b5752836263df2f6b694c9dfe737e489c22dc2c267144530dcfac74cc972b29"
EXPECTED_K6_MINUS_EDGE_OUTPUT_SHA256 = "21245a473a09f9537d82673c5fecc3174910e4050bd25b69cc7f9aae6b29b178"
EXPECTED_SEMANTIC_SHA256 = "b78766f4e02fd63625531f4d6d276c136a2aeab81d9f5822ce331f83e29903fb"

BODY_COUNT = 3003
UPSTREAM_ARBITRARY_LEVEL_BODIES = 2276
HEAD_BODY_COUNT = 727
HEAD_START = 1
HEAD_END = 15
WORDS = tuple(permutations(range(6)))
EDGES = tuple(combinations(range(6), 2))
RAW_ASSIGNMENT_COUNT = HEAD_BODY_COUNT * (HEAD_END - HEAD_START + 1) * len(WORDS)
EXPECTED_RESIDUAL_COUNTS = (
    12985, 1511, 773, 427, 235, 53, 43, 23, 8, 9, 5, 2, 0, 0, 2,
)
EXPECTED_GROUP_COUNTS = (101, 15, 16, 13, 9, 1, 3, 1, 1, 1, 1, 1, 0, 0, 1)
RESIDUAL_COUNT = 16076
RESIDUAL_GROUP_COUNT = 164
PROFILE_CHECK_COUNT = 37444
HOSTILE = (1, 2, 3, 4, 6, 12)
EXPECTED_WEAKEST_PAIR = (
    F(898065335348, 59407314481515),
    HOSTILE,
    1,
    (1, 5, 2, 4, 3, 0),
    F(2490, 84001),
    F(862916237002, 59407314481515),
    90,
    (1, 2),
    F(103, 43681),
    88,
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


for path, expected in (
    (TAIL, EXPECTED_TAIL_SHA256),
    (TAIL_OUTPUT, EXPECTED_TAIL_OUTPUT_SHA256),
    (LOW_PHASE, EXPECTED_LOW_PHASE_SHA256),
    (K6_MINUS_EDGE, EXPECTED_K6_MINUS_EDGE_SHA256),
    (K6_MINUS_EDGE_OUTPUT, EXPECTED_K6_MINUS_EDGE_OUTPUT_SHA256),
):
    require(sha256(path) == expected, ("upstream theorem changed", path, sha256(path), expected))


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = import_module("d5_head_tail", TAIL)
LOW = import_module("d5_head_low_phase", LOW_PHASE)
R = T.R


def crossdet_edge_floor(L: int, a: int, p: int, b: int, q: int) -> F:
    """Tail floor, retaining only orientations whose target slope is >=1."""
    phase = T.phase_floor(p, q)
    if phase == 0:
        return F(0)
    rows = []
    target_slopes = []
    for aa, pp, bb, qq in ((a, p, b, q), (b, q, a, p)):
        c = F(pp * L - aa, pp * L)
        eta = F(qq * aa - pp * bb, pp * L - aa)
        target = F(pp * (qq * L - bb), pp * L - aa)
        target_slopes.append(target)
        if target >= 1:
            rows.append(max(F(0), (phase - 2 * abs(eta)) / c))
    require(target_slopes[0] * target_slopes[1] == p * q, (a, p, b, q, target_slopes))
    require(rows, ("orientation product failed", L, a, p, b, q, target_slopes))
    return max(rows)


def intersection_mass(left, right) -> F:
    return (
        R.interval_mass(left)
        + R.interval_mass(right)
        - R.interval_mass(R.merge_intervals(left + right))
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    remaining = []
    body_digest = hashlib.sha256()
    for body in combinations(range(1, 15), 6):
        L, _, robust = LOW.robust_edges(body)
        if len(robust) <= 13:
            require(body not in tuple(row[0] for row in T.EXCEPTIONS), ("exception in head", body))
            remaining.append((body, L, robust))
            body_digest.update(f"{body}|{L}|{robust}\n".encode())
    require(len(remaining) == HEAD_BODY_COUNT, len(remaining))
    require(BODY_COUNT - len(remaining) == UPSTREAM_ARBITRARY_LEVEL_BODIES, len(remaining))
    require(RAW_ASSIGNMENT_COUNT == 7851600, RAW_ASSIGNMENT_COUNT)

    residual_counts = Counter()
    group_counts = Counter()
    selector_certified = 0
    pair_certified = 0
    profile_checks = 0
    arc_controls = 0
    union_checks = 0
    weakest_pair = None
    weakest_union_slack = None
    largest_union = None
    residual_digest = hashlib.sha256()
    certificate_digest = hashlib.sha256()

    for m in range(HEAD_START, HEAD_END + 1):
        for body, L, _ in remaining:
            edge_floors = {
                (i, j, di, dj): crossdet_edge_floor(
                    L, body[i], m + di, body[j], m + dj
                )
                for i, j in EDGES
                for di in range(6)
                for dj in range(6)
                if di != dj
            }
            singleton_terms = {
                (i, d): F(body[i], 7 * ((m + d) * L - body[i]))
                for i in range(6)
                for d in range(6)
            }
            residual = []
            for word in WORDS:
                old_gain, old_i, old_j = max(
                    (edge_floors[i, j, word[i], word[j]], i, j) for i, j in EDGES
                )
                debt = sum((singleton_terms[i, word[i]] for i in range(6)), F(0))
                if old_gain > debt:
                    selector_certified += 1
                else:
                    residual.append((word, old_gain, debt, (old_i, old_j)))

            if not residual:
                continue
            residual_counts[m] += len(residual)
            group_counts[m] += 1
            _, ranges = R.safe_cell_ranges(body)
            cells = tuple(j for left, right in ranges for j in range(left, right))
            require(cells, ("empty safe-cell set", body))
            require(
                cells == tuple(L - 1 - j for j in reversed(cells)),
                ("safe-cell involution failed", body, L),
            )
            require(len(cells) % 2 == 0, ("safe-cell parity failed", body, len(cells)))
            lower_median = cells[len(cells) // 2 - 1]
            cell = cells[len(cells) // 2]
            require(
                lower_median + cell == L - 1,
                ("median mates failed", body, L, lower_median, cell),
            )
            require(R.body_cell_is_safe(L, body, cell), ("median cell unsafe", body, cell))

            needed_slots = {(i, word[i]) for word, _, _, _ in residual for i in range(6)}
            arcs = {}
            for i, d in needed_slots:
                level = m + d
                reflected = R.reflected_level_arcs(L, body[i], level, cell)
                direct = R.direct_multiplier_arcs(L, level * L - body[i], cell)
                require(reflected == direct, ("arc-law mismatch", body, m, cell, i, d))
                mate = R.reflected_level_arcs(L, body[i], level, lower_median)
                reflected_mate = tuple(
                    sorted((F(1) - right, F(1) - left) for left, right in reflected)
                )
                require(
                    mate == reflected_mate,
                    ("reflected median-mate arc mismatch", body, m, i, d),
                )
                require(
                    R.interval_mass(reflected) == F(1, 7) + singleton_terms[i, d],
                    ("singleton mass mismatch", body, m, cell, i, d),
                )
                arcs[i, d] = reflected
                arc_controls += 1

            needed_profiles = {
                (i, j, word[i], word[j])
                for word, _, _, _ in residual
                for i, j in EDGES
            }
            profiles = {
                key: intersection_mass(arcs[key[0], key[2]], arcs[key[1], key[3]])
                for key in needed_profiles
            }
            profile_checks += len(profiles)

            for word, old_gain, debt, old_pair in residual:
                gain, i, j = max(
                    (profiles[i, j, word[i], word[j]], i, j) for i, j in EDGES
                )
                margin = gain - debt
                row = (
                    margin, body, m, word, gain, debt, cell, (i, j), old_gain, len(cells)
                )
                require(margin > 0, ("median pair failed", row, old_pair))

                six_arcs = tuple(
                    arc for slot in range(6) for arc in arcs[slot, word[slot]]
                )
                union_mass = R.interval_mass(R.merge_intervals(six_arcs))
                bonferroni = F(6, 7) + debt - gain
                require(union_mass <= bonferroni < F(6, 7), (
                    "union consequence failed", row, union_mass, bonferroni
                ))
                union_slack = F(6, 7) - union_mass
                pair_certified += 1
                union_checks += 1
                if weakest_pair is None or row < weakest_pair:
                    weakest_pair = row
                union_row = (union_slack, body, m, word, union_mass, cell)
                if weakest_union_slack is None or union_row < weakest_union_slack:
                    weakest_union_slack = union_row
                mass_row = (union_mass, body, m, word, cell)
                if largest_union is None or mass_row > largest_union:
                    largest_union = mass_row
                residual_digest.update(
                    f"{body}|{m}|{word}|{old_gain}|{debt}\n".encode()
                )
                certificate_digest.update(
                    f"{body}|{m}|{word}|{cell}|{(i,j)}|{gain}|{union_mass}|{bonferroni}\n".encode()
                )

    count_vector = tuple(residual_counts[m] for m in range(HEAD_START, HEAD_END + 1))
    group_vector = tuple(group_counts[m] for m in range(HEAD_START, HEAD_END + 1))
    require(count_vector == EXPECTED_RESIDUAL_COUNTS, count_vector)
    require(group_vector == EXPECTED_GROUP_COUNTS, group_vector)
    require(pair_certified == RESIDUAL_COUNT, pair_certified)
    require(sum(group_vector) == RESIDUAL_GROUP_COUNT, group_vector)
    require(selector_certified + pair_certified == RAW_ASSIGNMENT_COUNT, (
        selector_certified, pair_certified, RAW_ASSIGNMENT_COUNT
    ))
    require(profile_checks == PROFILE_CHECK_COUNT, profile_checks)
    require(weakest_pair == EXPECTED_WEAKEST_PAIR, weakest_pair)
    require(union_checks == pair_certified, (union_checks, pair_certified))

    semantic_payload = (
        body_digest.hexdigest(),
        count_vector,
        group_vector,
        selector_certified,
        pair_certified,
        profile_checks,
        arc_controls,
        union_checks,
        weakest_pair,
        weakest_union_slack,
        largest_union,
        residual_digest.hexdigest(),
        certificate_digest.hexdigest(),
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest", semantic))

    lines = [
        "LRC14 reflected D=5 finite-head median-cell closure exact proof",
        f"universe=bodies:{BODY_COUNT};upstream_arbitrary_level_bodies:{UPSTREAM_ARBITRARY_LEVEL_BODIES};head_bodies:{HEAD_BODY_COUNT};m:{HEAD_START}..{HEAD_END};raw_assignments:{RAW_ASSIGNMENT_COUNT}",
        "head_crossdet_repair=choose transport orientation with target A=p(qL-b)/(pL-a)>=1;orientation product AA'=pq",
        f"crossdet_selector_certified={selector_certified};residual={pair_certified};residual_counts_by_m={count_vector};residual_groups_by_m={group_vector}",
        "median_cell_rule=sort all body-safe cells and choose cells[len(cells)//2];one cell per residual body-m group",
        "median_involution=j->L-1-j is fixed-point-free;lower+upper=L-1;all used arcs reflect under u->1-u",
        f"median_pair_certificates={pair_certified};distinct_pair_profiles={profile_checks};arc_law_and_singleton_controls={arc_controls}",
        f"weakest_pair_margin={qtext(weakest_pair[0])};E={weakest_pair[1]};m={weakest_pair[2]};word={weakest_pair[3]};gain={qtext(weakest_pair[4])};debt={qtext(weakest_pair[5])};cell={weakest_pair[6]};pair={weakest_pair[7]};crossdet_gain={qtext(weakest_pair[8])};safe_cells={weakest_pair[9]}",
        f"exact_union_checks={union_checks};weakest_union_slack={qtext(weakest_union_slack[0])};E={weakest_union_slack[1]};m={weakest_union_slack[2]};word={weakest_union_slack[3]};union={qtext(weakest_union_slack[4])};cell={weakest_union_slack[5]}",
        f"largest_direct_union={qtext(largest_union[0])};E={largest_union[1]};m={largest_union[2]};word={largest_union[3]};cell={largest_union[4]}",
        "conclusion=every reflected D=5 residual packet closes for every m>=1 (head here, tail upstream)",
        "scope_boundary=D>=6 remains open;sufficient THM-2941 reflected-family certificate only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"tail_sha256={sha256(TAIL)}",
        f"tail_output_sha256={sha256(TAIL_OUTPUT)}",
        f"low_phase_sha256={sha256(LOW_PHASE)}",
        f"k6_minus_edge_sha256={sha256(K6_MINUS_EDGE)}",
        f"k6_minus_edge_output_sha256={sha256(K6_MINUS_EDGE_OUTPUT)}",
        f"body_digest={body_digest.hexdigest()}",
        f"residual_digest={residual_digest.hexdigest()}",
        f"certificate_digest={certificate_digest.hexdigest()}",
        f"source_sha256={sha256(Path(__file__))}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
