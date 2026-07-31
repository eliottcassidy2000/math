#!/usr/bin/env python3
r"""Exact finite referee for the K6 matching corollary of the pair theorem.

Import the labelled-pair overlap lower bound from the adjacent reflected
``k=1`` referee.  For a six-body packet let ``A_i`` be its six reflected
clauses and, for edge ``ij`` and positive integer base ``b``, put

    w_ij(b) = G(|v_ij|) - (4 S_ij(b)+|v_ij|/2)/b.

The pair theorem proves ``mu(A_i intersect A_j) >= w_ij(b)``.  If ``M`` is a
matching, the pointwise Bonferroni inequality

    1_(union A_i) <= sum_i 1_A_i - sum_(ij in M) 1_(A_i intersect A_j)

holds because a point lying in ``k`` clauses belongs to at most ``floor(k/2)``
selected matching edges.  Hence the strictly sufficient certificate is

    max_M sum_(ij in M) max(0,max_b w_ij(b)) > epsilon(E,q).       (M)

This is only a sufficient test.  Failure of (M) is not a survivor theorem and
does not assert that the six-clause union has mass at least ``6/7``.

For a fixed edge, ``v_ij`` is independent of ``b``.  Between the two real
kinks ``q_i-e_i/L`` and ``q_j-e_j/L``, the penalty divided by ``b`` is
fractional-linear and monotone.  Its positive-integer minimum is therefore at

    {max(1,q_i-1),q_i,max(1,q_j-1),q_j}.               (B)

All edge weights are nonnegative after truncation.  Every matching of ``K6``
extends to a perfect matching, so it suffices to test its 15 perfect matchings.

There is also a clean three-edge tariff.  Suppose a perfect matching pairs
unequal levels, with lower level ``p_a`` and gap ``delta_a`` on edge ``a``.
Taking ``b=p_a`` gives ``S=|v|<=delta_a+1/14``.  Since
``G>=35/2976`` and ``epsilon<=1/39``, the packet closes if

    sum_a (delta_a+1/14)/p_a < 373/174096.              (T)

Indeed ``3*(35/2976)-1/39=(9/2)*(373/174096)``.  An
unequal-level perfect matching exists exactly when no level multiplicity
exceeds three.  This dichotomy is structural, not a closure equivalence.

For the hostile multiplicity-four family below, same-level edges have
``|v|<5/7`` and hence ``G(|v|)=0``.  On a cross-level edge, every candidate
base is at most ``2Q``, while ``S>=|v|>=Q-11/168`` and ``G<=1/49``.  Thus its
penalty is at least ``9(Q-11/168)/(4Q)>1/49`` for every ``Q>=1``.  All edge
weights therefore vanish uniformly; this diagnoses failure of this
certificate, not failure of LRC.

Runtime gates enumerate all 76 matchings and all 15 perfect matchings of K6,
verify extension, exhaust the level-partition criterion on ``{1,2,3,4}^6``,
audit the finite optimizer (B) on three structurally distinct bodies and all
endpoint level pairs in ``{1,2,3}``, fire a positive tariff control, and retain
the multiplicity-four hostile family ``(2Q,Q,Q,Q,Q,2Q)`` where every matching
weight is zero.  The optimizer theorem itself is the preceding monotonicity
argument, not an inference from this hostile audit.
Assertions are explicit ``RuntimeError`` gates and remain live under ``-O``.
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from math import lcm
from pathlib import Path


PAIR = Path(__file__).with_name(
    "lrc14_j7_reflected_nearest_level_matching_tail_referee_thm2941.py"
)
ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_k6_matching_corollary_referee_thm2941.out"
)
EXPECTED_PAIR_SHA256 = "173b0edc01159be5ae7ec8f2c6a0d7d36bae347c67c8d1592c1f0976af6c1fb5"
EXPECTED_SEMANTIC_SHA256 = "a0ab3bf5344312f1bcc40a27c1649c376bc66687cffebc87392c8a323fbeb504"
TARIFF = Q(373, 174096)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


SPEC = spec_from_file_location("reflected_pair_referee", PAIR)
require(SPEC is not None and SPEC.loader is not None, "cannot load pair referee")
R = module_from_spec(SPEC)
SPEC.loader.exec_module(R)


Edge = tuple[int, int]
Matching = tuple[Edge, ...]


def canonical_edge(i: int, j: int) -> Edge:
    return (i, j) if i < j else (j, i)


def perfect_matchings(vertices: tuple[int, ...] = tuple(range(6))) -> tuple[Matching, ...]:
    if not vertices:
        return ((),)
    first = vertices[0]
    rows = []
    for position in range(1, len(vertices)):
        second = vertices[position]
        rest = vertices[1:position] + vertices[position + 1 :]
        for tail in perfect_matchings(rest):
            rows.append((canonical_edge(first, second),) + tail)
    return tuple(sorted(rows))


def all_matchings(vertices: tuple[int, ...] = tuple(range(6))) -> tuple[Matching, ...]:
    if not vertices:
        return ((),)
    first = vertices[0]
    rows = list(all_matchings(vertices[1:]))
    for position in range(1, len(vertices)):
        second = vertices[position]
        rest = vertices[1:position] + vertices[position + 1 :]
        for tail in all_matchings(rest):
            rows.append((canonical_edge(first, second),) + tail)
    return tuple(sorted(set(tuple(sorted(row)) for row in rows)))


PERFECT = perfect_matchings()
ALL = all_matchings()


def epsilon(E: tuple[int, ...], L: int, levels: tuple[int, ...]) -> Q:
    return sum((Q(e, 7 * (q * L - e)) for e, q in zip(E, levels)), Q(0))


def edge_candidates(levels: tuple[int, ...], i: int, j: int) -> tuple[int, ...]:
    return tuple(
        sorted(
            {
                max(1, levels[i] - 1),
                levels[i],
                max(1, levels[j] - 1),
                levels[j],
            }
        )
    )


def raw_edge_weight(
    E: tuple[int, ...], L: int, levels: tuple[int, ...], i: int, j: int, base: int
) -> Q:
    lam_i = Q(levels[i] - base) - Q(E[i], L)
    lam_j = Q(levels[j] - base) - Q(E[j], L)
    size = abs(lam_i) + abs(lam_j)
    velocity = lam_i - lam_j
    require(velocity != 0, ("zero labelled velocity", E, levels, i, j))
    return R.minimum_tent_average(abs(velocity)) - (4 * size + abs(velocity) / 2) / base


def edge_weight(E: tuple[int, ...], L: int, levels: tuple[int, ...], i: int, j: int) -> Q:
    return max(
        Q(0),
        *(raw_edge_weight(E, L, levels, i, j, base) for base in edge_candidates(levels, i, j)),
    )


def weight_table(E: tuple[int, ...], L: int, levels: tuple[int, ...]) -> dict[Edge, Q]:
    return {
        (i, j): edge_weight(E, L, levels, i, j)
        for i, j in combinations(range(6), 2)
    }


def matching_weight(weights: dict[Edge, Q], matching: Matching) -> Q:
    return sum((weights[edge] for edge in matching), Q(0))


def best_perfect_matching(weights: dict[Edge, Q]) -> tuple[Q, Matching]:
    return max((matching_weight(weights, row), row) for row in PERFECT)


def matching_closes(E: tuple[int, ...], levels: tuple[int, ...]) -> tuple[bool, Q, Matching, Q]:
    L = 14 * lcm(*E)
    best, matching = best_perfect_matching(weight_table(E, L, levels))
    debt = epsilon(E, L, levels)
    return best > debt, best, matching, debt


def unequal_matching_rows(levels: tuple[int, ...]) -> tuple[Matching, ...]:
    return tuple(
        row
        for row in PERFECT
        if all(levels[i] != levels[j] for i, j in row)
    )


def tariff_value(levels: tuple[int, ...], matching: Matching) -> Q:
    total = Q(0)
    for i, j in matching:
        low, high = sorted((levels[i], levels[j]))
        require(low < high, ("tariff edge has equal levels", levels, matching, i, j))
        total += Q(high - low, low) + Q(1, 14 * low)
    return total


def qtext(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    require(len(PERFECT) == 15, ("perfect matching count changed", len(PERFECT)))
    require(len(ALL) == 76, ("matching count changed", len(ALL)))
    extension_digest = hashlib.sha256()
    for matching in ALL:
        extensions = tuple(row for row in PERFECT if set(matching) <= set(row))
        require(extensions, ("matching has no perfect extension", matching))
        extension_digest.update(f"{matching}|{extensions}\n".encode())

    require(
        3 * R.TENT_FLOOR - Q(1, 39) == Q(9, 2) * TARIFF,
        "tariff identity changed",
    )
    body_rows = tuple(
        (14 * lcm(*E), E) for E in combinations(range(1, 15), 6)
    )
    sum_ratio_max = max((Q(sum(E), L), E, L) for L, E in body_rows)
    pair_ratio_max = max(
        (Q(abs(e - f), L), E, L, e, f)
        for L, E in body_rows
        for e, f in combinations(E, 2)
    )
    require(sum_ratio_max == (Q(1, 6), (1, 2, 3, 4, 6, 12), 168), ("sum ratio changed", sum_ratio_max))
    require(pair_ratio_max == (Q(11, 168), (1, 2, 3, 4, 6, 12), 168, 1, 12), ("pair ratio changed", pair_ratio_max))
    require(pair_ratio_max[0] < Q(1, 14), "tariff label allowance lost")
    require(Q(9, 4) * (1 - pair_ratio_max[0]) > Q(1, 49), "hostile cross-edge bound lost")
    partition_digest = hashlib.sha256()
    partition_rows = 0
    for levels in product(range(1, 5), repeat=6):
        expected = max(Counter(levels).values()) <= 3
        rows = unequal_matching_rows(levels)
        require(bool(rows) == expected, ("unequal matching criterion failed", levels, rows))
        partition_digest.update(f"{levels}|{len(rows)}\n".encode())
        partition_rows += 1

    optimizer_digest = hashlib.sha256()
    optimizer_rows = 0
    optimizer_bodies = (
        (1, 2, 3, 4, 6, 12),
        (1, 2, 4, 9, 11, 12),
        (3, 5, 7, 10, 13, 14),
    )
    for E in optimizer_bodies:
        L = 14 * lcm(*E)
        for i, j in combinations(range(6), 2):
            for qi, qj in product(range(1, 4), repeat=2):
                levels = (qi, qj)
                full = tuple(
                    (
                        raw_edge_weight(E, L, tuple(qi if k == i else qj if k == j else 1 for k in range(6)), i, j, base),
                        base,
                    )
                    for base in range(1, max(qi, qj) + 8)
                )
                candidates = tuple(
                    (
                        raw_edge_weight(E, L, tuple(qi if k == i else qj if k == j else 1 for k in range(6)), i, j, base),
                        base,
                    )
                    for base in edge_candidates(levels, 0, 1)
                )
                require(max(full)[0] == max(candidates)[0], ("edge optimizer failed", E, i, j, qi, qj))
                optimizer_digest.update(f"{E}|{i}|{j}|{qi}|{qj}|{max(candidates)}\n".encode())
                optimizer_rows += 1

    positive_E = (1, 2, 3, 4, 6, 12)
    positive_levels = (10000, 10000, 10001, 10001, 10002, 10002)
    positive_unequal = unequal_matching_rows(positive_levels)
    positive_tariff, positive_matching = min(
        (tariff_value(positive_levels, row), row) for row in positive_unequal
    )
    require(positive_tariff < TARIFF, ("positive tariff did not fire", positive_tariff))
    positive_result = matching_closes(positive_E, positive_levels)
    require(positive_result[0], ("positive weighted matching did not close", positive_result))

    hostile_E = (1, 2, 4, 9, 11, 12)
    hostile_rows = []
    for scale in (1, 2, 3, 7, 16, 64):
        hostile_levels = (2 * scale, scale, scale, scale, scale, 2 * scale)
        require(max(Counter(hostile_levels).values()) == 4, "hostile multiplicity changed")
        require(not unequal_matching_rows(hostile_levels), "hostile gained unequal perfect matching")
        hostile_L = 14 * lcm(*hostile_E)
        weights = weight_table(hostile_E, hostile_L, hostile_levels)
        require(all(value == 0 for value in weights.values()), ("hostile gained weight", scale, weights))
        result = matching_closes(hostile_E, hostile_levels)
        require(not result[0] and result[1] == 0, ("hostile matching unexpectedly closed", scale, result))
        hostile_rows.append((scale, result[1], result[3]))

    semantic_payload = (
        PERFECT,
        len(ALL),
        extension_digest.hexdigest(),
        TARIFF,
        sum_ratio_max,
        pair_ratio_max,
        partition_rows,
        partition_digest.hexdigest(),
        optimizer_rows,
        optimizer_digest.hexdigest(),
        positive_E,
        positive_levels,
        positive_tariff,
        positive_matching,
        positive_result,
        tuple(hostile_rows),
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_PAIR_SHA256 is not None:
        require(sha256(PAIR) == EXPECTED_PAIR_SHA256, ("pair source changed", sha256(PAIR)))
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic changed", semantic))

    lines = [
        "LRC14 reflected k=1 K6 maximum-weight-matching sufficient certificate",
        f"pair_referee_sha256={sha256(PAIR)}",
        "edge_weight=max(0,max_b(G(|v|)-(4S+|v|/2)/b))",
        "edge_base_candidates={max(1,q_i-1),q_i,max(1,q_j-1),q_j}",
        "closure_if=maximum matching weight>epsilon(E,q);sufficient only",
        "finite_reduction=nonnegative K6 weights;76 matchings extend to 15 perfect matchings",
        f"matchings={len(ALL)};perfect_matchings={len(PERFECT)};extension_digest_sha256={extension_digest.hexdigest()}",
        f"body_constant_max_sumE_over_L={sum_ratio_max};max_pair_label_gap_over_L={pair_ratio_max}",
        f"optimizer_rows={optimizer_rows};optimizer_digest_sha256={optimizer_digest.hexdigest()}",
        f"partition_rows={partition_rows};partition_digest_sha256={partition_digest.hexdigest()}",
        "unequal_perfect_matching_iff=maximum level multiplicity<=3",
        f"three_edge_tariff=sum(delta+1/14)/p<{qtext(TARIFF)}",
        f"positive_levels={positive_levels};tariff={qtext(positive_tariff)};matching={positive_matching};best_weight={qtext(positive_result[1])};epsilon={qtext(positive_result[3])};status=CLOSES",
        f"multiplicity_four_hostile=E={hostile_E};levels=(2Q,Q,Q,Q,Q,2Q);rows={tuple(hostile_rows)};all_edge_weights=0;status=NO_CERTIFICATE",
        "nonconsequence=NO_CERTIFICATE does not mean physical survival or union mass>=6/7",
        f"semantic_sha256={semantic}",
        "status=PROVED analytic sufficient corollary plus FINITE-EXACT referee controls",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
