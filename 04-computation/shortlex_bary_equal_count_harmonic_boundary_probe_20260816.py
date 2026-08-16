#!/usr/bin/env python3
"""Exact b-ary shortlex equal-count harmonic boundary companion.

This standard-library probe supports the all-b analytic theorem candidate:
for every b>=2, selecting the first or last 1/b block of every positive
shortlex level gives identical level counts b^(n-1), but distinct logarithmic
densities.  It also records the ternary depth-two Boolean K4 face and the
finite-harmonic-mass boundary for any bounded-width branch ray.
"""

from __future__ import annotations

import ast
from fractions import Fraction
import hashlib
from itertools import combinations, product
from pathlib import Path
import sys


sys.set_int_max_str_digits(0)


EXPECTED_SEMANTIC_SHA256 = "8889cd4d5d3da976a099824e6242fef1a9b9f80fe83e8ae1a4c6fba5dbdf517f"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def level_start(base: int, depth: int) -> int:
    return 1 + (base**depth - 1) // (base - 1)


def level_interval(base: int, depth: int) -> tuple[int, int]:
    start = level_start(base, depth)
    return start, start + base**depth - 1


def edge_blocks(base: int, depth: int) -> tuple[tuple[int, int], tuple[int, int]]:
    start, end = level_interval(base, depth)
    count = base ** (depth - 1)
    return (start, start + count - 1), (end - count + 1, end)


def harmonic_interval(left: int, right: int) -> Fraction:
    return sum((Fraction(1, value) for value in range(left, right + 1)), Fraction(0))


def target_ratios(base: int) -> tuple[Fraction, Fraction]:
    return Fraction(2 * base - 1, base), Fraction(base * base, base * base - base + 1)


def exact_level_ledger() -> tuple[object, ...]:
    rows = []
    for base in range(2, 9):
        left_ratio, right_ratio = target_ratios(base)
        require(left_ratio > right_ratio, ("target order", base))
        require(
            (2 * base - 1) * (base * base - base + 1) - base**3
            == (base - 1) ** 3,
            ("ratio gap identity", base),
        )
        for depth in range(1, 7):
            start, end = level_interval(base, depth)
            left, right = edge_blocks(base, depth)
            require(end + 1 == level_start(base, depth + 1), (base, depth, "seam"))
            require(left[1] - left[0] + 1 == right[1] - right[0] + 1 == base ** (depth - 1),
                    (base, depth, "count"))
            require(left[0] == start and right[1] == end, (base, depth, "boundary"))
            require(left[1] < right[0], (base, depth, "edge blocks overlap"))
            if base**depth <= 120_000:
                left_mass = harmonic_interval(*left)
                right_mass = harmonic_interval(*right)
                require(left_mass > right_mass, (base, depth, "mass order"))
                rows.append((base, depth, start, end, left, right, left_mass, right_mass))
    return tuple(rows)


def stage_domination_ledger() -> tuple[object, ...]:
    rows = []
    for base in range(2, 9):
        earlier = 0
        for stage in range(9):
            length = base ** (2**stage)
            ratio = Fraction(earlier, length)
            rows.append((base, stage, length, earlier, ratio))
            earlier += length
        require(rows[-1][-1] < Fraction(1, 10**20), ("stage domination", base, rows[-1]))
    return tuple(rows)


def ternary_boolean_face() -> tuple[object, ...]:
    words = tuple(product(range(3), repeat=2))
    require(len(words) == 9, len(words))
    t9_arcs = tuple((words[i], words[j]) for i in range(9) for j in range(i + 1, 9))
    require(len(t9_arcs) == 36, len(t9_arcs))
    face = tuple(word for word in words if set(word) <= {0, 2})
    require(face == ((0, 0), (0, 2), (2, 0), (2, 2)), face)
    t4_arcs = tuple((face[i], face[j]) for i in range(4) for j in range(i + 1, 4))
    require(len(t4_arcs) == 6, len(t4_arcs))
    binary = {word: (word[0] // 2, word[1] // 2) for word in face}
    covectors = ((1, 0), (0, 1), (1, 1))
    matchings = []
    covered = []
    for covector in covectors:
        fibres = {}
        for word in face:
            value = sum(a * b for a, b in zip(covector, binary[word])) % 2
            fibres.setdefault(value, []).append(word)
        matching = tuple(tuple(pair) for fibre in sorted(fibres) for pair in combinations(fibres[fibre], 2))
        require(len(matching) == 2, (covector, matching))
        matchings.append((covector, matching))
        covered.extend(tuple(sorted(edge)) for edge in matching)
    require(len(set(covered)) == 6, covered)
    require(len(tuple(combinations(t4_arcs, 2))) == 15, "six edge objects need 15 comparisons")
    return words, t9_arcs, face, t4_arcs, tuple(matchings), tuple(sorted(covered))


def ray_mass_bounds() -> tuple[object, ...]:
    rows = []
    for base in range(2, 13):
        for width in (1, 2, 5):
            # start_n >= b^n/(b-1), so width/start_n <= width*(b-1)/b^n.
            geometric_bound = Fraction(width * (base - 1), 1) * Fraction(1, base - 1)
            require(geometric_bound == width, (base, width, geometric_bound))
            partial = sum(
                Fraction(width, level_start(base, depth))
                for depth in range(1, 20)
            )
            require(partial < geometric_bound, (base, width, partial, geometric_bound))
            rows.append((base, width, partial, geometric_bound))
    return tuple(rows)


level_rows = exact_level_ledger()
stage_rows = stage_domination_ledger()
ternary = ternary_boolean_face()
ray_rows = ray_mass_bounds()

formula_rows = tuple(
    (
        base,
        target_ratios(base),
        (2 * base - 1) * (base * base - base + 1) - base**3,
    )
    for base in range(2, 101)
)
semantic = (
    formula_rows,
    level_rows,
    stage_rows,
    ternary,
    ray_rows,
    "all-b proof uses Riemann sums with c=1/(b-1), summable O(b^-n) errors, and superdominant stages",
)
semantic_hash = hashlib.sha256(repr(semantic).encode("utf-8")).hexdigest()
if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
    require(semantic_hash == EXPECTED_SEMANTIC_SHA256,
            ("semantic hash", semantic_hash))

tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")

print("== b-ary shortlex equal-count harmonic boundary probe ==")
print("level_n=[1+(b^n-1)/(b-1), next_start-1]; first/last blocks each have b^(n-1) indices")
print("limiting level masses=(log((2b-1)/b),log(b^2/(b^2-b+1)))")
print("strict gap certificate=(2b-1)(b^2-b+1)-b^3=(b-1)^3>0")
print("log-density candidates=level_mass/log(b); superdominant alternating stages destroy existence")
print(f"exact_level_rows={len(level_rows)}; level_ledger_sha256={hashlib.sha256(repr(level_rows).encode('utf-8')).hexdigest()}")
print(f"stage_rows={len(stage_rows)}; final_stage_ratio_gate=(bases=2..8,count=7,each<1/10^20)")
print(f"ternary_depth2=(vertices={len(ternary[0])},lex_T9_arcs={len(ternary[1])},boolean_face={ternary[2]},T4_arcs={len(ternary[3])},matchings={ternary[4]})")
print("six K4 edges are edge objects; a separate T6 needs 15 comparisons and an added orientation")
print("bounded-width branch ray has finite harmonic mass; in particular a one-word-per-level Fibonacci/Berggren ray has logarithmic density zero")
print(f"ray_bound_rows={len(ray_rows)}; semantic_sha256={semantic_hash}")
print("scope=shortlex/harmonic theorem candidate only; no ancestry transport, LRC, Keller, or automatic tournament conclusion")
print("all exact checks passed")
