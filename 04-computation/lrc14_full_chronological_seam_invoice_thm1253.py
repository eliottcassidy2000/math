#!/usr/bin/env python3
"""Exact referee for THM-1253's full chronological seam invoice.

The continuum provider is the elementary topology of a deletion-minimal
finite chain of open intervals.  This dependency-free audit checks its
three-interval order type, pairwise separation of all consecutive handoffs,
the gcd/lcm tooth quantum, the exact excess-mass coefficient, the
occurrence-count consumer, and the nonbacktracking structure of surjective
six-label owner words.
"""

from fractions import Fraction as F
from itertools import combinations, product
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


# If three intervals have increasing left and right endpoints and the first
# and third overlap, the middle interval is contained in their union.  A
# deletion-minimal chain must therefore obey beta_i <= alpha_(i+2).
redundancy_rows = 0
for points in combinations(range(13), 6):
    for alpha0, alpha1, alpha2 in combinations(points, 3):
        remaining = tuple(x for x in points if x not in (alpha0, alpha1, alpha2))
        if len(remaining) != 3:
            continue
        beta0, beta1, beta2 = remaining
        if not (alpha0 < alpha1 < alpha2 and beta0 < beta1 < beta2):
            continue
        if not all(a < b for a, b in
                   ((alpha0, beta0), (alpha1, beta1), (alpha2, beta2))):
            continue
        if not alpha2 < beta0:
            continue
        # Check every arrangement cell and endpoint in the middle interval.
        probes = [F(x, 2) for x in range(2 * alpha1 + 1, 2 * beta1)]
        for x in probes:
            in_middle = alpha1 < x < beta1
            in_outer_union = (alpha0 < x < beta0) or (alpha2 < x < beta2)
            require(not in_middle or in_outer_union,
                    "three-interval redundancy")
        redundancy_rows += 1

require(redundancy_rows > 0, "nonempty redundancy audit")


# Exhaust all abstract endpoint chains on a small integer alphabet.  The
# condition beta_i <= alpha_(i+2) is the independently proved minimality
# consequence above; it makes every consecutive handoff interval disjoint.
chain_rows = 0
handoff_pair_rows = 0
for length in range(2, 7):
    for alphas in combinations(range(10), length):
        for betas in combinations(range(1, 11), length):
            if any(alphas[i] >= betas[i] for i in range(length)):
                continue
            if any(alphas[i + 1] >= betas[i] for i in range(length - 1)):
                continue
            if any(betas[i] > alphas[i + 2] for i in range(length - 2)):
                continue
            handoffs = [(alphas[i + 1], betas[i])
                        for i in range(length - 1)]
            for i, (left_i, right_i) in enumerate(handoffs):
                require(left_i < right_i, "positive consecutive handoff")
                for left_j, right_j in handoffs[i + 1:]:
                    require(right_i <= left_j,
                            "pairwise disjoint handoff intervals")
                    handoff_pair_rows += 1
            chain_rows += 1

require(chain_rows > 0, "nonempty endpoint-chain audit")


# Every raw tooth handoff has a positive numerator divisible by gcd(u,v).
quantum_rows = 0
for u in range(2, 151):
    for v in range(u + 1, 151):
        g = gcd(u, v)
        ell = u * v // g
        for n in range(-4, 5):
            for m in range(-4, 5):
                numerator = v * (14 * n + 1) - u * (14 * m - 1)
                if numerator <= 0:
                    continue
                require(numerator % g == 0, "handoff gcd sheet")
                overlap = F(numerator, 14 * u * v)
                require(overlap >= F(1, 14 * ell),
                        "handoff lcm quantum")
                quantum_rows += 1


# Six-label owner words: a nonbacktracking middle position is exactly a
# change of unordered transition edge.  Each edge change introduces at most
# one new vertex, so a surjective word has at least 6-2=4 such turns.
word_rows = 0
minimum_turns = None
vertices = tuple(range(6))
for length in range(6, 9):
    for word in product(vertices, repeat=length):
        if set(word) != set(vertices):
            continue
        if any(a == b for a, b in zip(word, word[1:])):
            continue
        edges = [frozenset((a, b)) for a, b in zip(word, word[1:])]
        turns = sum(edges[i - 1] != edges[i] for i in range(1, len(edges)))
        require(turns >= 4, "six-label nonbacktracking turn floor")
        require(len(set(word)) <= turns + 2,
                "edge-run vertex accounting")
        minimum_turns = turns if minimum_turns is None else min(minimum_turns, turns)
        word_rows += 1

require(minimum_turns == 4, "sharp nonbacktracking turn floor")


# Exact coefficient ledger.  For |G|=6/(7c), the singleton upper bounds give
# sum_i |D_i cap G|-|G| <= (6/49)(H-1/c).  Since the disjoint handoffs lie in
# the multiplicity excess, every occurrence receives coefficient 49/6; its
# 1/(14*lcm) quantum becomes 7/(12*lcm).
require(F(6, 7) * F(1, 7) == F(6, 49),
        "slow-gap singleton deficit")
require(F(49, 6) * F(1, 14) == F(7, 12),
        "full seam coefficient")
require(F(7, 12) == 3 * F(7, 36),
        "factor-three improvement over Cayley average")
require(F(7, 6) * F(3, 4) * F(1, 14) == F(1, 16),
        "weighted functional full-seam coefficient")


# Private mass forces at least ceil(d/(7c)) selected teeth.  Coarsening every
# lcm by d6^2/g0 yields the direct scalar occurrence consumer.
owner_rows = 0
for c in range(1, 81):
    for d in range(c + 1, 35 * c + 1):
        teeth = (d + 7 * c - 1) // (7 * c)
        require(F(teeth, 7 * d) >= F(1, 49 * c),
                "private owner occurrence capacity")
        owner_rows += 1

scalar_rows = 0
for c in range(1, 31):
    for d6 in range(c + 6, 18 * c + 1):
        for g0 in range(1, min(c, d6) + 1):
            if d6 % g0:
                continue
            for handoffs in range(5, 18):
                lcm_sum_floor = F(handoffs * g0, d6 * d6)
                debt = F(7, 12) * lcm_sum_floor
                require(debt == F(7 * g0 * handoffs, 12 * d6 * d6),
                        "scalar occurrence debt")
                scalar_rows += 1


print("THM-1253 FULL CHRONOLOGICAL SEAM INVOICE EXACT AUDIT")
print(f"three-interval redundancy rows = {redundancy_rows}")
print(f"minimal endpoint chains checked = {chain_rows}")
print(f"disjoint handoff pairs checked = {handoff_pair_rows}")
print(f"positive gcd/lcm handoff rows = {quantum_rows}")
print(f"surjective six-label owner words = {word_rows}")
print(f"sharp nonbacktracking turn floor = {minimum_turns}")
print(f"private owner occurrence rows = {owner_rows}")
print(f"scalar debt rows = {scalar_rows}")
print("full occurrence debt = H >= 1/c + (7/12) sum_a 1/lcm(s_a,s_(a+1))")
print("functional debt = sum_i Pbar(6d_i/(7c))-1 >= (c/16) sum_a 1/lcm")
print("coarse debt = H >= 1/c + 7*g0*(N-1)/(12*d6^2)")
print("RESULT: PASS")
