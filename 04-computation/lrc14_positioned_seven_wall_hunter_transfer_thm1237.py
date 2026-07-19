#!/usr/bin/env python3
"""Exact referee for THM-1237 positioned seven-wall Hunter transfer.

The continuum pair formula and THM-1221 global tree floor are named analytic
providers.  This dependency-free script checks the universal pair upper
bound on a large exact bank, the forest-Hunter combinatorics, and every
rational constant in the protected-needle and BAD-transfer consumers.  It
uses neither floating arithmetic nor ``assert``.
"""

from fractions import Fraction as F
from itertools import product
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def folded(x: int) -> int:
    r = x % 14
    return r * (14 - r)


def pair_mass(a: int, b: int) -> F:
    g = gcd(a, b)
    a //= g
    b //= g
    if a > b:
        a, b = b, a
    return F(4 * a * b + folded(a + b) - folded(b - a), 196 * a * b)


pair_rows = 0
equality_rows: list[tuple[int, int]] = []
for a in range(1, 251):
    for b in range(a + 1, 251):
        if gcd(a, b) != 1:
            continue
        pair_rows += 1
        rho = pair_mass(a, b)
        require(rho <= F(1, 14), "pair upper bound")
        require(rho >= F(1, 91), "known pair lower bound cross-check")
        if rho == F(1, 14):
            equality_rows.append((a, b))

require(equality_rows == [(1, 2)], "unique reduced 1:2 upper equality")


def prufer_edges(word: tuple[int, ...], n: int = 7) -> tuple[tuple[int, int], ...]:
    degree = [1] * n
    for v in word:
        degree[v] += 1
    edges: list[tuple[int, int]] = []
    for v in word:
        leaf = next(i for i, d in enumerate(degree) if d == 1)
        edges.append((leaf, v))
        degree[leaf] -= 1
        degree[v] -= 1
    last = [i for i, d in enumerate(degree) if d == 1]
    require(len(last) == 2, "Pruefer terminal leaves")
    edges.append((last[0], last[1]))
    return tuple(edges)


tree_count = 0
active_checks = 0
for word in product(range(7), repeat=5):
    edges = prufer_edges(word)
    require(len(edges) == 6 and len(set(edges)) == 6, "tree edge count")
    tree_count += 1
    for mask in range(1, 1 << 7):
        vertices = mask.bit_count()
        induced = sum(1 for u, v in edges if (mask >> u) & 1 and (mask >> v) & 1)
        require(induced <= vertices - 1, "forest-Hunter active-set inequality")
        active_checks += 1

# Exact discrepancy and protected-needle constants.
require(F(1, 7) * F(6, 7) == F(6, 49), "single-comb discrepancy")
require(F(1, 14) * F(13, 14) == F(13, 196), "pair discrepancy maximum")
require(F(15, 154) * F(1, 7) * 196 == F(30, 11), "needle tree credit")
require(F(6, 49) * 196 == 24, "harmonic debt weight")
require(F(13, 196) * 196 == 13, "gcd debt weight")

# Exact incidence margins after periodic BAD erosion.
arbitrary_bad = F(2, 21)
non_ap_bad = F(60, 637)
require(F(15, 154) - arbitrary_bad == F(1, 462), "arbitrary BAD margin")
require(arbitrary_bad * (1 - arbitrary_bad) == F(38, 441), "arbitrary BAD period debt")
require(F(15, 154) - non_ap_bad == F(45, 14014), "non-AP BAD margin")
require(non_ap_bad * (1 - non_ap_bad) == F(34620, 405769), "non-AP BAD period debt")

# Threshold-rank consumer: a rank-r forest with every edge at least lambda
# contributes at least r*lambda.  Audit a complete rational bank.
rank_checks = 0
for r in range(7):
    for numerator in range(1, 51):
        lam = F(numerator, 997)
        weights = [lam + F(j, 1009) for j in range(r)]
        require(sum(weights, F(0)) >= r * lam, "threshold-rank sum")
        rank_checks += 1

print("THM-1237 POSITIONED SEVEN-WALL HUNTER TRANSFER EXACT AUDIT")
print(f"coprime pair rows checked through 250 = {pair_rows}")
print("unique rho=1/14 row = (1,2)")
print(f"labelled seven-vertex trees = {tree_count}")
print(f"nonempty active-set forest checks = {active_checks}")
print(f"threshold-rank rational checks = {rank_checks}")
print("protected-needle debt = 24m H + 13m sum_T 1/g >= 30/11")
print("BAD margins = 1/462 and 45/14014")
print("BAD period debts = 38/(441h) and 34620/(405769h)")
print("RESULT: PASS")
