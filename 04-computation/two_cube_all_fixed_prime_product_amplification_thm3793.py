#!/usr/bin/env python3
"""Exact probe for all-fixed-order inert-prime amplification of THM-3793.

This companion checks the elementary-symmetric collision inequality, finite
squarefree inert product rows, and the sharp exponent-three hostile.  The
all-scale conclusion is the displayed combinatorial identity plus the
reciprocal-prime bound already proved in THM-3730.
"""

from __future__ import annotations

import ast
import hashlib
import json
from fractions import Fraction
from itertools import combinations
from pathlib import Path


GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    if not condition:
        raise RuntimeError(label)
    GATES += 1


def primes_through(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            sieve[p * p:limit + 1:p] = b"\x00" * (((limit - p * p) // p) + 1)
    return [p for p in range(2, limit + 1) if sieve[p]]


def elementary_symmetric(weights: list[Fraction], degree: int) -> Fraction:
    values = [Fraction(0) for _ in range(degree + 1)]
    values[0] = Fraction(1)
    for weight in weights:
        for j in range(degree, 0, -1):
            values[j] += weight * values[j - 1]
    return values[degree]


inert = [p for p in primes_through(997) if p >= 5 and p % 3 == 2]
weights = [Fraction(1, p) for p in inert]
A = sum(weights, Fraction(0))
B = sum((weight * weight for weight in weights), Fraction(0))

rows: list[dict[str, str | int]] = []
for j in range(1, 9):
    e_j = elementary_symmetric(weights, j)
    ordered_distinct = e_j
    for factor in range(2, j + 1):
        ordered_distinct *= factor
    if j == 1:
        collision_lower = A
    else:
        collision_lower = A**j - Fraction(j * (j - 1), 2) * B * A ** (j - 2)
    gate(ordered_distinct >= collision_lower,
         f"ordered collision union bound j={j}")
    gate(e_j > 0, f"positive elementary symmetric mass j={j}")
    rows.append({
        "j": j,
        "prime_count": len(inert),
        "e_j_num_digits": len(str(e_j.numerator)),
        "e_j_den_digits": len(str(e_j.denominator)),
        "collision_margin_sign": 1 if ordered_distinct > collision_lower else 0,
    })

# Independent finite row-incidence control.  Squarefree inert products have
# pairwise-disjoint constructed values; this is the finite shadow of the
# singleton theorem's cross-row conclusion.
small_primes = [5, 11, 17, 23]
seen: dict[int, tuple[int, int, int]] = {}
finite_rows = 0
finite_values = 0
for j in range(1, 4):
    for subset in combinations(small_primes, j):
        d = 1
        for p in subset:
            d *= p
        gate(all(d % p == 0 and (d // p) % p != 0 for p in subset),
             f"squarefree row d={d}")
        gate(all(p % 3 == 2 for p in subset), f"inert row d={d}")
        row_values: set[int] = set()
        for x in range(1, (d + 1) // 2):
            y = d - x
            value = x**3 + y**3
            gate(value < d**3, f"row height d={d},x={x}")
            gate(value not in row_values, f"within-row injectivity d={d}")
            gate(value not in seen, f"cross-row injectivity d={d},x={x}")
            row_values.add(value)
            seen[value] = (d, x, y)
        gate(len(row_values) == (d - 1) // 2,
             f"complete odd row cardinality d={d}")
        row_mass_floor = Fraction(d - 1, 2 * d * d)
        gate(row_mass_floor >= Fraction(2, 5 * d),
             f"uniform 2/(5d) row floor d={d}")
        finite_rows += 1
        finite_values += len(row_values)

gate(54**3 + 71**3 == 15**3 + 80**3 == 515375,
     "primitive exponent-three collision hostile")
gate(54 + 71 == 125 and 15 + 80 == 95,
     "hostile has distinct pair sums")
gate(125 == 5**3, "hostile begins exactly at inert exponent three")

semantic = {
    "family": "d is a product of j distinct primes p=2 mod 3; all x<d-x rows",
    "height": "x^3+(d-x)^3<d^3, so j primes <=Z land below Z^(3j)",
    "row_mass": "sum row m^(-2/3)>=(d-1)/(2d^2)>=2/(5d)",
    "symmetric_bound": "j=1 separately; for j>=2, j!*e_j>=A^j-binomial(j,2)*B*A^(j-2)",
    "asymptotic": "liminf H(X)/(log log X)^j>=2/(5*2^j*j!) for every fixed j>=1",
    "consequence": "H(X) dominates every fixed power of log log X",
    "scope": "explicit singleton subfamily only; no support asymptotic or critical residue",
    "finite": {"rows": finite_rows, "values": finite_values, "j_max": 3},
}
blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "no inactive Python assert")

print("probe=two_cube_all_fixed_prime_product_amplification")
print(f"finite_prime_bank={small_primes};rows={finite_rows};values={finite_values}")
print(f"symmetric_bank=inert_primes_through_997:{len(inert)};orders=1..8")
print("bound_j1=H(Z^3)>=(2/5)*A(Z)")
print("bound_j_ge_2=H(Z^(3j))>=(2/(5*j!))*(A(Z)^j-C(j,2)B(Z)A(Z)^(j-2))")
print("asymptotic=for_every_fixed_j>=1:liminf_H(X)/(loglogX)^j>=2/(5*2^j*j!)")
print("growth=H(X)_dominates_every_fixed_power_of_loglogX")
print("scope=singleton_subfamily_lower_bound_only;support_asymptotic_and_residue_OPEN")
print(f"semantic_sha256={hashlib.sha256(blob).hexdigest()}")
print(f"GATES={GATES}")
print("PASS")
