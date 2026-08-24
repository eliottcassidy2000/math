#!/usr/bin/env python3
"""Exact discovery probe for adjoining the ramified prime 3 to THM-3793."""

from __future__ import annotations

import hashlib
import itertools
import math
from fractions import Fraction


MAX_SUM = 2500
COMP_SUM = math.floor((4 ** (1 / 3)) * MAX_SUM) + 2
GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def factor(n: int) -> dict[int, int]:
    ans: dict[int, int] = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            ans[p] = ans.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        ans[n] = ans.get(n, 0) + 1
    return ans


def valuation(n: int, p: int) -> int:
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def admissible(x: int, y: int, cap3: int) -> bool:
    g = math.gcd(x, y)
    d = x + y
    s = d // g
    for p in factor(d):
        if p == 3:
            if valuation(s, 3) > cap3:
                return False
        elif p % 3 == 2:
            if valuation(s, p) > 2:
                return False
        else:
            return False
    return True


def primitive_three_invoice(x: int, y: int) -> tuple[int, int]:
    g = math.gcd(x, y)
    x //= g
    y //= g
    q = x * x - x * y + y * y
    e = valuation(x + y, 3)
    return e, valuation(q, 3)


strict: dict[int, tuple[int, int]] = {}
relaxed: dict[int, tuple[int, int]] = {}
strict_rows = relaxed_rows = 0
for d in range(3, MAX_SUM + 1):
    row_strict = row_relaxed = False
    for x in range(1, (d - 1) // 2 + 1):
        y = d - x
        m = x**3 + y**3
        if admissible(x, y, 1):
            gate(m not in strict, f"strict construction collision before global scan m={m}")
            strict[m] = (x, y)
            row_strict = True
        if admissible(x, y, 2):
            relaxed.setdefault(m, (x, y))
            row_relaxed = True
    strict_rows += int(row_strict)
    relaxed_rows += int(row_relaxed)

# Any competing positive representation of a candidate with sum d<=MAX_SUM
# has d'<(4m)^(1/3)<4^(1/3)d, so this sum universe is complete.
strict_hostiles: list[tuple[int, tuple[int, int], tuple[int, int]]] = []
relaxed_hostiles: list[tuple[int, tuple[int, int], tuple[int, int]]] = []
strict_hits = {m: 0 for m in strict}
for d in range(3, COMP_SUM + 1):
    for x in range(1, (d - 1) // 2 + 1):
        y = d - x
        m = x**3 + y**3
        if m in strict:
            strict_hits[m] += 1
        if m in strict and strict[m] != (x, y):
            strict_hostiles.append((m, strict[m], (x, y)))
        if m in relaxed and relaxed[m] != (x, y):
            relaxed_hostiles.append((m, relaxed[m], (x, y)))

gate(not strict_hostiles, "strict cap-three-one singleton universe")
gate(all(count == 1 for count in strict_hits.values()), "every strict candidate encountered exactly once")
gate(bool(relaxed_hostiles), "relaxed cap-three-two has hostile")

# Check the exact primitive q valuation law independently on the full box.
for x in range(1, 1001):
    for y in range(x + 1, 1002):
        if math.gcd(x, y) != 1:
            continue
        e, q3 = primitive_three_invoice(x, y)
        gate(q3 == int(e > 0), f"primitive three invoice x={x},y={y}")

first = min(relaxed_hostiles)
m, pair1, pair2 = first
gate(
    m == 4104 and {pair1, pair2} == {(2, 16), (9, 15)},
    "sharp ramified exponent-two hostile",
)

# Exact finite Euler-band identity on an independently enumerated small bank.
bank = (5, 11)
main_mass = Fraction(0)
square_error = Fraction(0)
row_surrogate = Fraction(0)
band_values: set[int] = set()
for order in range(1, len(bank) + 1):
    for choice in itertools.combinations(bank, order):
        d = math.prod(choice)
        main_mass += Fraction(1, d)
        square_error += Fraction(1, d * d)
        row_surrogate += Fraction(d - 1, 2 * d * d)
        row_surrogate += Fraction(3 * d - 1, 18 * d * d)
        for pair_sum in (d, 3 * d):
            for x in range(1, (pair_sum - 1) // 2 + 1):
                y = pair_sum - x
                value = x**3 + y**3
                gate(value in strict, f"ramified Euler band candidate d={pair_sum},x={x}")
                gate(value not in band_values, f"ramified Euler band disjoint m={value}")
                band_values.add(value)
gate(
    row_surrogate == Fraction(2, 3) * main_mass - Fraction(5, 9) * square_error,
    "exact two-row Euler mass identity",
)
gate(
    square_error
    == math.prod((Fraction(1) + Fraction(1, p * p) for p in bank), start=Fraction(1))
    - 1,
    "square-weight error Euler product",
)

g1 = math.gcd(*pair1)
g2 = math.gcd(*pair2)
semantic = (
    "ramified prime 3 is admissible with primitive-sum exponent at most one",
    "inert primes retain primitive-sum exponent cap two",
    "every admissible candidate has a globally singleton positive distinct two-cube fibre",
    "the exponent-two relaxation at prime 3 has an exact hostile",
    "paired d and 3d rows change the Euler main coefficient from two-fifths to two-thirds",
)
print(
    "UNIVERSE",
    f"candidate_sum<={MAX_SUM};competitor_sum<={COMP_SUM};"
    f"strict_rows={strict_rows};strict_values={len(strict)};"
    f"relaxed_rows={relaxed_rows};relaxed_values={len(relaxed)}",
)
print(
    "FIRST_CAP3_TWO_HOSTILE",
    f"m={m};pair1={pair1};sum1={sum(pair1)};g1={g1};s1={sum(pair1)//g1};"
    f"pair2={pair2};sum2={sum(pair2)};g2={g2};s2={sum(pair2)//g2}",
)
print("RELAXED_HOSTILE_INCIDENCES", len(relaxed_hostiles))
print(
    "RAMIFIED_EULER_BAND",
    f"bank={bank};values={len(band_values)};surrogate={row_surrogate};"
    f"main={main_mass};square_error={square_error};coefficient_gain=5/3",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic).encode()).hexdigest())
print("GATES", GATES)
print("RESULT PASS")
