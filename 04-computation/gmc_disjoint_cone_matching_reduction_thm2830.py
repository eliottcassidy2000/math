#!/usr/bin/env python3
"""Exact referee for the THM-2830 matching reduction.

This companion verifies a bounded universe and exact algebraic controls.
It does not prove the universal matching inequality.
"""

from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import comb, factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def h(a, b):
    return comb(a + b, a)


def multinomial3(a, b, c):
    return factorial(a + b + c) // (factorial(a) * factorial(b) * factorial(c))


def triple(a, b, c):
    return (
        multinomial3(a + 1, b, c)
        + multinomial3(a, b + 1, c)
        + multinomial3(a, b, c + 1)
        - multinomial3(a, b, c)
    )


def tau(a, b, c):
    total = a + b + c
    return (
        Fraction(total + 1, a + 1)
        + Fraction(total + 1, b + 1)
        + Fraction(total + 1, c + 1)
        - 1
    )


def p_sum(i, values):
    ans = 0
    for slot in range(4):
        rest = values[:slot] + values[slot + 1 :]
        ans += h(i, values[slot]) * triple(*rest)
    return ans


def matching_gap(i, a, b, c, d):
    values = (a, b, c, d)
    return p_sum(i, values) - 3 * (
        triple(i, a, b) * h(c, d) + triple(i, c, d) * h(a, b)
    )


def normalized_matching_gap(i, a, b, c, d):
    return (
        factorial(i)
        * factorial(a)
        * factorial(b)
        * factorial(c)
        * factorial(d)
        * matching_gap(i, a, b, c, d)
    )


def all_pair_sum(i, values):
    total = 0
    for mask in range(16):
        if mask.bit_count() != 2:
            continue
        left = [values[j] for j in range(4) if (mask >> j) & 1]
        right = [values[j] for j in range(4) if not ((mask >> j) & 1)]
        total += triple(i, *left) * h(*right)
    return total


def polarized(i, values):
    return Fraction(p_sum(i, values) - all_pair_sum(i, values), 2)


def poly_mul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for j, x in enumerate(a):
        for k, y in enumerate(b):
            out[j + k] += x * y
    return out


def poly_add(a, b):
    out = [Fraction(0)] * max(len(a), len(b))
    for j in range(len(out)):
        out[j] = (a[j] if j < len(a) else 0) + (b[j] if j < len(b) else 0)
    return out


def poly_scale(scale, poly):
    return [scale * value for value in poly]


def moment(poly):
    return sum(x * factorial(j) for j, x in enumerate(poly))


def d_poly(index):
    out = [Fraction(0)] * (index + 2)
    out[index] = -Fraction(1, factorial(index))
    out[index + 1] = Fraction(1, factorial(index + 1))
    return out


def direct_triple(a, b, c):
    return moment(poly_mul(poly_mul(d_poly(a), d_poly(b)), d_poly(c)))


def direct_orientation(i, upper_atoms):
    v = [Fraction(0)]
    for atom in upper_atoms:
        v = poly_add(v, d_poly(atom))
    v2 = poly_mul(v, v)
    v3 = poly_mul(v2, v)
    di = d_poly(i)
    return (
        2 * moment(v3) * moment(poly_mul(di, v))
        - 3 * moment(poly_mul(di, v2)) * moment(v2)
    )


def direct_polarized(i, values):
    total = Fraction(0)
    for mask in range(16):
        chosen = [values[j] for j in range(4) if (mask >> j) & 1]
        sign = -1 if (4 - len(chosen)) % 2 else 1
        total += sign * direct_orientation(i, chosen)
    return total / factorial(4)


def pair_product_control(p, q):
    left = poly_mul(d_poly(p), d_poly(q))
    degree = p + q
    right = [Fraction(0)] * len(left)
    right[degree] += Fraction(comb(degree, p), factorial(degree))
    atom = d_poly(degree + 1)
    scale = comb(degree + 2, p + 1)
    for j, value in enumerate(atom):
        right[j] += scale * value
    return left == right


def falling(x, order):
    if x < order:
        return 0
    return factorial(x) // factorial(x - order)


def pair_factorial_hostile():
    i, y, z, order = 0, 3, 44, 4
    total = y + z
    rho = Fraction(196, 705)
    low = Fraction(falling(total - 1, order), 5)
    spike = falling(total + 1, order)
    phi = (low + rho * spike) / (1 + rho)
    endpoint = Fraction(falling(y, order) + falling(z, order), 2)
    return phi, endpoint, phi - endpoint


def laguerre_hostile(t):
    return 2 * t * (5640 * t * t + 371 * t - 3)


for a in range(7):
    for b in range(7):
        for c in range(7):
            require(triple(a, b, c) == direct_triple(a, b, c), "triple tensor")

for p in range(8):
    for q in range(8):
        require(pair_product_control(p, q), "pair product")

require(
    tau(3, 4, 5)
    == Fraction(
        triple(3, 4, 5) * factorial(3) * factorial(4) * factorial(5),
        factorial(12),
    ),
    "tau normalization",
)

polarization_controls = 0
for i in range(3):
    for values in combinations_with_replacement(range(i + 1, i + 5), 4):
        require(polarized(i, values) == direct_polarized(i, values), "polarization")
        polarization_controls += 1

matching_cells = 0
strict_cells = 0
equality_cells = 0
forward_cells = 0
minimum_positive = None

for i in range(5):
    values_range = range(i + 1, i + 9)
    for values in combinations_with_replacement(values_range, 4):
        a, b, c, d = values
        matchings = (
            (a, b, c, d),
            (a, c, b, d),
            (a, d, b, c),
        )
        gaps = []
        for matching in matchings:
            gap = matching_gap(i, *matching)
            matching_cells += 1
            require(gap >= 0, "matching gap")
            gaps.append(gap)
            if gap == 0:
                equality_cells += 1
                require(values == (i + 1,) * 4, "unexpected equality")
            else:
                strict_cells += 1
                if minimum_positive is None or gap < minimum_positive:
                    minimum_positive = gap
        require(polarized(i, values) >= 0, "polarized coefficient")
        require(sum(gaps) == 3 * (2 * polarized(i, values)), "matching average")

        for slot in range(4):
            raised = list(values)
            raised[slot] += 1
            before = normalized_matching_gap(i, a, b, c, d)
            after = normalized_matching_gap(i, *raised)
            forward_cells += 1
            require(after >= before, "normalized forward difference")

phi, endpoint, hostile_gap = pair_factorial_hostile()
require(phi == Fraction(1467522360, 901), "pair hostile phi")
require(endpoint == 1629012, "pair hostile endpoint")
require(hostile_gap == -Fraction(217452, 901), "pair hostile gap")

laguerre_value = laguerre_hostile(Fraction(1, 1000))
require(laguerre_value < 0, "Laguerre hostile")

cone_cells = 0
cone_equalities = 0
for boundary in range(1, 5):
    lower_vectors = product(range(3), repeat=boundary)
    for lower in lower_vectors:
        if not any(lower):
            continue
        for upper in product(range(3), repeat=3):
            if not any(upper):
                continue
            u = [Fraction(0)]
            v = [Fraction(0)]
            for index, coefficient in enumerate(lower):
                u = poly_add(u, poly_scale(coefficient, d_poly(index)))
            for offset, coefficient in enumerate(upper):
                v = poly_add(v, poly_scale(coefficient, d_poly(boundary + offset)))
            v2 = poly_mul(v, v)
            value = (
                2 * moment(poly_mul(v2, v)) * moment(poly_mul(u, v))
                - 3 * moment(poly_mul(u, v2)) * moment(v2)
            )
            cone_cells += 1
            require(value >= 0, "cone orientation")
            if value == 0:
                cone_equalities += 1
                require(
                    all(coefficient == 0 for coefficient in lower[:-1])
                    and lower[-1] > 0
                    and upper[0] > 0
                    and upper[1:] == (0, 0),
                    "cone equality",
                )

print("THM-2830 DISJOINT-CONE MATCHING REDUCTION -- exact referee")
print("status=FINITE-EXACT REDUCTION; UNIVERSAL MATCHING SIGN OPEN")
print("triple_tensor_cells=343")
print("pair_product_cells=64")
print(f"polarization_controls={polarization_controls}")
print(f"matching_cells={matching_cells}")
print(f"matching_strict={strict_cells}")
print(f"matching_equalities={equality_cells}")
print(f"normalized_forward_cells={forward_cells}")
print(f"minimum_positive_matching_gap={minimum_positive}")
print(f"pair_factorial_phi={phi}")
print(f"pair_factorial_endpoint={endpoint}")
print(f"pair_factorial_gap={hostile_gap}")
print(f"laguerre_t=1/1000_gap={laguerre_value}")
print(f"cone_cells={cone_cells}")
print(f"cone_equalities={cone_equalities}")
print("all_exact_controls=PASS")
