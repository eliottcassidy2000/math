#!/usr/bin/env python3
"""Exact companion for THM-2323.

All polynomial arithmetic is over Q (in fact Z).  No floating-point
evaluation of roots of unity is used.
"""

from fractions import Fraction
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def trim(poly):
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul(left, right):
    out = [Fraction(0) for _ in range(len(left) + len(right) - 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return trim(out)


def poly_divmod(numerator, denominator):
    den = trim([Fraction(x) for x in denominator])
    rem = trim([Fraction(x) for x in numerator])
    require(den != [0], "zero polynomial divisor")
    if len(rem) < len(den):
        return [Fraction(0)], rem
    quotient = [Fraction(0) for _ in range(len(rem) - len(den) + 1)]
    while len(rem) >= len(den) and rem != [0]:
        shift = len(rem) - len(den)
        scale = rem[-1] / den[-1]
        quotient[shift] += scale
        for j, value in enumerate(den):
            rem[shift + j] -= scale * value
        rem = trim(rem)
    return trim(quotient), trim(rem)


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def cyclotomic(n, cache=None):
    if cache is None:
        cache = {}
    if n in cache:
        return cache[n]
    numerator = [Fraction(-1)] + [Fraction(0)] * (n - 1) + [Fraction(1)]
    denominator = [Fraction(1)]
    for d in divisors(n):
        if d < n:
            denominator = poly_mul(denominator, cyclotomic(d, cache))
    quotient, remainder = poly_divmod(numerator, denominator)
    require(remainder == [0], f"cyclotomic division failed at n={n}")
    require(all(value.denominator == 1 for value in quotient),
            f"nonintegral cyclotomic coefficient at n={n}")
    result = [int(value) for value in quotient]
    cache[n] = result
    return result


def reduce_mod(poly, modulus):
    _, remainder = poly_divmod(poly, modulus)
    return trim(remainder)


def euler_phi(n):
    result = n
    p = 2
    remaining = n
    while p * p <= remaining:
        if remaining % p == 0:
            while remaining % p == 0:
                remaining //= p
            result -= result // p
        p += 1
    if remaining > 1:
        result -= result // remaining
    return result


def mobius_squarefree(n):
    sign = 1
    p = 2
    remaining = n
    while p * p <= remaining:
        if remaining % p == 0:
            remaining //= p
            sign = -sign
            if remaining % p == 0:
                return 0
            while remaining % p == 0:
                remaining //= p
        p += 1
    if remaining > 1:
        sign = -sign
    return sign


def ramanujan_sum(n, m):
    d = gcd(n, m)
    q = n // d
    return mobius_squarefree(q) * euler_phi(n) // euler_phi(q)


def field_polynomial(terms, n):
    """Return sum coefficient*zeta^exponent as a polynomial mod X^n-1."""
    out = [Fraction(0) for _ in range(n)]
    for exponent, coefficient in terms:
        out[exponent % n] += Fraction(coefficient)
    return trim(out)


N = 91
R = N // 7 - 1
PHI = cyclotomic(N)
UNITS = [k for k in range(1, N) if gcd(k, N) == 1]

require(R == 12, "wrong 91-needle correlation radius")
require(set(range(-R, R + 1))
        == {s - r for r in range(13) for s in range(13)},
        "thirteen-tooth difference support mismatch")
require(len(UNITS) == 72, "wrong primitive colour count")
require(len(PHI) - 1 == euler_phi(N) == 72, "wrong Phi_91 degree")
require(2 * R == 24 < euler_phi(N), "cyclotomic width inequality failed")

# A deterministic rational word/bare needle: 0 <= word <= bare.
bare = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9]
word_sites = {0, 2, 5, 7, 11}
word = [value if r in word_sites else 0 for r, value in enumerate(bare)]
require(all(0 <= a <= b for a, b in zip(word, bare)),
        "word/bare order failed")

correlation = {}
for d in range(-R, R + 1):
    correlation[d] = sum(
        word[r] * bare[r + d]
        for r in range(13)
        if 0 <= r + d < 13
    )

diagonal = correlation[0]
require(diagonal == sum(a * b for a, b in zip(word, bare)) > 0,
        "positive diagonal failed")

fixed_colour_nonzero = 0
for k in UNITS:
    terms = [(k * d, coefficient)
             for d, coefficient in correlation.items()]
    remainder = reduce_mod(field_polynomial(terms, N), PHI)
    require(remainder != [0],
            f"deterministic fixed-colour correlation vanished at k={k}")
    fixed_colour_nonzero += 1
require(fixed_colour_nonzero == 72, "not all primitive colours checked")

# THM-2319's aggregate-zero cross-energy control.
# Fixed colour: 1 + 12*zeta^(7k), nonzero for every primitive k.
hostile_fixed_nonzero = 0
hostile_sum = [Fraction(0)]
for k in UNITS:
    polynomial = field_polynomial([(0, 1), (7 * k, 12)], N)
    remainder = reduce_mod(polynomial, PHI)
    require(remainder != [0],
            f"aggregate hostile has a fixed-colour zero at k={k}")
    hostile_fixed_nonzero += 1
    if len(hostile_sum) < len(polynomial):
        hostile_sum += [Fraction(0)] * (len(polynomial) - len(hostile_sum))
    for exponent, coefficient in enumerate(polynomial):
        hostile_sum[exponent] += coefficient

require(hostile_fixed_nonzero == 72, "hostile fixed-colour count mismatch")
require(ramanujan_sum(91, 7) == -6, "wrong Ramanujan value c_91(7)")
require(72 + 12 * ramanujan_sum(91, 7) == 0,
        "aggregate hostile scalar cancellation failed")
require(reduce_mod(hostile_sum, PHI) == [0],
        "aggregate hostile field cancellation failed")

# The all-depth 7*13^a stalk.  The displayed difference is exact for all a;
# sample many levels as an independent arithmetic control.
for a in range(1, 21):
    n = 7 * 13**a
    radius = 13**a - 1
    phi_n = euler_phi(n)
    require(phi_n == 72 * 13 ** (a - 1),
            f"totient formula failed at a={a}")
    require(phi_n - 2 * radius == 46 * 13 ** (a - 1) + 2,
            f"stalk width identity failed at a={a}")
    require(phi_n > 2 * radius,
            f"stalk cyclotomic inequality failed at a={a}")

# LRC jump and height ledgers.
for S in range(1, 101):
    j_word = 6 * S
    j_bare = 2 * S
    product_bound = j_word * j_bare
    h_max = product_bound - 1
    n_max = 90 + 91 * h_max
    require(product_bound == 12 * S * S,
            f"jump-product bound failed at S={S}")
    require(h_max == 12 * S * S - 1,
            f"gauge bound failed at S={S}")
    require(n_max == 1092 * S * S - 1,
            f"frequency bound failed at S={S}")

# CRT fibre graphs used in the conditional incidence corollary.
def adjacent_91(left, right):
    return left != right and gcd((right - left) % 91, 91) == 1


full_vertices = list(range(91))
deleted_vertices = [z for z in full_vertices if z % 7 != 0]
full_edges = sum(
    adjacent_91(left, right)
    for left in full_vertices
    for right in full_vertices
) // 2
deleted_edges = sum(
    adjacent_91(left, right)
    for left in deleted_vertices
    for right in deleted_vertices
) // 2

require(all(
    adjacent_91(left, right)
    == (left % 7 != right % 7 and left % 13 != right % 13)
    for left in full_vertices
    for right in full_vertices
), "CRT product adjacency failed")
require(full_edges == 91 * 72 // 2, "full CRT fibre edge count failed")
require(deleted_edges == 78 * 60 // 2, "deleted CRT fibre edge count failed")
require(all(adjacent_91(left, right)
            for left in range(7) for right in range(left + 1, 7)),
        "K_7 clique control failed")
require(all(adjacent_91(left, right)
            for left in range(1, 7) for right in range(left + 1, 7)),
        "K_6 clique control failed")
require(all(
    not adjacent_91(left, right) or left % 7 != right % 7
    for left in full_vertices
    for right in full_vertices
), "seven-colouring control failed")
require(all(
    not adjacent_91(left, right) or left % 7 != right % 7
    for left in deleted_vertices
    for right in deleted_vertices
), "six-colouring control failed")


def coherent_lift_control(alpha, delta, r, a, z_left, z_right, h_left, w):
    require(delta >= 1 and gcd(r, 91) == 1, "bad cofactor data")
    g0 = 7**alpha * 13**delta
    d_prime = g0 * r
    modulus = 91 * g0
    require(euler_phi(modulus) - 2 * (modulus // 7 - 1)
            == 46 * g0 + 2, "conditional modulus gap failed")
    k0 = 1
    k_left = k0 + d_prime * z_left
    k_right = k0 + d_prime * z_right
    h_right = h_left + r * w
    q_left = k_left + modulus * h_left
    q_right = k_right + modulus * h_right
    bracket = (z_right - z_left) + 91 * w
    require(adjacent_91(z_left, z_right), "chosen fibre edge is not unit")
    require(q_right - q_left == d_prime * bracket,
            "coherent-lift factorization failed")
    require(gcd(bracket, 91) == 1, "bracket lost unit colour")
    require(gcd(a, 91) == 1 and gcd(a, d_prime) == 1,
            "physical multiplier hypotheses failed")
    g = 13**3
    c2 = g * a
    c3 = g * d_prime
    multiplier = a * bracket
    require(gcd(c2, c3) == g, "normalized gcd data failed")
    require(c2 * (q_right - q_left) == multiplier * c3,
            "physical c3-edge identity failed")
    require(gcd(multiplier, 91) == 1,
            "physical multiplier lost unit colour")


coherent_lift_control(
    alpha=1, delta=2, r=5, a=11,
    z_left=0, z_right=1, h_left=2, w=3,
)
coherent_lift_control(
    alpha=0, delta=2, r=5, a=11,
    z_left=0, z_right=1, h_left=2, w=3,
)

print("THM-2323 exact companion")
print(f"N={N}, R={R}, 2R={2*R}, phi(N)={euler_phi(N)}")
print(f"primitive_colours={len(UNITS)}")
print(f"deterministic_diagonal={diagonal}")
print(f"deterministic_fixed_colour_nonzero={fixed_colour_nonzero}")
print(f"c_91(7)={ramanujan_sum(91, 7)}")
print("aggregate_hostile_sum=0")
print(f"aggregate_hostile_fixed_colour_nonzero={hostile_fixed_nonzero}")
print("stalk_gap(a)=46*13^(a-1)+2 > 0 for all a>=1")
print("lrc_common_h<=12*S^2-1")
print("lrc_common_n<=1092*S^2-1")
print("full_fibre=K7xK13, vertices=91, degree=72, chi=7")
print("deleted_fibre=K6xK13, vertices=78, degree=60, chi=6")
print("conditional_thresholds: alpha>=1,r<=6; alpha=0,r<=5")
print("coherent_lift_controls=2")
print("all exact checks passed")
