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
require(4 * R < N, "91-needle is not in the acute Galois sector")

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

# The acute-sector inequality is independent of the prime factorization.
for a in range(1, 21):
    n = 7 * 13**a
    radius = 13**a - 1
    require(4 * radius < n,
            f"stalk acute-sector inequality failed at a={a}")

for cofactor in (1, 2, 3, 5, 6, 11, 30, 330, 2310):
    n = 7 * cofactor
    radius = n // 7 - 1
    require(4 * radius < n,
            f"arbitrary-prime acute sector failed at cofactor={cofactor}")

# N=210 is an exact control beyond the old degree comparison:
# 2R=58 >= phi(210)=48, yet every primitive fixed colour is nonzero.
N_WIDE = 210
R_WIDE = N_WIDE // 7 - 1
PHI_WIDE = cyclotomic(N_WIDE)
UNITS_WIDE = [k for k in range(1, N_WIDE) if gcd(k, N_WIDE) == 1]
bare_wide = [1 + (3 * r + r * r) % 11 for r in range(R_WIDE + 1)]
word_wide = [
    value if r % 4 in (0, 1) else 0
    for r, value in enumerate(bare_wide)
]
correlation_wide = {}
for d in range(-R_WIDE, R_WIDE + 1):
    correlation_wide[d] = sum(
        word_wide[r] * bare_wide[r + d]
        for r in range(R_WIDE + 1)
        if 0 <= r + d <= R_WIDE
    )
require(correlation_wide[0] > 0, "wide control lost its diagonal")
require(2 * R_WIDE >= euler_phi(N_WIDE),
        "wide control does not cross the old degree boundary")
wide_fixed_nonzero = 0
for k in UNITS_WIDE:
    terms = [(k * d, coefficient)
             for d, coefficient in correlation_wide.items()]
    remainder = reduce_mod(field_polynomial(terms, N_WIDE), PHI_WIDE)
    require(remainder != [0],
            f"wide fixed-colour correlation vanished at k={k}")
    wide_fixed_nonzero += 1
require(wide_fixed_nonzero == euler_phi(N_WIDE) == 48,
        "wide primitive-colour count mismatch")

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


def direct_incidence_control(d_prime, a, h_left, h_right):
    require(d_prime % 13 == 0, "d' must contain the deep 13 factor")
    require(gcd(a, d_prime) == 1 and gcd(a, 91) == 1,
            "physical multiplier hypotheses failed")
    modulus = 91 * d_prime
    k0 = 1
    require(1 <= k0 < d_prime and gcd(k0, d_prime) == 1,
            "bad K0 representative")
    retained = [
        z for z in range(91)
        if gcd(k0 + d_prime * z, modulus) == 1
    ]
    expected_count = 91 if d_prime % 7 == 0 else 78
    require(len(retained) == expected_count,
            "wrong primitive coherent-fibre size")
    require(all(1 <= k0 + d_prime * z < modulus for z in retained),
            "K_z left its canonical range")
    pair = next(
        (left, right)
        for left in retained
        for right in retained
        if adjacent_91(left, right)
    )
    z_left, z_right = pair
    k_left = k0 + d_prime * z_left
    k_right = k0 + d_prime * z_right
    q_left = k_left + modulus * h_left
    q_right = k_right + modulus * h_right
    bracket = (z_right - z_left) + 91 * (h_right - h_left)
    require(adjacent_91(z_left, z_right), "chosen fibre edge is not unit")
    require(q_right - q_left == d_prime * bracket,
            "direct-modulus factorization failed")
    require(gcd(bracket, 91) == 1, "bracket lost unit colour")
    g = 13**3
    c2 = g * a
    c3 = g * d_prime
    multiplier = a * bracket
    require(gcd(c2, c3) == g, "normalized gcd data failed")
    require(c2 * (q_right - q_left) == multiplier * c3,
            "physical c3-edge identity failed")
    require(gcd(multiplier, 91) == 1,
            "physical multiplier lost unit colour")
    S = 1
    jump_product = 12 * S * S
    require(0 <= h_left < jump_product and 0 <= h_right < jump_product,
            "sample gauges exceed the LRC product bound")
    require(1 <= q_left <= modulus * jump_product - 1,
            "left atom exceeds its positive bound")
    require(1 <= q_right <= modulus * jump_product - 1,
            "right atom exceeds its positive bound")
    require(abs(bracket) <= 91 * jump_product - 1,
            "normalized edge bound failed")
    require(abs(multiplier) <= a * (1092 * S * S - 1),
            "physical edge bound failed")
    return len(retained)


full_direct_count = direct_incidence_control(
    d_prime=7 * 13**2 * 330,
    a=17,
    h_left=2,
    h_right=5,
)
deleted_direct_count = direct_incidence_control(
    d_prime=13**2 * 330,
    a=17,
    h_left=2,
    h_right=5,
)
require(full_direct_count == 91 and deleted_direct_count == 78,
        "direct incidence fibre controls failed")

# Exact failure boundary: if seven divides a, every normalized multiplier
# a*B is seven-divisible, regardless of the unit bracket B.
obstructed_a = 7
unit_brackets = [b for b in range(-200, 201) if gcd(b, 91) == 1]
require(unit_brackets, "empty unit-bracket control")
require(all(gcd(obstructed_a * b, 91) > 1 for b in unit_brackets),
        "seven-divisible-a obstruction failed")

print("THM-2323 exact companion")
print(f"N={N}, R={R}, 2R={2*R}, phi(N)={euler_phi(N)}")
print(f"primitive_colours={len(UNITS)}")
print(f"deterministic_diagonal={diagonal}")
print(f"deterministic_fixed_colour_nonzero={fixed_colour_nonzero}")
print(f"c_91(7)={ramanujan_sum(91, 7)}")
print("aggregate_hostile_sum=0")
print(f"aggregate_hostile_fixed_colour_nonzero={hostile_fixed_nonzero}")
print("galois_sector: |d|<=N/7-1 gives |arg(zeta^d)|<2*pi/7<pi/2")
print("arbitrary_prime_control: N=210, 2R=58>=phi(N)=48")
print(f"arbitrary_prime_fixed_colour_nonzero={wide_fixed_nonzero}")
print("lrc_common_h<=12*S^2-1")
print("lrc_common_n<=1092*S^2-1")
print("full_fibre=K7xK13, vertices=91, degree=72")
print("deleted_fibre=K6xK13, vertices=78, degree=60")
print("direct_N=91*d_prime incidence controls=2")
print("normalized_unit_edge_bound<=1092*S^2-1")
print("seven_divisible_a_obstruction=exact")
print("all exact checks passed")
