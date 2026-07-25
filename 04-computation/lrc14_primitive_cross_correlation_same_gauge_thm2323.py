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

# The exact radius works without any seven-divisibility hypothesis.
for n in range(2, 501):
    radius = (n + 6) // 7 - 1
    require(0 <= radius and 7 * radius < n,
            f"arbitrary-modulus support radius failed at n={n}")
    require(4 * radius < n,
            f"arbitrary-modulus acute sector failed at n={n}")

# N=30 is a direct non-seven-divisible control, at the old degree
# boundary 2R=phi(N).
N_NONSEVEN = 30
R_NONSEVEN = (N_NONSEVEN + 6) // 7 - 1
PHI_NONSEVEN = cyclotomic(N_NONSEVEN)
UNITS_NONSEVEN = [
    k for k in range(1, N_NONSEVEN) if gcd(k, N_NONSEVEN) == 1
]
bare_nonseven = [2, 7, 1, 8, 2]
word_nonseven = [2, 0, 1, 0, 2]
require(len(bare_nonseven) == R_NONSEVEN + 1,
        "nonseven needle length mismatch")
correlation_nonseven = {}
for d in range(-R_NONSEVEN, R_NONSEVEN + 1):
    correlation_nonseven[d] = sum(
        word_nonseven[r] * bare_nonseven[r + d]
        for r in range(R_NONSEVEN + 1)
        if 0 <= r + d <= R_NONSEVEN
    )
nonseven_fixed_nonzero = 0
for k in UNITS_NONSEVEN:
    terms = [
        (k * d, coefficient)
        for d, coefficient in correlation_nonseven.items()
    ]
    remainder = reduce_mod(
        field_polynomial(terms, N_NONSEVEN), PHI_NONSEVEN
    )
    require(remainder != [0],
            f"nonseven fixed-colour correlation vanished at k={k}")
    nonseven_fixed_nonzero += 1
require(2 * R_NONSEVEN == euler_phi(N_NONSEVEN) == 8,
        "nonseven old-degree boundary mismatch")
require(nonseven_fixed_nonzero == euler_phi(N_NONSEVEN),
        "nonseven primitive-colour count mismatch")

# The D_a automorphic form at a!=1.  The sites are indexed by their
# short signed a-residues, so every surviving displacement obeys
# |signed(a*d)|<N/7.  Galois straightening is then exact over Phi_N.
def signed_residue(value, modulus):
    residue = value % modulus
    if 2 * residue > modulus:
        residue -= modulus
    return residue


N_ARITH = 65
A_ARITH = 3
PHI_ARITH = cyclotomic(N_ARITH)
UNITS_ARITH = [
    k for k in range(1, N_ARITH) if gcd(k, N_ARITH) == 1
]
INV_A_ARITH = pow(A_ARITH, -1, N_ARITH)
short_residues = list(range(-4, 5))
arith_sites = [(e * INV_A_ARITH) % N_ARITH for e in short_residues]
bare_arith = {
    site: 1 + (5 * index + index * index) % 11
    for index, site in enumerate(arith_sites)
}
word_arith = {
    site: value if index % 3 != 1 else 0
    for index, (site, value) in enumerate(bare_arith.items())
}
require(all(0 <= word_arith[r] <= bare_arith[r] for r in arith_sites),
        "arithmetic-comb word/bare order failed")
correlation_arith = {}
for d in range(N_ARITH):
    coefficient = sum(
        word_arith[r] * bare_arith.get((r + d) % N_ARITH, 0)
        for r in arith_sites
    )
    if coefficient:
        correlation_arith[d] = coefficient
        displacement = signed_residue(A_ARITH * d, N_ARITH)
        require(7 * abs(displacement) < N_ARITH,
                f"arithmetic-comb support escaped at d={d}")
require(correlation_arith.get(0, 0) > 0,
        "arithmetic-comb diagonal vanished")
arith_fixed_nonzero = 0
for k in UNITS_ARITH:
    terms = [
        (k * d, coefficient)
        for d, coefficient in correlation_arith.items()
    ]
    remainder = reduce_mod(field_polynomial(terms, N_ARITH), PHI_ARITH)
    require(remainder != [0],
            f"arithmetic-comb correlation vanished at k={k}")
    straightened = [
        (A_ARITH * d, coefficient)
        for d, coefficient in correlation_arith.items()
    ]
    automorphism = A_ARITH * pow(k, -1, N_ARITH)
    mapped_terms = [
        (automorphism * exponent, coefficient)
        for exponent, coefficient in terms
    ]
    require(field_polynomial(mapped_terms, N_ARITH)
            == field_polynomial(straightened, N_ARITH),
            f"arithmetic-comb Galois map failed at k={k}")
    straight_remainder = reduce_mod(
        field_polynomial(straightened, N_ARITH), PHI_ARITH
    )
    require(straight_remainder != [0],
            f"arithmetic-comb straightening vanished at k={k}")
    arith_fixed_nonzero += 1
require(arith_fixed_nonzero == euler_phi(N_ARITH) == 48,
        "arithmetic-comb primitive-colour count mismatch")

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

# Universal P_g carrier at N=13*d': two primitive colours separated by
# d' give two common word/bare atoms at a thirteen-primitive c3 distance.
# The residue enumeration also checks P_v P_(13^b)=P_g.
P_POWER = 13**2
V_UNIT = 17
composite_residues = {
    s + V_UNIT * r
    for s in range(V_UNIT)
    for r in range(P_POWER)
}
require(composite_residues == set(range(V_UNIT * P_POWER)),
        "P_v P_(13^b) residue enumeration failed")

D_CARRIER = 13 * 5 * 11
N_CARRIER = 13 * D_CARRIER
K0_CARRIER = 2
K1_CARRIER = K0_CARRIER + D_CARRIER
require(gcd(K0_CARRIER, D_CARRIER) == 1,
        "carrier K0 is not primitive over d'")
require(1 <= K0_CARRIER < K1_CARRIER < N_CARRIER,
        "carrier colours left their canonical range")
require(gcd(K0_CARRIER, N_CARRIER) == 1
        and gcd(K1_CARRIER, N_CARRIER) == 1,
        "carrier colours are not primitive")

G_CARRIER = 13**3 * V_UNIT
A_COFACTOR = 3
require(gcd(A_COFACTOR, D_CARRIER) == 1,
        "middle-owner normalized factors are not coprime")
C2_CARRIER = G_CARRIER * A_COFACTOR
C3_CARRIER = G_CARRIER * D_CARRIER
require(gcd(C2_CARRIER, C3_CARRIER) == G_CARRIER,
        "middle-owner common carrier mismatch")
H0_CARRIER = 2
H1_CARRIER = 5
Q0_CARRIER = K0_CARRIER + N_CARRIER * H0_CARRIER
Q1_CARRIER = K1_CARRIER + N_CARRIER * H1_CARRIER
T_CARRIER = 1 + 13 * (H1_CARRIER - H0_CARRIER)
require(Q1_CARRIER - Q0_CARRIER == D_CARRIER * T_CARRIER,
        "13-colour coherent-lift factorization failed")
require(G_CARRIER * (Q1_CARRIER - Q0_CARRIER)
        == T_CARRIER * C3_CARRIER,
        "physical P_g carrier edge failed")
require(T_CARRIER % 13 != 0,
        "P_g carrier multiplier lost thirteen-primitivity")
require((G_CARRIER // 13**3) * K0_CARRIER % 13
        == (G_CARRIER // 13**3) * K1_CARRIER % 13,
        "P_g carrier root character changed")
require(abs(T_CARRIER) <= 156 - 12,
        "P_g carrier multiplier bound failed at S=1")
for S in range(1, 101):
    L = 12 * S * S
    require(all(
        0 < abs(1 + 13 * delta) <= 13 * L - 12
        and (1 + 13 * delta) % 13 != 0
        for delta in range(-(L - 1), L)
    ), f"universal thirteen-unit bound failed at S={S}")

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
print("galois_sector: |d|<=ceil(N/7)-1 gives |arg(zeta^d)|<2*pi/7<pi/2")
print("arbitrary_modulus_control: N=30, 2R=phi(N)=8")
print(f"arbitrary_modulus_fixed_colour_nonzero={nonseven_fixed_nonzero}")
print("arithmetic_comb_control: N=65, a=3, surviving |ad|<N/7")
print(f"arithmetic_comb_fixed_colour_nonzero={arith_fixed_nonzero}")
print("old_degree_failure_control: N=210, 2R=58>=phi(N)=48")
print(f"arbitrary_prime_fixed_colour_nonzero={wide_fixed_nonzero}")
print("lrc_common_h<=12*S^2-1")
print("lrc_common_n<=1092*S^2-1")
print("P_g_carrier=N=13*d_prime, primitive_colours=2")
print("universal_thirteen_unit_edge_bound<=156*S^2-12")
print("full_fibre=K7xK13, vertices=91, degree=72")
print("deleted_fibre=K6xK13, vertices=78, degree=60")
print("direct_N=91*d_prime incidence controls=2")
print("normalized_unit_edge_bound<=1092*S^2-1")
print("seven_divisible_a_obstruction=exact")
print("all exact checks passed")
