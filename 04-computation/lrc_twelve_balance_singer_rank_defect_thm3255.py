#!/usr/bin/env python3
"""Exact companion for THM-3255's multiplicative Singer rank defect."""

import ast
from fractions import Fraction as F
from hashlib import sha256
from math import gcd, lcm
from pathlib import Path


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT / "01-canon/theorems/THM-3246-all-dilation-second-owner-seam-stabilization-and-sign-word.md":
        "6badc0c9aba09b56d3d055a96cb8ef8b619d8492508bf21476eba5f624b13055",
    ROOT / "01-canon/theorems/THM-3252-singer-compactified-owner-hodge-word-universal-charged-cyclicity.md":
        "1f8797de2d5fac74814fb78ca4f4d500de8c42eb14a6e1721e5f3e2a2810a873",
    ROOT / "01-canon/theorems/THM-3253-positive-owner-mass-newton-cyclicity-and-maximal-common-heisenberg-module.md":
        "b94aea11abe97a6cc1a3826a91fab59d4c04e15f0e6acd9c924c5463b7bd63e8",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


for dependency, expected in DEPENDENCIES.items():
    require(sha256(lf_bytes(dependency)).hexdigest() == expected,
            ("dependency hash drift", dependency.name))

tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(tree)
)
require(assert_nodes == 0, "optimization-sensitive assert")
require(float_literals == 0, "floating literal")


def trim(poly):
    poly = list(poly)
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return tuple(poly)


def poly_mul(left, right):
    answer = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] += a * b
    return trim(answer)


def poly_divmod_monic(dividend, divisor):
    dividend = list(trim(dividend))
    divisor = trim(divisor)
    require(divisor[-1] == 1, "monic divisor")
    if len(dividend) < len(divisor):
        return (0,), tuple(dividend)
    quotient = [0] * (len(dividend) - len(divisor) + 1)
    while len(dividend) >= len(divisor):
        shift = len(dividend) - len(divisor)
        coefficient = dividend[-1]
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            dividend[index + shift] -= coefficient * value
        dividend = list(trim(dividend))
    return trim(quotient), trim(dividend)


def exact_quotient(dividend, divisor):
    quotient, remainder = poly_divmod_monic(dividend, divisor)
    require(remainder == (0,), ("nonzero exact remainder", remainder))
    return quotient


def divisors(number):
    return tuple(value for value in range(1, number + 1)
                 if number % value == 0)


def prime_divisors(number):
    answer = []
    candidate = 2
    while candidate * candidate <= number:
        if number % candidate == 0:
            answer.append(candidate)
            while number % candidate == 0:
                number //= candidate
        candidate += 1
    if number > 1:
        answer.append(number)
    return tuple(answer)


# Build all cyclotomic polynomials needed below without a CAS.
cyclotomic = {}
for order in divisors(168):
    polynomial = tuple([-1] + [0] * (order - 1) + [1])
    for proper in divisors(order):
        if proper < order:
            polynomial = exact_quotient(polynomial, cyclotomic[proper])
    cyclotomic[order] = polynomial
require(poly_mul(poly_mul(poly_mul(poly_mul(cyclotomic[2], cyclotomic[3]),
                                           cyclotomic[4]), cyclotomic[6]),
                              cyclotomic[12]) == (1,) * 12,
        "C12 cyclotomic factorization")


def mass_coefficients(cell):
    """Return the (g^2,g,1) coefficients of THM-3246's cleared mass."""
    if cell <= 5 or cell >= 162:
        return 12096, -1032, 2
    if 6 <= cell <= 23 or 144 <= cell <= 161:
        return 12096, -24, 0
    if 24 <= cell <= 71:
        return 16044 - 168 * cell, 48, 0
    if 72 <= cell <= 95:
        return 4032, 96, 0
    if 96 <= cell <= 143:
        return 168 * cell - 12012, 48, 0
    raise RuntimeError(("cell", cell))


mass_triples = tuple(mass_coefficients(cell) for cell in range(168))
residue_triples = tuple(
    tuple(sum(mass_triples[cell][degree]
              for cell in range(residue, 168, 12))
          for degree in range(3))
    for residue in range(12)
)
require(set(residue_triples) == {(120960, -528, 2)},
        "uniform twelve-residue mass")
c12 = (1,) * 12
for degree in range(3):
    exact_quotient(tuple(row[degree] for row in mass_triples), c12)

missing_orders = (2, 3, 4, 6, 12)
missing_frequencies = tuple(range(14, 168, 14))
require(sum(len(cyclotomic[order]) - 1 for order in missing_orders) == 11,
        "missing degree")
require(tuple(168 // gcd(168, frequency)
              for frequency in missing_frequencies) ==
        (12, 6, 4, 3, 12, 2, 12, 3, 4, 6, 12),
        "missing character orders")


def vector_remainder(order):
    remainders = []
    for degree in range(3):
        _, remainder = poly_divmod_monic(
            tuple(row[degree] for row in mass_triples), cyclotomic[order])
        remainders.append(remainder)
    width = max(map(len, remainders))
    return tuple(tuple(remainders[degree][index]
                       if index < len(remainders[degree]) else 0
                       for degree in range(3))
                 for index in range(width))


# The two rational-root orders are handled over Z, not by a misleading
# finite-field root-free search.
special_eight = tuple((0, 1008 * sign, -2 * sign)
                      for sign in (1, 1, -1, -1))
special_twenty_four = tuple((0, 2016 * sign, -4 * sign)
                            for sign in (-1, -1, -1, -1, 0, 0, 1, 1))
require(vector_remainder(8) == special_eight,
        "Phi8 remainder -2(504g-1)(x-1)(x+1)^2")
require(vector_remainder(24) == special_twenty_four,
        "Phi24 remainder 4(504g-1)(x+1)(x6-x2-1)")


def is_prime(number):
    if number < 2:
        return False
    divisor = 2
    while divisor * divisor <= number:
        if number % divisor == 0:
            return False
        divisor += 1
    return True


def evaluate_component(component, value, prime):
    answer = 0
    power = 1
    for row in mass_triples:
        answer = (answer + row[component] * power) % prime
        power = power * value % prime
    return answer


# (order, prime, primitive generator, selected root, A, B, C, discriminant).
finite_certificates = (
    (7, 29, 2, 16, 14, 0, 16, 3),
    (14, 43, 3, 27, 33, 34, 18, 27),
    (21, 211, 2, 180, 173, 17, 163, 167),
    (28, 113, 3, 81, 8, 91, 46, 29),
    (42, 127, 3, 27, 47, 2, 76, 67),
    (56, 113, 3, 9, 98, 90, 71, 43),
    (84, 337, 10, 227, 162, 270, 146, 197),
    (168, 673, 5, 625, 147, 130, 250, 462),
)
for order, prime, generator, root, aa, bb, cc, discriminant in finite_certificates:
    require(is_prime(prime), ("certificate prime", order, prime))
    require(all(pow(generator, (prime - 1) // divisor, prime) != 1
                for divisor in prime_divisors(prime - 1)),
            ("primitive generator", order, prime))
    require(root == pow(generator, (prime - 1) // order, prime),
            ("selected primitive root", order))
    require(pow(root, order, prime) == 1 and
            all(pow(root, order // divisor, prime) != 1
                for divisor in prime_divisors(order)),
            ("root exact order", order))
    observed = tuple(evaluate_component(component, root, prime)
                     for component in range(3))
    require(observed == (aa, bb, cc),
            ("Fourier quadratic", order, observed))
    require(discriminant == (bb * bb - 4 * aa * cc) % prime,
            ("quadratic discriminant", order))
    require(pow(discriminant, (prime - 1) // 2, prime) == prime - 1,
            ("nonsquare discriminant", order))

covered_orders = (1,) + missing_orders + (8, 24) + tuple(
    row[0] for row in finite_certificates)
require(set(covered_orders) == set(divisors(168)), "all cyclotomic orders")


def mass_word(dilation):
    return tuple(a * dilation * dilation + b * dilation + c
                 for a, b, c in mass_triples)


for dilation in (1, 2, 18, 80):
    zero_orders = tuple(order for order in divisors(168)
                        if poly_divmod_monic(mass_word(dilation),
                                             cyclotomic[order])[1] == (0,))
    require(zero_orders == missing_orders,
            ("direct mass zero orders", dilation, zero_orders))

# The trivial character is positive: total mass numerator is
# 24(60480g^2-264g+1), increasing and positive from g=1 onward.
require(tuple(sum(row[degree] for row in mass_triples)
              for degree in range(3)) == (1451520, -6336, 24),
        "trivial character polynomial")
require(60480 - 264 + 1 > 0 and 2 * 60480 - 264 > 0,
        "trivial character positivity")


# Reconstruct THM-3252's signed Hodge word independently.
L0 = 168
THETA = F(1, 14)


def frac(value):
    return value - value.numerator // value.denominator


def b3bar(value):
    value = frac(value)
    return value ** 3 - F(3, 2) * value ** 2 + F(1, 2) * value


def comb(frequency, phase=F(0)):
    phase = frac(phase)
    intervals = []
    for integer in range(-2, frequency + 3):
        lower = max(F(0), (F(integer) + phase - THETA) / frequency)
        upper = min(F(1), (F(integer) + phase + THETA) / frequency)
        if lower < upper:
            intervals.append((lower, upper))
    return tuple(sorted(intervals))


def centered_overlap(left, right):
    i = j = 0
    total = F(0)
    while i < len(left) and j < len(right):
        lower = max(left[i][0], right[j][0])
        upper = min(left[i][1], right[j][1])
        if lower < upper:
            total += ((upper - F(1, 2)) ** 2
                      - (lower - F(1, 2)) ** 2) / 2
        if left[i][1] < right[j][1]:
            i += 1
        elif right[j][1] < left[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return total


def barycenter(p, alpha, q, beta):
    return centered_overlap(comb(p, alpha), comb(q, beta))


def correction(cell):
    return (barycenter(3, F(cell + 1, L0), 5, F(2 * (cell + 1), L0))
            - barycenter(3, F(cell, L0), 5, F(2 * cell, L0)))


def limit(cell):
    cross = -1
    rr, ss = cell % L0, 2 * cell % L0
    difference = 5 * rr - 3 * ss
    uu, vv = F(difference + cross, L0), F(-difference, L0)
    aa, bb = F(3, 14), F(5, 14)
    psi = (b3bar(uu + aa - bb) + b3bar(uu - aa + bb)
           + b3bar(vv + aa - bb) + b3bar(vv - aa + bb)
           - b3bar(uu + aa + bb) - b3bar(uu - aa - bb)
           - b3bar(vv + aa + bb) - b3bar(vv - aa - bb))
    return F(1, 49) + F(28, 15 * cross) * psi


owner_q = tuple(
    (F(2 if cell <= 5 or cell >= 162 else 0)
     - 2 * limit(cell) + 1848 * correction(cell)) / 423360
    for cell in range(L0)
)
owner_digest = sha256("\n".join(map(str, owner_q)).encode("ascii")).hexdigest()
require(owner_digest ==
        "b53d77a69f39a5f8c893b4cdceaaeecdd0dc70a16be02e6d28f3f6d7b520feef",
        "owner Hodge word digest")
scale = 1
for value in owner_q:
    scale = lcm(scale, value.denominator)
scaled_q = tuple(int(value * scale) for value in owner_q)
require(scale == 32006016000 and sum(scaled_q) == 1296000,
        "scaled Hodge word")
require(tuple(sum(scaled_q[residue::12]) for residue in range(12))
        == (108000,) * 12, "uniform Hodge residue sums")
exact_quotient(scaled_q, c12)

q_remainders = []
q_zero_orders = []
for order in divisors(168):
    remainder = poly_divmod_monic(scaled_q, cyclotomic[order])[1]
    q_remainders.append((order, remainder))
    if remainder == (0,):
        q_zero_orders.append(order)
require(tuple(q_zero_orders) == missing_orders, "Hodge exact zero orders")
q_remainder_digest = sha256("\n".join(
    "%d:%s" % (order, ",".join(map(str, remainder)))
    for order, remainder in q_remainders
).encode("ascii")).hexdigest()


# Unit decimations preserve the complete missing-character set.  A single
# labeled phase atom has coefficient one in every character and repairs it;
# a constant word does not repair any nontrivial missing character.
units = tuple(value for value in range(168) if gcd(value, 168) == 1)
require(all(tuple(sorted(multiplier * frequency % 168
                         for frequency in missing_frequencies))
            == missing_frequencies for multiplier in units),
        "unit invariance of missing characters")
delta_phase = (1,) + (0,) * 167
constant_word = (1,) * 168
require(all(poly_divmod_monic(delta_phase, cyclotomic[order])[1] != (0,)
            for order in missing_orders), "one phase atom repairs all modes")
require(all(poly_divmod_monic(constant_word, cyclotomic[order])[1] == (0,)
            for order in missing_orders), "constant correction does not repair")


print("THM-3255 TWELVE-BALANCE MULTIPLICATIVE SINGER DEFECT EXACT AUDIT")
print("dependency_hash_checks=%d" % len(DEPENDENCIES))
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("mass_residue_sums=12*(120960*g^2-528*g+2)")
print("common_factor=(x^12-1)/(x-1)=Phi2*Phi3*Phi4*Phi6*Phi12,degree=11")
print("special_orders=8:-2*(504g-1)*(x-1)*(x+1)^2;24:4*(504g-1)*(x+1)*(x^6-x^2-1)")
print("finite_root_free_certificates=%s" % (finite_certificates,))
print("direct_mass_zero_orders_g=(1,2,18,80):%s" % (missing_orders,))
print("all_integer_g>=1_mass_phase_rank=168-11=157")
print("hodge_digest=%s,scale=%d,residue_sums=12*108000" %
      (owner_digest, scale))
print("hodge_zero_orders=%s,phase_rank=157" % (missing_orders,))
print("hodge_remainder_digest=%s" % q_remainder_digest)
print("all_8064_translates_decimations_share_missing_frequencies=%s" %
      (missing_frequencies,))
print("joint_span_criterion=sidecar_nonzero_on_all_11_missing_characters")
print("single_labeled_phase_atom_repairs=PASS,constant_correction_repairs=FAIL")
print("scope=coefficient_plane_phase-translation-no-canonical-LRC-owner-marker")
print("all_exact_checks=PASS")
