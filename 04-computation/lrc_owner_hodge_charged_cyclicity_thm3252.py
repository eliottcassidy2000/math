#!/usr/bin/env python3
"""Exact companion for THM-3252's compactified owner cyclicity theorem."""

import ast
from collections import Counter
from fractions import Fraction as F
from functools import reduce
from hashlib import sha256
from math import gcd, lcm
from pathlib import Path


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT / "01-canon/theorems/THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate.md":
        "ef77a1f8fce16eb851eb38d5110a61ab73aa693f2d0ee9e11a912aa4fc302c87",
    ROOT / "01-canon/theorems/THM-3246-all-dilation-second-owner-seam-stabilization-and-sign-word.md":
        "6badc0c9aba09b56d3d055a96cb8ef8b619d8492508bf21476eba5f624b13055",
    ROOT / "04-computation/lrc_second_owner_all_dilation_seam_thm3246.py":
        "e23b098b38aa2199a348f48f8ab4ac0ce5913c870ead972bd31296494fc25a4b",
    ROOT / "05-knowledge/results/lrc_second_owner_all_dilation_seam_thm3246.out":
        "d7f7dd96b01c597113e78f903cad36246cb47b10e9a1758cb831aa0e83e8cebc",
    ROOT / "01-canon/theorems/THM-3250-charged-heisenberg-blowup-address-intertwiner-and-pointed-multiplicity-gate.md":
        "7e91a07e38d6869e6621b9057594cdf4745827a98de9618681d09427102b27ea",
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

# Independent reconstruction of THM-3246's q-word from the promoted closed
# numerator constants and THM-3224's exact Bernoulli endpoint formulas.
L0 = 168
THETA = F(1, 14)


def frac(x):
    return x - x.numerator // x.denominator


def b3bar(x):
    x = frac(x)
    return x ** 3 - F(3, 2) * x ** 2 + F(1, 2) * x


def comb(frequency, phase=F(0)):
    phase = frac(phase)
    ans = []
    for k in range(-2, frequency + 3):
        lo = max(F(0), (F(k) + phase - THETA) / frequency)
        hi = min(F(1), (F(k) + phase + THETA) / frequency)
        if lo < hi:
            ans.append((lo, hi))
    ans.sort()
    return ans


def centered_overlap(left, right):
    i = j = 0
    total = F(0)
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            total += ((hi - F(1, 2)) ** 2
                      - (lo - F(1, 2)) ** 2) / 2
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
    # Qe-Pf=5*1-3*2=-1 is load-bearing.
    cross = -1
    r, s = cell % L0, 2 * cell % L0
    d = 5 * r - 3 * s
    u, v = F(d + cross, L0), F(-d, L0)
    a, b = F(3, 14), F(5, 14)
    psi = (b3bar(u + a - b) + b3bar(u - a + b)
           + b3bar(v + a - b) + b3bar(v - a + b)
           - b3bar(u + a + b) - b3bar(u - a - b)
           - b3bar(v + a + b) - b3bar(v - a - b))
    return F(1, 49) + F(28, 15 * cross) * psi


def kappa(cell):
    return 2 if cell <= 5 or cell >= 162 else 0


owner_q = tuple(
    (F(kappa(cell)) - 2 * limit(cell) + 1848 * correction(cell)) / 423360
    for cell in range(L0)
)
owner_digest = sha256("\n".join(map(str, owner_q)).encode("ascii")).hexdigest()
require(owner_digest ==
        "b53d77a69f39a5f8c893b4cdceaaeecdd0dc70a16be02e6d28f3f6d7b520feef",
        "owner word digest")
require(sum(owner_q) == F(1, 24696), "owner Hodge sum")
require((sum(value > 0 for value in owner_q),
         sum(value < 0 for value in owner_q),
         sum(value == 0 for value in owner_q)) == (156, 12, 0),
        "owner sign census")
require(min(value for value in owner_q if value > 0) == F(1, 6174000),
        "positive owner minimum")
require(max(owner_q) == F(751, 666792000), "positive owner maximum")

scale = 1
for value in owner_q:
    scale = lcm(scale, value.denominator)
scaled_q = tuple(int(value * scale) for value in owner_q)
hodge_completion = -sum(owner_q)
scaled_completion = int(hodge_completion * scale)
require((scale, sum(scaled_q), scaled_completion)
        == (32006016000, 1296000, -1296000), "primitive scaling")
require(reduce(gcd, tuple(abs(value) for value in scaled_q)
               + (abs(scaled_completion),)) == 1, "scaled word primitive")

# THM-3234's deterministic Singer plane F_13[u]/(u^2-2), alpha=1+2u.
def field_mul(x, y):
    a, b = x
    c, d = y
    return ((a * c + 2 * b * d) % 13, (a * d + b * c) % 13)


alpha = (1, 2)
points = []
point = (1, 0)
for _ in range(168):
    points.append(point)
    point = field_mul(point, alpha)
require(point == (1, 0) and len(set(points)) == 168,
        "Singer orbit")


def determinant_mod(matrix, prime):
    a = [[entry % prime for entry in row] for row in matrix]
    determinant = 1
    for column in range(len(a)):
        pivot_row = next(
            (row for row in range(column, len(a)) if a[row][column]), None)
        if pivot_row is None:
            return 0
        if pivot_row != column:
            a[column], a[pivot_row] = a[pivot_row], a[column]
            determinant = -determinant
        pivot = a[column][column]
        determinant = determinant * pivot % prime
        inverse = pow(pivot, prime - 2, prime)
        for row in range(column + 1, len(a)):
            multiplier = a[row][column] * inverse % prime
            for index in range(column + 1, len(a)):
                a[row][index] = (
                    a[row][index] - multiplier * a[column][index]) % prime
    return determinant % prime


def determinant_bareiss(matrix):
    a = [row[:] for row in matrix]
    previous = 1
    sign = 1
    for column in range(len(a) - 1):
        pivot_row = next(
            row for row in range(column, len(a)) if a[row][column])
        if pivot_row != column:
            a[column], a[pivot_row] = a[pivot_row], a[column]
            sign = -sign
        pivot = a[column][column]
        for row in range(column + 1, len(a)):
            for index in range(column + 1, len(a)):
                numerator = (a[row][index] * pivot
                             - a[row][column] * a[column][index])
                require(numerator % previous == 0,
                        ("Bareiss divisibility", column, row, index))
                a[row][index] = numerator // previous
        previous = pivot
        for row in range(column + 1, len(a)):
            a[row][column] = 0
    return sign * a[-1][-1]


def gauge_matrix(multiplier, shift, completion):
    matrix = [[0] * 13 for _ in range(13)]
    matrix[0][0] = completion
    for owner, value in enumerate(scaled_q):
        x, y = points[(shift + multiplier * owner) % 168]
        matrix[x][y] = value
    return matrix


prime = 1000000007
units = tuple(a for a in range(168) if gcd(a, 168) == 1)
require(len(units) == 48 and scale % prime != 0, "gauge/modulus control")


def gauge_census(completion):
    determinants = tuple(
        determinant_mod(gauge_matrix(multiplier, shift, completion), prime)
        for multiplier in units
        for shift in range(168)
    )
    frequencies = Counter(determinants)
    require(len(determinants) == 8064, "Singer gauge census")
    require(all(determinants), "singular Singer gauge")
    require(len(frequencies) == 168
            and set(frequencies.values()) == {48}, "determinant orbit census")
    ordered_digest = sha256(
        "\n".join(map(str, determinants)).encode("ascii")).hexdigest()
    frequency_digest = sha256(
        "\n".join("%d:%d" % (key, frequencies[key])
                  for key in sorted(frequencies)).encode("ascii")).hexdigest()
    product = reduce(lambda x, y: x * y % prime, determinants, 1)
    return (len(frequencies), product, ordered_digest, frequency_digest,
            min(frequencies), max(frequencies))


zero_census = gauge_census(0)
hodge_census = gauge_census(scaled_completion)
require(zero_census == (
    168, 463433994,
    "e8391c1b3f0a4584f6df15953aa89ed49a8b7dedaf8707ccbfa3689b8fa229e2",
    "7e15bdc29de5647c3aae0a3e36125c3944e1fe38b58f038c3bb0ed45395e7c54",
    540717, 999459290), "zero-completion census")
require(hodge_census == (
    168, 565425212,
    "c4a504a86966c5dc7c0d375e962ae2419fe99271c79b8ae7ddc2f25a5d82d730",
    "bf31ba1bebb4664b6b929237ec57d328ee0b901595f61c6d02e111dc2d5ce85a",
    6636781, 993363226), "Hodge-completion census")

identity_zero = determinant_bareiss(gauge_matrix(1, 0, 0))
identity_hodge = determinant_bareiss(
    gauge_matrix(1, 0, scaled_completion))
require(identity_zero ==
        1306489535376679896012084927085361464824505576544950660174960053,
        "exact identity zero determinant")
require(identity_hodge ==
        -41129257721095723275557032050508281832453443425542278710444911947,
        "exact identity Hodge determinant")
require((identity_zero % prime, identity_hodge % prime)
        == (825265491, 852361250), "exact/modular determinant agreement")

# The rational central contrast has zero neutral Fourier mode and coefficient
# 13 in each of the twelve charged cyclotomic blocks.
contrast = (12,) + (-1,) * 12
require(sum(contrast) == 0, "central contrast mean")
charged_reductions = []
for kappa_value in range(1, 13):
    coefficients = [0] * 13
    for delta, value in enumerate(contrast):
        coefficients[(-kappa_value * delta) % 13] += value
    reduced = tuple(coefficients[index] - coefficients[12]
                    for index in range(12))
    charged_reductions.append(reduced)
require(set(charged_reductions) == {(13,) + (0,) * 11},
        "charged contrast transform")

print("THM-3252 SINGER-COMPACTIFIED OWNER CYCLICITY EXACT AUDIT")
print("dependency_hash_checks=%d" % len(DEPENDENCIES))
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("owner_digest=%s,signs=(156,12,0),sum=1/24696" % owner_digest)
print("scale=%d,scaled_sum=1296000,hodge_completion=-1296000" % scale)
print("singer_plane=F13[u]/(u^2-2),alpha=(1,2),orbit=168")
print("gauge_count=8064,prime=%d,all_determinants_nonzero=PASS" % prime)
print("zero_completion_census=%s" % (zero_census,))
print("hodge_completion_census=%s" % (hodge_census,))
print("identity_zero_determinant=%d" % identity_zero)
print("identity_hodge_determinant=%d" % identity_hodge)
print("central_contrast=(12,-1^12),neutral=0,charged_blocks=13")
print("charged_orbit_span_dimension=12*169=2028")
print("scope=signed-second-corrector-linear-frame-not-positive-physical-current")
print("all_exact_checks=PASS")
