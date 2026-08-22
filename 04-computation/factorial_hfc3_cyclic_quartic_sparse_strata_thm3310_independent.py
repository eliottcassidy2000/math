"""Exact sparse-stratum probe for THM-3304's cyclic quartic cell."""

from __future__ import annotations

import argparse
import ast
import hashlib
from fractions import Fraction
from itertools import combinations
from math import comb, factorial
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path)
args = parser.parse_args()
output_lines = []


# Restricted Fourier monomials (a exponent, b exponent) in THM-3304 order.
BASIS = ((0, 1), (2, 0), (1, 2), (0, 4), (3, 1))
NAMES = ("b", "a2", "ab2", "b4", "a3b")


def multinomial3(i, j, ell):
    return factorial(i + j + ell) // (factorial(i) * factorial(j)
                                      * factorial(ell))


def kernel_coefficient(r, s):
    value = 0
    for ell in range(min(r, s) + 1):
        if (r - ell) % 3 or (s - ell) % 3:
            continue
        value += multinomial3((r - ell) // 3, (s - ell) // 3, ell) * 3**ell
    return value


def moment_ab(r, s):
    coefficient = kernel_coefficient(r, s)
    return Fraction(2 * factorial(r) * factorial(s) * coefficient,
                    factorial(r + s + 2))


def pair_moment(left, right, exponent, t):
    expression = 0
    for j in range(exponent + 1):
        r = (exponent - j) * BASIS[left][0] + j * BASIS[right][0]
        s = (exponent - j) * BASIS[left][1] + j * BASIS[right][1]
        expression += sp.Rational(comb(exponent, j)) * sp.Rational(
            moment_ab(r, s).numerator, moment_ab(r, s).denominator) * t**j
    return sp.Poly(expression, t, domain=sp.QQ)


def weak_compositions(total, length, prefix=()):
    if length == 1:
        yield prefix + (total,)
        return
    for value in range(total + 1):
        yield from weak_compositions(total - value, length - 1,
                                     prefix + (value,))


def multinomial(parts):
    answer = factorial(sum(parts))
    for part in parts:
        answer //= factorial(part)
    return answer


def triple_moment(indices, exponent, y, z):
    expression = 0
    for counts in weak_compositions(exponent, 3):
        r = sum(count * BASIS[index][0]
                for count, index in zip(counts, indices))
        s = sum(count * BASIS[index][1]
                for count, index in zip(counts, indices))
        value = moment_ab(r, s)
        expression += (sp.Rational(multinomial(counts) * value.numerator,
                                   value.denominator)
                       * y**counts[1] * z**counts[2])
    return sp.Poly(expression, y, z, domain=sp.QQ)


t = sp.symbols("t")
rows = []
for left in range(len(BASIS)):
    for right in range(left + 1, len(BASIS)):
        m3 = pair_moment(left, right, 3, t)
        m6 = pair_moment(left, right, 6, t)
        gcd36 = sp.gcd(m3, m6)
        require(m3.degree() == 3 and m6.degree() == 6,
                "projective endpoint degree retained")
        require(m3.nth(0) != 0 and m3.LC() != 0,
                "both pure cubic endpoints are nonzero")
        rows.append((NAMES[left], NAMES[right], gcd36.degree(),
                     sp.factor(sp.resultant(m3.as_expr(), m6.as_expr(), t))))

require(len(rows) == 10, "all ten coordinate lines")
require(all(row[2] == 0 and row[3] != 0 for row in rows),
        "every coordinate-line resultant is nonzero")
output_lines.append(
    "universe=ten projective coordinate lines in THM-3304 cyclic quartic basis")
for row in rows:
    output_lines.append("pair=" + row[0] + "," + row[1]
                        + "; gcd_M3_M6_degree=" + str(row[2])
                        + "; resultant=" + str(row[3]))
output_lines.append(
    "nonzero_resultants=" + str(sum(row[3] != 0 for row in rows)) + "/10")

# Each coordinate P2 is covered by the affine chart c_first=1 plus the
# already excluded line c_first=0.  Test the first three surviving moments.
y, z = sp.symbols("y z")
triple_rows = []
for indices in combinations(range(len(BASIS)), 3):
    forms = [triple_moment(indices, exponent, y, z)
             for exponent in (3, 6, 9)]
    basis = sp.groebner([form.as_expr() for form in forms], y, z,
                        order="grevlex", domain=sp.QQ)
    contains_one = any(poly.as_expr().is_number and poly.as_expr() != 0
                       for poly in basis.polys)
    triple_rows.append((tuple(NAMES[index] for index in indices),
                        contains_one, len(basis.polys)))
    output_lines.append("triple=" + ",".join(triple_rows[-1][0])
                        + "; ideal_M3_M6_M9_is_unit=" + str(contains_one)
                        + "; groebner_size=" + str(len(basis.polys)))

require(len(triple_rows) == 10 and all(row[1] for row in triple_rows),
        "every coordinate-plane affine ideal is the unit ideal")
syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
require(not any(isinstance(node, ast.Assert) for node in ast.walk(syntax)),
        "no assertion-dependent truth gate")
require(not any(isinstance(node, ast.Constant) and isinstance(node.value, float)
                for node in ast.walk(syntax)),
        "no floating literals")
output_lines.extend([
    "unit_triple_charts="
    + str(sum(row[1] for row in triple_rows)) + "/10",
    "projective_coverage=affine first-coordinate chart plus excluded coordinate line",
    "consequence=no nonzero cyclic quartic of Fourier-basis support at most three kills M3,M6,M9",
    "scope=support four and five remain open; no HFC(3) or FC(3) conclusion",
    "controls=all pair resultants; all triple unit ideals; projective boundaries; no assert or float nodes",
])
payload = "\n".join(output_lines) + "\n"
payload += "payload_sha256=" + hashlib.sha256(payload.encode("ascii")).hexdigest() + "\n"
if args.output is not None:
    with args.output.open("w", encoding="ascii", newline="\n") as handle:
        handle.write(payload)
print(payload, end="")
