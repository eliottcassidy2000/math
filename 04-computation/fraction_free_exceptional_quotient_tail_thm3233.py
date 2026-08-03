#!/usr/bin/env python3
"""Exact controls for THM-3233's exceptional quotient-tail renormalization."""

from __future__ import annotations

import ast
import hashlib
from fractions import Fraction
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = (
    ROOT
    / "01-canon"
    / "theorems"
    / "THM-3231-fraction-free-simple-pivot-wall-second-order-blowup-carry.md"
)
EXPECTED_DEPENDENCY_SHA256 = (
    "95f193b2f1569ded324928ec4400f5e4555dd38cb925817c2d5394c3a4224174"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def trim(series: list[object]) -> list[object]:
    return list(series)


def multiply(left: list[object], right: list[object], length: int) -> list[object]:
    out = [0] * length
    for i, left_i in enumerate(left):
        for j, right_j in enumerate(right):
            if i + j < length:
                out[i + j] += left_i * right_j
    return [sp.expand(value) for value in out]


def pseudo(left: list[object], right: list[object], length: int | None = None) -> list[object]:
    available = min(len(left), len(right)) - 2
    if length is None:
        length = available
    require(0 <= length <= available, ("pseudo length", length, available))
    connection = left[0] * right[1] - left[1] * right[0]
    return [
        sp.expand(
            left[0] ** 2 * right[index + 2]
            - left[0] * right[0] * left[index + 2]
            - connection * left[index + 1]
        )
        for index in range(length)
    ]


def divide_series(numerator: list[object], denominator: list[object], length: int) -> list[object]:
    require(denominator[0] != 0, "series denominator")
    quotient: list[object] = []
    for index in range(length):
        value = numerator[index] if index < len(numerator) else 0
        value -= sum(
            denominator[j] * quotient[index - j]
            for j in range(1, min(index, len(denominator) - 1) + 1)
        )
        quotient.append(sp.cancel(value / denominator[0]))
    return quotient


# Scalar homogeneity and common-factor covariance, symbolically.
LENGTH = 5
c, d = sp.symbols("c d")
h = list(sp.symbols("h0:9"))
f = list(sp.symbols("f0:9"))
g = list(sp.symbols("g0:9"))
scaled = pseudo([c * x for x in f], [d * x for x in g], LENGTH)
base = pseudo(f, g, LENGTH)
for index in range(LENGTH):
    require(sp.expand(scaled[index] - c**2 * d * base[index]) == 0,
            ("scalar homogeneity", index))

hf = multiply(h, f, 9)
hg = multiply(h, g, 9)
common_left = pseudo(hf, hg, LENGTH)
common_right = multiply(h, base, LENGTH)
for index in range(LENGTH):
    require(
        sp.expand(common_left[index] - h[0] ** 2 * common_right[index]) == 0,
        ("common factor covariance", index),
    )


# Recheck the exceptional base and the first normalized successor.
a = sp.symbols("a")
v = list(sp.symbols("v0:11"))
t = list(sp.symbols("t0:12"))
w = [a] + v
x = pseudo(w, t, 8)
y = pseudo(x, w, 6)
q = divide_series(t, v, 10)
r = q[3:]
base_controls = 0
for index, coefficient in enumerate(y[:5]):
    polynomial = sp.Poly(coefficient, a)
    require(polynomial.nth(0) == polynomial.nth(1) == 0,
            ("exceptional square", index))
    expected = -t[0] * v[0] ** 3 * sum(
        v[j] * r[index - j] for j in range(index + 1)
    )
    require(sp.cancel(polynomial.nth(2) - expected) == 0,
            ("exceptional quotient tail", index))
    base_controls += 1

# R_7/a^4 on the exceptional divisor is the first new-tail row.
first_tail = pseudo(r, [1] + [0] * (len(r) - 1), 3)
expected_kappa = t[0] ** 3 * v[0] ** 9
next_row = pseudo(y, x, 3)
for index, coefficient in enumerate(next_row):
    polynomial = sp.Poly(coefficient, a)
    actual = polynomial.nth(4)
    expected = expected_kappa * sum(
        v[j] * first_tail[index - j] for j in range(index + 1)
    )
    require(sp.cancel(actual - expected) == 0,
            ("first recursive exceptional row", index))


# Twice-Pell normal exponents and the scalar gauge.
normal_exponents = [0, 2]
for _ in range(10):
    normal_exponents.append(2 * normal_exponents[-1] + normal_exponents[-2])
pell = [0, 1]
for _ in range(10):
    pell.append(2 * pell[-1] + pell[-2])
require(normal_exponents == [2 * value for value in pell],
        ("twice Pell", normal_exponents, pell))

kappas = [t[0] * v[0], -t[0] * v[0] ** 3]
for _ in range(5):
    kappas.append(sp.expand(kappas[-1] ** 2 * kappas[-2] * v[0] ** 2))
require(kappas[2] == t[0] ** 3 * v[0] ** 9, "first scalar gauge")


# Exact boundary controls for E(r,1).
def fraction_pseudo(left: list[Fraction], right: list[Fraction], length: int) -> list[Fraction]:
    connection = left[0] * right[1] - left[1] * right[0]
    return [
        left[0] ** 2 * right[index + 2]
        - left[0] * right[0] * left[index + 2]
        - connection * left[index + 1]
        for index in range(length)
    ]


one = [Fraction(1)] + [Fraction(0)] * 8
lam = Fraction(2)
geometric = [lam**index for index in range(9)]
delayed_simple = [Fraction(0)] + geometric[:-1]
delayed_double = [Fraction(0), Fraction(0), Fraction(1)] + [Fraction(0)] * 6
generic = [Fraction(1), Fraction(2), Fraction(5), Fraction(7), Fraction(11),
           Fraction(13), Fraction(17), Fraction(19), Fraction(23)]

boundary_rows = {
    "zero": fraction_pseudo([Fraction(0)] * 9, one, 7),
    "geometric": fraction_pseudo(geometric, one, 7),
    "delayed_simple": fraction_pseudo(delayed_simple, one, 7),
    "delayed_double": fraction_pseudo(delayed_double, one, 7),
    "generic": fraction_pseudo(generic, one, 7),
}
require(all(value == 0 for value in boundary_rows["zero"]),
        "quadratic quotient zero tail")
require(all(value == 0 for value in boundary_rows["geometric"]),
        "geometric terminal")
require(boundary_rows["delayed_simple"][0] == 1,
        "delayed simple return")
require(all(value == 0 for value in boundary_rows["delayed_double"]),
        "delayed double collapse")
require(any(value != 0 for value in boundary_rows["generic"]),
        "generic survivor")


source = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source)
require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "assert statement")
require(
    not any(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    ),
    "float literal",
)
require(lf_sha256(DEPENDENCY) == EXPECTED_DEPENDENCY_SHA256, "dependency drift")


print("THM-3233 exact exceptional quotient-tail renormalization")
print("dependency_pin=THM-3231")
print("scalar_homogeneity_symbolic_coefficients=5 PASS")
print("common_factor_covariance_symbolic_coefficients=5 PASS")
print("exceptional_base_quotient_tail_coefficients=5 PASS")
print("first_recursive_exceptional_row_coefficients=3 PASS")
print("normal_exponents=0,2,4,10,24,58,140 twice_Pell=PASS")
print("scalar_gauge_initial=t0*v0,-t0*v0^3 next=t0^3*v0^9 PASS")
print("tail_boundaries=quadratic_zero,geometric_terminal,delayed_simple_return,delayed_double_zero PASS")
print("scope=exceptional_special_fibre_not_generic_root_or_physical_selector")
print("all_exact_checks=PASS")
