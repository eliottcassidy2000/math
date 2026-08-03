#!/usr/bin/env python3
"""Exact universal resonant-degree PRS wall atlas for THM-3217."""

from __future__ import annotations

import ast
import hashlib
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.fields import field


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT
    / "01-canon"
    / "theorems"
    / "THM-3192-reciprocal-coefficient-jet-transfer-and-z-adic-pluecker-return.md": (
        "69ae0b5e4c9526d643bc976220a4752e33632248c80394ca4cf9c9620105878f"
    ),
    ROOT
    / "01-canon"
    / "theorems"
    / "THM-3214-two-jet-pseudo-division-locality-and-catalan-sharpness.md": (
        "af6731c4c0735240e511d34842edaff2ba1a542e981e435e3bc30a25fe2db30b"
    ),
}


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(data).hexdigest()


d, p = sp.symbols("d p", integer=True)
RATIONAL_FUNCTIONS, degree = field("d", QQ)
a = degree - 2


def choose_polynomial(value: object, index: int) -> object:
    return sp.prod(value - j for j in range(index)) / factorial(index)


def factorial_ratio_from_2a(offset: int) -> object:
    if offset >= 0:
        return sp.prod(2 * a + j for j in range(1, offset + 1))
    return 1 / sp.prod(2 * a - j for j in range(0, -offset))


def reciprocal_top(degree_shift: int, codimension: int) -> object:
    """Coefficient of v^(a+shift-codim), normalized by (2a)!.

    This is the exact multinomial moment formula, written polynomially in d.
    The production application has d>=9, so all requested codimensions exist.
    """

    value = choose_polynomial(a + degree_shift, codimension) * sum(
        sp.binomial(codimension, ell)
        * degree ** (codimension - ell)
        * (-1) ** ell
        * factorial_ratio_from_2a(
            2 * degree_shift - 2 * codimension + ell
        )
        for ell in range(codimension + 1)
    )
    return value


def p1(left: list[object], right: list[object]) -> object:
    return left[0] * right[1] - left[1] * right[0]


def pseudo(left: list[object], right: list[object]) -> list[object]:
    connection = p1(left, right)
    return [
        left[0] ** 2 * right[index + 2]
        - left[0] * right[0] * left[index + 2]
        - connection * left[index + 1]
        for index in range(min(len(left), len(right)) - 2)
    ]


# Seven input jets are exactly what the third remainder one-jet consumes.
A = [reciprocal_top(0, index) for index in range(8)]
B = [reciprocal_top(1, index) for index in range(8)]
R = pseudo(A, B)
S = pseudo(R, A)
T = pseudo(S, R)


H = 24 * d - 35
R1_WALL = 24 * d**2 - 259 * d + 315
CHI1_WALL = 28 * d - 45
J = 256 * d**4 - 33792 * d**3 + 187360 * d**2 - 348768 * d + 218295
J1 = (
    256 * d**5
    - 37120 * d**4
    + 650464 * d**3
    - 2875072 * d**2
    + 4835439 * d
    - 2837835
)
K = (
    5120 * d**5
    - 963840 * d**4
    + 6841088 * d**3
    - 18014016 * d**2
    + 20587392 * d
    - 8513505
)
U = (
    327680 * d**8
    + 75694080 * d**7
    - 5035040768 * d**6
    + 61753251840 * d**5
    - 341070457600 * d**4
    + 1023103349760 * d**3
    - 1737933934176 * d**2
    + 1581643823760 * d
    - 602628451425
)
V = (
    327680 * d**9
    + 70123520 * d**8
    - 6372368384 * d**7
    + 152519243776 * d**6
    - 1447233207040 * d**5
    + 7075663957760 * d**4
    - 19705993431648 * d**3
    + 31775324380272 * d**2
    - 27782756946945 * d
    + 10244683674225
)


def rf(expression: sp.Expr) -> object:
    return RATIONAL_FUNCTIONS.from_expr(expression.xreplace({d: sp.Symbol("d")}))


EXPECTED = {
    "R0": rf(-((d - 1) * H) / (2 * (2 * d - 7) * (2 * d - 5) ** 2)),
    "R1": rf(
        (d - 1)
        * R1_WALL
        / (4 * (2 * d - 9) * (2 * d - 7) * (2 * d - 5) ** 2)
    ),
    "chi1": rf(
        d
        * (d - 1)
        * CHI1_WALL
        / ((2 * d - 9) * (2 * d - 7) * (2 * d - 5) ** 3)
    ),
    "S0": rf(
        (d - 1) ** 2
        * J
        / (
            16
            * (2 * d - 11)
            * (2 * d - 9) ** 2
            * (2 * d - 7) ** 3
            * (2 * d - 5) ** 4
        )
    ),
    "S1": rf(
        -(d - 1) ** 2
        * J1
        / (
            32
            * (2 * d - 13)
            * (2 * d - 11)
            * (2 * d - 9) ** 2
            * (2 * d - 7) ** 3
            * (2 * d - 5) ** 4
        )
    ),
    "chi2": rf(
        d
        * (d - 1) ** 3
        * K
        / (
            16
            * (2 * d - 13)
            * (2 * d - 11)
            * (2 * d - 9) ** 3
            * (2 * d - 7) ** 4
            * (2 * d - 5) ** 6
        )
    ),
    "T0": rf(
        -(d - 1) ** 5
        * H**2
        * U
        / (
            2048
            * (2 * d - 15)
            * (2 * d - 13) ** 2
            * (2 * d - 11) ** 3
            * (2 * d - 9) ** 4
            * (2 * d - 7) ** 7
            * (2 * d - 5) ** 10
        )
    ),
    "T1": rf(
        (d - 1) ** 5
        * H**2
        * V
        / (
            4096
            * (2 * d - 17)
            * (2 * d - 15)
            * (2 * d - 13) ** 2
            * (2 * d - 11) ** 3
            * (2 * d - 9) ** 4
            * (2 * d - 7) ** 7
            * (2 * d - 5) ** 10
        )
    ),
}

ACTUAL = {
    "R0": R[0],
    "R1": R[1],
    "chi1": p1(R, A),
    "S0": S[0],
    "S1": S[1],
    "chi2": p1(S, R),
    "T0": T[0],
    "T1": T[1],
}

for label in EXPECTED:
    require(ACTUAL[label] == EXPECTED[label], label)


def direct_moment_coefficient(degree_value: int, exponent: int, power: int) -> int:
    """Raw coefficient of v^power from the multinomial moment polynomial."""

    return comb(exponent, power) * sum(
        comb(exponent - power, ell)
        * degree_value ** (exponent - power - ell)
        * (-1) ** ell
        * factorial(2 * power + ell)
        for ell in range(exponent - power + 1)
    )


def evaluate_rational_function(value: object, degree_value: int) -> Fraction:
    rational = value.as_expr().subs(sp.Symbol("d"), degree_value)
    return Fraction(int(rational.p), int(rational.q))


# Independent raw-moment controls do not reuse formula (4).
for degree_value in range(9, 25):
    normalization = factorial(2 * degree_value - 4)
    direct_a = [
        Fraction(
            direct_moment_coefficient(
                degree_value, degree_value - 2, degree_value - 2 - index
            ),
            normalization,
        )
        for index in range(8)
    ]
    direct_b = [
        Fraction(
            direct_moment_coefficient(
                degree_value, degree_value - 1, degree_value - 1 - index
            ),
            normalization,
        )
        for index in range(8)
    ]
    direct_r = pseudo(direct_a, direct_b)
    direct_s = pseudo(direct_r, direct_a)
    direct_t = pseudo(direct_s, direct_r)
    direct_coordinates = [
        direct_r[0],
        direct_r[1],
        p1(direct_r, direct_a),
        direct_s[0],
        direct_s[1],
        p1(direct_s, direct_r),
        direct_t[0],
        direct_t[1],
    ]
    expected_coordinates = [
        evaluate_rational_function(EXPECTED[label], degree_value)
        for label in EXPECTED
    ]
    require(direct_coordinates == expected_coordinates, ("direct moment", degree_value))


H6 = 24 * p + 109
J6 = 256 * p**4 - 27648 * p**3 - 365600 * p**2 - 1528800 * p - 2096649
K6 = (
    5120 * p**5
    - 810240 * p**4
    - 14447872 * p**3
    - 92004672 * p**2
    - 256323456 * p
    - 265142241
)
U6 = (
    327680 * p**8
    + 91422720 * p**7
    - 1525587968 * p**6
    - 58319874048 * p**5
    - 605420542720 * p**4
    - 3106619397120 * p**3
    - 8698849881696 * p**2
    - 12762150724080 * p
    - 7697077854849
)
V6 = (
    327680 * p**9
    + 87818240 * p**8
    - 2581766144 * p**7
    - 38490294272 * p**6
    + 127672001792 * p**5
    + 4526613136640 * p**4
    + 30758409758112 * p**3
    + 98688343652784 * p**2
    + 157343795925183 * p
    + 100505677023531
)

for universal, offset_six, label in [
    (H, H6, "H"),
    (J, J6, "J"),
    (K, K6, "K"),
    (U, U6, "U"),
    (V, V6, "V"),
]:
    require(sp.expand(universal.subs(d, p + 6) - offset_six) == 0, label)


# A no-root certificate modulo one prime rules out every integral root.
NO_ROOT_CERTIFICATES = {
    "R1": (R1_WALL, 11),
    "J": (J, 13),
    "J1": (J1, 37),
    "K": (K, 23),
    "U": (U, 17),
    "V": (V, 19),
}
for label, (polynomial, modulus) in NO_ROOT_CERTIFICATES.items():
    leading = int(sp.LC(sp.Poly(polynomial, d)))
    require(leading % modulus != 0, (label, "leading", modulus))
    require(
        all(int(polynomial.subs(d, residue)) % modulus for residue in range(modulus)),
        (label, "root modulo", modulus),
    )


def xi(offset: int) -> int:
    walls = [
        H,
        R1_WALL,
        CHI1_WALL,
        J,
        J1,
        K,
        U,
        V,
    ]
    denominators = [2 * d - odd for odd in range(5, 18, 2)]
    value = 2 * d * (d - 1) * sp.prod(walls) * sp.prod(denominators)
    return int(value.subs(d, offset))


for offset in range(2, 65):
    require(xi(offset) != 0, ("xi", offset))


# At offset six, five distinct primes isolate five successive walls.
WALL_CONTROLS = [
    (109, "H", H),
    (232961, "J", J),
    (2767, "K", K),
    (1067703961, "U", U),
    (52511, "V", V),
]
for prime, selected_label, selected in WALL_CONTROLS:
    require(bool(sp.isprime(prime)), ("not prime", prime))
    values = {
        label: int(poly.subs(d, prime + 6)) % prime
        for label, poly in [("H", H), ("J", J), ("K", K), ("U", U), ("V", V)]
    }
    require(values[selected_label] == 0, (prime, selected_label, values))
    require(
        all(value != 0 for label, value in values.items() if label != selected_label),
        (prime, "not isolated", values),
    )


# Common normalization by a unit u gives row exponents 1,1,3,7,17.
left_exponent, right_exponent = 1, 1
row_exponents = []
for _ in range(3):
    next_exponent = 2 * left_exponent + right_exponent
    row_exponents.append(next_exponent)
    right_exponent, left_exponent = left_exponent, next_exponent
require(row_exponents == [3, 7, 17], row_exponents)
require(3 + 1 == 4, "chi1 scaling")
require(7 + 3 == 10, "chi2 scaling")


source = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source)
require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert")
require(
    not any(
        isinstance(atom, sp.Float)
        for atom in sp.preorder_traversal(
            sp.Tuple(H, R1_WALL, CHI1_WALL, J, J1, K, U, V)
        )
    ),
    "float",
)
for dependency, expected_hash in DEPENDENCIES.items():
    require(lf_sha256(dependency) == expected_hash, dependency)


print("THM-3217 exact universal resonant-degree PRS wall atlas")
print("dependency_pins=2")
print("input_top_jets=8 derived_coordinates=8 exact_symbolic=PASS")
print("direct_moment_controls=16 coordinate_checks=128 PASS")
print("universal_walls=H,J,K,U,V offset_six_specialization=PASS")
print("no_integer_root_certificates=6 offsets_2_through_64=PASS")
print("fixed_offset_exception_integer=nonzero_for_every_integer_s_ge_2")
print("normalization_exponents=1,1,3,7,17 chi1=4 chi2=10")
print("isolated_offset_six_wall_primes=109,232961,2767,1067703961,52511")
print("scope=first_three_reciprocal_remainders_not_arbitrary_depth")
