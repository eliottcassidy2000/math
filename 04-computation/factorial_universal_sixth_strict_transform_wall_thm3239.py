#!/usr/bin/env python3
"""Exact universal sixth strict-transform wall controls for THM-3239."""

from __future__ import annotations

import ast
import contextlib
import hashlib
import io
import runpy
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
THM3223 = (
    ROOT
    / "01-canon/theorems/THM-3223-fourth-fifth-resonant-prs-primitive-walls-pell-content-and-pivot-resurrection.md"
)
THM3231 = (
    ROOT
    / "01-canon/theorems/THM-3231-fraction-free-simple-pivot-wall-second-order-blowup-carry.md"
)
SCRIPT3223 = (
    ROOT
    / "04-computation/factorial_fourth_fifth_prs_primitive_walls_thm3223.py"
)
EXPECTED_HASHES = {
    THM3223: "2db785e5169b1c8b0eb414b4070765b6b76f47c0d02e92680eb206139c54dfb3",
    THM3231: "95f193b2f1569ded324928ec4400f5e4555dd38cb925817c2d5394c3a4224174",
    SCRIPT3223: "2afe51afe4719922b739e4aa0ea43fc285a130864f6db1d1ba281430d4e93fdc",
}


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(data).hexdigest()


for path, wanted_hash in EXPECTED_HASHES.items():
    require(lf_sha256(path) == wanted_hash, ("dependency drift", path))

captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    namespace = runpy.run_path(SCRIPT3223)
require("scope=selected_pivots_not_full_row_or_arbitrary_depth" in captured.getvalue(), "dependency output")

reciprocal_top = namespace["reciprocal_top"]
pseudo = namespace["pseudo"]
rf = namespace["rf"]
d = namespace["d"]
H = namespace["H"]
J = namespace["J"]
U = namespace["U"]
W13 = namespace["W13"]

W28 = (
    169920813940688814080000 * d**28
    + 331141682207614360879104000 * d**27
    - 134162410011258685536379863040 * d**26
    + 5300802398351890701896525021184 * d**25
    - 684739297310367200983522643804160 * d**24
    + 130331904744222160513651684933632000 * d**23
    - 12943819225240695078283631960879267840 * d**22
    + 713250673128714558996347298753011515392 * d**21
    - 24834141237082462639059202312165067325440 * d**20
    + 598783632792278487019889843452862948966400 * d**19
    - 10626718017854951523648995520000467801210880 * d**18
    + 144696459098683056904745037387918164398964736 * d**17
    - 1556104532530719261165315051927889255805747200 * d**16
    + 13492893612749552869825791823397116965185126400 * d**15
    - 95734068097296356726622150832929787026834718720 * d**14
    + 561609982081860252702393234504629255248835248128 * d**13
    - 2743000014814555875946490646195207597676182896640 * d**12
    + 11199447544319326718253356230653228142356647116800 * d**11
    - 38279428000372095326340982875930990156643565568000 * d**10
    + 109404450823229543746513469845082648703787025301504 * d**9
    - 260447513492258947342576045172684040607636400701440 * d**8
    + 512815846815398540657575958180230914285371601715200 * d**7
    - 825986183003965800574110651732582518260652206080000 * d**6
    + 1070580568879390678018966410506355769337873203200000 * d**5
    - 1089659124692158529207880465753484653373043380800000 * d**4
    + 838910469802753648130728862980511764057677888000000 * d**3
    - 459389089498986441643744053260439590725236765000000 * d**2
    + 159486965048208251611930715459002067591134125000000 * d
    - 26393359344902929074282529688291558301498193359375
)

DQ = (
    2**112
    * (2 * d - 27)
    * (2 * d - 25) ** 2
    * (2 * d - 23) ** 3
    * (2 * d - 21) ** 4
    * (2 * d - 19) ** 5
    * (2 * d - 17) ** 6
    * (2 * d - 15) ** 11
    * (2 * d - 13) ** 16
    * (2 * d - 11) ** 27
    * (2 * d - 9) ** 38
    * (2 * d - 7) ** 65
    * (2 * d - 5) ** 92
)
D6 = (
    2**168
    * (2 * d - 27)
    * (2 * d - 25) ** 2
    * (2 * d - 23) ** 3
    * (2 * d - 21) ** 4
    * (2 * d - 19) ** 7
    * (2 * d - 17) ** 10
    * (2 * d - 15) ** 17
    * (2 * d - 13) ** 24
    * (2 * d - 11) ** 41
    * (2 * d - 9) ** 58
    * (2 * d - 7) ** 99
    * (2 * d - 5) ** 140
)

EXPECTED_QUOTIENT = rf(-(d - 1) ** 46 * H**16 * J**6 * U**4 * W28 / DQ)
EXPECTED_RHO6 = rf(
    -(d - 1) ** 70 * H**24 * J**10 * U**4 * W13**2 * W28 / D6
)

A = [reciprocal_top(0, index) for index in range(13)]
B = [reciprocal_top(1, index) for index in range(13)]
R = pseudo(A, B)
S = pseudo(R, A)
T = pseudo(S, R)
W = pseudo(T, S)
X = pseudo(W, T)
Y = pseudo(X, W)
require(tuple(map(len, (A, R, S, T, W, X, Y))) == (13, 11, 9, 7, 5, 3, 1), "lengths")
require(Y[0] / W[0] ** 2 == EXPECTED_QUOTIENT, "strict quotient")
require(Y[0] == EXPECTED_RHO6, "sixth pivot")


def direct_moment_coefficient(degree_value: int, exponent: int, power: int) -> int:
    return comb(exponent, power) * sum(
        comb(exponent - power, ell)
        * degree_value ** (exponent - power - ell)
        * (-1) ** ell
        * factorial(2 * power + ell)
        for ell in range(exponent - power + 1)
    )


def evaluate(value: object, degree_value: int) -> Fraction:
    rational = value.as_expr().subs(sp.Symbol("d"), degree_value)
    return Fraction(int(rational.p), int(rational.q))


for degree_value in range(14, 21):
    normalization = factorial(2 * degree_value - 4)
    direct_a = [
        Fraction(
            direct_moment_coefficient(
                degree_value, degree_value - 2, degree_value - 2 - index
            ),
            normalization,
        )
        for index in range(13)
    ]
    direct_b = [
        Fraction(
            direct_moment_coefficient(
                degree_value, degree_value - 1, degree_value - 1 - index
            ),
            normalization,
        )
        for index in range(13)
    ]
    direct_r = pseudo(direct_a, direct_b)
    direct_s = pseudo(direct_r, direct_a)
    direct_t = pseudo(direct_s, direct_r)
    direct_w = pseudo(direct_t, direct_s)
    direct_x = pseudo(direct_w, direct_t)
    direct_y = pseudo(direct_x, direct_w)
    require(
        direct_y[0] == evaluate(EXPECTED_RHO6, degree_value),
        ("direct sixth pivot", degree_value),
    )
    require(
        direct_y[0] / direct_w[0] ** 2
        == evaluate(EXPECTED_QUOTIENT, degree_value),
        ("direct strict quotient", degree_value),
    )


def polynomial_power_mod(
    base: sp.Poly, exponent: int, modulus_polynomial: sp.Poly
) -> sp.Poly:
    answer = sp.Poly(1, d, modulus=modulus_polynomial.get_modulus())
    factor = base.rem(modulus_polynomial)
    while exponent:
        if exponent & 1:
            answer = (answer * factor).rem(modulus_polynomial)
        factor = (factor * factor).rem(modulus_polynomial)
        exponent //= 2
    return answer


reduced = sp.Poly(W28, d, modulus=73)
variable = sp.Poly(d, d, modulus=73)
require(sp.Poly(W28, d).content() == 1, "W28 primitive")
require(reduced.degree() == 28, "W28 degree drop")
require(
    (polynomial_power_mod(variable, 73**28, reduced) - variable)
    .rem(reduced)
    .is_zero,
    "W28 Frobenius closure",
)
for exponent in (14, 4):
    test = polynomial_power_mod(variable, 73**exponent, reduced)
    require(sp.gcd(reduced, test - variable).degree() == 0, ("W28 gcd", exponent))
require(sp.Poly(W28, d, modulus=2) == sp.Poly(1, d, modulus=2), "W28 parity")

pell = [0, 1]
for _ in range(5):
    pell.append(2 * pell[-1] + pell[-2])
require(pell[1:] == [1, 2, 5, 12, 29, 70], pell)
gauge = [1, 1]
for _ in range(6):
    gauge.append(2 * gauge[-1] + gauge[-2])
require(gauge == [1, 1, 3, 7, 17, 41, 99, 239], gauge)

source = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source)
require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert")
require(
    not any(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    ),
    "float literal",
)

print("THM-3239 exact universal sixth strict-transform wall")
print("dependency_pins=THM-3223,THM-3231,THM-3223-script")
print("input_top_jets=13 row_lengths=13,11,9,7,5,3,1")
print("direct_moment_controls=7 quotient_and_full_checks=14 PASS")
print("strict_quotient=-(d-1)^46*H^16*J^6*U^4*W28/DQ")
print("full_rho6=-(d-1)^70*H^24*J^10*U^4*W13^2*W28/D6")
print("W28=primitive_degree28_irreducible_mod73_parity1")
print("DQ_2_power=112 exponents=1,2,3,4,5,6,11,16,27,38,65,92")
print("D6_2_power=168 exponents=1,2,3,4,7,10,17,24,41,58,99,140")
print("pell_pivot_exponents=1,2,5,12,29,70 gauge=1,1,3,7,17,41,99,239")
print("fixed_offset_extension=finite_exception_set_through_row6")
print("scope=selected_sixth_pivot_and_strict_transform_not_arbitrary_depth")
