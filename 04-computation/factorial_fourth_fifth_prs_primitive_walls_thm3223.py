#!/usr/bin/env python3
"""Exact fourth/fifth reciprocal PRS primitive-wall controls for THM-3223."""

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
DEPENDENCY = (
    ROOT
    / "01-canon"
    / "theorems"
    / "THM-3217-universal-resonant-degree-prs-wall-atlas-and-fixed-offset-exception-set.md"
)
EXPECTED_DEPENDENCY_SHA256 = (
    "2cdc04367a1ae8b86261912ee595aaf62423c32dd4e28d50933d8f9dc59221bb"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(data).hexdigest()


d = sp.symbols("d", integer=True)
RATIONAL_FUNCTIONS, degree = field("d", QQ)
a = degree - 2


def choose_polynomial(value: object, index: int) -> object:
    answer = RATIONAL_FUNCTIONS.one
    for j in range(index):
        answer *= value - j
    return answer / factorial(index)


def factorial_ratio_from_2a(offset: int) -> object:
    answer = RATIONAL_FUNCTIONS.one
    if offset >= 0:
        for j in range(1, offset + 1):
            answer *= 2 * a + j
        return answer
    for j in range(0, -offset):
        answer *= 2 * a - j
    return 1 / answer


def reciprocal_top(degree_shift: int, codimension: int) -> object:
    return choose_polynomial(a + degree_shift, codimension) * sum(
        comb(codimension, ell)
        * degree ** (codimension - ell)
        * (-1) ** ell
        * factorial_ratio_from_2a(
            2 * degree_shift - 2 * codimension + ell
        )
        for ell in range(codimension + 1)
    )


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


# Eleven input coefficients are exactly the Jet_10 budget of the fifth pivot.
A = [reciprocal_top(0, index) for index in range(11)]
B = [reciprocal_top(1, index) for index in range(11)]
R = pseudo(A, B)
S = pseudo(R, A)
T = pseudo(S, R)
W = pseudo(T, S)
X = pseudo(W, T)
require(tuple(map(len, (A, R, S, T, W, X))) == (11, 9, 7, 5, 3, 1), "lengths")


H = 24 * d - 35
G = 24 * d**2 - 259 * d + 315
I = 28 * d - 45
J = 256 * d**4 - 33792 * d**3 + 187360 * d**2 - 348768 * d + 218295
J_PLUS = (
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


W13 = (
    140928614400 * d**13
    - 4450810855424 * d**12
    + 2187898619166720 * d**11
    - 119479225810944000 * d**10
    + 2401299899413954560 * d**9
    - 26074123979998887936 * d**8
    + 177312698542783856640 * d**7
    - 809582018591761817600 * d**6
    + 2562896130982456504320 * d**5
    - 5668624022456998123776 * d**4
    + 8630773322305700966400 * d**3
    - 8649353474475612069600 * d**2
    + 5150611778276800836000 * d
    - 1384085705619177185625
)

W20 = (
    6734508720128000 * d**20
    - 9560308579093708800 * d**19
    - 104590212636279308288 * d**18
    + 77325440060392923463680 * d**17
    - 10990353099918947941089280 * d**16
    + 585398859813499740975267840 * d**15
    - 16234517220500771292093349888 * d**14
    + 280269110820744266239872860160 * d**13
    - 3319321835643659813075309035520 * d**12
    + 28555759263867535936704384860160 * d**11
    - 184855141179174118241589312618496 * d**10
    + 920280708402038646121095278100480 * d**9
    - 3567288992936033291861681133649920 * d**8
    + 10820035211743102835629862352322560 * d**7
    - 25630500379637406272477883155103744 * d**6
    + 46975014922828624287489671780597760 * d**5
    - 65372851083689767775383071292512000 * d**4
    + 66843975116236841444114398164960000 * d**3
    - 47400207255249718865601589483500000 * d**2
    + 20842506093742030140968421821925000 * d
    - 4284107145365556700983792046921875
)


D4 = (
    2**28
    * (2 * d - 19)
    * (2 * d - 17) ** 2
    * (2 * d - 15) ** 3
    * (2 * d - 13) ** 4
    * (2 * d - 11) ** 7
    * (2 * d - 9) ** 10
    * (2 * d - 7) ** 17
    * (2 * d - 5) ** 24
)
D5 = (
    2**69
    * (2 * d - 23)
    * (2 * d - 21) ** 2
    * (2 * d - 19) ** 3
    * (2 * d - 17) ** 4
    * (2 * d - 15) ** 7
    * (2 * d - 13) ** 10
    * (2 * d - 11) ** 17
    * (2 * d - 9) ** 24
    * (2 * d - 7) ** 41
    * (2 * d - 5) ** 58
)


def rf(expression: sp.Expr) -> object:
    return RATIONAL_FUNCTIONS.from_expr(expression.xreplace({d: sp.Symbol("d")}))


EXPECTED_W0 = rf(-(d - 1) ** 12 * H**4 * J**2 * W13 / D4)
EXPECTED_X0 = rf(-(d - 1) ** 29 * H**10 * J**4 * U**2 * W20 / D5)
require(W[0] == EXPECTED_W0, "fourth pivot")
require(X[0] == EXPECTED_X0, "fifth pivot")


def direct_moment_coefficient(degree_value: int, exponent: int, power: int) -> int:
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


# A second path expands the original moments before applying E.
for degree_value in range(12, 21):
    normalization = factorial(2 * degree_value - 4)
    direct_a = [
        Fraction(
            direct_moment_coefficient(
                degree_value, degree_value - 2, degree_value - 2 - index
            ),
            normalization,
        )
        for index in range(11)
    ]
    direct_b = [
        Fraction(
            direct_moment_coefficient(
                degree_value, degree_value - 1, degree_value - 1 - index
            ),
            normalization,
        )
        for index in range(11)
    ]
    direct_r = pseudo(direct_a, direct_b)
    direct_s = pseudo(direct_r, direct_a)
    direct_t = pseudo(direct_s, direct_r)
    direct_w = pseudo(direct_t, direct_s)
    direct_x = pseudo(direct_w, direct_t)
    require(
        (direct_w[0], direct_x[0])
        == (
            evaluate_rational_function(EXPECTED_W0, degree_value),
            evaluate_rational_function(EXPECTED_X0, degree_value),
        ),
        ("direct moment", degree_value),
    )


def polynomial_power_mod(
    base: sp.Poly, exponent: int, modulus_polynomial: sp.Poly
) -> sp.Poly:
    answer = sp.Poly(1, d, modulus=modulus_polynomial.get_modulus())
    factor = base.rem(modulus_polynomial)
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = (answer * factor).rem(modulus_polynomial)
        factor = (factor * factor).rem(modulus_polynomial)
        remaining //= 2
    return answer


def rabin_irreducible(polynomial: sp.Expr, prime: int) -> None:
    reduced = sp.Poly(polynomial, d, modulus=prime)
    n = sp.degree(polynomial, d)
    require(reduced.degree() == n, ("degree drop", n, prime))
    variable = sp.Poly(d, d, modulus=prime)
    require(
        (polynomial_power_mod(variable, prime**n, reduced) - variable)
        .rem(reduced)
        .is_zero,
        ("Frobenius closure", n, prime),
    )
    for divisor in sp.factorint(n):
        test = polynomial_power_mod(variable, prime ** (n // divisor), reduced)
        require(sp.gcd(reduced, test - variable).degree() == 0, (n, prime, divisor))


require(sp.Poly(W13, d).content() == 1, "W13 primitive")
require(sp.Poly(W20, d).content() == 1, "W20 primitive")
rabin_irreducible(W13, 47)
rabin_irreducible(W20, 29)
require(sp.Poly(W13, d, modulus=2) == sp.Poly(1, d, modulus=2), "W13 parity")
require(sp.Poly(W20, d, modulus=2) == sp.Poly(1, d, modulus=2), "W20 parity")


row_scalings = [1, 1]
for _ in range(5):
    row_scalings.append(2 * row_scalings[-1] + row_scalings[-2])
require(row_scalings == [1, 1, 3, 7, 17, 41, 99], row_scalings)
pell = [0, 1]
for _ in range(4):
    pell.append(2 * pell[-1] + pell[-2])
require(pell[1:] == [1, 2, 5, 12, 29], pell)


def modular_value(value: object, prime: int, degree_value: int) -> int:
    rational = value.as_expr().subs(sp.Symbol("d"), degree_value)
    numerator = int(rational.p) % prime
    denominator = int(rational.q) % prime
    require(denominator != 0, ("denominator wall", prime, degree_value))
    return numerator * pow(denominator, -1, prime) % prime


WALL_TABLE = {
    prime: tuple(
        modular_value(value, prime, prime + 2)
        for value in (T[0], W[0], W[1], X[0])
    )
    for prime in (41, 43)
}
require(WALL_TABLE[41] == (6, 15, 2, 0), WALL_TABLE[41])
require(WALL_TABLE[43] == (14, 0, 24, 23), WALL_TABLE[43])
require((24**2 * 14) % 43 == 23, "anchored-zero return")


# The wall transition is an identity of whole rows, not only constants.
# On f_0=0, E(f,g)=f_1*g_0*z^(-1)f.  The transported row then makes the
# following fraction-free row vanish identically.  A double anchored zero
# already collapses at the first wall step.
wall_clutch_controls = 0
for seed in range(1, 9):
    f = [0] + [seed * (index + 2) + index**2 + 1 for index in range(1, 8)]
    g = [2 * seed + 1] + [seed + 3 * index + index**2 for index in range(1, 8)]
    transported = pseudo(f, g)
    scalar = f[1] * g[0]
    require(
        transported
        == [scalar * f[index + 1] for index in range(len(transported))],
        ("whole-row wall clutch", seed),
    )
    require(
        all(coefficient == 0 for coefficient in pseudo(transported, f)),
        ("post-clutch terminal collapse", seed),
    )
    double_zero = [0, 0] + f[2:]
    require(
        all(coefficient == 0 for coefficient in pseudo(double_zero, g)),
        ("double-zero collapse", seed),
    )
    wall_clutch_controls += 3


# The two new walls are isolated from every THM-3217 numerator at s=2.
OLD_WALLS = (H, G, I, J, J_PLUS, K, U, V)
for prime, selected, other in ((41, W20, W13), (43, W13, W20)):
    require(int(selected.subs(d, 2)) % prime == 0, ("selected wall", prime))
    require(int(other.subs(d, 2)) % prime != 0, ("other new wall", prime))
    require(
        all(int(polynomial.subs(d, 2)) % prime for polynomial in OLD_WALLS),
        ("old wall collision", prime),
    )
    require(
        all((4 - odd) % prime for odd in range(5, 24, 2)),
        ("denominator collision", prime),
    )


source = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source)
require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert")
require(
    not any(
        isinstance(atom, sp.Float)
        for atom in sp.preorder_traversal(
            sp.Tuple(H, G, I, J, J_PLUS, K, U, V, W13, W20, D4, D5)
        )
    ),
    "float",
)
require(lf_sha256(DEPENDENCY) == EXPECTED_DEPENDENCY_SHA256, "dependency drift")


print("THM-3223 exact fourth/fifth resonant PRS wall tower")
print("dependency_pin=THM-3217")
print("input_top_jets=11 exact_rational_pivots=2 PASS")
print("direct_moment_controls=9 coordinate_checks=18 PASS")
print("rho4_factors=(d-1)^12*H^4*J^2*W13 denominator_2_power=28")
print("rho5_factors=(d-1)^29*H^10*J^4*U^2*W20 denominator_2_power=69")
print("primitive_degrees=13,20 irreducible_mod_primes=47,29 parity_mod2=1,1")
print("pell_d_minus_1_exponents=1,2,5,12,29")
print("normalization_scalings=1,1,3,7,17,41,99")
print("offset2_p41_table=(6,15,2,0) offset2_p43_table=(14,0,24,23)")
print("anchored_zero_return_p43=23")
print("wall_clutch_controls=24 whole_row=PASS next_row_zero=PASS")
print("fixed_offset_extension=finite_exception_set_through_row5")
print("scope=selected_pivots_not_full_row_or_arbitrary_depth")
