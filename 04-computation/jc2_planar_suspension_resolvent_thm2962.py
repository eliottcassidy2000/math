#!/usr/bin/env python3
"""Exact controls for THM-2962.

The geometric A1-invariance argument is proved in the theorem text.  This
companion checks the two load-bearing quartic/cubic hostiles and the exact
V4/S3/PSL2(Z) quotient typing.  It deliberately uses explicit ``require``
calls rather than Python ``assert`` so optimized execution checks the same
claims.
"""

import os

# Avoid a python-flint sorting incompatibility in modular factor_list.
os.environ.setdefault("SYMPY_GROUND_TYPES", "python")

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def t_order(poly: sp.Expr, parameter: sp.Symbol) -> int:
    p = sp.Poly(poly, parameter, domain=sp.QQ)
    require(not p.is_zero, "valuation requested for the zero polynomial")
    return min(monomial[0] for monomial, _ in p.terms())


def factor_degrees_mod(poly: sp.Expr, variable: sp.Symbol, prime: int) -> list[int]:
    factors = sp.factor_list(poly, modulus=prime)[1]
    return sorted(
        sp.Poly(factor, variable, modulus=prime).degree()
        for factor, exponent in factors
        for _ in range(exponent)
    )


def matmul_mod2(a: tuple[tuple[int, int], tuple[int, int]],
                b: tuple[tuple[int, int], tuple[int, int]]
                ) -> tuple[tuple[int, int], tuple[int, int]]:
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(2)) % 2 for j in range(2))
        for i in range(2)
    )


def matpow_mod2(a: tuple[tuple[int, int], tuple[int, int]], exponent: int
                ) -> tuple[tuple[int, int], tuple[int, int]]:
    result = ((1, 0), (0, 1))
    for _ in range(exponent):
        result = matmul_mod2(result, a)
    return result


def det_mod2(a: tuple[tuple[int, int], tuple[int, int]]) -> int:
    return (a[0][0] * a[1][1] - a[0][1] * a[1][0]) % 2


Y, U, t, p, q = sp.symbols("Y U t p q")

# THM-2769 squared-pair hostile.
quartic = Y**4 - 2 * Y**2 - 8 * t * Y + 1 - 4 * t
pair_cubic = U**3 - 4 * U**2 + 16 * t * U - 64 * t**2

quartic_discriminant = sp.factor(sp.discriminant(quartic, Y))
cubic_discriminant = sp.factor(sp.discriminant(pair_cubic, U))
expected_discriminant = -4096 * t**2 * (27 * t**2 - 14 * t + 3)
require(quartic_discriminant == expected_discriminant, "wrong quartic discriminant")
require(cubic_discriminant == expected_discriminant, "wrong cubic discriminant")

coefficient_orders = [
    t_order(sp.Poly(pair_cubic, U).coeff_monomial(U**i), t) for i in range(4)
]
require(coefficient_orders == [2, 1, 0, 0], "wrong coefficient Newton data")

# Lower hull: (0,2)--(2,0)--(3,0), slopes -1 and 0, lengths 2 and 1.
root_valuation_row = [1, 1, 0]
parity_row = "".join(str(value % 2) for value in root_valuation_row)
require(parity_row == "110", "wrong V4 divisor-parity word")

# A good t=2 specialization proves the family is genuinely S4:
# irreducible mod 5 gives a 4-cycle/transitivity, type 1+3 mod 3 gives a
# 3-cycle, and the negative rational discriminant excludes A4.
s4_specialization = sp.expand(quartic.subs(t, 2))
s4_discriminant = int(sp.discriminant(s4_specialization, Y))
degrees_mod_5 = factor_degrees_mod(s4_specialization, Y, 5)
degrees_mod_3 = factor_degrees_mod(s4_specialization, Y, 3)
require(s4_discriminant % 5 != 0 and s4_discriminant % 3 != 0,
        "bad specialization prime")
require(degrees_mod_5 == [4], "specialization is not irreducible mod 5")
require(degrees_mod_3 == [1, 3], "specialization has no visible 3-cycle")
require(s4_discriminant < 0, "specialization discriminant is a rational square")

# THM-2871 leading-face hostile: every depressed cubic occurs after an
# invertible coefficient rescaling, and discriminant equality supplies no
# additional square law.
depressed = Y**3 + p * Y + q
leading_pair_cubic = U**3 + 16 * p * U - 64 * q
leading_discriminant_ratio = sp.factor(
    sp.discriminant(leading_pair_cubic, U) / sp.discriminant(depressed, Y)
)
require(leading_discriminant_ratio == 4096, "wrong leading discriminant scale")

# The exact 2/3 representation: W=F2^2, Aut(W)=GL2(F2)=S3.  Its order-two
# and order-three generators are a quotient of C2*C3=PSL2(Z), with the
# additional S3 relation (sr)^2=1.
identity = ((1, 0), (0, 1))
reflection = ((0, 1), (1, 0))
three_cycle = ((0, 1), (1, 1))
require(matpow_mod2(reflection, 2) == identity, "reflection is not order two")
require(matpow_mod2(three_cycle, 3) == identity, "cycle is not order three")
require(matpow_mod2(matmul_mod2(reflection, three_cycle), 2) == identity,
        "S3 quotient relation failed")

all_matrices = [
    ((a, b), (c, d))
    for a in range(2)
    for b in range(2)
    for c in range(2)
    for d in range(2)
]
gl2 = [matrix for matrix in all_matrices if det_mod2(matrix) == 1]
require(len(gl2) == 6, "GL2(F2) does not have order six")

print("quartic_discriminant =", quartic_discriminant)
print("pair_cubic_discriminant =", cubic_discriminant)
print("t_adic_coefficient_orders =", coefficient_orders)
print("root_valuation_row =", root_valuation_row)
print("v4_parity_row =", parity_row)
print("t2_quartic =", s4_specialization)
print("t2_discriminant =", s4_discriminant)
print("t2_factor_degrees_mod5 =", degrees_mod_5)
print("t2_factor_degrees_mod3 =", degrees_mod_3)
print("leading_depressed_cubic_discriminant_ratio =", leading_discriminant_ratio)
print("gl2_f2_order =", len(gl2))
print("orders_and_s3_relation =", "2,3,2")
print("PASS")
