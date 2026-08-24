#!/usr/bin/env python3
"""Independent exact hostile audit for THM-3989.

This file deliberately imports no canonical companion.  It attacks the six
load-bearing steps with different finite models: the birational chart, the
full Laurent coefficient cone, a bounded intersection equality (rather than
only residue visibility), the complete bounded syzygy kernel, convolution
signs and scalar moment, the UFD multiplicity law, and shear descent.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from math import gcd

import sympy as sp


GATES = 0


def gate(condition: object, label: str) -> None:
    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)


def vanishes(expression: sp.Expr, label: str) -> None:
    gate(sp.factor(sp.cancel(expression)) == 0, label)


def laurent_terms(expression: sp.Expr, variables: tuple[sp.Symbol, ...]) -> dict[tuple[int, ...], sp.Expr]:
    """Return an exact sparse exponent dictionary after full expansion."""
    result: dict[tuple[int, ...], sp.Expr] = {}
    for term in sp.Add.make_args(sp.expand(expression)):
        powers = term.as_powers_dict()
        exponents = tuple(int(powers.get(variable, 0)) for variable in variables)
        coefficient = term
        for variable, exponent in zip(variables, exponents):
            coefficient /= variable**exponent
        coefficient = sp.cancel(coefficient)
        result[exponents] = sp.expand(result.get(exponents, 0) + coefficient)
    return {key: value for key, value in result.items() if value != 0}


def bracket(first: sp.Expr, second: sp.Expr, q: sp.Symbol, r: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(first, q) * sp.diff(second, r)
        - sp.diff(first, r) * sp.diff(second, q)
    )


s, tau, p, y = sp.symbols("s tau p y")
H = p**3 - y**2
P = s**2 + tau
Y = s * P
X = s / tau
U = s**2 / tau


# I. Audit both directions of the logarithmic chart and its orientation.
vanishes(H.subs({p: P, y: Y}) - tau * P**2, "chart: H=tau*p^2")
vanishes(Y / P - s, "chart: s=y/p")
vanishes((P**3 - Y**2) / P**2 - tau, "chart: tau=H/p^2")
vanishes(Y * P / (P**3 - Y**2) - X, "chart: x=yp/H")
vanishes(Y**2 / (P**3 - Y**2) - U, "chart: u=y^2/H")
vanishes(P.subs({s: y / p, tau: H / p**2}) - p, "inverse chart: p")
vanishes(Y.subs({s: y / p, tau: H / p**2}) - y, "inverse chart: y")
vanishes(
    bracket(X, tau, s, tau) - 1 / tau,
    "orientation: dx wedge dtau=ds wedge dtau/tau",
)


# II. Exhaust a substantially larger monomial cone than the canonical grid.
# The test checks every coefficient, not only the negative ones.
monomial_rows = 0
negative_rows = 0
for total in range(0, 12):
    for a, b, c, e in itertools.product(range(total + 1), repeat=4):
        if a + b + c + e != total:
            continue
        actual = sp.expand(X**a * U**b * P**c * Y**e)
        closed = sp.expand(
            s ** (a + 2 * b + e)
            * tau ** (-a - b)
            * (s**2 + tau) ** (c + e)
        )
        vanishes(actual - closed, f"cone closed form {a,b,c,e}")
        terms = laurent_terms(actual, (s, tau))
        for (s_exp, tau_exp), coefficient in terms.items():
            gate(coefficient != 0, f"cone nonzero coefficient {a,b,c,e,s_exp,tau_exp}")
            if tau_exp < 0:
                j = -tau_exp
                gate(s_exp >= j, f"cone conductor {a,b,c,e,s_exp,tau_exp}")
                negative_rows += 1
        monomial_rows += 1

# Independently span every depth symbol s^(d+n), including the exceptional
# semigroup hole n=1, and verify that no deeper term occurs.
symbol_rows = 0
for d in range(1, 13):
    for n in range(0, 31):
        if n == 0:
            witness = X**d
        elif n == 1:
            witness = X ** (d - 1) * U
        elif n % 2 == 0:
            witness = X**d * P ** (n // 2)
        else:
            witness = X**d * P ** ((n - 3) // 2) * Y
        terms = laurent_terms(witness, (s, tau))
        lowest_tau = min(q for _, q in terms)
        gate(lowest_tau == -d, f"symbol depth d={d},n={n}")
        gate(terms.get((d + n, -d), 0) == 1, f"symbol value d={d},n={n}")
        symbol_rows += 1


# III. Bounded full intersection equality.  For generator degree <=5, form
# the exact Laurent vector space, impose all negative-coefficient equations,
# and compare the resulting *whole polynomials* with the p,y-monomial span.
def vector_matrix(expressions: list[sp.Expr], basis: list[tuple[int, int]]) -> sp.Matrix:
    sparse = [laurent_terms(expression, (s, tau)) for expression in expressions]
    return sp.Matrix([
        [entry.get(exponent, 0) for entry in sparse]
        for exponent in basis
    ])


intersection_rows: dict[int, dict[str, int]] = {}
for cutoff in range(1, 6):
    b2_monomials: list[sp.Expr] = []
    for exponents in itertools.product(range(cutoff + 1), repeat=4):
        if sum(exponents) <= cutoff:
            a, b, c, e = exponents
            b2_monomials.append(sp.expand(X**a * U**b * P**c * Y**e))
    cusp_monomials = [
        sp.expand(P**c * Y**e)
        for c in range(cutoff + 1)
        for e in range(cutoff + 1 - c)
    ]
    all_expressions = b2_monomials + cusp_monomials
    all_support = sorted({
        exponent
        for expression in all_expressions
        for exponent in laurent_terms(expression, (s, tau))
    })
    negative_support = [exponent for exponent in all_support if exponent[1] < 0]
    nonnegative_support = [exponent for exponent in all_support if exponent[1] >= 0]
    negative_matrix = vector_matrix(b2_monomials, negative_support)
    kernel = negative_matrix.nullspace()
    full_b2_matrix = vector_matrix(b2_monomials, nonnegative_support)
    intersection_matrix = (
        full_b2_matrix * sp.Matrix.hstack(*kernel)
        if kernel else sp.zeros(len(nonnegative_support), 0)
    )
    cusp_matrix = vector_matrix(cusp_monomials, nonnegative_support)
    intersection_rank = intersection_matrix.rank()
    cusp_rank = cusp_matrix.rank()
    joined_rank = intersection_matrix.row_join(cusp_matrix).rank()
    gate(intersection_rank == cusp_rank == joined_rank,
         f"full bounded intersection cutoff={cutoff}")
    # The normalization parameter s must remain absent in every polynomial
    # intersection vector; this is a particularly cheap hostile direction.
    hostile = sp.zeros(len(nonnegative_support), 1)
    if (1, 0) in nonnegative_support:
        hostile[nonnegative_support.index((1, 0)), 0] = 1
        gate(intersection_matrix.row_join(hostile).rank() > intersection_rank,
             f"normalization parameter absent cutoff={cutoff}")
    intersection_rows[cutoff] = {
        "b2_columns": len(b2_monomials),
        "kernel": len(kernel),
        "intersection_rank": intersection_rank,
    }


# IV. Solve the whole bounded pG+yK in (H) kernel and reconstruct L,M for
# every basis vector.  This tests necessity, not just parametrized positives.
syzygy_rows: dict[int, int] = {}
for cutoff in range(0, 6):
    py_basis = [(a, b) for a in range(cutoff + 1)
                for b in range(cutoff + 1 - a)]
    g_symbols = sp.symbols(f"g0:{len(py_basis)}")
    k_symbols = sp.symbols(f"k0:{len(py_basis)}")
    G = sp.Add(*(coefficient * p**a * y**b
                 for coefficient, (a, b) in zip(g_symbols, py_basis)))
    K = sp.Add(*(coefficient * p**a * y**b
                 for coefficient, (a, b) in zip(k_symbols, py_basis)))
    quotient, remainder = sp.div(sp.expand(p * G + y * K), H, y)
    remainder_poly = sp.Poly(sp.expand(remainder), p, y)
    equations = [coefficient for _, coefficient in remainder_poly.terms()]
    matrix, _ = sp.linear_eq_to_matrix(equations, g_symbols + k_symbols)
    nullspace = matrix.nullspace()
    for index, vector in enumerate(nullspace):
        substitution = dict(zip(g_symbols + k_symbols, vector))
        G0 = sp.expand(G.subs(substitution))
        K0 = sp.expand(K.subs(substitution))
        L0 = sp.cancel((p * G0 + y * K0) / H)
        gate(L0.is_polynomial(p, y), f"syzygy L cutoff={cutoff},row={index}")
        M0 = sp.cancel((G0 - p**2 * L0) / y)
        gate(M0.is_polynomial(p, y), f"syzygy M cutoff={cutoff},row={index}")
        vanishes(G0 - (p**2 * L0 + y * M0),
                 f"syzygy G reconstruction cutoff={cutoff},row={index}")
        vanishes(K0 - (-y * L0 - p * M0),
                 f"syzygy K reconstruction cutoff={cutoff},row={index}")
        vanishes(
            (Y * P / (P**3 - Y**2)) * G0.subs({p: P, y: Y})
            + (Y**2 / (P**3 - Y**2)) * K0.subs({p: P, y: Y})
            - Y * L0.subs({p: P, y: Y}),
            f"syzygy Laurent descent cutoff={cutoff},row={index}",
        )
    syzygy_rows[cutoff] = len(nullspace)


# V. Recompute the Laurent signs for unrelated coefficient polynomials over
# a wider support.  Then test moment and shallow jet on exact Keller hostiles.
a_coeff = {
    i: (i * i + 7) + (i + 4) * s + (2 * i * i + 3) * s**2
       + (i**3 + 11) * s**4
    for i in range(-6, 6)
}
c_coeff = {
    j: (j * j + 13) + (3 * j + 8) * s + (j**2 - j + 5) * s**3
       + (j**4 + 2) * s**5
    for j in range(-5, 7)
}
A = sp.Add(*(coefficient * tau**index for index, coefficient in a_coeff.items()))
C = sp.Add(*(coefficient * tau**index for index, coefficient in c_coeff.items()))
direct = sp.expand(tau * bracket(A, C, s, tau))
for weight in range(-11, 12):
    expected = sp.expand(sum(
        j * sp.diff(a_coeff[i], s) * c_coeff[j]
        - i * a_coeff[i] * sp.diff(c_coeff[j], s)
        for i in a_coeff for j in c_coeff if i + j == weight
    ))
    vanishes(direct.coeff(tau, weight) - expected,
             f"convolution sign weight={weight}")
moment = sp.expand(sum(i * a_coeff[i] * c_coeff.get(-i, 0) for i in a_coeff))
vanishes(direct.coeff(tau, 0) + sp.diff(moment, s), "moment derivative sign")

for q in range(1, 10):
    A0 = X
    C0 = sp.expand(tau + X**q)
    vanishes(tau * bracket(A0, C0, s, tau) - 1, f"Keller hostile q={q}")
    coefficients_a = {index: sp.expand(A0).coeff(tau, index)
                      for index in range(-q - 1, q + 2)}
    coefficients_c = {index: sp.expand(C0).coeff(tau, index)
                      for index in range(-q - 1, q + 2)}
    moment0 = sp.expand(sum(
        index * coefficients_a[index] * coefficients_c.get(-index, 0)
        for index in coefficients_a
    ))
    vanishes(moment0 + s, f"scalar moment q={q}")
    shallow = (
        coefficients_a.get(1, 0).subs(s, 0)
        * sp.diff(coefficients_c.get(-1, 0), s).subs(s, 0)
        - sp.diff(coefficients_a.get(-1, 0), s).subs(s, 0)
        * coefficients_c.get(1, 0).subs(s, 0)
    )
    vanishes(shallow + 1, f"shallow jet q={q}")


# VI. Independent UFD multiplicity audit and target-shear descent.
ufd_rows = 0
for d, e in ((1, 1), (1, 4), (2, 3), (2, 5), (3, 4), (4, 6), (6, 10), (6, 15), (8, 12)):
    g = gcd(d, e)
    for m, r in ((1, 1), (2, -2), (3, 3)):
        h = s**g * (s - r)**m * (s**2 + s + 1)
        aa = sp.Rational(2, 3) * h ** (d // g)
        cc = sp.Rational(-5, 7) * h ** (e // g)
        vanishes(d * aa * sp.diff(cc, s) - e * sp.diff(aa, s) * cc,
                 f"UFD positive d={d},e={e},m={m}")
        gate(sp.rem(aa, s**d, s) == 0, f"UFD a conductor d={d},e={e}")
        gate(sp.rem(cc, s**e, s) == 0, f"UFD c conductor d={d},e={e}")
        # A one-root perturbation is the minimal hostile to common powers.
        hostile = sp.expand(cc + s**e * (s - r + 1))
        gate(sp.factor(d * aa * sp.diff(hostile, s)
                       - e * sp.diff(aa, s) * hostile) != 0,
             f"UFD hostile d={d},e={e},m={m}")
        ufd_rows += 1


def depth(expression: sp.Expr) -> int:
    terms = laurent_terms(expression, (s, tau))
    return max(0, -min(exponent[1] for exponent in terms))


shear_rows = 0
for d, q in ((1, 2), (1, 5), (2, 2), (2, 3), (3, 4), (5, 3)):
    e = q * d
    h = s**d * (s + 2)
    alpha = sp.Rational(2, 3)
    beta = sp.Rational(-7, 5)
    A1 = sp.expand(alpha * h * tau**(-d) + (s + 1) * tau**(-d + 1) + tau**2)
    C1 = sp.expand(beta * h**q * tau**(-e) + (s**2 + 3) * tau**(-e + 1) + tau)
    before = depth(C1)
    C2 = sp.expand(C1 - beta / alpha**q * A1**q)
    after = depth(C2)
    gate(before == e and after < e, f"shear lowers d={d},e={e}")
    vanishes(
        tau * bracket(A1, C2, s, tau) - tau * bracket(A1, C1, s, tau),
        f"shear preserves bracket d={d},e={e}",
    )
    gate(depth(A1) + after < d + e, f"shear sum descent d={d},e={e}")
    shear_rows += 1


summary = {
    "gates": GATES,
    "monomial_rows": monomial_rows,
    "negative_rows": negative_rows,
    "symbol_rows": symbol_rows,
    "intersection_rows": intersection_rows,
    "syzygy_rows": syzygy_rows,
    "convolution_weights": 23,
    "moment_hostiles": 9,
    "ufd_rows": ufd_rows,
    "shear_rows": shear_rows,
    "verdict": "PASS; theorem quantifiers and signs survive independent hostile audit",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3989 independent every-line hostile audit")
print(f"GATES={GATES}")
print(f"MONOMIAL_ROWS={monomial_rows};NEGATIVE_ROWS={negative_rows};SYMBOL_ROWS={symbol_rows}")
print(f"INTERSECTION={intersection_rows}")
print(f"SYZYGY_KERNELS={syzygy_rows}")
print("CONVOLUTION_SIGNS=PASS;SCALAR_MOMENT=-S;SHALLOW_JET=-1")
print(f"UFD_ROWS={ufd_rows};SHEAR_ROWS={shear_rows}")
print("VERDICT=PASS;NON-DIVIDING_CELLS_AND_JC2_REMAIN_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
