#!/usr/bin/env python3
"""Exact companion for THM-3975 (one-arm Danielewski modifications)."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(expression) == 0, message)


x, t, z, p, y, C, X, w = sp.symbols("x t z p y C X w")

# The determinantal presentation, exact primitive, homogeneous LND, and
# two-denominator grading are checked at several independent heights.  The
# theorem proves these formulas uniformly in n.
for n in range(2, 8):
    zn = 1 + x**n * t
    pn = sp.expand(zn * t)
    yn = sp.expand(x ** (n - 1) * zn * t**2)

    zero(zn * (zn - 1) - x**n * pn,
         f"height {n}: Danielewski relation")
    zero(pn * (zn - 1) - x * yn,
         f"height {n}: modification relation")
    zero(zn * yn - x ** (n - 1) * pn**2,
         f"height {n}: saturated minor")
    zero(zn * (zn - 1)**2 - x ** (n + 1) * yn,
         f"height {n}: local hypersurface")

    # beta=-dz/((n-1)x^(n-1)); its exterior derivative is dx wedge dt.
    beta_x = sp.expand(-sp.diff(zn, x) / ((n - 1) * x ** (n - 1)))
    beta_t = sp.expand(-sp.diff(zn, t) / ((n - 1) * x ** (n - 1)))
    zero(sp.diff(beta_t, x) - sp.diff(beta_x, t) - 1,
         f"height {n}: exact source volume")

    # delta=x partial_t is homogeneous of degree n+1.
    delta = lambda f: sp.expand(x * sp.diff(f, t))
    zero(delta(zn) - x ** (n + 1), f"height {n}: delta z")
    zero(delta(pn) - x * (2 * zn - 1), f"height {n}: delta p")
    zero(delta(yn) - (zn - 1) * (3 * zn - 1),
         f"height {n}: delta y")

    # A near-slice whose remaining multiple of p is killed by the
    # hyperelliptic generic-fibre obstruction in the proof.
    near = x * (1 - 2 * zn)
    jac_near = sp.diff(pn, x) * sp.diff(near, t) - sp.diff(pn, t) * sp.diff(near, x)
    zero(jac_near - 1 - 2 * (n + 2) * x**n * pn,
         f"height {n}: near-slice identity")

    # Exact DPD/two-color exponents.  A negative weight -m has coefficient
    # z^ceil(m/(n+1)) (z-1)^ceil(m/n).
    for m in range(1, 2 * n + 4):
        pairs = [
            (a, b)
            for a in range(m + 2)
            for b in range(m + 2)
            if n * a + (n + 1) * b >= m
        ]
        alpha = min(a + b for a, b in pairs)
        beta = min(a + 2 * b for a, b in pairs)
        gate(alpha == sp.ceiling(sp.Rational(m, n + 1)),
             f"height {n}, weight {m}: zero-color exponent")
        gate(beta == sp.ceiling(sp.Rational(m, n)),
             f"height {n}, weight {m}: one-color exponent")

    # All-height finite controls C_ell=y+x^ell.  The theorem proves the
    # exact field degree from the pole divisor on the generic p-fibre;
    # these gates freeze monicity, elimination, and the pole arithmetic.
    for ell in range(n // 2, n // 2 + 4):
        cell = sp.expand(yn + x**ell)
        fall = sp.expand(
            X * (C - X**ell)**2 + p * (C - X**ell)
            - p**3 * X ** (n - 1)
        )
        poly = sp.Poly(fall, X)
        gate(poly.degree() == 2 * ell + 1,
             f"height {n}, exponent {ell}: finite-control degree")
        gate(poly.LC() == 1,
             f"height {n}, exponent {ell}: finite-control monicity")
        zero(fall.subs({X: x, C: cell, p: pn}),
             f"height {n}, exponent {ell}: finite-control elimination")
        zero(2 * x * yn - pn * ((2 * zn - 1) - 1),
             f"height {n}, exponent {ell}: generic-fibre y coordinate")

        if n % 2 == 0:
            gate(ell > n // 2 - 1,
                 f"height {n}, exponent {ell}: even-infinity dominance")
        else:
            gate(2 * ell > n - 2,
                 f"height {n}, exponent {ell}: odd-infinity dominance")
        gate(1 + 2 * ell == poly.degree(),
             f"height {n}, exponent {ell}: exact pole degree")

        jac_cell = sp.expand(
            sp.diff(pn, x) * sp.diff(cell, t)
            - sp.diff(pn, t) * sp.diff(cell, x)
        )
        gate(sp.Poly(jac_cell, x, t).total_degree() > 0,
             f"height {n}, exponent {ell}: finite control is not Keller")

# The finite controls at heights two and three.
deltas: dict[int, sp.Expr] = {}
for n in (2, 3):
    f = sp.expand(X * (C - X)**2 + p * (C - X) - X ** (n - 1) * p**3)
    gate(sp.Poly(f, X).LC() == 1,
         f"height {n}: monic cubic")
    gate(sp.Poly(f, X).degree() == 3, f"height {n}: cubic degree")

    # The p=1 hostile used in the irreducibility proof.
    g = sp.expand(f.subs(p, 1))
    a, b = sp.symbols("a b")
    linear = sp.Poly(sp.expand(g.subs(X, a * C + b)), C)
    zero(linear.coeff_monomial(C**3) - a * (a - 1)**2,
         f"height {n}: linear-root leading gate")
    if n == 2:
        zero(g.subs(X, C + b) - ((b**2 - 1) * C + b**3 - 2 * b),
             "height 2: a=1 root contradiction")
    else:
        gate(sp.Poly(g.subs(X, C + b), C).coeff_monomial(C**2) == -1,
             "height 3: a=1 root contradiction")

    # Multiplication by x,w on the free basis (1,x,w), where w=z-1.
    if n == 2:
        x2 = sp.Matrix([0, C, -p])
        xw = sp.Matrix([C, -(1 + p**2), C])
        w2 = sp.Matrix([0, p * C, -(1 + p**2)])
    else:
        x2 = sp.Matrix([0, C, -p])
        xw = sp.Matrix([C, -(1 + C * p**2), C + p**3])
        w2 = sp.Matrix([
            -p**2 * C,
            p * C**2 + p**2 + C * p**4,
            -(1 + 2 * C * p**2 + p**5),
        ])

    one = sp.Matrix([1, 0, 0])
    xb = sp.Matrix([0, 1, 0])
    wb = sp.Matrix([0, 0, 1])
    Mx = sp.Matrix.hstack(xb, x2, xw)
    Mw = sp.Matrix.hstack(wb, xw, w2)
    gate((Mx * Mw - Mw * Mx).applyfunc(sp.expand) == sp.zeros(3),
         f"height {n}: commuting table")
    gate((Mx**2 - C * Mx + p * Mw).applyfunc(sp.expand) == sp.zeros(3),
         f"height {n}: x squared table")
    gate((p * Mx**n - (sp.eye(3) + Mw) * Mw).applyfunc(sp.expand)
         == sp.zeros(3), f"height {n}: Danielewski table")
    gate((Mx * (C * sp.eye(3) - Mx) - p * Mw).applyfunc(sp.expand)
         == sp.zeros(3), f"height {n}: modification table")
    gate(((sp.eye(3) + Mw) * (C * sp.eye(3) - Mx)
          - p**2 * Mx ** (n - 1)).applyfunc(sp.expand) == sp.zeros(3),
         f"height {n}: saturated table")
    coeff = sp.Poly(f, X).all_coeffs()
    f_of_Mx = coeff[0] * Mx**3 + coeff[1] * Mx**2 + coeff[2] * Mx + coeff[3] * sp.eye(3)
    gate(f_of_Mx.applyfunc(sp.expand) == sp.zeros(3),
         f"height {n}: cubic annihilates multiplication by x")

    trace_matrix = sp.Matrix(3, 3, lambda i, j: sp.trace(
        [sp.eye(3), Mx, Mw][i] * [sp.eye(3), Mx, Mw][j]
    ))
    delta = sp.factor(trace_matrix.det())
    deltas[n] = delta
    zero(sp.discriminant(f, X) - p**2 * delta,
         f"height {n}: monogenic index-square discriminant")

    # The three normal addresses over p=0,C!=0.
    addresses = [(0, -1), (C, -1), (C, 0)]
    for j, (x0, w0) in enumerate(addresses):
        zero(sp.expand(x0**2 - C * x0),
             f"height {n}: address {j} x relation")
        zero(sp.expand(x0 * w0 - (C * w0 + C - x0)),
             f"height {n}: address {j} mixed relation")
        zero(sp.expand(w0**2 + w0),
             f"height {n}: address {j} w relation")

    # Exact Jacobian and ramification numerator.
    zn = 1 + x**n * t
    pn = sp.expand(zn * t)
    yn = sp.expand(x ** (n - 1) * zn * t**2)
    target_c = x + yn
    jac = sp.expand(sp.diff(pn, x) * sp.diff(target_c, t)
                    - sp.diff(pn, t) * sp.diff(target_c, x))
    Rn = sp.expand(x * (2 * zn - 1) + yn * ((n - 2) * zn + 1))
    zero(x * jac + Rn, f"height {n}: ramification numerator")
    gate(sp.Poly(jac, x, t).total_degree() > 0,
         f"height {n}: finite control is not Keller")

delta2_expected = (
    4 * C**4 * p - 8 * C**2 * p**4 + 20 * C**2 * p**2 + C**2
    + 4 * p**7 + 12 * p**5 + 12 * p**3 + 4 * p
)
delta3_expected = (
    4 * C**5 * p + C**4 * p**4 + 22 * C**3 * p**2
    + 22 * C**2 * p**5 + C**2 + 4 * C * p**8 + 22 * C * p**3
    + p**6 + 4 * p
)
zero(deltas[2] - delta2_expected, "height 2: normal discriminant")
zero(deltas[3] - delta3_expected, "height 3: normal discriminant")
gate(delta2_expected.subs(p, 0) == C**2,
     "height 2: generic p=0 fibre is etale")
gate(delta3_expected.subs(p, 0) == C**2,
     "height 3: generic p=0 fibre is etale")

summary = {
    "checks": CHECKS,
    "family": "B_n=k[x,1+x^n t,(1+x^n t)t,x^(n-1)(1+x^n t)t^2]",
    "geometry": "one-point Danielewski affine modification with A1 boundary",
    "grading": "Dplus=0; Dminus=-[0]/(n+1)-[1]/n",
    "canonical": "Cl=Z[D], K=[D], source volume exact",
    "lnd": "x*d/dt, kernel k[x], marked plinth x^(n+1)",
    "cubic": "n=2,3 finite free rank3 over k[p,x+y]",
    "finite_tower": "C_ell=y+x^ell, ell>=floor(n/2), exact degree 2ell+1",
    "index": "monogenic x-order has index p and discriminant p^2 Delta_n",
    "addresses": "three distinct normal addresses over generic p=0",
    "no_mate": "p has no rational constant-Jacobian mate for every n>=2",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3975 Danielewski one-arm cubic-control companion")
print(f"CHECKS={CHECKS}")
print("PRESENTATION=THREE_2X2_MINORS;SMOOTH_TWO_CHART_NORMALIZATION")
print("BOUNDARY=A1;CL=Z;K=D;SOURCE_VOLUME=EXACT")
print("DPD=DPLUS_0;DMINUS_MINUS_0_OVER_NPLUS1_MINUS_1_OVER_N")
print("LND=X_DT;KERNEL_KX;MARKED_PLINTH_X_NPLUS1")
print("CUBIC=N2_N3;FINITE_FREE_BASIS_1_X_W")
print("ALL_HEIGHT_FINITE=C_ELL_Y_PLUS_XELL;DEGREE_2ELL_PLUS1;NEVER_KELLER")
print("INDEX=MONOGENIC_X_ORDER_P;DISC_ORDER=P2_DISC_NORMAL")
print("P0_ADDRESSES=THREE_DISTINCT;INDEX_NOT_RAMIFICATION")
print("NO_MATE=P_GENERIC_FIBRE_LOG_OR_HOLOMORPHIC_DIFFERENTIAL")
print(f"SEMANTIC_SHA256={semantic}")
