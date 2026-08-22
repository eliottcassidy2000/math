#!/usr/bin/env python3
"""Exact controls for THM-3551's one-invariant planar mate no-go.

The all-degree statements are proved symbolically in the theorem.  This
companion checks the sector identities, rational cages, gradient hostiles,
bounded response systems, and the proposed high-genus ray factorization.
Every gate uses explicit exceptions so normal and optimized runs agree.
"""

from __future__ import annotations

import hashlib
from math import gcd

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def jac(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(f, x) * sp.diff(g, y) - sp.diff(f, y) * sp.diff(g, x))


def bounded_response(P: sp.Expr, cap: int) -> tuple[int, int, int, int]:
    mons = [x**i * y**j for i in range(cap + 1) for j in range(cap + 1 - i)]
    images = [sp.Poly(jac(P, mon), x, y) for mon in mons]
    rhs = sp.Poly(1, x, y)
    output = sorted(set().union(*(poly.monoms() for poly in images), rhs.monoms()))
    matrix = sp.Matrix([[poly.coeff_monomial(mon) for poly in images] for mon in output])
    column = sp.Matrix([rhs.coeff_monomial(mon) for mon in output])
    return cap, len(mons), matrix.rank(), matrix.row_join(column).rank()


x, y, t, u, w, z = sp.symbols("x y t u w z")
records: list[object] = []


# Multiplicative primitive ray P=x*phi(x^p y^q).
sector_rows = []
for p in range(1, 6):
    for q in range(1, 6):
        if gcd(p, q) != 1:
            continue
        phi = 1 + 2 * t + 3 * t**2
        for degree in range(7):
            coeffs = sp.symbols(f"a_{p}_{q}_{degree}_0:{degree + 1}")
            f = sum(c * t**i for i, c in enumerate(coeffs))
            response = sp.expand(
                (phi + p * t * sp.diff(phi, t)) * f
                + q * t * phi * sp.diff(f, t)
            )
            poly = sp.Poly(response - 1, t)
            solution = sp.linsolve(poly.all_coeffs(), coeffs)
            require(solution == sp.EmptySet, "multiplicative-ray response unexpectedly solved")
            top = sp.Poly(response, t).coeff_monomial(t ** (degree + 2))
            expected = 3 * (1 + 2 * p + q * degree) * coeffs[-1]
            require(sp.expand(top - expected) == 0, "multiplicative top multiplier")
        sector_rows.append((p, q, 1 + 2 * p))

T = x * y
phi_xy = 1 + T + T**2
P_cage = sp.expand(x * phi_xy)
Q_cage = y / phi_xy
require(sp.cancel(jac(P_cage, Q_cage)) == 1, "rational cage Jacobian")
phi_uw = 1 + u * w + (u * w) ** 2
inverse = (u / phi_uw, w * phi_uw)
require(
    sp.cancel(P_cage.subs({x: inverse[0], y: inverse[1]}, simultaneous=True)) == u
    and sp.cancel(Q_cage.subs({x: inverse[0], y: inverse[1]}, simultaneous=True)) == w,
    "rational cage inverse",
)
require(list(sp.groebner([sp.diff(P_cage, x), sp.diff(P_cage, y)], x, y)) == [1],
        "multiplicative hostile is not a submersion")
records.append(("multiplicative", tuple(sector_rows), sp.factor(P_cage)))


# Additive diagonal P=x+h(xy): the generic rational response has simple residues.
h = t**2 + t**4
P_diag = sp.expand(x + h.subs(t, x * y))
resultant = sp.factor(sp.resultant(z - h, sp.diff(h, t), t))
require(resultant != 0, "diagonal generic denominator is inseparable")
require(list(sp.groebner([sp.diff(P_diag, x), sp.diff(P_diag, y)], x, y)) == [1],
        "diagonal hostile is not a submersion")
diag_ranks = tuple(bounded_response(P_diag, cap) for cap in (1, 2, 3, 4, 6, 8, 10))
require(all(rank + 1 == augmented for _, _, rank, augmented in diag_ranks),
        "bounded diagonal response unexpectedly solved")
require(jac(P_diag, x * y) == x, "diagonal birational-coordinate Jacobian")
records.append(("diagonal", resultant, diag_ranks))


# Thickened primitive ray P=x+h(x^p y^q), q>=2: residue-one top coefficient.
residue_rows = []
for p, q in ((1, 2), (2, 3), (3, 4), (3, 5)):
    require(gcd(p, q) == 1, "nonprimitive residue control")
    for N in range(1, 4 * q + 2, q):
        coeffs = sp.symbols(f"b_{p}_{q}_{N}_0:{p * N + q + 2}")
        A = sum(c * x**i for i, c in enumerate(coeffs))
        obstruction = sp.Poly(p * N * A - q * x * sp.diff(A, x), x)
        solution = sp.linsolve(obstruction.all_coeffs(), coeffs)
        require(solution == sp.FiniteSet(tuple(0 for _ in coeffs)),
                "residue-one top kernel unexpectedly nonzero")
        residue_rows.append((p, q, N))

V = x * y**2
P_residue = sp.expand(x + V**2 + V**3)
require(list(sp.groebner([sp.diff(P_residue, x), sp.diff(P_residue, y)], x, y)) == [1],
        "residue hostile is not a submersion")
records.append(("residue", tuple(residue_rows), sp.factor(P_residue)))


# The apparent hyperelliptic escape is one multiplicative ray after refactoring.
a, lam = 4, sp.Integer(1)
T_hyp = x ** (a - 1) * y**a
P_hyp = sp.expand(x + (x * y) ** a + lam * x ** (2 * a - 1) * y ** (2 * a))
require(sp.expand(P_hyp - x * (1 + T_hyp + lam * T_hyp**2)) == 0,
        "hyperelliptic scaffold did not collapse to one ray")
require(list(sp.groebner([sp.diff(P_hyp, x), sp.diff(P_hyp, y)], x, y)) == [1],
        "squarefree hyperelliptic hostile is not a submersion")
records.append(("hyperelliptic-ray", a, sp.factor(P_hyp)))


# Positive constant-family controls.
P_positive = 3 * x + 7
Q_positive = sp.Rational(1, 3) * y + P_positive**2
require(jac(P_positive, Q_positive) == 1, "positive affine/shear control")
records.append(("positive", P_positive, Q_positive))


semantic = hashlib.sha256(repr(tuple(records)).encode("utf-8")).hexdigest()
print("THM3551 ONE-INVARIANT MATE NO-GO EXACT CONTROLS")
print(f"MULTIPLICATIVE primitive_pairs={len(sector_rows)} degree_caps=0..6 rational_cage=PASS")
print(f"DIAGONAL resultant={resultant} bounded={diag_ranks} rational_residue=PASS")
print(f"RESIDUE top_cells={len(residue_rows)} gradient_unit=PASS")
print("HYPERELLIPTIC proposed_escape=ONE_RAY lambda=1 a=4 NO_MATE")
print("POSITIVE constant_family=PASS")
print(f"SEMANTIC_SHA256={semantic}")
print("PASS")
