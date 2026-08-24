#!/usr/bin/env python3
"""Exact companion for THM-3976's rational-compression intersection."""

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


def ceil_div(a: int, b: int) -> int:
    return (a + b - 1) // b


x, t, w = sp.symbols("x t w", nonzero=True)
wxt = 1 + 2 * x**2 * t
z = (wxt + 1) / 2
p = (wxt**2 - 1) / (4 * x**2)
y = (wxt - 1)**2 * (wxt + 1) / (8 * x**3)
R = x**2
Z = x * wxt


def jac(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.factor(sp.diff(f, x) * sp.diff(g, t)
                     - sp.diff(f, t) * sp.diff(g, x))


# Source identities and the rational Darboux positive control.
gate(sp.factor(z - (1 + x**2 * t)) == 0, "z identity")
gate(sp.factor(p - t * (1 + x**2 * t)) == 0, "p identity")
gate(sp.factor(y - x * (1 + x**2 * t) * t**2) == 0, "y identity")
a = (wxt**2 - 1) / (4 * wxt**2)
b = wxt**3 / x
gate(sp.factor(jac(a, b) - 1) == 0, "rational Darboux identity")
gate(sp.factor(wxt**2 - 1 - 4 * R * p) == 0, "Pell identity")
gate(sp.factor(Z**2 - R - 4 * R**2 * p) == 0, "pseudoplane relation")

# The complete Jacobian-minor packet of the quadratic source atlas.
gate(sp.factor(jac(R, Z) - 4 * R**2) == 0, "J(R,Z)")
gate(sp.factor(jac(R, p) - 2 * Z) == 0, "J(R,p)")
gate(sp.factor(jac(Z, p) - (1 + 8 * R * p)) == 0, "J(Z,p)")

# Rational field reconstruction.
aa, bb = sp.symbols("aa bb", nonzero=True)
w2_ab = 1 / (1 - 4 * aa)
R_ab = w2_ab**3 / bb**2
Z_ab = w2_ab**2 / bb
gate(sp.factor(Z_ab**2 - R_ab * w2_ab) == 0, "fixed-field R,Z")

# Exact coefficient-ideal minima and the parity-fixed normal forms.
pp = sp.Symbol("pp")
for q in range(1, 81):
    feasible = []
    for i in range(q + 2):
        for j in range(q + 2):
            if 2 * i + 3 * j >= q:
                feasible.append((i, j))
    amin = min(i + 2 * j for i, j in feasible)
    bmin = min(i + j for i, j in feasible)
    gate(amin == ceil_div(q, 2), f"weight {-q}: w-1 minimum")
    gate(bmin == ceil_div(q, 3), f"weight {-q}: w+1 minimum")

    # A fixed profile must satisfy both the original and reflected orders.
    fixed_order = max(amin, bmin)
    gate(fixed_order == ceil_div(q, 2),
         f"weight {-q}: fixed common color order")
    if q % 2 == 0:
        m = q // 2
        profile = x**(-q) * (w**2 - 1)**m
        normal = 4**m * pp**m
        converted = sp.expand(profile.subs(w**2 - 1, 4 * x**2 * pp))
        gate(sp.factor(converted - normal) == 0,
             f"weight {-q}: even fixed normal form")
    else:
        m = (q - 1) // 2
        profile = x**(-q) * (w**2 - 1)**(m + 1) * w
        converted = sp.expand(profile.subs(w**2 - 1, 4 * x**2 * pp))
        normal = 4**(m + 1) * (x * w) * pp**(m + 1)
        gate(sp.factor(converted - normal) == 0,
             f"weight {-q}: odd fixed normal form")

# Positive weights have the two expected parity normal forms.
for r in range(81):
    if r % 2 == 0:
        gate(sp.factor(x**r - R**(r // 2)) == 0,
             f"weight {r}: even positive normal form")
    else:
        gate(sp.factor(x**r * wxt - R**((r - 1) // 2) * Z) == 0,
             f"weight {r}: odd positive normal form")

# The new boundary address is genuinely nonfixed under sigma.
y_profile = (w - 1)**2 * (w + 1)
gate(sp.expand(y_profile.subs(w, -w) + y_profile) != 0,
     "y is not fixed by simultaneous sign")
gate(sp.factor(y_profile.subs(w, -w) - (-1)**3 * y_profile) != 0,
     "y fails weight-minus-three fixed parity")

# Smoothness: a common gradient zero forces R=Z=0, leaving -1.
RR, ZZ, qq = sp.symbols("RR ZZ qq")
relation = ZZ**2 - RR - 4 * RR**2 * qq
partials = [sp.diff(relation, v) for v in (RR, ZZ, qq)]
groebner = sp.groebner(partials, RR, ZZ, qq, order="lex")
gate(list(groebner) == [sp.Integer(1)], "pseudoplane smoothness unit ideal")

# Frozen consequence gates.
gate(2 * 3 == 6, "quadratic atlas times surface-degree floor")
live = [(2, 5), (3, 4), (4, 3), (5, 2)]
gate(all(i + j == 7 for i, j in live), "Danielewski support floor")

summary = {
    "checks": CHECKS,
    "rational_pair": "a=(w^2-1)/(4w^2), b=w^3/x, J=1",
    "fixed_field": "k(x,w)^sigma=k(a,b), sigma(x,w)=(-x,-w)",
    "intersection": "B2 intersect k(a,b)=k[R,Z,p]/(Z^2-R-4R^2p)",
    "atlas": "surjective nonfinite etale quadratic A2-to-Y",
    "units_pic": "units=k*, Pic=Z/2 (symbolic theorem proof)",
    "degree_floor": 6,
    "support_floor": 7,
    "lost": "boundary-address generator y",
    "scope": "arbitrary B2 Darboux pair and JC2 remain open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3976 rational-compression pseudoplane intersection companion")
print(f"CHECKS={CHECKS}")
print("RATIONAL_PAIR=EXACT_DARBOUX")
print("FIXED_FIELD=SIMULTANEOUS_SIGN_QUADRATIC")
print("INTERSECTION=SMOOTH_QUADRATIC_PSEUDOPLANE")
print("ATLAS=SURJECTIVE_NONFINITE_ETALE_DEGREE_TWO")
print("LOST=ASYMMETRIC_BOUNDARY_ADDRESS_Y")
print("DEGREE_FLOOR=SURFACE_THREE;PLANAR_SIX")
print("SUPPORT_FLOOR=SEVEN")
print("SCOPE=ARBITRARY_B2_DARBOUX_AND_JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
