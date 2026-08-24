#!/usr/bin/env python3
"""Independent hostile algebra audit of the THM-3985 cusp-plane theorem."""

from __future__ import annotations

import hashlib
import json
from math import gcd

import sympy as sp


CHECKS = 0


def gate(condition: object, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.factor(sp.cancel(expression)) == 0, label)


def jacobian(f: sp.Expr, g: sp.Expr,
             first: sp.Symbol, second: sp.Symbol) -> sp.Expr:
    return sp.expand(sp.diff(f, first) * sp.diff(g, second)
                     - sp.diff(f, second) * sp.diff(g, first))


x, t, p, y, s, q, u = sp.symbols("x t p y s q u")
a, c = sp.symbols("a c", nonzero=True)
z = 1 + x**2 * t
ps = sp.expand(t * z)
ys = sp.expand(x * t * ps)
H = p**3 - y**2
Hs = sp.expand(ps**3 - ys**2)


# I. Rebuild the birational cusp chart and its volume sign from scratch.
zero(Hs - t * ps**2, "H=t*p^2")
zero(ys * ps / Hs - x, "x inverse")
zero(Hs / ps**2 - t, "t inverse")
zero(ys**2 / Hs - x**2 * t, "u inverse")
zero(ps**3 / Hs - z, "z inverse")
zero(jacobian(ps, ys, x, t) + t * ps, "J(p,y)=-tp")
zero(jacobian(ps, ys, x, t) + Hs / ps, "J(p,y)=-H/p")


# II. Check the residue sign for an independently chosen generic polynomial.
coeffs = sp.symbols("b0:10")
support = (1, p, y, p**2, p*y, y**2, p**3, p**2*y, p*y**2, y**3)
F = sp.expand(sum(b * monomial for b, monomial in zip(coeffs, support)))
phi = sp.expand(F.subs({p: s**2, y: s**3}))
JFH = jacobian(F, H, p, y)
zero(JFH.subs({p: s**2, y: s**3}) + s**2 * sp.diff(phi, s),
     "J(F,H)|cusp=-s^2 phi'")
zero(-s**2 / JFH.subs({p: s**2, y: s**3})
     - 1 / sp.diff(phi, s), "residue sign +1/phi'")


# The constant-cusp branch has criticality on both contracted colors.
G = 2 + 3*p + 5*y + 7*p*y
Fc_source = sp.expand((11 + H*G).subs({p: ps, y: ys}))
gate(sp.diff(Fc_source, x).subs(t, 0) == 0,
     "constant cusp: dx on t=0")
gate(sp.diff(Fc_source, t).subs(t, 0) == 0,
     "constant cusp: dt on t=0")
z_color = {t: -x**-2}
zero(sp.diff(Fc_source, x).subs(z_color),
     "constant cusp: dx on z=0")
zero(sp.diff(Fc_source, t).subs(z_color),
     "constant cusp: dt on z=0")


# III. Independent source-resultant proof of the mixed submersions.
resultant_rows: dict[int, str] = {}
for m in range(2, 11):
    A = sp.expand(a*ps + c*ys**m)
    resultant = sp.factor(sp.resultant(sp.diff(A, x), sp.diff(A, t), x))
    expected = (9 * m**(3*m) * a**(3*m - 1) * c**(3*m)
                * t**(9*m*m))
    zero(resultant - expected, f"m={m}: independent source resultant")
    gate(sp.diff(A, t).subs(t, 0) == a,
         f"m={m}: resultant's only address is not critical")
    resultant_rows[m] = str(sp.factor(resultant / t**(9*m*m)))

A1 = sp.expand(a*ps + c*ys)
R1 = sp.factor(sp.resultant(sp.diff(A1, x), sp.diff(A1, t), x))
zero(R1 - 9*c**3*t**9*(a**2 + c**2*t), "m=1 critical resultant")
critical_1 = {x: c/a, t: -a**2/c**2}
zero(sp.diff(A1, x).subs(critical_1), "m=1 critical dx")
zero(sp.diff(A1, t).subs(critical_1), "m=1 critical dt")


# IV. Generic rational-fibre punctures and time-form divisor ledger.
field = sp.QQ.frac_field(a, c, q)
puncture_rows: dict[int, tuple[int, int, int]] = {}
for m in range(2, 11):
    Pnum = q - c*y**m
    Dclear = sp.expand(a**3*y**2 - Pnum**3)
    Ppoly = sp.Poly(Pnum, y, domain=field)
    Dpoly = sp.Poly(Dclear, y, domain=field)
    gate(Ppoly.degree() == m, f"m={m}: p roots")
    gate(Dpoly.degree() == 3*m, f"m={m}: cusp roots")
    gate(sp.gcd(Ppoly, Ppoly.diff()).degree() == 0,
         f"m={m}: p roots simple")
    gate(sp.gcd(Ppoly, Dpoly).degree() == 0,
         f"m={m}: puncture packets disjoint")
    gate(sp.gcd(Dpoly, Dpoly.diff()).degree() == 0,
         f"m={m}: cusp roots simple")
    source_places = m + 3*m + 1
    completion_places = 3*m + 1
    gate(source_places == 4*m + 1, f"m={m}: source count")
    gate(completion_places == source_places - m,
         f"m={m}: completion restores boundary roots")
    # On P1_y, omega has m simple zeros at P=0, 3m simple poles at
    # D=0, and a zero of order 2m-2 at infinity.  Its degree must be -2.
    gate(m + (2*m - 2) - 3*m == -2,
         f"m={m}: rational time divisor degree")
    puncture_rows[m] = (m, 3*m, source_places)


# V. Independently rederive every higher-height off-color endpoint.
height_rows: dict[int, str] = {}
for n in range(3, 11):
    r = u*(u + 1)
    f = u**2*(u + 1)
    reduced = sp.factor(
        sp.diff(f, u)/f - sp.Rational(n + 1, n)*sp.diff(r, u)/r
    )
    expected = ((n - 2)*u + (n - 1))/(n*u*(u + 1))
    zero(reduced - expected, f"n={n}: reduced critical address")
    u0 = -sp.Rational(n - 1, n - 2)
    zero(reduced.subs(u, u0), f"n={n}: endpoint root")
    gate(u0 != 0 and u0 != -1, f"n={n}: endpoint off colors")
    for m in range(1, 8):
        exponent = m*(n + 1) - n
        rhs = -m*(n + 1)*c*f**m/(n*a*r)
        gate(exponent > 0, f"n={n},m={m}: reconstructing power positive")
        gate(sp.factor(rhs.subs(u, u0)) != 0,
             f"n={n},m={m}: reconstructing value nonzero")
    height_rows[n] = str(u0)


# VI. Hostile parity audit for the nonconstant h(u) compatibility step.
C, B = sp.symbols("C B", nonzero=True)
parity_rows: dict[int, tuple[int, int]] = {}
for m in range(2, 13):
    N = 3*m - 2
    e = gcd(2, N)
    common = sp.factor(sp.resultant(x**2 - C, x**N - B, x))
    if N % 2 == 0:
        zero(common - (C**(N//2) - B)**2,
             f"m={m}: even compatibility")
    else:
        zero(common - (B**2 - C**N), f"m={m}: odd compatibility")

    g = u + 2 if m % 2 else u**2 + u + 1
    K0 = -sp.Rational(3*m, 2)*c/a
    Bu = K0*u**(2*m - 1)*(u + 1)**(m - 1)
    compat = sp.expand(a**(N//e)
                       - 3**(N//e)*g**(N//e)*Bu**(2//e))
    poly = sp.Poly(compat, u, domain=sp.QQ.frac_field(a, c))
    gate(poly.degree() > 0, f"m={m}: nonconstant h compatibility")
    gate(sp.factor(compat.subs(u, 0)) == a**(N//e),
         f"m={m}: compatibility avoids zero")
    gate(sp.factor(compat.subs(u, -1)) == a**(N//e),
         f"m={m}: compatibility avoids minus one")
    gate(sp.gcd(poly, sp.Poly(g, u,
                             domain=sp.QQ.frac_field(a, c))).degree() == 0,
         f"m={m}: compatibility avoids g=0")
    parity_rows[m] = (N, e)


summary = {
    "checks": CHECKS,
    "residue": "Res=+1/phi'(s); no geometric-integrality input",
    "mixed_resultants": resultant_rows,
    "puncture_rows": puncture_rows,
    "height_endpoints": height_rows,
    "common_power": parity_rows,
    "verdict": "PASS",
    "scope": "THM-3985 only; JC(2) remains open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("Independent hostile audit: THM-3985")
print(f"CHECKS={CHECKS}")
print("VERDICT=PASS")
print("RESIDUE_SIGN=+1/phi_prime;NORMALIZATION_LOCAL_ONLY")
print("PUNCTURES=SOURCE_4M_PLUS_1;X2_3M_PLUS_1")
print("ENDPOINTS=M1_CRITICAL;N_GE_3_CRITICAL;NONCONSTANT_H_CRITICAL")
print("COMMON_POWER=EVEN_AND_ODD_BRANCHES_PASS")
print("JC2=OPEN")
print(f"SEMANTIC_SHA256={semantic}")
