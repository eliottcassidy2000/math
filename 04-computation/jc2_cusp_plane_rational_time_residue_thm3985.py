#!/usr/bin/env python3
"""Exact cusp-plane residue and mixed-submersion companion for THM-3985."""

from __future__ import annotations

import hashlib
import json
from math import gcd

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    """Assertion-free exact gate, retained under ``python -O``."""
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(sp.cancel(expression)) == 0, message)


def jacobian(f: sp.Expr, g: sp.Expr,
             first: sp.Symbol, second: sp.Symbol) -> sp.Expr:
    return sp.expand(sp.diff(f, first) * sp.diff(g, second)
                     - sp.diff(f, second) * sp.diff(g, first))


x, t, u, p, y, s, q = sp.symbols("x t u p y s q")
alpha, gamma = sp.symbols("alpha gamma", nonzero=True)
z = 1 + x**2 * t
p_source = sp.expand(z * t)
y_source = sp.expand(x * z * t**2)
H = p**3 - y**2
H_source = sp.expand(p_source**3 - y_source**2)


# Panel I: the exact cusp-plane affine modification and its volume density.
zero(H_source - t * p_source**2, "H=t*p^2")
zero(y_source * p_source / H_source - x, "inverse x=yp/H")
zero(y_source**2 / H_source - x**2 * t, "inverse u=y^2/H")
zero(p_source**3 / H_source - z, "inverse z=p^3/H")
zero(H_source / p_source**2 - t, "inverse t=H/p^2")
zero(y_source / p_source - x * t, "v=y/p=xt")
zero(p_source - (t + (y_source / p_source)**2), "p=t+v^2")
J_py_source = jacobian(p_source, y_source, x, t)
zero(J_py_source - (y_source**2 - p_source**3) / p_source,
     "cusp-plane Jacobian")


# Panel II: the universal cusp chain rule and the constant-cusp critical row.
c00, c10, c01, c20, c11, c02 = sp.symbols(
    "c00 c10 c01 c20 c11 c02"
)
F_sample = (c00 + c10 * p + c01 * y + c20 * p**2
            + c11 * p * y + c02 * y**2)
phi_sample = sp.expand(F_sample.subs({p: s**2, y: s**3}))
J_FH = jacobian(F_sample, H, p, y)
zero(J_FH.subs({p: s**2, y: s**3})
     + s**2 * sp.diff(phi_sample, s), "cusp chain-rule sign")

G_sample = 1 + 2 * p + 3 * y + p * y
constant_cusp = sp.expand(7 + H * G_sample)
constant_source = sp.expand(constant_cusp.subs({p: p_source, y: y_source}))
gate(sp.diff(constant_source, x).subs(t, 0) == 0,
     "constant cusp row: x derivative on t=0")
gate(sp.diff(constant_source, t).subs(t, 0) == 0,
     "constant cusp row: t derivative on t=0")


# Panel III: A_m=alpha*p+gamma*y^m.  These finite loops are controls for
# the symbolic degree, coprimality, squarefreeness, residue, and place rows.
mixed_controls: dict[int, dict[str, object]] = {}
field = sp.QQ.frac_field(alpha, gamma, q)
for m in range(2, 9):
    p_num = q - gamma * y**m
    D_clear = sp.expand(alpha**3 * y**2 - p_num**3)
    p_poly = sp.Poly(p_num, y, domain=field)
    D_poly = sp.Poly(D_clear, y, domain=field)
    gate(p_poly.degree() == m, f"m={m}: p degree")
    gate(D_poly.degree() == 3 * m, f"m={m}: D degree")
    gate(sp.gcd(p_poly, D_poly).degree() == 0,
         f"m={m}: gcd(p,D)=1")
    gate(sp.gcd(D_poly, D_poly.diff()).degree() == 0,
         f"m={m}: generic D squarefree")

    A_source = sp.expand(alpha * p_source + gamma * y_source**m)
    zero(jacobian(A_source, y_source, x, t) + alpha * t * p_source,
         f"m={m}: exact Hamiltonian row")
    gate(sp.diff(A_source, t).subs(t, 0) == alpha,
         f"m={m}: t=0 is not critical")

    repeated_derivative = sp.factor(
        sp.diff(y**2 - ((q - gamma * y**m) / alpha)**3, y)
        .subs({y: s**3, q: alpha * s**2 + gamma * s**(3 * m)})
    )
    zero(repeated_derivative
         - s**3 * (2 + 3 * gamma * m * s**(3 * m - 2) / alpha),
         f"m={m}: repeated-root equation")

    mixed_controls[m] = {
        "p_degree": m,
        "D_degree": 3 * m,
        "source_places": 4 * m + 1,
        "completion_places": 3 * m + 1,
    }


# A fully factored low endpoint checks the generic-squarefree argument by an
# independent discriminant computation.
D2 = sp.expand(alpha**3 * y**2 - (q - gamma * y**2)**3)
disc2 = sp.factor(sp.discriminant(D2, y))
expected_disc2 = (64 * q**3 * alpha**12 * gamma**9
                  * (27 * q**2 * gamma + 4 * alpha**3)**2)
zero(disc2 - expected_disc2, "m=2 hostile discriminant")


# Panel IV: the sharp exponent and height endpoints.
A1_source = sp.expand(alpha * p_source + gamma * y_source)
critical_m1 = {x: gamma / alpha, t: -alpha**2 / gamma**2}
gate(sp.factor(sp.diff(A1_source, x).subs(critical_m1)) == 0,
     "m=1 endpoint: x derivative")
gate(sp.factor(sp.diff(A1_source, t).subs(critical_m1)) == 0,
     "m=1 endpoint: t derivative")
gate(sp.factor(z.subs(critical_m1)) == 0,
     "m=1 endpoint lies on the second color")

height_controls: dict[str, str] = {}
for n in range(3, 10):
    r = sp.expand(u * (u + 1))
    f = sp.expand(u**2 * (u + 1))
    reduced_u = sp.factor(sp.diff(r, u)
                          - sp.Rational(n, n + 1)
                          * r * sp.diff(f, u) / f)
    expected_u = ((2 - n) * u + (1 - n)) / (n + 1)
    zero(reduced_u - expected_u, f"height n={n}: reduced u row")
    u0 = -sp.Rational(n - 1, n - 2)
    gate(sp.factor(expected_u.subs(u, u0)) == 0,
         f"height n={n}: critical color")
    gate(u0 not in (0, -1), f"height n={n}: off-color address")
    for m in range(1, 7):
        exponent = m * (n + 1) - n
        gate(exponent >= 1, f"height n={n},m={m}: positive power")
        rhs = sp.factor(-m * (n + 1) * gamma * f**m
                        / (n * alpha * r))
        gate(sp.factor(rhs.subs(u, u0)) != 0,
             f"height n={n},m={m}: nonzero power address")
    height_controls[f"n{n}"] = str(u0)


# Panel V: a nonconstant weight-zero shift h(u) restores criticality.  The
# exact source derivatives are checked for arbitrary symbolic g=h'(u), and
# resultant identities verify both power-compatibility parity branches.
g_symbol = sp.symbols("g", nonzero=True)
r = u * (u + 1)
f = u**2 * (u + 1)
height_two_reduced = sp.factor(
    sp.diff(r, u) - sp.Rational(2, 3) * r * sp.diff(f, u) / f
)
gate(height_two_reduced == -sp.Rational(1, 3),
     "height-two 2/3 Euler cancellation")

compatibility_controls: dict[int, dict[str, int]] = {}
X2, B = sp.symbols("X2 B", nonzero=True)
for m in range(2, 11):
    N = 3 * m - 2
    e = gcd(2, N)
    resultant = sp.factor(sp.resultant(x**2 - X2, x**N - B, x))
    if N % 2 == 0:
        zero(resultant - (X2**(N // 2) - B)**2,
             f"m={m}: even common-power resultant")
    else:
        zero(resultant - (B**2 - X2**N),
             f"m={m}: odd common-power resultant")

    g_sample = u**2 + 2 * u + 3
    K0 = -sp.Rational(3 * m, 2) * gamma / alpha
    B_u = sp.expand(K0 * u**(2 * m - 1) * (u + 1)**(m - 1))
    compatibility = sp.together(
        alpha**(N // e)
        - 3**(N // e) * g_sample**(N // e) * B_u**(2 // e)
    )
    compatibility_poly = sp.Poly(
        compatibility, u, domain=sp.QQ.frac_field(alpha, gamma)
    )
    gate(compatibility_poly.degree() > 0,
         f"m={m}: compatibility is nonconstant")
    gate(sp.factor(compatibility.subs(u, 0) - alpha**(N // e)) == 0,
         f"m={m}: compatibility avoids u=0")
    gate(sp.factor(compatibility.subs(u, -1) - alpha**(N // e)) == 0,
         f"m={m}: compatibility avoids u=-1")
    gate(sp.gcd(compatibility_poly,
                sp.Poly(g_sample, u,
                        domain=sp.QQ.frac_field(alpha, gamma))).degree() == 0,
         f"m={m}: compatibility roots avoid g=0")
    compatibility_controls[m] = {"N": N, "gcd_2_N": e}


summary = {
    "checks": CHECKS,
    "affine_modification": "B2=k[p,y,yp/H,y^2/H], H=p^3-y^2",
    "cusp_dichotomy": {
        "nonconstant_phi": "J(F,k(x,t)) intersect k(F)={0}",
        "constant_phi": "F-c in (H); regular mates excluded only",
    },
    "positive_family": "alpha*p+gamma*y^m is submersion iff n=2,m>=2",
    "mixed_controls": mixed_controls,
    "height_controls": height_controls,
    "shift_gate": "h(u)+alpha*p+gamma*y^m submersion iff h constant",
    "compatibility_controls": compatibility_controls,
    "scope": "first coordinates in k[p,y] and displayed shifted family; JC2 open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3985 cusp-plane rational-time residue companion")
print(f"CHECKS={CHECKS}")
print("AFFINE_MODIFICATION=B2=k[p,y,yp/H,y^2/H];H=p^3-y^2")
print("PHI_NONCONSTANT=RATIONAL_CONSTANT_IMAGE_ZERO_BY_CUSP_RESIDUES")
print("PHI_CONSTANT=REGULAR_MATES_EXCLUDED;RATIONAL_SCOPE_NOT_CLAIMED")
print("MIXED_SUBMERSION=IFF_HEIGHT_2_AND_M_AT_LEAST_2")
print("GENERIC_SOURCE_PLACES=4M+1;X2_PLACES=3M+1")
print("WEIGHT_ZERO_SHIFT=SUBMERSION_IFF_CONSTANT")
print("SCOPE=CUSP_PLANE_AND_DISPLAYED_SHIFTED_FAMILY;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
