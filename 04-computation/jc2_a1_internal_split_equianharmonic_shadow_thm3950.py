#!/usr/bin/env python3
"""Exact companion for THM-3950's A1 split debt and j=0 shadow."""

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


d = sp.symbols("delta")


def dreduce(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.together(sp.cancel(expression)).as_numer_denom()
    num = sp.rem(sp.Poly(sp.expand(numerator), d), sp.Poly(d**2 + 3, d))
    den = sp.rem(sp.Poly(sp.expand(denominator), d), sp.Poly(d**2 + 3, d))
    return sp.cancel(num.as_expr() / den.as_expr())


def dzero(expression: sp.Expr, message: str) -> None:
    gate(dreduce(expression) == 0, message)


omega = (-1 + d) / 2
omega2 = (-1 - d) / 2
dzero(omega**2 + omega + 1, "omega relation")
dzero(omega - omega2 - d, "delta relation")
dzero(d**2 + 3, "delta square")


# ---------------------------------------------------------------------------
# Universal normalized-A1 pullback.  The abstract proof first extracts
# r,s from p_i=r_i^2,q_i=2r_i^3 and then uses gcd(r,s)=g.  Here R,S are
# the coprime cores and c is the remaining common multiplier.
# ---------------------------------------------------------------------------

R, S, c, h = sp.symbols("R S c h")
U = S + omega2 * R
V = S - omega * R
g = dreduce(c * U * V)
r = dreduce(g * R)
s = dreduce(g * S)
a = dreduce(c * (S - R) * V**2)
b = dreduce(c * (S + R) * U**2)

dzero(U * V - (S**2 - d * R * S - R**2),
      "two denominator factors multiply correctly")
dzero(a * b - (s**2 - r**2), "assigned factors multiply to p1-p0")
dzero(a * (s**2 - omega * r**2) - (s**3 - r**3),
      "first internally split row")
dzero(b * (s**2 - omega2 * r**2) - (s**3 + r**3),
      "second internally split row")

# Before inserting g=cUV, cancellation gives aU=g(S-R)V and
# bV=g(S+R)U.  These freeze the divisibility step used in the proof.
dzero(a * U - g * (S - R) * V, "first denominator-debt equation")
dzero(b * V - g * (S + R) * U, "second denominator-debt equation")


# ---------------------------------------------------------------------------
# Companion cubic and its irreducible quadratic shadow.
# ---------------------------------------------------------------------------

D = dreduce((1 - omega2) * b - (1 - omega) * a)
E = dreduce((b - a) * a * b)
C = dreduce(E + D * h**2 - 2 * h**3)
B0 = dreduce(D - 2 * r)
Qh = dreduce(-2 * h**2 + B0 * h + r * B0)

dzero(E + D * r**2 - 2 * r**3, "r is a companion-cubic root")
dzero(C - (h - r) * Qh, "companion cubic factors linear times quadratic")
dzero(D - (3 * r + d * c * S**3), "compressed D identity")
dzero(D - 2 * r + c * (R - S) * (R + S) * (R + d * S),
      "first quadratic-discriminant factor")
dzero(D + 6 * r + c * (3 * R + d * S) ** 3 / 3,
      "second quadratic-discriminant factor")

disc_h = dreduce(sp.discriminant(Qh, h))
F4 = (R - S) * (R + S) * (R + d * S) * (3 * R + d * S)
dzero(disc_h - (D - 2 * r) * (D + 6 * r),
      "quadratic shadow discriminant")
dzero(
    disc_h
    - c**2 * (R - S) * (R + S) * (R + d * S)
    * (3 * R + d * S) ** 3 / 3,
    "exact four-color discriminant invoice",
)
dzero(disc_h / (c**2 * (3 * R + d * S) ** 2 / 3) - F4,
      "quadratic shadow square class")


# ---------------------------------------------------------------------------
# The fixed four-color cover is nonsingular and equianharmonic.
# ---------------------------------------------------------------------------

x, y, z = sp.symbols("x y z")
f4 = dreduce((x - 1) * (x + 1) * (x + d) * (3 * x + d))
coeff = sp.Poly(f4, x).all_coeffs()
qa, qb, qc, qd, qe = coeff
I4 = dreduce(12 * qa * qe - 3 * qb * qd + qc**2)
J4 = dreduce(
    72 * qa * qc * qe + 9 * qb * qc * qd
    - 27 * qa * qd**2 - 27 * qb**2 * qe - 2 * qc**3
)

gate(sp.degree(f4, x) == 4, "fixed cover is quartic")
gate(dreduce(sp.discriminant(f4, x)) == -110592,
     "fixed quartic is squarefree")
gate(I4 == 0, "binary quartic invariant I vanishes")
gate(J4 == 1728, "binary quartic invariant J is nonzero")
gate(dreduce(4 * I4**3 - J4**2) == -2985984,
     "binary quartic discriminant invariant is nonzero")
dzero(F4 - S**4 * f4.subs(x, R / S),
      "four-color packet is pullback of fixed quartic")


# ---------------------------------------------------------------------------
# Degree-three ratio map and its S3 Galois-closure curve.
# ---------------------------------------------------------------------------

Nphi = sp.expand((1 - x) * (1 - omega * x) ** 2)
Dphi = sp.expand((1 + x) * (1 + omega2 * x) ** 2)
phi = sp.cancel(Nphi / Dphi)
Fz = dreduce(Nphi - z * Dphi)
disc_Fz = dreduce(sp.discriminant(Fz, x))

gate(max(sp.degree(Nphi, x), sp.degree(Dphi, x)) == 3,
     "universal factor-ratio map has degree three")
dzero(a / b - phi.subs(x, R / S),
      "assigned-factor ratio is the degree-three map")
dzero(
    disc_Fz + 6 * z * (z - 1) * ((d + 3) * z + d - 3),
    "degree-three map branch-value discriminant",
)
dzero(((3 - d) / (3 + d)) + omega,
      "fourth finite expression is the branch value -omega")

offdiag_num = sp.expand(Nphi.subs(x, y) * Dphi - Nphi * Dphi.subs(x, y))
offdiag = dreduce(sp.cancel(offdiag_num / (y - x)))
offdiag_expected = (
    x**2 + y**2 + 4 * x * y
    - d * (x**2 * y + x * y**2) + d * (x + y)
)
dzero(offdiag - offdiag_expected, "off-diagonal fiber quadratic")
dzero(sp.discriminant(offdiag, y) + f4,
      "off-diagonal discriminant is the j=0 quartic")

# Exact degree-three isogeny sidecar between the two j=0 models.
isogeny_multiplier = x * (x - omega2) / ((x + 1) ** 2 * (x + omega) ** 3)
dzero(
    phi * (phi - 1) * (phi + omega)
    - (3 + d) * isogeny_multiplier**2 * f4 / 6,
    "explicit j=0 three-isogeny square identity",
)


# ---------------------------------------------------------------------------
# Smallest positive survivor: one reduced A1 graph, irreducible elliptic
# residual, normal quadratic surface, and two local A2 character witnesses.
# ---------------------------------------------------------------------------

P, Y, W = sp.symbols("P Y W")
gY = dreduce((Y + omega2) * (Y - omega))
AY = dreduce((Y - 1) * (Y - omega) ** 2)
BY = dreduce((Y + 1) * (Y + omega2) ** 2)
p0 = P
p1 = sp.expand(P + AY * BY)
L1 = p1 - omega * P
L2 = p1 - omega2 * P
q0 = dreduce(BY * L2 - AY * L1)
q1 = dreduce(BY * L2 + AY * L1)
H = dreduce(q0**2 - 4 * P**3)
DY = dreduce((1 - omega2) * BY - (1 - omega) * AY)
NY = dreduce(
    4 * P**2 - (DY**2 - 4 * gY**2) * P
    + gY**2 * (DY - 2 * gY) ** 2
)

dzero(gY - (Y**2 - d * Y - 1), "explicit common-root debt")
dzero(p1.subs(P, gY**2) - (Y * gY) ** 2,
      "second cusp square on the graph")
dzero(q0.subs(P, gY**2) - 2 * gY**3,
      "first cusp cube on the graph")
dzero(q1.subs(P, gY**2) - 2 * (Y * gY) ** 3,
      "second cusp cube on the graph")
dzero(q1**2 - 4 * p1**3 - H, "explicit common discriminant")
dzero(H - (gY**2 - P) * NY, "A1 graph plus residual factorization")
gate(sp.Poly(H, P, Y).total_degree() == 14,
     "explicit full discriminant has degree fourteen")
gate(sp.Poly(NY, P, Y).total_degree() == 10,
     "explicit residual component has degree ten")
gate(sp.degree(NY, P) == 2, "residual is quadratic over the Y-line")
dzero(
    sp.diff(H, P).subs(P, gY**2) - 4 * d * gY**3 * Y**3,
    "A1 graph occurs with multiplicity one",
)

disc_NY = dreduce(sp.discriminant(NY, P))
dzero(disc_NY - (DY - 2 * gY) ** 3 * (DY + 6 * gY),
      "residual P-quadratic discriminant")
explicit_f4 = (1 - Y) * (1 + Y) * (1 + d * Y) * (3 + d * Y)
dzero(
    disc_NY / ((DY - 2 * gY) ** 2 * (3 + d * Y) ** 2 / 3)
    - explicit_f4,
    "explicit residual square class",
)
gate(dreduce(sp.discriminant(explicit_f4, Y)) != 0,
     "explicit residual square class has four simple roots")

# The two exclusive local A2 witnesses used for Cardano independence.
dzero(BY - AY + d * Y + 1, "first exclusive witness factor")
Jlin = (1 + d) * Y + 3 - d
dzero(AY + omega * BY - Y**2 * Jlin / 2,
      "second exclusive witness factor")
gate(dreduce(sp.resultant(BY - AY, AY, Y)) == 2,
     "first witness avoids A=0")
gate(dreduce(sp.resultant(BY - AY, BY, Y)) == 2,
     "first witness avoids B=0")
gate(dreduce(sp.resultant(Jlin, AY, Y)) == -8 * (d + 1),
     "second simple witness avoids A=0")
gate(dreduce(sp.resultant(Jlin, BY, Y)) == -8 * (d - 1),
     "second simple witness avoids B=0")
dzero(q0.subs(P, 0) - (BY - AY) * AY * BY,
      "first A2 local cusp coordinate")
dzero(
    q1.subs(P, -AY * BY) - omega * AY * BY * (AY + omega * BY),
    "second A2 local cusp coordinate",
)


summary = {
    "checks": CHECKS,
    "pullback": "g=c(S+omega2R)(S-omegaR); exact A,B debt formulas",
    "shadow": "quadratic residual pulls back fixed nonsingular j=0 quartic",
    "ratio_map": "degree3 phi; branch values 0,1,-omega,infinity",
    "positive_survivor": "reduced A1 graph plus irreducible genus1 residual",
    "quadratic_surface": "normal with two independent local A2 Cardano classes",
    "jc2": "open; no affine-plane cubic atlas or Keller pair",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3950 A1 internal-split equianharmonic-shadow companion")
print(f"CHECKS={CHECKS}")
print("DEBT=g=c(S+omega^2R)(S-omegaR)")
print("A=c(S-R)(S-omegaR)^2;B=c(S+R)(S+omega^2R)^2")
print("SHADOW_DISC_SQUARECLASS=(R-S)(R+S)(R+delta*S)(3R+delta*S)")
print("FIXED_COVER=GENUS1;j=0;I=0;J=1728")
print("RATIO_MAP=DEGREE3;BRANCH_VALUES=0,1,-omega,infinity")
print("EXPLICIT_H=REDUCED_A1_GRAPH_PLUS_IRREDUCIBLE_GENUS1_RESIDUAL")
print("EXPLICIT_SURFACE=NORMAL;CARDANO_CHARACTER_RANK_AT_LEAST_2")
print("SCOPE=NO_AFFINE_PLANE_CUBIC_ATLAS;NO_KELLER_MAP;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
