#!/usr/bin/env python3
"""Exact algebraic companion for THM-3951's boundary-incidence obstruction."""

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
# Universal minimal-denominator-debt packet on an A1 normalization.
# R,S stand for coprime polynomial values at the generic parameter.
# ---------------------------------------------------------------------------

R, S, c, h, P, T = sp.symbols("R S c h P T")
U = S + omega2 * R
V = S - omega * R
g = dreduce(c * U * V)
r = dreduce(g * R)
s1 = dreduce(g * S)
A = dreduce(c * (S - R) * V**2)
B = dreduce(c * (S + R) * U**2)
E = dreduce((B - A) * A * B)
C = dreduce((1 - omega2) * B - (1 - omega) * A)
q0 = dreduce(E + C * P)
F = dreduce(T**3 - 3 * P * T - q0)

dzero(A * B - (s1**2 - r**2), "universal square-difference split")
dzero(A * (s1**2 - omega * r**2) - (s1**3 - r**3),
      "universal first cube row")
dzero(B * (s1**2 - omega2 * r**2) - (s1**3 + r**3),
      "universal second cube row")
dzero(C - (3 * r + d * c * S**3), "compressed coefficient C")
dzero(E + r**2 * (C - 2 * r), "graph-root constant identity")
dzero(q0.subs(P, r**2) - 2 * r**3, "graph lies on discriminant")
dzero(F.subs({P: r**2, T: -r}), "graph lies on cubic ramification")

# A polynomial root of the monic depressed cubic would force
# C^3+27E=0.  Dividing by B^3 gives a genuinely nonzero cubic equation
# for the nonconstant ratio A/B.
z = sp.symbols("z")
root_relation = dreduce(
    ((1 - omega2) - (1 - omega) * z) ** 3 + 27 * z * (1 - z)
)
gate(sp.degree(root_relation, z) == 3,
     "cubic root obstruction is a nonzero algebraic equation")
gate(dreduce(sp.LC(sp.Poly(root_relation, z), z)) != 0,
     "cubic root obstruction has nonzero leading coefficient")
dzero((C**3 + 27 * E) / B**3 - root_relation.subs(z, A / B),
      "universal cubic-domain root criterion")


# ---------------------------------------------------------------------------
# Exact source-ramification factorization and the three clean color fibres.
# ---------------------------------------------------------------------------

Q = dreduce(-2 * h**2 + (C - 2 * r) * h + r * (C - 2 * r))
ramification_hidden = dreduce(E + C * h**2 - 2 * h**3)
dzero(ramification_hidden - (h - r) * Q,
      "ramification factors as graph times residual")
dzero(sp.discriminant(Q, h) - (C - 2 * r) * (C + 6 * r),
      "residual discriminant")
dzero(Q.subs(h, r) - 2 * d * c * r * S**3,
      "universal graph-residual incidence identity")
dzero(sp.diff(F, T) - 3 * (T**2 - P), "relative derivative")
dzero(sp.diff(F, P).subs({P: r**2, T: -r}) + d * c * S**3,
      "generic graph prime is smooth and ramified")

# At each of the three colors R=0, U=0, V=0 one has r=0 and the same
# nonzero smoothness coefficient delta*c*S^3.  The substitutions for U=0
# and V=0 retain R as a generic nonzero color-fibre value.
color_substitutions = {
    "R=0": {R: 0},
    "U=0": {S: -omega2 * R},
    "V=0": {S: omega * R},
}
for label, substitution in color_substitutions.items():
    r_color = dreduce(r.subs(substitution))
    C_color = dreduce(C.subs(substitution))
    S_color = dreduce(S.subs(substitution))
    dzero(r_color, f"{label}: graph root vanishes")
    dzero(C_color - d * c * S_color**3,
          f"{label}: common smoothness coefficient")
    dzero(Q.subs(substitution).subs(h, 0),
          f"{label}: residual meets graph")
    dzero(sp.diff(Q, h).subs(substitution).subs(h, 0)
          - d * c * S_color**3,
          f"{label}: residual germ is smooth")


# ---------------------------------------------------------------------------
# Explicit R=1,S=Y,c=1 survivor.  This freezes irreducibility and two
# transverse affine intersections without relying on a genericity argument.
# ---------------------------------------------------------------------------

Y = sp.symbols("Y")
s = dreduce((Y + omega2) * (Y - omega))
AY = dreduce((Y - 1) * (Y - omega)**2)
BY = dreduce((Y + 1) * (Y + omega2)**2)
EY = dreduce((BY - AY) * AY * BY)
CY = dreduce((1 - omega2) * BY - (1 - omega) * AY)
q0Y = dreduce(EY + CY * P)
FY = dreduce(T**3 - 3 * P * T - q0Y)
QY = dreduce(-2 * h**2 + (s + d * Y**3) * h + s * (s + d * Y**3))
HY = dreduce(q0Y**2 - 4 * P**3)
NY = dreduce(
    4 * P**2 - (CY**2 - 4 * s**2) * P + s**2 * (CY - 2 * s)**2
)

dzero(s - (Y**2 - d * Y - 1), "explicit common debt")
dzero(CY - (3 * s + d * Y**3), "explicit compressed C")
dzero(EY + s**2 * (s + d * Y**3), "explicit compressed E")
dzero(q0Y.subs(P, s**2) - 2 * s**3, "explicit graph branch")
dzero(FY.subs({P: s**2, T: -s}), "explicit graph ramification")
dzero(HY - (s**2 - P) * NY, "explicit graph plus residual branch")
dzero(
    EY + CY * h**2 - 2 * h**3 - (h - s) * QY,
    "explicit ramification graph plus residual",
)
dzero(QY.subs(h, s) - 2 * d * s * Y**3,
      "explicit incidence polynomial")

disc_QY = dreduce(sp.discriminant(QY, h))
disc_factor = d * (Y - 1) * (Y + 1) * (d * Y + 1) * (Y - d)**3
dzero(disc_QY - (s + d * Y**3) * (9 * s + d * Y**3),
      "explicit residual discriminant product")
dzero(disc_QY - disc_factor, "explicit residual discriminant factorization")
square_class = sp.expand((Y - 1) * (Y + 1) * (d * Y + 1) * (Y - d))
gate(dreduce(sp.resultant(square_class, sp.diff(square_class, Y), Y)) != 0,
     "residual square class has four distinct roots")
gate(sp.degree(QY, h) == 2, "residual remains quadratic")

# The domain test is independent of the factorization of the discriminant.
gate(dreduce(EY.subs(Y, 1)) == 0, "domain hostile point has E=0")
gate(dreduce(CY.subs(Y, 1) + 2 * d) == 0,
     "domain hostile point has C=-2delta")
gate(dreduce((CY**3 + 27 * EY).subs(Y, 1)) != 0,
     "explicit cubic has no polynomial root")

y_points = (omega, -omega2)
gate(dreduce(y_points[0] - y_points[1]) != 0,
     "the two clean affine colors are distinct")
incidence_restriction = dreduce(QY.subs(h, s))
for index, y0 in enumerate(y_points, start=1):
    dzero(s.subs(Y, y0), f"intersection {index}: s vanishes")
    dzero(FY.subs({P: 0, T: 0, Y: y0}),
          f"intersection {index}: point lies on cubic")
    dzero(sp.diff(FY, T).subs({P: 0, T: 0, Y: y0}),
          f"intersection {index}: point is ramified")
    gate(dreduce(sp.diff(FY, P).subs({P: 0, T: 0, Y: y0})) != 0,
         f"intersection {index}: cubic surface is smooth")
    dzero(QY.subs({h: 0, Y: y0}),
          f"intersection {index}: residual meets graph")
    gate(dreduce(sp.diff(QY, h).subs({h: 0, Y: y0})) != 0,
         f"intersection {index}: residual germ is smooth")
    gate(dreduce(sp.diff(incidence_restriction, Y).subs(Y, y0)) != 0,
         f"intersection {index}: graph and residual are transverse")

# The residual prime is generically in the smooth locus of X0: a quadratic
# in h cannot divide the nonzero linear relative-P derivative.  The exact
# resultant supplies a hostile algebraic control.
Fp_on_ramification = dreduce(3 * h - CY)
gate(dreduce(sp.resultant(QY, Fp_on_ramification, h)) != 0,
     "residual ramification prime is generically smooth")


summary = {
    "checks": CHECKS,
    "forest": "two boundary primes cannot have two distinct smooth incidences",
    "universal": "three clean colors force at least two finite incidences",
    "cubic": "natural depressed cubic is a domain",
    "explicit": "R=1,S=Y has irreducible residual and two transverse points",
    "zmt": "a same-field Keller source would be an A2 open in the normalization",
    "scope": "extra common debt and other factor distributions remain open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3951 affine-plane boundary-incidence forest companion")
print(f"CHECKS={CHECKS}")
print("FOREST=TWO_BOUNDARY_PRIMES_CANNOT_MEET_TWICE")
print("UNIVERSAL_COLORS=0,-omega,omega^2;AT_LEAST_TWO_FINITE")
print("RAMIFICATION=(h-r)*Q;Q_IRREDUCIBLE_EQUANHARMONIC")
print("EXPLICIT_INTERSECTIONS=(0,0,omega),(0,0,-omega^2);TRANSVERSE")
print("ZMT=KELLER_SOURCE_A2_OPEN_IN_FINITE_NORMALIZATION")
print("CONCLUSION=NO_SAME_FUNCTION_FIELD_PLANAR_KELLER_CHART")
print("SCOPE=EXTRA_COMMON_DEBT_AND_OTHER_FACTOR_DISTRIBUTIONS_OPEN;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
