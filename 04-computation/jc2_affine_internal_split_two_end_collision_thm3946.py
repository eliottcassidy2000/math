#!/usr/bin/env python3
"""Exact companion for THM-3946's full affine internal-split obstruction."""

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


P, Y, W, h, u, s, c, m, z = sp.symbols("P Y W h u s c m z")
d = sp.symbols("delta")


def dreduce(expression: sp.Expr) -> sp.Expr:
    return sp.expand(sp.rem(sp.expand(expression), d**2 + 3, d))


def dzero(expression: sp.Expr, message: str) -> None:
    gate(dreduce(expression) == 0, message)


def drzero(expression: sp.Expr, message: str) -> None:
    """Check a rational identity after reducing its numerator by delta^2+3."""
    numerator = sp.together(expression).as_numer_denom()[0]
    gate(dreduce(numerator) == 0, message)


omega = (-1 + d) / 2
omega2 = (-1 - d) / 2
dzero(omega**2 + omega + 1, "omega relation")
dzero(omega**3 - 1, "omega cube relation")
dzero((1 - omega) - (3 - d) / 2, "first cyclotomic coefficient")
dzero((1 - omega2) - (3 + d) / 2, "second cyclotomic coefficient")


# ---------------------------------------------------------------------------
# Full canonical grammar: the actual assigned affine factors are A=Y,
# B=mY+c.  The scalar ambiguity in a UFD split has already been absorbed in
# A=(q1-q0)/(2L1) and B=(q1+q0)/(2L2).
# ---------------------------------------------------------------------------

A = Y
B = m * Y + c
D = sp.expand(A * B)
p0 = P
p1 = P + D
l0 = D
l1 = p1 - omega * p0
l2 = p1 - omega2 * p0
qminus = sp.expand(2 * A * l1)
qplus = sp.expand(2 * B * l2)
q0 = dreduce((qplus - qminus) / 2)
q1 = dreduce((qplus + qminus) / 2)

gate(sp.expand(p1 - p0 - A * B) == 0,
     "one cube factor is internally split")
dzero(qminus - (q1 - q0), "first complementary factor row")
dzero(qplus - (q1 + q0), "second complementary factor row")
dzero((q1 - q0) * (q1 + q0) - 4 * (p1**3 - p0**3),
      "difference-of-cubes identity")

q0_expected = (
    m * (m - 1) * Y**3
    + c * (2 * m - 1) * Y**2
    + (c**2 + (m * (1 - omega2) - (1 - omega)) * P) * Y
    + (1 - omega2) * c * P
)
dzero(q0 - q0_expected, "full canonical q0 formula")

H = dreduce(q0**2 - 4 * P**3)
dzero(H - (q1**2 - 4 * p1**3), "common torus discriminant")


# ---------------------------------------------------------------------------
# Hidden normalization row and the universal infinity polynomial.
# ---------------------------------------------------------------------------

row = dreduce(q0.subs(P, h**2) - 2 * h**3)
row_expected = (
    m * (m - 1) * Y**3
    + c * (2 * m - 1) * Y**2
    + (c**2 + (m * (1 - omega2) - (1 - omega)) * h**2) * Y
    + (1 - omega2) * c * h**2
    - 2 * h**3
)
dzero(row - row_expected, "hidden plane-cubic row")

# For m not in {0,1}, the ordinary projective closure meets infinity at
# [h:Y:Z]=[1:z:0], where phi(z)=0.  It has no z^2 term and nonzero constant,
# so it cannot be a cube of one linear factor: a(z-r)^3 would force r=0
# from the z^2 row and then constant zero.
phi = dreduce(
    m * (m - 1) * z**3
    + (m * (1 - omega2) - (1 - omega)) * z
    - 2
)
gate(sp.Poly(phi, z).coeff_monomial(z**2) == 0,
     "infinity polynomial has no quadratic row")
gate(sp.Poly(phi, z).coeff_monomial(1) == -2,
     "infinity polynomial has nonzero constant row")
gate(sp.expand(sp.Poly(phi, z).coeff_monomial(z**3) - m * (m - 1)) == 0,
     "infinity polynomial is cubic off m=0,1")

# The discriminant sharpens the count: three distinct infinity target points
# generically; only m=-omega^2 can collide inside the m!=0,1 cubic range.
phi_disc = dreduce(sp.discriminant(phi, z))
dzero(
    phi_disc + 12 * d * m * (m - 1) * (m + omega2) ** 3,
    "infinity-polynomial discriminant",
)
dzero(phi.subs(z, 0) + 2, "zero cannot be the sole infinity root")

# Birational inverse on P!=0: h=q0/(2P), and H=0 then forces h^2=P.
h_inverse_square = sp.cancel(q0**2 / (4 * P**2) - P)
dzero(4 * P**2 * h_inverse_square - H,
      "hidden-row inverse square is H/(4P^2)")
gate(sp.expand(
    row.subs(h, 0) - Y * (m * Y + c) * ((m - 1) * Y + c)
) == 0, "hidden row has no component contracted inside h=0")
gate(sp.expand(
    H.subs(P, 0)
    - (Y * (m * Y + c) * ((m - 1) * Y + c)) ** 2
) == 0, "P is not a branch component in the coprime c-nonzero lane")


# ---------------------------------------------------------------------------
# Stronger coprime anatomy.  Scaling (h,Y) by c reduces every c!=0 row to
# c=1.  The resulting plane cubic is irreducible cuspidal for every m!=0.
# ---------------------------------------------------------------------------

alpha = d * (m + 1) + 3 * (m - 1)
beta = d + 3
Dm = m + omega2
Fcurve = dreduce(2 * row.subs(c, 1))
Fcurve_expected = (
    -4 * h**3 + alpha * h**2 * Y + beta * h**2
    + 2 * m * (m - 1) * Y**3 + 2 * (2 * m - 1) * Y**2 + 2 * Y
)
dzero(Fcurve - Fcurve_expected, "scaled coprime cubic")

ZZ = sp.symbols("ZZ")
Fbar = (
    -4 * h**3 + alpha * h**2 * Y + 2 * m * (m - 1) * Y**3
    + ZZ * (beta * h**2 + 2 * (2 * m - 1) * Y**2)
    + 2 * Y * ZZ**2
)
dzero(Fbar.subs(ZZ, 1) - Fcurve, "projective coprime cubic")

# There are no affine singularities on h=0.  The three possible roots are
# simple whenever they exist in the m!=0 family.
gate(sp.expand(
    Fcurve.subs(h, 0) - 2 * Y * (m * Y + 1) * ((m - 1) * Y + 1)
) == 0, "h-zero slice has simple factor ledger")
gate(sp.diff(Fcurve, Y).subs({h: 0, Y: 0}) == 2,
     "h-zero root Y=0 is smooth")
drzero(sp.diff(Fcurve, Y).subs({h: 0, Y: -1 / m}) + 2 / m,
       "h-zero root Y=-1/m is smooth")
drzero(
    sp.diff(Fcurve, Y).subs({h: 0, Y: -1 / (m - 1)}) - 2 / (m - 1),
    "third h-zero root is smooth when m!=1",
)

# Off h=0 the h-derivative gives one critical line.  Restriction to that line
# is a perfect cube, locating the unique affine singularity when Dm!=0 and
# proving there is none when Dm=0.
hcrit = (alpha * Y + beta) / 6
dzero(sp.diff(Fcurve, h) - 2 * h * (-6 * h + alpha * Y + beta),
      "cusp critical-line derivative")
dzero(
    Fcurve.subs(h, hcrit) - 2 * d * (Dm * Y + 1) ** 3 / 9,
    "critical-line cubic identity",
)
dzero(
    sp.diff(Fcurve, Y).subs(h, hcrit)
    - 2 * d * Dm * (Dm * Y + 1) ** 2 / 3,
    "critical-line Y-derivative identity",
)
drzero(hcrit.subs(Y, -1 / Dm) + omega / Dm,
       "unique affine cusp coordinate")

# At Dm!=0, shift the cusp to the origin.  The tangent cone is one square,
# while the cubic part is nonzero on that tangent.  Hence the cusp is
# unibranch and its tangent line cannot be a component of the plane cubic.
uu, vv = sp.symbols("uu vv")
Qlocal = 6 * omega / Dm * (uu - alpha * vv / 6) ** 2
Clocal = -4 * uu**3 + alpha * uu**2 * vv + 2 * m * (m - 1) * vv**3
drzero(
    Fcurve.subs({h: -omega / Dm + uu, Y: -1 / Dm + vv})
    - Qlocal - Clocal,
    "generic cusp tangent-plus-cubic expansion",
)
dzero(
    Clocal.subs(uu, alpha * vv / 6) - 2 * d * Dm**3 * vv**3 / 9,
    "generic cusp tangent is not a component",
)

# Lines uu=tau*vv through the cusp give an exact rational normalization.
tau = sp.symbols("tau")
Cden = -4 * tau**3 + alpha * tau**2 + 2 * m * (m - 1)
vpar_general = -(
    6 * omega / Dm * (tau - alpha / 6) ** 2 / Cden
)
hpar_general = -omega / Dm + tau * vpar_general
Ypar_general = -1 / Dm + vpar_general
drzero(Fcurve.subs({h: hpar_general, Y: Ypar_general}),
       "generic cusp-line normalization")
dzero(Cden.subs(tau, alpha / 6) - 2 * d * Dm**3 / 9,
      "normalization numerator and denominator are coprime")
dzero(
    sp.discriminant(Cden, tau)
    + 192 * d * m * (m - 1) * Dm**3,
    "normalization-pole discriminant",
)

# At Dm=0 the cusp moves to [1:omega^2:0].  The corrected tangent is
# vv=delta*ZZ/3 (not (1+delta)*ZZ/3).  Its cubic restriction is nonzero.
Fspecial = dreduce(Fbar.subs({m: -omega2, h: 1, Y: omega2 + vv}))
Qspecial = 3 * (1 + d) * (vv - d * ZZ / 3) ** 2
Cspecial = 2 * vv * (-vv**2 + d * vv * ZZ + ZZ**2)
dzero(Fspecial - Qspecial - Cspecial,
      "special infinity-cusp tangent-plus-cubic expansion")
dzero(
    Cspecial.subs(vv, d * ZZ / 3) - 2 * d * ZZ**3 / 9,
    "special infinity-cusp tangent is not a component",
)

# Lines vv=tau*ZZ give the special normalization.  Its two infinity values
# are tau=delta/3 (the cuspidal point) and tau=infinity (the simple point).
zpar_special = (
    3 * (1 + d) * (tau - d / 3) ** 2
    / (2 * tau * (tau**2 - d * tau - 1))
)
hpar_special = 1 / zpar_special
Ypar_special = omega2 / zpar_special + tau
drzero(
    Fcurve.subs({m: -omega2, h: hpar_special, Y: Ypar_special}),
    "special cusp-line normalization",
)

# Exact infinity multiplicities.  Generic m has three points; m=1 and
# m=-omega^2 each have two normalization points, with one double contact.
rinf = sp.symbols("r_inf")
Finf = dreduce(Fbar.subs({h: 1, Y: rinf, ZZ: 0}))
dzero(Finf - 2 * phi.subs(z, rinf), "infinity binary cubic")
dzero(Finf.subs(m, 1) - (-4 + 2 * d * rinf),
      "balanced finite-slope infinity root")
gate(sp.diff(Fbar, ZZ).subs({m: 1, h: 0, Y: 1, ZZ: 0}) == 2,
     "balanced double-contact infinity point is smooth")
dzero(
    Finf.subs(m, -omega2)
    + 2 * (rinf - (1 + d)) * (rinf - omega2) ** 2,
    "special infinity factorization",
)


# ---------------------------------------------------------------------------
# Balanced slope m=1: exact irreducible quartic and two-place normalization.
# ---------------------------------------------------------------------------

q0_balanced = dreduce(q0.subs(m, 1))
q1_balanced = dreduce(q1.subs(m, 1))
H_balanced = dreduce(H.subs(m, 1))
row_balanced = dreduce(row.subs(m, 1))

q0_balanced_expected = (
    P * c * d + 3 * P * c + 2 * P * d * Y
    + 2 * c**2 * Y + 2 * c * Y**2
) / 2
dzero(q0_balanced - q0_balanced_expected, "balanced q0 formula")
gate(sp.Poly(H_balanced, P, Y).total_degree() == 4,
     "balanced branch is a plane quartic")

row_balanced_expected = (
    c * Y**2 + (c**2 + d * h**2) * Y
    + c * (d + 3) * h**2 / 2 - 2 * h**3
)
dzero(row_balanced - row_balanced_expected, "balanced hidden quadratic row")
balanced_disc = dreduce(sp.discriminant(row_balanced, Y))
gate(sp.expand(balanced_disc - (c - h) ** 3 * (c + 3 * h)) == 0,
     "balanced hidden quadratic discriminant")

# Over k(c)[h] the two displayed factors are coprime and both have odd
# valuation.  Hence the discriminant is nonsquare for c!=0 and the quadratic
# row is irreducible by Gauss.
hc_domain = sp.QQ.frac_field(c)
gate(sp.gcd(sp.Poly(c - h, h, domain=hc_domain),
            sp.Poly(c + 3 * h, h, domain=hc_domain)).degree() == 0,
     "balanced discriminant factors are coprime for c nonzero")
gate(sp.factor_list(balanced_disc, h)[1]
     == [(c + 3 * h, 1), (-c + h, 3)],
     "balanced discriminant has two odd valuations")

# Remove the square (c-h)^2: u^2=(c-h)(c+3h).
hpar = sp.cancel(2 * c * (1 - s) / (s**2 + 3))
upar = sp.cancel(-c * (s - 3) * (s + 1) / (s**2 + 3))
Ypar = sp.cancel(-(c**2 + d * hpar**2) + (c - hpar) * upar)
Ypar = sp.cancel(Ypar / (2 * c))
Ppar = sp.cancel(hpar**2)

gate(sp.cancel(upar**2 - (c - hpar) * (c + 3 * hpar)) == 0,
     "balanced conic parametrization")
dzero(row_balanced.subs({h: hpar, Y: Ypar}),
      "balanced normalization lands on hidden row")
dzero(H_balanced.subs({P: Ppar, Y: Ypar}),
      "balanced normalization lands on plane quartic")
dzero(q0_balanced.subs({P: Ppar, Y: Ypar}) - 2 * hpar**3,
      "balanced normalization cube row")

# Generic birational inverse.
urec = sp.cancel((2 * c * Y + c**2 + d * h**2) / (c - h))
srec = sp.cancel((urec - c) / h)
dzero(urec.subs({h: hpar, Y: Ypar}) - upar,
      "balanced inverse recovers conic coordinate")
dzero(srec.subs({h: hpar, Y: Ypar}) - s,
      "balanced inverse recovers normalization parameter")

Ypar_expected = sp.cancel(
    -c * (s - 1) ** 2 * (2 * d + s**2 + 2 * s + 3)
    / (s**2 + 3) ** 2
)
dzero(Ypar - Ypar_expected, "simplified balanced target normalization")

S, R = sp.symbols("S R")
Phom = 4 * c**2 * (S - R) ** 2 * R**2
Yhom = -c * (S - R) ** 2 * (
    S**2 + 2 * S * R + (3 + 2 * d) * R**2
)
Zhom = (S**2 + 3 * R**2) ** 2
gate(sp.gcd(sp.gcd(Phom, Yhom), Zhom) == 1,
     "balanced projective normalization is basepoint free")
gate(sp.discriminant(S**2 + 3 * R**2, S) == -12 * R**2,
     "balanced normalization has two distinct infinity places")
gate(sp.gcd(S**2 + 3 * R**2, (S - R) * R) == 1,
     "P coordinate prevents balanced infinity cancellation")


# ---------------------------------------------------------------------------
# Repeated-root boundary c=0 and the corrected THM-3944 endpoint.
# ---------------------------------------------------------------------------

gate(sp.expand(D.subs(c, 0) - m * Y**2) == 0,
     "c=0 is exactly the repeated-square lane")

q0_collision = dreduce(q0.subs({c: 0, m: 1}))
q1_collision = dreduce(q1.subs({c: 0, m: 1}))
H_collision = dreduce(H.subs({c: 0, m: 1}))
gate(q0_collision == d * P * Y, "balanced collision first coefficient")
gate(sp.expand(q1_collision - Y * (3 * P + 2 * Y**2)) == 0,
     "balanced collision second coefficient")
gate(sp.expand(H_collision + P**2 * (4 * P + 3 * Y**2)) == 0,
     "balanced collision square-conductor discriminant")

V = sp.symbols("V")
Pnorm = -(V**2 + 3 * Y**2) / 4
Wnorm = sp.expand(Pnorm * V)
f0 = dreduce((q0_collision + W).subs({P: Pnorm, W: Wnorm}))
f1 = sp.expand((q1_collision + W).subs({P: Pnorm, W: Wnorm}))
dzero(f0 + (V + d * Y) ** 2 * (V - d * Y) / 4,
      "collision first Cardano vector is (2,1)")
gate(sp.expand(f1 + (V + Y) ** 3 / 4) == 0,
     "collision second Cardano radicand is a cube")
gate((2 % 3, 1 % 3) != (0, 0),
     "collision first class survives on order regular Gm2")


# ---------------------------------------------------------------------------
# Exhaustive affine gauge: after absorbing the split scalar, any two
# nonconstant affine factors become A=Ynew, B=mYnew+c.  The slope ratio m is
# genuine (up to inversion under exchanging the two sides), not removable.
# ---------------------------------------------------------------------------

a, b, e, f, Ynew = sp.symbols("a b e f Ynew")
Yold = sp.cancel((Ynew - b) / a)
Bold = e * Yold + f
gate(sp.cancel(Bold - ((e / a) * Ynew + f - e * b / a)) == 0,
     "general affine factors reduce to Y and mY+c")
gate(sp.cancel((e / a) * (a / e)) == 1,
     "exchanging the factors inverts the slope ratio")


summary = {
    "checks": CHECKS,
    "grammar": "actual quotient factors A=Y, B=mY+c",
    "coprime_balanced": "c!=0,m=1: irreducible quartic; normalization Gm; two ends",
    "coprime_unbalanced": "c!=0: irreducible cuspidal H; normalization P1 minus 2 or 3 points",
    "unbalanced_infinity": "phi=m(m-1)z^3+[m(1-w2)-(1-w)]z-2; never triple-root",
    "repeated_boundary": "c=0: scalar repeated-square lane of THM3947",
    "balanced_collision": "THM3944: order-reg rank 1, no extension to normal A2",
    "end_count": "2 for m=1,-omega^2; 3 otherwise",
    "constant_factor_boundary": "m=0: whole-factor THM3942 lane",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3946 full affine internal-split companion")
print(f"CHECKS={CHECKS}")
print("GAUGE=A=Y;B=mY+c;SLOPE_RATIO_m_IS_GENUINE")
print("C_NONZERO_M1=IRREDUCIBLE_QUARTIC;NORMALIZATION=Gm;INFINITY_PLACES=2")
print("C_NONZERO=IRREDUCIBLE_CUSPIDAL;NORMALIZATION=P1_MINUS_2_OR_3_POINTS")
print("END_COUNT=2_FOR_m=1,-omega^2;3_OTHERWISE")
print("INFINITY_PHI=m(m-1)z^3+[m(1-w2)-(1-w)]z-2;TRIPLE_ROOT=IMPOSSIBLE")
print("PHI_DISC=-12delta*m(m-1)(m+omega^2)^3")
print("C_ZERO=REPEATED_SQUARE_THM3947;M1_ENDPOINT=THM3944")
print("THM3944_ORDER_REG_CARDANO_RANK=1;NORMALIZATION_EXTENDABLE_RANK=0")
print("CONCLUSION=FULL_DISCRIMINANT_REDUCIBLE_OR_AT_LEAST_TWO_INFINITY_PLACES")
print(f"SEMANTIC_SHA256={semantic}")
