#!/usr/bin/env python3
"""Exact companion for THM-3820's universal residual-discriminant identity."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


u, e, G, D = sp.symbols("u e G D")
Y, Z, t = sp.symbols("Y Z t")

# The coefficient-free two-equation compression used by the pure-r nodal
# carrier papers.  The first equation is the z-Hamiltonian after r=u/e; the
# second is the compatibility of z^2=(1+2u)/(9G) and z^3=u(1+u)/e.
K_source = 1 + 2 * u
P = sp.expand(G * u**2 - K_source * (2 * e**3 + u * e * D))
Q = sp.expand(e**2 * K_source**3 - 729 * G**3 * u**2 * (1 + u) ** 2)

resultant = sp.factor(sp.resultant(P, Q, u))
H = sp.factor(resultant / (G**3 * e**4))
gate(sp.Poly(H, e, G, D).is_zero is False, "nonzero universal residual")
zero(resultant - G**3 * e**4 * H, "universal resultant factorization")

H_expected = (
    -729 * D**4 * e**2
    + 1458 * D**3 * G * e
    + 8748 * D**3 * e**4
    + 2125764 * D**2 * G**3 * e**4
    - 729 * D**2 * G**2
    - 17496 * D**2 * G * e**3
    - 34992 * D**2 * e**6
    - 4251528 * D * G**4 * e**3
    - 8503056 * D * G**3 * e**6
    + 11664 * D * G**2 * e**2
    + 52488 * D * G * e**5
    + 46656 * D * e**8
    - 2 * D * e
    + 2125764 * G**5 * e**2
    + 8503056 * G**4 * e**5
    + 8503056 * G**3 * e**8
    - 2916 * G**3 * e
    - 17496 * G**2 * e**4
    - 23328 * G * e**7
    + G
)
zero(H - H_expected, "explicit universal residual")
gate(len(sp.Poly(H, e, G, D).terms()) == 20, "twenty universal residual terms")

# Euler normalization.  Equivalently G=e^3Y and eD-G=e^3Z, so
# D=e^2(Y+Z).  Every resulting e-exponent is 3 modulo 7.
euler_substitution = {G: e**3 * Y, D: e**2 * (Y + Z)}
H_euler = sp.Poly(sp.expand(H.subs(euler_substitution)), e)
gate(
    {degree[0] % 7 for degree, _ in H_euler.terms()} == {3},
    "one mod-seven Euler residue class",
)
gate({degree[0] for degree, _ in H_euler.terms()} == {3, 10, 17},
     "Euler residual has exactly three layers")

K = Z * (Y + Z) ** 2 + 4 * (2 - Z) * (Y + 2 * Z)
R = (Y + Z) ** 2 - 8 * Z + 16
W = (Y + Z) * (Z - 4) ** 2 + 4 * Y * (Z - 2)
F = (
    -(Y + 2 * Z)
    - 729 * t * (Z - 4) * K
    + 2125764 * t**2 * Y**3 * (Z - 2) ** 2
)
zero(H_euler.as_expr() - e**3 * F.subs(t, e**7),
     "universal Euler/mod-seven quadratic")
gate(sp.Poly(F, t).degree() == 2, "generic residual degree two in t")
zero(sp.Poly(F, t).LC() - 4 * 3**12 * Y**3 * (Z - 2) ** 2,
     "quadratic leading coefficient")

# The residual discriminant is the source discriminant times one square.
disc_t = sp.factor(sp.discriminant(F, t))
zero(
    (Z - 4) ** 2 * K**2
    + 16 * Y**3 * (Z - 2) ** 2 * (Y + 2 * Z)
    - R * W**2,
    "core discriminant square identity",
)
zero(disc_t - 3**12 * R * W**2, "residual discriminant factorization")

p = sp.factor(P.subs(euler_substitution) / e**3)
p_expected = -(Y + 2 * Z) * u**2 - (4 + Y + Z) * u - 2
zero(p - p_expected, "Euler-normalized source quadratic")
disc_u = sp.factor(sp.discriminant(p, u))
zero(disc_u - R, "normalized source discriminant")
disc_P = sp.factor(sp.discriminant(P, u))
zero(disc_P - e**2 * ((D - 4 * e**2) ** 2 + 8 * e * G),
     "unnormalized source discriminant")
zero(
    disc_P.subs(euler_substitution) - e**6 * R,
    "Euler-normalized unscaled source discriminant",
)
zero(
    disc_t - 3**12 * W**2 * disc_P.subs(euler_substitution) / e**6,
    "source-to-residual discriminant identity",
)

# Direct elimination after normalization.  This separately checks that F is
# the two-valued image of the source quadratic under the rational t-map.
q_normalized = (1 + 2 * u) ** 3 - 729 * t * Y**3 * u**2 * (1 + u) ** 2
zero(Q.subs(euler_substitution) / e**2 - q_normalized.subs(t, e**7),
     "normalized compatibility equation")
zero(sp.resultant(p, q_normalized, u) - Y**3 * F,
     "normalized direct resultant")

# In the generic source quadratic algebra, reduce
# phi(u)=(1+2u)^3/(729Y^3u^2(1+u)^2).  Its affine u-coefficient is the
# collision factor W.  Hence W=0 identifies two distinct source sheets even
# when R is nonzero; R=0 instead means the source sheets themselves collide.
fraction_field = sp.QQ.frac_field(Y, Z)
p_poly = sp.Poly(p, u, domain=fraction_field)
denominator = 729 * Y**3 * u**2 * (1 + u) ** 2
inverse = sp.invert(sp.Poly(denominator, u, domain=fraction_field), p_poly)
phi_remainder = sp.rem(
    sp.Poly(sp.expand((1 + 2 * u) ** 3 * inverse.as_expr()),
            u, domain=fraction_field),
    p_poly,
).as_expr()
gate(sp.Poly(phi_remainder, u).degree() == 1, "generic t-map is affine on source algebra")
phi_slope = sp.factor(sp.Poly(phi_remainder, u).coeff_monomial(u))
zero(
    phi_slope - (Y + 2 * Z) * W / (2916 * Y**3 * (Z - 2) ** 2),
    "source-sheet collision slope",
)
zero(
    sp.rem(
        sp.Poly(sp.together(denominator * phi_remainder - (1 + 2 * u) ** 3),
                u, domain=fraction_field),
        p_poly,
    ).as_expr(),
    "affine t-map reduction",
)

# Sharp boundary/hostile controls.
zero(F.subs(Y, 0) + Z * (2 + 729 * t * (Z - 4) ** 3),
     "G=0 arm continuation")
gate(sp.Poly(F.subs(Z, 2), t).degree() == 1,
     "Z=2 sends one residual value to infinity")
zero(F.subs(Z, 2) - (-(Y + 4) + 2916 * t * (Y + 2) ** 2),
     "Z=2 exact degree-drop formula")

collision_point = {Y: -sp.Rational(9, 5), Z: 1}
zero(W.subs(collision_point), "finite W collision control")
gate(R.subs(collision_point) == sp.Rational(216, 25),
     "W collision has distinct source sheets")
zero(disc_t.subs(collision_point), "W collision residual double root")

source_double_point = {Y: -sp.Rational(1, 2), Z: sp.Rational(5, 2)}
zero(R.subs(source_double_point), "finite source double-root control")
gate(W.subs(source_double_point) != 0, "source double root avoids W collision")
zero(disc_t.subs(source_double_point), "source double root residual double root")

dense_point = {Y: 2, Z: 3}
gate(R.subs(dense_point) == 17, "dense source discriminant control")
gate(W.subs(dense_point) == 13, "dense collision-factor control")
gate(disc_t.subs(dense_point) == 3**12 * 17 * 13**2,
     "dense residual discriminant control")

# The original G=0 formula cross-checks the arm continuation before Euler
# normalization.  It is not used to divide by G.
zero(
    H.subs(G, 0) - D * e * (-729 * e * (D - 4 * e**2) ** 3 - 2),
    "unnormalized G=0 arm law",
)

identity_packet = "\n".join(
    sp.srepr(sp.expand(item)) for item in (H, F, K, R, W, disc_t, p, phi_slope)
).encode()
semantic = {
    "source": "P=G*u^2-(1+2u)(2e^3+u*e*D); Q=e^2(1+2u)^3-729G^3u^2(1+u)^2",
    "resultant": "Res_u(P,Q)=G^3*e^4*H",
    "euler": "G=e^3Y;D=e^2(Y+Z);t=e^7;H=e^3F",
    "quadratic": "F=-(Y+2Z)-729t(Z-4)K+2125764t^2Y^3(Z-2)^2",
    "discriminant": "Disc_t(F)=3^12*R*W^2;Disc_u(P/e^3)=R",
    "collision": "W is the affine slope of the t-map on the generic source quadratic algebra",
    "degree_drop": "Y=0 is divided-arm continuation;Z=2 drops the t-quadratic to a line",
    "open": "does not close degree>=6 profiles or mixed carrier corrections",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3820-universal-euler-mod-seven-residual-quadratic-discriminant")
print("universal_resultant=Res_u(P,Q)=G^3*e^4*H")
print("euler_chart=G:e^3Y;D:e^2(Y+Z);t:e^7")
print("layers=H_after_euler_has_e_degrees_3,10,17")
print("F=-(Y+2Z)-729*t*(Z-4)*K+2125764*t^2*Y^3*(Z-2)^2")
print("K=Z*(Y+Z)^2+4*(2-Z)*(Y+2Z)")
print("R=(Y+Z)^2-8Z+16")
print("W=(Y+Z)*(Z-4)^2+4Y*(Z-2)")
print("discriminant=Disc_t(F)=3^12*R*W^2")
print("source=Disc_u(P/e^3)=R;Disc_u(P)=e^6*R_after_euler")
print("collision=W_is_t_map_slope_up_to_a_unit_on_generic_source_quadratic")
print("boundaries=Y0_arm_continuation;Z2_degree_drop;W0_sheet_identification;R0_source_collision")
print("open=degree_ge_6_pure_r;mixed_carriers;Darboux_pair")
print(f"identity_packet_sha256={hashlib.sha256(identity_packet).hexdigest()}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
