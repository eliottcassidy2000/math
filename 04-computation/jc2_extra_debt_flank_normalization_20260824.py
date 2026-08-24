#!/usr/bin/env python3
"""Exact scout for extra common debt in the universal THM-3950 family.

The script keeps the common multiplier symbolic for the global identities,
then specializes to the smallest symmetric extra debt c=UV.  It verifies the
normal-form identity which makes each shared-color point an A2 surface
singularity and the two residual normalization addresses which do not
separate in that normal surface.
"""

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


delta = sp.symbols("delta")


def dreduce(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.together(sp.cancel(expression)).as_numer_denom()
    num = sp.rem(
        sp.Poly(sp.expand(numerator), delta),
        sp.Poly(delta**2 + 3, delta),
    )
    den = sp.rem(
        sp.Poly(sp.expand(denominator), delta),
        sp.Poly(delta**2 + 3, delta),
    )
    return sp.cancel(num.as_expr() / den.as_expr())


def dzero(expression: sp.Expr, message: str) -> None:
    gate(dreduce(expression) == 0, message)


omega = (-1 + delta) / 2
omega2 = (-1 - delta) / 2
dzero(omega**2 + omega + 1, "omega relation")
dzero(omega - omega2 - delta, "delta relation")
dzero(delta**2 + 3, "delta square")


# ---------------------------------------------------------------------------
# Universal THM-3950 packet.  R and S are treated as independent values;
# after substitution they may be arbitrary coprime polynomials in t.
# ---------------------------------------------------------------------------

RR, SS, cc, PP, hh, xx = sp.symbols("R S c_univ P_univ h xi_univ")
UU = SS + omega2 * RR
VV = SS - omega * RR
r0_univ = dreduce(RR * UU * VV)
A0_univ = dreduce((SS - RR) * VV**2)
B0_univ = dreduce((SS + RR) * UU**2)
D0_univ = dreduce(
    (1 - omega2) * B0_univ - (1 - omega) * A0_univ
)
E0_univ = dreduce((B0_univ - A0_univ) * A0_univ * B0_univ)

dzero(
    D0_univ - (3 * r0_univ + delta * SS**3),
    "universal compressed D0 identity",
)
dzero(
    E0_univ - r0_univ**2 * (2 * r0_univ - D0_univ),
    "universal graph-root E0 identity",
)

q_univ = dreduce(cc * D0_univ * PP + cc**3 * E0_univ)
H_univ = dreduce(q_univ**2 - 4 * PP**3)
N_univ = dreduce(
    4 * PP**2
    - cc**2 * (D0_univ**2 - 4 * r0_univ**2) * PP
    + cc**4 * r0_univ**2 * (D0_univ - 2 * r0_univ) ** 2
)
dzero(
    H_univ - (cc**2 * r0_univ**2 - PP) * N_univ,
    "universal graph plus residual factorization",
)
dzero(
    N_univ.subs(PP, cc**2 * r0_univ**2)
    + 4 * delta * cc**4 * r0_univ**3 * SS**3,
    "universal graph-residual intersection divisor",
)

F_univ = dreduce(hh**3 - 3 * PP * hh + q_univ)
W_univ = dreduce(3 * hh - cc * D0_univ)
X_univ = dreduce(
    PP
    - (
        W_univ**2
        + 3 * cc * D0_univ * W_univ
        + 3 * cc**2 * D0_univ**2
    )
    / 27
)
K_univ = dreduce(D0_univ**3 / 27 + E0_univ)
dzero(
    F_univ - (-X_univ * W_univ + cc**3 * K_univ),
    "universal A-type local normal form",
)

Q0_univ = dreduce(
    -2 * xx**2
    + (D0_univ - 2 * r0_univ) * xx
    + r0_univ * (D0_univ - 2 * r0_univ)
)
disc_Q0_univ = dreduce(sp.discriminant(Q0_univ, xx))
dzero(
    disc_Q0_univ - (D0_univ - 2 * r0_univ) * (D0_univ + 6 * r0_univ),
    "universal residual discriminant",
)
dzero(
    disc_Q0_univ
    - (RR - SS)
    * (RR + SS)
    * (RR + delta * SS)
    * (3 * RR + delta * SS) ** 3
    / 3,
    "universal four-color residual discriminant",
)

# At each of the three collision colors R=0, U=0, V=0, coprimality makes S
# a unit.  The following substitutions freeze r0=E0=0, D0=delta*S^3 and
# the two distinct residual roots 0,D0/2.  The U/V substitutions are written
# in terms of S so that no division remains.
color_substitutions = {
    "R": {RR: 0},
    "U": {RR: -omega * SS},
    "V": {RR: omega2 * SS},
}
for name, substitution in color_substitutions.items():
    dzero(r0_univ.subs(substitution), f"universal {name}-color r0")
    dzero(E0_univ.subs(substitution), f"universal {name}-color E0")
    dzero(
        D0_univ.subs(substitution) - delta * SS**3,
        f"universal {name}-color D0",
    )
    dzero(
        K_univ.subs(substitution) - (delta * SS**3) ** 3 / 27,
        f"universal {name}-color K",
    )
    dzero(
        Q0_univ.subs(substitution)
        - xx * (delta * SS**3 - 2 * xx),
        f"universal {name}-color residual addresses",
    )


# ---------------------------------------------------------------------------
# Degree-one root-ratio packet, with a symbolic common multiplier c.
# ---------------------------------------------------------------------------

Y, P, Z, c = sp.symbols("Y P Z c")
U = Y + omega2
V = Y - omega
g0 = dreduce(U * V)
A0 = dreduce((Y - 1) * V**2)
B0 = dreduce((Y + 1) * U**2)
D0 = dreduce((1 - omega2) * B0 - (1 - omega) * A0)
E0 = dreduce((B0 - A0) * A0 * B0)

dzero(g0 - (Y**2 - delta * Y - 1), "minimal denominator product")
dzero(D0 - (3 * g0 + delta * Y**3), "compressed D0 identity")
dzero(E0 - g0**2 * (2 * g0 - D0), "graph-root E0 identity")

A = dreduce(c * A0)
B = dreduce(c * B0)
g = dreduce(c * g0)
q0 = dreduce(c * D0 * P + c**3 * E0)
H = dreduce(q0**2 - 4 * P**3)
N = dreduce(
    4 * P**2
    - c**2 * (D0**2 - 4 * g0**2) * P
    + c**4 * g0**2 * (D0 - 2 * g0) ** 2
)

dzero(A * B - ((Y * g) ** 2 - g**2), "assigned factors")
dzero(q0.subs(P, g**2) - 2 * g**3, "graph cusp root")
dzero(H - (g**2 - P) * N, "graph plus residual factorization")
dzero(
    N.subs(P, g**2) + 4 * delta * c**4 * g0**3 * Y**3,
    "graph-residual intersection divisor",
)

disc_N = dreduce(sp.discriminant(N, P))
dzero(
    disc_N - c**4 * (D0 - 2 * g0) ** 3 * (D0 + 6 * g0),
    "residual quadratic discriminant",
)
fixed_four_color = (1 - Y) * (1 + Y) * (1 + delta * Y) * (3 + delta * Y)
dzero(
    disc_N
    / (c**4 * (D0 - 2 * g0) ** 2 * (3 + delta * Y) ** 2 / 3)
    - fixed_four_color,
    "residual square class is independent of extra debt",
)
gate(
    dreduce(sp.discriminant(fixed_four_color, Y)) != 0,
    "four residual branch colors are distinct",
)


# ---------------------------------------------------------------------------
# The natural depressed-cubic surface and its exact local normal form.
# ---------------------------------------------------------------------------

F = dreduce(Z**3 - 3 * P * Z + q0)
gate(sp.degree(F, Z) == 3, "natural cover has rank three")
dzero(sp.discriminant(F, Z) + 27 * H, "natural cubic discriminant")

W = dreduce(3 * Z - c * D0)
X = dreduce(P - (W**2 + 3 * c * D0 * W + 3 * c**2 * D0**2) / 27)
K = dreduce(D0**3 / 27 + E0)
dzero(F - (-X * W + c**3 * K), "exact A-type surface normal form")

# The coefficient of P and the P-free term are coprime: division of
# Z^3+E0 by D0-3Z leaves K, which is not identically zero.  This freezes the
# elementary irreducibility argument for the natural cubic.
dzero(
    (Z**3 + E0).subs(Z, D0 / 3) - K,
    "linear-in-P primitive remainder",
)
gate(K != 0, "natural cubic primitive remainder is nonzero")


# ---------------------------------------------------------------------------
# Residual ramification prime.  Away from c=0 put xi=Z/c; its normalization
# is independent of c and is the same j=0 double cover as in THM-3950.
# ---------------------------------------------------------------------------

xi = sp.symbols("xi")
Q0 = dreduce(
    -2 * xi**2
    + (D0 - 2 * g0) * xi
    + g0 * (D0 - 2 * g0)
)
Qc = dreduce(
    -2 * Z**2
    + c * (D0 - 2 * g0) * Z
    + c**2 * g0 * (D0 - 2 * g0)
)
Cram = dreduce(E0 * c**3 + D0 * c * Z**2 - 2 * Z**3)

dzero(Cram - (Z - c * g0) * Qc, "ramification graph-residual split")
dzero(Qc.subs(Z, c * xi) - c**2 * Q0, "residual normalization scaling")
dzero(
    Qc * Qc.subs(Z, -Z) - N.subs(P, Z**2),
    "residual ramification prime maps onto N",
)
dzero(
    Qc.subs(Z, c * g0) - 2 * delta * c**2 * g0 * Y**3,
    "graph-residual ramification intersection divisor",
)
dzero(
    sp.discriminant(Q0, xi) - (D0 - 2 * g0) * (D0 + 6 * g0),
    "residual normalization discriminant",
)


# ---------------------------------------------------------------------------
# Exact flank values and the smallest symmetric extra debt c=UV.
# ---------------------------------------------------------------------------

alpha_U = -omega2
alpha_V = omega

flanks = {
    "U": (alpha_U, -delta, delta / 9, -xi * (delta + 2 * xi), 1),
    "V": (alpha_V, delta, -delta / 9, -xi * (-delta + 2 * xi), -1),
}

for name, (alpha, d_value, k_value, q_value, g_derivative) in flanks.items():
    dzero(g0.subs(Y, alpha), f"{name}: g0 vanishes")
    dzero(D0.subs(Y, alpha) - d_value, f"{name}: D0 is a unit")
    dzero(E0.subs(Y, alpha), f"{name}: E0 vanishes")
    dzero(K.subs(Y, alpha) - k_value, f"{name}: K is a unit")
    dzero(Q0.subs(Y, alpha) - q_value, f"{name}: two residual addresses")
    gate(
        dreduce(sp.discriminant(q_value, xi)) == -3,
        f"{name}: residual addresses are distinct",
    )
    dzero(
        sp.diff(g0, Y).subs(Y, alpha) - g_derivative,
        f"{name}: flank is simple",
    )

# At a unit multiplier, the graph root xi=g0 and the residual root through
# xi=0 have opposite nonzero slopes.  The implicit derivative of Q0 at the
# root is -g0', freezing transversality of the two ramification curves.
for name, (alpha, _, _, _, g_derivative) in flanks.items():
    q_xi = dreduce(sp.diff(Q0, xi).subs({Y: alpha, xi: 0}))
    q_y = dreduce(sp.diff(Q0, Y).subs({Y: alpha, xi: 0}))
    gate(q_xi != 0, f"{name}: residual small branch is smooth")
    dzero(-q_y / q_xi + g_derivative, f"{name}: opposite residual slope")

c_extra = g0
F_extra = dreduce(F.subs(c, c_extra))
Qc_extra = dreduce(Qc.subs(c, c_extra))
H_extra = dreduce(H.subs(c, c_extra))
N_extra = dreduce(N.subs(c, c_extra))

gate(sp.Poly(H_extra, P, Y).total_degree() == 26, "extra-debt H degree")
gate(sp.Poly(N_extra, P, Y).total_degree() == 18, "extra-debt residual degree")
dzero(
    H_extra - (c_extra**2 * g0**2 - P) * N_extra,
    "extra-debt graph plus residual",
)

# Both color factors now vanish at both flanks, with orders (1,3) and (3,1).
for name, (alpha, _, _, _, _) in flanks.items():
    dzero(A.subs(c, c_extra).subs(Y, alpha), f"{name}: A shared zero")
    dzero(B.subs(c, c_extra).subs(Y, alpha), f"{name}: B shared zero")
    dzero(F_extra.subs({Y: alpha, P: 0, Z: 0}), f"{name}: triple-root point")
    dzero(sp.diff(F_extra, P).subs({Y: alpha, P: 0, Z: 0}),
          f"{name}: surface P derivative vanishes")
    dzero(sp.diff(F_extra, Y).subs({Y: alpha, P: 0, Z: 0}),
          f"{name}: surface Y derivative vanishes")
    dzero(sp.diff(F_extra, Z).subs({Y: alpha, P: 0, Z: 0}),
          f"{name}: surface Z derivative vanishes")

# The exact coordinate identity already proves that the completion at either
# flank is XW=t^3 times a unit: c=UV has a simple zero and K is a unit.  The
# two roots displayed above map under (Y,xi)->(Y,c^2 xi^2,c xi) to the same
# point (Y,0,0), so the irreducible residual ramification curve is
# non-unibranch there.  Q0's nonsquare four-color discriminant freezes its
# irreducibility; the checks below freeze the two image coordinates.
for name, (alpha, d_value, _, _, _) in flanks.items():
    residual_roots = (sp.Integer(0), d_value / 2)
    gate(residual_roots[0] != residual_roots[1],
         f"{name}: two normalization points")
    for root in residual_roots:
        dzero(Q0.subs({Y: alpha, xi: root}),
              f"{name}: residual normalization root")
        dzero((c_extra * xi).subs({Y: alpha, xi: root}),
              f"{name}: residual Z images coalesce")
        dzero((c_extra**2 * xi**2).subs({Y: alpha, xi: root}),
              f"{name}: residual P images coalesce")


summary = {
    "checks": CHECKS,
    "family": "arbitrary coprime R,S and nonzero common multiplier c(t)",
    "coprime_debt": "THM-3951 gives two graph/residual incidences",
    "shared_color": "completed natural cubic is A_(3m-1)",
    "symmetric_minimal_extra_debt": "c=(Y+omega^2)(Y-omega)",
    "shared_zero_result": "normal surface; two residual addresses coalesce",
    "consequence": "all c excluded in the universal nonconstant-ratio packet",
    "jc2": "open; other internal-split grammars remain",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("JC2 universal extra-debt color-normalization exact scout")
print(f"CHECKS={CHECKS}")
print("F=Z^3-3PZ+cD0P+c^3E0")
print("LOCAL_NORMAL_FORM=XW=c^3(D0^3/27+E0)")
print("ORD_FINITE_R,U,V_COLOR(c)=m=>COMPLETION=A_(3m-1)")
print("GCD(c,RUV)=1=>THM3951_TWO-INCIDENCE_FOREST")
print("GCD(c,RUV)!=1=>ONE_RESIDUAL_PRIME_IS_NONUNIBRANCH")
print("c=UV=>TWO_A2_POINTS;AMBIENT_SURFACE_ALREADY_NORMAL")
print("E_RES_HAS_TWO_NORMALIZATION_ADDRESSES_AT_EACH_SHARED_ZERO")
print("MINIMAL_LITERAL_CLEAN_GATE_ESCAPE=U_OR_V;SYMMETRIC=UV")
print("OUTCOME=ALL_c_CLOSED_IN_UNIVERSAL_NONCONSTANT-RATIO_A1_PACKET")
print("SCOPE=OTHER_SPLIT_GRAMMARS_AND_JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
