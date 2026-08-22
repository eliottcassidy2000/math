#!/usr/bin/env python3
"""Fast exact structural checks for the level-three Keller norm divisor.

This companion does not expand the global numerator J.  It verifies the
quotient-algebra formulas and the Newton-face calculation which prove the
order-seven pole of Norm(H) along L.  It also supplies a hostile finite-sheet
point.  All arithmetic is exact over QQ; ``require`` remains active under
``python -O``.
"""

from __future__ import annotations

import hashlib
import pickle
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
H_PATH = ROOT / "05-knowledge/results/keller_L2_polynomial_opus_20260728.pkl"
H_RAW_SHA256 = "5a9459b3149e500c1b00b67bd804aa7e607de06bf4610c7cdf5fa26d41d74ce9"
H_LEDGER_SHA256 = "b146c11f33e895c08f72303d282e2668d955e0a58d9268a1b445d4d5202016c2"

raw = H_PATH.read_bytes()
require(hashlib.sha256(raw).hexdigest() == H_RAW_SHA256, "transported H pickle changed")
H = pickle.loads(raw)

a, b, c, w, u = sp.symbols("a b c w u")
H_poly = sp.Poly(H, a, b, c)
H_ledger = "\n".join(
    f"{monomial}:{coefficient}" for monomial, coefficient in H_poly.terms()
)
require(
    hashlib.sha256(H_ledger.encode("ascii")).hexdigest() == H_LEDGER_SHA256,
    "transported H coefficient ledger changed",
)

L = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
T = 4 - 3 * b * c
S = 27 * a * c**2 - 9 * b * c + 8
K = 9 * a * c - b
D = 18 * a * c - 3 * b**2 * c + 2 * b
M = 27 * a * c**2 - 9 * b * c + 26

Y0 = 81 * a * b * c**2 - 72 * a * c - 15 * b**2 * c + 16 * b
A1 = 27 * a * b * c**2 + 54 * a * c - 9 * b**2 * c + 2 * b
A2 = 27 * a * b**2 * c**2 + 18 * a * b * c - 48 * a - 9 * b**3 * c + 10 * b**2
Z0 = (
    -2916 * a**3 * c**4
    + 2916 * a**2 * b * c**3
    - 4536 * a**2 * c**2
    + 621 * a * b**3 * c**3
    - 1026 * a * b**2 * c**2
    + 504 * a * b * c
    + 64 * a
    - 207 * b**4 * c**2
    + 454 * b**3 * c
    - 256 * b**2
)

# Reduced inverse coordinates in QQ(a,b,c)[w]/(E):
#     2S*q_y=Y,  8S*q_z=Z.
Y = Y0 + 6 * L * w - 3 * K * L * w**2
Z = Z0 + 6 * L * A1 * w - 9 * L * A2 * w**2
E = L * w**3 + T * w - 2 * c

# Re-derive those compact formulas from the rational inverse chart rather
# than merely checking their internal factorization.
Dchart = (12 * a - b**2) * w**2 + b * w + 2
y_numerator = 3 * a * w * (K * w + 2)
qy_numerator = b * Dchart - y_numerator
qz_numerator = (2 * w - c) * Dchart - 3 * w**2 * qy_numerator
coefficient_field = sp.QQ.frac_field(a, b, c)


def quotient_remainder(expression) -> sp.Poly:
    return sp.Poly(expression, w, domain=coefficient_field).rem(
        sp.Poly(E, w, domain=coefficient_field)
    )


require(
    quotient_remainder(Dchart * Y - 2 * S * qy_numerator).is_zero,
    "the reduced q_y formula does not match the inverse chart",
)
require(
    quotient_remainder(8 * S * qz_numerator - Z * w**3 * Dchart).is_zero,
    "the reduced q_z formula does not match the inverse chart",
)

# These identities expose the first nonzero u=1/w terms on the two divergent
# Newton branches.  The cubic gives L=-T*u^2+2c*u^3 there.
require(sp.expand(Y0 + 3 * K * T - 2 * D) == 0, "Y leading identity failed")
require(
    sp.expand(-6 * A1 * T - 18 * A2 * c + 24 * D) == 0,
    "Z linear identity failed",
)
require(
    sp.expand(Z0 + 9 * A2 * T + 4 * M * L) == 0,
    "Z constant cancellation failed",
)

L_on_branch = -T * u**2 + 2 * c * u**3
Y_on_branch = sp.expand(Y0 + 6 * L_on_branch / u - 3 * K * L_on_branch / u**2)
Z_on_branch = sp.expand(
    -9 * A2 * T
    - 4 * M * L_on_branch
    + 6 * A1 * L_on_branch / u
    - 9 * A2 * L_on_branch / u**2
)
expected_Y = 2 * D - 6 * (T + K * c) * u + 12 * c * u**2
expected_Z = -24 * D * u + (12 * A1 * c + 4 * M * T) * u**2 - 8 * M * c * u**3
require(sp.expand(Y_on_branch - expected_Y) == 0, "u-expansion of q_y failed")
require(sp.expand(Z_on_branch - expected_Z) == 0, "u-expansion of q_z failed")

# The H-face for the valuation (v(w),v(q_y),v(q_z))=(-1/2,0,+1/2).
face_weight = max(i - k for (i, _j, k), _coefficient in H_poly.terms())
face = sum(
    coefficient * a**i * b**j * c**k
    for (i, j, k), coefficient in H_poly.terms()
    if i - k == face_weight
)
expected_face = -63078912 * a**7 * (3 * a * c - 2 * b) ** 3
require(face_weight == 7, "the transported H Newton face has the wrong weight")
require(sp.expand(face - expected_face) == 0, "the transported H Newton face changed")

# One point on L=0 makes all generic-DVR unit hypotheses simultaneous and
# also controls the finite Newton root.  Here the finite root is w=2.
hostile_target = {a: sp.Rational(2, 27), b: sp.Integer(1), c: sp.Integer(1)}
require(L.subs(hostile_target) == 0, "hostile target is not on L")
require(
    (
        c.subs(hostile_target),
        T.subs(hostile_target),
        S.subs(hostile_target),
        D.subs(hostile_target),
    )
    == (1, 1, 1, sp.Rational(1, 3)),
    "c,T,S,D unit witness changed",
)
hostile_w = sp.Integer(2)
hostile_Dchart = Dchart.subs({**hostile_target, w: hostile_w})
hostile_y = sp.cancel(qy_numerator.subs({**hostile_target, w: hostile_w}) / hostile_Dchart)
hostile_z = sp.cancel(qz_numerator.subs({**hostile_target, w: hostile_w}) / (hostile_w**3 * hostile_Dchart))
hostile_q = (hostile_w, hostile_y, hostile_z)
require(
    hostile_q == (sp.Integer(2), sp.Rational(5, 6), sp.Rational(-7, 8)),
    "finite inverse point changed",
)

xq, yq, zq = hostile_q
xy = xq * yq
image_q = (
    (1 + xy) ** 3 * zq + yq**2 * (1 + xy) * (4 + 3 * xy),
    yq + 3 * xq * (1 + xy) ** 2 * zq + 3 * xq * yq**2 * (4 + 3 * xy),
    2 * xq - 3 * xq**2 * yq - xq**3 * zq,
)
require(image_q == tuple(hostile_target[symbol] for symbol in (a, b, c)), "finite inverse point failed F(q)=target")
hostile_H = sp.cancel(H.subs({a: xq, b: yq, c: zq}))
expected_hostile_H = sp.Rational(
    3393794313700733412883215882425216567,
    359414999291950792704,
)
require(hostile_H == expected_hostile_H and hostile_H != 0, "finite-sheet H control changed")

print("== exact level-three norm-divisor structure ==")
print(f"transported H raw/ledger sha256={H_RAW_SHA256}/{H_LEDGER_SHA256}")
print("inverse-chart quotient formulas: PASS")
print("generic DVR at (L): c,T,S,D are units (simultaneous witness (2/27,1,1))")
print("Newton roots of E: one v_L(w)=0 and two v_L(w)=-1/2")
print("divergent sheets: q_y~D/S, q_z~-3(D/S)u, u=1/w")
print("H face: max(i-k)=7; in(H)=-63078912*x^7*(3*x*z-2*y)^3")
print("divergent sheets: v_L(H(q))=-7/2 each; leading 3*w*q_z-2*q_y=-11D/S")
print(f"finite sheet: q={hostile_q}; H(q)={hostile_H} != 0")
print("therefore v_L(Norm(H))=-7")
print("scope: structural pole proof only; global numerator and JC implications are not inferred here")
print("all exact checks passed")
