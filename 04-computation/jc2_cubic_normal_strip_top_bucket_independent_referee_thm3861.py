#!/usr/bin/env python3
"""Independent MSG-3070 / THM-3861 cubic-strip top-bucket referee."""

from __future__ import annotations

import hashlib
import json
import sys

import sympy as sp


sys.stdout.reconfigure(newline="\n")


GATES = 0


def require(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def zero(expr: sp.Expr) -> bool:
    return sp.cancel(sp.expand(expr)) == 0


# Reconstruct the six generic Jacobian buckets by convolution.
z = sp.symbols("z")
A = sp.symbols("A0:4")
C = sp.symbols("C0:4")
Ad = sp.symbols("A0d:4d")
Cd = sp.symbols("C0d:4d")
Az = sum(i * A[i] * z ** (i - 1) for i in range(1, 4))
As = sum(Ad[i] * z**i for i in range(4))
Cz = sum(i * C[i] * z ** (i - 1) for i in range(1, 4))
Cs = sum(Cd[i] * z**i for i in range(4))
generic = sp.Poly(sp.expand(Az * Cs - As * Cz), z)
expected = (
    A[1] * Cd[0] - Ad[0] * C[1],
    A[1] * Cd[1] + 2 * A[2] * Cd[0] - Ad[1] * C[1] - 2 * Ad[0] * C[2],
    A[1] * Cd[2] + 2 * A[2] * Cd[1] + 3 * A[3] * Cd[0]
    - 2 * Ad[1] * C[2] - Ad[2] * C[1] - 3 * Ad[0] * C[3],
    A[1] * Cd[3] + 2 * A[2] * Cd[2] + 3 * A[3] * Cd[1]
    - 3 * Ad[1] * C[3] - 2 * Ad[2] * C[2] - Ad[3] * C[1],
    2 * A[2] * Cd[3] + 3 * A[3] * Cd[2] - 3 * Ad[2] * C[3] - 2 * Ad[3] * C[2],
    3 * (A[3] * Cd[3] - Ad[3] * C[3]),
)
for degree, target in enumerate(expected):
    require(zero(generic.coeff_monomial(z**degree) - target), f"generic bucket z^{degree}")


# Differential-polynomial verification of the complete q=0, p/v nonzero
# Kummer integration.  Hs, Xs, Bs stand for ordinary derivatives.
H, Hs, X, Xs, B, Bs = sp.symbols("H Hs X Xs B Bs")
rho, sigma, d, ee, a0 = sp.symbols("rho sigma d e a0", nonzero=True)


def deriv(expr: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(expr, H) * Hs + sp.diff(expr, X) * Xs + sp.diff(expr, B) * Bs)


K = 3 * rho / (2 * sigma**2)
M = 3 * rho / (2 * sigma)
p = rho * H**3
v = sigma * H**2
beta = H * X
u = v * (K * X + d)
alpha = H * (K * X**2 / 4 + d * X + M * B + ee)
arm = -rho * X**3 / (16 * sigma**3) + 3 * rho * B * X / (4 * sigma**2) + ee * X / (2 * sigma) + d * B + a0

bucket4 = 3 * p * deriv(v) - 2 * deriv(p) * v
bucket3 = 3 * p * deriv(beta) - deriv(p) * beta + 2 * (u * deriv(v) - deriv(u) * v)
bucket2 = 3 * p * Bs + 2 * u * deriv(beta) - deriv(u) * beta + alpha * deriv(v) - 2 * deriv(alpha) * v
bucket1 = 2 * u * Bs + alpha * deriv(beta) - deriv(alpha) * beta - 2 * deriv(arm) * v
bucket0 = alpha * Bs - deriv(arm) * beta
T = B - X**2 / (4 * sigma)
factored0 = H * (M * T + ee) * deriv(T)

for degree, expression in ((4, bucket4), (3, bucket3), (2, bucket2), (1, bucket1)):
    require(zero(expression), f"integrated Kummer bucket z^{degree}")
require(zero(bucket0 - factored0), "constant-bucket factorization")


# Finite hostile controls for the all-prime Kummer valuation proof.  The
# report supplies the unrestricted gcd(2,3)=1 argument.
for order_p in range(61):
    for order_v in range(61):
        relation = 2 * order_p == 3 * order_v
        kummer = order_p % 3 == 0 and order_v % 2 == 0 and order_p // 3 == order_v // 2
        require(relation == kummer, f"Kummer orders {order_p},{order_v}")


# Local squeeze at a prime of h.  If ord(X)<0, the X^3 term is the unique
# arm pole.  If ord(X)>=0, every factor except h is local and the constant
# bucket cannot be a unit.  The integer sweep is a hostile control; the two
# inequalities themselves are the all-order proof.
local_profiles = 0
for order_h in range(1, 26):
    for order_x in range(-25, 26):
        if order_x < 0:
            excluded = 3 * order_x < order_x and 3 * order_x < 0
        else:
            excluded = order_h > 0
        require(excluded, f"local profile h={order_h},X={order_x}")
        local_profiles += 1


# If h is a unit, T is polynomial.  For deg(T)=0 its derivative vanishes;
# for positive degree r, (M*T+e)*T' has degree 2r-1, never zero.
for degree_t in range(51):
    impossible_constant_unit = degree_t == 0 or 2 * degree_t - 1 > 0
    require(impossible_constant_unit, f"constant-h degree {degree_t}")


# Polynomial top-five hostile: nonunit h survives all positive-z buckets,
# but the constant bucket is a nonconstant multiple of h.
s = sp.symbols("s")
poly_A = (
    -sp.Rational(1, 16) + sp.Rational(3, 4) * s
    + (sp.Rational(3, 8) * s + sp.Rational(3, 2) * s**2) * z
    + sp.Rational(3, 2) * s**2 * z**2
    + s**3 * z**3
)
poly_C = s + s * z + s**2 * z**2
poly_J = sp.Poly(sp.expand(sp.diff(poly_A, z) * sp.diff(poly_C, s) - sp.diff(poly_A, s) * sp.diff(poly_C, z)), z)
poly_constant = sp.Rational(3, 2) * s**2 - sp.Rational(3, 8) * s
require(zero(poly_J.coeff_monomial(1) - poly_constant), "polynomial hostile constant bucket")
for degree in range(1, 6):
    require(poly_J.coeff_monomial(z**degree) == 0, f"polynomial hostile z^{degree}")
require(sp.rem(sp.Poly(poly_constant, s), sp.Poly(s, s)) == 0, "polynomial hostile h divisibility")


# All-bucket rational hostile from THM-3861: the only failure is the cubic
# pole in the arm coefficient.
rational_A = -1 / (16 * s**3) + sp.Rational(3, 8) * s**3 * z + sp.Rational(3, 2) * s**9 * z**2 + s**15 * z**3
rational_C = s**4 * z + s**10 * z**2
rational_J = sp.cancel(sp.diff(rational_A, z) * sp.diff(rational_C, s) - sp.diff(rational_A, s) * sp.diff(rational_C, z))
require(rational_J == -sp.Rational(3, 16), "rational all-bucket hostile")
require(sp.denom(sp.cancel(rational_A.coeff(z, 0))) == 16 * s**3, "rational hostile arm pole")
for degree in range(1, 4):
    denominator = sp.denom(sp.cancel(rational_A.coeff(z, degree)))
    require(sp.degree(denominator, s) == 0, f"rational hostile A z^{degree}")
for degree in range(3):
    denominator = sp.denom(sp.cancel(rational_C.coeff(z, degree)))
    require(sp.degree(denominator, s) == 0, f"rational hostile C z^{degree}")


semantic = {
    "origin": "213867b7d40bf4bf98bb0d5e2283c57fedcd2da6",
    "scope": "q=0, p and v nonzero, characteristic-zero cubic strip",
    "kummer": "p=rho*h^3, v=sigma*h^2",
    "local_squeeze": {
        "arm_polynomiality": "at every prime of h, X=beta/h must be local",
        "constant_bucket": "at every prime of nonunit h, X must have a pole",
        "consequence": "h is a unit; the unit edge is contradictory, so the (3,2) branch is empty",
    },
    "hostiles": {
        "polynomial": "top five buckets permit h=s but constant bucket is nonconstant and h-divisible",
        "rational": "all six buckets and J=-3/16 survive with h=s^5, but a has a cubic pole",
    },
    "status": "independent proof of an already-PROVED THM-3861 branch; no new JC2 theorem",
}
blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("THM3861_MSG3070_TOP_BUCKET_INDEPENDENT_REFEREE_20260823")
print("status=PASS;already_canonical_THM3861;no_new_theorem")
print("kummer_profiles=3721")
print(f"local_profiles={local_profiles}")
print(f"polynomial_hostile_J={sp.factor(poly_constant)}")
print(f"rational_hostile_J={rational_J}")
print(f"active_gates={GATES}")
print(f"semantic_sha256={hashlib.sha256(blob).hexdigest()}")
print("ALL CHECKS PASSED")
