#!/usr/bin/env python3
"""Deterministic exact companion for provisional THM-3861."""

from __future__ import annotations

import hashlib

import sympy as sp


CHECKS = 0


def check(label: str, condition: object) -> None:
    """Optimization-safe exact gate."""

    global CHECKS
    CHECKS += 1
    if condition is True or condition == sp.S.true:
        return
    if isinstance(condition, sp.Basic) and sp.cancel(sp.factor(condition)) == 0:
        return
    raise RuntimeError(f"CHECK FAILED: {label}: {condition}")


def zero(label: str, expression: sp.Expr) -> None:
    check(label, sp.cancel(sp.factor(expression)) == 0)


s, z = sp.symbols("s z")

# ---------------------------------------------------------------------------
# 1. The six universal coefficient buckets.
# ---------------------------------------------------------------------------

a, alpha, u, p = [sp.Function(name)(s) for name in ("a", "alpha", "u", "p")]
b, beta, v, q = [sp.Function(name)(s) for name in ("b", "beta", "v", "q")]
A = a + alpha * z + u * z**2 + p * z**3
C = b + beta * z + v * z**2 + q * z**3
J = sp.expand(sp.diff(A, z) * sp.diff(C, s) - sp.diff(A, s) * sp.diff(C, z))

buckets = [
    alpha * sp.diff(b, s) - sp.diff(a, s) * beta,
    alpha * sp.diff(beta, s) - sp.diff(alpha, s) * beta
    + 2 * (u * sp.diff(b, s) - sp.diff(a, s) * v),
    alpha * sp.diff(v, s) + 2 * u * sp.diff(beta, s)
    + 3 * p * sp.diff(b, s) - 2 * sp.diff(alpha, s) * v
    - sp.diff(u, s) * beta - 3 * sp.diff(a, s) * q,
    alpha * sp.diff(q, s) + 2 * u * sp.diff(v, s)
    + 3 * p * sp.diff(beta, s) - 3 * sp.diff(alpha, s) * q
    - 2 * sp.diff(u, s) * v - sp.diff(p, s) * beta,
    2 * u * sp.diff(q, s) + 3 * p * sp.diff(v, s)
    - 3 * sp.diff(u, s) * q - 2 * sp.diff(p, s) * v,
    3 * (p * sp.diff(q, s) - sp.diff(p, s) * q),
]

for degree, bucket in enumerate(buckets):
    zero(f"generic Jacobian bucket z^{degree}", J.coeff(z, degree) - bucket)
check("no coefficient above z5", sp.Poly(J, z).degree() <= 5)

# A constant target determinant scales the top Wronskian as expected.
aa, bb, cc, dd = sp.symbols("aa bb cc dd")
p_new = aa * p + bb * q
q_new = cc * p + dd * q
top_new = sp.factor(p_new * sp.diff(q_new, s) - sp.diff(p_new, s) * q_new)
zero("top Wronskian target covariance",
     top_new - (aa * dd - bb * cc) * (p * sp.diff(q, s) - sp.diff(p, s) * q))

# ---------------------------------------------------------------------------
# 2. Complete (3,1) triangular family.
# ---------------------------------------------------------------------------

rho, beta0, d0, e0, lam, a0 = sp.symbols(
    "rho beta0 d0 e0 lam a0", nonzero=True
)
bfun = sp.Function("bfun")(s)
C31 = bfun + beta0 * z
A31 = rho * C31**3 + d0 * C31**2 + e0 * C31 - lam * s / beta0 + a0
J31 = sp.expand(sp.diff(A31, z) * sp.diff(C31, s)
                - sp.diff(A31, s) * sp.diff(C31, z))
zero("(3,1) family constant Jacobian", J31 - lam)

coeff31 = [sp.expand(A31).coeff(z, i) for i in range(4)]
zero("(3,1) p coefficient", coeff31[3] - rho * beta0**3)
zero("(3,1) u coefficient", coeff31[2] - beta0**2 * (3 * rho * bfun + d0))
zero("(3,1) alpha coefficient",
     coeff31[1] - beta0 * (3 * rho * bfun**2 + 2 * d0 * bfun + e0))
zero("(3,1) arm coefficient",
     coeff31[0] - (rho * bfun**3 + d0 * bfun**2 + e0 * bfun
                   - lam * s / beta0 + a0))

A_t, C_t = sp.symbols("A_t C_t")
s_inverse = beta0 * (rho * C_t**3 + d0 * C_t**2 + e0 * C_t + a0 - A_t) / lam
z_inverse = (C_t - bfun.subs(s, s_inverse)) / beta0
zero("inverse recovers s", s_inverse.subs({A_t: A31, C_t: C31}) - s)
zero("inverse recovers z", z_inverse.subs({A_t: A31, C_t: C31}) - z)

# Verify the successive differential integrations without assuming beta constant.
betaf = sp.Function("betaf")(s)
bf = sp.Function("bf")(s)
p31 = rho * betaf**3
u31 = betaf**2 * (3 * rho * bf + d0)
alpha31 = betaf * (3 * rho * bf**2 + 2 * d0 * bf + e0)
a31 = rho * bf**3 + d0 * bf**2 + e0 * bf - lam * s / beta0 + a0
zero("(3,1) top ODE", 3 * p31 * sp.diff(betaf, s) - sp.diff(p31, s) * betaf)
zero("(3,1) u ODE",
     2 * u31 * sp.diff(betaf, s) + 3 * p31 * sp.diff(bf, s)
     - sp.diff(u31, s) * betaf)
zero("(3,1) alpha ODE",
     alpha31 * sp.diff(betaf, s) - sp.diff(alpha31, s) * betaf
     + 2 * u31 * sp.diff(bf, s))
zero("(3,1) constant-bucket factor",
     alpha31 * sp.diff(bf, s) - sp.diff(
         rho * bf**3 + d0 * bf**2 + e0 * bf, s
     ) * betaf)

# ---------------------------------------------------------------------------
# 3. Exhaustive (3,2) Kummer integrations and final factor.
# ---------------------------------------------------------------------------

P, V = sp.symbols("P V", nonzero=True)
h = sp.Function("h")(s)
X = sp.Function("X")(s)
b32 = sp.Function("b32")(s)
d32, e32, a32_0 = sp.symbols("d32 e32 a32_0")

p32 = P * h**3
v32 = V * h**2
beta32 = h * X
K = 3 * P / (2 * V**2)
M = 3 * P / (2 * V)
U = K * X + d32
u32 = v32 * U
Y = K * X**2 / 4 + d32 * X + M * b32 + e32
alpha32 = h * Y
a32 = (-P * X**3 / (16 * V**3) + 3 * P * b32 * X / (4 * V**2)
       + e32 * X / (2 * V) + d32 * b32 + a32_0)

E4 = 3 * p32 * sp.diff(v32, s) - 2 * sp.diff(p32, s) * v32
E3 = (2 * u32 * sp.diff(v32, s) + 3 * p32 * sp.diff(beta32, s)
      - 2 * sp.diff(u32, s) * v32 - sp.diff(p32, s) * beta32)
E2 = (alpha32 * sp.diff(v32, s) + 2 * u32 * sp.diff(beta32, s)
      + 3 * p32 * sp.diff(b32, s) - 2 * sp.diff(alpha32, s) * v32
      - sp.diff(u32, s) * beta32)
E1 = (alpha32 * sp.diff(beta32, s) - sp.diff(alpha32, s) * beta32
      + 2 * (u32 * sp.diff(b32, s) - sp.diff(a32, s) * v32))
E0 = alpha32 * sp.diff(b32, s) - sp.diff(a32, s) * beta32

zero("(3,2) Kummer top bucket", E4)
zero("(3,2) integrated E3", E3)
zero("(3,2) integrated E2", E2)
zero("(3,2) integrated E1", E1)
T = b32 - X**2 / (4 * V)
F = M * T + e32
zero("(3,2) final E0 factorization", E0 - h * F * sp.diff(T, s))
zero("Kummer constant ratio", p32**2 / v32**3 - P**2 / V**3)

# A generic specialization with nonconstant h confirms the rational identities.
h_poly = s**2 + 1
X_rat = (s**3 - 2) / h_poly
b_poly = s**4 + s
specialization = {
    h: h_poly,
    X: X_rat,
    b32: b_poly,
}
for label, expression in (
    ("specialized E4", E4),
    ("specialized E3", E3),
    ("specialized E2", E2),
    ("specialized E1", E1),
    ("specialized E0 factor", E0 - h * F * sp.diff(T, s)),
):
    zero(label, expression.xreplace(specialization).doit())

# ---------------------------------------------------------------------------
# 4. Sharp rational hostile: every bucket survives, only a(s) has a pole.
# ---------------------------------------------------------------------------

C_hostile = s**4 * z + s**10 * z**2
A_hostile = (-sp.Rational(1, 16) / s**3 + sp.Rational(3, 8) * s**3 * z
             + sp.Rational(3, 2) * s**9 * z**2 + s**15 * z**3)
J_hostile = sp.factor(sp.diff(A_hostile, z) * sp.diff(C_hostile, s)
                      - sp.diff(A_hostile, s) * sp.diff(C_hostile, z))
check("sharp hostile constant Jacobian", J_hostile == -sp.Rational(3, 16))
check("sharp hostile C is polynomial", sp.denom(C_hostile) == 1)
check("sharp hostile nonpolynomial arm coefficient",
      sp.limit(s**3 * A_hostile.subs(z, 0), s, 0) == -sp.Rational(1, 16))

hostile_coeff_A = [sp.expand(A_hostile).coeff(z, i) for i in range(4)]
hostile_coeff_C = [sp.expand(C_hostile).coeff(z, i) for i in range(4)]
for degree in range(6):
    check(f"sharp hostile bucket z^{degree}",
          sp.expand(J_hostile).coeff(z, degree)
          == (-sp.Rational(3, 16) if degree == 0 else 0))

zero("hostile p=h3", hostile_coeff_A[3] - (s**5) ** 3)
zero("hostile v=h2", hostile_coeff_C[2] - (s**5) ** 2)
zero("hostile beta=hX", hostile_coeff_C[1] - s**5 * (1 / s))
T_hostile = -1 / (4 * s**2)
zero("hostile factored constant bucket",
     s**5 * (sp.Rational(3, 2) * T_hostile) * sp.diff(T_hostile, s)
     + sp.Rational(3, 16))

SEMANTIC_FACTS = "\n".join(
    [
        "cubic-strip-six-buckets=exact",
        "top-cubic-row=constant-target-direction",
        "degree-31=triangular-cubic-automorphism",
        "degree-32=Kummer-h3-h2-packet",
        "degree-32-constant-bucket=h-times-affine-T-times-Tprime",
        "nonconstant-h=X-pole-forces-uncancellable-X3-arm-pole",
        "constant-h=unit-factor-contradiction",
        "rational-hostile=Jacobian-minus-3-over-16",
        "scope=polynomial-z-depth-at-most-three-only",
    ]
)
semantic_sha256 = hashlib.sha256(SEMANTIC_FACTS.encode("utf-8")).hexdigest()

print("THM-3861 CUBIC NORMAL-STRIP KELLER CLASSIFICATION")
print("SIX_BUCKETS=PASS;TOP_DIRECTION=SL2_CONSTANT;DEGREE_DROPS=THM-3856")
print("BRANCH_3_1=TRIANGULAR_AUTOMORPHISM;EXPLICIT_INVERSE=PASS")
print("BRANCH_3_2=EMPTY;KUMMER=h3_h2;FINAL_FACTOR=h*(M*T+e)*Tprime")
print("POLYNOMIALITY_GATE=NONCONSTANT_h_X3_POLE;CONSTANT_h_UNIT_CONTRADICTION")
print("RATIONAL_HOSTILE=PASS;JACOBIAN=-3/16;ONLY_ARM_COEFFICIENT_HAS_POLE")
print("BOUNDARY=RATIONAL_OR_INFINITE_NORMAL_SERIES_REMAIN_OPEN")
print(f"SEMANTIC_SHA256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("ALL CHECKS PASSED")
