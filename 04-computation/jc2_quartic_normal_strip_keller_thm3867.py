#!/usr/bin/env python3
"""Deterministic exact companion for provisional THM-3867."""

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


s, z, x = sp.symbols("s z x")

# ---------------------------------------------------------------------------
# 1. The eight universal quartic coefficient buckets.
# ---------------------------------------------------------------------------

a, alpha, u, p, r = [
    sp.Function(name)(s) for name in ("a", "alpha", "u", "p", "r")
]
b, beta, v, q, t = [
    sp.Function(name)(s) for name in ("b", "beta", "v", "q", "t")
]
A = a + alpha * z + u * z**2 + p * z**3 + r * z**4
C = b + beta * z + v * z**2 + q * z**3 + t * z**4
Jac = sp.expand(sp.diff(A, z) * sp.diff(C, s) - sp.diff(A, s) * sp.diff(C, z))

buckets = [
    alpha * sp.diff(b, s) - sp.diff(a, s) * beta,
    alpha * sp.diff(beta, s) - sp.diff(alpha, s) * beta
    + 2 * u * sp.diff(b, s) - 2 * sp.diff(a, s) * v,
    alpha * sp.diff(v, s) - 2 * sp.diff(alpha, s) * v
    + 2 * u * sp.diff(beta, s) - sp.diff(u, s) * beta
    + 3 * p * sp.diff(b, s) - 3 * sp.diff(a, s) * q,
    alpha * sp.diff(q, s) - 3 * sp.diff(alpha, s) * q
    + 2 * u * sp.diff(v, s) - 2 * sp.diff(u, s) * v
    + 3 * p * sp.diff(beta, s) - sp.diff(p, s) * beta
    + 4 * r * sp.diff(b, s) - 4 * sp.diff(a, s) * t,
    alpha * sp.diff(t, s) - 4 * sp.diff(alpha, s) * t
    + 2 * u * sp.diff(q, s) - 3 * sp.diff(u, s) * q
    + 3 * p * sp.diff(v, s) - 2 * sp.diff(p, s) * v
    + 4 * r * sp.diff(beta, s) - sp.diff(r, s) * beta,
    2 * u * sp.diff(t, s) - 4 * sp.diff(u, s) * t
    + 3 * p * sp.diff(q, s) - 3 * sp.diff(p, s) * q
    + 4 * r * sp.diff(v, s) - 2 * sp.diff(r, s) * v,
    3 * p * sp.diff(t, s) - 4 * sp.diff(p, s) * t
    + 4 * r * sp.diff(q, s) - 3 * sp.diff(r, s) * q,
    4 * (r * sp.diff(t, s) - sp.diff(r, s) * t),
]

for degree, bucket in enumerate(buckets):
    zero(f"generic Jacobian bucket z^{degree}", Jac.coeff(z, degree) - bucket)
check("no coefficient above z7", sp.Poly(Jac, z).degree() <= 7)

aa, bb, cc, dd = sp.symbols("aa bb cc dd")
r_new = aa * r + bb * t
t_new = cc * r + dd * t
top_new = sp.factor(r_new * sp.diff(t_new, s) - sp.diff(r_new, s) * t_new)
zero(
    "top Wronskian target covariance",
    top_new - (aa * dd - bb * cc) * (r * sp.diff(t, s) - sp.diff(r, s) * t),
)

# ---------------------------------------------------------------------------
# 2. The removable (4,1) and (4,2) rows.
# ---------------------------------------------------------------------------

rho, lam, beta0 = sp.symbols("rho lam beta0", nonzero=True)
c3, c2, c1, c0 = sp.symbols("c3 c2 c1 c0")
bfun = sp.Function("bfun")(s)
C41 = bfun + beta0 * z
A41 = (
    rho * C41**4 + c3 * C41**3 + c2 * C41**2 + c1 * C41
    - lam * s / beta0 + c0
)
J41 = sp.expand(
    sp.diff(A41, z) * sp.diff(C41, s) - sp.diff(A41, s) * sp.diff(C41, z)
)
zero("(4,1) triangular family constant Jacobian", J41 - lam)
zero("(4,1) top coefficient", sp.expand(A41).coeff(z, 4) - rho * beta0**4)

At, Ct = sp.symbols("At Ct")
s_inv = beta0 * (
    rho * Ct**4 + c3 * Ct**3 + c2 * Ct**2 + c1 * Ct + c0 - At
) / lam
z_inv = (Ct - bfun.subs(s, s_inv)) / beta0
zero("(4,1) inverse recovers s", s_inv.subs({At: A41, Ct: C41}) - s)
zero("(4,1) inverse recovers z", z_inv.subs({At: A41, Ct: C41}) - z)

# Target shears preserve the bracket and erase the forced top coefficient.
rho_shear = sp.symbols("rho_shear")
generic_shear_2 = sp.expand(A - rho_shear * C**2)
generic_shear_4 = sp.expand(A - rho_shear * C**4)
J_shear_2 = sp.expand(
    sp.diff(generic_shear_2, z) * sp.diff(C, s)
    - sp.diff(generic_shear_2, s) * sp.diff(C, z)
)
J_shear_4 = sp.expand(
    sp.diff(generic_shear_4, z) * sp.diff(C, s)
    - sp.diff(generic_shear_4, s) * sp.diff(C, z)
)
zero("quadratic target shear preserves Jacobian", J_shear_2 - Jac)
zero("quartic target shear preserves Jacobian", J_shear_4 - Jac)

r2, v2 = [sp.Function(name)(s) for name in ("r2", "v2")]
zero(
    "(4,2) top ODE under r=rho*v^2",
    (4 * r2 * sp.diff(v2, s) - 2 * sp.diff(r2, s) * v2).subs(
        {r2: rho * v2**2, sp.diff(r2, s): sp.diff(rho * v2**2, s)}
    ),
)
r1, beta1 = [sp.Function(name)(s) for name in ("r1", "beta1")]
zero(
    "(4,1) top ODE under r=rho*beta^4",
    (4 * r1 * sp.diff(beta1, s) - sp.diff(r1, s) * beta1).subs(
        {r1: rho * beta1**4, sp.diff(r1, s): sp.diff(rho * beta1**4, s)}
    ),
)

# ---------------------------------------------------------------------------
# 3. Kummer scaling and the complete depressed (4,3) packet.
# ---------------------------------------------------------------------------

R, Q = sp.symbols("R Q", nonzero=True)
h = sp.Function("h")(s)
r43 = R * h**4
q43 = Q * h**3
zero("(4,3) Kummer top ODE", 4 * r43 * sp.diff(q43, s) - 3 * sp.diff(r43, s) * q43)

# Verify the moving-scale bracket cancellation directly on generic profiles.
y = sp.symbols("y")
f0, f1, f2, f3 = [sp.Function(f"f{i}")(s) for i in range(4)]
g0, g1, g2 = [sp.Function(f"g{i}")(s) for i in range(3)]
Fbar = R * y**4 + f3 * y**3 + f2 * y**2 + f1 * y + f0
Gbar = Q * y**3 + g2 * y**2 + g1 * y + g0
Fzs = Fbar.subs(y, h * z)
Gzs = Gbar.subs(y, h * z)
Jzs = sp.expand(
    sp.diff(Fzs, z) * sp.diff(Gzs, s) - sp.diff(Fzs, s) * sp.diff(Gzs, z)
)
Jys = sp.expand(sp.diff(Fbar, y) * sp.diff(Gbar, s) - sp.diff(Fbar, s) * sp.diff(Gbar, y))
zero("moving scale gives J_zs=h*J_ys", Jzs - h * Jys.subs(y, h * z))

D, E, Fconst, Gconst = sp.symbols("D E Fconst Gconst")
Bfun = sp.Function("Bfun")(s)
bdep = sp.Function("bdep")(s)
K = 4 * R / (3 * Q)
H0 = D / Q
I = 2 * R / (9 * Q**2)
Jc = 2 * E / (3 * Q)
Udep = K * Bfun + E
Ldep = K * bdep + H0 * Bfun + Fconst
adep = H0 * bdep + I * Bfun**2 + Jc * Bfun + Gconst
Adep = R * x**4 + D * x**3 + Udep * x**2 + Ldep * x + adep
Cdep = Q * x**3 + Bfun * x + bdep
Jdep = sp.expand(
    sp.diff(Adep, x) * sp.diff(Cdep, s)
    - sp.diff(Adep, s) * sp.diff(Cdep, x)
)

for degree in (4, 3, 2):
    zero(f"depressed integrated bucket x^{degree}", Jdep.coeff(x, degree))
W = K * Bfun + 2 * E
x1_expected = W * sp.diff(bdep, s) + (K * bdep + Fconst) * sp.diff(Bfun, s)
x0_expected = (
    (K * bdep + Fconst) * sp.diff(bdep, s)
    - (2 * I * Bfun**2 + Jc * Bfun) * sp.diff(Bfun, s)
)
zero("depressed x1 exact differential", Jdep.coeff(x, 1) - x1_expected)
zero("depressed x0 exact factor", Jdep.coeff(x, 0) - x0_expected)
check("depressed degree at most four", sp.Poly(Jdep, x).degree() <= 4)

# Nonzero-W integration and autonomous one-variable equation.
N = sp.symbols("N")
b_packet = (N - Fconst * Bfun) / W
Sconst = K * N + 2 * E * Fconst
zero("nonzero-W integrated x1", x1_expected.subs(bdep, b_packet).doit())
zero(
    "Kb+F=S/W",
    (K * b_packet + Fconst) - Sconst / W,
)
x0_packet = sp.cancel(x0_expected.subs(bdep, b_packet).doit())
Wprime = sp.diff(W, s)
autonomous = -3 * Q / (16 * R**2) * (
    W**2 - 2 * E * W + 4 * R * Sconst**2 / W**3
) * Wprime
zero("nonzero-W autonomous constant bucket", x0_packet - autonomous)

Phi = W**3 / 3 - E * W**2 - 2 * R * Sconst**2 / W**2
zero(
    "autonomous equation is a rational first integral",
    autonomous + 3 * Q * sp.diff(Phi, s) / (16 * R**2),
)

# Identically-zero W is a distinct packet; no division by W is used.
bzero = sp.Function("bzero")(s)
Bzero = -2 * E / K
Uzero = -E
Lzero = K * bzero + H0 * Bzero + Fconst
azero = H0 * bzero + I * Bzero**2 + Jc * Bzero + Gconst
A_W0 = R * x**4 + D * x**3 + Uzero * x**2 + Lzero * x + azero
C_W0 = Q * x**3 + Bzero * x + bzero
J_W0 = sp.expand(
    sp.diff(A_W0, x) * sp.diff(C_W0, s)
    - sp.diff(A_W0, s) * sp.diff(C_W0, x)
)
for degree in range(1, 5):
    zero(f"W-identically-zero bucket x^{degree}", J_W0.coeff(x, degree))
zero(
    "W-identically-zero constant bucket",
    J_W0.coeff(x, 0) - (K * bzero + Fconst) * sp.diff(bzero, s),
)

# ---------------------------------------------------------------------------
# 4. Exact local residual constants after the best possible arm cancellation.
# ---------------------------------------------------------------------------

Bsym, gsym = sp.symbols("Bsym gsym", nonzero=True)
pole_A_lead = I * Bsym**2 + K * Bsym * gsym**2 + R * gsym**4
pole_residual = sp.factor(pole_A_lead.subs(gsym**2, -Bsym / Q))
zero("W-pole residual after C-arm cancellation", pole_residual + R * Bsym**2 / (9 * Q**2))
zero("W-zero-channel residual coefficient", R - K * Q + R / 3)

# Valuation arithmetic for the three nonconstant-h channels.
m, n = sp.symbols("m n", positive=True, integer=True)
check("W-pole valuation H=3m+1", (-2 * m) + (-m - 1) == -3 * m - 1)
check("W-zero valuation H=2n+1", (-3 * n) + (n - 1) == -2 * n - 1)
check("W-identically-zero valuation H=2n+1", (-n) + (-n - 1) == -2 * n - 1)

# ---------------------------------------------------------------------------
# 5. Two sharp rational hostiles, including W identically zero.
# ---------------------------------------------------------------------------

# First hostile: W=s^-2, h=s^7.  The C arm is canceled by d^2=-3/4;
# every remaining coefficient is polynomial except the forced A arm.
dalg = sp.symbols("dalg")
x_host = s**7 * z + dalg / s
A_host = sp.expand(x_host**4 + x_host**2 / s**2 + sp.Rational(1, 8) / s**4)
C_host = sp.expand(x_host**3 + sp.Rational(3, 4) * x_host / s**2)
J_host = sp.factor(
    sp.diff(A_host, z) * sp.diff(C_host, s)
    - sp.diff(A_host, s) * sp.diff(C_host, z)
)
zero("W-pole hostile constant Jacobian", J_host - sp.Rational(3, 8))

quadratic_relation = sp.Poly(dalg**2 + sp.Rational(3, 4), dalg, domain="QQ[s,z]")


def reduce_dalg(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.together(expression).as_numer_denom()
    reduced = sp.Poly(numerator, dalg, domain="QQ[s,z]").rem(quadratic_relation).as_expr()
    return sp.cancel(reduced / denominator)


A_host_red = sp.expand(reduce_dalg(A_host))
C_host_red = sp.expand(reduce_dalg(C_host))
zero("W-pole hostile C arm cancels", C_host_red.coeff(z, 0))
zero("W-pole hostile A residual", A_host_red.coeff(z, 0) + 1 / (16 * s**4))
for degree in range(1, 5):
    check(
        f"W-pole hostile A non-arm coefficient z^{degree} polynomial",
        not sp.denom(A_host_red.coeff(z, degree)).has(s),
    )
for degree in range(4):
    check(
        f"W-pole hostile C coefficient z^{degree} polynomial",
        not sp.denom(C_host_red.coeff(z, degree)).has(s),
    )

# Second hostile: W is identically zero, h=s^7, and C's arm again cancels.
x_host_W0 = s**7 * z + 1 / s
A_host_W0 = sp.expand(x_host_W0**4 - sp.Rational(4, 3) * x_host_W0 / s**3)
C_host_W0 = sp.expand(x_host_W0**3 - 1 / s**3)
J_host_W0 = sp.factor(
    sp.diff(A_host_W0, z) * sp.diff(C_host_W0, s)
    - sp.diff(A_host_W0, s) * sp.diff(C_host_W0, z)
)
zero("W-identically-zero hostile Jacobian", J_host_W0 + 4)
zero("W-identically-zero hostile C arm cancels", C_host_W0.coeff(z, 0))
zero("W-identically-zero hostile A residual", A_host_W0.coeff(z, 0) + 1 / (3 * s**4))
for degree in range(1, 5):
    check(
        f"W-identically-zero hostile A non-arm z^{degree} polynomial",
        not sp.denom(A_host_W0.coeff(z, degree)).has(s),
    )
for degree in range(4):
    check(
        f"W-identically-zero hostile C z^{degree} polynomial",
        not sp.denom(C_host_W0.coeff(z, degree)).has(s),
    )

# Every hostile Jacobian bucket is checked independently of its factor call.
for label, hostile_jac, constant in (
    ("W-pole", sp.expand(J_host), sp.Rational(3, 8)),
    ("W-identically-zero", sp.expand(J_host_W0), -4),
):
    for degree in range(8):
        zero(
            f"{label} hostile bucket z^{degree}",
            hostile_jac.coeff(z, degree) - (constant if degree == 0 else 0),
        )

SEMANTIC_FACTS = "\n".join(
    [
        "quartic-strip-eight-buckets=exact",
        "top-quartic-row=constant-target-direction",
        "degree-41=quartic-target-shear-to-cubic",
        "degree-42=quadratic-target-shear-to-cubic",
        "degree-43=Kummer-h4-h3-depressed-packet",
        "degree-43-nonzero-W=autonomous-rational-first-integral",
        "degree-43-W-pole=residual-minus-R-B2-over-9Q2",
        "degree-43-W-zero=residual-minus-R-g4-over-3",
        "degree-43-W-identically-zero=separate-two-n-plus-one-channel",
        "constant-h=polynomial-unit-contradiction",
        "rational-hostile-W-pole=Jacobian-3-over-8",
        "rational-hostile-W-identically-zero=Jacobian-minus-4",
        "scope=polynomial-z-depth-at-most-four-only",
    ]
)
semantic_sha256 = hashlib.sha256(SEMANTIC_FACTS.encode("utf-8")).hexdigest()

print("THM-3867 QUARTIC NORMAL-STRIP KELLER CLASSIFICATION")
print("EIGHT_BUCKETS=PASS;TOP_DIRECTION=SL2_CONSTANT")
print("BRANCH_4_1=SHEAR_C4_TO_CUBIC;BRANCH_4_2=SHEAR_C2_TO_CUBIC")
print("BRANCH_4_3=EMPTY;KUMMER=h4_h3;DEPRESSED_PACKET=COMPLETE")
print("NONZERO_W=RATIONAL_FIRST_INTEGRAL;W_IDENTICALLY_ZERO=SPLIT_EXPLICITLY")
print("LOCAL_CHANNELS=W_POLE_3m+1,W_ZERO_2n+1,W_IDZERO_2n+1")
print("SIMULTANEOUS_ARM_GATE=RESIDUALS_-R_B2/(9Q2),_-R_g4/3")
print("RATIONAL_HOSTILES=PASS;JACOBIANS=3/8,-4;ONLY_A_ARM_HAS_POLE")
print("BOUNDARY=TRANSVERSE_DEGREE_FIVE_OR_NONPOLYNOMIAL_SERIES_REMAIN_OPEN")
print(f"SEMANTIC_SHA256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("ALL CHECKS PASSED")
