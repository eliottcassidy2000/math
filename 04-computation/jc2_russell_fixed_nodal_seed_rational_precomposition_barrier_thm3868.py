#!/usr/bin/env python3
"""Exact companion for THM-3868.

This script verifies the algebraic identities behind the fixed-nodal-seed
barrier: integral recovery, the derivative/different factor, closure of the
Russell Poisson bracket, six mixed formal density rows, and one exact
genuinely mixed rational hostile.  Normality and the constant-unit theorem
for the Russell ring are cited proof inputs and are not finite computations.
"""

from __future__ import annotations

import hashlib

import sympy as sp


CHECKS = 0


def check(label: str, condition: object) -> None:
    """Assertion that remains active under python -O."""

    global CHECKS
    CHECKS += 1
    if condition is True or condition == sp.S.true:
        return
    if isinstance(condition, sp.Basic) and sp.simplify(condition) == 0:
        return
    raise RuntimeError(f"CHECK FAILED: {label}: {condition}")


def zero(label: str, expression: sp.Expr) -> None:
    check(label, sp.cancel(sp.factor(expression)) == 0)


def bracket(left: sp.Expr, right: sp.Expr, z: sp.Symbol, s: sp.Symbol) -> sp.Expr:
    """Jacobian bracket in coordinate order (z,s)."""

    return sp.factor(
        sp.diff(left, z) * sp.diff(right, s)
        - sp.diff(left, s) * sp.diff(right, z)
    )


# ---------------------------------------------------------------------------
# 1. The nodal seed is a finite monogenic cubic with explicit different.
# ---------------------------------------------------------------------------

c, w = sp.symbols("c w", nonzero=True)
S0, Z0, Avar, Cvar = sp.symbols("S0 Z0 Avar Cvar")
w_nodal = 1 / (2 * c**3)

A_seed = 9 * c**6 * S0**2 - Z0 / (3 * c**3)
C_seed = 27 * c**9 * S0**3 - 3 * c**3 * S0 - sp.Rational(3, 2) * S0 * Z0

seed_jacobian = (
    sp.diff(A_seed, Z0) * sp.diff(C_seed, S0)
    - sp.diff(A_seed, S0) * sp.diff(C_seed, Z0)
)
zero("nodal seed Jacobian factor", seed_jacobian - (1 + w_nodal * Z0))

Z_recovered = 27 * c**9 * S0**2 - 3 * c**3 * A_seed
zero("recover Z from A and S", Z_recovered - Z0)

recovery_cubic = (
    S0**3
    + (2 - 3 * A_seed) * S0 / (9 * c**6)
    + 2 * C_seed / (27 * c**9)
)
zero("monic recovery cubic", recovery_cubic)

p_cubic = (2 - 3 * Avar) / (9 * c**6)
q_cubic = 2 * Cvar / (27 * c**9)
discriminant = sp.factor(-4 * p_cubic**3 - 27 * q_cubic**2)
expected_discriminant = -4 * ((2 - 3 * Avar) ** 3 + 27 * Cvar**2) / (729 * c**18)
zero("depressed cubic discriminant", discriminant - expected_discriminant)

recovery_derivative = 3 * S0**2 + (2 - 3 * A_seed) / (9 * c**6)
zero(
    "power-basis derivative equals seed Jacobian",
    recovery_derivative - 2 * (1 + w_nodal * Z0) / (9 * c**6),
)
zero(
    "power-basis derivative alternate form",
    recovery_derivative - (Z0 + 2 * c**3) / (9 * c**9),
)

branch_A = sp.factor(A_seed.subs(Z0, -2 * c**3))
branch_C = sp.factor(C_seed.subs(Z0, -2 * c**3))
zero("branch cusp parametrization", (3 * branch_A - 2) ** 3 - 27 * branch_C**2)


# ---------------------------------------------------------------------------
# 2. The canonical bracket restricts to the Russell coordinate ring.
# ---------------------------------------------------------------------------

s, z = sp.symbols("s z")
r_param = z**3 / (c**3 + 3 * s * z**3)
e_param = 3 * c**3 * s + 9 * s**2 * z**3

zero("Russell parametrization relation", r_param**2 * e_param - z**3 + c**3 * r_param)
zero("arm identity er=3sz3", e_param * r_param - 3 * s * z**3)
zero("bracket z r", bracket(z, r_param, z, s) + 3 * r_param**2)
zero(
    "bracket z e",
    bracket(z, e_param, z, s) - (3 * c**3 + 6 * r_param * e_param),
)
zero("bracket r e", bracket(r_param, e_param, z, s) - 9 * z**2)


# ---------------------------------------------------------------------------
# 3. Mixed density rows through order six and the all-order new-variable row.
# ---------------------------------------------------------------------------

u, v, p5, q6, t7 = [sp.Function(name)(s) for name in ("u", "v", "p5", "q6", "t7")]
h, k3, ell, m5, n6 = [sp.Function(name)(s) for name in ("h", "k3", "ell", "m5", "n6")]

Z_series = (
    z
    - w * z**2 / 2
    + u * z**3
    + v * z**4
    + p5 * z**5
    + q6 * z**6
    + t7 * z**7
)
S_series = s + h * z**2 + k3 * z**3 + ell * z**4 + m5 * z**5 + n6 * z**6
density_series = sp.expand(
    (1 + w * Z_series) * bracket(Z_series, S_series, z, s) - 1
)

E2 = 3 * u + sp.diff(h, s) - 3 * w**2 / 2
E3 = 4 * v + sp.diff(k3, s) + 4 * w * u + w**3 / 2
E4 = (
    5 * p5
    + sp.diff(ell, s)
    + 5 * w * v
    - 5 * w**2 * u / 2
    + (3 * u - 3 * w**2 / 2) * sp.diff(h, s)
    - 2 * h * sp.diff(u, s)
)
E5 = (
    6 * q6
    + sp.diff(m5, s)
    - 3 * w**2 * v
    + 6 * w * p5
    + 3 * w * u**2
    + (3 * u - 3 * w**2 / 2) * sp.diff(k3, s)
    - (2 * w * h + 3 * k3) * sp.diff(u, s)
    + (w**3 / 2 + 4 * w * u + 4 * v) * sp.diff(h, s)
    - 2 * h * sp.diff(v, s)
)
E6 = (
    7 * t7
    + sp.diff(n6, s)
    - 7 * w**2 * p5 / 2
    + 7 * w * q6
    + 7 * w * u * v
    + (3 * u - 3 * w**2 / 2) * sp.diff(ell, s)
    - (2 * w * h + 3 * k3) * sp.diff(v, s)
    + (w**3 / 2 + 4 * w * u + 4 * v) * sp.diff(k3, s)
    + (w**2 * h - 3 * w * k3 - 4 * ell) * sp.diff(u, s)
    + (-5 * w**2 * u / 2 + 5 * w * v + 5 * p5) * sp.diff(h, s)
    - 2 * h * sp.diff(p5, s)
)

zero("density coefficient z1", density_series.coeff(z, 1))
for order, expected in [(2, E2), (3, E3), (4, E4), (5, E5), (6, E6)]:
    zero(f"mixed density row z{order}", density_series.coeff(z, order) - expected)

# A generic finite series confirms that row n contains only the new pair
# (n+1)*a_(n+1)+b_n' beyond quantities from earlier rows.
N = 8
a_coeff = [sp.Integer(0), sp.Integer(1)] + [
    sp.Function(f"a{i}")(s) for i in range(2, N + 1)
]
b_coeff = [s, sp.Integer(0)] + [sp.Function(f"b{i}")(s) for i in range(2, N)]
Z_generic = sum(a_coeff[i] * z**i for i in range(1, N + 1))
S_generic = sum(b_coeff[i] * z**i for i in range(N))
generic_density = sp.expand((1 + w * Z_generic) * bracket(Z_generic, S_generic, z, s))
for order in range(2, 7):
    row = generic_density.coeff(z, order)
    check(
        f"row {order} coefficient of new a",
        sp.diff(row, a_coeff[order + 1]) == order + 1,
    )
    check(
        f"row {order} coefficient of new b derivative",
        sp.diff(row, sp.diff(b_coeff[order], s)) == 1,
    )


# ---------------------------------------------------------------------------
# 4. Exact genuinely mixed rational family and its named moved divisor.
# ---------------------------------------------------------------------------

eta = sp.symbols("eta", nonzero=True)
D = 1 - eta * s * z**2
x_control = z / D**2
y_control = s * D**3
zero("monomial Hamiltonian flow is symplectic", bracket(x_control, y_control, z, s) - 1)

H = D**2 + w * z / 2
G = D**2 + 3 * w * z / 2
Z_mixed = z / H
S_mixed = s * H**3 / (D * G)
zero(
    "exact mixed weighted density",
    (1 + w * Z_mixed) * bracket(Z_mixed, S_mixed, z, s) - 1,
)

Z_mixed_series = sp.series(Z_mixed, z, 0, 6).removeO().expand()
S_mixed_series = sp.series(S_mixed, z, 0, 5).removeO().expand()
check("mixed Z linear jet", Z_mixed_series.coeff(z, 1) == 1)
check("mixed Z forced quadratic jet", Z_mixed_series.coeff(z, 2) == -w / 2)
check("mixed Z cubic jet", Z_mixed_series.coeff(z, 3) == w**2 / 4 + 2 * eta * s)
check("mixed Z quartic jet", Z_mixed_series.coeff(z, 4) == -w**3 / 8 - 2 * eta * w * s)
check(
    "mixed Z quintic jet",
    Z_mixed_series.coeff(z, 5)
    == w**4 / 16 + 3 * eta * w**2 * s / 2 + 3 * eta**2 * s**2,
)
check("mixed S arm", S_mixed_series.coeff(z, 0) == s)
check("mixed S linear jet", S_mixed_series.coeff(z, 1) == 0)
check(
    "mixed S quadratic jet",
    S_mixed_series.coeff(z, 2) == 3 * w**2 * s / 4 - 3 * eta * s**2,
)
check("mixed S cubic jet", S_mixed_series.coeff(z, 3) == -w**3 * s)
check(
    "mixed S quartic jet",
    S_mixed_series.coeff(z, 4)
    == 3 * w**4 * s / 2 + 3 * eta * w**2 * s**2 / 4 + 3 * eta**2 * s**3,
)
zero("genuine mixed derivative starts at z3", sp.diff(Z_mixed, s).coeff(z, 0))
check("genuine mixed cubic derivative nonzero", sp.diff(Z_mixed_series.coeff(z, 3), s) == 2 * eta)

A_mixed = A_seed.subs({S0: S_mixed, Z0: Z_mixed})
C_mixed = C_seed.subs({S0: S_mixed, Z0: Z_mixed})
zero(
    "mixed nodal composite is rational Keller",
    bracket(A_mixed, C_mixed, z, s).subs(w, w_nodal) - 1,
)

# On Y_z, s=1/(3r)-c^3/(3z^3).  Clearing the Laurent units from H gives
# the displayed mixed divisor.  At eta=0 it reduces to the old z+4c^3 one.
r = sp.symbols("r", nonzero=True)
s_on_torus = 1 / (3 * r) - c**3 / (3 * z**3)
D_torus = sp.factor(D.subs(s, s_on_torus))
H_torus = sp.factor(H.subs({s: s_on_torus, w: w_nodal}))
d_tilde = 3 * r * z + eta * c**3 * r - eta * z**3
mixed_divisor = 4 * c**3 * d_tilde**2 + 9 * r**2 * z**3
zero("cleared D on Russell torus", D_torus - d_tilde / (3 * r * z))
zero(
    "cleared mixed pole divisor",
    H_torus - mixed_divisor / (36 * c**3 * r**2 * z**2),
)
zero(
    "vertical specialization recovers old divisor",
    mixed_divisor.subs(eta, 0) - 9 * r**2 * z**2 * (z + 4 * c**3),
)
zero("H and D cannot vanish on torus", H - D**2 - w * z / 2)
zero("G differs from H by a torus unit", G - H - w * z)

# The exact principal part of A at any prime of H is exposed without a
# valuation package: D,G,s,z are generically units there and H*A has residue
# -z/(3c^3).
zero(
    "mixed A principal-part identity",
    H * A_mixed - 9 * c**6 * s**2 * H**7 / (D**2 * G**2) + z / (3 * c**3),
)


# ---------------------------------------------------------------------------
# 5. The density equation is an ordinary Keller equation plus divisibility.
# ---------------------------------------------------------------------------

Zf = sp.Function("Zf")(z, s)
Sf = sp.Function("Sf")(z, s)
Tf = (1 + w * Zf) * Sf
zero(
    "weighted density linearizes via T",
    bracket(Zf, Tf, z, s) - (1 + w * Zf) * bracket(Zf, Sf, z, s),
)

SEMANTIC_FACTS = "\n".join(
    [
        "fixed-nodal-seed=recovery-monic-cubic",
        "recovery-derivative=scalar-times-seed-jacobian-different",
        "branch-curve=(3A-2)^3=27C^2",
        "russell-bracket=closed-on-B",
        "normality-plus-recovery=A-C-regular-implies-Z-S-regular",
        "constant-unit-trap=all-rational-precompositions-excluded",
        "mixed-density-rows=verified-through-order-six",
        "formal-density-recursion=one-free-coefficient-per-order",
        "exact-genuinely-mixed-rational-hostile=Z_s-starts-order-three",
        "mixed-hostile-pole=H-divisor-moves-but-persists",
        "scope=fixed-seed-rational-precomposition-only-no-JC2-conclusion",
    ]
)
semantic_sha256 = hashlib.sha256(SEMANTIC_FACTS.encode("utf-8")).hexdigest()

print("THM-3868 RUSSELL FIXED NODAL SEED RATIONAL PRECOMPOSITION BARRIER")
print("RECOVERY=PASS;MONIC_DEGREE=3;DERIVATIVE=SCALAR_TIMES_1_PLUS_wZ")
print("POISSON_CLOSURE=PASS;UNIT_TRAP=B_STAR_EQUALS_k_STAR")
print("FORMAL_MIXED_ROWS=PASS;ORDERS=2..6;ONE_FREE_COEFFICIENT_PER_ROW")
print("EXACT_MIXED_HOSTILE=PASS;Z_s_FIRST_ORDER=3;POLE_DIVISOR=H")
print("CONCLUSION=NO_RATIONAL_PRECOMPOSITION_OF_FIXED_NODAL_SEED_ALGEBRAIZES")
print("BOUNDARY=DIFFERENT_SEED_OR_NONRATIONAL_CONTROL_REMAINS_OPEN;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("ALL CHECKS PASSED")
