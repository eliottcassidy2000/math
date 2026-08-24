#!/usr/bin/env python3
"""Clean-room audit of the THM-3905 equality-color third response."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


e, a, L, kappa = sp.symbols("epsilon a L kappa")
u = sp.symbols("u0:4")
v = sp.symbols("v0:4")

# (label, coefficient, T exponent, f exponent, K exponent)
support = (
    ("L4", L**4, 0, 0, 0),
    ("f", 6 * a * L**4, 0, 1, 0),
    ("f2", 12 * a**2 * L**4, 0, 2, 0),
    ("f3", 8 * a**3 * L**4, 0, 3, 0),
    ("KTf", 6 * L**2, 1, 1, 1),
    ("KTf2", 12 * a * L**2, 1, 2, 1),
    ("KTf3", 6 * a**2 * L**2, 1, 3, 1),
    ("T2", -6 * a * L**2, 2, 0, 0),
    ("T2f", -6 * a**2 * L**2, 2, 1, 0),
    ("T2f2", 3 * a**3 * L**2, 2, 2, 0),
    ("KT3", -8, 3, 0, 1),
    ("KT3f", -6 * a, 3, 1, 1),
    ("T4", -3 * a**2, 4, 0, 0),
    ("K2f3", 2 * L**2, 0, 3, 2),
    ("K2f4", 3 * a * L**2, 0, 4, 2),
    ("K2T2f2", -3, 2, 2, 2),
)


def deficit(n: int, p: int, q: int, r: int) -> int:
    return (4 - p - q) * n + 4 - 2 * r


included = []
for label, _coefficient, p, q, r in support:
    rows = tuple(deficit(n, p, q, r) for n in range(1, 8))
    if min(rows) <= 3:
        included.append((label, p, q, r, rows[:4]))
    else:
        require(all(value >= 4 for value in rows), f"excluded cutoff: {label}")

require(
    [row[0] for row in included]
    == ["KTf2", "KTf3", "KT3", "KT3f", "K2f3", "K2f4", "K2T2f2"],
    "deficit <=3 support list",
)


def truncated_series(coefficients: tuple[sp.Symbol, ...], n: int) -> sp.Expr:
    return sum(coefficients[j] * e**j for j in range(min(n, 3) + 1))


checks = 0
for n in range(1, 8):
    U = truncated_series(u, n)
    V = truncated_series(v, n)
    T = e ** (-n) * U
    f = e ** (-n) * V
    K = e ** (-2) * (1 + kappa * e**2)
    P = a * L**2
    r = a * T + K * f
    A = K * T + a * P * f
    S = sp.expand(
        L**4
        + 2 * (3 * A + 3 * P + r**2) * L**2 * f
        + (8 * A + 6 * P + 3 * r**2) * (P * f**2 - T**2)
    )

    # Multiply the proposed C_-C_+ identity by V^2 so no formal division is
    # hidden in the audit.
    proposed_cleared = sp.expand(
        3 * a * L**2 * V**4
        + 2 * L**2 * e**n * (1 + 2 * kappa * e**2) * V**3
        + 6 * e**2 * (kappa * V**2 + a * U * V) * (a * L**2 * V**2 - U**2)
        + (1 if n == 1 else 0)
        * e**3
        * (12 * a * L**2 * U * V**2 - 8 * U**3)
    )
    actual_cleared = sp.expand(e ** (4 * n + 4) * S + 3 * U**2 * V**2)
    difference = sp.expand(actual_cleared - proposed_cleared)
    for order in range(4):
        require(sp.expand(difference).coeff(e, order) == 0, f"n={n}, eps^{order}")
        checks += 1

# Freeze the coefficient-epsilon^3 boundary in normalized form.
u0, u1, _u2, _u3 = u
v0, v1, v2, v3 = v
E0 = a * L**2 * v0**2 - u0**2
E1 = 2 * a * L**2 * v0 * v1 - 2 * u0 * u1
ratio1 = (u1 * v0 - u0 * v1) / v0**2
d1, d2, d3 = sp.symbols("delta1 delta2 delta3")
R3 = sp.expand(
    6 * a * L**2 * (v0 * v3 + v1 * v2)
    + 2 * L**2 * (d1 * v2 + d2 * v1 + d3 * v0)
    + 4 * kappa * L**2 * d1 * v0
    + 6 * ((kappa + a * u0 / v0) * E1 + a * ratio1 * E0)
    + d1 * (12 * a * L**2 * u0 - 8 * u0**3 / v0**2)
)
require(sp.denom(sp.cancel(v0**3 * R3)) == 1, "epsilon^3 denominator clearing")

semantic = {
    "support": [(label, p, q, r, rows) for label, p, q, r, rows in included],
    "law": (
        "3aL2V2+2L2 eps^n(1+2kappa eps^2)V+"
        "6eps^2(kappa+aU/V)(aL2V2-U2)+"
        "[n=1]eps^3(12aL2U-8U3/V2) mod eps^4"
    ),
    "boundaries": {
        "n1": "marked eps source plus K2 curvature and KTf2/KT3 source",
        "n2": "marked eps^2 source",
        "n3": "marked eps^3 source",
        "n_ge_4": "no marked source through eps^3",
    },
    "scope": "necessary equality-color third jet only; no square existence or JC2 closure",
}
blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("THM3905_EQUALITY_COLOR_THIRD_RESPONSE_INDEPENDENT_AUDIT")
print("status=PASS;candidate_exact_mod_epsilon4;necessity_only")
for label, p, q, r, rows in included:
    print(f"support={label};tuple=({p},{q},{r});deficits_n1_to_n4={rows}")
print(f"coefficient_checks={checks}")
print(f"epsilon3_cleared={sp.factor(v0**3 * R3)}")
print(f"semantic_sha256={hashlib.sha256(blob).hexdigest()}")
print("ALL CHECKS PASSED")
