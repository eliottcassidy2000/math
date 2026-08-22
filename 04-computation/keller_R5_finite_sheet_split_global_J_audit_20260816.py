#!/usr/bin/env python3
"""Independent split-prime audit of the fixed R_5 finite-sheet value.

Unlike the recursive L/H companion, this audit reconstructs the frozen
66,146-term polynomial J.  At p=71 the outer inverse cubic above
q=(2,5/6,-7/8) splits into three rational roots.  The audit evaluates
N(J) separately at those three inverse points, forms G=L^43 N(J), and
multiplies the three G-values to obtain R_5(q)=L(q)^271 N(G)(q).

This is intentionally slower but representation-disjoint from the primary
tower computation.  Scope is the same fixed-map finite-sheet gate only.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import runpy
from fractions import Fraction
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT / "04-computation/keller_tropical_norm_face_recurrence_probe_20260816.py"
PARENT_SHA256 = "fe1b03de9061c997a1abba6b88753e589c0a01f9f92762e49fe7ea0504ce9797"
require(hashlib.sha256(PARENT.read_bytes()).hexdigest() == PARENT_SHA256, "parent probe changed")

captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    namespace = runpy.run_path(str(PARENT))
require(captured.getvalue().rstrip().endswith("all exact checks passed"), "parent replay failed")

prime = 71
q = (Fraction(2), Fraction(5, 6), Fraction(-7, 8))
fraction_mod = namespace["fraction_mod"]
invariants = namespace["invariants"]
inv_mod = namespace["inv_mod"]
evaluate_j_norm_at_target = namespace["evaluate_j_norm_at_target"]


def fixed_map_mod(point):
    x, y, z = point
    u = (1 + x * y) % prime
    return (
        (u**3 * z + y**2 * u * (4 + 3 * x * y)) % prime,
        (y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y)) % prime,
        (2 * x - 3 * x**2 * y - x**3 * z) % prime,
    )


target = tuple(fraction_mod(value, prime) for value in q)
outer = invariants(*target, prime)
require(outer["L"] != 0 and outer["S"] != 0, "outer inverse denominator vanished")
outer_roots = [
    w
    for w in range(prime)
    if (outer["L"] * w**3 + outer["T"] * w - 2 * target[2]) % prime == 0
]
require(outer_roots == [10, 23, 38], "outer split-root ledger changed")

rows = []
norm_g = 1
for w in outer_roots:
    y = (
        outer["Y0"]
        + 6 * outer["L"] * w
        - 3 * outer["K"] * outer["L"] * w * w
    ) * inv_mod(2 * outer["S"], prime) % prime
    z = (
        outer["Z0"]
        + 6 * outer["L"] * outer["A1"] * w
        - 9 * outer["L"] * outer["A2"] * w * w
    ) * inv_mod(8 * outer["S"], prime) % prime
    point = (w, y, z)
    require(fixed_map_mod(point) == target, f"w={w}: outer inverse graph changed")
    norm_j, l_value, inner_discriminant = evaluate_j_norm_at_target(
        tuple(Fraction(value) for value in point), prime
    )
    g_value = pow(l_value, 43, prime) * norm_j % prime
    require(norm_j != 0 and l_value != 0 and inner_discriminant != 0 and g_value != 0, f"w={w}: bad inner gate")
    rows.append((w, y, z, norm_j, l_value, inner_discriminant, g_value))
    norm_g = norm_g * g_value % prime

require(
    rows
    == [
        (10, 0, 2, 68, 26, 4, 64),
        (23, 1, 36, 14, 57, 68, 66),
        (38, 18, 60, 40, 20, 44, 19),
    ],
    "split global-J row ledger changed",
)
r5_value = pow(outer["L"], 271, prime) * norm_g % prime
require(r5_value == 43, "split global-J R5 value changed")

semantic = "\n".join(
    [f"p={prime};outerL={outer['L']};R5={r5_value}"]
    + [":".join(map(str, row)) for row in rows]
)
semantic_sha256 = hashlib.sha256(semantic.encode("ascii")).hexdigest()

print("== fixed Keller R5 split global-J audit ==")
print(f"parent sha256={PARENT_SHA256}")
print(f"p={prime}; outer roots={outer_roots}; rows={rows}")
print(f"R5(2,5/6,-7/8) mod {prime}={r5_value}")
print(f"semantic sha256={semantic_sha256}")
print("scope: independent finite-sheet audit only; no renewal, image-prime, degree-243, or general JC claim")
print("all exact checks passed")
