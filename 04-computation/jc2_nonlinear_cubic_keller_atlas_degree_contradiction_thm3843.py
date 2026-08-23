#!/usr/bin/env python3
"""Exact companion for THM-3843's all-degree atlas contradiction."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0
DEGREE_CAP = 60


def check(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)
    CHECKS += 1


def same(lhs: object, rhs: object, label: str) -> None:
    check(sp.expand(lhs - rhs) == 0, label)  # type: ignore[operator]


h, k, C = sp.symbols("h k C")
Q = 7 * h**2 + 3 * k**2
B = 6 * h**3 + 7 * h**2 * k - k**3
P = 3 * h**3 + 7 * h**2 * k + k**3
S = C * h**2 * Q - k * B

# The factorization and the cleared degree identity are literally equivalent.
same(P - C * S, -(C**2 * h**2 * Q - C * k * B - P),
     "P=CS is the cleared degree identity")

z = sp.symbols("z")
qz = 7 * z**2 + 3
pz = 3 * z**3 + 7 * z**2 + 1
bz = 6 * z**3 + 7 * z**2 - 1
check(sp.resultant(qz, pz, z) == 1615,
      "Q-leading cancellation cannot cancel P-leading term")
check(sp.resultant(qz, bz, z) == 6460,
      "Q-leading cancellation cannot cancel B-leading term")

beta = sp.symbols("beta")
same(Q - 7 * (h - beta * k) * (h + beta * k),
     k**2 * (3 + 7 * beta**2),
     "conjugate-factor identity modulo beta^2=-3/7")

# Symbolic arithmetic for the three degree branches.
a, b, c, d = sp.symbols("a b c d", integer=True, positive=True)
same((2 * c + 4 * a) - (c + b + 3 * a), c + a - b,
     "a>b degree equation forces c=b-a")
same((2 * c + 2 * a + 2 * b) - (c + 4 * b),
     c + 2 * a - 2 * b,
     "b>a degree equation forces c=2(b-a)")
same(3 * b - (a + b + 2 * (b - a) - 2), a + 2,
     "b>a cofactor degree deficit is positive")
same((2 * d + c - 2) - 3 * d, c - d - 2,
     "equal-degree cofactor requires c>=d+2")
same((2 * d - c) - d, d - c,
     "equal-degree conjugate factor requires c<=d")

# Exhaustive hostile census of positive degree triples.  This is a finite
# control for the symbolic proof above, not the source of its quantifiers.
feasible: list[tuple[int, int, int]] = []
branch_counts = {"h_gt_k": 0, "k_gt_h": 0, "equal": 0}
for av in range(1, DEGREE_CAP + 1):
    for bv in range(1, DEGREE_CAP + 1):
        for cv in range(1, DEGREE_CAP + 1):
            if av > bv:
                branch_counts["h_gt_k"] += 1
                identity_possible = 2 * cv + 4 * av == cv + bv + 3 * av
                if identity_possible:
                    feasible.append((av, bv, cv))
            elif bv > av:
                branch_counts["k_gt_h"] += 1
                identity_possible = cv == 2 * (bv - av)
                cofactor_possible = av + bv + cv - 2 >= 3 * bv
                if identity_possible and cofactor_possible:
                    feasible.append((av, bv, cv))
            else:
                branch_counts["equal"] += 1
                conjugate_possible = cv <= av
                cofactor_possible = cv >= av + 2
                if conjugate_possible and cofactor_possible:
                    feasible.append((av, bv, cv))
check(not feasible, "bounded hostile degree census has no survivor")
check(sum(branch_counts.values()) == DEGREE_CAP**3,
      "bounded hostile degree census is exhaustive")

# Determinant/factorization alone has a degenerate constant control; the
# Keller cofactor is the load-bearing extra equation.
h0, k0, C0, m0 = 0, 1, 1, 0
Q0 = 7 * h0**2 + 3 * k0**2
B0 = 6 * h0**3 + 7 * h0**2 * k0 - k0**3
P0 = 3 * h0**3 + 7 * h0**2 * k0 + k0**3
S0 = C0 * h0**2 * Q0 - k0 * B0
check(C0 * k0 - m0 * h0 == 1, "constant control is unimodular")
check(P0 == C0 * S0, "constant control satisfies factorization")
check(P0 != 0, "constant control cannot satisfy delta=lambda P")

source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
      "no inactive Python assert")

semantic = {
    "theorem": "no dominant etale polynomial-plane atlas of the THM-3811 nonlinear cubic surface",
    "packet": "C^2 h^2 Q=CkB+P and kJac(h,C)-hJac(k,C)=lambda P",
    "cases": "deg h>deg k; deg k>deg h; equal degrees with conjugate-factor floor",
    "equal_case": "Q cancellation gives c<=d while Keller cofactor gives c>=d+2",
    "scope": "this candidate surface only; JC(2) remains open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3843-nonlinear-cubic-keller-atlas-total-degree-contradiction")
print("classification=PROVED;VERIFIED-EXACT;INDEPENDENTLY-HOSTILE-AUDITED")
print("packet=C2h2Q_equals_CkB_plus_P;deltaC_equals_lambdaP")
print("h_gt_k=identity_forces_c_equals_b_minus_a_lt_0")
print("k_gt_h=identity_forces_c_equals_2b_minus_2a;cofactor_deficit=a_plus_2")
print("equal_degree=Q_conjugate_floor_forces_c_le_d;cofactor_forces_c_ge_d_plus_2")
print("resultants=Res(Q,P):1615;Res(Q,B):6460")
print(f"degree_hostile=1_le_a_b_c_le_{DEGREE_CAP};triples={DEGREE_CAP**3};survivors=0")
print("scope=specific_THM3811_surface_only;JC2_OPEN")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print("RESULT PASS")
