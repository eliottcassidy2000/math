#!/usr/bin/env python3
"""Exact controls for the common-AP step-three 2x4 cubic-pseudoplane peel.

This is a research candidate, not canon.  It verifies the four exhaustive
sign rows at d=3, the endpoint-power and adjacent-collision identities that
force the scalar bucket to retain the arm divisor, and the first d=4 hostile
rows.  No arbitrary d or multiple-collision conclusion is encoded.
"""

from __future__ import annotations

import ast
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def same(left: sp.Expr, right: sp.Expr, message: str) -> None:
    gate(sp.simplify(sp.together(left - right)) == 0, message)


def polynomial_factor(left: sp.Expr, divisor: sp.Expr,
                      quotient: sp.Expr, message: str) -> None:
    """Certify a displayed factor with a denominator-free quotient."""
    same(left, divisor * quotient, message)
    gate(sp.denom(sp.together(quotient)) == 1, f"{message}: quotient not polynomial")


def bracket(u: int, f: sp.Expr, v: int, g: sp.Expr, w: sp.Symbol) -> sp.Expr:
    return sp.expand(u * f * sp.diff(g, w) - v * sp.diff(f, w) * g)


def sign_possible(u: int, v: int) -> bool:
    return u * v > 0 or (u == 0 and v == 0)


def ap_rows(d: int) -> dict[int, tuple[tuple[int, int], ...]]:
    offsets = (0, d, 2 * d, 3 * d)
    buckets: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for i, u in enumerate((0, d)):
        for j, v in enumerate(offsets):
            buckets[u + v].append((i, j))
    result: dict[int, tuple[tuple[int, int], ...]] = {}
    for scalar_q, scalar_edges in buckets.items():
        if len(scalar_edges) != 2:
            continue
        rows: list[tuple[int, int]] = []
        for a in range(-scalar_q - 1, 0):
            b = -scalar_q - 2 - a
            if b >= 0:
                continue
            fw = (a, a + d)
            gw = tuple(b + v for v in offsets)
            if all(
                sign_possible(fw[i], gw[j])
                for q, edges in buckets.items()
                if q != scalar_q and len(edges) == 1
                for i, j in edges
            ):
                rows.append((a, b))
        result[scalar_q] = tuple(rows)
    return result


rows3 = ap_rows(3)
rows4 = ap_rows(4)
gate(
    rows3 == {
        3: ((-2, -3), (-1, -4)),
        6: ((-2, -6), (-1, -7)),
        9: (),
    },
    f"step-three row table changed: {rows3}",
)
gate(
    rows4 == {
        4: ((-3, -3), (-2, -4), (-1, -5)),
        8: ((-3, -7), (-2, -8), (-1, -9)),
        12: ((-3, -11),),
    },
    f"step-four hostile table changed: {rows4}",
)


w, c = sp.symbols("w c", nonzero=True)
Delta = w**3 - c**3
T = w**3
lam, mu = sp.symbols("lam mu", nonzero=True)
nu, kap = sp.symbols("nu kap")


# Scalar q=3, row (a,b)=(-2,-3).
# Bottom commutation gives p=t^2, g0=lam*t^3.  Since the negative profile
# p is divisible by squarefree Delta, polynomiality forces Delta|t.
U = sp.Function("U")(w)
t = Delta * U
p = t**2
q = w * sp.Function("Q1")(T)
g0 = lam * t**3
g1 = sp.Function("G0")(T)
s31 = bracket(-2, p, 0, g1, w) + bracket(1, q, -3, g0, w)
s31_q = U**2 * (
    -2 * sp.diff(g1, w)
    + 3 * lam * q * sp.diff(t, w)
    + 3 * lam * sp.diff(q, w) * t
)
polynomial_factor(s31, Delta**2, s31_q, "q=3 row (-2,-3) Delta^2")


# Scalar q=3, row (a,b)=(-1,-4).
# Bottom commutation gives g0=lam*p^4.  Both weight -1 profiles share Delta.
P = sp.Function("P")(T)
N = sp.Function("N")(T)
p = w**2 * Delta * P
g1 = w**2 * Delta * N
q = w**2 * sp.Function("Q2")(T)
g0 = lam * p**4
s32 = bracket(-1, p, -1, g1, w) + bracket(2, q, -4, g0, w)
s32_q = (
    w**4 * (sp.diff(P, w) * N - P * sp.diff(N, w))
    + 8 * lam * q * w**6 * Delta * P**3 * sp.diff(p, w)
    + 4 * lam * sp.diff(q, w) * w**8 * Delta**2 * P**4
)
polynomial_factor(s32, Delta**2, s32_q, "q=3 row (-1,-4) Delta^2")


# Scalar q=6, row (a,b)=(-2,-6).
# Endpoints give g0=lam*p^3 and g3=mu*q^3.  The adjacent equations have
# g1=3*lam*p^2*q+H with [-2,p;-3,H]=0 and
# g2=3*mu*p*q^2+nu.  First audit the exact reductions.
p = w * Delta * sp.Function("P2")(T)
q = w * sp.Function("Q3")(T)
H = sp.Function("H")(w)
g0 = lam * p**3
g3 = mu * q**3
g1 = 3 * lam * p**2 * q + H
g2 = 3 * mu * p * q**2 + nu
b1 = bracket(-2, p, -3, g1, w) + bracket(1, q, -6, g0, w)
b3 = bracket(-2, p, 3, g3, w) + bracket(1, q, 0, g2, w)
same(b1, bracket(-2, p, -3, H, w), "q=6 lower residual equation")
same(b3, 0, "q=6 upper collision primitive")
s33 = bracket(-2, p, 0, g2, w) + bracket(1, q, -3, g1, w)
expected33 = 6 * (lam - mu) * p * sp.diff(p * q**2, w) + bracket(1, q, -3, H, w)
same(s33, expected33, "q=6 scalar decomposition")

# If H=0, the scalar retains p and hence Delta.
s33_zero = sp.expand(s33.subs(H, 0))
s33_zero_q = (
    6 * (lam - mu) * w * sp.Function("P2")(T)
    * sp.diff(p * q**2, w)
)
polynomial_factor(s33_zero, Delta, s33_zero_q, "q=6 zero-residual Delta")
same(s33_zero.subs(nu, 0), (Delta * s33_zero_q).subs(nu, 0),
     "q=6 zero residual and zero integration constant")

# If H is nonzero, commutation forces p=t^2,H=kap*t^3.  Since Delta|p and
# Delta is squarefree, Delta|t; this branch actually retains Delta^2.
U2 = sp.Function("U2")(w)
t2 = Delta * U2
p2 = t2**2
H2 = kap * t2**3
g1_2 = 3 * lam * p2**2 * q + H2
g2_2 = 3 * mu * p2 * q**2 + nu
s33_nonzero = bracket(-2, p2, 0, g2_2, w) + bracket(1, q, -3, g1_2, w)
s33_nonzero_q = (
    6 * (lam - mu) * U2**2 * sp.diff(p2 * q**2, w)
    + 3 * kap * q * U2**2 * sp.diff(t2, w)
    + 3 * kap * sp.diff(q, w) * Delta * U2**3
)
polynomial_factor(
    s33_nonzero, Delta**2, s33_nonzero_q,
    "q=6 nonzero-residual Delta^2",
)


# Scalar q=6, row (a,b)=(-1,-7).
# Bottom commutation gives g0=lam*p^7.  Solving the adjacent collision gives
# g1=p^4(7*lam*p^2*q+nu); the other scalar profile has weight -1.
p = w**2 * Delta * sp.Function("P3")(T)
q = w**2 * sp.Function("Q4")(T)
g0 = lam * p**7
g1 = p**4 * (7 * lam * p**2 * q + nu)
g2 = w**2 * Delta * sp.Function("N3")(T)
b1 = bracket(-1, p, -4, g1, w) + bracket(2, q, -7, g0, w)
same(b1, 0, "q=6 row (-1,-7) adjacent primitive")
s34 = bracket(-1, p, -1, g2, w) + bracket(2, q, -4, g1, w)
M = 7 * lam * w**4 * Delta**2 * sp.Function("P3")(T)**2 * q + nu
P3 = sp.Function("P3")(T)
N3 = sp.Function("N3")(T)
J = w**8 * P3**4 * M
s34_q = (
    w**4 * (sp.diff(P3, w) * N3 - P3 * sp.diff(N3, w))
    + 2 * q * (
        4 * Delta * sp.diff(Delta, w) * J
        + Delta**2 * sp.diff(J, w)
    )
    + 4 * sp.diff(q, w) * Delta**2 * J
)
polynomial_factor(s34, Delta**2, s34_q, "q=6 row (-1,-7) Delta^2")
same(s34.subs(nu, 0), (Delta**2 * s34_q).subs(nu, 0),
     "q=6 row (-1,-7) zero integration constant")


# Squarefreeness and endpoint-power controls.
same(sp.gcd(Delta, sp.diff(Delta, w)), 1, "Delta not squarefree")
z = sp.Function("z")(w)
same(bracket(-2, z**2, -3, z**3, w), 0, "(-2,-3) endpoint power")
same(bracket(-1, z, -4, z**4, w), 0, "(-1,-4) endpoint power")
same(bracket(-2, z, -6, z**3, w), 0, "(-2,-6) endpoint power")
same(bracket(-1, z, -7, z**7, w), 0, "(-1,-7) endpoint power")

semantic = {
    "candidate": "all common-AP exact 2x4 and 4x2 cells with d=3 are empty",
    "rows": rows3,
    "mechanism": "endpoint powers plus adjacent collision ODE retain Delta in scalar bucket",
    "open": "d>=4 AP, two-disjoint-pair, other multiple-collision cells",
    "not_claimed": "arbitrary 2x4, arbitrary Darboux pair, JC2",
}
semantic_hash = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()
source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))), "bare assert")

print("CANDIDATE=cubic-pseudoplane-common-AP-step-three-peel")
print(f"D3_ROWS={rows3};ALL_CLOSED=Delta")
print("D3_q3=(-2,-3):Delta^2;(-1,-4):Delta^2")
print("D3_q6=(-2,-6):Delta_or_Delta^2;(-1,-7):Delta^2")
print(f"OPEN_CONTROL_D4={rows4}")
print("NEXT_OPEN=common_AP_d>=4;two_disjoint_pairs;three_chain_plus_isolated")
print(f"SEMANTIC_SHA256={semantic_hash}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
