#!/usr/bin/env python3
"""Exact hostile companion for THM-3787's cubic pseudo-plane support census."""

from __future__ import annotations

import ast
import hashlib
import json
from collections import defaultdict
from itertools import product
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


def bracket(a: int | sp.Expr, f: sp.Expr,
            b: int | sp.Expr, g: sp.Expr, w: sp.Symbol) -> sp.Expr:
    return sp.expand(a * f * sp.diff(g, w) - b * sp.diff(f, w) * g)


# ---------------------------------------------------------------------------
# Exact contribution-bucket census for 2x3 supports.
# ---------------------------------------------------------------------------


def buckets(r: int, s: int, t: int) -> dict[int, list[tuple[int, int]]]:
    answer: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for i, u in enumerate((0, r)):
        for j, v in enumerate((0, s, t)):
            answer[u + v].append((i, j))
    return dict(answer)


def collision_offsets(r: int, s: int, t: int) -> tuple[int, ...]:
    return tuple(sorted(offset for offset, rows in buckets(r, s, t).items()
                        if len(rows) > 1))


collision_counts: dict[str, int] = defaultdict(int)
collision_packet: list[str] = []
for r in range(1, 25):
    for s in range(1, 25):
        for t in range(s + 1, 33):
            actual = collision_offsets(r, s, t)
            predicted = set()
            if r == s:
                predicted.add(r)
            if r == t:
                predicted.add(r)
            if t == r + s:
                predicted.add(t)
            gate(actual == tuple(sorted(predicted)),
                 f"collision grammar r={r},s={s},t={t}")
            gate((len(actual) == 2) == (r == s and t == 2 * r),
                 f"two collisions iff common AP r={r},s={s},t={t}")
            if len(actual) == 0:
                label = "none"
            elif len(actual) == 2:
                label = "common_ap"
            elif r == s:
                label = "lower"
            elif r == t:
                label = "upper"
            else:
                label = "diagonal"
            collision_counts[label] += 1
            collision_packet.append(f"{r},{s},{t}:{actual}")

gate(collision_offsets(3, 1, 7) == (), "no-collision hostile")
gate(collision_offsets(3, 3, 8) == (3,), "lower hostile")
gate(collision_offsets(4, 2, 4) == (4,), "upper hostile")
gate(collision_offsets(3, 2, 5) == (5,), "diagonal hostile")
gate(collision_offsets(3, 3, 6) == (3, 6), "common AP hostile")


# Necessary exact-support commutation rule after removable scalar terms are
# discarded: singleton commuting weights have the same strict sign.  A (0,0)
# endpoint is retained as a seam; a zero/nonzero singleton is removable.
def singleton_commuting_possible(u: int, v: int) -> bool:
    return u * v > 0 or (u == 0 and v == 0)


irregular_counts: dict[str, int] = defaultdict(int)
for r in range(1, 17):
    for s in range(1, 17):
        for t in range(s + 1, 25):
            kind = None
            scalar_offset = None
            if r == s and t != 2 * r:
                kind, scalar_offset = "lower", r
            elif r == t:
                kind, scalar_offset = "upper", r
            elif t == r + s and r != s:
                kind, scalar_offset = "diagonal", t
            if kind is None:
                continue
            survivors = 0
            for a in range(-64, 65):
                b = -2 - scalar_offset - a
                weights = (a, a + r, b, b + s, b + t)
                addresses = {
                    (0, 2): 0,
                    (0, 3): s,
                    (0, 4): t,
                    (1, 2): r,
                    (1, 3): r + s,
                    (1, 4): r + t,
                }
                singleton_edges = [edge for edge, offset in addresses.items()
                                   if offset != scalar_offset]
                if all(singleton_commuting_possible(weights[i], weights[j])
                       for i, j in singleton_edges):
                    survivors += 1
            gate(survivors == 0,
                 f"irregular sign/zero closure {kind} {r},{s},{t}")
            irregular_counts[kind] += 1


# ---------------------------------------------------------------------------
# Common-step 2x3 resonance arithmetic and exact ODE factorizations.
# ---------------------------------------------------------------------------

lower_resonances: list[tuple[int, int]] = []
lower_exception_controls = 0
for d in range(2, 97):
    for A in range(2, d):
        resonance = 3 * A - d - 2
        if resonance == 0:
            lower_resonances.append((d, A))
        # A hypothetical higher nominal P,M degree would require resonance<0
        # and the leading cancellation 2A=d+2.  The two signs are incompatible.
        gate(not (resonance < 0 and 2 * A == d + 2),
             f"lower exceptional branch d={d},A={A}")
        lower_exception_controls += 1

upper_resonances: list[tuple[int, int]] = []
upper_exception_controls = 0
for d in range(4, 97):
    for A in range(3, d):
        resonance = 3 * A - 2 * d - 2
        if resonance == 0:
            upper_resonances.append((d, A))
            gate(A % 2 == 0, f"upper resonance parity d={d},A={A}")
            C0 = (A - 2) // 2
            gate((A, d) == (2 * C0 + 2, 3 * C0 + 2),
                 f"upper resonance parametrization d={d},A={A}")
        # A hypothetical higher nominal first bracket requires resonance>0;
        # leading cancellation would force 2A=d+2 and hence resonance=2-A<0.
        gate(not (resonance > 0 and 2 * A == d + 2),
             f"upper exceptional branch d={d},A={A}")
        upper_exception_controls += 1

w = sp.symbols("w")
p = sp.Function("p")(w)
q = sp.Function("q")(w)
H = sp.Function("H")(w)
M = sp.Function("M")(w)

# Lower A=1 seam: after normalizing endpoint constants, M=p is the unique
# polynomial solution modulo the Delta-divisible homogeneous remainder.
for d in range(2, 18):
    C = d - 1
    nonscalar = bracket(-1, p, C, q, w) + bracket(C, q, -1, p + H, w)
    same(nonscalar, C * q * sp.diff(H, w) + sp.diff(q, w) * H,
         f"lower A=1 ODE d={d}")
    scalar = bracket(-1, p, -1, p, w) + bracket(
        C, q, -(d + 1), p ** (d + 1), w
    )
    same(scalar, (d + 1) * p**d *
         (C * q * sp.diff(p, w) + sp.diff(q, w) * p),
         f"lower A=1 scalar d={d}")

# Lower resonance d=3A-2.  Keep independent endpoint constants l0,m0 so the
# factor checks both a nonzero p branch and the pure homogeneous H branch.
l0, m0 = sp.symbols("l0 m0", nonzero=True)
for A in range(2, 15):
    C = 2 * A - 2
    T = A - 2
    middle = 2 * m0 * p * q + H
    nonscalar = bracket(-A, p, 2 * C, m0 * q**2, w) + bracket(
        C, q, T, middle, w
    )
    same(nonscalar, C * q * sp.diff(H, w) - T * sp.diff(q, w) * H,
         f"lower resonant ODE A={A}")
    scalar = bracket(-A, p, T, middle, w) + bracket(
        C, q, -2 * A, l0 * p**2, w
    )
    scalar_reduced = scalar.subs(
        sp.diff(H, w), sp.Rational(T, C) * sp.diff(q, w) * H / q
    )
    expected = (C * q * sp.diff(p, w) + A * p * sp.diff(q, w)) * (
        2 * (l0 - m0) * p - sp.Rational(T, C) * H / q
    )
    same(scalar_reduced, expected, f"lower scalar factor A={A}")

# Upper resonance A=2C0+2,d=3C0+2.
for C0 in range(1, 15):
    A = 2 * C0 + 2
    middle = 2 * l0 * p * q + H
    nonscalar = bracket(-A, p, -(C0 + 2), middle, w) + bracket(
        C0, q, -2 * A, l0 * p**2, w
    )
    same(nonscalar, -A * p * sp.diff(H, w) +
         (C0 + 2) * sp.diff(p, w) * H,
         f"upper resonant ODE C0={C0}")
    scalar = bracket(-A, p, 2 * C0, m0 * q**2, w) + bracket(
        C0, q, -(C0 + 2), middle, w
    )
    scalar_reduced = scalar.subs(
        sp.diff(H, w), sp.Rational(C0 + 2, A) * sp.diff(p, w) * H / p
    )
    expected = (C0 * q * sp.diff(p, w) + A * p * sp.diff(q, w)) * (
        2 * (l0 - m0) * q + sp.Rational(C0 + 2, A) * H / p
    )
    same(scalar_reduced, expected, f"upper scalar factor C0={C0}")


# ---------------------------------------------------------------------------
# Scalar-centred 3x3 AP: adjacent integrations, elimination, and a=1 seam.
# ---------------------------------------------------------------------------

K = sp.Function("K")(w)
N = sp.Function("N")(w)
alpha = sp.Function("alpha")(w)
beta = sp.Function("beta")(w)
Cconst, Dconst = sp.symbols("Cconst Dconst")

ap_rows = 0
for d in range(2, 15):
    for a in range(1, d):
        low_bucket = bracket(
            -(d + a), K ** (d + a), a - 2, beta, w
        ) + bracket(
            -a, alpha, -(d - a + 2), K ** (d - a + 2), w
        )
        low_factor = (d - a + 2) * K ** (d + 2) * sp.diff(
            alpha / K**a, w
        ) - (d + a) * K ** (d + 2) * sp.diff(
            K ** (a - 2) * beta, w
        )
        same(low_bucket, low_factor, f"3x3 lower integration d={d},a={a}")

        high_bucket = bracket(
            -a, alpha, d + a - 2, N ** (d + a - 2), w
        ) + bracket(
            d - a, N ** (d - a), a - 2, beta, w
        )
        high_factor = (d - a) * N ** (d - 2) * sp.diff(
            beta / N ** (a - 2), w
        ) - (d + a - 2) * N ** (d - 2) * sp.diff(
            alpha * N**a, w
        )
        same(high_bucket, high_factor,
             f"3x3 upper integration d={d},a={a}")
        ap_rows += 1

# Elimination from the two integrated equations is independent of d.
lam, mu = sp.symbols("lam mu", nonzero=True)
for a in range(2, 15):
    beta_expr = mu * alpha * N ** (2 * a - 2) + Dconst * N ** (a - 2)
    eliminated = sp.expand(
        alpha - lam * K ** (2 * a - 2) * beta_expr - Cconst * K**a
    )
    expected = sp.expand(
        alpha * (1 - lam * mu * (K * N) ** (2 * a - 2))
        - K**a * (lam * Dconst * (K * N) ** (a - 2) + Cconst)
    )
    same(eliminated, expected, f"3x3 elimination a={a}")

# a=1 exact Wronskian.  These formulas retain all four nonzero endpoint
# constants and therefore also check the nu=1 zero branch.
d_sym = sp.symbols("d", integer=True, positive=True)
n_sym = d_sym + 1
m_sym = d_sym - 1
A0, B0, L0, M0 = sp.symbols("A0 B0 L0 M0", nonzero=True)
nu = A0 * M0 / (B0 * L0)
cc, dd = sp.symbols("cc dd")
alpha_expr = (cc * K + (A0 / B0) * dd / N) / (1 - nu)
beta_expr = ((M0 / L0) * cc * K + dd / N) / (1 - nu)
middle_wronskian = bracket(-1, alpha_expr, -1, beta_expr, w)
expected_wronskian = cc * dd / (1 - nu) * sp.diff(K * N, w) / N**2
same(middle_wronskian, expected_wronskian, "3x3 a=1 Wronskian")

endpoint_scalar = bracket(
    -n_sym, A0 * K**n_sym, m_sym, M0 * N**m_sym, w
) + bracket(
    m_sym, L0 * N**m_sym, -n_sym, B0 * K**n_sym, w
)
expected_endpoint = B0 * L0 * (1 - nu) * n_sym * m_sym * (
    K ** (n_sym - 1) * N ** (m_sym - 1) * sp.diff(K * N, w)
)
same(endpoint_scalar, expected_endpoint, "3x3 a=1 endpoint scalar")

total_a1 = sp.factor(endpoint_scalar + middle_wronskian)
expected_a1 = sp.diff(K * N, w) * (
    B0 * L0 * (1 - nu) * n_sym * m_sym *
    K ** (n_sym - 1) * N ** (m_sym - 1)
    + cc * dd / ((1 - nu) * N**2)
)
same(total_a1, expected_a1, "3x3 a=1 complete factor")


# Negative-divisor and zero-weight hostile controls.
Delta = w**3 - sp.symbols("c", nonzero=True)**3
u0, v0, u1, v1 = sp.symbols("u0 v0 u1 v1")
delta_prime = sp.diff(Delta, w)
for left, right in [(-1, -1), (-2, -1), (-4, -7)]:
    value = sp.expand(
        left * Delta * u0 * (delta_prime * v0 + Delta * v1)
        - right * (delta_prime * u0 + Delta * u1) * Delta * v0
    )
    same(sp.rem(value, Delta, w), 0,
         f"negative arm divisibility {left},{right}")
for negative in (-1, -2, -5):
    zero_left = sp.expand(negative * Delta * u0 * v1)
    zero_right = sp.expand(negative * v1 * Delta * u0)
    same(sp.rem(zero_left, Delta, w), 0,
         f"zero seam left {negative}")
    same(sp.rem(zero_right, Delta, w), 0,
         f"zero seam right {negative}")


semantic = {
    "3x3": "all common-step AP cells with central scalar bucket are empty; other 3x3 cells remain open",
    "a1": "lower 2x3 ODE seam and 3x3 radical Wronskian seam closed exactly",
    "collision": "2x3 collisions iff r=s or r=t or t=r+s; double iff common AP",
    "common_2x3": "both scalar placements empty; resonances d=3A-2 and 3A=2d+2 also fail",
    "field": "algebraically closed characteristic zero",
    "grading": "B_u=x^u*w^rho(u)*(w^3-c^3)^max(0,ceil(-u/3))*k[w^3]; bracket shift +2",
    "irregular_2x3": "all three single-collision grammars empty",
    "open": "2x4 and noncentral or non-AP 3x3",
    "surface": "R^2E-Z^3+c^3R=0",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
packet_blob = "\n".join(collision_packet).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3787-cubic-pseudoplane-complete-low-support-darboux-nonentry")
print("surface=R^2E-Z^3+c^3R;weights=(3,1,-3);bracket_shift=2")
print("homogeneous=B_u=x^u*w^rho*Delta^ceil_negative_thirds*k[w^3]")
print("collision_counts=" + ",".join(
    f"{key}:{collision_counts[key]}" for key in sorted(collision_counts)
))
print("irregular_controls=" + ",".join(
    f"{key}:{irregular_counts[key]}" for key in sorted(irregular_counts)
))
print(f"lower_resonances={len(lower_resonances)};first={lower_resonances[:4]}")
print(f"upper_resonances={len(upper_resonances)};first={upper_resonances[:4]}")
print(f"exception_controls=lower:{lower_exception_controls},upper:{upper_exception_controls}")
print(f"three_by_three_integration_rows={ap_rows}")
print("closed=all_2x3_and_3x2;central_common_AP_3x3")
print("open=2x4;noncentral_or_non_AP_3x3;arbitrary_Darboux")
print(f"packet_sha256={hashlib.sha256(packet_blob).hexdigest()}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
