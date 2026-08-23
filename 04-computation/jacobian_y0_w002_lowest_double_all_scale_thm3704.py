#!/usr/bin/env python3
"""Exact symbolic identity audit for THM-3704."""

from collections import defaultdict

import sympy as sp


def require(condition, message):
    if not condition:
        raise AssertionError(message)


# The W002 fibre word is scale-independent; record sums in units of n.
A_SUPPORT = (0, 1, 2)
B_SUPPORT = (0, 1, 2, 4)
fibres = defaultdict(list)
for i, left in enumerate(A_SUPPORT):
    for j, right in enumerate(B_SUPPORT):
        fibres[left + right].append((i, j))
FIBRES = tuple(tuple(fibres[key]) for key in sorted(fibres))
EXPECTED_FIBRES = (
    ((0, 0),),
    ((0, 1), (1, 0)),
    ((0, 2), (1, 1), (2, 0)),
    ((1, 2), (2, 1)),
    ((0, 3), (2, 2)),
    ((1, 3),),
    ((2, 3),),
)
require(FIBRES == EXPECTED_FIBRES, "wrong W002 fibre word")


def same_sign_or_both_zero(r, s):
    return (r == 0 and s == 0) or (r > 0 and s > 0) or (r < 0 and s < 0)


def weights_a(n):
    return (-2, n - 2, 2 * n - 2), (1 - n, 1, n + 1, 3 * n + 1)


def weights_b(n):
    return (1 - n, 1, n + 1), (-2, n - 2, 2 * n - 2, 4 * n - 2)


singleton_cells = ((0, 0), (1, 3), (2, 3))
pa1, qa1 = weights_a(1)
pb1, qb1 = weights_b(1)
pa2, qa2 = weights_a(2)
pb2, qb2 = weights_b(2)
require(not same_sign_or_both_zero(pa1[0], qa1[0]), "A n=1 should fail 00")
require(not same_sign_or_both_zero(pb1[0], qb1[0]), "B n=1 should fail 00")
require(not same_sign_or_both_zero(pa2[1], qa2[3]), "A n=2 should fail 13")
require(
    all(same_sign_or_both_zero(pb2[i], qb2[j]) for i, j in singleton_cells),
    "B n=2 should survive singleton signs",
)


# Work in the differential polynomial ring C[H,K,H',K']; the derivation below
# is exact on every monomial used by the proof.
ell = sp.symbols("ell", integer=True, positive=True)
H, K, Hp, Kp = sp.symbols("H K Hp Kp", nonzero=True)
a, b0, c0, d0, t0, kappa, lam = sp.symbols(
    "a b0 c0 d0 t0 kappa lam", nonzero=True
)


def derivative(expression):
    return sp.diff(expression, H) * Hp + sp.diff(expression, K) * Kp


def wedge(r, left, s, right):
    return sp.expand(s * derivative(left) * right - r * left * derivative(right))


def zero(expression, label):
    reduced = sp.simplify(sp.powsimp(sp.factor(expression), force=True))
    require(reduced == 0, label)


parity_branches = (
    (1, 2 * ell, 2, 2 * ell - 1),
    (2, 2 * ell + 1, 1, ell),
)


def audit_orientation_a(delta, n, epsilon, m):
    F = H**epsilon
    G0 = H**m
    gamma = a * t0 * (3 * n + 1) / (d0 * (2 * n - 2))
    beta = c0 * (n - 2) * gamma / (d0 * (2 * n - 2))
    M = kappa * K ** (n + 1) + gamma * F * K ** (n + 3)
    L = lam * K - beta * F * K**3

    top = wedge(-2, a * F, 3 * n + 1, t0 * K ** (3 * n + 1))
    top += wedge(2 * n - 2, d0 * K ** (2 * n - 2), n + 1, M)
    next_row = wedge(n - 2, c0 * K ** (n - 2), n + 1, M)
    next_row += wedge(2 * n - 2, d0 * K ** (2 * n - 2), 1, L)

    euler = Hp * K + delta * H * Kp
    scalar = wedge(-2, a * F, 1, L)
    scalar += wedge(n - 2, c0 * K ** (n - 2), 1 - n, b0 * G0)
    quotient = a * epsilon * H ** (epsilon - 1) * (
        lam - 3 * beta * F * K**2
    )
    quotient -= b0 * c0 * (n - 2) * m * H ** (m - 1) * K ** (n - 3)

    zero(top, f"orientation A top row delta={delta}")
    zero(next_row, f"orientation A next row delta={delta}")
    zero(scalar - euler * quotient, f"orientation A scalar delta={delta}")


def audit_orientation_b(delta, n, epsilon, m):
    F = H**m
    G0 = H**epsilon
    gamma = a * t0 * (4 * n - 2) / (d0 * (n + 1))
    beta = c0 * gamma / (d0 * (n + 1))
    M = kappa * K ** (2 * n - 2) + gamma * F * K ** (3 * n - 3)
    L = lam * K ** (n - 2) - beta * F * K ** (2 * n - 3)

    top = wedge(1 - n, a * F, 4 * n - 2, t0 * K ** (4 * n - 2))
    top += wedge(n + 1, d0 * K ** (n + 1), 2 * n - 2, M)
    next_row = wedge(1, c0 * K, 2 * n - 2, M)
    next_row += wedge(n + 1, d0 * K ** (n + 1), n - 2, L)

    euler = Hp * K + delta * H * Kp
    scalar = wedge(1 - n, a * F, n - 2, L)
    scalar += wedge(1, c0 * K, -2, b0 * G0)
    quotient = a * m * H ** (m - 1) * (
        lam * (n - 2) * K ** (n - 3)
        - beta * (2 * n - 3) * F * K ** (2 * n - 4)
    )
    quotient -= b0 * c0 * epsilon * H ** (epsilon - 1)

    zero(top, f"orientation B top row delta={delta}")
    zero(next_row, f"orientation B next row delta={delta}")
    zero(scalar - euler * quotient, f"orientation B scalar delta={delta}")


for branch in parity_branches:
    audit_orientation_a(*branch)
    audit_orientation_b(*branch)


# The active small boundary B,n=2 avoids the formal zero*K^-1 expression.
n = sp.Integer(2)
gamma = a * t0 * (4 * n - 2) / (d0 * (n + 1))
beta = c0 * gamma / (d0 * (n + 1))
M = kappa * K**2 + gamma * H * K**3
L = lam - beta * H * K
euler = Hp * K + H * Kp
top = wedge(-1, a * H, 6, t0 * K**6) + wedge(3, d0 * K**3, 2, M)
next_row = wedge(1, c0 * K, 2, M) + wedge(3, d0 * K**3, 0, L)
scalar = wedge(-1, a * H, 0, L) + wedge(1, c0 * K, -2, b0 * H**2)
zero(top, "orientation B n=2 top row")
zero(next_row, "orientation B n=2 next row")
zero(scalar + H * euler * (a * beta + 2 * b0 * c0), "orientation B n=2 scalar")


print("THM-3704 exact W002 audit")
print("fibres = 00|01=10|02=11=20|12=21|03=22|13|23")
print("A n=1 gate = 00 zero/nonzero FAIL")
print("B n=1 gate = 00 zero/nonzero FAIL")
print("A n=2 gate = 13 zero/nonzero FAIL")
print("B n=2 singleton gate = PASS; Euler factor identity = PASS")
print("orientation A delta=1,2 integrated rows and scalar factors = PASS")
print("orientation B delta=1,2 integrated rows and scalar factors = PASS")
print("common factor = H'K + delta H K'; h|H and b|K make it nonunit")
print("scope = scalar fibre 01=10 in both orientations, not all W002 placements")
print("ALL CHECKS PASSED")
