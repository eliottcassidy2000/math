#!/usr/bin/env python3
"""Exact controls for THM-3552's two-channel Kummer obstruction.

The paper proof gives the all-parameter valuation argument.  This companion
freezes the first C5 example, its gradient ideal, Kummer branch/divisor
ledger, bounded mate systems, a broken-gradient hostile, and a tame positive
response control.  Explicit exceptions remain active under ``python -O``.
"""

from __future__ import annotations

import hashlib
from math import gcd

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def jac(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(f, x) * sp.diff(g, y) - sp.diff(f, y) * sp.diff(g, x))


def convex_hull(points: list[tuple[int, int]]) -> list[tuple[int, int]]:
    points = sorted(set(points))

    def cross(o: tuple[int, int], a: tuple[int, int], b: tuple[int, int]) -> int:
        return ((a[0] - o[0]) * (b[1] - o[1])
                - (a[1] - o[1]) * (b[0] - o[0]))

    lower: list[tuple[int, int]] = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[int, int]] = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def rank_packet(P: sp.Expr, cap: int) -> tuple[int, int, int, int, int, int]:
    mons = [x**i * y**j for i in range(cap + 1) for j in range(cap + 1 - i)]
    images = [sp.Poly(jac(P, mon), x, y) for mon in mons]
    rhs = sp.Poly(1, x, y)
    output = sorted(set().union(*(poly.monoms() for poly in images), rhs.monoms()))
    matrix = sp.MutableSparseMatrix(
        len(output), len(mons),
        {(i, j): poly.coeff_monomial(mon)
         for j, poly in enumerate(images)
         for i, mon in enumerate(output)
         if poly.coeff_monomial(mon) != 0},
    )
    column = sp.Matrix([rhs.coeff_monomial(mon) for mon in output])
    rank = matrix.rank()
    augmented = matrix.row_join(column).rank()
    return cap, len(mons), len(output), rank, augmented, len(mons) - rank


x, y, T, Z = sp.symbols("x y T Z")
p, q = 2, 5
TT = x**p * y**q
phi = 1 + TT + TT**2
Psi = TT + sp.Rational(1, 2) * TT**2 + sp.Rational(1, 3) * TT**3
P = sp.expand(x * phi + Psi)

# Differentiate in an independent symbol to avoid treating TT as a generator.
phi_T = 1 + T + T**2
Psi_T = T + sp.Rational(1, 2) * T**2 + sp.Rational(1, 3) * T**3
require(sp.diff(Psi_T, T) == phi_T, "Psi prime does not equal phi")

gradient_basis = sp.groebner([sp.diff(P, x), sp.diff(P, y)], x, y, domain=sp.QQ)
require(list(gradient_basis) == [1], "C5 first coordinate is not a submersion")
require(sp.expand(jac(P, TT) - q * phi * TT / y) == 0,
        "generic-fibre response identity")

support = sorted(mon for mon, _ in sp.Poly(P, x, y).terms())
hull = convex_hull(support)
twice_area = abs(sum(
    hull[i][0] * hull[(i + 1) % len(hull)][1]
    - hull[i][1] * hull[(i + 1) % len(hull)][0]
    for i in range(len(hull))
))
require(hull == [(1, 0), (5, 10), (6, 15), (2, 5)], "Newton hull")
require(twice_area == 20, "Newton area")

phi_disc = sp.discriminant(phi_T, T)
moving_resultant = sp.factor(sp.resultant(Z - Psi_T, phi_T, T))
moving_discriminant = sp.factor(sp.discriminant(Z - Psi_T, T))
require(phi_disc == -3, "phi squarefreeness")
require(moving_resultant != 0 and moving_discriminant != 0,
        "generic branch separation")

# y^5=T*phi(T)^2/(Z-Psi(T))^2.  All seven branch orders are coprime to five.
branch_orders = (1, 2, 2, -2, -2, -2, 1)
require(all(gcd(q, abs(order)) == 1 for order in branch_orders),
        "branch is not fully ramified")
genus = (2 - 2 * q + len(branch_orders) * (q - 1)) // 2
require(genus == 10, "Kummer genus")

omega_orders = (("T=0", 0, 1), ("phi=0", 1, 2),
                ("Z-Psi=0", 2, 3), ("infinity", 10, 1))
require(all(order >= 0 for _, order, _ in omega_orders), "omega has a pole")
canonical_degree = sum(order * multiplicity for _, order, multiplicity in omega_orders)
require(canonical_degree == 2 * genus - 2 == 18, "omega divisor degree")

caps = (1, 2, 3, 4, 5, 6, 8, 10, 12, 15)
rank_packets = tuple(rank_packet(P, cap) for cap in caps)
require(all(rank + 1 == augmented for _, _, _, rank, augmented, _ in rank_packets),
        "bounded mate unexpectedly exists")

# Hostile: if Psi'=1 is not divisible by phi, the gradient ideal is proper.
P_bad = sp.expand(x * phi + TT)
bad_basis = sp.groebner([sp.diff(P_bad, x), sp.diff(P_bad, y)], x, y, domain=sp.QQ)
require(list(bad_basis) != [1], "broken divisibility hostile remained a submersion")

# Positive response control for the linear datum.
P_pos = 2 * x + 3
Q_pos = sp.Rational(1, 2) * y + P_pos**3
require(jac(P_pos, Q_pos) == 1, "positive mate control")

general_ledgers = []
for pp, qq, d, e in ((2, 3, 1, 2), (2, 5, 2, 3), (3, 5, 1, 3), (4, 7, 3, 5)):
    require(gcd(pp, qq) == 1 and 2 <= pp < qq and e >= d + 1,
            "general ledger parameters")
    g_inf = gcd(qq, pp * (e - d) - 1)
    orders = (
        0,
        pp - 1,
        qq - 1 - pp,
        (pp * (e - d) - 1 + d * qq) // g_inf - 1,
    )
    require(all(order >= 0 for order in orders), "general omega valuation")
    general_ledgers.append((pp, qq, d, e, g_inf, orders))

packet = (
    sp.factor(P), tuple(support), tuple(hull), twice_area, phi_disc,
    moving_resultant, moving_discriminant, branch_orders, genus,
    omega_orders, rank_packets, tuple(general_ledgers), tuple(bad_basis),
)
semantic = hashlib.sha256(repr(packet).encode("utf-8")).hexdigest()

print("THM3552 TWO-CHANNEL KUMMER OBSTRUCTION EXACT CONTROLS")
print(f"P={P}")
print(f"SUPPORT={support} HULL={hull} NORMALIZED_AREA={twice_area}")
print("GRADIENT ideal=(1) JAC_PT=5*phi*T/y")
print(f"SEPARATION disc_phi={phi_disc} resultant={moving_resultant}")
print(f"KUMMER branch_orders={branch_orders} genus={genus}")
print(f"OMEGA orders={omega_orders} divisor_degree={canonical_degree} HOLomorphic=PASS")
print(f"BOUNDED packets={rank_packets}")
print(f"GENERAL ledgers={tuple(general_ledgers)}")
print("HOSTILE Psi_prime_not_divisible_by_phi=CRITICAL POSITIVE linear_mate=PASS")
print(f"SEMANTIC_SHA256={semantic}")
print("PASS")
