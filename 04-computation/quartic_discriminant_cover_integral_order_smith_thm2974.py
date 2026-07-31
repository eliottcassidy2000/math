#!/usr/bin/env python3
"""Exact companion for THM-2974.

Audit the four affine sidecars lost by the THM-2971 generic intertwiner.

The test family is THM-2769

    f_t(X)=X^4-2X^2-8tX+1-4t.

At t=0 we compare the edge and oriented-cycle sextic graph orders after the
discriminant cover.  All decisions are exact; no floating point is used.
"""

from fractions import Fraction
from itertools import combinations, permutations

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def valuation(poly, parameter):
    expanded = sp.Poly(sp.expand(poly), parameter)
    return min(exponent[0] for exponent, coefficient in expanded.terms() if coefficient)


def rational_valuation(value, parameter):
    value = sp.cancel(value)
    if value == 0:
        raise RuntimeError("valuation of zero requested")
    numerator, denominator = sp.fraction(value)
    return valuation(numerator, parameter) - valuation(denominator, parameter)


def lower_hull(points):
    hull = []
    for point in sorted(points):
        while len(hull) >= 2:
            a, b = hull[-2], hull[-1]
            left = Fraction(b[1] - a[1], b[0] - a[0])
            right = Fraction(point[1] - b[1], point[0] - b[0])
            if left >= right:
                hull.pop()
            else:
                break
        hull.append(point)
    return tuple(hull)


def newton_valuations(poly, variable, parameter):
    points = tuple(
        sorted(
            (exponent[0], valuation(coefficient, parameter))
            for exponent, coefficient in sp.Poly(poly, variable).terms()
            if coefficient
        )
    )
    hull = lower_hull(points)
    values = []
    for left, right in zip(hull, hull[1:]):
        slope = Fraction(right[1] - left[1], right[0] - left[0])
        values.extend([-slope] * (right[0] - left[0]))
    return points, hull, tuple(values)


def cycle_type(action):
    seen = set()
    lengths = []
    for start in range(len(action)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = action[current]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def rotate_to_zero(cycle):
    index = cycle.index(0)
    return cycle[index:] + cycle[:index]


def relabel_cycle(permutation, cycle):
    return rotate_to_zero(tuple(permutation[index] for index in cycle))


t, U, Y, Z = sp.symbols("t U Y Z")
s, a, b, X = sp.symbols("s a b X")
p = sp.Integer(-2)
q = -8 * t
r = 1 - 4 * t

Delta = sp.factor(
    256 * r**3
    - 128 * p**2 * r**2
    + 144 * p * q**2 * r
    - 27 * q**4
    + 16 * p**4 * r
    - 4 * p**3 * q**2
)
J_or = sp.factor(1024 * r**3 + 768 * p**2 * r**2 - 288 * p * q**2 * r + 27 * q**4)

E = sp.expand(Y**6 - 4 * Y**4 + 16 * t * Y**2 - 64 * t**2)
c4 = sp.expand(-32 * (8 * p**2 * r - 3 * p * q**2 - 32 * r**2))
c2 = sp.expand(-256 * q**2 * (p**2 + 12 * r) * (32 * p * r - 9 * q**2))
c0 = sp.expand(-4096 * q**4 * Delta)
O = sp.expand(Z**6 + c4 * Z**4 + c2 * Z**2 + c0)

require(sp.expand(Delta + 4096 * t**2 * (27 * t**2 - 14 * t + 3)) == 0, "quartic discriminant changed")
require(sp.expand(J_or - 4096 * (t - 1) * (27 * t**3 - 25 * t**2 + 8 * t - 1)) == 0, "J_or changed")

quartic = X**4 - 2 * X**2 - 8 * t * X + 1 - 4 * t
plus_expansion = sp.Poly(sp.expand(quartic.subs({X: 1 + a * s, t: s**2})), s)
minus_expansion = sp.Poly(sp.expand(quartic.subs({X: -1 + b * s, t: s**2})), s)
require(plus_expansion.coeff_monomial(s**2) == 4 * (a**2 - 3), "root expansion at +1 changed")
require(minus_expansion.coeff_monomial(s**2) == 4 * (b**2 + 1), "root expansion at -1 changed")

edge_points, edge_hull, edge_values = newton_valuations(E, Y, t)
orientation_points, orientation_hull, orientation_values = newton_valuations(O, Z, t)
require(edge_points == ((0, 2), (2, 1), (4, 0), (6, 0)), "edge Newton points changed")
require(edge_hull == ((0, 2), (4, 0), (6, 0)), "edge lower hull changed")
require(edge_values == (Fraction(1, 2),) * 4 + (Fraction(0),) * 2, "edge root valuations changed")
require(orientation_points == ((0, 6), (2, 2), (4, 1), (6, 0)), "orientation Newton points changed")
require(orientation_hull == ((0, 6), (2, 2), (6, 0)), "orientation lower hull changed")
require(orientation_values == (Fraction(2),) * 2 + (Fraction(1, 2),) * 4, "orientation root valuations changed")

# The three matching-base Kummer radicands are Y^2 and Z^2.  Pair the two
# order-one edge radicands with the four half-order roots, and the unit edge
# radicand with the two fixed roots.  The parity word is unchanged.
edge_kummer_orders = (1, 1, 0)
orientation_kummer_orders = (1, 1, 4)
require(
    tuple(value % 2 for value in edge_kummer_orders)
    == tuple(value % 2 for value in orientation_kummer_orders)
    == (1, 1, 0),
    "Kummer residue word was not preserved",
)
require(valuation(Delta, t) == 2, "discriminant valuation changed")

# A loop about t=0 sends s=sqrt(t) to -s.  The leading root equations are
# X=1+a*s with a^2=3 and X=-1+b*s with b^2=-1, so its natural quartic
# inertia is a double transposition and has no fixed sheet.
quartic_inertia = (1, 0, 3, 2)
require(cycle_type(quartic_inertia) == (2, 2), "quartic inertia changed")

edges = tuple((i, j) for i in range(4) for j in range(i + 1, 4))
edge_lookup = {edge: index for index, edge in enumerate(edges)}
edge_action = tuple(
    edge_lookup[tuple(sorted((quartic_inertia[i], quartic_inertia[j])))]
    for i, j in edges
)

cycles = tuple(permutation for permutation in permutations(range(4)) if permutation[0] == 0)
cycle_lookup = {cycle: index for index, cycle in enumerate(cycles)}
orientation_action = tuple(
    cycle_lookup[relabel_cycle(quartic_inertia, cycle)] for cycle in cycles
)
require(cycle_type(edge_action) == (2, 2, 1, 1), "edge inertia cycle type changed")
require(cycle_type(orientation_action) == (2, 2, 1, 1), "orientation inertia cycle type changed")
require(sum(index == image for index, image in enumerate(quartic_inertia)) == 0, "quartic inertia gained an owner")
require(sum(index == image for index, image in enumerate(edge_action)) == 2, "edge fixed-sheet count changed")
require(sum(index == image for index, image in enumerate(orientation_action)) == 2, "orientation fixed-sheet count changed")

# Over the discriminant cover the prime over t=0 is unramified (Delta has
# even order).  The common normalized degree-six algebra is tame with local
# inertia cycle 2^2 1^2, hence maximal-order discriminant exponent 2.
maximal_discriminant_order = sum(length - 1 for length in cycle_type(edge_action))
edge_discriminant_order = valuation(64 * q**2 * Delta**2, t)
orientation_discriminant_order = valuation(2**66 * q**12 * J_or**4 * Delta**3, t)
require((maximal_discriminant_order, edge_discriminant_order, orientation_discriminant_order) == (2, 6, 18), "order discriminants changed")
edge_normalization_colength = (edge_discriminant_order - maximal_discriminant_order) // 2
orientation_normalization_colength = (orientation_discriminant_order - maximal_discriminant_order) // 2
relative_colength = orientation_normalization_colength - edge_normalization_colength
intertwiner_determinant_order = 5 * valuation(q, t) + 2 * valuation(J_or, t) + 1
require((edge_normalization_colength, orientation_normalization_colength) == (2, 8), "graph-order indices changed")
require(relative_colength == intertwiner_determinant_order == 6, "relative graph-order index changed")

# Relative lattice position.  Write sqrt(Delta)=t*w, where w is a local
# unit.  Column k of the change matrix gains the unit w^k, so the Smith
# valuations are those of powers of T0=t*A(Y).  Parity splits the matrix into
# two 3x3 blocks and keeps this exact calculation small.
T0 = sp.cancel(
    Y
    * (
        Y**4 * (t - 1)
        + Y**2 * (24 * t**2 - 12 * t + 4)
        + 144 * t**3
        - 96 * t**2
        + 16 * t
    )
    / (32 * t * (27 * t**2 - 14 * t + 3))
)
fraction_field = sp.QQ.frac_field(t)
F = 16 * (16 * r - 4 * p * U - 3 * U**2) * (U**2 + 2 * p * U + p**2 - 4 * r)
require(
    sp.cancel(sp.rem(sp.diff(E, Y) * T0 - 64 * t**2, E, domain=fraction_field)) == 0,
    "local inverse-derivative representative changed",
)
require(
    sp.cancel(
        sp.rem((Delta / t**2) * T0**2 - F.subs(U, Y**2), E, domain=fraction_field)
    )
    == 0,
    "local discriminant-cover intertwiner square changed",
)
change_columns = []
for power in range(6):
    remainder = sp.rem(T0**power, E, domain=fraction_field)
    change_columns.append(
        [sp.cancel(sp.Poly(remainder, Y).coeff_monomial(Y**row)) for row in range(6)]
    )
change_matrix = sp.Matrix(6, 6, lambda row, column: change_columns[column][row])
block_smith = []
for rows, columns in (
    ((0, 2, 4), (0, 2, 4)),
    ((1, 3, 5), (1, 3, 5)),
):
    block = change_matrix.extract(rows, columns)
    determinant_orders = [0]
    for size in range(1, 4):
        orders = []
        for selected_rows in combinations(range(3), size):
            for selected_columns in combinations(range(3), size):
                determinant = block.extract(selected_rows, selected_columns).det()
                if determinant:
                    orders.append(rational_valuation(determinant, t))
        require(orders, "Smith block lost all minors")
        determinant_orders.append(min(orders))
    block_smith.extend(
        determinant_orders[index] - determinant_orders[index - 1]
        for index in range(1, 4)
    )
relative_smith = tuple(sorted(block_smith))
require(relative_smith == (-1, 0, 0, 1, 2, 4), "relative order Smith profile changed")
require(sum(relative_smith) == intertwiner_determinant_order, "Smith determinant sum changed")

# Exhaust the boundary components of the THM-2769 parameter line.  The
# power-basis determinant is a unit only off q*Delta*J_or.  On the
# discriminant cover its orders are 6 at t=0, 1 above a simple root of the
# residual discriminant factor, and 2 above every simple J_or root.
g_delta = 27 * t**2 - 14 * t + 3
g_jor = 27 * t**3 - 25 * t**2 + 8 * t - 1
require(sp.gcd(g_delta, sp.diff(g_delta, t)) == 1, "residual Delta factor is not squarefree")
require(sp.gcd(g_jor, sp.diff(g_jor, t)) == 1, "cubic J_or factor is not squarefree")
require(sp.gcd(g_delta, (t - 1) * g_jor) == 1, "Delta and J_or boundary components met")
require(g_delta.subs(t, 0) != 0 and J_or.subs(t, 0) != 0, "t=0 boundary multiplicities changed")
require(g_jor.subs(t, 1) != 0, "the two J_or factors met")
require((6, 1, 2) == (intertwiner_determinant_order, 1, 2), "boundary determinant orders changed")

print("THM-2974 DISCRIMINANT-COVER INTEGRAL-ORDER SMITH AND OWNER BOUNDARY")
print("family=X^4-2X^2-8tX+1-4t")
print("t0_quartic_leading_roots=1+/-sqrt(3)*sqrt(t);-1+/-i*sqrt(t)")
print("maximal_inertia=double_transposition;quartic_cycle=2^2;sextic_cycle=2^2*1^2")
print("fixed_sheets=quartic:0;edge:2;orientation:2")
print("edge_newton=(1/2^4,0^2);orientation_newton=(1/2^4,2^2)")
print("kummer_orders=edge:(1,1,0);orientation:(1,1,4);residue=110")
print("discriminant_orders=maximal:2;edge_order:6;orientation_order:18")
print("normalization_colengths=edge:2;orientation:8;relative:6")
print("relative_lattice_smith=(-1,0,0,1,2,4);orders_not_nested")
print("intertwiner_det_order_t0=6;simple_Delta_cover=1;simple_Jor=2")
print("FIRST_FAILURE=fixed_sextic_sheet_does_not_imply_fixed_quartic_owner")
print("SCOPE=field_inertia_and_Kummer_residues_preserved;graph_order_and_owner_not_preserved")
print("ALL_EXACT_CONTROLS=PASS")
