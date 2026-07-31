#!/usr/bin/env python3
"""Exact companion for THM-2971.

Verify the explicit edge/orientation sextic algebra isomorphism on the
quartic discriminant cover, its exact S4/A4 descent typing, and its sharp
separability and base-field boundaries.
"""

from itertools import permutations

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def parity(perm):
    return sum(
        perm[i] > perm[j]
        for i in range(len(perm))
        for j in range(i + 1, len(perm))
    ) % 2


def rotate_to_zero(cycle):
    j = cycle.index(0)
    return cycle[j:] + cycle[:j]


def relabel_cycle(g, cycle):
    return rotate_to_zero(tuple(g[i] for i in cycle))


def relabel_edge(g, edge):
    return tuple(sorted((g[edge[0]], g[edge[1]])))


def signed_pair_action(pairs, relabel, g):
    lookup = {}
    for j, pair in enumerate(pairs):
        lookup[pair[0]] = (j, 0)
        lookup[pair[1]] = (j, 1)
    sigma = [None] * 3
    delta = [None] * 3
    for i, pair in enumerate(pairs):
        target = relabel(g, pair[0])
        require(target in lookup, "signed-pair relabelling left the carrier")
        j, bit = lookup[target]
        require(
            relabel(g, pair[1]) == pairs[j][bit ^ 1],
            "signed-pair inverse mate changed",
        )
        sigma[i] = j
        delta[i] = bit
    require(sorted(sigma) == [0, 1, 2], "block action is not a permutation")
    return tuple(delta), tuple(sigma)


def central_twist(action):
    delta, sigma = action
    return tuple(bit ^ 1 for bit in delta), sigma


def permutation_order(g):
    power = tuple(range(len(g)))
    for order in range(1, 25):
        power = tuple(g[power[i]] for i in range(len(g)))
        if power == tuple(range(len(g))):
            return order
    raise RuntimeError("permutation order exceeded S4 bound")


S4 = tuple(permutations(range(4)))
A4 = tuple(g for g in S4 if parity(g) == 0)
edge_pairs = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)
orientation_pairs = (
    ((0, 2, 1, 3), (0, 3, 1, 2)),
    ((0, 3, 2, 1), (0, 1, 2, 3)),
    ((0, 1, 3, 2), (0, 2, 3, 1)),
)

for g in S4:
    edge_action = signed_pair_action(edge_pairs, relabel_edge, g)
    orientation_action = signed_pair_action(
        orientation_pairs, relabel_cycle, g
    )
    require(
        orientation_action
        == (edge_action if parity(g) == 0 else central_twist(edge_action)),
        "central-sign graph changed",
    )

edge0 = edge_pairs[0][0]
orientation0 = orientation_pairs[0][0]
edge_stabilizer = {
    g for g in S4 if relabel_edge(g, edge0) == edge0
}
orientation_stabilizer = {
    g for g in S4 if relabel_cycle(g, orientation0) == orientation0
}
require(len(edge_stabilizer) == len(orientation_stabilizer) == 4, "stabilizer order changed")
require(
    sorted(permutation_order(g) for g in edge_stabilizer) == [1, 2, 2, 2],
    "edge stabilizer is not V4",
)
require(
    sorted(permutation_order(g) for g in orientation_stabilizer)
    == [1, 2, 4, 4],
    "orientation stabilizer is not C4",
)
require(
    edge_stabilizer.intersection(A4)
    == orientation_stabilizer.intersection(A4),
    "the two selected stabilizers do not agree over A4",
)
require(
    len(edge_stabilizer.intersection(A4)) == 2,
    "common A4 stabilizer is not C2",
)

transposition = (0, 1, 3, 2)
edge_fixed = sum(
    relabel_edge(transposition, edge) == edge
    for pair in edge_pairs
    for edge in pair
)
orientation_fixed = sum(
    relabel_cycle(transposition, cycle) == cycle
    for pair in orientation_pairs
    for cycle in pair
)
require((edge_fixed, orientation_fixed) == (2, 0), "base-field character hostile changed")


p, q, r, U, Y, Z = sp.symbols("p q r U Y Z")
Delta = (
    256 * r**3
    - 128 * p**2 * r**2
    + 144 * p * q**2 * r
    - 27 * q**4
    + 16 * p**4 * r
    - 4 * p**3 * q**2
)
J_or = 1024 * r**3 + 768 * p**2 * r**2 - 288 * p * q**2 * r + 27 * q**4
S = U**3 + 2 * p * U**2 + (p**2 - 4 * r) * U - q**2
C = sp.diff(S, U)
F = 16 * (16 * r - 4 * p * U - 3 * U**2) * (
    U**2 + 2 * p * U + p**2 - 4 * r
)
E = sp.expand(S.subs(U, Y**2))
Eprime = sp.diff(E, Y)
K = sp.QQ.frac_field(p, q, r)

numerator = (
    16 * Y**4 * p**2 * r
    - 6 * Y**4 * p * q**2
    - 64 * Y**4 * r**2
    + 32 * Y**2 * p**3 * r
    - 10 * Y**2 * p**2 * q**2
    - 128 * Y**2 * p * r**2
    + 24 * Y**2 * q**2 * r
    + 16 * p**4 * r
    - 4 * p**3 * q**2
    - 128 * p**2 * r**2
    + 80 * p * q**2 * r
    - 9 * q**4
    + 256 * r**3
)
A = sp.cancel(-4 * Y * numerator / (q * Delta))
require(
    sp.cancel(sp.rem(Eprime * A + 8 * q, E, domain=K)) == 0,
    "explicit inverse-derivative formula failed",
)
require(
    sp.cancel(
        sp.rem(Delta * A**2 - F.subs(U, Y**2), E, domain=K)
    )
    == 0,
    "intertwiner square did not equal F(s^2)",
)
require(sp.expand(A.subs(Y, -Y) + A) == 0, "intertwiner is not odd")
B_in_Y2 = sp.cancel(A / Y)
require(
    sp.expand(B_in_Y2.subs(Y, -Y) - B_in_Y2) == 0,
    "odd intertwiner quotient is not even",
)
B = sum(
    coefficient * U ** (exponent[0] // 2)
    for exponent, coefficient in sp.Poly(B_in_Y2, Y).terms()
)
require(sp.degree(B, U) == 2, "cubic-base Kummer multiplier degree changed")
require(
    sp.cancel(sp.rem(Delta * U * B**2 - F, S, domain=K)) == 0,
    "cubic-base Kummer ratio changed",
)

orientation_coefficient_four = -32 * (8 * p**2 * r - 3 * p * q**2 - 32 * r**2)
orientation_coefficient_two = (
    -256 * q**2 * (p**2 + 12 * r) * (32 * p * r - 9 * q**2)
)
orientation_coefficient_zero = -4096 * q**4 * Delta
O = (
    Z**6
    + orientation_coefficient_four * Z**4
    + orientation_coefficient_two * Z**2
    + orientation_coefficient_zero
)
O_poly = sp.Poly(O, Z)
require(O_poly.degree() == 6, "orientation resultant degree changed")
require(
    all(exponent[0] % 2 == 0 for exponent, _ in O_poly.terms()),
    "orientation sextic is not even",
)
require(
    sp.rem(
        F**3
        + orientation_coefficient_four * F**2
        + orientation_coefficient_two * F
        + orientation_coefficient_zero,
        S,
        domain=K,
    )
    == 0,
    "orientation norm polynomial did not vanish at F(u)",
)
disc_edge = 64 * q**2 * Delta**2
disc_orientation = 2**66 * q**12 * J_or**4 * Delta**3
change_determinant_square = (2**30 * q**5 * J_or**2) ** 2 * Delta
require(
    sp.expand(disc_orientation - disc_edge * change_determinant_square)
    == 0,
    "power-basis discriminant ratio changed",
)

a0, a1, a2 = sp.symbols("a0 a1 a2")
a3 = -a0 - a1 - a2
symbolic_roots = (a0, a1, a2, a3)
p_from_roots = sum(
    symbolic_roots[i] * symbolic_roots[j]
    for i in range(4)
    for j in range(i + 1, 4)
)
q_from_roots = -sum(
    symbolic_roots[i] * symbolic_roots[j] * symbolic_roots[k]
    for i in range(4)
    for j in range(i + 1, 4)
    for k in range(j + 1, 4)
)
r_from_roots = sp.prod(symbolic_roots)
s_symbolic = a0 + a1
d1_symbolic = a0 - a1
d2_symbolic = a2 - a3
vandermonde_symbolic = sp.prod(
    symbolic_roots[i] - symbolic_roots[j]
    for i in range(4)
    for j in range(i + 1, 4)
)
Sprime_at_edge = C.subs(
    {p: p_from_roots, q: q_from_roots, r: r_from_roots, U: s_symbolic**2}
)
require(
    sp.expand(
        s_symbolic * (d1_symbolic**2 - d2_symbolic**2)
        + 4 * q_from_roots
    )
    == 0,
    "root-level q identity changed",
)
require(
    sp.expand(
        vandermonde_symbolic
        - d1_symbolic * d2_symbolic * Sprime_at_edge
    )
    == 0,
    "root-level Vandermonde identity changed",
)


def omega(cycle, roots):
    i, j, k, ell = cycle
    d1 = roots[i] - roots[k]
    d2 = roots[j] - roots[ell]
    return d1 * d2 * (d1**2 - d2**2)


roots = (-4, -1, 2, 3)
p0 = sum(roots[i] * roots[j] for i in range(4) for j in range(i + 1, 4))
e3 = sum(
    roots[i] * roots[j] * roots[k]
    for i in range(4)
    for j in range(i + 1, 4)
    for k in range(j + 1, 4)
)
q0 = -e3
r0 = sp.prod(roots)
V0 = sp.prod(
    roots[i] - roots[j]
    for i in range(4)
    for j in range(i + 1, 4)
)
require(V0**2 == Delta.subs({p: p0, q: q0, r: r0}), "Vandermonde square failed")
A0 = sp.cancel(A.subs({p: p0, q: q0, r: r0}))
root_controls = []
for edge_pair, orientation_pair in zip(edge_pairs, orientation_pairs):
    for edge, cycle in zip(edge_pair, orientation_pair):
        s0 = roots[edge[0]] + roots[edge[1]]
        value = sp.cancel(V0 * A0.subs(Y, s0))
        expected = omega(cycle, roots)
        require(value == expected, "common-gauge root formula changed")
        root_controls.append((s0, expected))

E0 = sp.Poly(E.subs({p: p0, q: q0, r: r0}), Y)
T0 = sp.Poly(V0 * A0, Y)
change_columns = []
for exponent in range(6):
    remainder = sp.Poly(
        sp.rem(T0.as_expr() ** exponent, E0.as_expr(), Y), Y
    )
    change_columns.append([remainder.nth(i) for i in range(6)])
change_matrix = sp.Matrix(
    6, 6, lambda row, column: change_columns[column][row]
)
change_determinant_control = change_matrix.det()
change_determinant_expected = (
    2**30
    * q0**5
    * J_or.subs({p: p0, q: q0, r: r0}) ** 2
    * V0
)
require(
    change_determinant_control == change_determinant_expected,
    "power-basis determinant sign or normalization changed",
)

boundary_O = sp.factor(O.subs({p: 1, q: 4, r: -3}))
require(
    boundary_O == (Z**2 - 1280) ** 2 * (Z**2 + 14080),
    "J_or=0 boundary factorization changed",
)
require(
    Delta.subs({p: 1, q: 4, r: -3}) == -22000
    and J_or.subs({p: 1, q: 4, r: -3}) == 0,
    "J_or boundary typing changed",
)

print("THM-2971 discriminant-cover sextic algebra intertwiner")
print("group_orders=|S4|:24;|A4|:12;edge_stabilizer:V4;orientation_stabilizer:C4")
print("A4_common_stabilizer=C2;S4_transposition_fixed_counts=edge:2/orientation:0")
print(f"symbolic_degrees=S:{sp.degree(S,U)};E:{sp.degree(E,Y)};O:{sp.degree(O,Z)};A:{sp.degree(A,Y)}")
print("intertwiner=T(Y)=sqrt(Delta)*A(Y)=-8*q*sqrt(Delta)/E'(Y) mod E")
print("kummer_gate=F(U)=Delta*U*B(U)^2 mod S;degree_B:2")
print("root_formula_symbolic=s(d1^2-d2^2):-4q;V:d1*d2*Sprime(s^2);PASS")
print("power_basis_det=2^30*q^5*J_or^2*sqrt(Delta);discriminant_ratio:PASS")
print("symbolic_gates=Eprime_inverse:PASS;T_square:PASS;O_of_T:PASS;oddness:PASS")
print(f"numeric_control=(p,q,r,Delta,V)=({p0},{q0},{r0},{V0**2},{V0})")
print(f"six_root_controls={tuple(root_controls)}")
print("boundary=(p,q,r)=(1,4,-3);Delta=-22000;J_or=0;O=(Z^2-1280)^2(Z^2+14080)")
print("all_exact_controls=PASS")
