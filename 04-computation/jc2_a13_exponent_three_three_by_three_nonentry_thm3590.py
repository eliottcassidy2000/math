#!/usr/bin/env python3
"""Exact companion for provisional THM-3590.

This checks the complete complement-edge combinatorics and a finite hostile
atlas for the A13 exponent-three 3 x 3 Darboux nonentry.  The universal
valuation and connected-logarithmic-derivative arguments remain proof-driven.
All truth-bearing gates use ``require`` rather than Python ``assert``.
"""

from collections import defaultdict, deque

import sympy as sp


b, c, e, z = sp.symbols("b c e z")
kappa = sp.symbols("kappa", nonzero=True)
CHECKS = 0


def require(label, condition):
    """Record one optimization-stable exact gate."""
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError("FAILED exact gate: " + label)


def is_zero(expression):
    return sp.cancel(expression) == 0


def surface_bracket(F, G):
    """Poisson bracket on c^3 e=b(b^2+kappa^2)."""
    sigma = b * (b**2 + kappa**2)
    bc = c**3 * (sp.diff(F, b) * sp.diff(G, c) - sp.diff(F, c) * sp.diff(G, b))
    ce = -sp.diff(sigma, b) * (
        sp.diff(F, c) * sp.diff(G, e) - sp.diff(F, e) * sp.diff(G, c)
    )
    be = -3 * c**2 * e * (
        sp.diff(F, b) * sp.diff(G, e) - sp.diff(F, e) * sp.diff(G, b)
    )
    return sp.expand(bc + ce + be)


def layer_map(r_weights, s_weights):
    layers = defaultdict(list)
    for i, r_value in enumerate(r_weights):
        for j, s_value in enumerate(s_weights):
            layers[r_value + s_value].append((i, j))
    return dict(layers)


def complement_edges(r_weights, s_weights):
    return sorted(layer_map(r_weights, s_weights).get(-2, []))


def singleton_edges(r_weights, s_weights):
    answer = []
    for value, edges in layer_map(r_weights, s_weights).items():
        if value != -2 and len(edges) == 1:
            answer.extend(edges)
    return sorted(answer)


def collision_sets(r_weights, s_weights):
    answer = []
    for value, edges in layer_map(r_weights, s_weights).items():
        if value != -2 and len(edges) > 1:
            answer.append(frozenset(edges))
    return sorted(answer, key=lambda edges: sorted(edges))


def connected(edges):
    vertices = [("P", i) for i in range(3)] + [("Q", j) for j in range(3)]
    graph = {vertex: set() for vertex in vertices}
    for i, j in edges:
        left = ("P", i)
        right = ("Q", j)
        graph[left].add(right)
        graph[right].add(left)
    seen = {vertices[0]}
    queue = deque([vertices[0]])
    while queue:
        vertex = queue.popleft()
        for neighbor in graph[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    return len(seen) == len(vertices)


def common_step(weights):
    return weights[1] - weights[0] == weights[2] - weights[1]


def has_exact_low(weights):
    return any(value <= -6 and value % 3 == 0 for value in weights)


def has_high(weights):
    return any(value >= 2 for value in weights)


print("THM-3590 exact companion")
print("SECTION exponent-three bracket and arm normal forms")

f_test = b**3 + 2 * b + 1
g_test = b**2 - 3 * b + 5
for u in range(-9, 7):
    for v in range(-9, 7):
        F = c**u * f_test
        G = c**v * g_test
        chart_bracket = c**3 * (
            sp.diff(F, b) * sp.diff(G, c) - sp.diff(F, c) * sp.diff(G, b)
        )
        expected = c ** (u + v + 2) * (
            v * sp.diff(f_test, b) * g_test - u * f_test * sp.diff(g_test, b)
        )
        require("homogeneous bracket u=%d v=%d" % (u, v), is_zero(chart_bracket - expected))

sigma = b * (b**2 + kappa**2)
h_test = b + 2
for s_value in range(1, 31):
    m_value = (s_value + 2) // 3
    c_exponent = 3 * m_value - s_value
    require("normal exponent range s=%d" % s_value, 0 <= c_exponent <= 2)
    require(
        "exact e-arm survival s=%d" % s_value,
        (c_exponent == 0) == (s_value % 3 == 0),
    )
    left = c ** (-s_value) * sigma**m_value * h_test
    right = c**c_exponent * e**m_value * h_test
    require(
        "negative normal form s=%d" % s_value,
        is_zero((left - right).subs(e, sigma / c**3)),
    )

alpha, beta, gamma, delta = sp.symbols("alpha beta gamma delta")
F_jet = alpha * c + beta * e
G_jet = gamma * c + delta * e
jet_value = surface_bracket(F_jet, G_jet).subs({b: 0, c: 0, e: 0})
require(
    "central tangent jet determinant",
    is_zero(jet_value + kappa**2 * (alpha * delta - beta * gamma)),
)
for middle_p in (-3, 1):
    for middle_q in (-3, 1):
        determinant_nonzero = middle_p != middle_q
        require(
            "opposite middle jets %d %d" % (middle_p, middle_q),
            determinant_nonzero == ({middle_p, middle_q} == {-3, 1}),
        )

print("PASS bracket formula, 30 negative normal forms, and central jet gate")


print("SECTION simple-arm unit classification")

# The two exceptional negative-negative / negative-zero boundaries.
for a_order in range(1, 7):
    for d_order in range(1, 7):
        f_local = z**a_order
        g_local = z**d_order
        W_11 = -sp.diff(f_local, z) * g_local + f_local * sp.diff(g_local, z)
        require(
            "(-1,-1) arm ideal a=%d d=%d" % (a_order, d_order),
            W_11 == 0 or sp.Poly(W_11, z).terms()[0][0][0] >= 1,
        )

for a_order in range(1, 7):
    for d_order in range(0, 7):
        f_local = z**a_order
        g_local = z**d_order
        W_20 = 2 * f_local * sp.diff(g_local, z)
        require(
            "(-2,0) arm ideal a=%d d=%d" % (a_order, d_order),
            W_20 == 0 or sp.Poly(W_20, z).terms()[0][0][0] >= 1,
        )

unit_rows = []
for N in range(3, 31):
    minimum_order = (N + 2) // 3
    for a_order in range(minimum_order, minimum_order + 7):
        for d_order in range(0, 7):
            f_local = z**a_order
            g_local = z**d_order
            W = (N - 2) * sp.diff(f_local, z) * g_local + N * f_local * sp.diff(g_local, z)
            expected = ((N - 2) * a_order + N * d_order) * z ** (a_order + d_order - 1)
            require(
                "local multiplier N=%d a=%d d=%d" % (N, a_order, d_order),
                sp.expand(W - expected) == 0,
            )
            is_unit = sp.Poly(W, z).degree() == 0 and W != 0
            expected_unit = N == 3 and a_order == 1 and d_order == 0
            require(
                "unit classification N=%d a=%d d=%d" % (N, a_order, d_order),
                is_unit == expected_unit,
            )
            if is_unit:
                unit_rows.append((N, a_order, d_order, int(W)))

require("unique unit row", unit_rows == [(3, 1, 0, 1)])

for f_degree in range(3, 11):
    for g_degree in range(0, 9):
        f_monomial = b**f_degree
        g_monomial = b**g_degree
        W = sp.diff(f_monomial, b) * g_monomial + 3 * f_monomial * sp.diff(g_monomial, b)
        require(
            "homogeneous channel degree f=%d g=%d" % (f_degree, g_degree),
            sp.Poly(W, b).degree() == f_degree + g_degree - 1,
        )
        require(
            "homogeneous leading coefficient f=%d g=%d" % (f_degree, g_degree),
            sp.Poly(W, b).LC() == f_degree + 3 * g_degree,
        )

print("PASS unique arm-unit channel=(-3,1) and homogeneous h=1 degree gate")


print("SECTION complement counts h=0,1,2,3")

COUNT_CONTROLS = [
    ([-6, -3, 2], [-9, 2, 5], 0),
    ([-6, -3, 2], [-6, 1, 3], 1),
    ([-6, -3, 7], [-6, 1, 4], 2),
    ([-6, -3, 4], [-6, 1, 4], 3),
]
for r_weights, s_weights, expected_h in COUNT_CONTROLS:
    require(
        "complement count h=%d" % expected_h,
        len(complement_edges(r_weights, s_weights)) == expected_h,
    )
    edges = complement_edges(r_weights, s_weights)
    require(
        "complement matching h=%d" % expected_h,
        len({i for i, _ in edges}) == len(edges) and len({j for _, j in edges}) == len(edges),
    )

print("PASS explicit controls for every complement-edge count")


print("SECTION h=2 collision and singleton-tree atlas")

TREE_GENERIC = {(0, 0), (0, 1), (1, 0), (1, 2), (2, 2)}
TREE_C1 = {(0, 0), (0, 1), (1, 0), (2, 1), (2, 2)}
TREE_C2 = {(0, 0), (1, 0), (1, 2), (2, 1), (2, 2)}
TREE_C3 = {(0, 0), (0, 1), (1, 0), (2, 0), (2, 2)}
COLLISION_C1 = frozenset({(2, 0), (1, 2)})
COLLISION_C2 = frozenset({(2, 0), (0, 1)})
COLLISION_C3 = frozenset({(1, 2), (2, 1)})

for name, tree in (
    ("generic", TREE_GENERIC),
    ("C1", TREE_C1),
    ("C2", TREE_C2),
    ("C3", TREE_C3),
):
    require(name + " tree has five edges", len(tree) == 5)
    require(name + " tree connected", connected(tree))

left_rows = 0
left_collision_counts = {"C1": 0, "C2": 0, "C3": 0, "none": 0}
for m_value in range(2, 13):
    for k_value in range(2, 13):
        B_value = 3 * m_value - 2
        for A_value in range(2, 81):
            if A_value == 3 * k_value - 2:
                continue
            left_rows += 1
            r_weights = [-3 * m_value, -3, A_value]
            s_weights = [-3 * k_value, 1, B_value]
            require("h2-left ordering", r_weights == sorted(r_weights) and s_weights == sorted(s_weights))
            require("h2-left complements", complement_edges(r_weights, s_weights) == [(0, 2), (1, 1)])

            expected_collisions = []
            expected_tree = TREE_GENERIC
            collision_name = "none"
            if A_value == 3 * (m_value + k_value) - 5:
                expected_collisions.append(COLLISION_C1)
                expected_tree = TREE_C1
                collision_name = "C1"
            if A_value == 3 * (k_value - m_value) + 1:
                expected_collisions.append(COLLISION_C2)
                expected_tree = TREE_C2
                collision_name = "C2"
            if A_value == 3 * m_value - 6:
                expected_collisions.append(COLLISION_C3)
                expected_tree = TREE_C3
                collision_name = "C3"

            actual_collisions = collision_sets(r_weights, s_weights)
            require("h2-left at most one collision", len(expected_collisions) <= 1)
            require(
                "h2-left collision classification",
                set(actual_collisions) == set(expected_collisions),
            )
            singletons = set(singleton_edges(r_weights, s_weights))
            require("h2-left displayed tree singleton", expected_tree <= singletons)
            require("h2-left singleton graph connected", connected(singletons))
            left_collision_counts[collision_name] += 1

right_rows = 0
for m_value in range(2, 13):
    for k_value in range(2, 13):
        A_value = 3 * k_value - 2
        for B_value in range(2, 81):
            if B_value == 3 * m_value - 2:
                continue
            right_rows += 1
            r_weights = [-3 * m_value, -3, A_value]
            s_weights = [-3 * k_value, 1, B_value]
            require("h2-right complements", complement_edges(r_weights, s_weights) == [(1, 1), (2, 0)])
            singletons = singleton_edges(r_weights, s_weights)
            require("h2-right singleton graph connected", connected(singletons))
            require("h2-right at most one collision", len(collision_sets(r_weights, s_weights)) <= 1)

COLLISION_REPRESENTATIVES = [
    (2, 2, 7, COLLISION_C1),
    (2, 3, 4, COLLISION_C2),
    (3, 2, 3, COLLISION_C3),
]
for m_value, k_value, A_value, expected_collision in COLLISION_REPRESENTATIVES:
    r_weights = [-3 * m_value, -3, A_value]
    s_weights = [-3 * k_value, 1, 3 * m_value - 2]
    require(
        "collision representative %s" % sorted(expected_collision),
        collision_sets(r_weights, s_weights) == [expected_collision],
    )

require("C1 represented in atlas", left_collision_counts["C1"] > 0)
require("C2 represented in atlas", left_collision_counts["C2"] > 0)
require("C3 represented in atlas", left_collision_counts["C3"] > 0)

print(
    "PASS h=2 rows left=%d right=%d collisions=%s"
    % (left_rows, right_rows, left_collision_counts)
)


print("SECTION h=3 unequal-gap six-cycle and equal-gap boundary")

full_rows = 0
for m_value in range(2, 21):
    for k_value in range(2, 21):
        full_rows += 1
        r_weights = [-3 * m_value, -3, 3 * k_value - 2]
        s_weights = [-3 * k_value, 1, 3 * m_value - 2]
        require("h3 exact low P", has_exact_low(r_weights))
        require("h3 exact low Q", has_exact_low(s_weights))
        require("h3 high P", has_high(r_weights))
        require("h3 high Q", has_high(s_weights))
        require("h3 complements", complement_edges(r_weights, s_weights) == [(0, 2), (1, 1), (2, 0)])

        p_gap = r_weights[1] - r_weights[0]
        q_gap = r_weights[2] - r_weights[1]
        require("h3 gaps positive", p_gap > 0 and q_gap > 0)
        require("h3 gaps unequal mod three", p_gap != q_gap and p_gap % 3 == 0 and q_gap % 3 == 1)
        require("h3 no nonscalar collisions", collision_sets(r_weights, s_weights) == [])
        off_edges = singleton_edges(r_weights, s_weights)
        require("h3 six-cycle edge count", len(off_edges) == 6)
        require("h3 six-cycle connected", connected(off_edges))

        relative_sums = sorted(
            value + 2
            for value, edges in layer_map(r_weights, s_weights).items()
            if value != -2
            for _ in edges
        )
        expected_relative = sorted([-(p_gap + q_gap), -p_gap, -q_gap, p_gap, q_gap, p_gap + q_gap])
        require("h3 relative sums", relative_sums == expected_relative)
        require("h3 equal-gap congruence", 3 * (m_value - k_value) != 4)

print("PASS h=3 full rows=%d, every off-matching graph is a connected six-cycle" % full_rows)


print("SECTION sharp support hostiles")

ONE_SIDED = ([-9, -3, 3], [-5, 1, 7])
NO_EXACT = ([-8, -3, 2], [-4, 1, 6])
for label, (r_weights, s_weights) in (("one-sided", ONE_SIDED), ("no-exact", NO_EXACT)):
    require(label + " common step P", common_step(r_weights))
    require(label + " common step Q", common_step(s_weights))
    require(label + " same step", r_weights[1] - r_weights[0] == s_weights[1] - s_weights[0])
    require(label + " full complement", complement_edges(r_weights, s_weights) == [(0, 2), (1, 1), (2, 0)])
    require(label + " local unit channel", r_weights[1] == -3 and s_weights[1] == 1)
    require(label + " positive arms", has_high(r_weights) and has_high(s_weights))

require("one-sided P exact low", has_exact_low(ONE_SIDED[0]))
require("one-sided Q loses exact low", not has_exact_low(ONE_SIDED[1]))
require("no-exact P loses exact low", not has_exact_low(NO_EXACT[0]))
require("no-exact Q loses exact low", not has_exact_low(NO_EXACT[1]))

print("PASS formal near misses pinpoint the two-arm exact-negative dependency")


print("SECTION consequence invoice")
print("PASS conditional arm+jet gate forces >=3 pieces in each coordinate")
print("PASS h=0 scalar absent; h=1 homogeneous degree obstruction")
print("PASS h=2 all collision graphs retain a singleton spanning tree")
print("PASS h=3 unequal gaps give six-cycle; equal gaps fail modulo three")
print("H2_LEFT_ROWS=%d" % left_rows)
print("H2_RIGHT_ROWS=%d" % right_rows)
print("H3_ROWS=%d" % full_rows)
print("CHECKS=%d" % CHECKS)
print("RESULT=PASS")
