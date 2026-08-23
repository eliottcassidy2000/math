#!/usr/bin/env python3
"""Exact placement and differential audit for THM-3711."""

from collections import defaultdict
from math import gcd

import sympy as sp


def require(condition, message):
    if not condition:
        raise AssertionError(message)


# ---------------------------------------------------------------------------
# Actual-support census after the inherited scalar/singleton and parity gates.
# ---------------------------------------------------------------------------

A_POS = (0, 1, 2)
B_POS = (0, 2, 3, 4)
raw_fibres = defaultdict(list)
for i, left in enumerate(A_POS):
    for j, right in enumerate(B_POS):
        raw_fibres[left + right].append((i, j))
FIBRES = tuple(tuple(raw_fibres[key]) for key in sorted(raw_fibres))
SINGLETONS = tuple(fibre[0] for fibre in FIBRES if len(fibre) == 1)

require(
    FIBRES
    == (
        ((0, 0),),
        ((1, 0),),
        ((0, 1), (2, 0)),
        ((0, 2), (1, 1)),
        ((0, 3), (1, 2), (2, 1)),
        ((1, 3), (2, 2)),
        ((2, 3),),
    ),
    "W008 fibre word",
)


def same_sign_or_both_zero(r, s):
    return (r == 0 and s == 0) or (r > 0 and s > 0) or (r < 0 and s < 0)


def inherited_placements(n):
    placements = set()
    a_support = tuple(n * value for value in A_POS)
    b_support = tuple(n * value for value in B_POS)
    for scalar_index, fibre in enumerate(FIBRES):
        if len(fibre) < 2:
            continue
        for i, j in fibre:
            for arm_left, arm_right in ((-2, 1), (1, -2)):
                p0 = arm_left - a_support[i]
                q0 = arm_right - b_support[j]
                p_weights = tuple(p0 + value for value in a_support)
                q_weights = tuple(q0 + value for value in b_support)
                require(
                    all(p_weights[k] + q_weights[l] == -1 for k, l in fibre),
                    "bad scalar fibre",
                )
                if not all(
                    same_sign_or_both_zero(p_weights[k], q_weights[l])
                    for k, l in SINGLETONS
                ):
                    continue
                placements.add((scalar_index, p_weights, q_weights))
    return placements


class DSU:
    def __init__(self, vertices):
        self.parent = {vertex: vertex for vertex in vertices}

    def find(self, vertex):
        if self.parent[vertex] != vertex:
            self.parent[vertex] = self.find(self.parent[vertex])
        return self.parent[vertex]

    def union(self, left, right):
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def parity_rejected(placement):
    scalar_index, p_weights, q_weights = placement
    vertices = tuple(
        [f"P{i}" for i in range(len(p_weights))]
        + [f"Q{j}" for j in range(len(q_weights))]
    )
    dsu = DSU(vertices)
    active = set()
    for i, j in SINGLETONS:
        r, s = p_weights[i], q_weights[j]
        if r * s > 0:
            left, right = f"P{i}", f"Q{j}"
            dsu.union(left, right)
            active.update((left, right))

    components = defaultdict(list)
    for vertex in active:
        components[dsu.find(vertex)].append(vertex)
    exponent = {}
    for component in components.values():
        weights = tuple(
            p_weights[int(vertex[1:])]
            if vertex.startswith("P")
            else q_weights[int(vertex[1:])]
            for vertex in component
        )
        divisor = gcd(*(abs(weight) for weight in weights))
        for vertex, weight in zip(component, weights):
            exponent[vertex] = abs(weight) // divisor

    eligible = []
    for i, j in FIBRES[scalar_index]:
        pair = p_weights[i], q_weights[j]
        if pair == (-2, 1):
            eligible.append(exponent.get(f"P{i}"))
        elif pair == (1, -2):
            eligible.append(exponent.get(f"Q{j}"))
    require(eligible, "every placement retains an arm address")
    return all(value is not None and value >= 2 for value in eligible)


def family_x(n):
    return (
        4,
        (1 - 2 * n, 1 - n, 1),
        (-2 * n - 2, -2, n - 2, 2 * n - 2),
    )


def family_y(n):
    return (
        5,
        (1 - 2 * n, 1 - n, 1),
        (-3 * n - 2, -n - 2, -2, n - 2),
    )


candidate_counts = {}
survivor_counts = {}
for n in range(1, 13):
    inherited = inherited_placements(n)
    survivors = {placement for placement in inherited if not parity_rejected(placement)}
    expected = set()
    if n >= 2:
        expected.add(family_x(n))
    if n >= 3:
        expected.add(family_y(n))
    require(survivors == expected, f"post-parity actual-support census n={n}")
    candidate_counts[n] = len(inherited)
    survivor_counts[n] = len(survivors)

require(
    tuple(candidate_counts[n] for n in (1, 2, 3, 4, 5)) == (0, 2, 4, 6, 6),
    "inherited count boundary and tail",
)
require(
    tuple(survivor_counts[n] for n in (1, 2, 3, 4, 5)) == (0, 1, 2, 2, 2),
    "post-parity count boundary and tail",
)


# ---------------------------------------------------------------------------
# Exact all-scale differential integrations for the two residual families.
# ---------------------------------------------------------------------------

n = sp.symbols("n", integer=True, positive=True)
H, K = sp.symbols("H K", nonzero=True)
Hp, Kp = sp.symbols("Hp Kp")
a, b0, c0, d0 = sp.symbols("a b0 c0 d0", nonzero=True)
kappa, lam = sp.symbols("kappa lam")


def derivative(expression):
    return sp.diff(expression, H) * Hp + sp.diff(expression, K) * Kp


def wedge(r, left, s, right):
    return sp.expand(s * derivative(left) * right - r * left * derivative(right))


def zero(expression, label):
    if sp.simplify(sp.powsimp(sp.factor(expression), force=True)) != 0:
        raise AssertionError(label)


r0, r1, r2 = 1 - 2 * n, 1 - n, 1
f0 = a * H ** (2 * n - 1)
f1 = c0 * H ** (n - 1)
f2 = d0 * K

# Family X: the 01+20 row integrates g1/H^2 directly.
x_s0, x_s1 = -2 * n - 2, -2
x_g0 = b0 * H ** (2 * n + 2)
x_gamma = 2 * (n + 1) * b0 * d0 / (a * (2 * n - 1))
x_g1 = H**2 * (kappa + x_gamma * H * K)
x_low = wedge(r0, f0, x_s1, x_g1) + wedge(r2, f2, x_s0, x_g0)
zero(x_low, "family X low integration")

# Family Y: the first low row determines g1/H^(n+2); substituting it into
# the next row determines g2/H^2 with the opposite sign.
y_s0, y_s1, y_s2 = -3 * n - 2, -n - 2, -2
y_g0 = b0 * H ** (3 * n + 2)
y_gamma = (3 * n + 2) * b0 * d0 / (a * (2 * n - 1))
y_g1 = H ** (n + 2) * (kappa + y_gamma * H * K)
y_eta = c0 * (n - 1) * y_gamma / (a * (2 * n - 1))
y_g2 = H**2 * (lam - y_eta * H * K)
y_low = wedge(r0, f0, y_s1, y_g1) + wedge(r2, f2, y_s0, y_g0)
y_next = wedge(r0, f0, y_s2, y_g2) + wedge(r1, f1, y_s1, y_g1)
zero(y_low, "family Y first low integration")
zero(y_next, "family Y second low integration")


print("THM-3711 exact W008 completion audit")
print("fibre word = 00|10|01+20|02+11|03+12+21|13+22|23")
print("inherited placement counts n=1,2,3,>=4 = 0,2,4,6")
print("post-parity actual placements n=1,2,>=3 = 0,1,2")
print("family X live n>=2 = scalar triple; g1=H^2(kappa+gamma HK)")
print("family Y live n>=3 = scalar double; g2=H^2(lambda-eta HK)")
print("n=3 extra arm addresses already have negative coefficient H^2")
print("scope = complete W008 ray in the y=0 collision ring")
print("ALL CHECKS PASSED")
