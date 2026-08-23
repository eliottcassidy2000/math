#!/usr/bin/env python3
"""Exact placement and differential audit for THM-3712."""

from collections import defaultdict
from math import gcd

import sympy as sp


def require(condition, message):
    if not condition:
        raise AssertionError(message)


# ---------------------------------------------------------------------------
# Actual-support census after the inherited scalar/singleton and parity gates.
# ---------------------------------------------------------------------------

A_POS = (0, 2, 3)
B_POS = (0, 1, 2, 3)
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
        ((0, 1),),
        ((0, 2), (1, 0)),
        ((0, 3), (1, 1), (2, 0)),
        ((1, 2), (2, 1)),
        ((1, 3), (2, 2)),
        ((2, 3),),
    ),
    "W007 fibre word",
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


def family_a(n):
    return (
        4,
        (-2 * n - 2, -2, n - 2),
        (1 - 2 * n, 1 - n, 1, n + 1),
    )


def family_b(n):
    return (
        4,
        (1 - 2 * n, 1, n + 1),
        (-2 * n - 2, -n - 2, -2, n - 2),
    )


def family_c(n):
    return (
        5,
        (1 - 3 * n, 1 - n, 1),
        (-2 * n - 2, -n - 2, -2, n - 2),
    )


def family_d(n):
    return (
        5,
        (-2 * n - 2, -2, n - 2),
        (1 - 3 * n, 1 - 2 * n, 1 - n, 1),
    )


candidate_counts = {}
survivor_counts = {}
for scale in range(1, 13):
    inherited = inherited_placements(scale)
    survivors = {
        placement for placement in inherited if not parity_rejected(placement)
    }
    if scale <= 2:
        expected = set()
    elif scale == 3:
        # C=D as an actual support pair; two arm anchors are not two systems.
        expected = {family_a(scale), family_b(scale), family_c(scale)}
        require(family_c(scale) == family_d(scale), "n=3 C/D merger")
    else:
        expected = {
            family_a(scale), family_b(scale),
            family_c(scale), family_d(scale),
        }
    require(survivors == expected, f"post-parity actual-support census n={scale}")
    candidate_counts[scale] = len(inherited)
    survivor_counts[scale] = len(survivors)

require(
    tuple(candidate_counts[scale] for scale in (1, 2, 3, 4, 5))
    == (2, 4, 6, 8, 8),
    "inherited count boundary and tail",
)
require(
    tuple(survivor_counts[scale] for scale in (1, 2, 3, 4, 5))
    == (0, 0, 3, 4, 4),
    "post-parity count boundary and tail",
)


# ---------------------------------------------------------------------------
# Exact all-scale differential integrations for the four residual families.
# ---------------------------------------------------------------------------

n = sp.symbols("n", integer=True, positive=True)
H, K = sp.symbols("H K", nonzero=True)
Hp, Kp = sp.symbols("Hp Kp")
a, b0, c0, d0, t0 = sp.symbols("a b0 c0 d0 t0", nonzero=True)
kappa, lam = sp.symbols("kappa lam")


def derivative(expression):
    return sp.diff(expression, H) * Hp + sp.diff(expression, K) * Kp


def wedge(r, left, s, right):
    return sp.expand(s * derivative(left) * right - r * left * derivative(right))


def zero(expression, label):
    if sp.simplify(sp.powsimp(sp.factor(expression), force=True)) != 0:
        raise AssertionError(label)


# Family A: scalar 12+21 with persistent arm 12=(-2,1).
a_f0 = a * H ** (2 * n + 2)
a_g0 = b0 * H ** (2 * n - 1)
a_g2 = d0 * K
a_gamma = 2 * (n + 1) * a * d0 / ((2 * n - 1) * b0)
a_f1 = H**2 * (kappa + a_gamma * H * K)
a_low = wedge(-2 * n - 2, a_f0, 1, a_g2)
a_low += wedge(-2, a_f1, 1 - 2 * n, a_g0)
zero(a_low, "family A arm transport")

# Family B: scalar 12+21 with persistent arm 12=(1,-2).
b_f0 = a * H ** (2 * n - 1)
b_g0 = b0 * H ** (2 * n + 2)
b_f1 = d0 * K
b_gamma = 2 * (n + 1) * b0 * d0 / (a * (2 * n - 1))
b_g2 = H**2 * (kappa + b_gamma * H * K)
b_low = wedge(1 - 2 * n, b_f0, -2, b_g2)
b_low += wedge(1, b_f1, -2 * n - 2, b_g0)
zero(b_low, "family B arm transport")

# Family C: scalar 13+22 with persistent arm 22=(1,-2).  The triple first
# makes f1/H^(n-1) polynomial; the lower double then transfers it into g2/H^2.
c_f0 = a * H ** (3 * n - 1)
c_g0 = b0 * H ** (2 * n + 2)
c_g1 = c0 * H ** (n + 2)
c_f2 = d0 * K
c_g3 = t0 * K ** (n - 2)
c_alpha = a * t0 * (3 * n - 1) / (c0 * (n + 2))
c_beta = 2 * (n + 1) * b0 * d0 / (c0 * (n + 2))
c_base = kappa + c_alpha * (H * K) ** (n - 2) - c_beta * H * K
c_f1 = H ** (n - 1) * c_base
c_theta = 2 * (n + 1) * b0 / (a * (3 * n - 1))
c_g2 = H**2 * (lam + c_theta * c_base)
c_triple = wedge(1 - 3 * n, c_f0, n - 2, c_g3)
c_triple += wedge(1 - n, c_f1, -n - 2, c_g1)
c_triple += wedge(1, c_f2, -2 * n - 2, c_g0)
c_low = wedge(1 - 3 * n, c_f0, -2, c_g2)
c_low += wedge(1 - n, c_f1, -2 * n - 2, c_g0)
zero(c_triple, "family C triple integration")
zero(c_low, "family C arm transport")

# Family D: scalar 13+22 with persistent arm 13=(-2,1).  The triple itself
# integrates f1/H^2.  At n=3, C=D and C also supplies the second arm square.
d_f0 = a * H ** (2 * n + 2)
d_g0 = b0 * H ** (3 * n - 1)
d_g1 = c0 * H ** (2 * n - 1)
d_f2 = d0 * K ** (n - 2)
d_g3 = t0 * K
d_alpha = 2 * (n + 1) * a * t0 / ((2 * n - 1) * c0)
d_beta = (3 * n - 1) * b0 * d0 / ((2 * n - 1) * c0)
d_f1 = H**2 * (
    kappa + d_alpha * H * K - d_beta * (H * K) ** (n - 2)
)
d_triple = wedge(-2 * n - 2, d_f0, 1, d_g3)
d_triple += wedge(-2, d_f1, 1 - 2 * n, d_g1)
d_triple += wedge(n - 2, d_f2, 1 - 3 * n, d_g0)
zero(d_triple, "family D triple arm transport")


print("THM-3712 exact W007 completion audit")
print("fibre word = 00|01|02+10|03+11+20|12+21|13+22|23")
print("inherited placement counts n=1,2,3,>=4 = 2,4,6,8")
print("post-parity actual placements n<=2,n=3,n>=4 = 0,3,4")
print("n=3 actual-support merger = C=D with both arm addresses retained")
print("families A,B = one-row H^2 arm transport")
print("family C = triple polynomiality then lower-row H^2 transport")
print("family D = triple-row H^2 transport")
print("scope = complete W007 ray in the y=0 collision ring")
print("ALL CHECKS PASSED")
