#!/usr/bin/env python3
"""Exact odd-scale W003 lowest-scalar audit for THM-3715."""

from collections import defaultdict
from math import gcd

import sympy as sp


def require(condition, message):
    if not condition:
        raise AssertionError(message)


# ---------------------------------------------------------------------------
# Exact post-parity actual-support census on scalar fibre 01+10.
# ---------------------------------------------------------------------------

A_POS = (0, 1, 3)
B_POS = (0, 1, 2, 3)
raw_fibres = defaultdict(list)
for i, left in enumerate(A_POS):
    for j, right in enumerate(B_POS):
        raw_fibres[left + right].append((i, j))
FIBRES = tuple(tuple(raw_fibres[key]) for key in sorted(raw_fibres))
SINGLETONS = tuple(fibre[0] for fibre in FIBRES if len(fibre) == 1)

require(FIBRES[1] == ((0, 1), (1, 0)), "W003 lowest double")


def same_sign_or_both_zero(r, s):
    return (r == 0 and s == 0) or (r > 0 and s > 0) or (r < 0 and s < 0)


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


def inherited_lowest(n):
    scalar_index = 1
    fibre = FIBRES[scalar_index]
    a_support = tuple(n * value for value in A_POS)
    b_support = tuple(n * value for value in B_POS)
    placements = set()
    for i, j in fibre:
        for arm_left, arm_right in ((-2, 1), (1, -2)):
            p0 = arm_left - a_support[i]
            q0 = arm_right - b_support[j]
            p_weights = tuple(p0 + value for value in a_support)
            q_weights = tuple(q0 + value for value in b_support)
            if not all(
                same_sign_or_both_zero(p_weights[k], q_weights[l])
                for k, l in SINGLETONS
            ):
                continue
            placements.add((scalar_index, p_weights, q_weights))
    return placements


def parity_rejected(placement):
    scalar_index, p_weights, q_weights = placement
    vertices = tuple(
        [f"P{i}" for i in range(3)] + [f"Q{j}" for j in range(4)]
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
        if (p_weights[i], q_weights[j]) == (-2, 1):
            eligible.append(exponent.get(f"P{i}"))
        elif (p_weights[i], q_weights[j]) == (1, -2):
            eligible.append(exponent.get(f"Q{j}"))
    require(eligible, "every placement retains an arm address")
    return all(value is not None and value >= 2 for value in eligible)


def family_o1(n):
    return (
        1,
        (1 - n, 1, 2 * n + 1),
        (-2, n - 2, 2 * n - 2, 3 * n - 2),
    )


def family_o2(n):
    return (
        1,
        (-2, n - 2, 3 * n - 2),
        (1 - n, 1, n + 1, 2 * n + 1),
    )


for scale in range(1, 14):
    survivors = {
        placement
        for placement in inherited_lowest(scale)
        if not parity_rejected(placement)
    }
    if scale < 3 or scale % 2 == 0:
        expected = set()
    elif scale == 3:
        require(family_o1(scale) == family_o2(scale), "n=3 O1/O2 merger")
        expected = {family_o1(scale)}
    else:
        expected = {family_o1(scale), family_o2(scale)}
    require(survivors == expected, f"odd lowest-scalar census n={scale}")


# ---------------------------------------------------------------------------
# Exact all-odd-scale half-charge identities.
# ---------------------------------------------------------------------------

m = sp.symbols("m", integer=True, positive=True)
n = 2 * m + 1
H, K = sp.symbols("H K", nonzero=True)
Hp, Kp = sp.symbols("Hp Kp")
a, b0, d0, e0, t0 = sp.symbols("a b0 d0 e0 t0", nonzero=True)
kappa, lam = sp.symbols("kappa lam")


def derivative(expression):
    return sp.diff(expression, H) * Hp + sp.diff(expression, K) * Kp


def wedge(r, left, s, right):
    return sp.expand(s * derivative(left) * right - r * left * derivative(right))


def zero(expression, label):
    if sp.simplify(sp.powsimp(sp.factor(expression), force=True)) != 0:
        raise AssertionError(label)


euler = Hp * K + 2 * H * Kp


# Orientation O1: arm 10=(1,-2).
o1_f0 = a * H**m
o1_g0 = b0 * H
o1_f2 = d0 * K ** (2 * n + 1)
o1_g2 = e0 * K ** (2 * n - 2)
o1_g3 = t0 * K ** (3 * n - 2)
o1_alpha = (3 * n - 2) * a * t0 / (2 * (n - 1) * e0)
o1_beta = (2 * n + 1) * b0 * d0 / (2 * (n - 1) * e0)
o1_f1_terms = (
    (kappa, 0, 1),
    (-o1_alpha, m, n),
    (o1_beta, 1, 3),
)
o1_f1 = sum(coefficient * H**p * K**q for coefficient, p, q in o1_f1_terms)
o1_rho = (3 * n - 2) * t0 / ((2 * n + 1) * d0)
o1_g1 = lam * K ** (n - 2) + o1_rho * o1_f1 * K ** (n - 3)
o1_top = wedge(1, o1_f1, 3 * n - 2, o1_g3)
o1_top += wedge(2 * n + 1, o1_f2, n - 2, o1_g1)
o1_triple = wedge(1 - n, o1_f0, 3 * n - 2, o1_g3)
o1_triple += wedge(1, o1_f1, 2 * n - 2, o1_g2)
o1_triple += wedge(2 * n + 1, o1_f2, -2, o1_g0)
o1_scalar = wedge(1 - n, o1_f0, n - 2, o1_g1)
o1_scalar += wedge(1, o1_f1, -2, o1_g0)

o1_g1_terms = [(lam, 0, n - 2)]
o1_g1_terms += [
    (o1_rho * coefficient, p, q + n - 3)
    for coefficient, p, q in o1_f1_terms
]
o1_scalar_quotient = sum(
    a * coefficient * m * q * H ** (m + p - 1) * K ** (q - 1)
    for coefficient, p, q in o1_g1_terms
)
o1_scalar_quotient += sum(
    -b0 * coefficient * (2 * p + 1) * H**p * K ** (q - 1)
    for coefficient, p, q in o1_f1_terms
)
require(
    all(sp.expand(q - (2 * p + 1)) == 0 for _coefficient, p, q in o1_f1_terms),
    "O1 f1 half-charge",
)
require(
    all(
        sp.expand(q - (n - 2 + 2 * p)) == 0
        for _coefficient, p, q in o1_g1_terms
    ),
    "O1 g1 half-charge",
)
zero(o1_top, "O1 top")
zero(o1_triple, "O1 triple")
zero(o1_scalar - euler * o1_scalar_quotient, "O1 scalar Euler")


# Orientation O2: arm 01=(-2,1).
o2_f0 = a * H
o2_g0 = b0 * H**m
o2_f2 = d0 * K ** (3 * n - 2)
o2_g2 = e0 * K ** (n + 1)
o2_g3 = t0 * K ** (2 * n + 1)
o2_alpha = (2 * n + 1) * a * t0 / ((n + 1) * e0)
o2_beta = (3 * n - 2) * b0 * d0 / ((n + 1) * e0)
o2_u_terms = (
    (kappa, 0, 0),
    (-o2_alpha, 1, 2),
    (o2_beta, m, n - 1),
)
o2_u = sum(coefficient * H**p * K**q for coefficient, p, q in o2_u_terms)
o2_f1_terms = tuple(
    (coefficient, p, q + n - 2) for coefficient, p, q in o2_u_terms
)
o2_f1 = K ** (n - 2) * o2_u
o2_rho = (2 * n + 1) * t0 / ((3 * n - 2) * d0)
o2_g1_terms = ((lam, 0, 1),) + tuple(
    (o2_rho * coefficient, p, q + 1)
    for coefficient, p, q in o2_u_terms
)
o2_g1 = K * (lam + o2_rho * o2_u)
o2_top = wedge(n - 2, o2_f1, 2 * n + 1, o2_g3)
o2_top += wedge(3 * n - 2, o2_f2, 1, o2_g1)
o2_triple = wedge(-2, o2_f0, 2 * n + 1, o2_g3)
o2_triple += wedge(n - 2, o2_f1, n + 1, o2_g2)
o2_triple += wedge(3 * n - 2, o2_f2, 1 - n, o2_g0)
o2_scalar = wedge(-2, o2_f0, 1, o2_g1)
o2_scalar += wedge(n - 2, o2_f1, 1 - n, o2_g0)
o2_scalar_quotient = sum(
    a * coefficient * (2 * p + 1) * H**p * K ** (q - 1)
    for coefficient, p, q in o2_g1_terms
)
o2_scalar_quotient += sum(
    -b0 * coefficient * m * q * H ** (m + p - 1) * K ** (q - 1)
    for coefficient, p, q in o2_f1_terms
)
require(
    all(
        sp.expand(q - (n - 2 + 2 * p)) == 0
        for _coefficient, p, q in o2_f1_terms
    ),
    "O2 f1 half-charge",
)
require(
    all(sp.expand(q - (2 * p + 1)) == 0 for _coefficient, p, q in o2_g1_terms),
    "O2 g1 half-charge",
)
zero(o2_top, "O2 top")
zero(o2_triple, "O2 triple")
zero(o2_scalar - euler * o2_scalar_quotient, "O2 scalar Euler")


print("THM-3715 exact odd W003 lowest-scalar audit")
print("post-parity actual placements = none for even n or n=1")
print("odd actual placements n=3,n>=5 = 1,2")
print("n=3 merger = O1=O2 with both arm addresses retained")
print("O1 and O2 upper-row integrations PASS")
print("half-charge laws in both scalar addresses PASS")
print("common scalar factor = H'K+2HK'")
print("scope = complete W003 scalar fibre 01+10")
print("ALL CHECKS PASSED")
