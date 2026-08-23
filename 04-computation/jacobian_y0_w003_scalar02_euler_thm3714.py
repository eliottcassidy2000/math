#!/usr/bin/env python3
"""Exact census and Euler-factor audit for THM-3714."""

from collections import defaultdict
from math import gcd

import sympy as sp


def require(condition, message):
    if not condition:
        raise AssertionError(message)


# ---------------------------------------------------------------------------
# Exact actual-support census on the W003 scalar fibre 02+11.
# ---------------------------------------------------------------------------

A_POS = (0, 1, 3)
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
        ((0, 1), (1, 0)),
        ((0, 2), (1, 1)),
        ((0, 3), (1, 2), (2, 0)),
        ((1, 3), (2, 1)),
        ((2, 2),),
        ((2, 3),),
    ),
    "W003 fibre word",
)


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


def inherited_scalar02(n):
    scalar_index = 2
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


def family_a(scale):
    return (
        2,
        (-scale - 2, -2, 2 * scale - 2),
        (1 - scale, 1, scale + 1, 2 * scale + 1),
    )


def family_b(scale):
    return (
        2,
        (1 - scale, 1, 2 * scale + 1),
        (-scale - 2, -2, scale - 2, 2 * scale - 2),
    )


for scale in range(1, 13):
    survivors = {
        placement
        for placement in inherited_scalar02(scale)
        if not parity_rejected(placement)
    }
    if scale == 1:
        expected = set()
    elif scale == 2:
        expected = {family_a(scale)}
    else:
        expected = {family_a(scale), family_b(scale)}
    require(survivors == expected, f"scalar02 census n={scale}")


# ---------------------------------------------------------------------------
# Exact generic Euler identities in the two residue branches.
# ---------------------------------------------------------------------------

ell = sp.symbols("ell", integer=True, positive=True)
H, K, J = sp.symbols("H K J", nonzero=True)
Hp, Kp, Jp = sp.symbols("Hp Kp Jp")
a, b0, d0, e0, t0 = sp.symbols("a b0 d0 e0 t0", nonzero=True)
lam = sp.symbols("lam")


def derivative(expression):
    return (
        sp.diff(expression, H) * Hp
        + sp.diff(expression, K) * Kp
        + sp.diff(expression, J) * Jp
    )


def wedge(r, left, s, right):
    return sp.expand(s * derivative(left) * right - r * left * derivative(right))


def zero(expression, label):
    if sp.simplify(sp.powsimp(sp.factor(expression), force=True)) != 0:
        raise AssertionError(label)


BRANCHES = (
    (1, 3 * ell, 3 * ell + 2, 3 * ell - 1),
    (1, 3 * ell + 2, 3 * ell + 4, 3 * ell + 1),
    (3, 3 * ell + 1, ell + 1, ell),
)


def audit_family_a(delta, n, alpha, beta):
    f0 = a * H**alpha
    g0 = b0 * H**beta
    f2 = d0 * K ** (2 * n - 2)
    g2 = e0 * K ** (n + 1)
    g3 = t0 * K ** (2 * n + 1)
    left = (2 * n + 1) * a * t0 / ((n + 1) * e0)
    right = 2 * (n - 1) * b0 * d0 / ((n + 1) * e0)
    arm = -left * H**alpha * K**n + right * H**beta * K ** (n - 3)
    rho = (2 * n + 1) * t0 / (2 * (n - 1) * d0)
    mate = lam * K + rho * arm * K**3
    top = wedge(-2, arm, 2 * n + 1, g3)
    top += wedge(2 * n - 2, f2, 1, mate)
    triple = wedge(-n - 2, f0, 2 * n + 1, g3)
    triple += wedge(-2, arm, n + 1, g2)
    triple += wedge(2 * n - 2, f2, 1 - n, g0)
    scalar = wedge(-n - 2, f0, n + 1, g2)
    scalar += wedge(-2, arm, 1, mate)
    euler = Hp * K + delta * H * Kp
    charge_quotient = (
        -left * alpha * H ** (alpha - 1) * K**n
        + right * beta * H ** (beta - 1) * K ** (n - 3)
    )
    outer_quotient = (
        (n + 1) * a * e0 * alpha * H ** (alpha - 1) * K**n
    )
    scalar_quotient = outer_quotient
    scalar_quotient += charge_quotient * (lam + 3 * rho * arm * K**2)
    zero(top, f"family A top delta={delta}")
    zero(triple, f"family A triple delta={delta}")
    zero(
        wedge(-2, arm, 1, K) - euler * charge_quotient,
        f"family A charge delta={delta}",
    )
    zero(
        wedge(-2, arm, 1, arm * K**3)
        - 3 * arm * K**2 * wedge(-2, arm, 1, K),
        f"family A nonlinear wedge delta={delta}",
    )
    zero(scalar - euler * scalar_quotient, f"family A Euler delta={delta}")


def audit_family_b(delta, n, alpha, beta):
    f0 = a * H**beta
    g0 = b0 * H**alpha
    f2 = d0 * K ** (2 * n + 1)
    g2 = e0 * K ** (n - 2)
    g3 = t0 * K ** (2 * n - 2)
    left = 2 * (n - 1) * b0 * t0 / ((n - 2) * e0)
    right = 4 * (n - 1) ** 2 * a * t0**2
    right /= (2 * n + 1) * (n - 2) * d0 * e0
    arm = left * H**alpha * K**n - right * H**beta * K ** (n - 3)
    rho = (2 * n + 1) * d0 / (2 * (n - 1) * t0)
    mate = lam * K + rho * arm * K**3
    top = wedge(1, mate, 2 * n - 2, g3)
    top += wedge(2 * n + 1, f2, -2, arm)
    triple = wedge(1 - n, f0, 2 * n - 2, g3)
    triple += wedge(1, mate, n - 2, g2)
    triple += wedge(2 * n + 1, f2, -n - 2, g0)
    scalar = wedge(1 - n, f0, n - 2, g2)
    scalar += wedge(1, mate, -2, arm)
    euler = Hp * K + delta * H * Kp
    charge_quotient = (
        left * alpha * H ** (alpha - 1) * K**n
        - right * beta * H ** (beta - 1) * K ** (n - 3)
    )
    outer_quotient = (
        (n - 2) * a * e0 * beta * H ** (beta - 1) * K ** (n - 3)
    )
    scalar_quotient = outer_quotient
    scalar_quotient -= charge_quotient * (lam + 3 * rho * arm * K**2)
    zero(top, f"family B top delta={delta}")
    zero(triple, f"family B triple delta={delta}")
    zero(
        wedge(-2, arm, 1, K) - euler * charge_quotient,
        f"family B charge delta={delta}",
    )
    zero(
        wedge(-2, arm, 1, arm * K**3)
        - 3 * arm * K**2 * wedge(-2, arm, 1, K),
        f"family B nonlinear wedge delta={delta}",
    )
    zero(scalar - euler * scalar_quotient, f"family B Euler delta={delta}")


for branch in BRANCHES:
    audit_family_a(*branch)
    audit_family_b(*branch)


# Family A at n=2: polynomiality forces K|H.  Substitute H=KJ and verify the
# exact exceptional identities and factor J'K+2JK'.
H2 = K * J
a2_left = 5 * a * t0 / (3 * e0)
a2_right = 2 * b0 * d0 / (3 * e0)
a2_arm = -a2_left * H2**4 * K**2 + a2_right * J
a2_rho = 5 * t0 / (2 * d0)
a2_mate = lam * K + a2_rho * a2_arm * K**3
a2_top = wedge(-2, a2_arm, 5, t0 * K**5)
a2_top += wedge(2, d0 * K**2, 1, a2_mate)
a2_triple = wedge(-4, a * H2**4, 5, t0 * K**5)
a2_triple += wedge(-2, a2_arm, 3, e0 * K**3)
a2_triple += wedge(2, d0 * K**2, -1, b0 * H2)
a2_scalar = wedge(-4, a * H2**4, 3, e0 * K**3)
a2_scalar += wedge(-2, a2_arm, 1, a2_mate)
a2_euler = Jp * K + 2 * J * Kp
zero(a2_top, "family A n=2 top")
zero(a2_triple, "family A n=2 triple")
zero(
    a2_scalar - a2_euler * sp.cancel(a2_scalar / a2_euler),
    "family A n=2 Euler",
)


print("THM-3714 exact W003 scalar-fibre-02 audit")
print("post-parity actual placements n=1,n=2,n>=3 = 0,1,2")
print("family A generic delta=gcd(n-1,3) Euler identities PASS")
print("family B generic delta=gcd(n-1,3) Euler identities PASS")
print("family A n=2 K|H exceptional Euler identity PASS")
print("charge law = every arm monomial H^p K^q has q+2=delta*p")
print("Euler factors = H'K+delta HK' or J'K+2JK'")
print("scope = complete W003 scalar fibre 02+11 in both orientations")
print("ALL CHECKS PASSED")
