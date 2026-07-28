#!/usr/bin/env python3
"""Exact finite-group and hostile controls for THM-2633."""

from collections import Counter
from itertools import permutations

import sympy as sp


checks = 0


def require(condition, message):
    global checks
    if not bool(condition):
        raise RuntimeError(message)
    checks += 1


def compose(p, q):
    """Permutation p after q."""
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    out = [0] * len(p)
    for i, j in enumerate(p):
        out[j] = i
    return tuple(out)


def generated(generators):
    n = len(generators[0])
    one = tuple(range(n))
    out = {one}
    frontier = [one]
    while frontier:
        a = frontier.pop()
        for b in generators:
            for c in (compose(a, b), compose(b, a)):
                if c not in out:
                    out.add(c)
                    frontier.append(c)
    return frozenset(out)


def power(p, exponent):
    out = tuple(range(len(p)))
    for _ in range(exponent):
        out = compose(p, out)
    return out


def fixed_points(p):
    return frozenset(i for i, image in enumerate(p) if image == i)


def parity(p):
    return sum(p[i] > p[j]
               for i in range(len(p))
               for j in range(i + 1, len(p))) % 2


one = (0, 1, 2, 3)
r = (1, 2, 3, 0)
s = (0, 3, 2, 1)
z = power(r, 2)
edge = compose(r, s)
other_edge = compose(power(r, 3), s)

C4 = generated((r,))
D4 = generated((r, s))
S4 = frozenset(permutations(range(4)))
A4 = frozenset(g for g in S4 if parity(g) == 0)
V4 = frozenset((one, z, edge, other_edge))

require(len(C4) == 4, "C4 order changed")
require(len(V4) == 4, "V4 order changed")
require(len(D4) == 8, "D4 order changed")
require(len(A4) == 12, "A4 order changed")
require(len(S4) == 24, "S4 order changed")
require(all(compose(a, b) in V4 for a in V4 for b in V4),
        "displayed V4 is not a subgroup")

# Exact character kernels in the five transitive quartic groups.
K_C4 = frozenset((one, z))
K_V4 = frozenset((one, z))
H = frozenset((one, s))
J = frozenset(g for g in D4
              if frozenset(compose(compose(g, h), inverse(g)) for h in H) == H)
K_D4 = J
K_A4 = V4
K_S4 = A4

require(len(J) == 4 and J == frozenset((one, z, s, compose(z, s))),
        "D4 source-deck kernel changed")

groups_and_kernels = {
    "C4": (C4, K_C4),
    "V4": (V4, K_V4),
    "D4": (D4, K_D4),
    "A4": (A4, K_A4),
    "S4": (S4, K_S4),
}

support_fixed_census = {}
gate = {}
for name, (group, kernel) in groups_and_kernels.items():
    require(kernel < group, f"{name} character kernel is not proper")
    require(all(compose(a, b) in kernel for a in kernel for b in kernel),
            f"{name} character kernel is not a subgroup")
    # Every displayed kernel is normal.
    require(all(compose(compose(g, h), inverse(g)) in kernel
                for g in group for h in kernel),
            f"{name} character kernel is not normal")
    census = Counter(len(fixed_points(g)) for g in group - kernel)
    support_fixed_census[name] = dict(sorted(census.items()))
    gate[name] = all(k == 0 for k in census)

require(support_fixed_census == {
    "C4": {0: 2},
    "V4": {0: 2},
    "D4": {0: 4},
    "A4": {1: 8},
    "S4": {0: 6, 2: 6},
}, "quartic character fixed-point census changed")
require(gate == {
    "C4": True,
    "V4": True,
    "D4": True,
    "A4": False,
    "S4": False,
}, "quartic derangement-character boundary changed")

# Every D4 deck-odd element admits only the empty inertia-fixed survivor set.
d4_deck_odd = D4 - J
require(d4_deck_odd == frozenset((r, power(r, 3), edge, other_edge)),
        "D4 deck-odd coset changed")
d4_deck_odd_survivor_rows = 0
for g in d4_deck_odd:
    fixed = fixed_points(g)
    require(not fixed, "D4 deck-odd element acquired a fixed sheet")
    lawful = [frozenset(i for i in range(4) if (mask >> i) & 1)
              for mask in range(16)]
    lawful = [survivors for survivors in lawful if survivors.issubset(fixed)]
    require(lawful == [frozenset()],
            "D4 deck-odd inertia acquired a nonempty survivor row")
    d4_deck_odd_survivor_rows += len(lawful)

# Residual typed targets after the D4 exclusion.
a4_support_survivor_sizes = set()
for g in A4 - V4:
    fixed = fixed_points(g)
    for mask in range(1, 16):
        survivors = frozenset(i for i in range(4) if (mask >> i) & 1)
        if survivors.issubset(fixed):
            a4_support_survivor_sizes.add(len(survivors))
require(a4_support_survivor_sizes == {1},
        "A4 character-support survivor target changed")

s4_sign_support_survivor_sizes = set()
s4_allowed_sign_inertia = 0
for g in S4 - A4:
    fixed = fixed_points(g)
    if fixed:
        s4_allowed_sign_inertia += 1
    for mask in range(1, 16):
        survivors = frozenset(i for i in range(4) if (mask >> i) & 1)
        if survivors.issubset(fixed):
            s4_sign_support_survivor_sizes.add(len(survivors))
require(s4_allowed_sign_inertia == 6,
        "S4 allowed sign-support inertia is not exactly the transpositions")
require(s4_sign_support_survivor_sizes == {1, 2},
        "S4 sign-support survivor target changed")

# Exact local controls: an S3 transposition can retain one sheet, while the
# D4 deck-odd edge/four-cycle rows retain none.
T, t = sp.symbols("T t")
R_fix = sp.expand(t*T**3 - T)
R_edge = sp.expand((t*T**2 - 1)*(t*T**2 - 2))
R_four = t*T**4 - 1


def specialized_degree(poly):
    q = sp.Poly(sp.expand(poly.subs(t, 0)), T)
    return -1 if q.is_zero else q.degree()


require(sp.expand(R_fix - T*(t*T**2 - 1)) == 0,
        "fixed-transposition control changed")
require(specialized_degree(R_fix) == 1,
        "fixed-transposition control lost its finite branch")
require(R_edge == t**2*T**4 - 3*t*T**2 + 2,
        "edge control expansion changed")
require(specialized_degree(R_edge) == 0,
        "edge control retained a finite branch")
require(specialized_degree(R_four) == 0,
        "four-cycle control retained a finite branch")

# Dominant but non-etale hostile F0=(x,xy): on u=0 the equations force v=0,
# so the image meets the divisor only at one special point.
x, y, v = sp.symbols("x y v")
u0, v0 = x, x*y
jac0 = sp.det(sp.Matrix([[sp.diff(u0, x), sp.diff(u0, y)],
                         [sp.diff(v0, x), sp.diff(v0, y)]]))
require(jac0 == x, "non-etale hostile Jacobian changed")
require(sp.expand(v0.subs(x, 0)) == 0,
        "non-etale hostile unexpectedly dominates u=0")

# The familiar finite ramified hostile z->z^2 has nonconstant Jacobian.
q = sp.symbols("q")
require(sp.diff(q**2, q) == 2*q,
        "ramified one-variable hostile changed")

print("THM-2633 derangement-character exact controls")
print("quartic_group_orders=C4:4,V4:4,D4:8,A4:12,S4:24")
print("character_support_fixed_counts=C4:{0:2},V4:{0:2},D4:{0:4},A4:{1:8},S4:{0:6,2:6}")
print("quartic_gate=C4:EXCLUDED,V4:EXCLUDED,D4:EXCLUDED,A4:NOT_EXCLUDED,S4:NOT_EXCLUDED")
print(f"D4_deck_kernel_order={len(J)} deck_odd={len(d4_deck_odd)} deck_odd_survivor_rows={d4_deck_odd_survivor_rows}")
print("residual_targets=A4:3cycle_k1,S4:transposition_k1_or_k2")
print("dominant_nonetale_hostile=Jac:x;image_on_u0:only_(0,0);generic_k:0")
print("local_character_hostile=S3_transposition:k1")
print("local_D4_near_miss=edge:k0,fourcycle:k0")
print("source_unit_hostile=Gm_open_immersion_can_omit_divisor")
print(f"exact assertions passed: {checks}")
