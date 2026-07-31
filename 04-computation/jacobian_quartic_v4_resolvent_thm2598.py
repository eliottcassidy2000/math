# -*- coding: utf-8 -*-
"""Exact companion for THM-2598: quartic V4 resolvent transfer boundary."""

from itertools import permutations

import sympy as sp


failures = []


def require(label, condition):
    """Optimization-safe exact check."""
    ok = bool(condition)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        failures.append(label)
        raise RuntimeError(f"failed exact check: {label}")


def zero(label, expression):
    require(label, sp.expand(expression) == 0)


T, W, U, Z = sp.symbols("T W U Z")
A, B, C, D, E = sp.symbols("A B C D E")
p, q, r, s, t = sp.symbols("p q r s t")

print("THM-2598 QUARTIC V4 RESOLVENT AUDIT")

# General integral resolvent and invariant depression.
f = A*T**4 + B*T**3 + C*T**2 + D*T + E
R = W**3 - C*W**2 + (B*D - 4*A*E)*W + 4*A*C*E - B**2*E - A*D**2
disc_f = sp.discriminant(f, T)
disc_R = sp.discriminant(R, W)
zero("general integral resolvent has identical discriminant", disc_R - disc_f)

I = C**2 - 3*B*D + 12*A*E
J = 2*C**3 - 9*B*C*D - 72*A*C*E + 27*B**2*E + 27*A*D**2
zero("resolvent depression is Z^3-I Z/3-J/27",
     R.subs(W, Z + C/sp.Integer(3)) - (Z**3 - I*Z/3 - J/27))
zero("universal cusp identity 27 Disc=4I^3-J^2", 27*disc_f - 4*I**3 + J**2)

# Root provenance under Vieta, including the root-difference proof.
x1, x2, x3, x4 = sp.symbols("x1 x2 x3 x4")
roots = (x1, x2, x3, x4)
e1 = sum(roots)
e2 = sum(roots[i]*roots[j] for i in range(4) for j in range(i + 1, 4))
e3 = sum(roots[i]*roots[j]*roots[k]
         for i in range(4) for j in range(i + 1, 4) for k in range(j + 1, 4))
e4 = sp.prod(roots)
vieta = {B: -A*e1, C: A*e2, D: -A*e3, E: A*e4}
betas = (
    A*(x1*x2 + x3*x4),
    A*(x1*x3 + x2*x4),
    A*(x1*x4 + x2*x3),
)
root_poly = sp.prod(W - beta for beta in betas)
zero("resolvent roots are the three complementary pair products",
     root_poly - R.subs(vieta))
zero("pairing difference factors into two quartic differences",
     betas[0] - betas[1] - A*(x1 - x4)*(x2 - x3))

# Depressed squared-pairing resolvent and homogeneous leading drop.
fA = A*T**4 + p*T**2 + q*T + r
SA = U**3 + 2*p*U**2 + (p**2 - 4*A*r)*U - A*q**2
zero("homogeneous squared-pairing resolvent has quartic discriminant",
     sp.discriminant(SA, U) - sp.discriminant(fA, T))
zero("leading drop is U(U+p)^2", SA.subs(A, 0) - U*(U + p)**2)
disc_expansion = A*(4*p**3*(4*p*r - q**2)
                    + A*(-128*p**2*r**2 + 144*p*q**2*r - 27*q**4)
                    + 256*A**2*r**3)
zero("exact leading-factor discriminant expansion",
     sp.discriminant(fA, T) - disc_expansion)

# V4 reconstruction on a fully explicit nondegenerate quartic.
Sf = U**3 - 14*U**2 + 49*U - 36
zero("f+ and f- have the same squared-pairing resolvent",
     Sf - (U - 1)*(U - 4)*(U - 9))
f_plus = T**4 - 7*T**2 + 6*T
f_minus = T**4 - 7*T**2 - 6*T
zero("f+ root factorization", f_plus - T*(T - 1)*(T - 2)*(T + 3))
zero("f- is the marked reversal f+(-T)", f_minus - f_plus.subs(T, -T))

sign_triples = ((1, 2, -3), (1, -2, 3), (-1, 2, 3), (-1, -2, -3))
reconstructed = []
for a, b, c in sign_triples:
    require("sign triple has fixed squares and product -q",
            (a*a, b*b, c*c, a*b*c) == (1, 4, 9, -6))
    reconstructed.append(((a+b+c)//2, (a-b-c)//2,
                          (-a+b-c)//2, (-a-b+c)//2))
expected = ((0, 1, 2, -3), (1, 0, -3, 2),
            (2, -3, 0, 1), (-3, 2, 1, 0))
require("four reconstruction sections are exactly the V4 relabellings",
        tuple(reconstructed) == expected)
zero("first biquadratic quotient algebra is split",
     (U**3 - 20*U**2 + 96*U) - U*(U - 8)*(U - 12))
zero("second biquadratic quotient algebra is split",
     (U**3 - 28*U**2 + 160*U) - U*(U - 8)*(U - 20))

pairings = (
    frozenset((frozenset((0, 1)), frozenset((2, 3)))),
    frozenset((frozenset((0, 2)), frozenset((1, 3)))),
    frozenset((frozenset((0, 3)), frozenset((1, 2)))),
)


def apply_pairing(perm, pairing):
    return frozenset(frozenset((perm[a], perm[b])) for a, b in pairing)


# Exhaustive transitive-subgroup and field/deck typing.  Permutations are
# tuples mapping old labels to new labels; composition is left after right.
ID4 = (0, 1, 2, 3)
ID3 = (0, 1, 2)
ALL_S4 = set(permutations(range(4)))


def compose(left, right):
    return tuple(left[right[i]] for i in range(len(right)))


def inverse(perm):
    answer = [0] * len(perm)
    for old, new in enumerate(perm):
        answer[new] = old
    return tuple(answer)


def generated(*generators):
    group = {ID4}
    frontier = [ID4]
    while frontier:
        old = frontier.pop()
        for generator in generators:
            new = compose(old, generator)
            if new not in group:
                group.add(new)
                frontier.append(new)
    return group


def orbit_sizes(group, degree, action=lambda g: g):
    unseen = set(range(degree))
    sizes = []
    while unseen:
        seed = min(unseen)
        orbit = {action(g)[seed] for g in group}
        sizes.append(len(orbit))
        unseen -= orbit
    return tuple(sorted(sizes))


def matching_action(perm):
    return tuple(pairings.index(apply_pairing(perm, pairing)) for pairing in pairings)


def point_stabilizer(group, point=0):
    return {g for g in group if g[point] == point}


def normalizer(group, subgroup):
    answer = set()
    for g in group:
        gi = inverse(g)
        if {compose(compose(g, h), gi) for h in subgroup} == subgroup:
            answer.add(g)
    return answer


def conjugate_group(group, by):
    byi = inverse(by)
    return {compose(compose(by, g), byi) for g in group}


def all_subgroups_of_s4():
    trivial = frozenset((ID4,))
    known = {trivial}
    frontier = [trivial]
    while frontier:
        subgroup = frontier.pop()
        for element in ALL_S4 - set(subgroup):
            enlarged = frozenset(generated(*(tuple(subgroup) + (element,))))
            if enlarged not in known:
                known.add(enlarged)
                frontier.append(enlarged)
    return known


def parity(perm):
    return sum(perm[i] > perm[j] for i in range(4) for j in range(i + 1, 4)) % 2


r4 = (1, 2, 3, 0)
s_vertex = (0, 3, 2, 1)
dt1 = (1, 0, 3, 2)
dt2 = (2, 3, 0, 1)
C4 = generated(r4)
V4 = generated(dt1, dt2)
D4 = generated(r4, s_vertex)
A4 = {perm for perm in ALL_S4 if parity(perm) == 0}
S4 = ALL_S4

group_rows = (
    ("C4", C4, 2, (1, 2), 1, 4),
    ("V4", V4, 1, (1, 1, 1), 1, 4),
    ("D4", D4, 2, (1, 2), 2, 2),
    ("A4", A4, 3, (3,), 3, 1),
    ("S4", S4, 6, (3,), 6, 1),
)
for name, group, image_order, matching_orbits, stabilizer_order, deck_order in group_rows:
    image = {matching_action(g) for g in group}
    H = point_stabilizer(group)
    deck = normalizer(group, H)
    require(f"transitive branch {name}: image/orbits/stabilizer/deck",
            orbit_sizes(group, 4) == (4,)
            and len(image) == image_order
            and orbit_sizes(group, 3, matching_action) == matching_orbits
            and len(H) == stabilizer_order
            and len(deck) // len(H) == deck_order)

subgroups = all_subgroups_of_s4()
transitive_subgroups = [set(group) for group in subgroups
                        if orbit_sizes(set(group), 4) == (4,)]
types_seen = set()
classification_ok = len(subgroups) == 30
for subgroup in transitive_subgroups:
    matches = [name for name, named_group, *_ in group_rows
               if any(conjugate_group(named_group, by) == subgroup for by in S4)]
    classification_ok = classification_ok and len(matches) == 1
    types_seen.update(matches)
require("all 30 S4 subgroups give exactly five transitive conjugacy types",
        classification_ok and len(transitive_subgroups) == 9
        and types_seen == {"C4", "V4", "D4", "A4", "S4"})

kernel = {g for g in S4 if matching_action(g) == ID3}
require("matching kernel is exactly V4 and every S3 label has four lifts",
        kernel == V4 and len({matching_action(g) for g in S4}) == 6
        and all(sum(matching_action(g) == image for g in S4) == 4
                for image in {matching_action(g) for g in S4}))

# In the D4 branch, the matching quadratic and the quadratic intermediate
# inside the quartic root field are different subfields of the closure.
H_D4 = point_stabilizer(D4)
root_intermediate_stabilizer = normalizer(D4, H_D4)
require("D4 generic deck C2 needs polynomial extension",
        len(H_D4) == 2 and len(root_intermediate_stabilizer) == 4)
require("D4 matching quadratic is not the root-field quadratic intermediate",
        H_D4.isdisjoint(V4 - {ID4})
        and generated(*(tuple(H_D4 | V4))) == D4
        and root_intermediate_stabilizer != V4)

# The unique rational D4 matching is the antipodal matching of its canonical
# square, and the generic deck C2 is its central half-turn.  This identifies
# the exact owner pair that a non-polynomial deck involution must split.
D4_fixed_pairings = [index for index in range(3)
                     if all(matching_action(g)[index] == index for g in D4)]
D4_center = {g for g in D4 if all(compose(g, h) == compose(h, g) for h in D4)}
D4_half_turns = D4_center - {ID4}
D4_channel_stabilizer = {g for g in S4
                         if apply_pairing(g, pairings[D4_fixed_pairings[0]])
                         == pairings[D4_fixed_pairings[0]]}
D4_half_turn = next(iter(D4_half_turns))
D4_antipodes = frozenset(
    frozenset((vertex, D4_half_turn[vertex])) for vertex in range(4)
)
require("D4 rational matching is the canonical square's antipodal matching",
        len(D4_fixed_pairings) == 1
        and len(D4_center) == 2
        and D4_channel_stabilizer == D4
        and D4_antipodes == pairings[D4_fixed_pairings[0]])
require("D4 generic deck C2 is the central antipodal half-turn",
        D4_half_turn in root_intermediate_stabilizer - H_D4
        and D4_half_turn in V4
        and matching_action(D4_half_turn) == ID3
        and all(compose(D4_half_turn, g) == compose(g, D4_half_turn)
                for g in D4))

# In the irreducible A4/S4 rows the quartic and matching fields are likewise
# incomparable: their stabilizers are incomparable and generate the closure.
for name, group in (("A4", A4), ("S4", S4)):
    H_root = point_stabilizer(group)
    H_matching = {g for g in group if matching_action(g)[0] == 0}
    require(f"{name} quartic and matching fields are incomparable",
            not H_root <= H_matching and not H_matching <= H_root
            and generated(*(tuple(H_root | H_matching))) == group)

# Tame inertia ledger: action on four roots and on three pairings.
def cycle_count(perm):
    seen = set()
    count = 0
    for start in range(len(perm)):
        if start not in seen:
            count += 1
            cur = start
            while cur not in seen:
                seen.add(cur)
                cur = perm[cur]
    return count


representatives = {
    "identity": (0, 1, 2, 3),
    "transposition": (1, 0, 2, 3),
    "double_transposition": (1, 0, 3, 2),
    "three_cycle": (1, 2, 0, 3),
    "four_cycle": (1, 2, 3, 0),
}
expected_ledger = {
    "identity": (0, 0, 0),
    "transposition": (1, 1, 0),
    "double_transposition": (2, 0, 1),
    "three_cycle": (2, 2, 0),
    "four_cycle": (3, 1, 1),
}
for name, perm in representatives.items():
    induced = tuple(pairings.index(apply_pairing(perm, pairing)) for pairing in pairings)
    d4 = 4 - cycle_count(perm)
    d3 = 3 - cycle_count(induced)
    require(f"inertia ledger {name}: (d4,d3,index tax)",
            (d4, d3, (d4 - d3)//2) == expected_ledger[name])

# Sharp double-transposition family and generic S4 specialization.
fst = (T**2 - 1)**2 - t*(T + s)
Rst = W**3 + 2*W**2 + (4*s*t - 4)*W + 8*s*t - t**2 - 8
Dst = -t**2*(256*s**3*t - 256*s**2 - 288*s*t + 27*t**2 + 256)
zero("double-transposition family quartic discriminant", sp.discriminant(fst, T) - Dst)
zero("double-transposition family resolvent discriminant", sp.discriminant(Rst, W) - Dst)
zero("t=0 resolvent is one simple plus one double root",
     Rst.subs(t, 0) - (W + 2)**2*(W - 2))
zero("normalized double-transposition chart has separable residual quadratic",
     sp.expand(Rst.subs(W, -2 + t*Z)/t**2)
     - (t*Z**3 - 4*Z**2 + 4*s*Z - 1))

special = sp.Poly(fst.subs({s: 2, t: 1}), T)
degrees_mod2 = sorted(poly.degree() for poly, exponent in sp.factor_list(special, modulus=2)[1]
                      for _ in range(exponent))
degrees_mod3 = sorted(poly.degree() for poly, exponent in sp.factor_list(special, modulus=3)[1]
                      for _ in range(exponent))
require("S4 witness has Frobenius type 4 mod 2", degrees_mod2 == [4])
require("S4 witness has Frobenius type 1+3 mod 3", degrees_mod3 == [1, 3])

# Minimal inertia and finite-cusp controls.
zero("four-cycle control resolvent", (W**3 + 4*t*W) - W*(W**2 + 4*t))
zero("four-cycle control discriminant", sp.discriminant(T**4 - t, T) + 256*t**3)
R3 = W**3 - 3*t*W - t - t**2
zero("three-cycle control discriminant",
     sp.discriminant((T**3 - t)*(T - 1), T) + 27*t**2*(t - 1)**2)
zero("three-cycle resolvent discriminant",
     sp.discriminant(R3, W) + 27*t**2*(t - 1)**2)
u, v = sp.symbols("u v")
zero("finite cusp control has common discriminant",
     sp.discriminant(T**4 + u*T + v, T) - sp.discriminant(W**3 - 4*v*W - u**2, W))
zero("finite cusp equation", sp.discriminant(T**4 + u*T + v, T) - (256*v**3 - 27*u**4))

print("FAILED CHECKS:", failures if failures else "NONE")
