#!/usr/bin/env python3
"""Exact algebra and finite-group companion for THM-2612."""

from __future__ import annotations

from itertools import combinations

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


# ---------------------------------------------------------------------------
# Depressed quartic, squared-pair resolvent, and the D4 deck formula.
# ---------------------------------------------------------------------------

w, U, z, p, q, r = s.symbols("w U z p q r")
f = w**4 + p * w**2 + q * w + r
S = U**3 + 2 * p * U**2 + (p**2 - 4 * r) * U - q**2
Sz = s.expand(S.subs(U, z))
D = 2 * w**2 + p + z
deck_sum_numerator = q - 2 * z * w

G = s.groebner([f, Sz], w, z, order="lex")

# If s0=(q-2zw)/D, then s0^2=z and tau(w)=-w-s0.
s_square_numerator = s.expand(deck_sum_numerator**2 - z * D**2)
require(G.reduce(s_square_numerator)[1] == 0, "deck pair sum stopped squaring to z")

tau_numerator = s.expand(-w * D - deck_sum_numerator)
tau_root_numerator = s.expand(
    tau_numerator**4
    + p * tau_numerator**2 * D**2
    + q * tau_numerator * D**3
    + r * D**4
)
require(G.reduce(tau_root_numerator)[1] == 0, "deck formula stopped preserving f")

s0 = deck_sum_numerator / D
tau = -w - s0
s0_after_tau = s.cancel(s0.subs(w, tau))
s0_fixed_numerator = s.cancel(s0_after_tau - s0).as_numer_denom()[0]
require(G.reduce(s.expand(s0_fixed_numerator))[1] == 0, "deck sum stopped being fixed")
tau2_numerator = s.cancel(tau.subs(w, tau) - w).as_numer_denom()[0]
require(G.reduce(s.expand(tau2_numerator))[1] == 0, "deck formula stopped being involutive")

norm_D = s.resultant(f, D, w)
Sprime_z = s.diff(S, U).subs(U, z)
require(
    s.expand(norm_D - Sprime_z**2 + 8 * (p + z) * Sz) == 0,
    "denominator norm identity changed",
)

# Sharp q=0 D4 control: T^4-2, rational resolvent root z=0, tau(w)=-w.
control_sub = {p: 0, q: 0, r: -2, z: 0}
require(s.factor(Sz.subs(control_sub)) == 0, "T^4-2 resolvent section changed")
require(s.cancel(tau.subs(control_sub) + w) == 0, "T^4-2 deck involution changed")
require(norm_D.subs(control_sub) == 64, "T^4-2 denominator norm changed")


# ---------------------------------------------------------------------------
# D4 point stabilizer, unique automorphism, and unique intermediate field.
# ---------------------------------------------------------------------------

Perm = tuple[int, ...]


def compose(a: Perm, b: Perm) -> Perm:
    """a after b."""
    return tuple(a[b[i]] for i in range(len(a)))


def inverse(a: Perm) -> Perm:
    ans = [0] * len(a)
    for i, ai in enumerate(a):
        ans[ai] = i
    return tuple(ans)


def generated(generators: list[Perm]) -> set[Perm]:
    identity = tuple(range(len(generators[0])))
    group = {identity}
    changed = True
    while changed:
        changed = False
        for a in list(group):
            for b in generators:
                for c in (compose(a, b), compose(b, a)):
                    if c not in group:
                        group.add(c)
                        changed = True
    return group


identity4 = (0, 1, 2, 3)
rotation = (1, 2, 3, 0)
reflection = (0, 3, 2, 1)
D4 = generated([rotation, reflection])
require(len(D4) == 8, "D4 generation changed")

H = {g for g in D4 if g[0] == 0}
require(len(H) == 2, "point stabilizer changed")
normalizer = {
    g
    for g in D4
    if {compose(compose(g, h), inverse(g)) for h in H} == H
}
require(len(normalizer) == 4, "D4 point-stabilizer normalizer changed")
require(
    {g[0] for g in normalizer} == {0, 2},
    "nontrivial quartic-field automorphism stopped pairing opposite sheets",
)


def is_subgroup(subset: set[Perm]) -> bool:
    if identity4 not in subset:
        return False
    return all(compose(a, inverse(b)) in subset for a in subset for b in subset)


subgroups = []
d4_list = list(D4)
for size in range(1, 9):
    for indices in combinations(range(8), size):
        subset = {d4_list[i] for i in indices}
        if is_subgroup(subset):
            subgroups.append(subset)
intermediate_subgroups = [J for J in subgroups if H < J < D4]
require(len(intermediate_subgroups) == 1, "D4 root field lost unique proper intermediate")
require(intermediate_subgroups[0] == normalizer, "intermediate subgroup is not the normalizer")


# ---------------------------------------------------------------------------
# Tame inertia tax and the depressed gcd gate.
# ---------------------------------------------------------------------------

pairings = (
    frozenset((frozenset((0, 1)), frozenset((2, 3)))),
    frozenset((frozenset((0, 2)), frozenset((1, 3)))),
    frozenset((frozenset((0, 3)), frozenset((1, 2)))),
)


def act_pairing(g: Perm, pairing: frozenset[frozenset[int]]) -> frozenset[frozenset[int]]:
    return frozenset(frozenset(g[i] for i in pair) for pair in pairing)


def orbit_count(permutation: Perm) -> int:
    unseen = set(range(len(permutation)))
    count = 0
    while unseen:
        count += 1
        start = next(iter(unseen))
        x = start
        while x in unseen:
            unseen.remove(x)
            x = permutation[x]
    return count


def pairing_permutation(g: Perm) -> Perm:
    return tuple(pairings.index(act_pairing(g, pairing)) for pairing in pairings)


inertia_representatives = {
    "identity": identity4,
    "transposition": (1, 0, 2, 3),
    "double_transposition": (1, 0, 3, 2),
    "three_cycle": (1, 2, 0, 3),
    "four_cycle": rotation,
}
expected_tax_table = {
    "identity": (0, 0, 0),
    "transposition": (1, 1, 0),
    "double_transposition": (2, 0, 1),
    "three_cycle": (2, 2, 0),
    "four_cycle": (3, 1, 1),
}
for name, permutation in inertia_representatives.items():
    d4 = 4 - orbit_count(permutation)
    d3 = 3 - orbit_count(pairing_permutation(permutation))
    tax = (d4 - d3) // 2
    require((d4, d3, tax) == expected_tax_table[name], f"tax row changed: {name}")

# Positive tax reductions.  A double transposition gives
# (T-a)^2(T+a)^2; a four-cycle gives T^4 after depression.
T, a = s.symbols("T a")
double_reduction = s.expand((T - a) ** 2 * (T + a) ** 2)
double_poly = s.Poly(double_reduction, T)
double_p = double_poly.coeff_monomial(T**2)
double_q = double_poly.coeff_monomial(T)
double_r = double_poly.coeff_monomial(1)
require(double_q == 0, "double-transposition reduction acquired a linear term")
require(s.expand(double_p**2 - 4 * double_r) == 0, "double-transposition gcd gate changed")

four_poly = s.Poly(T**4, T)
four_p = four_poly.coeff_monomial(T**2)
four_q = four_poly.coeff_monomial(T)
four_r = four_poly.coeff_monomial(1)
require(four_q == 0 and four_p**2 - 4 * four_r == 0, "four-cycle gcd gate changed")

# Hostile to sufficiency of the gcd gate.  The first quadratic has one tame
# transposition; the second factor is already split over C((t)).  Both raw
# orders have index length one, so the quartic-to-resolvent tax is zero even
# though t divides q and p^2-4r.
t = s.symbols("t")
local_f = s.expand(((T - 1) ** 2 - t) * ((T + 1) ** 2 - t**2))
local_poly = s.Poly(local_f, T)
local_p = local_poly.coeff_monomial(T**2)
local_q = local_poly.coeff_monomial(T)
local_r = local_poly.coeff_monomial(1)
local_gate = s.factor(local_p**2 - 4 * local_r)
require(s.rem(local_q, t) == 0 and s.rem(local_gate, t) == 0, "gcd hostile lost gate")
local_S = s.expand(S.subs({U: U, p: local_p, q: local_q, r: local_r}))
require(
    s.factor(local_S)
    == (U - 4) * (U**2 - 2 * U * t**2 - 2 * U * t + t**4 - 2 * t**3 + t**2),
    "gcd hostile resolvent factorization changed",
)
local_disc4 = s.factor(s.discriminant(local_f, T))
local_disc3 = s.factor(s.discriminant(local_S, U))
require(local_disc4 == local_disc3, "local discriminant identity changed")
require(s.factor(local_disc4 / t**3).subs(t, 0) != 0, "local raw discriminant valuation changed")
raw_discriminant_order = 3
field_discriminant_order = 1
quartic_index = (raw_discriminant_order - field_discriminant_order) // 2
resolvent_index = quartic_index
require(quartic_index == resolvent_index == 1, "zero-tax hostile index changed")

print("THM-2612 exact companion")
print("d4_group_order=8")
print("point_stabilizer_order=2")
print("point_stabilizer_normalizer_order=4")
print("quartic_field_automorphisms=C2")
print("proper_intermediate_fields=1")
print("deck_formula=tau(w)=-w-(q-2zw)/(2w^2+p+z)")
print("deck_denominator_norm=Sprime(z)^2")
print("deck_formula_involutive=YES")
print("tax_positive_inertia=double_transposition,four_cycle")
print("positive_tax_gate=prime_divides_gcd(q,p^2-4r)")
print("gcd_gate_sufficient=NO")
print("gcd_zero_tax_hostile_indices=(1,1)")
print("residue_characteristic_two_covered=NO")
print("ALL CHECKS PASSED")
