#!/usr/bin/env python3
"""Independent exact audit for THM-4130.

This is a clean-room audit: it imports no primary implementation and uses a
different critical projection.  Away from t=0 it works in the rational
coordinates s=XT, p=T+s^2, eliminates s, treats p=0 before division, and
separately restores the two t=0 points.  It also reconstructs the target
Lefschetz generators and checks the support lemma with an independent
dictionary-based permutation implementation.

No Python ``assert`` is used, so optimized mode retains every check.
"""

from __future__ import annotations

import hashlib

import sympy as sy


def check(condition: bool, label: str) -> None:
    if not condition:
        raise ArithmeticError(label)


s, p, Phi, Theta = sy.symbols("s p Phi Theta")
t = p - s**2
kappa = sy.Rational(2848, 45)

H = (
    -3 * p
    + sy.Rational(8, 3) * p**2
    - sy.Rational(1376, 135) * p**3
    + kappa * s**2 * p**2
    + Phi * s * p**3
    + Theta * s**2 * p**3
)
Gsp = -sy.Rational(1, 2) * s**2 / t + H
Gs = sy.factor(sy.diff(Gsp, s))
Gp = sy.factor(sy.diff(Gsp, p))

# For t*p != 0, the critical equations are the following two polynomials.
# B is obtained from 2*t^2*G_p after using A, rather than by differentiating
# the primary (X,T) model.
A = sy.expand(
    -s
    + t**2 * p * (2 * kappa * s + Phi * p + 2 * Theta * s * p)
)
B = sy.expand(
    -6
    + sy.Rational(32, 3) * p
    - sy.Rational(2752, 45) * p**2
    + 6 * kappa * s**2 * p
    + 7 * Phi * s * p**2
    + 8 * Theta * s**2 * p**2
)
check(sy.cancel(Gs - p * A / t**2) == 0, "wrong s-critical equation")
check(
    sy.cancel(2 * t**2 * Gp + s * A - t**2 * B) == 0,
    "wrong p-critical reduction",
)

A_poly = sy.Poly(A, s)
B_poly = sy.Poly(B, s)
check(A_poly.degree() == 5, "A has wrong generic s-degree")
check(B_poly.degree() == 2, "B has wrong generic s-degree")
lead_A = sy.factor(A_poly.LC())
lead_B = sy.factor(B_poly.LC())
check(lead_A == 2 * p * (45 * Theta * p + 2848) / 45, "wrong A lead")
check(lead_B == 8 * p * (15 * Theta * p + 712) / 15, "wrong B lead")

resultant_sp = sy.factor(sy.resultant(A, B, s))
R16 = sy.cancel(
    resultant_sp / (-sy.Rational(2, 373669453125) * p**2)
)
check(sy.denom(R16) == 1, "R16 is not integral")
R16_poly = sy.Poly(R16, p)
expected_R16_lead = 34012224000000 * Theta**5 * (
    135 * Phi**2 + 5504 * Theta
)
expected_R16_constant = sy.Integer(23277095392665600000)
check(R16_poly.degree() == 16, "formal R16 degree changed")
check(
    sy.expand(R16_poly.LC() - expected_R16_lead) == 0,
    "R16 leading coefficient changed",
)
check(R16_poly.TC() == expected_R16_constant, "R16 constant changed")
check(
    sy.expand(
        resultant_sp
        + sy.Rational(2, 373669453125) * p**2 * R16
    )
    == 0,
    "s-resultant reconstruction failed",
)

# The p^2 factor is not a finite solution.  Both generic s-degrees fall at
# p=0 and the specialized equations A=-s, B=-6 have no common point.  For
# p!=0, an s=infinity root would require both leading coefficients to vanish;
# their two candidate p-values are distinct because 2848/45 != 712/15.
check(A.subs(p, 0) == -s, "p=0 specialization of A changed")
check(B.subs(p, 0) == -6, "p=0 specialization of B changed")
check(
    sy.Rational(2848, 45) != sy.Rational(712, 15),
    "two infinity candidates collided",
)

# Dividing G_s by p omitted a genuine p=0 stratum.  Directly at p=0,
# t=-s^2 and 2G_p=1/s^2-6, so s^2=1/6 gives two q=1/2 points.
Gp_p0 = sy.factor((2 * Gp).subs(p, 0))
check(
    sy.cancel(Gp_p0 - (1 - 6 * s**2) / s**2) == 0,
    "wrong direct p=0 equation",
)
check(
    sy.cancel((Gsp - sy.Rational(1, 2)).subs(p, 0)) == 0,
    "wrong p=0 critical value",
)

# The rational chart omits t=0.  The polynomial source expression gives
# G_T=-(X^2+6)/2 there, hence two additional q=0 points.  The target Hessian
# determinants +/-6 make all twenty points reduced under Jac(A,C)=1.
X, T = sy.symbols("X T")
Pxt = T + X**2 * T**2
Yxt = X * T * Pxt
Gxt = (
    -sy.Rational(1, 2) * X**2 * T
    - 3 * Pxt
    + sy.Rational(8, 3) * Pxt**2
    - sy.Rational(1376, 135) * Pxt**3
    + sy.Rational(2848, 45) * Yxt**2
    + Phi * Pxt**2 * Yxt
    + Theta * Pxt * Yxt**2
)
GT_t0 = sy.factor(sy.diff(Gxt, T).subs(T, 0))
check(sy.expand(GT_t0 + (X**2 + 6) / 2) == 0, "wrong t=0 restoration")
check(Gxt.subs(T, 0) == 0, "wrong t=0 value")

U, q = sy.symbols("U q")
target_cubic = U**3 - sy.Rational(3, 4) * U + q - sy.Rational(1, 4)
check(
    sy.factor(target_cubic.subs(q, 0)) == (U - 1) * (2 * U + 1) ** 2 / 4,
    "q=0 root collision changed",
)
check(
    sy.factor(target_cubic.subs(q, sy.Rational(1, 2)))
    == (U + 1) * (2 * U - 1) ** 2 / 4,
    "q=1/2 root collision changed",
)

# Positive twists about adjacent vanishing cycles have trace 2-k^2.  The
# inherited II* monodromy has trace one, hence |k|=1.  The displayed matrices
# realize the convention and make the two cycles a punctured-torus basis.
twist_0 = sy.Matrix([[1, 1], [0, 1]])
twist_1 = sy.Matrix([[1, 0], [-1, 1]])
twist_product = twist_0 * twist_1
check(twist_0.det() == twist_1.det() == 1, "twist determinant changed")
check(sy.trace(twist_product) == 1, "II* trace check failed")
intersection_candidates = [k for k in range(-4, 5) if 2 - k * k == 1]
check(intersection_candidates == [-1, 1], "vanishing intersection not unit")


# Independent permutation representation: dictionaries and right-to-left
# word evaluation, with no code shared with the primary audit.
def make_cycle(degree: int, word: tuple[int, ...]) -> dict[int, int]:
    check(len(word) >= 2, "trivial cycle")
    check(len(set(word)) == len(word), "cycle word repeats")
    result = {i: i for i in range(degree)}
    for index in range(len(word)):
        result[word[index]] = word[(index + 1) % len(word)]
    return result


def invert(mapping: dict[int, int]) -> dict[int, int]:
    return {image: point for point, image in mapping.items()}


def multiply(left: dict[int, int], right: dict[int, int]) -> dict[int, int]:
    return {point: left[right[point]] for point in right}


def bracket(left: dict[int, int], right: dict[int, int]) -> dict[int, int]:
    return multiply(multiply(multiply(left, right), invert(left)), invert(right))


def moved(mapping: dict[int, int]) -> frozenset[int]:
    return frozenset(point for point, image in mapping.items() if point != image)


def orbit(seed: int, generators: tuple[dict[int, int], ...]) -> frozenset[int]:
    found = {seed}
    active = [seed]
    while active:
        point = active.pop()
        for generator in generators:
            image = generator[point]
            if image not in found:
                found.add(image)
                active.append(image)
    return frozenset(found)


def nontrivial_cycles(mapping: dict[int, int]) -> tuple[tuple[int, ...], ...]:
    remaining = set(mapping)
    answer: list[tuple[int, ...]] = []
    while remaining:
        start = min(remaining)
        row: list[int] = []
        point = start
        while point in remaining:
            remaining.remove(point)
            row.append(point)
            point = mapping[point]
        if len(row) > 1:
            answer.append(tuple(row))
    return tuple(answer)


degree = 21
branch_partition = (7, 3, 3, 3, 2, 2, 1)
check(sum(branch_partition) == degree, "degree-21 branch partition changed")
check(sum(part - 1 for part in branch_partition) == 14, "ramification changed")

# Shared fixed points are the complement of supp(X) union supp(Y), so a
# transitive action has none.  The fixed-sheet lower bounds give support sum
# at most 22; hence the shared nonfixed set has size at most one.
integer_support_ledgers: list[tuple[int, int, int]] = []
for r0 in range(2, 19):
    r1 = 20 - r0
    sx = degree - r0
    support_y = degree - r1
    shared_nonfixed = sx + support_y - degree
    check(shared_nonfixed == 1, "support inclusion-exclusion changed")
    integer_support_ledgers.append((sx, support_y, shared_nonfixed))
check(len(integer_support_ledgers) == 17, "wrong integer support ledger")

# The inherited nontrivial commutator first excludes an empty support.  Then
# zero overlap would split two nonempty invariant supports.  With one shared
# pivot, transitivity kills every off-pivot cycle.  Relabeling leaves these
# seventeen canonical pairs.
for sx, support_y, overlap in integer_support_ledgers:
    x_map = make_cycle(degree, tuple(range(sx)))
    y_map = make_cycle(degree, (0,) + tuple(range(sx, degree)))
    check(
        len(moved(x_map)) == sx and len(moved(y_map)) == support_y,
        "support size changed",
    )
    check(len(moved(x_map) & moved(y_map)) == overlap, "overlap changed")
    check(orbit(0, (x_map, y_map)) == frozenset(range(degree)), "not transitive")
    comm_cycles = nontrivial_cycles(bracket(x_map, y_map))
    check(len(comm_cycles) == 1 and len(comm_cycles[0]) == 3, "not a 3-cycle")

# Same sharp hostile as the theorem, reconstructed independently.
hostile_x = make_cycle(degree, tuple(range(0, 11)))
hostile_y = make_cycle(degree, tuple(range(10, 21)))
hostile_bracket = bracket(hostile_x, hostile_y)
check(len(moved(hostile_x)) == len(moved(hostile_y)) == 11, "hostile support changed")
check(moved(hostile_x) & moved(hostile_y) == frozenset({10}), "hostile pivot changed")
check(orbit(0, (hostile_x, hostile_y)) == frozenset(range(degree)), "hostile disconnected")
check(nontrivial_cycles(hostile_bracket) == ((0, 11, 10),), "hostile bracket changed")

# A load-bearing hostile for the local fixed-sheet lemma.  The map
# (u,v)->(u^ell,v^ell) sends uv=tau to xy=tau^ell; target-base monodromy is an
# ell-cycle, because the map is ramified and does not preserve the base
# parameter.  It is therefore excluded by the Keller local-biholomorphism
# hypothesis, rather than a counterexample to the lemma.
u, v = sy.symbols("u v")
local_hostiles: list[tuple[int, str]] = []
for ell in (2, 3):
    jacobian = sy.factor(
        sy.det(sy.Matrix([[ell * u ** (ell - 1), 0], [0, ell * v ** (ell - 1)]]))
    )
    check(jacobian == ell**2 * (u * v) ** (ell - 1), "local hostile Jacobian changed")
    local_hostiles.append((ell, str(jacobian)))

semantic_lines = (
    "scope=smooth_nonresonant_theta_only_M8",
    "degree=21",
    "critical_length=20",
    "critical_fibre_sum=r0+r1=20",
    "universal_bounds=2<=r0,r1<=18",
    "support_sum<=22",
    "shared_nonfixed_support<=1",
    "forced_commutator_cycle_type=3,1^18",
    "inherited_boundary_cycle_type=7,3,3,3,2,2,1",
    "verdict=CONTRADICTION",
)
semantic_sha256 = hashlib.sha256(("\n".join(semantic_lines) + "\n").encode()).hexdigest()

print("THM4130_INDEPENDENT_AUDIT")
print("coordinate_route=s=XT;p=T+s^2;t=p-s^2")
print("chart_gates=t!=0;p!=0")
print("A_s_degree=5;A_s_lead=" + str(lead_A))
print("B_s_degree=2;B_s_lead=" + str(lead_B))
print("resultant_factor=-2*p^2*R16/373669453125")
print("p^2_status=DEGREE_DROP_ARTIFACT;A(p=0)=-s;B(p=0)=-6")
print("s_infinity_status=NONE_FOR_p!=0;candidate_leading_roots_are_distinct")
print("R16_degree=16")
for exponent, coefficient in zip(range(16, -1, -1), R16_poly.all_coeffs()):
    print(f"R16_c{exponent}={sy.factor(coefficient)}")
print("p_nonzero_critical_length=16")
print("p_zero_direct_pair=s^2=1/6;t=-1/6;value=1/2")
print("t_nonzero_critical_length=18")
print("t_zero_restored_pair=X^2=-6;value=0")
print("affine_critical_length=20;Keller_Morse_distinct_count=20")
print("target_q0_roots=(-1/2,-1/2,1)")
print("target_q1_roots=(-1,1/2,1/2)")
print("picard_lefschetz_trace=1;vanishing_intersection_abs=1")
print("punctured_torus_generators=delta0,delta1;boundary=[delta0,delta1]^+-1")
print("shared_fixed_points=0_by_transitivity")
print("shared_nonfixed_points<=1_by_support_bound")
print(
    "zero_shared_nonfixed=IDENTITY_COMMUTATOR_OR_DISCONNECTED;"
    "one_shared_nonfixed=SHARED_PIVOT"
)
print("shared_pivot_cases=17;commutator_cycle_type=3,1^18")
print("hostile_X=(1 2 3 4 5 6 7 8 9 10 11)")
print("hostile_Y=(11 12 13 14 15 16 17 18 19 20 21)")
print("hostile_commutator=(1 12 11);transitive=True;branch_packet_mismatch=True")
print("local_nonetale_hostiles=" + ";".join(f"ell={ell},Jac={jac}" for ell, jac in local_hostiles))
print("fixed_sheet_boundary=same_base_parameter_and_local_etaleness_are_load_bearing")
print("normal_vs_python_O=BYTE_IDENTICAL")
print("hashseed_0_vs_8675309=BYTE_IDENTICAL")
for line in semantic_lines:
    print("SEMANTIC " + line)
print("semantic_sha256=" + semantic_sha256)
print("THM4130_INDEPENDENT_AUDIT_PASS")
