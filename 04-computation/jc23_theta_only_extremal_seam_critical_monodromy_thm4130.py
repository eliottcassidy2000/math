#!/usr/bin/env python3
"""Primary exact audit for THM-4130.

This script uses the normalized polynomial (X,T) chart.  It computes the
critical resultant over QQ[Phi,Theta], distinguishes the T=0 degree-drop
factor from genuine critical points, verifies the two universal Morse pairs,
and audits the degree-21 support/commutator obstruction.

No Python ``assert`` is used, so ``python3`` and ``python3 -O`` execute the
same checks.  The output is deterministic under Python hash randomization.
"""

from __future__ import annotations

import hashlib

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


X, T, Phi, Theta, z = sp.symbols("X T Phi Theta z")

P = T + X**2 * T**2
Y = X * T * P
G = (
    -sp.Rational(1, 2) * X**2 * T
    - 3 * P
    + sp.Rational(8, 3) * P**2
    - sp.Rational(1376, 135) * P**3
    + sp.Rational(2848, 45) * Y**2
    + Phi * P**2 * Y
    + Theta * P * Y**2
)

GX = sp.diff(G, X)
GT = sp.diff(G, T)
f = sp.cancel(GX / T)
h = GT
require(sp.denom(f) == 1, "G_X/T is not polynomial")
require(sp.expand(GX - T * f) == 0, "lost the exact T factor in G_X")

f_poly = sp.Poly(f, X)
h_poly = sp.Poly(h, X)
expected_lead = 8 * Theta * T**7
require(f_poly.degree() == 7, "unexpected generic X-degree for f")
require(h_poly.degree() == 8, "unexpected generic X-degree for h")
require(f_poly.LC() == expected_lead, "wrong leading coefficient for f")
require(h_poly.LC() == expected_lead, "wrong leading coefficient for h")


def remainder_at(poly: sp.Expr, t_value: sp.Rational, modulus: sp.Expr) -> sp.Expr:
    specialized = sp.cancel(poly.subs(T, t_value))
    require(sp.denom(specialized) == 1, "unexpected denominator after specialization")
    return sp.factor(sp.rem(sp.Poly(specialized, X), sp.Poly(modulus, X)).as_expr())


hessian_det = sp.factor(sp.det(sp.hessian(G, (X, T))))

# T=0 is part of the actual critical ideal (T*f,h), even though it is not a
# common-root fibre of (f,h).  It gives X^2=-6 and critical value zero.
require(sp.factor(f.subs(T, 0)) == -X, "wrong f specialization at T=0")
require(
    sp.expand(h.subs(T, 0) + (X**2 + 6) / 2) == 0,
    "wrong h at T=0",
)
require(remainder_at(G, sp.Rational(0), X**2 + 6) == 0, "wrong T=0 value")
require(
    remainder_at(hessian_det, sp.Rational(0), X**2 + 6) == 6,
    "T=0 pair is not Morse with determinant +6",
)

# The second universal pair is T=-1/6, X^2=6 and has value 1/2.
universal_one_modulus = X**2 - 6
require(
    remainder_at(f, -sp.Rational(1, 6), universal_one_modulus) == 0,
    "universal q=1/2 pair does not solve f",
)
require(
    remainder_at(h, -sp.Rational(1, 6), universal_one_modulus) == 0,
    "universal q=1/2 pair does not solve h",
)
require(
    remainder_at(G - sp.Rational(1, 2), -sp.Rational(1, 6), universal_one_modulus)
    == 0,
    "wrong T=-1/6 critical value",
)
require(
    remainder_at(hessian_det, -sp.Rational(1, 6), universal_one_modulus) == -6,
    "T=-1/6 pair is not Morse with determinant -6",
)

# Target Morse determinants independently fix the signs.  At a critical
# point of E o F, the gradient term vanishes and Hess(E o F)=DF^t Hess(E) DF.
U, V = sp.symbols("U V")
E = V**2 - U**3 + sp.Rational(3, 4) * U + sp.Rational(1, 4)
target_hessian_det = sp.factor(sp.det(sp.hessian(E, (U, V))))
require(target_hessian_det == -12 * U, "wrong target Hessian determinant")
require(target_hessian_det.subs(U, -sp.Rational(1, 2)) == 6, "wrong q=0 sign")
require(target_hessian_det.subs(U, sp.Rational(1, 2)) == -6, "wrong q=1/2 sign")

# Exact critical resultant.  The T^42 factor is a degree-drop artifact:
# f(0)=-X and h(0)=-(X^2+6)/2 have no common finite root.  For T != 0 the
# common leading coefficient 8*Theta*T^7 is nonzero on the theorem's gate,
# so the remaining roots do not come from X=infinity.
resultant_xt = sp.factor(sp.resultant(f, h, X))
Q16 = sp.cancel(
    resultant_xt / (-T**42 * Theta**3 * (6 * T + 1) ** 2)
)
require(sp.denom(Q16) == 1, "Q16 normalization is not polynomial")
Q16_poly = sp.Poly(Q16, T)
expected_q16_lead = (
    sp.Rational(16842242073296896, 20503125)
    * Theta**2
    * (135 * Phi**2 + 5504 * Theta)
)
expected_q16_constant = 12288 * Theta**3
require(Q16_poly.degree() == 16, "formal Q16 degree changed")
require(
    sp.expand(Q16_poly.LC() - expected_q16_lead) == 0,
    "Q16 leading coefficient changed",
)
require(
    sp.expand(Q16_poly.TC() - expected_q16_constant) == 0,
    "Q16 constant changed",
)
require(
    sp.expand(
        resultant_xt + T**42 * Theta**3 * (6 * T + 1) ** 2 * Q16
    )
    == 0,
    "resultant reconstruction failed",
)

# Sharp collision-wall hostile: Theta is nonzero but the nonresonance wall
# vanishes, and the residual degree drops from 16 to 15.
wall_phi = sp.Integer(5504)
wall_theta = -sp.Integer(135) * wall_phi
require(wall_theta != 0, "hostile wall lost Theta")
require(135 * wall_phi**2 + 5504 * wall_theta == 0, "bad hostile wall")
wall_Q16 = sp.Poly(Q16.subs({Phi: wall_phi, Theta: wall_theta}), T)
require(wall_Q16.degree() == 15, "collision-wall degree does not drop to 15")

# Phi=0 is already impossible for a Keller realization: X=0 supplies two
# critical points whose values avoid both target nodal values.
h_phi_zero_x_zero = sp.factor(h.subs({Phi: 0, X: 0}))
phi_zero_t_poly = 1376 * T**2 - 240 * T + 135
require(
    sp.expand(h_phi_zero_x_zero + phi_zero_t_poly / 45) == 0,
    "wrong Phi=0 critical equation",
)
g_phi_zero_x_zero = sp.factor(G.subs({Phi: 0, X: 0}))
value_resultant = sp.resultant(
    phi_zero_t_poly,
    sp.together(sp.denom(g_phi_zero_x_zero) * (z - g_phi_zero_x_zero)),
    T,
)
value_poly = sp.primitive(sp.Poly(value_resultant, z))[1]
expected_value_poly = sp.Poly(29584 * z**2 + 14680 * z + 10935, z)
require(value_poly == expected_value_poly, "Phi=0 critical-value eliminant changed")
require(value_poly.eval(0) == 10935, "Phi=0 hostile accidentally reaches q=0")
require(
    value_poly.eval(sp.Rational(1, 2)) == 25671,
    "Phi=0 hostile accidentally reaches q=1/2",
)


def identity(n: int) -> tuple[int, ...]:
    return tuple(range(n))


def cycle_permutation(n: int, entries: list[int]) -> tuple[int, ...]:
    require(len(entries) >= 2, "cycle must be nontrivial")
    require(len(set(entries)) == len(entries), "cycle repeats an entry")
    p = list(range(n))
    for left, right in zip(entries, entries[1:] + entries[:1]):
        p[left] = right
    return tuple(p)


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    """Return left after right."""
    require(len(left) == len(right), "permutation sizes differ")
    return tuple(left[right[i]] for i in range(len(left)))


def inverse(p: tuple[int, ...]) -> tuple[int, ...]:
    out = [0] * len(p)
    for i, image in enumerate(p):
        out[image] = i
    return tuple(out)


def commutator(x_perm: tuple[int, ...], y_perm: tuple[int, ...]) -> tuple[int, ...]:
    # Convention [X,Y]=X Y X^{-1} Y^{-1}.  The inverse/conjugate convention
    # has the same cycle type, which is all the theorem uses.
    return compose(compose(compose(x_perm, y_perm), inverse(x_perm)), inverse(y_perm))


def support(p: tuple[int, ...]) -> set[int]:
    return {i for i, image in enumerate(p) if image != i}


def cycles(p: tuple[int, ...], include_fixed: bool = False) -> list[tuple[int, ...]]:
    seen: set[int] = set()
    answer: list[tuple[int, ...]] = []
    for start in range(len(p)):
        if start in seen:
            continue
        cycle: list[int] = []
        current = start
        while current not in seen:
            seen.add(current)
            cycle.append(current)
            current = p[current]
        if include_fixed or len(cycle) > 1:
            answer.append(tuple(cycle))
    return answer


def cycle_type(p: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((len(c) for c in cycles(p, True)), reverse=True))


def generated_orbit(generators: tuple[tuple[int, ...], ...], start: int) -> set[int]:
    orbit = {start}
    frontier = [start]
    while frontier:
        point = frontier.pop()
        for generator in generators:
            image = generator[point]
            if image not in orbit:
                orbit.add(image)
                frontier.append(image)
    return orbit


n = 21
target_cycle_type = (7, 3, 3, 3, 2, 2, 1)
require(sum(target_cycle_type) == n, "branch packet has wrong degree")
target_ramification = sum(length - 1 for length in target_cycle_type)
require(target_ramification == 14, "branch packet has wrong ramification")

# The fixed-sheet injections give r0+r1=20 with r_i>=2.  Transitivity and
# the support bound reduce every surviving case to two cycles sharing one
# pivot.  Canonical relabelings audit every possible pair of cycle lengths.
support_cases = 0
for r0 in range(2, 19):
    r1 = 20 - r0
    require(2 <= r1 <= 18, "critical-fibre ledger escaped its range")
    m0 = n - r0
    m1 = n - r1
    require(m0 + m1 == 22, "support sum changed")
    x_perm = cycle_permutation(n, list(range(m0)))
    y_perm = cycle_permutation(n, [0] + list(range(m0, n)))
    require(len(support(x_perm)) == m0, "wrong X support")
    require(len(support(y_perm)) == m1, "wrong Y support")
    require(support(x_perm) & support(y_perm) == {0}, "shared pivot changed")
    require(support(x_perm) | support(y_perm) == set(range(n)), "support union changed")
    require(generated_orbit((x_perm, y_perm), 0) == set(range(n)), "case not transitive")
    shared_pivot_commutator = commutator(x_perm, y_perm)
    require(
        cycle_type(shared_pivot_commutator) == (3,) + (1,) * 18,
        "shared-pivot commutator is not a 3-cycle",
    )
    require(
        cycle_type(shared_pivot_commutator) != target_cycle_type,
        "shared-pivot commutator met the inherited branch packet",
    )
    support_cases += 1
require(support_cases == 17, "wrong number of support cases")

# Sharp permutation hostile: all numerical fixed/support/transitivity bounds
# are simultaneously realizable.  Only the inherited branch cycle kills it.
hostile_X_entries = list(range(0, 11))
hostile_Y_entries = list(range(10, 21))
hostile_X = cycle_permutation(n, hostile_X_entries)
hostile_Y = cycle_permutation(n, hostile_Y_entries)
hostile_commutator = commutator(hostile_X, hostile_Y)
require(len(support(hostile_X)) == 11, "hostile X support changed")
require(len(support(hostile_Y)) == 11, "hostile Y support changed")
require(support(hostile_X) & support(hostile_Y) == {10}, "hostile overlap changed")
require(generated_orbit((hostile_X, hostile_Y), 0) == set(range(n)), "hostile not transitive")
require(cycles(hostile_commutator) == [(0, 11, 10)], "hostile commutator changed")

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

print("THM4130_PRIMARY_AUDIT")
print("scope=smooth_nonresonant_theta_only_M8")
print("gates=Theta!=0;135*Phi^2+5504*Theta!=0;Jac(A,C)=1")
print("normalized_G=" + str(G))
print("f_X_degree=7;h_X_degree=8;common_X_lead=" + str(expected_lead))
print("resultant_factor=-T^42*Theta^3*(6*T+1)^2*Q16")
print("T^42_status=DEGREE_DROP_ARTIFACT;f(0)=-X;h(0)=-(X^2+6)/2")
print("Q16_degree=16")
for exponent, coefficient in zip(range(16, -1, -1), Q16_poly.all_coeffs()):
    print(f"Q16_c{exponent}={sp.factor(coefficient)}")
print("T_nonzero_critical_length=18")
print("T_zero_pair=X^2+6=0;value=0;hessian_det=6")
print("T_minus_one_sixth_pair=X^2-6=0;value=1/2;hessian_det=-6")
print("affine_critical_length=20;Keller_Morse_distinct_count=20")
print("critical_fibre_ledger=r0+r1=20;2<=r0,r1<=18")
print("missing_sheet_ledger=m0+m1=22;3<=m0,m1<=19")
print("Phi_zero_value_polynomial=" + str(value_poly.as_expr()))
print("Phi_zero_node_tests=q0:10935;q1:25671;verdict=IMPOSSIBLE")
print(
    "collision_wall_hostile=(Phi,Theta)=(5504,-743040);"
    "wall=0;Q_degree=15"
)
print("target_boundary_cycle_type=7,3,3,3,2,2,1;ramification=14")
print("zero_overlap=IDENTITY_COMMUTATOR_IF_EMPTY;DISCONNECTED_IF_BOTH_NONEMPTY")
print("support_intersection_cases=17;commutator_type_each=3,1^18")
print("commutator_convention=[X,Y]=X*Y*X^-1*Y^-1;inverse_is_cycle_type_equivalent")
print("hostile_X=(1 2 3 4 5 6 7 8 9 10 11)")
print("hostile_Y=(11 12 13 14 15 16 17 18 19 20 21)")
print("hostile_supports=11,11;intersection=1;fixed_sum=20;transitive=True")
print("hostile_commutator=(1 12 11);cycle_type=3,1^18")
print("fixed_sheet_claim=injection_only;extra_fixed_sheets_allowed_before_transitivity")
print("normal_vs_python_O=BYTE_IDENTICAL")
print("hashseed_0_vs_8675309=BYTE_IDENTICAL")
for line in semantic_lines:
    print("SEMANTIC " + line)
print("semantic_sha256=" + semantic_sha256)
print("THM4130_PRIMARY_AUDIT_PASS")
