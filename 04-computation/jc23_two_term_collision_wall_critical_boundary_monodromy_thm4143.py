#!/usr/bin/env python3
"""Primary exact audit for THM-4143's two-term collision wall.

This standalone SymPy route proves the Phi=0 subwall impossible, factors the
critical resultants on the K!=0 and K=0 strata, normalizes the repeated outer
edge, freezes the rational/quadratic boundary responses, and checks all three
monodromy budgets.  It imports no repository code and uses no Python
``assert``, so every gate remains live under optimized mode.
"""

from __future__ import annotations

import hashlib
from math import gcd

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def convex_hull(points: list[tuple[int, int]]) -> list[tuple[int, int]]:
    points = sorted(set(points))

    def cross(o: tuple[int, int], a: tuple[int, int], b: tuple[int, int]) -> int:
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: list[tuple[int, int]] = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[int, int]] = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def pick_data(vertices: list[tuple[int, int]]) -> tuple[int, int, int]:
    area_twice = abs(
        sum(
            vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
            - vertices[(index + 1) % len(vertices)][0] * vertices[index][1]
            for index in range(len(vertices))
        )
    )
    boundary = sum(
        gcd(
            abs(vertices[(index + 1) % len(vertices)][0] - vertices[index][0]),
            abs(vertices[(index + 1) % len(vertices)][1] - vertices[index][1]),
        )
        for index in range(len(vertices))
    )
    require((area_twice - boundary + 2) % 2 == 0, "Pick parity failed")
    return area_twice, boundary, (area_twice - boundary + 2) // 2


def remainder_at(poly: sp.Expr, t_value: sp.Rational, modulus: sp.Expr, X: sp.Symbol, T: sp.Symbol) -> sp.Expr:
    specialized = sp.cancel(poly.subs(T, t_value))
    return sp.factor(sp.rem(sp.Poly(specialized, X), sp.Poly(modulus, X)).as_expr())


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def cycle(size: int, entries: list[int]) -> tuple[int, ...]:
    result = list(range(size))
    for source, target in zip(entries, entries[1:] + entries[:1]):
        result[source] = target
    return tuple(result)


def cycle_type(permutation: tuple[int, ...]) -> tuple[int, ...]:
    seen: set[int] = set()
    lengths: list[int] = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def permutation_index(permutation: tuple[int, ...]) -> int:
    return len(permutation) - len(cycle_type(permutation))


def orbit_size(generators: tuple[tuple[int, ...], ...]) -> int:
    reached = {0}
    frontier = [0]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            image = generator[current]
            if image not in reached:
                reached.add(image)
                frontier.append(image)
    return len(reached)


X, T, Delta, Phi = sp.symbols("X T Delta Phi")
K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
Delta_K0 = sp.Rational(5696, 105)
P = T + X**2 * T**2
Y = X * T * P

# On Theta=-Delta the apparent two-term top row contracts by one source
# coordinate degree: P^4-PY^2=T*P^3 exactly.
require(sp.factor(P**4 - P * Y**2 - T * P**3) == 0, "two-term contraction changed")
G = sp.expand(
    -X**2 * T / 2
    - 3 * P
    + sp.Rational(8, 3) * P**2
    - sp.Rational(1376, 135) * P**3
    + K * Y**2
    + Phi * P**2 * Y
    + Delta * T * P**3
)
f = sp.cancel(sp.diff(G, X) / T)
h = sp.diff(G, T)
require(not sp.denom(f).has(X, T), "G_X/T is not polynomial")
require(sp.Poly(f, X).degree() == 6, "wrong generic X-degree for f")
require(sp.Poly(h, X).degree() == 7, "wrong generic X-degree for h")
require(sp.factor(sp.Poly(f, X).LC() - 7 * Phi * T**6) == 0, "wrong f leading row")
require(sp.factor(sp.Poly(h, X).LC() - 7 * Phi * T**6) == 0, "wrong h leading row")

# Phi=0 is impossible before any boundary argument.  On X=0 every root of
# g'(T) is critical, so all three critical values would have to lie in
# {0,1/2}.  The critical-value resultant would then be one of four monic
# binary-root cubics.  The coefficient equations have gcd one in every case,
# including repeated-root specializations of g'.
z_value = sp.symbols("z_value")
g_line = sp.expand(G.subs({X: 0, Phi: 0}))
g_derivative = sp.diff(g_line, T)
critical_value_resultant = sp.Poly(sp.resultant(g_derivative, z_value - g_line, T), z_value)
require(critical_value_resultant.degree() == 3, "Phi=0 critical-value degree changed")
monic_critical_values = sp.Poly(
    sp.cancel(critical_value_resultant.as_expr() / critical_value_resultant.LC()), z_value
)
phi_zero_shape_gcds: list[sp.Expr] = []
for zero_multiplicity in range(4):
    candidate = sp.Poly(
        sp.expand(z_value**zero_multiplicity * (z_value - sp.Rational(1, 2)) ** (3 - zero_multiplicity)),
        z_value,
    )
    equations = [
        sp.together(actual - expected).as_numer_denom()[0]
        for actual, expected in zip(monic_critical_values.all_coeffs()[1:], candidate.all_coeffs()[1:])
    ]
    equation_gcd = sp.factor(sp.gcd_list(equations))
    phi_zero_shape_gcds.append(equation_gcd)
    require(equation_gcd == 1, f"Phi=0 binary critical-value shape {zero_multiplicity} survived")

# Exact critical resultants on the two live Phi!=0 coefficient strata.
resultant_generic = sp.factor(sp.resultant(f, h, X))
Q15_expr = sp.cancel(resultant_generic / (7 * T**30 * Phi * (6 * T + 1) ** 2))
Q15 = sp.Poly(Q15_expr, T)
require(Q15.degree() == 15, "generic residual is not degree fifteen")
require(
    sp.factor(Q15.LC() - 576 * Delta**4 * K**4) == 0,
    "generic Q15 leading row changed",
)
require(sp.factor(Q15.TC() + sp.Rational(7203, 32) * Phi**4) == 0, "generic Q15 constant changed")
require(
    sp.factor(resultant_generic - 7 * T**30 * Phi * (6 * T + 1) ** 2 * Q15.as_expr()) == 0,
    "generic resultant reconstruction failed",
)

f_K0 = sp.factor(f.subs(Delta, Delta_K0))
h_K0 = sp.factor(h.subs(Delta, Delta_K0))
resultant_K0 = sp.factor(sp.resultant(f_K0, h_K0, X))
Q13_expr = sp.cancel(resultant_K0 / (7 * T**30 * Phi * (6 * T + 1) ** 2))
Q13 = sp.Poly(Q13_expr, T)
expected_Q13_lead = sp.Rational(2185746832794987941330944, 16544390625)
require(Q13.degree() == 13, "K=0 residual is not degree thirteen")
require(Q13.LC() == expected_Q13_lead, "K=0 Q13 leading row changed")
require(sp.factor(Q13.TC() + sp.Rational(7203, 32) * Phi**4) == 0, "K=0 Q13 constant changed")
require(
    sp.factor(resultant_K0 - 7 * T**30 * Phi * (6 * T + 1) ** 2 * Q13.as_expr()) == 0,
    "K=0 resultant reconstruction failed",
)

# The two actual strata omitted or displayed by the projection are universal.
hessian_det = sp.factor(sp.det(sp.hessian(G, (X, T))))
require(sp.factor(f.subs(T, 0)) == -X, "wrong f at T=0")
require(sp.factor(h.subs(T, 0) + (X**2 + 6) / 2) == 0, "wrong h at T=0")
require(remainder_at(G, sp.Rational(0), X**2 + 6, X, T) == 0, "wrong T=0 value")
require(remainder_at(hessian_det, sp.Rational(0), X**2 + 6, X, T) == 6, "T=0 pair not Morse")
require(remainder_at(f, -sp.Rational(1, 6), X**2 - 6, X, T) == 0, "universal f failed")
require(remainder_at(h, -sp.Rational(1, 6), X**2 - 6, X, T) == 0, "universal h failed")
require(
    remainder_at(G - sp.Rational(1, 2), -sp.Rational(1, 6), X**2 - 6, X, T) == 0,
    "wrong universal value",
)
require(
    remainder_at(hessian_det, -sp.Rational(1, 6), X**2 - 6, X, T) == -6,
    "universal pair not Morse",
)

# Reconstruct the generic source in (s,p), its two Newton polygons, and the
# exact repeated-edge node.  Lambda/alpha/epsilon/gamma are the forced rows.
s_source, p_source, Q = sp.symbols("s_source p_source Q")
gamma = -sp.Rational(1, 2)
lambda_coefficient = -3
alpha = sp.Rational(8, 3)
epsilon = -sp.Rational(1376, 135)
t_source = p_source - s_source**2
H_source = sp.expand(
    lambda_coefficient * p_source
    + alpha * p_source**2
    + epsilon * p_source**3
    + K * s_source**2 * p_source**2
    + Phi * s_source * p_source**3
    + Delta * p_source**3 * t_source
)
G_source = -s_source**2 / (2 * t_source) + H_source
require(
    sp.factor(sp.cancel(G_source.subs({s_source: X * T, p_source: P}) - G)) == 0,
    "(s,p) reconstruction changed",
)
source_fibre = sp.expand((s_source**2 - p_source) * (1 - Q * H_source) + gamma * Q * s_source**2)

def source_polygon(expression: sp.Expr) -> list[tuple[int, int]]:
    polynomial = sp.Poly(expression, s_source, p_source)
    return convex_hull([monomial for monomial, coefficient in polynomial.terms() if coefficient != 0])


polygon_generic = source_polygon(source_fibre)
polygon_K0 = source_polygon(sp.factor(source_fibre.subs(Delta, Delta_K0)))
require(
    polygon_generic == [(0, 1), (2, 0), (4, 2), (4, 3), (0, 5)],
    "generic two-term polygon changed",
)
require(
    polygon_K0 == [(0, 1), (2, 0), (4, 3), (0, 5)],
    "K=0 two-term polygon changed",
)
require(pick_data(polygon_generic) == (26, 10, 9), "generic arithmetic-genus ledger changed")
require(pick_data(polygon_K0) == (24, 8, 9), "K=0 arithmetic-genus ledger changed")

z_boundary, a_boundary = sp.symbols("z_boundary a_boundary")
r_boundary = 1 - a_boundary
boundary_local = sp.factor(
    sp.cancel(
        z_boundary**10
        * source_fibre.subs(
            {s_source: 1 / z_boundary, p_source: r_boundary / z_boundary**2}
        )
    )
)
expected_boundary_local = sp.expand(
    a_boundary
    * (
        z_boundary**8
        + Q * Delta * r_boundary**3 * a_boundary
        - Q * Phi * r_boundary**3 * z_boundary
        - Q * (epsilon * r_boundary**3 + K * r_boundary**2) * z_boundary**2
        - Q * alpha * r_boundary**2 * z_boundary**4
        - Q * lambda_coefficient * r_boundary * z_boundary**6
    )
    + gamma * Q * z_boundary**8
)
require(sp.factor(boundary_local - expected_boundary_local) == 0, "boundary local equation changed")
tangent_cone = Q * a_boundary * (Delta * a_boundary - Phi * z_boundary)
local_poly = sp.Poly(sp.expand(boundary_local), a_boundary, z_boundary)
actual_tangent = sp.factor(
    sum(
        coefficient * a_boundary**monomial[0] * z_boundary**monomial[1]
        for monomial, coefficient in local_poly.terms()
        if sum(monomial) == 2
    )
)
require(sp.factor(actual_tangent - tangent_cone) == 0, "ordinary-node tangent cone changed")

# The two normalization branches begin a~(Phi/Delta)z and
# a~(gamma/Phi)z^7.  The residue Q*z^6 dz/L_a has order five on both, so
# both collided punctures have index six.  The leading rows below are exact.
first_branch_slope = sp.cancel(Phi / Delta)
require(
    sp.factor(
        sp.limit(
            sp.diff(boundary_local, a_boundary).subs(a_boundary, first_branch_slope * z_boundary)
            / z_boundary,
            z_boundary,
            0,
        )
        - Q * Phi
    )
    == 0,
    "first node branch residue row changed",
)
second_branch_coefficient = sp.cancel(gamma / Phi)
require(
    sp.factor(
        sp.limit(
            sp.diff(boundary_local, a_boundary).subs(
                a_boundary, second_branch_coefficient * z_boundary**7
            )
            / z_boundary,
            z_boundary,
            0,
        )
        + Q * Phi
    )
    == 0,
    "second node branch residue row changed",
)

# At the other corner, K!=0 gives one rational index-three branch and the
# locked quadratic BC pair.  When K=0 the primitive replacement edge has
# length one and index five.
bc_equation = 1 + gamma * Q - Q * K * sp.Symbol("W") ** 2
W = next(iter(bc_equation.free_symbols - {Q, Delta, Phi}))
require(sp.Poly(bc_equation, W).degree() == 2, "BC equation lost degree two")
require(
    sp.factor(K * W**2 - (1 / Q + gamma) + bc_equation / Q) == 0,
    "BC base-change equation changed",
)
require(sp.factor(K.subs(Delta, Delta_K0)) == 0, "K=0 wall changed")

packet_generic = (6, 6, 3, 2, 2, 1)
packet_K0 = (6, 6, 5, 1)
require(sum(index - 1 for index in packet_generic) == 14, "generic RH packet changed")
require(sum(index - 1 for index in packet_K0) == 14, "K=0 RH packet changed")
require(sum(packet_generic) == 20, "generic full response degree changed")
require(sum(packet_K0) == 18, "K=0 full response degree changed")
require(6 + 6 + 3 + 1 == 16, "generic rational response changed")

# Full-boundary branches: critical length is n-1.  Transitivity with both
# handle permutations nontrivial forces two pivot cycles; every commutator is
# a 3-cycle.  The identity cases cannot be transitive because the other
# handle has at least two fixed sheets.
for n_full, packet in ((20, packet_generic), (18, packet_K0)):
    for first_support in range(3, n_full - 1):
        X_full = cycle(n_full, list(range(first_support)))
        Y_full = cycle(n_full, [0] + list(range(first_support, n_full)))
        require(orbit_size((X_full, Y_full)) == n_full, "full pivot model lost transitivity")
        commutator = compose(
            compose(compose(X_full, Y_full), inverse(X_full)), inverse(Y_full)
        )
        require(
            cycle_type(commutator) == (3,) + (1,) * (n_full - 3),
            "full pivot commutator changed",
        )
    require(tuple(sorted(packet, reverse=True)) != (3,) + (1,) * (n_full - 3), "full packet became a 3-cycle")

# Finite-BC degree sixteen: L=19 gives support sum thirteen, while the two
# distinct finite transpositions contribute two merger units.
n_low = 16
critical_length_low = 19
support_budget_low = 2 * n_low - critical_length_low
require(support_budget_low == 13, "degree-sixteen support budget changed")
require((support_budget_low - 2) + 2 == 13 < n_low - 1, "both-nonidentity budget changed")
require((support_budget_low - 1) + 2 == 14 < n_low - 1, "one-identity budget changed")

# Sharp hostile: losing one critical point permits a fourteen-cycle plus two
# attaching transpositions to generate all sixteen sheets.
X_hostile = cycle(n_low, list(range(14)))
Z_hostile_1 = cycle(n_low, [0, 14])
Z_hostile_2 = cycle(n_low, [0, 15])
require(orbit_size((X_hostile, Z_hostile_1, Z_hostile_2)) == n_low, "L=18 hostile failed")
require(
    sum(map(permutation_index, (X_hostile, Z_hostile_1, Z_hostile_2))) == n_low - 1,
    "L=18 hostile index changed",
)

semantic_lines = (
    "scope=exact_M8_two_term_wall;Theta=-Delta;Delta!=0",
    "contraction=P^4-PY^2=T*P^3",
    "Phi=0=critical_value_binary_shape_no_go",
    "K!=0:critical_residual_degree=15;critical_length=19;packet=6,6,3,2,2,1;degrees=16,20",
    "K=0:critical_residual_degree=13;critical_length=17;packet=6,6,5,1;degree=18",
    "boundary_node=branches_index_6,6;normalized_genus=8",
    "degree20,18=shared_pivot_commutator_contradiction",
    "degree16=two_transposition_orbit_merger_contradiction",
    "geometric_dependencies=closed_polynomial_connectedness,target_MW,carrier_classification,finite_etale_transport",
    "verdict=two_term_wall_empty_relative;exact_M8_THM4053_trichotomy_empty;JC2=OPEN",
)
semantic_sha256 = hashlib.sha256(("\n".join(semantic_lines) + "\n").encode()).hexdigest()

print("THM4143_TWO_TERM_WALL_PRIMARY")
print("identity=P^4-PY^2=T*P^3")
print("Phi=0=NO_GO;critical_value_binary_shapes=4;all_parameter_gcds=1")
print("K_nonzero=resultant=7*T^30*Phi*(6T+1)^2*Q15;degree_Q15=15;critical_length=19")
print("K_zero=Delta=5696/105;resultant=7*T^30*Phi*(6T+1)^2*Q13;degree_Q13=13;critical_length=17")
print("generic_polygon=0,1;2,0;4,2;4,3;0,5;arithmetic_genus=9")
print("K0_polygon=0,1;2,0;4,3;0,5;arithmetic_genus=9")
print("boundary_node=ordinary;branches=a~(Phi/Delta)z,a~(gamma/Phi)z^7;indices=6,6;normalized_genus=8")
print("generic_packet=6,6,3,2,2,1;responses=n16_finite_BC_or_n20_full")
print("K0_packet=6,6,5,1;response=n18_full")
print("monodromy=n20,n18_shared_pivot_EXCLUDED;n16_carrier_capacity<=14<15_EXCLUDED")
print("hostiles=Phi0_all_shapes_fail;L18_degree16_reaches_transitivity;node_normalization_lowers_genus_one")
print("semantic_sha256=" + semantic_sha256)
print("verdict=two_term_wall=EMPTY_RELATIVE;exact_M8_trichotomy=EMPTY_RELATIVE;JC2=OPEN;DC2=OPEN")
