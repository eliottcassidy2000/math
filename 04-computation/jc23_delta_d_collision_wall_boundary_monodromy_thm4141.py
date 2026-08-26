#!/usr/bin/env python3
"""Primary exact audit for THM-4141's Delta_D wall exclusion.

The calculation rebuilds the delta-only generic source polynomial, resolves
the repeated outer-edge root in a toric chart, computes its residue order,
and freezes the two resulting monodromy budgets.  It imports no repository
code and uses no Python ``assert``, so optimized mode retains every gate.
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


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def cycle(permutation_size: int, entries: list[int]) -> tuple[int, ...]:
    result = list(range(permutation_size))
    if entries:
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


s, p, r, Q = sp.symbols("s p r Q")
Delta = sp.Rational(5696, 15) / (6 * r**2 + 7)
K = sp.cancel(Delta * r**2)
Phi = sp.cancel(2 * Delta * r)

H = (
    -3 * p
    + sp.Rational(8, 3) * p**2
    - sp.Rational(1376, 135) * p**3
    + K * s**2 * p**2
    + Phi * s * p**3
    + Delta * p**4
)
F = sp.expand((s**2 - p) * (1 - Q * H) - sp.Rational(1, 2) * Q * s**2)

require(sp.factor(Phi**2 - 4 * K * Delta) == 0, "Delta_D wall changed")
require(sp.factor(K / Delta - r**2) == 0, "wall parameter changed")

# F is linear and primitive in Q.  Its two Q-coefficients are coprime for
# every coefficient specialization, since the Q-row reduces to -s^2/2
# modulo s^2-p.  Gauss therefore excludes a factorization over C[Q,s,p].
# Geometric connectedness additionally uses the finite-critical-locus/Stein
# gate proved in the theorem; irreducibility of the total Q-equation alone
# would not suffice after an algebraic base extension.
F_as_Q = sp.Poly(F, Q)
require(F_as_Q.degree() == 1, "generic-fibre equation is not linear in Q")
q_constant, q_linear = F_as_Q.nth(0), F_as_Q.nth(1)
require(sp.factor(q_constant - (s**2 - p)) == 0, "Q-constant row changed")
require(
    sp.factor(sp.rem(sp.Poly(q_linear, p), sp.Poly(s**2 - p, p)).as_expr() + s**2 / 2) == 0,
    "Q-linear row has wrong reduction modulo s^2-p",
)
require(sp.gcd(sp.Poly(q_constant, s, p), sp.Poly(q_linear, s, p)).total_degree() == 0, "generic fibre lost primitivity")

# The generic Newton polygon is unchanged by the optional interior
# cancellation K=epsilon.  Its repeated CD edge is handled below.
support = [(i, j) for (i, j), coefficient in sp.Poly(F, s, p).terms() if coefficient != 0]
polygon = convex_hull(support)
expected_polygon = [(0, 1), (2, 0), (4, 2), (2, 4), (0, 5)]
require(polygon == expected_polygon, "delta-only Newton polygon changed")
require(
    convex_hull([point for point in support if point != (2, 3)]) == expected_polygon,
    "optional K=epsilon cancellation changed the polygon",
)
area_twice, boundary_steps, interior_points = pick_data(polygon)
require((area_twice, boundary_steps, interior_points) == (24, 10, 8), "genus ledger changed")

# (primitive inward normal, supporting constant, lattice length,
#  ordinary residue order) for the five counterclockwise polygon edges.
edge_ledger = (
    ((1, 2), 2, 1, 0),
    ((-1, 1), -2, 2, 1),
    ((-1, -1), -6, 2, 3),
    ((-1, -2), -10, 1, 6),
    ((1, 0), 0, 4, 0),
)
for index, ((nx, ny), height, length, residue_order) in enumerate(edge_ledger):
    first = polygon[index]
    second = polygon[(index + 1) % len(polygon)]
    require(nx * first[0] + ny * first[1] == height, f"edge {index} first endpoint failed")
    require(nx * second[0] + ny * second[1] == height, f"edge {index} second endpoint failed")
    distance = nx + ny - height
    require(distance - 1 == residue_order, f"edge {index} residue order failed")
    require(
        gcd(abs(second[0] - first[0]), abs(second[1] - first[1])) == length,
        f"edge {index} length failed",
    )

# The CD edge polynomial has one double nonzero root.  The BC roots instead
# form one quadratic orbit over C(r)(Q), because their equation has odd
# Q-valuation.  The DE and AB roots are rational.
w = sp.symbols("w")
cd_edge = sp.factor(K + Phi * w + Delta * w**2)
require(sp.factor(cd_edge - Delta * (w + r) ** 2) == 0, "CD repeated edge changed")
bc_edge = 1 - Q / 2 - Q * K * w**2
require(sp.Poly(bc_edge, w).degree() == 2, "BC edge is not quadratic")
require(
    sp.factor(sp.discriminant(bc_edge, w) - 4 * Q * K * (1 - Q / 2)) == 0,
    "BC discriminant changed",
)
require(
    sp.factor(K * w**2 - (1 / Q - sp.Rational(1, 2))) == sp.factor(-bc_edge / Q),
    "BC quadratic base change changed",
)

# The length-four s=0 edge gives four nonzero affine points generically.
# A rational live specialization proves that its discriminant is not the
# zero function; the complete generic smoothness statement is inherited.
s_zero_edge = sp.Poly(1 - Q * H.subs(s, 0), p)
require(s_zero_edge.degree() == 4, "s=0 edge lost degree four")
require(
    sp.discriminant(s_zero_edge, p).subs({r: 1, Q: 1}) != 0,
    "s=0 generic squarefreeness control failed",
)

# Resolve the repeated root T=p/s=-r.  This is an exact chart of the full
# curve, not merely of its leading face.
u, v = sp.symbols("u v")
local = sp.cancel(u**6 * F.subs({s: 1 / u, p: (v - r) / u}))
require(not sp.denom(local).has(u, v), "local chart is not polynomial")
local = sp.expand(local)
require(
    sp.factor(local.subs(u, 0) + Delta * Q * v**2 * (v - r) ** 2) == 0,
    "boundary edge specialization changed",
)
linear_u = sp.factor(sp.diff(local, u).subs({u: 0, v: 0}))
quadratic_v = sp.factor(sp.diff(local, v, 2).subs({u: 0, v: 0}) / 2)
require(linear_u == -sp.Rational(1376, 135) * Q * r**3, "forced transverse row changed")
require(sp.factor(quadratic_v + Delta * Q * r**2) == 0, "tangent quadratic row changed")

# Since linear_u is nonzero on the live wall, the boundary point is smooth
# and u=c*v^2+O(v^3).  The pullback residue is
#   Q ds/F_p = -Q*u^3*du/local_v,
# so its first coefficient is nonzero at order v^6.
tangent_coefficient = sp.factor(-quadratic_v / linear_u)
expected_tangent = -sp.Rational(135, 1376) * Delta / r
require(sp.factor(tangent_coefficient - expected_tangent) == 0, "tangent coefficient changed")
residue_lead = sp.factor(-Q * tangent_coefficient**4 / quadratic_v)
require(
    sp.factor(residue_lead - tangent_coefficient**4 / (Delta * r**2)) == 0,
    "residue leading coefficient changed",
)
# The coefficient -1376/135 is load-bearing.  Removing the p^3 row makes
# the transverse derivative vanish; this is a hostile model, not a live
# specialization.
H_without_cubic = H + sp.Rational(1376, 135) * p**3
F_without_cubic = sp.expand(
    (s**2 - p) * (1 - Q * H_without_cubic) - sp.Rational(1, 2) * Q * s**2
)
local_without_cubic = sp.cancel(
    u**6 * F_without_cubic.subs({s: 1 / u, p: (v - r) / u})
)
require(
    sp.factor(sp.diff(local_without_cubic, u).subs({u: 0, v: 0})) == 0,
    "deleted-cubic hostile did not lose transversality",
)

infinity_packet = (7, 7, 2, 2, 1)
require(sum(index - 1 for index in infinity_packet) == 14, "ramification budget changed")
require(sum(infinity_packet[index] for index in (0, 1, 4)) == 15, "rational packet changed")
require(sum(infinity_packet) == 19, "full packet degree changed")

# High branch n=19.  If both handle generators are nonidentity, transitivity
# and the n+1 support ceiling force two cycles meeting in one pivot; every
# such commutator is a 3-cycle.  The identity case has trivial commutator.
n_high = 19
identity_high = tuple(range(n_high))
for x_support in range(2, n_high):
    X = cycle(n_high, list(range(x_support)))
    Y = cycle(n_high, [0] + list(range(x_support, n_high)))
    require(orbit_size((X, Y)) == n_high, "shared-pivot model lost transitivity")
    commutator = compose(compose(compose(X, Y), inverse(X)), inverse(Y))
    require(cycle_type(commutator) == (3,) + (1,) * 16, "shared-pivot commutator changed")
require(cycle_type(identity_high) != tuple(sorted(infinity_packet, reverse=True)), "identity hostile failed")
require((3,) + (1,) * 16 != tuple(sorted(infinity_packet, reverse=True)), "high packet became a 3-cycle")

# Low branch n=15.  Critical length 18 gives handle support at most n-3,
# hence handle index at most n-4.  The quadratic BC carrier contributes two
# more units, still one short of transitivity.
n_low = 15
critical_length = 18
support_budget = 2 * n_low - critical_length
handle_index_budget = support_budget - 1
carrier_index = 2
require(support_budget == n_low - 3, "low support budget changed")
require(handle_index_budget + carrier_index == n_low - 2, "low defect budget changed")
require(handle_index_budget + carrier_index < n_low - 1, "low branch no longer contradicts transitivity")

# Two sharp controls show that both one-unit margins are load-bearing.
X_actual = cycle(n_low, list(range(12)))
Z_actual_1 = cycle(n_low, [0, 12])
Z_actual_2 = cycle(n_low, [0, 13])
require(orbit_size((X_actual, Z_actual_1, Z_actual_2)) == 14, "actual sharp control changed")
require(
    sum(map(permutation_index, (X_actual, Z_actual_1, Z_actual_2))) == n_low - 2,
    "actual index sharpness changed",
)

X_one_less_fixed = cycle(n_low, list(range(13)))
Z_attach_1 = cycle(n_low, [0, 13])
Z_attach_2 = cycle(n_low, [0, 14])
require(orbit_size((X_one_less_fixed, Z_attach_1, Z_attach_2)) == n_low, "fixed-sheet hostile failed")
require(
    sum(map(permutation_index, (X_one_less_fixed, Z_attach_1, Z_attach_2))) == n_low - 1,
    "fixed-sheet hostile index changed",
)

Z_third = cycle(n_low, [0, 14])
require(
    orbit_size((X_actual, Z_actual_1, Z_actual_2, Z_third)) == n_low,
    "carrier-three hostile failed",
)
require(
    sum(map(permutation_index, (X_actual, Z_actual_1, Z_actual_2, Z_third))) == n_low - 1,
    "carrier-three hostile index changed",
)

semantic_lines = (
    "scope=delta_only_exact_M8_DeltaD_wall",
    "polygon=(0,1),(2,0),(4,2),(2,4),(0,5);genus=8",
    "collision_chart=L(0,v)=-Delta*Q*v^2*(v-r)^2;L_u(0,0)=-1376*Q*r^3/135",
    "collision_branch=u~-(135*Delta/(1376*r))*v^2;residue_order=6;index=7",
    "generic_source_Q_equation=primitive_and_irreducible;geometric_connectedness=finite_critical_Stein_gate",
    "infinity_packet=7,7,2,2,1;defect=14",
    "response_degrees=15,19",
    "critical_length=18",
    "degree19=shared_pivot_commutator_contradiction",
    "degree15=orbit_merger_index_at_most_13<14",
    "verdict=DeltaD_wall_empty_relative_to_inherited_Keller_and_transport_gates;JC2=OPEN",
)
semantic_sha256 = hashlib.sha256(("\n".join(semantic_lines) + "\n").encode()).hexdigest()

print("THM4141_DELTA_D_BOUNDARY_MONODROMY_PRIMARY")
print("polygon=(0,1),(2,0),(4,2),(2,4),(0,5);pick=area2:24,boundary:10,genus:8")
print("collision=L(0,v)=-Delta*Q*v^2*(v-r)^2;L_u=-1376*Q*r^3/135;smooth=YES")
print("collision_branch=u~-(135*Delta/(1376*r))*v^2;residue_order=6;index=7")
print("infinity_packet=7,7,2,2,1;ramification_defect=14")
print("affine_s0_packet=1,1,1,1")
print("response_degrees=15,19;critical_length=18")
print("degree19=boundary_cycle_type_7,7,2,2,1!=identity_or_3cycle;EXCLUDED")
print("degree15=handle_index<=11;carrier_index=2;total<=13<14;EXCLUDED")
print("hostiles=delete_cubic_loses_transversality;fixed_sum_n+2_allows_transitivity;carrier_index_3_allows_transitivity")
print("semantic_sha256=" + semantic_sha256)
print("verdict=DeltaD_wall=EMPTY_RELATIVE;JC2=OPEN;DC2=OPEN")
