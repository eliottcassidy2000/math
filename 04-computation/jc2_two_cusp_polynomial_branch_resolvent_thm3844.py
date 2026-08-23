#!/usr/bin/env python3
"""Exact companion for THM-3844's two-cusp resolvent design gate."""

from __future__ import annotations

import ast
import hashlib
import json
from itertools import combinations, product
from math import gcd
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


t, X, Y, W, Z, T = sp.symbols("t X Y W Z T")

# Polynomial normalization candidate.
X_t = t**2 * (2 * t - 3)
Y_t = t**3 * (3 * t - 4)
delta = (
    27 * X**4
    + 16 * X**3
    + 6 * X**2 * Y
    + 48 * X * Y**2
    - 16 * Y**3
    + 27 * Y**2
)
zero(delta.subs({X: X_t, Y: Y_t}), "quartic vanishes on normalization")
gate(sp.factor(delta) == delta, "quartic is irreducible over Q")
gate(sp.Poly(delta, X, Y).total_degree() == 4, "quartic total degree")

# Elimination and rational inverse prove that k[t] is the normalization.
resultant = sp.factor(
    sp.resultant(2 * t**3 - 3 * t**2 - X, 3 * t**4 - 4 * t**3 - Y, t)
)
zero(resultant - delta, "normalization eliminant")
inverse_numerator = -X**2 + 4 * X * Y + 3 * Y
inverse_denominator = 6 * X**2 + 4 * X + 2 * Y
zero(
    inverse_numerator.subs({X: X_t, Y: Y_t})
    - t * inverse_denominator.subs({X: X_t, Y: Y_t}),
    "birational inverse",
)
zero(t**3 - sp.Rational(3, 2) * t**2 - X_t / 2,
     "normalization parameter is integral")

# Normalization-conductor sidecar.  The two cusp addresses contribute order
# two and the two node branches order one.  The three displayed contractions
# also verify directly that c*k[t] lies in the branch ring, since 1,t,t^2
# generate k[t] over it by the integral cubic above.
node_factor = 2 * t**2 - 2 * t - 1
conductor = t**2 * (t - 1) ** 2 * node_factor
zero(
    conductor - (3 * X_t**2 + 2 * X_t + Y_t) / 6,
    "conductor generator contracts to the branch ring",
)
zero(
    t * conductor + (X_t**2 - 4 * X_t * Y_t - 3 * Y_t) / 12,
    "t times conductor contracts to the branch ring",
)
zero(
    t**2 * conductor + (X_t - Y_t) * (X_t + 2 * Y_t) / 9,
    "t squared times conductor contracts to the branch ring",
)

# Exact singular locus.
delta_X = sp.diff(delta, X)
delta_Y = sp.diff(delta, Y)
singular_groebner = sp.groebner([delta, delta_X, delta_Y], Y, X, order="lex")
expected_singular_groebner = sp.groebner(
    [Y - 8 * X**4 - 16 * X**3 - 7 * X**2,
     X**2 * (X + 1) ** 2 * (2 * X + 1)],
    Y,
    X,
    order="lex",
)
gate(singular_groebner == expected_singular_groebner, "complete singular locus")
for point in ((0, 0), (-1, -1), (sp.Rational(-1, 2), sp.Rational(1, 4))):
    gate(
        delta.subs({X: point[0], Y: point[1]}) == 0
        and delta_X.subs({X: point[0], Y: point[1]}) == 0
        and delta_Y.subs({X: point[0], Y: point[1]}) == 0,
        f"singular point {point}",
    )

# The two cusps and the node are visible directly on the normalization.
s = sp.symbols("s")
zero(sp.expand(X_t.subs(t, s)) - (-3 * s**2 + 2 * s**3),
     "first cusp X expansion")
zero(sp.expand(Y_t.subs(t, s)) - (-4 * s**3 + 3 * s**4),
     "first cusp tangent expansion")
zero(
    sp.expand(X_t.subs(t, 1 + s) + 1) - (3 * s**2 + 2 * s**3),
    "second cusp X expansion",
)
zero(
    sp.expand((Y_t - 2 * X_t - 1).subs(t, 1 + s))
    - (4 * s**3 + 3 * s**4),
    "second cusp tangent expansion",
)
node_polynomial = 2 * t**2 - 2 * t - 1
zero(sp.rem(X_t + sp.Rational(1, 2), node_polynomial, t),
     "node has two X preimages")
zero(sp.rem(Y_t - sp.Rational(1, 4), node_polynomial, t),
     "node has two Y preimages")
gate(sp.discriminant(node_polynomial, t) == 12, "node preimages are distinct")
slope = sp.cancel(sp.diff(Y_t, t) / sp.diff(X_t, t))
zero(slope - 2 * t, "normalization tangent slope")
gate(
    sp.rem(2 * t, node_polynomial, t) != sp.rem(2 * (1 - t), node_polynomial, t),
    "node tangent slopes differ",
)

# One projective place at infinity, and it is smooth.
delta_h = sp.expand(Z**4 * delta.subs({X: X / Z, Y: Y / Z}))
zero(delta_h.subs(Z, 0) - 27 * X**4, "unique projective infinity support")
infinity_point = {X: 0, Y: 1, Z: 0}
gate(delta_h.subs(infinity_point) == 0, "infinity point lies on closure")
gate(sp.diff(delta_h, Z).subs(infinity_point) == -16,
     "projective infinity point is smooth")

# Tangent contacts at the two cusps and their sole residual points.
zero(Y_t - t**3 * (3 * t - 4), "first cusp tangent pullback")
zero(Y_t - 2 * X_t - 1 - (t - 1) ** 3 * (3 * t + 1),
     "second cusp tangent pullback")


# Weak degree-two del Pezzo lattice.  Coordinates are H,E1,...,E7 with
# intersection diag(1,-1,...,-1).  The finite exact enumeration checks that
# every 2A2+A1 root configuration orthogonal to one boundary (-1)-curve has
# the same affine quotient torsion order three.
def lattice_dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return left[0] * right[0] - sum(a * b for a, b in zip(left[1:], right[1:]))


def lattice_neg(vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(-entry for entry in vector)


def canonical_sign(vector: tuple[int, ...]) -> tuple[int, ...]:
    return min(vector, lattice_neg(vector))


roots: set[tuple[int, ...]] = set()
for i in range(7):
    for j in range(7):
        if i != j:
            vector = [0] * 8
            vector[i + 1] = 1
            vector[j + 1] = -1
            roots.add(tuple(vector))
for indices in combinations(range(7), 3):
    vector = [0] * 8
    vector[0] = 1
    for index in indices:
        vector[index + 1] = -1
    roots.add(tuple(vector))
    roots.add(lattice_neg(tuple(vector)))
for omitted in range(7):
    vector = [2] + [-1] * 7
    vector[omitted + 1] = 0
    roots.add(tuple(vector))
    roots.add(lattice_neg(tuple(vector)))
gate(len(roots) == 126, "complete E7 root list")
gate(
    all(lattice_dot(root, root) == -2 for root in roots),
    "every listed E7 root has square minus two",
)

boundary_plus = (0, 0, 0, 0, 0, 0, 0, 1)
anticanonical = (3, -1, -1, -1, -1, -1, -1, -1)
boundary_minus = tuple(
    anticanonical[index] - boundary_plus[index] for index in range(8)
)
gate(lattice_dot(boundary_plus, boundary_plus) == -1, "first boundary square")
gate(lattice_dot(boundary_minus, boundary_minus) == -1, "second boundary square")
gate(lattice_dot(boundary_plus, boundary_minus) == 2, "hyperflex boundary contact")
gate(
    sp.Matrix([boundary_plus, boundary_minus]).rank() == 2,
    "the two boundary classes are independent",
)

orthogonal_roots = sorted(
    root for root in roots if lattice_dot(root, boundary_plus) == 0
)
gate(len(orthogonal_roots) == 72, "boundary orthogonal E6 root list")
root_lookup = set(orthogonal_roots)

# Deduplicate each A2 by its full six-root subsystem.
a2_systems: dict[
    tuple[tuple[int, ...], ...], tuple[tuple[int, ...], tuple[int, ...]]
] = {}
for left, right in combinations(orthogonal_roots, 2):
    if lattice_dot(left, right) != 1:
        continue
    root_sum = tuple(left[index] + right[index] for index in range(8))
    six_roots = tuple(
        sorted({left, right, lattice_neg(left), lattice_neg(right),
                root_sum, lattice_neg(root_sum)})
    )
    gate(all(root in root_lookup for root in six_roots), "A2 root closure")
    a2_systems[six_roots] = (left, right)
gate(len(a2_systems) == 120, "complete boundary-orthogonal A2 list")

# A configuration key contains two unordered orthogonal A2 systems and one
# orthogonal A1 sign-pair.
configurations: dict[
    tuple[object, ...],
    tuple[
        tuple[tuple[int, ...], tuple[int, ...]],
        tuple[tuple[int, ...], tuple[int, ...]],
        tuple[int, ...],
    ],
] = {}
a2_items = list(a2_systems.items())
for first_index, (first_key, first_pair) in enumerate(a2_items):
    for second_key, second_pair in a2_items[first_index + 1:]:
        if any(
            lattice_dot(left, right) != 0
            for left in first_key
            for right in second_key
        ):
            continue
        for a1_root in orthogonal_roots:
            if a1_root != canonical_sign(a1_root):
                continue
            if any(
                lattice_dot(a1_root, root) != 0
                for root in first_key + second_key
            ):
                continue
            key = (tuple(sorted((first_key, second_key))), a1_root)
            configurations.setdefault(
                key, (first_pair, second_pair, a1_root)
            )
gate(len(configurations) == 360, "complete 2A2+A1 boundary configurations")

torsion_orders: dict[int, int] = {}
local_torsion_patterns: dict[tuple[bool, bool, int], int] = {}
for first_pair, second_pair, a1_root in configurations.values():
    relation_rows = [
        *first_pair, *second_pair, a1_root, boundary_plus, boundary_minus
    ]
    matrix = sp.Matrix(relation_rows)
    minor_gcd = 0
    for omitted_column in range(8):
        columns = [index for index in range(8) if index != omitted_column]
        determinant = abs(int(matrix[:, columns].det(method="domain-ge")))
        minor_gcd = gcd(minor_gcd, determinant)
    gate(minor_gcd != 0, "configuration quotient has rank one")
    torsion_orders[minor_gcd] = torsion_orders.get(minor_gcd, 0) + 1

    # The unique nonzero coefficient line in the left kernel modulo three
    # constructs the saturated-lattice lift D with 3D in the relation
    # lattice.  Its two A2 discriminant characters are both primitive; its
    # A1 character vanishes.  This rules out a marking-dependent local-class
    # accident, rather than checking only the convenient model below.
    coefficient_lines = [
        coefficients
        for coefficients in product(range(3), repeat=7)
        if coefficients != (0,) * 7
        and all(
            sum(
                coefficients[row] * relation_rows[row][column]
                for row in range(7)
            ) % 3 == 0
            for column in range(8)
        )
    ]
    gate(len(coefficient_lines) == 2, "unique mod-three saturation line")
    torsion_numerator = tuple(
        sum(
            coefficient_lines[0][row] * relation_rows[row][column]
            for row in range(7)
        )
        for column in range(8)
    )
    gate(
        all(entry % 3 == 0 for entry in torsion_numerator),
        "integral saturated torsion lift",
    )
    torsion_lift = tuple(entry // 3 for entry in torsion_numerator)
    local_pattern = (
        (
            lattice_dot(torsion_lift, first_pair[0])
            + 2 * lattice_dot(torsion_lift, first_pair[1])
        ) % 3 != 0,
        (
            lattice_dot(torsion_lift, second_pair[0])
            + 2 * lattice_dot(torsion_lift, second_pair[1])
        ) % 3 != 0,
        lattice_dot(torsion_lift, a1_root) % 2,
    )
    local_torsion_patterns[local_pattern] = (
        local_torsion_patterns.get(local_pattern, 0) + 1
    )
gate(torsion_orders == {3: 360}, "every affine lattice quotient has torsion Z/3")
gate(
    local_torsion_patterns == {(True, True, 0): 360},
    "the global three-class is primitive at both A2 points and trivial at A1",
)

# One explicit marking makes the free and torsion generators transparent.
r1 = (0, 1, -1, 0, 0, 0, 0, 0)
r2 = (0, 0, 1, -1, 0, 0, 0, 0)
r3 = (0, 0, 0, 0, 1, -1, 0, 0)
r4 = (0, 0, 0, 0, 0, 1, -1, 0)
r5 = (1, -1, -1, -1, 0, 0, 0, 0)
explicit_matrix = sp.Matrix([r1, r2, r3, r4, r5, boundary_plus, boundary_minus])
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

smith = smith_normal_form(explicit_matrix, domain=ZZ)
gate(
    [smith[index, index] for index in range(7)] == [1, 1, 1, 1, 1, 1, 3],
    "explicit affine class quotient is Z plus Z/3",
)
torsion_vector = (0, 2, 0, 0, -1, 0, 0, 0)  # 2E1-E4
torsion_relation_coefficients = (4, 2, -2, -1, -3, 2)
torsion_relation = tuple(
    boundary_minus[index]
    + sum(
        coefficient * root[index]
        for coefficient, root in zip(
            torsion_relation_coefficients,
            (r1, r2, r3, r4, r5, boundary_plus),
        )
    )
    for index in range(8)
)
gate(
    torsion_relation == tuple(3 * entry for entry in torsion_vector),
    "explicit order-three divisor relation",
)
canonical = lattice_neg(anticanonical)
deck_torsion = tuple(
    lattice_dot(torsion_vector, canonical) * canonical[index]
    - torsion_vector[index]
    for index in range(8)
)
gate(
    tuple(deck_torsion[index] + torsion_vector[index] for index in range(8))
    == anticanonical,
    "Geiser involution acts by inversion after deleting both boundaries",
)
gate(lattice_dot(torsion_vector, r1) == -2,
     "torsion class is nontrivial at the first A2 cusp")
gate(lattice_dot(torsion_vector, r3) == 1,
     "torsion class is nontrivial at the second A2 cusp")
gate(lattice_dot(torsion_vector, r5) % 2 == 0,
     "torsion class is trivial at the A1 node")
first_local_pairing = (
    lattice_dot(torsion_vector, r1),
    lattice_dot(torsion_vector, r2),
)
# Orient the second A2 oppositely, by the ordered basis (r4,r3).  For an A2
# Cartan block the quotient character is (a,b) |-> a+2b mod 3.
second_local_pairing = (
    lattice_dot(torsion_vector, r4),
    lattice_dot(torsion_vector, r3),
)
gate(
    (first_local_pairing[0] + 2 * first_local_pairing[1]) % 3 == 1,
    "first A2 local class is primitive",
)
gate(
    (second_local_pairing[0] + 2 * second_local_pairing[1]) % 3 == 2,
    "second A2 local class is the oppositely oriented primitive class",
)


# The apparent mixed 3-class is exactly a pullback of the depressed-cubic
# cusp, so its S3 cubic is globally monogenic.
p_cubic = X - Y
q_cubic = (X**2 + Y) / 2
h_parameter = t * (t - 1)
zero(p_cubic.subs({X: X_t, Y: Y_t}) + 3 * h_parameter**2,
     "cusp pullback p coordinate")
zero(q_cubic.subs({X: X_t, Y: Y_t}) - 2 * h_parameter**3,
     "cusp pullback q coordinate")
zero(-4 * p_cubic**3 - 27 * q_cubic**2 + delta / 4,
     "quartic is a depressed-cusp pullback")
coefficient_jacobian = sp.det(
    sp.Matrix(
        [
            [sp.diff(p_cubic, X), sp.diff(p_cubic, Y)],
            [sp.diff(q_cubic, X), sp.diff(q_cubic, Y)],
        ]
    )
)
zero(coefficient_jacobian - (X + sp.Rational(1, 2)),
     "coefficient map Jacobian")
gate(
    coefficient_jacobian.subs({X: sp.Rational(-1, 2), Y: sp.Rational(1, 4)}) == 0,
    "coefficient map ramifies at the node",
)

monogenic = T**3 + p_cubic * T + q_cubic
monogenic_discriminant = sp.discriminant(monogenic, T)
zero(monogenic_discriminant + delta / 4, "monogenic discriminant")

# Over the quadratic resolvent W^2=delta, Cardano gives the cyclic cubic
# Kummer layer explicitly.  The product identity shows that its conjugate is
# exactly the cube paired with u^3 under uv=-p/3.
cardano_plus = -q_cubic / 2 + W / (12 * sp.sqrt(3))
cardano_minus = -q_cubic / 2 - W / (12 * sp.sqrt(3))
zero(
    sp.expand(cardano_plus * cardano_minus + p_cubic**3 / 27).subs(W**2, delta),
    "Cardano Kummer norm identity",
)

# Irreducibility hostile: at X=0 the monic cubic is Eisenstein at Y.
specialized_twice = sp.expand(2 * monogenic.subs(X, 0))
gate(sp.Poly(specialized_twice, T).LC() == 2, "specialized leading unit")
zero(sp.Poly(specialized_twice, T).coeff_monomial(T) + 2 * Y,
     "specialized linear coefficient is Y-divisible")
zero(sp.Poly(specialized_twice, T).coeff_monomial(1) - Y,
     "specialized constant has exact Y-order one")

# Different/companion identity: deleting the ramified sheet makes a
# nonconstant polynomial a unit, so the one-place branch is not a
# constant-unit Keller source.
different = 3 * T**2 + p_cubic
companion = -3 * T**2 - 4 * p_cubic
zero(
    sp.rem(different**2 * companion - monogenic_discriminant, monogenic, T),
    "different times companion discriminant identity",
)
gate(sp.Poly(different, T).total_degree() == 2, "different is nonconstant")


semantic = {
    "branch": "delta=27X4+16X3+6X2Y+48XY2-16Y3+27Y2",
    "normalization": "X=t2(2t-3);Y=t3(3t-4);conductor=t2(t-1)2(2t2-2t-1);A1;one projective infinity",
    "singularities": "A2 cusps at (0,0),(-1,-1);A1 node at (-1/2,1/4)",
    "resolvent": "W2=delta;projective ADE=2A2+A1;affine Cl=Z plus Z/3",
    "torsion": "all 360 marked root-boundary configurations have torsion order 3",
    "cusp_pullback": "p=X-Y;q=(X2+Y)/2;delta=-4*disc(p,q)",
    "monogenic": "T3+pT+q;Cardano cyclic layer;different=3T2+p becomes a nonconstant unit after deletion",
    "scope": "exact design hostile;no Keller atlas or JC counterexample",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate")
print("branch=quartic;normalization=A1;infinity_places=1")
print("finite_singularities=2A2_plus_A1")
print("quadratic_resolvent_affine_class_group=Z_plus_Z3")
print("root_boundary_configurations=360;torsion_order=3_uniform")
print("cusp_pullback=p_X_minus_Y;q_half_X2_plus_Y")
print("cubic=globally_monogenic;deleted_different_nonconstant_unit")
print("scope=design_hostile_no_plane_atlas_no_JC_counterexample")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
