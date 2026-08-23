"""Exact probe for the reserved THM-3851 tricuspidal tradeoff.

This is a research artifact, not a theorem promotion.  It freezes three
independent packets:

1. the standard tricuspidal quartic has a split bitangent line but no line
   whose pullback to the normalization is supported at one point; and
2. the degree-two del Pezzo root/boundary lattice for 3A2 has affine quotient
   (Z/3)^2, with the three local A2 characters subject to one full-support
   linear relation; and
3. the quartic is the discriminant of the reciprocal cubic, whose root sheet
   is A1 x G_m and whose ordered-root cover is the A2 torus (G_m)^2.

The geometric bridge from a double quartic to the weak degree-two del Pezzo
resolution/class-group quotient remains a human input, exactly as in
THM-3844.  Run with --verify-frozen PATH to compare the semantic transcript
against an immutable output file.
"""

from __future__ import annotations

import hashlib
import sys
from collections import Counter
from itertools import combinations, product
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

sys.stdout.reconfigure(newline="\n")


GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"gate failed: {label}")


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.expand(expression) == 0, label)


def lattice_dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return left[0] * right[0] - sum(
        a * b for a, b in zip(left[1:], right[1:])
    )


def lattice_neg(vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(-entry for entry in vector)


# ---------------------------------------------------------------------------
# I. The standard tricuspidal quartic and its normalization.
# ---------------------------------------------------------------------------
s, t, U, V, Z = sp.symbols("s t U V Z")

delta = U**2 * V**2 - 4 * (U**3 + V**3) * Z + 18 * U * V * Z**2 - 27 * Z**4
U_st = t * (2 * s**3 + t**3)
V_st = s * (s**3 + 2 * t**3)
Z_st = s**2 * t**2

zero(delta.subs({U: U_st, V: V_st, Z: Z_st}), "quartic parametrization")
gate(sp.Poly(delta, U, V, Z).total_degree() == 4, "quartic degree")
factor_unit, factor_packet = sp.factor_list(delta)
gate(
    factor_unit in (-1, 1)
    and len(factor_packet) == 1
    and factor_packet[0][1] == 1
    and sp.Poly(factor_packet[0][0], U, V, Z).total_degree() == 4,
    "quartic irreducibility control",
)

# A dense rational inverse.  On Z=1 and q=s/t it is
# q=(2u^2-6v)/(u^2 v+3u-4v^2).
q = sp.symbols("q")
u_q = 2 * q + q**-2
v_q = q**2 + 2 * q**-1
q_inverse = sp.cancel(
    (2 * u_q**2 - 6 * v_q) / (u_q**2 * v_q + 3 * u_q - 4 * v_q**2)
)
gate(q_inverse == q, "birational inverse on dense normalization chart")

# Exact Fourier equivalence with the independently promoted THM-3851 model
# Delta_3=(AB+AC+BC)^2-4ABC(A+B+C).  For omega^2+omega+1=0,
# U=A+omega B+omega^2 C, V=A+omega^2 B+omega C, Z=(A+B+C)/3
# sends D to 9 Delta_3 and the canonical bitangent A+B+C=0 to Z=0.
A0, B0, C0, omega = sp.symbols("A0 B0 C0 omega")
delta_three = (
    (A0 * B0 + A0 * C0 + B0 * C0) ** 2
    - 4 * A0 * B0 * C0 * (A0 + B0 + C0)
)
U_fourier = A0 + omega * B0 + omega**2 * C0
V_fourier = A0 + omega**2 * B0 + omega * C0
Z_fourier = (A0 + B0 + C0) / 3
fourier_difference = sp.rem(
    sp.Poly(
        sp.expand(
            delta.subs({U: U_fourier, V: V_fourier, Z: Z_fourier})
            - 9 * delta_three
        ),
        omega,
    ),
    sp.Poly(omega**2 + omega + 1, omega),
).as_expr()
zero(fourier_difference, "Fourier equivalence with canonical tricuspidal model")
fourier_matrix = sp.Matrix(
    [
        [1, omega, omega**2],
        [1, omega**2, omega],
        [sp.Rational(1, 3), sp.Rational(1, 3), sp.Rational(1, 3)],
    ]
)
fourier_determinant = sp.rem(
    sp.Poly(sp.expand(fourier_matrix.det()), omega),
    sp.Poly(omega**2 + omega + 1, omega),
).as_expr()
gate(fourier_determinant != 0, "Fourier projective change is invertible")

# Exact singular support in the affine chart Z=1.
delta_aff = sp.expand(delta.subs(Z, 1))
singular_gb = sp.groebner(
    [delta_aff, sp.diff(delta_aff, U), sp.diff(delta_aff, V)],
    V,
    U,
    domain=sp.QQ,
)
gate(
    sp.expand(singular_gb.polys[1].as_expr()) == U**6 - 54 * U**3 + 729,
    "singular eliminant is (U^3-27)^2",
)
gate(
    sp.expand(U**6 - 54 * U**3 + 729) == sp.expand((U**3 - 27) ** 2),
    "three reduced cusp addresses",
)
cusp_relation_numerator = sp.cancel(
    q * (3 * V - U**2).subs({U: 3 * q, V: 3 / q})
)
zero(cusp_relation_numerator.subs(q**3, 1), "cusp image relation")

# At q=1, the first derivative vanishes, the second derivative is nonzero,
# and the second/third derivative vectors are independent: an A2 cusp.
u2 = sp.diff(u_q, q, 2).subs(q, 1)
v2 = sp.diff(v_q, q, 2).subs(q, 1)
u3 = sp.diff(u_q, q, 3).subs(q, 1)
v3 = sp.diff(v_q, q, 3).subs(q, 1)
gate(sp.diff(u_q, q).subs(q, 1) == 0, "cusp first derivative U")
gate(sp.diff(v_q, q).subs(q, 1) == 0, "cusp first derivative V")
gate((u2, v2) == (6, 6), "cusp quadratic tangent vector")
gate(sp.det(sp.Matrix([[u2, v2], [u3, v3]])) != 0, "cusp A2 cubic escape")

# The line Z=0 is a split bitangent: two smooth normalization addresses,
# each with contact two.  The target points are smooth because Delta_Z != 0.
gate(Z_st == (s * t) ** 2, "split bitangent pullback")
gate(sp.diff(delta, Z).subs({U: 1, V: 0, Z: 0}) == -4, "first infinity smooth")
gate(sp.diff(delta, Z).subs({U: 0, V: 1, Z: 0}) == -4, "second infinity smooth")

# A cusp tangent is the sharp hostile: it has contact partition 3+1, not 4.
gate(sp.factor(U_st - V_st) == -(s - t) ** 3 * (s + t), "cusp tangent 3+1")

# Uniform no-fourfold-line gate.  A line aU+bV+cZ pulls back with coefficient
# vector (b,2a,c,2b,a).  Saturate the equality with lambda(alpha s+beta t)^4
# on the three projective root charts.  Every saturated ideal is the unit
# ideal, so no line section is supported at one normalization address.
a, b, c, lam, alpha, beta, inv = sp.symbols(
    "a b c lam alpha beta inv"
)
pull_coefficients = (b, 2 * a, c, 2 * b, a)
fourth_coefficients = (
    lam * alpha**4,
    4 * lam * alpha**3 * beta,
    6 * lam * alpha**2 * beta**2,
    4 * lam * alpha * beta**3,
    lam * beta**4,
)
coefficient_equations = [
    left - right for left, right in zip(pull_coefficients, fourth_coefficients)
]

interior_gb = sp.groebner(
    coefficient_equations + [inv * lam * alpha * beta - 1],
    a,
    b,
    c,
    lam,
    alpha,
    beta,
    inv,
    domain=sp.QQ,
)
alpha_zero_gb = sp.groebner(
    coefficient_equations + [alpha, inv * lam * beta - 1],
    a,
    b,
    c,
    lam,
    alpha,
    beta,
    inv,
    domain=sp.QQ,
)
beta_zero_gb = sp.groebner(
    coefficient_equations + [beta, inv * lam * alpha - 1],
    a,
    b,
    c,
    lam,
    alpha,
    beta,
    inv,
    domain=sp.QQ,
)
gate(
    len(interior_gb.polys) == 1 and interior_gb.polys[0].as_expr() == 1,
    "no interior fourth power",
)
gate(
    len(alpha_zero_gb.polys) == 1 and alpha_zero_gb.polys[0].as_expr() == 1,
    "no alpha-zero fourth power",
)
gate(
    len(beta_zero_gb.polys) == 1 and beta_zero_gb.polys[0].as_expr() == 1,
    "no beta-zero fourth power",
)


# ---------------------------------------------------------------------------
# II. Reciprocal-cubic and A2-torus explanation.
# ---------------------------------------------------------------------------
T = sp.symbols("T")
reciprocal_cubic = T**3 - U * T**2 + V * T - 1
zero(
    sp.discriminant(reciprocal_cubic, T) - delta_aff,
    "tricuspidal quartic is the reciprocal-cubic discriminant",
)
zero(
    reciprocal_cubic - (T * (T**2 - U * T + V) - 1),
    "the cubic root is a global unit",
)
V_source = T**-1 - T**2 + U * T
zero(
    sp.cancel(reciprocal_cubic.subs(V, V_source)),
    "cubic source is A1_U times G_m_T",
)

x, y = sp.symbols("x y", nonzero=True)
z = 1 / (x * y)
U_xy = x + y + z
V_xy = x * y + 1 / x + 1 / y
W_xy = (x - y) * (y - z) * (z - x)
zero(
    sp.cancel(delta_aff.subs({U: U_xy, V: V_xy}) - W_xy**2),
    "ordered-root Vandermonde squares to the discriminant",
)

# The logarithmic Jacobian is the alternating Vandermonde.  This identifies
# the same divisor as the torus quotient branch, not merely as a polynomial
# discriminant with matching support.
log_jacobian = sp.det(
    sp.Matrix(
        [
            [x * sp.diff(U_xy, x), y * sp.diff(U_xy, y)],
            [x * sp.diff(V_xy, x), y * sp.diff(V_xy, y)],
        ]
    )
)
zero(sp.cancel(log_jacobian + W_xy), "torus quotient logarithmic Jacobian")

# A transposition negates W; the 3-cycle (x,y,z)->(y,z,x) fixes it.  These
# exact substitutions are the A3/S3 quadratic-resolvent dictionary.
zero(
    sp.cancel(U_xy.subs({x: y, y: x}, simultaneous=True) - U_xy),
    "transposition fixes U",
)
zero(
    sp.cancel(V_xy.subs({x: y, y: x}, simultaneous=True) - V_xy),
    "transposition fixes V",
)
zero(
    sp.cancel(W_xy.subs({x: y, y: x}, simultaneous=True) + W_xy),
    "transposition negates W",
)
cycle_substitution = {x: y, y: z}
zero(
    sp.cancel(U_xy.subs(cycle_substitution, simultaneous=True) - U_xy),
    "three-cycle fixes U",
)
zero(
    sp.cancel(V_xy.subs(cycle_substitution, simultaneous=True) - V_xy),
    "three-cycle fixes V",
)
zero(
    sp.cancel(W_xy.subs(cycle_substitution, simultaneous=True) - W_xy),
    "three-cycle fixes W",
)


# ---------------------------------------------------------------------------
# III. The 3A2 degree-two del Pezzo root/boundary quotient.
# ---------------------------------------------------------------------------
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
gate(all(lattice_dot(root, root) == -2 for root in roots), "E7 root squares")

boundary_plus = (0, 0, 0, 0, 0, 0, 0, 1)
anticanonical = (3, -1, -1, -1, -1, -1, -1, -1)
canonical = lattice_neg(anticanonical)
boundary_minus = tuple(
    anticanonical[index] - boundary_plus[index] for index in range(8)
)
gate(lattice_dot(boundary_plus, boundary_plus) == -1, "first boundary square")
gate(lattice_dot(boundary_minus, boundary_minus) == -1, "second boundary square")
gate(lattice_dot(boundary_plus, boundary_minus) == 2, "boundary total contact")
gate(sp.Matrix([boundary_plus, boundary_minus]).rank() == 2, "boundary classes independent")

orthogonal_roots = sorted(
    root for root in roots if lattice_dot(root, boundary_plus) == 0
)
root_lookup = set(orthogonal_roots)
gate(len(orthogonal_roots) == 72, "boundary-orthogonal E6 roots")

a2_systems: dict[
    tuple[tuple[int, ...], ...], tuple[tuple[int, ...], tuple[int, ...]]
] = {}
for left, right in combinations(orthogonal_roots, 2):
    if lattice_dot(left, right) != 1:
        continue
    root_sum = tuple(left[index] + right[index] for index in range(8))
    subsystem = tuple(
        sorted(
            {
                left,
                right,
                lattice_neg(left),
                lattice_neg(right),
                root_sum,
                lattice_neg(root_sum),
            }
        )
    )
    gate(all(root in root_lookup for root in subsystem), "A2 root closure")
    a2_systems[subsystem] = (left, right)
gate(len(a2_systems) == 120, "complete boundary-orthogonal A2 systems")

a2_items = list(a2_systems.items())
compatible: dict[int, set[int]] = {index: set() for index in range(len(a2_items))}
for first in range(len(a2_items)):
    first_roots = a2_items[first][0]
    for second in range(first + 1, len(a2_items)):
        second_roots = a2_items[second][0]
        if all(
            lattice_dot(left, right) == 0
            for left in first_roots
            for right in second_roots
        ):
            compatible[first].add(second)
            compatible[second].add(first)

configurations: list[
    tuple[
        tuple[tuple[int, ...], tuple[int, ...]],
        tuple[tuple[int, ...], tuple[int, ...]],
        tuple[tuple[int, ...], tuple[int, ...]],
    ]
] = []
for first in range(len(a2_items)):
    for second in sorted(index for index in compatible[first] if index > first):
        for third in sorted(
            index
            for index in compatible[first].intersection(compatible[second])
            if index > second
        ):
            configurations.append(
                (
                    a2_items[first][1],
                    a2_items[second][1],
                    a2_items[third][1],
                )
            )
gate(len(configurations) == 40, "complete 3A2 configurations in E6")

smith_histogram: Counter[tuple[int, ...]] = Counter()
normal_histogram: Counter[tuple[int, int, int]] = Counter()
local_images: list[set[tuple[int, int, int]]] = []

for pairs in configurations:
    relation_rows = [
        pairs[0][0],
        pairs[0][1],
        pairs[1][0],
        pairs[1][1],
        pairs[2][0],
        pairs[2][1],
        boundary_plus,
        boundary_minus,
    ]
    relation_matrix = sp.Matrix(relation_rows)
    gate(relation_matrix.rank() == 8, "3A2 boundary quotient is finite")
    smith = smith_normal_form(relation_matrix, domain=ZZ)
    diagonal = tuple(abs(int(smith[index, index])) for index in range(8))
    smith_histogram[diagonal] += 1

    kernel_coefficients = [
        coefficients
        for coefficients in product(range(3), repeat=8)
        if all(
            sum(
                coefficients[row] * relation_rows[row][column]
                for row in range(8)
            )
            % 3
            == 0
            for column in range(8)
        )
    ]
    gate(len(kernel_coefficients) == 9, "rank-two mod-three saturation kernel")

    image: set[tuple[int, int, int]] = set()
    for coefficients in kernel_coefficients:
        numerator = tuple(
            sum(
                coefficients[row] * relation_rows[row][column]
                for row in range(8)
            )
            for column in range(8)
        )
        gate(all(entry % 3 == 0 for entry in numerator), "integral 3-torsion lift")
        lift = tuple(entry // 3 for entry in numerator)
        image.add(
            tuple(
                (
                    lattice_dot(lift, pair[0])
                    + 2 * lattice_dot(lift, pair[1])
                )
                % 3
                for pair in pairs
            )
        )
    gate(len(image) == 9, "local restriction injects the rank-two torsion")
    gate(
        all({entry[coordinate] for entry in image} == {0, 1, 2} for coordinate in range(3)),
        "every cusp projection is surjective",
    )

    annihilators = [
        normal
        for normal in product(range(3), repeat=3)
        if normal != (0, 0, 0)
        and all(sum(a * b for a, b in zip(normal, point)) % 3 == 0 for point in image)
    ]
    gate(len(annihilators) == 2, "one projective local relation")
    gate(all(all(coordinate != 0 for coordinate in normal) for normal in annihilators), "local relation has full cusp support")
    canonical_normal = min(annihilators)
    normal_histogram[canonical_normal] += 1
    local_images.append(image)

expected_smith = (1, 1, 1, 1, 1, 1, 3, 3)
gate(smith_histogram == Counter({expected_smith: 40}), "uniform affine class group (Z/3)^2")

# An explicit compatible marking, useful for a human proof candidate.
r1 = (0, 1, -1, 0, 0, 0, 0, 0)
r2 = (0, 0, 1, -1, 0, 0, 0, 0)
r3 = (0, 0, 0, 0, 1, -1, 0, 0)
r4 = (0, 0, 0, 0, 0, 1, -1, 0)
r5 = (1, -1, -1, -1, 0, 0, 0, 0)
r6 = (1, 0, 0, 0, -1, -1, -1, 0)
explicit_pairs = ((r1, r2), (r3, r4), (r5, r6))
gate(
    all(lattice_dot(pair[0], pair[1]) == 1 for pair in explicit_pairs),
    "explicit A2 simple pairs",
)
gate(
    all(
        lattice_dot(left, right) == 0
        for first, second in combinations(explicit_pairs, 2)
        for left in first
        for right in second
    ),
    "explicit A2 systems orthogonal",
)
explicit_matrix = sp.Matrix(
    [r1, r2, r3, r4, r5, r6, boundary_plus, boundary_minus]
)
explicit_smith = smith_normal_form(explicit_matrix, domain=ZZ)
gate(
    tuple(abs(int(explicit_smith[index, index])) for index in range(8))
    == expected_smith,
    "explicit quotient Smith form",
)

# Two transparent order-three lifts.  Reverse the orientation of the second
# A2 pair; their local character rows then span x+y+z=0.
torsion_one = (1, 0, 0, -1, 0, 0, -1, 0)
torsion_two = (0, 0, 0, -1, 0, 0, 0, 1)
relation_one_coefficients = (1, 2, 1, 2, 0, 0, 2, 1)
relation_two_coefficients = (1, 2, 0, 0, 2, 1, 1, -1)
gate(
    tuple(
        sum(
            relation_one_coefficients[row] * explicit_matrix[row, column]
            for row in range(8)
        )
        for column in range(8)
    )
    == tuple(3 * entry for entry in torsion_one),
    "first explicit order-three relation",
)
gate(
    tuple(
        sum(
            relation_two_coefficients[row] * explicit_matrix[row, column]
            for row in range(8)
        )
        for column in range(8)
    )
    == tuple(3 * entry for entry in torsion_two),
    "second explicit order-three relation",
)
oriented_pairs = ((r1, r2), (r4, r3), (r5, r6))


def local_character(divisor: tuple[int, ...]) -> tuple[int, int, int]:
    return tuple(
        (
            lattice_dot(divisor, pair[0])
            + 2 * lattice_dot(divisor, pair[1])
        )
        % 3
        for pair in oriented_pairs
    )


explicit_local_generators = (
    local_character(torsion_one),
    local_character(torsion_two),
)
gate(
    explicit_local_generators == ((1, 2, 0), (1, 0, 2)),
    "explicit local hyperplane generators",
)
gate(
    all(sum(character) % 3 == 0 for character in explicit_local_generators),
    "explicit local relation x+y+z=0",
)

# Geiser acts as -1 after both boundary classes are killed: sigma(D)+D is
# (D.K)K, and K=-(B_++B_-) vanishes in the affine quotient.
for basis_vector in [tuple(1 if index == coordinate else 0 for index in range(8)) for coordinate in range(8)]:
    geiser = tuple(
        lattice_dot(basis_vector, canonical) * canonical[index]
        - basis_vector[index]
        for index in range(8)
    )
    geiser_plus_identity = tuple(
        geiser[index] + basis_vector[index] for index in range(8)
    )
    multiple = lattice_dot(basis_vector, canonical)
    gate(
        geiser_plus_identity
        == tuple(-multiple * (boundary_plus[index] + boundary_minus[index]) for index in range(8)),
        "Geiser inversion modulo boundary",
    )


# ---------------------------------------------------------------------------
# IV. Deterministic transcript.
# ---------------------------------------------------------------------------
normal_histogram_text = ",".join(
    f"{normal}:{count}" for normal, count in sorted(normal_histogram.items())
)
semantic_lines = [
    "THM3851_TRICUSPIDAL_TRADEOFF_EXACT PASS",
    "CURVE U^2V^2-4(U^3+V^3)Z+18UVZ^2-27Z^4=0",
    "NORMALIZATION [s:t]->[t(2s^3+t^3):s(s^3+2t^3):s^2t^2]",
    "FOURIER_EQUIVALENCE U=A+wB+w^2C,V=A+w^2B+wC,Z=(A+B+C)/3;D=9Delta_3",
    "SINGULARITIES 3A2 at q^3=1",
    "SPLIT_BITANGENT Z=0 pulls back to s^2t^2 (two smooth places)",
    "NO_FOURFOLD_LINE every line pullback has support size at least two",
    "SHARP_HOSTILE cusp tangent U-V=0 pulls back to -(s-t)^3(s+t)",
    "RECIPROCAL_CUBIC disc(T^3-UT^2+VT-1)=delta and T is a global unit",
    "CUBIC_SOURCE Spec=k[U,T,T^-1]=A1_x_Gm (finite flat rank three over A2)",
    "ORDERED_ROOT_COVER (x,y,1/(xy)) in Gm^2;Vandermonde^2=delta",
    "RESOLVENT_DICTIONARY Gm^2/A3=Q and Gm^2/S3=A2 at function-field level",
    "E7_ROOTS 126;BOUNDARY_E6_ROOTS 72;A2_SYSTEMS 120;3A2_CONFIGS 40",
    "AFFINE_CLASS_GROUP every 3A2+split-boundary quotient has Smith 1,1,1,1,1,1,3,3",
    "AFFINE_UNITS k* because the two boundary classes are independent in Pic(S)",
    "LOCAL_RESTRICTION (Z/3)^2 injects into (Z/3)^3 as a full-support hyperplane",
    "EXPLICIT_LOCAL_GENERATORS (1,2,0),(1,0,2), hence x+y+z=0 after orientation",
    f"LOCAL_NORMAL_HISTOGRAM {normal_histogram_text}",
    "DECK_ACTION Geiser acts by inversion modulo the two boundary classes",
]
semantic_sha = hashlib.sha256("\n".join(semantic_lines).encode("utf-8")).hexdigest()
transcript = "\n".join(semantic_lines + [f"GATES {GATES}", f"SEMANTIC_SHA256 {semantic_sha}"]) + "\n"

if len(sys.argv) == 3 and sys.argv[1] == "--verify-frozen":
    frozen = Path(sys.argv[2]).read_bytes()
    gate(frozen == transcript.encode("utf-8"), "frozen transcript byte match")
    print("FROZEN_REPLAY PASS")
    print(f"FROZEN_SHA256 {hashlib.sha256(frozen).hexdigest()}")
elif len(sys.argv) == 1:
    sys.stdout.write(transcript)
else:
    raise SystemExit("usage: tricuspidal_tradeoff_exact.py [--verify-frozen PATH]")
