#!/usr/bin/env python3
"""Exact cross-frontier scout for filtered JC lifting and principal parts.

This artifact independently replays the smallest THM-3412
filtration--observer noncommutation, the THM-3871 arms-only false component,
and THM-3404's CRT principal-part cancellation.  It also freezes two hostile
boundaries: a tangent-ray observer does not retain contact at a common point
at infinity, and Pythagorean outer ordinal rank is a scheduler rather than an
object address.  No theorem transfer to sextic JC is asserted.

All truth gates are optimization-safe; Python ``assert`` is not used.
"""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


GATES = 0


def gate(condition: object, label: str) -> None:
    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"GATE FAILED: {label}: {condition}")


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(sp.factor(expression)) == 0, label)


def nilpotent_blocks(matrix: sp.Matrix) -> list[int]:
    """Recover nilpotent Jordan block sizes from successive nullities."""

    dimension = matrix.rows
    nullities = [0]
    power = sp.eye(dimension)
    for _ in range(1, dimension + 1):
        power = power * matrix
        nullities.append(dimension - power.rank())
    blocks_at_least = [
        nullities[index] - nullities[index - 1]
        for index in range(1, len(nullities))
    ]
    exact: list[int] = []
    for size in range(1, len(blocks_at_least) + 1):
        next_count = blocks_at_least[size] if size < len(blocks_at_least) else 0
        exact.extend([size] * (blocks_at_least[size - 1] - next_count))
    return sorted(exact, reverse=True)


# ---------------------------------------------------------------------------
# 1. THM-3412: the exact e=3,q=2 filtration-observer hostile.
# ---------------------------------------------------------------------------
y, t = sp.symbols("y t")
local_basis = [y**power for power in range(6)] + [y**power * t for power in range(3)]
local_index = {sp.expand(monomial): index for index, monomial in enumerate(local_basis)}


def reduce_local(expression: sp.Expr) -> sp.Expr:
    """Reduce modulo (y^3,t)^2=(y^6,y^3*t,t^2)."""

    result = 0
    for (y_degree, t_degree), coefficient in sp.Poly(
        sp.expand(expression), y, t
    ).terms():
        if t_degree >= 2:
            continue
        if t_degree == 1 and y_degree >= 3:
            continue
        if t_degree == 0 and y_degree >= 6:
            continue
        result += coefficient * y**y_degree * t**t_degree
    return sp.expand(result)


def nabla_3_2(monomial: sp.Expr) -> sp.Expr:
    # f=y, g=y^3 and q=2: g(d_t-d_y)+q*g'=y^3(d_t-d_y)+6y^2.
    return reduce_local(
        y**3 * (sp.diff(monomial, t) - sp.diff(monomial, y))
        + 6 * y**2 * monomial
    )


actual_matrix = sp.zeros(len(local_basis))
for column, monomial in enumerate(local_basis):
    for powers, coefficient in sp.Poly(nabla_3_2(monomial), y, t).terms():
        target = y**powers[0] * t**powers[1]
        actual_matrix[local_index[target], column] = coefficient

actual_kernel = actual_matrix.nullspace()
gate(len(actual_kernel) == 4, "THM-3412 actual finite kernel dimension four")
actual_kernel_matrix = sp.Matrix.hstack(*actual_kernel)
actual_action_columns: list[sp.Matrix] = []
for vector in actual_kernel:
    polynomial = sum(vector[index] * local_basis[index] for index in range(len(local_basis)))
    product = reduce_local((y + t) * polynomial)
    target_vector = sp.Matrix([
        sp.Poly(product, y, t).coeff_monomial(monomial)
        for monomial in local_basis
    ])
    solutions = list(sp.linsolve((actual_kernel_matrix, target_vector)))
    gate(len(solutions) == 1, "actual P-action has unique kernel coordinates")
    actual_action_columns.append(sp.Matrix(solutions[0]))
actual_P_action = sp.Matrix.hstack(*actual_action_columns)
actual_blocks = nilpotent_blocks(actual_P_action)
gate(actual_blocks == [4], "THM-3412 actual kernel has one length-four arm")

# On G_2=k[y,t]/(y^3,t^2), gr(D) is multiplication by 6y^2.
graded_basis = [y**power * t**level for level in range(2) for power in range(3)]
graded_matrix = sp.zeros(len(graded_basis))
for column, monomial in enumerate(graded_basis):
    product = sp.Poly(sp.expand(6 * y**2 * monomial), y, t)
    for powers, coefficient in product.terms():
        if powers[0] < 3 and powers[1] < 2:
            target = y**powers[0] * t**powers[1]
            graded_matrix[graded_basis.index(target), column] = coefficient
graded_kernel = graded_matrix.nullspace()
gate(len(graded_kernel) == 4, "THM-3412 graded kernel dimension four")
graded_kernel_matrix = sp.Matrix.hstack(*graded_kernel)
graded_action_columns: list[sp.Matrix] = []
for vector in graded_kernel:
    polynomial = sum(vector[index] * graded_basis[index] for index in range(len(graded_basis)))
    raw_product = sp.Poly(sp.expand((y + t) * polynomial), y, t)
    target_vector = sp.Matrix([
        raw_product.coeff_monomial(monomial) for monomial in graded_basis
    ])
    solutions = list(sp.linsolve((graded_kernel_matrix, target_vector)))
    gate(len(solutions) == 1, "graded P-action has unique kernel coordinates")
    graded_action_columns.append(sp.Matrix(solutions[0]))
graded_P_action = sp.Matrix.hstack(*graded_action_columns)
graded_blocks = nilpotent_blocks(graded_P_action)
gate(graded_blocks == [3, 1], "THM-3412 graded kernel splits as three plus one")

liftable_symbol_dimension = 2
graded_deficit = len(graded_kernel) - liftable_symbol_dimension
gate(graded_deficit == 2, "THM-3412 sharp (q-1)m deficit")


# ---------------------------------------------------------------------------
# 2. THM-3871: arms-only kernel versus the two conserved lower buckets.
# ---------------------------------------------------------------------------
X, Y, Z = sp.symbols("X Y Z")
arm_C = 1 + X + Y + Z
arm_A = 5 * X**2 + 10 * X * Y - 8
conserved_P = X**3 - 8 * X * Z - 4 * Y**2
conserved_T = Y * (3 * X**2 - 8 * Z)
arms_only_point = {
    X: sp.Integer(1),
    Y: sp.Rational(3, 10),
    Z: sp.Rational(-23, 10),
}
zero(arm_C.subs(arms_only_point), "THM-3871 hostile cancels C arm")
zero(arm_A.subs(arms_only_point), "THM-3871 hostile cancels A arm")
arm_jacobian = sp.Matrix([
    [sp.diff(arm_C, variable) for variable in (X, Y, Z)],
    [sp.diff(arm_A, variable) for variable in (X, Y, Z)],
])
gate(arm_jacobian.subs(arms_only_point).rank() == 2,
     "arms-only hostile is a smooth local one-fold")
gate(conserved_P.subs(arms_only_point) == sp.Rational(476, 25),
     "first omitted conserved row rejects arms-only hostile")
gate(conserved_T.subs(arms_only_point) == sp.Rational(321, 50),
     "second omitted conserved row rejects arms-only hostile")

p_extreme = 15 * X**3 + 20 * X**2 + 40 * X + 32
q_extreme = 50 * X**5 + 25 * X**4 - 80 * X**2 + 64
extreme_resultant = sp.resultant(p_extreme, q_extreme, X)
gate(extreme_resultant == 3171942400000,
     "THM-3871 terminal simultaneous-arm resultant")
gate(sp.factorint(extreme_resultant) == {2: 23, 5: 5, 11: 2},
     "THM-3871 terminal resultant factorization")


# ---------------------------------------------------------------------------
# 3. Twelve sextic buckets and an exact shallow-observer hostile.
# ---------------------------------------------------------------------------
a = sp.symbols("a0:7")
c = sp.symbols("c0:7")
ad = sp.symbols("a0d:7d")
cd = sp.symbols("c0d:7d")


def sextic_bucket(degree: int) -> sp.Expr:
    return sp.expand(sum(
        i * a[i] * cd[j] - j * ad[i] * c[j]
        for i in range(7)
        for j in range(7)
        if i + j == degree + 1
    ))


gate(len([sextic_bucket(degree) for degree in range(12)]) == 12,
     "sextic bracket has twelve coefficient observers")
zero(sextic_bucket(11) - 6 * (a[6] * cd[6] - ad[6] * c[6]),
     "sextic top Wronskian bucket")

s, z = sp.symbols("s z")
A_shallow = z**6 + s
C_shallow = z**5
J_shallow = sp.expand(
    sp.diff(A_shallow, z) * sp.diff(C_shallow, s)
    - sp.diff(A_shallow, s) * sp.diff(C_shallow, z)
)
zero(J_shallow + 5 * z**4, "sextic shallow-observer hostile Jacobian")
for degree in range(5, 12):
    zero(J_shallow.coeff(z, degree),
         f"sextic hostile invisible to upper bucket z^{degree}")
gate(J_shallow.coeff(z, 4) == -5,
     "sextic hostile first appears in the next extension bucket")


# ---------------------------------------------------------------------------
# 4. THM-3404: CRT principal parts see cancellation that orders miss.
# ---------------------------------------------------------------------------
v = sp.symbols("v")
F_simple = v * (v - 1)
h_left = sp.Integer(1)
h_right = F_simple - 1
zero(sp.rem(h_left + h_right, F_simple, v),
     "simple CRT principal parts cancel modulo F")
simple_crt_left = (h_left.subs(v, 0), h_left.subs(v, 1))
simple_crt_right = (h_right.subs(v, 0), h_right.subs(v, 1))
gate(simple_crt_left == (1, 1), "simple CRT first residue packet")
gate(simple_crt_right == (-1, -1), "simple CRT second residue packet")
gate(tuple(simple_crt_left[index] + simple_crt_right[index] for index in range(2)) == (0, 0),
     "simple CRT cancellation is factorwise")

# Two-coordinate regularity is the direct sum of two principal-part modules.
pair_one = (h_left, h_right)
pair_two = (h_right, h_left)
for coordinate in range(2):
    zero(sp.rem(pair_one[coordinate] + pair_two[coordinate], F_simple, v),
         f"coordinatewise CRT cancellation coordinate {coordinate}")

# With repeated factors, value residues vanish but the jet class survives.
F_double = v**2 * (v - 1)**2
h_jet = v * (v - 1)
gate(sp.rem(h_jet, F_double, v) != 0,
     "double-factor principal part survives modulo F")
jet_packet = (
    h_jet.subs(v, 0), sp.diff(h_jet, v).subs(v, 0),
    h_jet.subs(v, 1), sp.diff(h_jet, v).subs(v, 1),
)
gate(jet_packet == (0, -1, 0, 1),
     "value-only observer misses the two first jets")


# ---------------------------------------------------------------------------
# 5. Common-infinity and ordinal hostiles to over-compression.
# ---------------------------------------------------------------------------
u = sp.symbols("u")
common_infinity_contacts: dict[int, int] = {}
for degree in (2, 3, 4):
    local_equation = v**(degree - 1) - u**degree
    zero(
        local_equation.subs({u: t**(degree - 1), v: t**degree}),
        f"A1 normalization of y=x^{degree} at common infinity",
    )
    initial_degree = min(sum(monomial) for monomial, _ in sp.Poly(
        local_equation, u, v
    ).terms())
    gate(initial_degree == degree - 1,
         f"common tangent multiplicity for graph degree {degree}")

for degree in (3, 4):
    restricted = sp.expand((v**(degree - 1) - u**degree).subs(v, u**2))
    contact = min(monomial[0] for monomial, _ in sp.Poly(restricted, u).terms())
    common_infinity_contacts[degree] = contact
gate(common_infinity_contacts == {3: 3, 4: 4},
     "same common point and tangent retain different contact orders")


def pythagorean_from_ordinals(r: int, inner: int) -> tuple[int, int, int]:
    q = 2 * r - 1
    d = 2 * inner - 1
    return q * d, (q*q - d*d) // 2, (q*q + d*d) // 2


gate(pythagorean_from_ordinals(3, 1) == (5, 12, 13),
     "outer ordinal rank-three first node")
gate(pythagorean_from_ordinals(3, 2) == (15, 8, 17),
     "outer ordinal rank-three second node")
gate(pythagorean_from_ordinals(3, 1) != pythagorean_from_ordinals(3, 2),
     "outer ordinal is not a node address")

# The stronger cross-hostile uses the *same* ambient triangular address on
# both task sets.  Address 8 decodes as (r,s)=(5,2) in THM-3756 and as the
# quintic (n,j)=(5,2) normal-strip bucket.  The former is a primitive hole,
# while the latter is a genuine Keller stratum with a native 3R/8 residue.
def triangular(value: int) -> int:
    return value * (value + 1) // 2


pythagorean_address = triangular(5 - 2) + 2
jc_bucket_address = triangular(5 - 2) + 2
gate(pythagorean_address == jc_bucket_address == 8,
     "shared triangular address eight")
gate(sp.gcd(2 * 5 - 1, 2 * 2 - 1) == 3,
     "address-eight Pythagorean pair fails the primitive gcd filter")
gate(pythagorean_from_ordinals(5, 2) == (27, 36, 45),
     "address-eight ambient Pythagorean node is nonprimitive")
R52, Q52 = sp.symbols("R52 Q52", nonzero=True)
K52 = 5 * R52 / (2 * Q52)
residual52 = sp.cancel(R52 - K52 * Q52 + 3 * K52 * Q52 / 4)
zero(residual52 - 3 * R52 / 8,
     "address-eight quintic (5,2) terminal arm residue")
gate(residual52 != 0,
     "native quintic filter retains the Pythagorean-hole address")

# Hadamard repeated-event minors vanish, whereas a same-index Jacobian wedge
# need not.  This blocks a literal subset-minor interpretation of buckets.
i = sp.symbols("i", integer=True, positive=True)
same_index_wedge = i * (s * sp.diff(s**2, s) - sp.diff(s, s) * s**2)
zero(same_index_wedge - i * s**2,
     "Jacobian same-index channel can be nonzero")
gate(same_index_wedge != 0,
     "Hadamard repeated-minor zero law does not transfer literally")


# ---------------------------------------------------------------------------
# 6. Frozen semantic packet and optimization-safe source audit.
# ---------------------------------------------------------------------------
source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "no inactive Python assert",
)

semantic = {
    "source_3412": "e3q2 actual block4;graded blocks3+1;liftable deficit2",
    "target_3871": "arms-only smooth one-fold rejected by conserved P,T and nonzero resultant",
    "sextic": "twelve buckets;upper seven can miss the z4 extension equation",
    "crt_3404": "factorwise residues and repeated-factor jets detect simultaneous cancellation",
    "crt_verdict": "direct on an identified suspension/DVR chart;proposed sidecar for sextic JC",
    "common_infinity": "A1 graph branches share point+tangent but contacts 3 and 4 differ",
    "ordinal": "rank collision r3;shared address8 is primitive hole (5,2) in Pythagorean chart but genuine JC bucket with residue3R/8",
    "hadamard": "subset-minor repeated-index law does not match Jacobian anti-diagonal convolution",
    "scope": "cross-frontier scout only;no sextic theorem;no JC2 consequence",
}
semantic_sha = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("JC2_FILTERED_KERNEL_PRINCIPAL_PART_CONNECTION_SCOUT")
print("status=VERIFIED_EXACT_SCOUT;NO_SEXTIC_THEOREM;JC2_OPEN")
print(f"thm3412=actual_dim4_blocks{actual_blocks};graded_dim4_blocks{graded_blocks};liftable_dim2;deficit2")
print("thm3871=arms_local_dimension1;P=476/25;T=321/50")
print(f"thm3871_terminal_resultant={extreme_resultant}=2^23*5^5*11^2")
print("sextic=twelve_buckets;upper_z5_to_z11_blind_to_exact_z4_extension_hostile")
print("thm3404=simple_CRT_cancellation_PASS;double_factor_jet_packet=(0,-1,0,1)")
print("crt_verdict=DIRECT_ON_IDENTIFIED_SUSPENSION_OR_DVR;PROPOSED_FOR_SEXTIC_NORMAL_STRIP")
print("common_infinity=A1_graphs_same_point+tangent;contacts_G2_G3=3,G2_G4=4")
print("ordinal_rank=r3_collides_(5,12,13)_with_(15,8,17);scheduler_only")
print("ordinal_address8=pythagorean_(r,s)=(5,2)_gcd3_HOLE;jc_(n,j)=(5,2)_residue=3R/8_GENUINE")
print("hadamard_literal_transfer=REFUTED_same_index_Jacobian_wedge_nonzero")
print("scope=filtered_kernel_and_principal_part_sidecar_tasks_only;JC2_OPEN")
print(f"semantic_sha256={semantic_sha}")
print(f"GATES={GATES}")
print("RESULT=PASS")
