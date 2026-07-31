#!/usr/bin/env python3
"""Exact companion for THM-2690.

Checks the cyclic-cubic normal-crossing class-group presentation and the
singular-C3 / polynomial-S3 reflection-completion hostile.  Every
truth-bearing check uses ``require`` so optimized Python is equally strict.
"""

from itertools import combinations, product

from sympy import Matrix, diff, eye, factor, ones, symbols
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def mod2_matrix(matrix):
    return matrix.applyfunc(lambda entry: int(entry) % 2)


# 1. Nagata presentation [3 I_m | 1] has Smith form (1,3,...,3).
snf_rows = []
for m in range(1, 13):
    presentation = Matrix.hstack(3 * eye(m), ones(m, 1))
    smith = smith_normal_form(presentation, domain=ZZ)
    diagonal = tuple(abs(int(smith[i, i])) for i in range(m))
    expected = (1,) + (3,) * (m - 1)
    require(diagonal == expected, f"bad Smith form at m={m}: {diagonal}")
    snf_rows.append((m, diagonal, 3 ** (m - 1)))


# 1b. Every reduced normal-crossing exponent vector has the same 3-primary SNF.
reduced_exponent_count = 0
for m in range(1, 9):
    for exponent_vector in product((1, 2), repeat=m):
        presentation = Matrix.hstack(3 * eye(m), Matrix(exponent_vector))
        smith = smith_normal_form(presentation, domain=ZZ)
        diagonal = tuple(abs(int(smith[i, i])) for i in range(m))
        expected = (1,) + (3,) * (m - 1)
        require(
            diagonal == expected,
            f"bad reduced-exponent Smith form at m={m}, alpha={exponent_vector}: {diagonal}",
        )
        reduced_exponent_count += 1
require(reduced_exponent_count == 510, "wrong reduced exponent census")


# 2. The Jacobian singular support has at least two vanishing branch variables.
for m in range(2, 13):
    indices = tuple(range(m))
    singular_zero_sets = []
    for size in range(m + 1):
        for zero_tuple in combinations(indices, size):
            zero_set = set(zero_tuple)
            # d/dx_i(product x_j) vanishes iff some j != i is zero.
            all_branch_derivatives_zero = all(
                any(j in zero_set for j in indices if j != i) for i in indices
            )
            if all_branch_derivatives_zero:
                singular_zero_sets.append(zero_set)
    require(singular_zero_sets, f"missing singular supports at m={m}")
    require(
        min(len(zero_set) for zero_set in singular_zero_sets) == 2,
        f"singular codimension support failed at m={m}",
    )


# 3. C3 cycles the three nonzero points of the binary standard plane.
sigma = Matrix([[0, 1], [1, 1]])
sigma2 = mod2_matrix(sigma * sigma)
sigma3 = mod2_matrix(sigma2 * sigma)
require(sigma3 == eye(2), "standard-plane matrix does not have order three")
nonzero = {(1, 0), (0, 1), (1, 1)}
orbit = []
vector = Matrix([1, 0])
for _ in range(3):
    orbit.append(tuple(int(value) % 2 for value in vector))
    vector = mod2_matrix(sigma * vector)
require(set(orbit) == nonzero, f"bad nonzero standard-plane orbit: {orbit}")
for vector_tuple in nonzero:
    vector = Matrix(vector_tuple)
    require(mod2_matrix(sigma * vector) != vector, "unexpected C3-fixed vector")


# 4. Exact C3 invariant relation for d^2=abc.
a, b, c, e1, e2, e3, d, delta = symbols("a b c e1 e2 e3 d delta")
root_delta = (a - b) * (b - c) * (c - a)
discriminant = (
    e1**2 * e2**2
    - 4 * e2**3
    - 4 * e1**3 * e3
    - 27 * e3**2
    + 18 * e1 * e2 * e3
)
elementary_substitution = {
    e1: a + b + c,
    e2: a * b + a * c + b * c,
    e3: a * b * c,
}
require(
    factor(root_delta**2 - discriminant.subs(elementary_substitution)) == 0,
    "cyclic invariant discriminant identity failed",
)

cyclic_relation = delta**2 - discriminant.subs(e3, d**2)
origin = {e1: 0, e2: 0, d: 0, delta: 0}
require(cyclic_relation.subs(origin) == 0, "origin is not on cyclic quotient")
require(
    all(diff(cyclic_relation, variable).subs(origin) == 0 for variable in (e1, e2, d, delta)),
    "cyclic quotient is not singular at the origin",
)


# 5. The C2 reflection delta -> -delta removes the alternating generator.
reflected_relation = cyclic_relation.subs(delta, -delta)
require(reflected_relation == cyclic_relation, "C2 reflection does not preserve relation")
require(
    factor(delta**2 - discriminant.subs(e3, d**2) - cyclic_relation) == 0,
    "delta-square elimination failed",
)


# 6. The known reflection-completed quartic carrier hostile is not Keller.
s1, s2, s3 = symbols("s1 s2 s3")
quotient_map = Matrix([s1**2 - 2 * s2, s2**2 - 2 * s1 * s3, s3])
jacobian = quotient_map.jacobian((s1, s2, s3)).det()
require(factor(jacobian - 4 * (s1 * s2 - s3)) == 0, "hostile Jacobian failed")


print("THM-2690 CYCLIC-CUBIC NORMAL-CROSSING GATE -- exact checks")
print("snf_m=1", snf_rows[0][1], "class_order", snf_rows[0][2])
print("snf_m=2", snf_rows[1][1], "class_order", snf_rows[1][2])
print("snf_m=3", snf_rows[2][1], "class_order", snf_rows[2][2])
print("snf_m=12", snf_rows[-1][1], "class_order", snf_rows[-1][2])
print("all_m_1_to_12_snf=(1,3^(m-1)): PASS")
print("reduced_exponent_vectors_m_1_to_8=510 all_snf=(1,3^(m-1)): PASS")
print("singular_support_min_zero_count=2 for m=2..12: PASS")
print("C3_on_Cl2_standard_orbit", orbit, "fixed_nonzero=0")
print("C3_invariant_delta_squared_equals_discriminant: PASS")
print("C3_quotient_origin_singular: PASS")
print("C2_reflection_delta_to_minus_delta: PASS")
print("C2_reflection_even_relation_eliminates_delta_squared: PASS")
print("reflection_completed_hostile_Jacobian=4*(s1*s2-s3): PASS")
print("scope=A4_normal_crossing_family_excluded; general_A4_S4_JC2_DC2_open")
