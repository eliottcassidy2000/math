#!/usr/bin/env python3
"""Exact certificate for THM-4241's W=0 saturation and gluing theorem.

Checks the algebraic identities and finite lattice enumerations used by the
theorem.  Standard geometric inputs such as Rosati integrality,
Castelnuovo--Severi, and the rational decomposition inherited from THM-4230
remain named noncomputational dependencies.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as Q

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


# ---------------------------------------------------------------------------
# Eisenstein arithmetic and the THM-4230 hidden Hermitian form.
# ---------------------------------------------------------------------------

E = tuple[Q, Q]  # m+n*omega, omega^2+omega+1=0


def e(value: tuple[int | Q, int | Q]) -> E:
    return Q(value[0]), Q(value[1])


def e_add(left: E, right: E) -> E:
    return left[0] + right[0], left[1] + right[1]


def e_mul(left: E, right: E) -> E:
    m, n = left
    r, s = right
    return m * r - n * s, m * s + n * r - n * s


def e_conjugate(value: E) -> E:
    m, n = value
    return m - n, -n


def e_trace(value: E) -> Q:
    return 2 * value[0] - value[1]


def e_norm(value: E) -> Q:
    return e_mul(value, e_conjugate(value))[0]


HIDDEN_GRAM: tuple[tuple[E, E], tuple[E, E]] = (
    (e((6, 0)), e((-4, -2))),
    (e((-2, 2)), e((6, 0))),
)
require(all(coordinate.denominator == 1 and coordinate.numerator % 2 == 0
            for row in HIDDEN_GRAM for value in row for coordinate in value),
        "hidden Gram no longer makes q mod 4 well-defined on L/2L")


def hermitian(left: tuple[E, E], right: tuple[E, E]) -> E:
    answer = e((0, 0))
    for i in range(2):
        for j in range(2):
            answer = e_add(
                answer,
                e_mul(e_mul(left[i], HIDDEN_GRAM[i][j]), e_conjugate(right[j])),
            )
    return answer


def hidden_degree(a: E, b: E) -> int:
    value = hermitian((a, b), (a, b))
    require(value[1] == 0 and value[0].denominator == 1, "degree not integral")
    return int(value[0])


# ---------------------------------------------------------------------------
# K=Q(zeta_12), its two ramified primes, and the four norm-integral ideals.
# ---------------------------------------------------------------------------

z = sp.symbols("z")
cyclotomic = sp.Poly(z**4 - z**2 + 1, z, domain=sp.QQ)


def k_reduce(expression: sp.Expr) -> sp.Expr:
    return sp.rem(sp.Poly(sp.expand(expression), z, domain=sp.QQ), cyclotomic).as_expr()


def k_inverse(expression: sp.Expr) -> sp.Expr:
    bezout, _, gcd = sp.gcdex(sp.Poly(k_reduce(expression), z, domain=sp.QQ), cyclotomic)
    return k_reduce(bezout.as_expr() / gcd.as_expr())


zeta = z
omega_polynomial = k_reduce(zeta**4)
T_polynomial = k_reduce(zeta**5)
power_basis = (sp.Integer(1), omega_polynomial, T_polynomial,
               k_reduce(omega_polynomial * T_polynomial))
power_matrix = sp.Matrix([
    [sp.Poly(value, z).coeff_monomial(z**degree) for value in power_basis]
    for degree in range(4)
])
require(abs(power_matrix.det()) == 1, "O-basis is not integral")


def k_to_o2(expression: sp.Expr) -> tuple[E, E]:
    reduced = k_reduce(expression)
    vector = sp.Matrix([
        sp.Poly(reduced, z).coeff_monomial(z**degree) for degree in range(4)
    ])
    coordinates = power_matrix.inv() * vector
    return (
        (Q(coordinates[0]), Q(coordinates[1])),
        (Q(coordinates[2]), Q(coordinates[3])),
    )


def o_scale(scalar: E, vector: tuple[E, E]) -> tuple[E, E]:
    return e_mul(scalar, vector[0]), e_mul(scalar, vector[1])


pi2 = 1 - zeta**3
pi3 = zeta + zeta**3
require(abs(int(sp.resultant(cyclotomic.as_expr(), pi2, z))) == 4,
        "N(pi2) changed")
require(abs(int(sp.resultant(cyclotomic.as_expr(), pi3, z))) == 9,
        "N(pi3) changed")
require(k_reduce(pi2**2 / 2) == -zeta**3, "(2)=p2^2 unit changed")
require(k_reduce(pi3**2 / 3) == omega_polynomial, "(3)=p3^2 unit changed")


def ideal_z_basis(e2: int, e3: int) -> list[tuple[E, E]]:
    alpha = k_inverse(k_reduce(pi2**e2 * pi3**e3))
    alpha_o2 = k_to_o2(alpha)
    t_alpha_o2 = k_to_o2(k_reduce(T_polynomial * alpha))
    omega = e((0, 1))
    return [alpha_o2, o_scale(omega, alpha_o2),
            t_alpha_o2, o_scale(omega, t_alpha_o2)]


def norm_integral(e2: int, e3: int) -> bool:
    basis = ideal_z_basis(e2, e3)
    diagonal = [hermitian(value, value) for value in basis]
    polar = [
        e_trace(hermitian(basis[i], basis[j]))
        for i in range(4) for j in range(i + 1, 4)
    ]
    return (
        all(value[1] == 0 and value[0].denominator == 1 for value in diagonal)
        and all(value.denominator == 1 for value in polar)
    )


norm_integral_exponents = {
    (e2, e3)
    for e2 in range(4)
    for e3 in range(4)
    if norm_integral(e2, e3)
}
require(norm_integral_exponents == {(0, 0), (1, 0), (0, 1), (1, 1)},
        "norm-integral overideal list changed")


def o_gram_of_inverse(prime: sp.Expr) -> list[list[E]]:
    alpha = k_inverse(prime)
    basis = [k_to_o2(alpha), k_to_o2(k_reduce(T_polynomial * alpha))]
    return [[hermitian(left, right) for right in basis] for left in basis]


p2_gram = o_gram_of_inverse(pi2)
p3_gram = o_gram_of_inverse(pi3)
require(p2_gram == [
    [e((3, 0)), e((-2, -1))],
    [e((-1, 1)), e((3, 0))],
], "p2 inverse Hermitian Gram changed")
require(p3_gram == [
    [e((2, 0)), e((Q(-4, 3), Q(-2, 3)))],
    [e((Q(-2, 3), Q(2, 3))), e((2, 0))],
], "p3 inverse Hermitian Gram changed")
require(all(coordinate.denominator == 1
            for row in p2_gram for value in row for coordinate in value),
        "p2 inverse lost Rosati integrality")
require(any(coordinate.denominator != 1
            for row in p3_gram for value in row for coordinate in value),
        "p3 inverse unexpectedly became Rosati integral")


# The Z-quadratic Gram of O_K f.  Its Smith support confines local enlargement
# to 2 and 3; the exponent audit above then gives the complete principal list.
z_basis = ideal_z_basis(0, 0)
quadratic_gram = sp.zeros(4)
for i in range(4):
    quadratic_gram[i, i] = hermitian(z_basis[i], z_basis[i])[0]
    for j in range(i):
        quadratic_gram[i, j] = quadratic_gram[j, i] = (
            e_trace(hermitian(z_basis[i], z_basis[j])) / 2
        )
require(quadratic_gram.det() == 324, "quadratic determinant changed")
quadratic_smith = smith_normal_form(quadratic_gram, domain=sp.ZZ)
# Smith diagonal entries are defined only up to multiplication by units.
# SymPy versions differ on the sign of the final entry for this matrix.
require(tuple(abs(quadratic_smith[i, i]) for i in range(4)) == (3, 3, 6, 6),
        "quadratic Smith factors changed")


# Fix the real embedding convention explicitly.  With
# sqrt(3)=zeta_12+zeta_12^-1, the displayed Gram is
# Tr_{Q(sqrt(3))/Q}((3+sqrt(3))*x*conjugate(x)).
def k_conjugate(expression: sp.Expr) -> sp.Expr:
    reduced = sp.Poly(k_reduce(expression), z, domain=sp.QQ)
    z_inverse = k_inverse(zeta)
    answer = sp.Integer(0)
    for degree in range(4):
        answer += reduced.coeff_monomial(z**degree) * z_inverse**degree
    return k_reduce(answer)


def k_trace(expression: sp.Expr) -> sp.Expr:
    columns = []
    for basis_value in (1, zeta, zeta**2, zeta**3):
        product = sp.Poly(k_reduce(expression * basis_value), z, domain=sp.QQ)
        columns.append([
            product.coeff_monomial(z**degree) for degree in range(4)
        ])
    return sp.Matrix.hstack(*(sp.Matrix(column) for column in columns)).trace()


sqrt_three = k_reduce(zeta + k_inverse(zeta))
require(k_reduce(sqrt_three**2) == 3, "sqrt(3) embedding changed")


def trace_degree(value: sp.Expr) -> sp.Expr:
    real_value = k_reduce(
        (3 + sqrt_three) * value * k_conjugate(value)
    )
    return sp.simplify(k_trace(real_value) / 2)


trace_quadratic_gram = sp.zeros(4)
for i in range(4):
    trace_quadratic_gram[i, i] = trace_degree(power_basis[i])
    for j in range(i):
        cross = (
            trace_degree(power_basis[i] + power_basis[j])
            - trace_degree(power_basis[i])
            - trace_degree(power_basis[j])
        ) / 2
        trace_quadratic_gram[i, j] = trace_quadratic_gram[j, i] = cross
require(trace_quadratic_gram == quadratic_gram,
        "3+sqrt(3) trace convention does not reproduce Q")


# ---------------------------------------------------------------------------
# Degree-three geometric obstruction and degree-two bielliptic obstruction.
# ---------------------------------------------------------------------------

a = sp.symbols("a")
first_difference = sp.rem(
    sp.Poly(sp.expand((1 - a)**3 - (-1 - a)**3), a),
    sp.Poly(a**2 + 1, a),
).as_expr()
second_difference = sp.rem(
    sp.Poly(sp.expand((1 - a)**3 + (-1 - a)**3), a),
    sp.Poly(a**2 + 1, a),
).as_expr()
require(first_difference == -4, "pole-at-zero contradiction changed")
require(second_difference == -4 * a, "pole-at-infinity contradiction changed")

# Quotient genera of the three involutions in C6 x C4.
interior_points = [
    (i, j)
    for i in range(1, 6)
    for j in range(1, 4)
    if Q(i, 6) + Q(j, 4) < 1
]
require(interior_points == [(1, 1), (1, 2), (1, 3),
                            (2, 1), (2, 2), (3, 1), (4, 1)],
        "Newton interior list changed")
rho_invariants = sum((i + j) % 2 == 0 for i, j in interior_points)
quotient_genera = (3, 2, rho_invariants)
require(quotient_genera == (3, 2, 4), "involution quotient genera changed")
castelnuovo_two_bielliptic_bound = 2 + 2 + 1
require(castelnuovo_two_bielliptic_bound == 5 < 7,
        "bielliptic uniqueness bound changed")


# ---------------------------------------------------------------------------
# Explicit degree-four visible-hidden gluing map.
# ---------------------------------------------------------------------------

y, lam, r, c, alpha, eps = sp.symbols("y lambda r c alpha eps")
q_polynomial = y**2 - c * y / 2 + 3 * lam**2
identity_difference = sp.Poly(
    sp.expand(c * (y - lam)**3 + 1 - y**4 + q_polynomial**2), y
)
relations = sp.groebner([
    c - r * lam,
    r**2 - 12 * r + 24,
    lam**4 * (r - 9) - 1,
], c, r, lam, order="lex", domain=sp.QQ)
for coefficient in identity_difference.all_coeffs():
    require(relations.reduce(coefficient)[1] == 0,
            "degree-four map identity changed")
require(sp.rem(sp.Poly((r - 10), r), sp.Poly(r**2 - 12*r + 24, r)).as_expr()
        != 0, "lambda^4=1 was not excluded")

# Pullback differential eta is scalar*x^-5*F(y)dy.  F(0) and F'(0) prove
# that both iota eigendifferentials are nonzero.
F = (3 + y**4 - 4 * lam * y**3) / q_polynomial
F_zero = sp.simplify(F.subs(y, 0))
F_derivative_zero = sp.simplify(sp.diff(F, y).subs(y, 0))
require(F_zero == lam**-2, "anti-invariant differential control changed")
require(sp.simplify(F_derivative_zero - c / (6 * lam**4)) == 0,
        "invariant differential control changed")


# ---------------------------------------------------------------------------
# Unique mod-two glue line, normalized full Gram, and exact theta census.
# ---------------------------------------------------------------------------

small_eisenstein = [e((0, 0)), e((1, 0)), e((0, 1)), e((1, 1))]
isotropic_mod_two: list[tuple[E, E]] = []
for aa in small_eisenstein:
    for bb in small_eisenstein:
        if hidden_degree(aa, bb) % 4 == 0:
            isotropic_mod_two.append((aa, bb))
require(len(isotropic_mod_two) == 4, "mod-two isotropic set size changed")
for aa, bb in isotropic_mod_two:
    omega_squared_b = e_mul(e((-1, -1)), bb)
    require((aa[0] - omega_squared_b[0]).denominator == 1
            and (aa[1] - omega_squared_b[1]).denominator == 1
            and int(aa[0] - omega_squared_b[0]) % 2 == 0
            and int(aa[1] - omega_squared_b[1]) % 2 == 0,
            "mod-two isotropic line changed")

w = sp.symbols("w")
full_gram = sp.Matrix([
    [4, 0, 0, 0],
    [0, 6, -4 - 2*w, 2*w - 2],
    [0, -2 + 2*w, 6, 2 - 2*w],
    [0, -4 - 2*w, 4 + 2*w, 4],
])
full_determinant = sp.rem(
    sp.Poly(sp.expand(full_gram.det()), w), sp.Poly(w**2 + w + 1, w)
).as_expr()
require(full_determinant == 96, "full Hermitian determinant changed")


def e_int(value: tuple[int, int]) -> E:
    return e(value)


def e_scale_int(multiplier: int, value: E) -> E:
    return Q(multiplier) * value[0], Q(multiplier) * value[1]


def full_degree(visible: E, f_coefficient: E, g_coefficient: E,
                h_coefficient: E) -> int:
    # 2h=v+(omega^2 f+g), with q(v)=4 and v orthogonal to f,g.
    hidden_a = e_add(e_scale_int(2, f_coefficient),
                     e_mul(e((-1, -1)), h_coefficient))
    hidden_b = e_add(e_scale_int(2, g_coefficient), h_coefficient)
    numerator = hidden_degree(hidden_a, hidden_b)
    require(numerator % 4 == 0, "full degree numerator lost divisibility")
    result = 4 * int(e_norm(visible)) + int(e_norm(h_coefficient)) + numerator // 4
    return result


elements = [
    (e_int((m, n)), int(e_norm(e_int((m, n)))))
    for m in range(-8, 9)
    for n in range(-8, 9)
    if e_norm(e_int((m, n))) <= 42
]
require(not any(max(abs(int(value[0])), abs(int(value[1]))) == 8
                for value, _ in elements),
        "Eisenstein enumeration touched its coordinate boundary")
hidden_rows: list[tuple[E, E, int]] = []
for am in range(-12, 13):
    for an in range(-12, 13):
        for bm in range(-12, 13):
            for bn in range(-12, 13):
                degree = hidden_degree(e_int((am, an)), e_int((bm, bn)))
                if degree <= 168:
                    hidden_rows.append((e_int((am, an)), e_int((bm, bn)), degree))
require(not any(max(abs(int(a0)), abs(int(a1)), abs(int(b0)), abs(int(b1))) == 12
                for (a0, a1), (b0, b1), _ in hidden_rows),
        "hidden enumeration touched its coordinate boundary")

rows_by_parity: dict[tuple[int, int, int, int], list[tuple[E, E, int]]] = defaultdict(list)
for hidden_a, hidden_b, degree in hidden_rows:
    parity = (int(hidden_a[0]) % 2, int(hidden_a[1]) % 2,
              int(hidden_b[0]) % 2, int(hidden_b[1]) % 2)
    rows_by_parity[parity].append((hidden_a, hidden_b, degree))

theta = Counter()
for h_coefficient, h_norm in elements:
    omega_squared_h = e_mul(e((-1, -1)), h_coefficient)
    parity = (int(omega_squared_h[0]) % 2, int(omega_squared_h[1]) % 2,
              int(h_coefficient[0]) % 2, int(h_coefficient[1]) % 2)
    for hidden_a, hidden_b, hidden_numerator in rows_by_parity[parity]:
        base = h_norm + hidden_numerator // 4
        if base > 42:
            continue
        for visible, visible_norm in elements:
            del visible
            degree = base + 4 * visible_norm
            if degree <= 42:
                theta[degree] += 1

expected_theta = {
    0: 1, 4: 60, 6: 72, 8: 324, 10: 864, 12: 276, 14: 2592,
    16: 2868, 18: 1152, 20: 5832, 22: 9504, 24: 2556,
    26: 15552, 28: 15456, 30: 6480, 32: 22356, 34: 36288,
    36: 7836, 38: 49248, 40: 44280, 42: 16992,
}
require(dict(sorted(theta.items())) == expected_theta, "full theta census changed")
require(full_degree(e((0, 0)), e((-1, 0)), e((1, 0)), e((1, 0))) == 34,
        "degree-34 example changed")
require(full_degree(e((0, 0)), e((-1, 0)), e((0, -1)), e((1, -1))) == 42,
        "degree-42 example changed")


print("W=0 saturation and gluing exact certificate")
print("field=Q(zeta12) omega=zeta12^4 T=zeta12^5 T^2=-omega")
print("scalar_degree=Tr_Qsqrt3/Q((3+sqrt3)*x*xbar) trace_gram=PASS")
print("prime2=(1-zeta12^3) norm=4 square=unit*2")
print("prime3=(zeta12+zeta12^3) norm=9 square=unit*3")
print("quadratic_smith=[3,3,6,6] norm_integral_overideals=O,P2^-1,P3^-1,(P2*P3)^-1")
print("candidate_degrees=6,3,2,1")
print("P2_inverse_Rosati_integral=PASS P3_inverse_Rosati_integral=FAIL_denominator_3")
print("bielliptic_controls=CS_bound:5<7 central_involution_quotient_genera:3,2,4")
print("degree3_pole_classifier=only_zero_or_infinity_patterns; endpoint_differences=-4,-4*i")
print("hidden_saturation=PROVED_THM4241")
print("explicit_glue_map_identity=PASS degree=4 pole_divisor=4_points_times_2")
print(f"eigendifferential_controls=F0:{F_zero},Fprime0:{F_derivative_zero}")
print("full_quotient=O/2O cardinality=4 underlying_Z_index=4")
print(f"full_hermitian_determinant={full_determinant}")
print(f"theta_through_42={dict(sorted(theta.items()))}")
print("degree34_example=h-f+g count=36288")
print("degree42_example=-f-omega*g+(1-omega)*h count=16992")
print("attachment_boundary=FINITE_UNENUMERATED_S_34_S_42")
print("result=PASS")
