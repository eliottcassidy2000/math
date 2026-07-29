#!/usr/bin/env python3
"""Exact companion for THM-2843 four-slot projective resolvent reduction.

This script verifies the algebraic construction and the explicit
support-{0,1,2,3} cell.  It deliberately does not claim the full 378-cell
SFC(4) atlas.
"""

from itertools import product
from math import comb, factorial

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


S = sp.symbols("S")
X0, X1, X2, X3 = sp.symbols("X0 X1 X2 X3")
COEFFICIENTS = (X0, X1, X2, X3)


def factorial_readout(polynomial):
    """Apply L(S^n)=n! exactly."""

    expanded = sp.Poly(sp.expand(polynomial), S)
    return sp.expand(
        sum(
            coefficient * factorial(exponent[0])
            for exponent, coefficient in expanded.terms()
        )
    )


def weak_compositions(total, slots):
    if slots == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, slots - 1):
            yield (first,) + tail


def normalized_moment_formula(support, power):
    """Coefficient formula for H=sum X_j S^a_j/a_j!."""

    result = 0
    for multiplicities in weak_compositions(power, len(support)):
        multinomial = factorial(power)
        weighted_degree = 0
        monomial = 1
        denominator = 1
        for index, multiplicity in enumerate(multiplicities):
            multinomial //= factorial(multiplicity)
            weighted_degree += multiplicity * support[index]
            monomial *= COEFFICIENTS[index] ** multiplicity
            denominator *= factorial(support[index]) ** multiplicity
        result += (
            sp.Rational(multinomial * factorial(weighted_degree), denominator)
            * monomial
        )
    return sp.expand(result)


def normalized_moment_direct(support, power):
    polynomial = sum(
        COEFFICIENTS[index] * S**degree / factorial(degree)
        for index, degree in enumerate(support)
    )
    return factorial_readout(polynomial**power)


# General construction controls on several nonsymmetric supports and all
# moment degrees needed by the k<=2 four-slot window.
support_bank = (
    (0, 1, 2, 3),
    (0, 1, 3, 5),
    (0, 2, 5, 7),
    (1, 3, 4, 6),
)
moment_formula_cells = 0
positive_gram_cells = 0
for support in support_bank:
    for power in range(1, 7):
        formula = normalized_moment_formula(support, power)
        direct = normalized_moment_direct(support, power)
        require(formula == direct, f"moment formula {support=} {power=}")
        if power == 1:
            require(
                formula == sum(COEFFICIENTS),
                f"normalized first moment {support=}",
            )
        moment_formula_cells += 1

    quadratic = normalized_moment_formula(support, 2)
    gram = sp.hessian(quadratic, COEFFICIENTS) / 2
    for size in range(1, 5):
        require(
            sp.det(gram[:size, :size]) > 0,
            f"positive factorial Gram principal minor {support=} {size=}",
        )
    positive_gram_cells += 1


# Explicit support {0,1,2,3}, after eliminating M_1 by adjacent differences.
x, y, z = sp.symbols("x y z")
divided_powers = tuple(S**degree / factorial(degree) for degree in range(4))
adjacent = tuple(
    divided_powers[index + 1] - divided_powers[index] for index in range(3)
)
H = x * adjacent[0] + y * adjacent[1] + z * adjacent[2]
Q = factorial_readout(H**2)
C = factorial_readout(H**3)
F = factorial_readout(H**4)

expected_Q = x**2 + 2 * x * y + 2 * x * z + 2 * y**2 + 6 * y * z + 6 * z**2
sum_of_squares_Q = (x + y + z) ** 2 + (y + 2 * z) ** 2 + z**2
require(sp.expand(Q - expected_Q) == 0, "explicit quadratic")
require(sp.expand(Q - sum_of_squares_Q) == 0, "quadratic sum of squares")
require(sp.det(sp.hessian(Q, (x, y, z)) / 2) == 1, "quadratic Gram determinant")


# Parametrize X^2+Y^2+Z^2=0, where
# X=x+y+z, Y=y+2z, Z=z.
u, v = sp.symbols("u v")
I = sp.I
capital_x = u**2 - v**2
capital_y = I * (u**2 + v**2)
capital_z = 2 * u * v
conic_substitution = {
    z: capital_z,
    y: capital_y - 2 * capital_z,
    x: capital_x - capital_y + capital_z,
}
require(
    sp.expand(Q.subs(conic_substitution)) == 0,
    "conic parametrization",
)

c6 = sp.expand(C.subs(conic_substitution) / 2)
f8 = sp.expand(F.subs(conic_substitution) / 3)
expected_c6 = (
    (-5 - 2 * I) * u**6
    + (-36 + 18 * I) * u**5 * v
    + (27 + 126 * I) * u**4 * v**2
    + 152 * u**3 * v**3
    + (-27 + 126 * I) * u**2 * v**4
    + (-36 - 18 * I) * u * v**5
    + (5 - 2 * I) * v**6
)
expected_f8 = (
    (44 - 96 * I) * u**8
    + (-992 - 1264 * I) * u**7 * v
    + (-10368 + 4160 * I) * u**6 * v**2
    + (6624 + 39216 * I) * u**5 * v**3
    + 51048 * u**4 * v**4
    + (-6624 + 39216 * I) * u**3 * v**5
    + (-10368 - 4160 * I) * u**2 * v**6
    + (992 - 1264 * I) * u * v**7
    + (44 + 96 * I) * v**8
)
require(sp.expand(c6 - expected_c6) == 0, "binary sextic")
require(sp.expand(f8 - expected_f8) == 0, "binary octic")

c6_affine = sp.Poly(c6.subs(v, 1), u, domain=sp.QQ_I)
f8_affine = sp.Poly(f8.subs(v, 1), u, domain=sp.QQ_I)
require(c6_affine.degree() == 6, "sextic degree")
require(f8_affine.degree() == 8, "octic degree")
require(
    sp.gcd(c6_affine, sp.Poly(sp.diff(c6_affine.as_expr(), u), u, domain=sp.QQ_I))
    == 1,
    "sextic squarefree",
)
require(sp.gcd(c6_affine, f8_affine) == 1, "sextic-octic coprime")

explicit_resultant = sp.resultant(
    c6_affine.as_expr(),
    f8_affine.as_expr(),
    u,
)
expected_resultant = sp.Integer(208741470184115575361509867388928)
expected_factorization = {
    sp.Integer(2): 36,
    sp.Integer(3): 7,
    sp.Integer(67): 1,
    sp.Integer(11702701): 1,
    sp.Integer(1771410437): 1,
}
require(explicit_resultant == expected_resultant, "explicit binary resultant")
require(
    sp.factorint(explicit_resultant) == expected_factorization,
    "resultant factorization",
)

sextic_discriminant = sp.discriminant(c6_affine.as_expr(), u)
expected_discriminant = sp.Integer(6291863054003994624)
require(sextic_discriminant == expected_discriminant, "sextic discriminant")
require(
    sp.factorint(sextic_discriminant)
    == {
        sp.Integer(2): 20,
        sp.Integer(3): 12,
        sp.Integer(31): 3,
        sp.Integer(379): 1,
    },
    "sextic discriminant factorization",
)


# The affine z=1 quotient has dimension six.  Its fourth-moment norm can be
# obtained from the shape Groebner basis and a univariate resultant.
groebner_basis = sp.groebner(
    [Q.subs(z, 1), C.subs(z, 1)],
    y,
    x,
    order="lex",
)
require(len(groebner_basis.polys) == 2, "shape Groebner basis size")
shape_linear = sp.Poly(groebner_basis.polys[0].as_expr(), y, x)
shape_sextic = sp.Poly(groebner_basis.polys[1].as_expr(), x)
require(shape_linear.degree(y) == 1, "shape basis linear y equation")
require(shape_sextic.degree() == 6, "shape basis quotient dimension")
linear_coefficient = shape_linear.coeff_monomial(y)
linear_remainder = sp.expand(shape_linear.as_expr() - linear_coefficient * y)
require(not linear_remainder.has(y), "shape linear remainder")

fourth_numerator = sp.Poly(
    sp.cancel(
        F.subs(z, 1).subs(y, -linear_remainder / linear_coefficient)
        * linear_coefficient**4
    ),
    x,
)
shape_resultant = sp.resultant(
    shape_sextic.as_expr(),
    fourth_numerator.as_expr(),
    x,
)
quotient_norm = sp.cancel(
    shape_resultant
    / (
        shape_sextic.LC() ** fourth_numerator.degree()
        * linear_coefficient ** (4 * 6)
    )
)
expected_quotient_norm = sp.Rational(
    9070189700378194715889733632,
    707281,
)
require(quotient_norm == expected_quotient_norm, "six-dimensional quotient norm")
require(quotient_norm > 0, "positive explicit quotient norm")

shape_sextic_qq = sp.Poly(shape_sextic.as_expr(), x, domain=sp.QQ)
fourth_reduced_qq = sp.Poly(
    fourth_numerator.as_expr() / linear_coefficient**4,
    x,
    domain=sp.QQ,
)
norm_columns = []
for degree in range(6):
    remainder = (
        fourth_reduced_qq * sp.Poly(x**degree, x, domain=sp.QQ)
    ).rem(shape_sextic_qq)
    norm_columns.append([remainder.nth(row) for row in range(6)])
direct_norm_matrix = sp.Matrix(
    6,
    6,
    lambda row, column: norm_columns[column][row],
)
require(
    sp.det(direct_norm_matrix) == quotient_norm,
    "direct six-dimensional multiplication determinant",
)

face_resultant = sp.factor(
    sp.resultant(Q.subs(z, 0), C.subs(z, 0), x)
)
require(face_resultant == 116 * y**6, "three-slot dehomogenizing face")


# Exact moving-plane divisibility and the sharp abstract facewise hostile.
hidden_Q = x**2 + y**2 + z**2
hidden_C = 2 * x**3 + 2 * y**3 + (x + y) * z**2
hidden_F = (x - y) * x**3
plane_lambda = x - y
hidden_B1 = x + y
hidden_A2 = (x + y) * (x - y)
hidden_B2 = 0
hidden_A3 = x**3
require(
    sp.expand(hidden_C - hidden_Q * hidden_B1 - plane_lambda * hidden_A2) == 0,
    "cubic plane-ideal divisibility",
)
require(
    sp.expand(hidden_F - hidden_Q * hidden_B2 - plane_lambda * hidden_A3) == 0,
    "quartic plane-ideal divisibility",
)
require(
    sp.expand(hidden_C.subs(y, x) - 2 * x * hidden_Q.subs(y, x)) == 0,
    "binary q divides cubic",
)
require(hidden_F.subs(y, x) == 0, "binary q divides quartic")

hidden_point = {x: 1, y: 1, z: I * sp.sqrt(2)}
require(sp.simplify(hidden_Q.subs(hidden_point)) == 0, "hidden point quadratic")
require(sp.simplify(hidden_C.subs(hidden_point)) == 0, "hidden point cubic")
require(sp.simplify(hidden_F.subs(hidden_point)) == 0, "hidden point quartic")

coordinate_face_resultants = (
    sp.factor(sp.resultant(hidden_Q.subs(z, 0), hidden_C.subs(z, 0), x)),
    sp.factor(sp.resultant(hidden_Q.subs(y, 0), hidden_C.subs(y, 0), x)),
    sp.factor(sp.resultant(hidden_Q.subs(x, 0), hidden_C.subs(x, 0), y)),
)
require(
    coordinate_face_resultants == (8 * y**6, z**6, z**6),
    "coordinate faces miss hidden common point",
)


# The no-real-points hypothesis is exactly what makes the real norm positive.
a_symbol, b_symbol = sp.symbols("a b", real=True)
complex_pair_matrix = sp.Matrix(
    [
        [a_symbol, -b_symbol],
        [b_symbol, a_symbol],
    ]
)
require(
    sp.det(complex_pair_matrix) == a_symbol**2 + b_symbol**2,
    "conjugate-pair positive norm",
)
nonreduced_pair_matrix = sp.Matrix(
    [
        [0, 0, 0, -1],
        [1, 0, 0, 0],
        [0, 1, 0, -2],
        [0, 0, 1, 0],
    ]
)
require(sp.det(nonreduced_pair_matrix) == 1, "nonreduced conjugate-pair norm")
real_point_matrix = sp.Matrix([[0, 1], [1, 0]])
require(sp.det(real_point_matrix) == -1, "real-point sign hostile")


# Exact failure of the standard PSL_2(Z)=C_2*C_3 generators to preserve c6.
sigma_matrix = sp.Matrix([[0, -1], [1, 0]])
tau_matrix = sp.Matrix([[0, -1], [1, 1]])
require(sigma_matrix**2 == -sp.eye(2), "sigma projective order two")
require(tau_matrix**3 == -sp.eye(2), "tau projective order three")

sigma_c6 = sp.expand(c6.subs({u: -v, v: u}, simultaneous=True))
tau_c6 = sp.expand(c6.subs({u: -v, v: u + v}, simultaneous=True))
modular_lambda = sp.Rational(-21, 29) + sp.Rational(20, 29) * I
sigma_defect = sp.Poly(
    sp.expand(sigma_c6 - modular_lambda * c6),
    u,
    v,
)
tau_defect = sp.Poly(
    sp.expand(tau_c6 - modular_lambda * c6),
    u,
    v,
)
expected_sigma_defect = sp.Rational(648, 29) + sp.Rational(1620, 29) * I
expected_tau_defect = sp.Rational(1518, 29) + sp.Rational(1272, 29) * I
require(
    sigma_defect.coeff_monomial(u**5 * v) == expected_sigma_defect,
    "sigma sextic nonaction",
)
require(
    tau_defect.coeff_monomial(u**5 * v) == expected_tau_defect,
    "tau sextic nonaction",
)


# Optional finite base-field control: the same sextic has trivial stabilizer
# in PGL_2(F_17).  This is a finite control, not a proof about every geometric
# automorphism in characteristic zero.
def normalize_matrix(matrix_entries, prime):
    for entry in matrix_entries:
        if entry % prime:
            inverse = pow(entry, -1, prime)
            return tuple(value * inverse % prime for value in matrix_entries)
    raise RuntimeError("zero projective matrix")


def transform_binary_coefficients(coefficients, matrix_entries, prime):
    """Coefficients are ordered U^n, U^(n-1)V, ..., V^n."""

    aa, bb, cc, dd = matrix_entries
    degree = len(coefficients) - 1
    transformed = [0] * (degree + 1)
    for second_degree, coefficient in enumerate(coefficients):
        first_degree = degree - second_degree
        for first_v_degree in range(first_degree + 1):
            first_term = (
                comb(first_degree, first_v_degree)
                * pow(aa, first_degree - first_v_degree, prime)
                * pow(bb, first_v_degree, prime)
            )
            for second_v_degree in range(second_degree + 1):
                second_term = (
                    comb(second_degree, second_v_degree)
                    * pow(cc, second_degree - second_v_degree, prime)
                    * pow(dd, second_v_degree, prime)
                )
                target_degree = first_v_degree + second_v_degree
                transformed[target_degree] = (
                    transformed[target_degree]
                    + coefficient * first_term * second_term
                ) % prime
    return tuple(transformed)


prime = 17
square_root_minus_one = 4
gaussian_coefficients = (
    (-5, -2),
    (-36, 18),
    (27, 126),
    (152, 0),
    (-27, 126),
    (-36, -18),
    (5, -2),
)
finite_coefficients = tuple(
    (real + imaginary * square_root_minus_one) % prime
    for real, imaginary in gaussian_coefficients
)
projective_matrices = set()
for aa, bb, cc, dd in product(range(prime), repeat=4):
    if (aa * dd - bb * cc) % prime:
        projective_matrices.add(normalize_matrix((aa, bb, cc, dd), prime))
require(len(projective_matrices) == 4896, "PGL2(F17) order")

finite_stabilizer = []
first_nonzero_index = next(
    index for index, coefficient in enumerate(finite_coefficients) if coefficient
)
for matrix_entries in projective_matrices:
    transformed = transform_binary_coefficients(
        finite_coefficients,
        matrix_entries,
        prime,
    )
    scalar = (
        transformed[first_nonzero_index]
        * pow(finite_coefficients[first_nonzero_index], -1, prime)
    ) % prime
    if all(
        (transformed[index] - scalar * finite_coefficients[index]) % prime == 0
        for index in range(7)
    ):
        finite_stabilizer.append(matrix_entries + (scalar,))
finite_stabilizer.sort()
require(
    finite_stabilizer == [(1, 0, 0, 1, 1)],
    "trivial PGL2(F17) sextic stabilizer",
)


print("THM-2843 FOUR-SLOT PROJECTIVE RESOLVENT EXACT COMPANION")
print("status=STRUCTURAL-REDUCTION CONTROLS; NOT FULL SFC4 ATLAS")
print(f"moment_formula_cells={moment_formula_cells}")
print(f"positive_gram_cells={positive_gram_cells}")
print("explicit_support=0,1,2,3")
print("quadratic_gram_determinant=1")
print("quotient_dimension=6")
print(f"quotient_norm={quotient_norm}")
print(f"binary_resultant={explicit_resultant}")
print("binary_resultant_factorization=2^36*3^7*67*11702701*1771410437")
print(f"sextic_discriminant={sextic_discriminant}")
print("moving_plane_divisibility=PASS")
print("abstract_hidden_plane_hostile=PASS")
print("conjugate_pair_norm_controls=PASS")
print(f"sigma_u5v_defect={expected_sigma_defect}")
print(f"tau_u5v_defect={expected_tau_defect}")
print(f"pgl2_f17_order={len(projective_matrices)}")
print(f"pgl2_f17_sextic_stabilizer={len(finite_stabilizer)}")
print("full_378_cell_atlas_claimed=NO")
print("all_exact_controls=PASS")
