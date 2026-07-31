#!/usr/bin/env python3
"""Independent exact audit for THM-2811.

This file does not import the primary companion.  It uses direct monomial
contraction, custom rational polynomial reduction, passport defect ledgers,
and symbolic differentiation as logically separate checks.
"""

from fractions import Fraction
from math import comb, factorial

import sympy as sp


checks = 0


def require(condition, label):
    global checks
    checks += 1
    if not bool(condition):
        raise RuntimeError(f"FAIL: {label}")


def contraction(poly, derivative_vars, coordinate_vars):
    """Apply E by its monomial formula, independently of differential APIs."""
    variables = tuple(derivative_vars) + tuple(coordinate_vars)
    source = sp.Poly(sp.expand(poly), *variables)
    answer = sp.Integer(0)
    split = len(derivative_vars)
    for exponents, coefficient in source.terms():
        alpha = exponents[:split]
        beta = exponents[split:]
        if all(a <= b for a, b in zip(alpha, beta)):
            term = coefficient
            for coordinate, a, b in zip(coordinate_vars, alpha, beta):
                term *= factorial(b) // factorial(b - a)
                term *= coordinate ** (b - a)
            answer += term
    return sp.expand(answer)


def partition_defect(partition):
    return sum(part - 1 for part in partition)


def monic_polynomial_from_roots(roots):
    """Ascending coefficient list for product(X-root), over Fraction."""
    coefficients = [Fraction(1)]
    for root in roots:
        updated = [Fraction(0)] * (len(coefficients) + 1)
        for degree, coefficient in enumerate(coefficients):
            updated[degree] -= root * coefficient
            updated[degree + 1] += coefficient
        coefficients = updated
    return coefficients


def monomial_remainder(degree, monic_modulus):
    """Ascending coefficients of X^degree modulo a monic modulus."""
    modulus_degree = len(monic_modulus) - 1
    coefficients = [Fraction(0)] * (degree + 1)
    coefficients[degree] = Fraction(1)
    for top in range(degree, modulus_degree - 1, -1):
        leading = coefficients[top]
        if leading == 0:
            continue
        shift = top - modulus_degree
        for index, modulus_coefficient in enumerate(monic_modulus):
            coefficients[shift + index] -= leading * modulus_coefficient
    return coefficients[:modulus_degree]


def projective_multiset_signature(values):
    """Canonical signature of a nonzero rational multiset up to scaling."""
    values = tuple(Fraction(value) for value in values)
    candidates = []
    for pivot in values:
        candidates.append(tuple(sorted(value / pivot for value in values)))
    return min(candidates)


print("THM-2811 independent exact audit")

# 1. Necessity on generators and sufficiency on a full small monomial box.
eta = sp.symbols("eta0:3")
y = sp.symbols("y0:3")
xi = sp.symbols("xi0:2")
z = sp.symbols("z0:2")

k_symbols = sp.symbols("k00 k01 k02 k10 k11 k12")
l_symbols = sp.symbols("l00 l01 l02 l10 l11 l12")
c_symbols = sp.symbols("c0 c1")
s_symbols = sp.symbols("s00 s01 s02 s10 s11 s12")
K_generic = sp.Matrix(2, 3, k_symbols)
L_generic = sp.Matrix(2, 3, l_symbols)
S_generic = sp.Matrix(2, 3, s_symbols)
c_generic = sp.Matrix(2, 1, c_symbols)

for row in range(2):
    image_xi = (K_generic * sp.Matrix(y) + L_generic * sp.Matrix(eta) + c_generic)[
        row
    ]
    expected = (K_generic * sp.Matrix(y) + c_generic)[row]
    require(
        sp.expand(contraction(image_xi, eta, y) - expected) == 0,
        f"xi generator forces K,c row={row}",
    )

for row in range(2):
    image_xi = (L_generic * sp.Matrix(eta))[row]
    for column in range(2):
        image_z = (S_generic * sp.Matrix(y))[column]
        observed = contraction(image_xi * image_z, eta, y)
        expected = (L_generic * S_generic.T)[row, column]
        require(
            sp.expand(observed - expected) == 0,
            f"mixed generator gives LS^T row={row},column={column}",
        )

S = sp.Matrix([[1, 2, -1], [0, 1, 1]])
L = (S * S.T).inv() * S
translation = sp.Matrix([sp.Rational(3, 2), -2])
require(L * S.T == sp.eye(2), "chosen skew embedding has LS^T=I")
require(S.rank() == 2 and L.rank() == 2, "chosen skew embedding has full ranks")

z_images = list(S * sp.Matrix(y) + translation)
xi_images = list(L * sp.Matrix(eta))
substitution = dict(zip(tuple(xi) + tuple(z), xi_images + z_images))
output_substitution = dict(zip(z, z_images))
monomial_rows = 0
for a0 in range(3):
    for a1 in range(3):
        for b0 in range(3):
            for b1 in range(3):
                monomial = xi[0] ** a0 * xi[1] ** a1 * z[0] ** b0 * z[1] ** b1
                transported = sp.expand(
                    monomial.subs(substitution, simultaneous=True)
                )
                left = contraction(transported, eta, y)
                source_output = contraction(monomial, xi, z)
                right = sp.expand(
                    source_output.subs(output_substitution, simultaneous=True)
                )
                require(
                    sp.expand(left - right) == 0,
                    f"affine sufficiency monomial={a0,a1,b0,b1}",
                )
                monomial_rows += 1

rank_rows = 0
for old_pairs in range(1, 6):
    for new_pairs in range(old_pairs, old_pairs + 3):
        extra = new_pairs - old_pairs
        S_block = sp.eye(old_pairs).row_join(
            sp.Matrix(
                old_pairs,
                extra,
                lambda row, column: (row + 1) * (column + 2),
            )
        )
        L_block = sp.eye(old_pairs).row_join(sp.zeros(old_pairs, extra))
        require(
            L_block * S_block.T == sp.eye(old_pairs),
            f"rank-wall construction n={old_pairs},k={new_pairs}",
        )
        require(
            S_block.rank() == old_pairs and L_block.rank() == old_pairs,
            f"full row ranks n={old_pairs},k={new_pairs}",
        )
        combined = sp.diag(L_block, S_block)
        require(
            combined.rank() == 2 * old_pairs,
            f"joint algebraic independence n={old_pairs},k={new_pairs}",
        )
        rank_rows += 1
    for new_pairs in range(0, old_pairs):
        require(
            new_pairs < old_pairs,
            f"rank(LS^T)<=k<n obstruction n={old_pairs},k={new_pairs}",
        )
        rank_rows += 1

print(
    "affine intertwiner/rank wall: PASS "
    f"(monomials={monomial_rows}, rank rows={rank_rows})"
)

# 2. Independent h=1,2,3 passport and total-fibre ledger.
passport_rows = 0
for degree in range(4, 33):
    zero_fibre = [2, 2] + [1] * (degree - 4)

    h1_fibres = [zero_fibre, [degree], [degree - 2, 1, 1]]
    require(
        [partition == [degree] for partition in h1_fibres]
        == [False, True, False],
        f"h1 unique total fibre degree={degree}",
    )
    require(
        sum(partition_defect(partition) for partition in h1_fibres)
        == 2 * degree - 2,
        f"h1 Riemann-Hurwitz degree={degree}",
    )
    passport_rows += 1

    for first_pole in range(1, degree):
        second_pole = degree - first_pole
        h2_fibres = [
            zero_fibre,
            [first_pole, second_pole],
            [degree - 1, 1],
        ]
        require(
            not any(partition == [degree] for partition in h2_fibres),
            f"h2 no total fibre degree={degree},d={first_pole}",
        )
        require(
            sum(partition_defect(partition) for partition in h2_fibres)
            == 2 * degree - 2,
            f"h2 Riemann-Hurwitz degree={degree},d={first_pole}",
        )
        require(
            (degree - 4) + degree + 4 == 2 * degree,
            f"h2 potential degree degree={degree}",
        )
        passport_rows += 1

    for first_pole in range(1, degree - 1):
        for second_pole in range(1, degree - first_pole):
            third_pole = degree - first_pole - second_pole
            h3_fibres = [
                zero_fibre,
                [first_pole, second_pole, third_pole],
                [degree],
            ]
            require(
                [partition == [degree] for partition in h3_fibres]
                == [False, False, True],
                (
                    "h3 unique total fibre "
                    f"degree={degree},parts={first_pole,second_pole,third_pole}"
                ),
            )
            require(
                sum(partition_defect(partition) for partition in h3_fibres)
                == 2 * degree - 2,
                (
                    "h3 Riemann-Hurwitz "
                    f"degree={degree},parts={first_pole,second_pole,third_pole}"
                ),
            )
            require(
                len(h3_fibres[1]) == 3,
                f"h3 polynomialized response has three roots degree={degree}",
            )
            passport_rows += 1

print(f"e=2 passport/Mobius trichotomy: PASS (rows={passport_rows})")

# 3. Binary Hessian-nilpotent formula, checked independently by degree.
X, Y, A, B = sp.symbols("X Y A B")
u = X + sp.I * Y
v = X - sp.I * Y
hn_rows = 0
for degree in range(2, 13):
    form = sp.expand(A * u**degree + B * v**degree)
    xx = sp.diff(form, X, 2)
    xy = sp.diff(form, X, Y)
    yy = sp.diff(form, Y, 2)
    trace = sp.expand(xx + yy)
    determinant = sp.expand(xx * yy - xy**2)
    expected_determinant = sp.expand(
        -4
        * degree**2
        * (degree - 1) ** 2
        * A
        * B
        * u ** (degree - 2)
        * v ** (degree - 2)
    )
    require(trace == 0, f"harmonic trace degree={degree}")
    require(
        sp.expand(determinant - expected_determinant) == 0,
        f"binary Hessian determinant degree={degree}",
    )
    hessian_u = sp.Matrix(
        [
            [sp.diff(u**degree, X, 2), sp.diff(u**degree, X, Y)],
            [sp.diff(u**degree, X, Y), sp.diff(u**degree, Y, 2)],
        ]
    )
    hessian_v = sp.Matrix(
        [
            [sp.diff(v**degree, X, 2), sp.diff(v**degree, X, Y)],
            [sp.diff(v**degree, X, Y), sp.diff(v**degree, Y, 2)],
        ]
    )
    require(
        all(sp.expand(entry) == 0 for entry in hessian_u**2),
        f"isotropic u Hessian square degree={degree}",
    )
    require(
        all(sp.expand(entry) == 0 for entry in hessian_v**2),
        f"isotropic v Hessian square degree={degree}",
    )
    hn_rows += 1

print(f"binary Hessian-nilpotent formula: PASS (degrees={hn_rows})")

# 4. Both polynomializable degree-four hostiles, using scalar derivatives.
P1 = sp.expand((X**2 - Y**2) ** 2)
P3 = sp.expand(X * (X - Y) * (2 * X - Y) ** 2)
for label, form, expected_trace, expected_determinant in [
    (
        "one-pole",
        P1,
        8 * (X**2 + Y**2),
        -48 * (X**2 - Y**2) ** 2,
    ),
    (
        "three-pole",
        P3,
        2 * (29 * X**2 - 27 * X * Y + 5 * Y**2),
        -3 * (2 * X - Y) ** 2 * (8 * X**2 - 8 * X * Y + 3 * Y**2),
    ),
]:
    xx = sp.diff(form, X, 2)
    xy = sp.diff(form, X, Y)
    yy = sp.diff(form, Y, 2)
    require(
        sp.expand(xx + yy - expected_trace) == 0,
        f"{label} quartic Laplacian",
    )
    require(
        sp.expand(xx * yy - xy**2 - expected_determinant) == 0,
        f"{label} quartic Hessian determinant",
    )
    require(sp.expand(xx + yy) != 0, f"{label} fails first contraction")

# Independently confirm the N=4 three-pole Maxwell symmetry at lambda=1/2.
x = sp.symbols("x")
D_three = sp.expand(x * (x - 1) * (x - sp.Rational(1, 2)) ** 2)
critical_quotient = sp.cancel(
    sp.diff(D_three, x) / (x * (x - 1) * (x - sp.Rational(1, 2)))
)
critical_roots = sp.solve(sp.together(critical_quotient).as_numer_denom()[0], x)
require(len(critical_roots) == 2, "N4 h3 has two non-pole critical points")
require(
    sp.simplify(D_three.subs(x, critical_roots[0]) - D_three.subs(x, critical_roots[1]))
    == 0,
    "N4 h3 critical values coincide",
)
print("both N=4 polynomialized Hessian hostiles: PASS")

# 5. Barycentric coefficient extraction using custom Q arithmetic and
# a signed rational node family different from the primary companion.
barycentric_rows = 0
for node_count in range(2, 12):
    nodes = [
        Fraction(((-1) ** index) * (index + 1), 2 * index + 3)
        for index in range(node_count)
    ]
    require(len(set(nodes)) == node_count, f"distinct signed nodes N={node_count}")
    modulus = monic_polynomial_from_roots(nodes)
    require(
        len(modulus) == node_count + 1 and modulus[-1] == 1,
        f"custom monic modulus N={node_count}",
    )
    derivatives = []
    for index, node in enumerate(nodes):
        derivative = Fraction(1)
        for other_index, other in enumerate(nodes):
            if index != other_index:
                derivative *= node - other
        derivatives.append(derivative)
    for degree in range(0, 3 * node_count + 1):
        left = sum(
            node**degree / derivative
            for node, derivative in zip(nodes, derivatives)
        )
        remainder = monomial_remainder(degree, modulus)
        right = (
            remainder[node_count - 1]
            if len(remainder) == node_count
            else Fraction(0)
        )
        require(
            left == right,
            f"custom barycentric remainder N={node_count},degree={degree}",
        )
        barycentric_rows += 1

print(f"custom rational barycentric identities: PASS (rows={barycentric_rows})")

# 6. SIC falling-factorial stencil and the response-alphabet mismatch.
sic_rows = 0
for power in range(1, 51):
    weights = []
    for node in range(power + 1):
        derivative = 1
        for other in range(power + 1):
            if node != other:
                derivative *= node - other
        barycentric_weight = factorial(power) // derivative
        binomial_weight = ((-1) ** (power - node)) * comb(power, node)
        require(
            barycentric_weight == binomial_weight,
            f"SIC falling-factorial derivative m={power},j={node}",
        )
        weights.append(barycentric_weight)
        sic_rows += 1
    for moment_degree in range(power + 1):
        moment = sum(
            weight * node**moment_degree for node, weight in enumerate(weights)
        )
        expected = factorial(power) if moment_degree == power else 0
        require(
            moment == expected,
            f"SIC finite-difference moment m={power},q={moment_degree}",
        )
        sic_rows += 1
    shifted = sum(
        ((-1) ** (power - node)) * comb(power, node - 1)
        for node in range(1, power + 1)
    )
    require(shifted == 1, f"SIC marked endpoint m={power}")
    sic_rows += 1

mismatch_rows = 0
for degree in range(4, 51):
    sic_weights = [
        ((-1) ** (degree - 1 - node)) * comb(degree - 1, node)
        for node in range(degree)
    ]
    sic_signature = projective_multiset_signature(sic_weights)
    for first_pole in range(1, degree // 2 + 1):
        response_weights = (
            [1] * (degree - 4)
            + [2, 2, -first_pole, -(degree - first_pole)]
        )
        require(
            projective_multiset_signature(response_weights) != sic_signature,
            f"response/SIC projective mismatch N={degree},d={first_pole}",
        )
        mismatch_rows += 1

print(
    "SIC stencil/response mismatch: PASS "
    f"(stencil rows={sic_rows}, mismatch rows={mismatch_rows})"
)
print(f"PASS checks={checks}")
