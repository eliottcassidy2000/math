#!/usr/bin/env python3
"""Exact hostile probe for prescribed-coordinate orders at the newest Keller prime.

The script has two independent lanes.

1.  Over Q it uses the rational point q0=(2/27,1,1) on V(L), factors the
    other two points in F^{-1}(F(q0)) as one quadratic algebra, and computes
    the three last-step coordinate blocks by exact resultants.  For each of
    y, z, and u=1/x, the degree-nine special-fibre polynomial has exactly one
    repeated linear factor: the double shadow on the unique L=0 block.  Thus
    no cross-block residue collision is forced at the level-two newest prime.

2.  Over finite fields it searches for complete split inverse trees rooted
    at a point q0 on V(L).  It then multiplies all last-step y-, z-, and
    reciprocal-x cubics and checks that gcd(P,P') has degree exactly one.
    These are finite exact witnesses against an identically vanishing
    cross-block resultant at the corresponding characteristic-zero divisor.

The geometric all-level reduction and its scope belong in the companion
reflection/theorem candidate.  Finite-field witnesses are not promoted here
to an all-level noncollision theorem.
"""

from __future__ import annotations

from hashlib import sha256
import json

import sympy as sp


EXPECTED_SEMANTIC_SHA256 = "54e24e2b22edbc0882be21fc1c64f8b451a6c8767af779be52a3e16adcb6db34"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def l_value(a, b, c):
    return 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2


def f_map(point):
    x, y, z = point
    unit = 1 + x * y
    return (
        unit**3 * z + y**2 * unit * (4 + 3 * x * y),
        y + 3 * x * unit**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
        2 * x - 3 * x**2 * y - x**3 * z,
    )


def inverse_section(target, x):
    """The rational inverse section from THM-3519, with x left algebraic."""
    a, b, c = target
    denominator = (12 * a - b**2) * x**2 + b * x + 2
    y = b - 3 * a * x * ((9 * a * c - b) * x + 2) / denominator
    z = (2 * x - c - 3 * x**2 * y) / x**3
    return x, sp.cancel(y), sp.cancel(z)


def y_block(target, variable):
    """Constant-leading cubic for the last-step source y-coordinate."""
    a, b, c = target
    r = b - variable
    return sp.expand(27 * a**2 * c - 18 * a * r + 3 * b * r**2 - 2 * r**3)


def z_block(target, variable):
    """Constant-leading cubic Q_z of THM-2546."""
    a, b, c = target
    L = l_value(a, b, c)
    q2 = 324 * a**2 * c**2 - 216 * a * b * c + 408 * a - 15 * b**3 * c + 6 * b**2
    Sz = 27 * a**2 * c**2 - 18 * a * b * c + 52 * a + b**3 * c + 14 * b**2
    Tz = (
        729 * a**4 * c**4
        - 972 * a**3 * b * c**3
        + 2322 * a**3 * c**2
        + 54 * a**2 * b**3 * c**3
        + 270 * a**2 * b**2 * c**2
        - 3735 * a**2 * b * c
        - 338 * a**2
        - 36 * a * b**4 * c**2
        + 122 * a * b**3 * c
        + 1372 * a * b**2
        + b**6 * c**2
        - 2 * b**5 * c
        - 80 * b**4
    )
    return sp.expand(8 * variable**3 + q2 * variable**2 + 6 * L * Sz * variable + L * Tz)


def reciprocal_x_block(target, variable):
    """Integral reciprocal cubic L+(4-3bc)u^2-2cu^3."""
    a, b, c = target
    return sp.expand(l_value(a, b, c) + (4 - 3 * b * c) * variable**2 - 2 * c * variable**3)


def primitive_qq_poly(expression, variable):
    polynomial = sp.Poly(sp.cancel(expression), variable, domain=sp.QQ)
    return polynomial.monic()


def rational_level_two_witness():
    X, Y, Z, U = sp.symbols("X Y Z U")
    q0 = (sp.Rational(2, 27), sp.Integer(1), sp.Integer(1))
    require(l_value(*q0) == 0, "q0 must lie on L")
    q0_square_factors = {
        "h_y": sp.Rational(27, 2) * q0[0],
        "h_z": sp.Rational(27, 32) * z_square_factor_value(q0),
        "h_u": sp.Rational(27 * q0[0] * q0[2] ** 2 - 9 * q0[1] * q0[2] + 8,
                            2 * q0[2] ** 2),
    }
    require(q0_square_factors == {"h_y": 1, "h_z": sp.Rational(49, 32), "h_u": sp.Rational(1, 2)},
            q0_square_factors)
    eta = tuple(sp.factor(value) for value in f_map(q0))
    outer = sp.factor(l_value(*eta) * X**3 + (4 - 3 * eta[1] * eta[2]) * X - 2 * eta[2])
    quotient, remainder = sp.div(sp.Poly(outer, X), sp.Poly(X - q0[0], X))
    require(remainder.is_zero and quotient.degree() == 2, (outer, quotient, remainder))
    h = quotient.monic().as_expr()

    qx = inverse_section(eta, X)
    require(all(sp.cancel(qx[i].subs(X, q0[0]) - q0[i]) == 0 for i in range(3)), qx)

    rows = {}
    for name, variable, constructor, expected_gcd in (
        ("y", Y, y_block, Y - sp.Rational(1, 3)),
        ("z", Z, z_block, Z),
        ("u", U, reciprocal_x_block, U),
    ):
        block0 = primitive_qq_poly(constructor(q0, variable), variable)
        generic_numerator = sp.together(constructor(qx, variable)).as_numer_denom()[0]
        other_product_raw = sp.resultant(h, generic_numerator, X)
        other_product = primitive_qq_poly(other_product_raw, variable)
        total = primitive_qq_poly(block0.as_expr() * other_product.as_expr(), variable)
        gcd = sp.gcd(total, total.diff()).monic()
        require(total.degree() == 9, (name, "degree", total.degree()))
        require(gcd == sp.Poly(expected_gcd, variable, domain=sp.QQ).monic(), (name, gcd))
        require(sp.gcd(block0, other_product).degree() == 0, (name, "cross-block collision"))
        require(sp.gcd(other_product, other_product.diff()).degree() == 0, (name, "other blocks collide"))
        rows[name] = {
            "block0_factorization": str(sp.factor(block0.as_expr())),
            "other_degree": other_product.degree(),
            "total_degree": total.degree(),
            "derivative_gcd": str(gcd.as_expr()),
            "coefficient_digest": digest([str(c) for c in total.all_coeffs()]),
        }

    return {
        "q0": tuple(str(v) for v in q0),
        "eta": tuple(str(v) for v in eta),
        "outer_factorization": str(outer),
        "quadratic_factor": str(sp.factor(h)),
        "ramified_block_square_factors": {name: str(value) for name, value in q0_square_factors.items()},
        "rows": rows,
    }


def f_mod(point, prime):
    x, y, z = point
    unit = (1 + x * y) % prime
    return (
        (unit**3 * z + y**2 * unit * (4 + 3 * x * y)) % prime,
        (y + 3 * x * unit**2 * z + 3 * x * y**2 * (4 + 3 * x * y)) % prime,
        (2 * x - 3 * x**2 * y - x**3 * z) % prime,
    )


def l_mod(point, prime):
    a, b, c = point
    return (27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2) % prime


def s_mod(point, prime):
    a, b, c = point
    return (27 * a * c**2 - 9 * b * c + 8) % prime


def z_square_factor_value(point):
    a, b, c = point
    return (
        729 * a**3 * b**3 * c**4
        + 4374 * a**3 * b**2 * c**3
        + 8748 * a**3 * b * c**2
        + 5832 * a**3 * c
        - 729 * a**2 * b**4 * c**3
        - 2754 * a**2 * b**3 * c**2
        - 2268 * a**2 * b**2 * c
        + 648 * a**2 * b
        + 27 * a * b**6 * c**3
        + 297 * a * b**5 * c**2
        + 360 * a * b**4 * c
        - 268 * a * b**3
        - 9 * b**7 * c**2
        - 28 * b**6 * c
        + 28 * b**5
    )


def z_square_factor_mod(point, prime):
    return z_square_factor_value(point) % prime


def l_smooth_mod(point, prime):
    a, b, c = point
    gradient = (
        54 * a * c**2 - 18 * b * c + 16,
        -18 * a * c + 3 * b**2 * c - 2 * b,
        54 * a**2 * c - 18 * a * b + b**3,
    )
    return any(value % prime for value in gradient)


def poly_mod(expression, variable, prime):
    return sp.Poly(expression, variable, modulus=prime).monic()


def resultant_mod(left, right_expression, variable, prime):
    right = sp.Poly(right_expression, variable, modulus=prime)
    return int(left.resultant(right)) % prime


def interpolate_mod(values, variable, prime):
    """Return the degree < len(values) interpolant over F_prime.

    ``values[j]`` is the value at j.  Coefficients are stored low first while
    the Lagrange basis is assembled; degrees here are at most nine.
    """
    count = len(values)
    answer = [0] * count
    for i, value in enumerate(values):
        basis = [1]
        denominator = 1
        for j in range(count):
            if i == j:
                continue
            new_basis = [0] * (len(basis) + 1)
            for degree, coefficient in enumerate(basis):
                new_basis[degree] = (new_basis[degree] - j * coefficient) % prime
                new_basis[degree + 1] = (new_basis[degree + 1] + coefficient) % prime
            basis = new_basis
            denominator = denominator * (i - j) % prime
        scale = value * pow(denominator, -1, prime) % prime
        for degree, coefficient in enumerate(basis):
            answer[degree] = (answer[degree] + scale * coefficient) % prime
    return sp.Poly(sum(coefficient * variable**degree for degree, coefficient in enumerate(answer)),
                   variable, modulus=prime)


def c_coerce(value, prime):
    if isinstance(value, tuple):
        return tuple(entry % prime for entry in value)
    return (int(value) % prime, 0, 0)


def c_add(left, right, prime):
    a, b = c_coerce(left, prime), c_coerce(right, prime)
    return tuple((a[i] + b[i]) % prime for i in range(3))


def c_neg(value, prime):
    return tuple(-entry % prime for entry in c_coerce(value, prime))


def c_sub(left, right, prime):
    return c_add(left, c_neg(right, prime), prime)


def c_mul(left, right, relation, prime):
    a, b = c_coerce(left, prime), c_coerce(right, prime)
    raw = [0] * 5
    for i in range(3):
        for j in range(3):
            raw[i + j] = (raw[i + j] + a[i] * b[j]) % prime
    # X^3 = relation[0] + relation[1] X + relation[2] X^2.
    for degree in (4, 3):
        coefficient = raw[degree]
        raw[degree] = 0
        for j in range(3):
            raw[degree - 3 + j] = (
                raw[degree - 3 + j] + coefficient * relation[j]
            ) % prime
    return tuple(raw[:3])


def c_pow(value, exponent, relation, prime):
    require(exponent >= 0, exponent)
    answer = c_coerce(1, prime)
    base = c_coerce(value, prime)
    while exponent:
        if exponent & 1:
            answer = c_mul(answer, base, relation, prime)
        base = c_mul(base, base, relation, prime)
        exponent >>= 1
    return answer


def c_matrix(value, relation, prime):
    basis = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    columns = [c_mul(value, item, relation, prime) for item in basis]
    return [[columns[column][row] for column in range(3)] for row in range(3)]


def det3(matrix, prime):
    a = matrix
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    ) % prime


def c_inverse(value, relation, prime):
    matrix = [row + [1 if index == 0 else 0] for index, row in enumerate(c_matrix(value, relation, prime))]
    for column in range(3):
        pivot = next((row for row in range(column, 3) if matrix[row][column] % prime), None)
        require(pivot is not None, ("nonunit cubic element", value, relation, prime))
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        scale = pow(matrix[column][column] % prime, -1, prime)
        matrix[column] = [(entry * scale) % prime for entry in matrix[column]]
        for row in range(3):
            if row == column:
                continue
            scale = matrix[row][column] % prime
            matrix[row] = [
                (matrix[row][j] - scale * matrix[column][j]) % prime for j in range(4)
            ]
    answer = tuple(matrix[row][3] for row in range(3))
    require(c_mul(value, answer, relation, prime) == c_coerce(1, prime), "cubic inverse")
    return answer


def c_divide(left, right, relation, prime):
    return c_mul(left, c_inverse(right, relation, prime), relation, prime)


def c_norm(value, relation, prime):
    return det3(c_matrix(value, relation, prime), prime)


def c_l_value(point, relation, prime):
    a, b, c = point
    return c_add(
        c_add(
            c_sub(
                c_add(c_mul(27, c_mul(c_pow(a, 2, relation, prime), c_pow(c, 2, relation, prime), relation, prime), relation, prime),
                      c_mul(16, a, relation, prime), prime),
                c_mul(18, c_mul(c_mul(a, b, relation, prime), c, relation, prime), relation, prime), prime),
            c_mul(c_pow(b, 3, relation, prime), c, relation, prime), prime),
        c_neg(c_pow(b, 2, relation, prime), prime), prime)


def c_f_map(point, relation, prime):
    x, y, z = point
    xy = c_mul(x, y, relation, prime)
    unit = c_add(1, xy, prime)
    first = c_add(
        c_mul(c_pow(unit, 3, relation, prime), z, relation, prime),
        c_mul(c_mul(c_pow(y, 2, relation, prime), unit, relation, prime),
              c_add(4, c_mul(3, xy, relation, prime), prime), relation, prime),
        prime)
    second = c_add(
        y,
        c_add(
            c_mul(3, c_mul(c_mul(x, c_pow(unit, 2, relation, prime), relation, prime), z, relation, prime), relation, prime),
            c_mul(3, c_mul(c_mul(x, c_pow(y, 2, relation, prime), relation, prime),
                           c_add(4, c_mul(3, xy, relation, prime), prime), relation, prime), relation, prime),
            prime),
        prime)
    third = c_sub(
        c_sub(c_mul(2, x, relation, prime),
              c_mul(3, c_mul(c_pow(x, 2, relation, prime), y, relation, prime), relation, prime), prime),
        c_mul(c_pow(x, 3, relation, prime), z, relation, prime), prime)
    return first, second, third


def cubic_inverse_point(parent, prime):
    a, b, c = parent
    lead = l_mod(parent, prime)
    require(lead != 0 and c % prime != 0, ("cubic chart", parent, prime))
    lead_inverse = pow(lead, -1, prime)
    # X^3 + alpha X + beta = 0.
    alpha = (4 - 3 * b * c) * lead_inverse % prime
    beta = (-2 * c) * lead_inverse % prime
    relation = ((-beta) % prime, (-alpha) % prime, 0)
    x = (0, 1, 0)
    denominator = c_add(
        c_add(c_mul(12 * a - b * b, c_pow(x, 2, relation, prime), relation, prime),
              c_mul(b, x, relation, prime), prime),
        2, prime)
    numerator_factor = c_add(c_mul(9 * a * c - b, x, relation, prime), 2, prime)
    y = c_sub(
        b,
        c_divide(c_mul(3 * a, c_mul(x, numerator_factor, relation, prime), relation, prime),
                 denominator, relation, prime),
        prime)
    z_numerator = c_sub(
        c_sub(c_mul(2, x, relation, prime), c, prime),
        c_mul(3, c_mul(c_pow(x, 2, relation, prime), y, relation, prime), relation, prime),
        prime)
    z = c_divide(z_numerator, c_pow(x, 3, relation, prime), relation, prime)
    point = (x, y, z)
    require(
        c_f_map(point, relation, prime) == tuple(c_coerce(value, prime) for value in parent),
        ("cubic inverse reconstruction", parent, prime),
    )
    return relation, point


def c_y_block_value(target, scalar, relation, prime):
    a, b, c = target
    r = c_sub(b, scalar, prime)
    return c_sub(
        c_add(c_mul(27, c_mul(c_pow(a, 2, relation, prime), c, relation, prime), relation, prime),
              c_mul(3, c_mul(b, c_pow(r, 2, relation, prime), relation, prime), relation, prime), prime),
        c_add(c_mul(18, c_mul(a, r, relation, prime), relation, prime),
              c_mul(2, c_pow(r, 3, relation, prime), relation, prime), prime),
        prime)


def c_z_block_value(target, scalar, relation, prime):
    a, b, c = target
    L = c_l_value(target, relation, prime)
    a2, a3, a4 = (c_pow(a, exponent, relation, prime) for exponent in (2, 3, 4))
    b2, b3, b4, b5, b6 = (c_pow(b, exponent, relation, prime) for exponent in range(2, 7))
    c2, c3, c4 = (c_pow(c, exponent, relation, prime) for exponent in (2, 3, 4))
    def prod(*values):
        answer = c_coerce(1, prime)
        for value in values:
            answer = c_mul(answer, value, relation, prime)
        return answer
    q2 = c_add(c_sub(c_add(c_mul(324, prod(a2, c2), relation, prime), c_mul(408, a, relation, prime), prime),
                           c_mul(216, prod(a, b, c), relation, prime), prime),
               c_sub(c_mul(6, b2, relation, prime), c_mul(15, prod(b3, c), relation, prime), prime), prime)
    Sz = c_add(c_sub(c_add(c_mul(27, prod(a2, c2), relation, prime), c_mul(52, a, relation, prime), prime),
                           c_mul(18, prod(a, b, c), relation, prime), prime),
               c_add(prod(b3, c), c_mul(14, b2, relation, prime), prime), prime)
    terms = (
        (729, prod(a4, c4)), (-972, prod(a3, b, c3)), (2322, prod(a3, c2)),
        (54, prod(a2, b3, c3)), (270, prod(a2, b2, c2)), (-3735, prod(a2, b, c)),
        (-338, a2), (-36, prod(a, b4, c2)), (122, prod(a, b3, c)), (1372, prod(a, b2)),
        (1, prod(b6, c2)), (-2, prod(b5, c)), (-80, b4),
    )
    Tz = c_coerce(0, prime)
    for coefficient, value in terms:
        Tz = c_add(Tz, c_mul(coefficient, value, relation, prime), prime)
    w = c_coerce(scalar, prime)
    return c_add(
        c_add(c_mul(8, c_pow(w, 3, relation, prime), relation, prime),
              c_mul(q2, c_pow(w, 2, relation, prime), relation, prime), prime),
        c_add(c_mul(6, c_mul(L, c_mul(Sz, w, relation, prime), relation, prime), relation, prime),
              c_mul(L, Tz, relation, prime), prime), prime)


def c_u_block_value(target, scalar, relation, prime):
    _a, b, c = target
    u = c_coerce(scalar, prime)
    return c_add(
        c_l_value(target, relation, prime),
        c_sub(c_mul(4, c_pow(u, 2, relation, prime), relation, prime),
              c_add(c_mul(3, c_mul(b, c_mul(c, c_pow(u, 2, relation, prime), relation, prime), relation, prime), relation, prime),
                    c_mul(2, c_mul(c, c_pow(u, 3, relation, prime), relation, prime), relation, prime), prime), prime),
        prime)


def quotient_norm_poly(parent, name, variable, prime):
    evaluators = {
        "y": c_y_block_value,
        "z": c_z_block_value,
        "u": c_u_block_value,
    }
    relation, inverse_point = cubic_inverse_point(parent, prime)
    values = [
        c_norm(evaluators[name](inverse_point, value, relation, prime), relation, prime)
        for value in range(10)
    ]
    polynomial = interpolate_mod(values, variable, prime).monic()
    require(polynomial.degree() == 9, ("quotient norm degree", parent, name, polynomial))
    return polynomial


def split_quotient_crosscheck(blocks, prime):
    Y, Z, U = sp.symbols("Y Z U")
    grouped = {}
    for block in blocks:
        grouped.setdefault(f_mod(block, prime), []).append(block)
    if any(len(fibre) != 3 for fibre in grouped.values()):
        return False
    for parent, fibre in grouped.items():
        for name, variable, constructor in (
            ("y", Y, y_block),
            ("z", Z, z_block),
            ("u", U, reciprocal_x_block),
        ):
            if parent[2] == 0:
                return False
            direct = sp.Poly(1, variable, modulus=prime)
            for block in fibre:
                direct *= poly_mod(constructor(block, variable), variable, prime)
            try:
                quotient = quotient_norm_poly(parent, name, variable, prime)
            except RuntimeError:
                return False
            if direct.monic() != quotient:
                return False
    return True


def coordinate_collision_row(blocks, prime):
    Y, Z, U = sp.symbols("Y Z U")
    result = {}
    for name, variable, constructor in (
        ("y", Y, y_block),
        ("z", Z, z_block),
        ("u", U, reciprocal_x_block),
    ):
        total = sp.Poly(1, variable, modulus=prime)
        block_polynomials = []
        for block in blocks:
            polynomial = poly_mod(constructor(block, variable), variable, prime)
            require(polynomial.degree() == 3, (prime, name, block, polynomial))
            block_polynomials.append(polynomial)
            total *= polynomial
        total = total.monic()
        gcd = sp.gcd(total, total.diff()).monic()
        pairwise = True
        for i in range(len(block_polynomials)):
            for j in range(i):
                if sp.gcd(block_polynomials[i], block_polynomials[j]).degree() != 0:
                    pairwise = False
                    break
            if not pairwise:
                break
        result[name] = {
            "degree": total.degree(),
            "derivative_gcd_degree": gcd.degree(),
            "derivative_gcd": str(gcd.as_expr()),
            "pairwise_coprime": pairwise,
            "coefficient_digest": digest([int(c) % prime for c in total.all_coeffs()]),
        }
    return result


def quotient_depth_three_collision_row(parents, q0, prime):
    """Test level four through nine exact cubic quotient norms.

    Here ``parents`` is the completely split depth-two fibre above F^3(q0).
    Each of its nine points has a cubic inverse algebra.  Taking a resultant
    over that cubic gives the product of the three last-step coordinate
    blocks without requiring those 27 predecessor points to be rational.
    """
    X, Y, Z, U = sp.symbols("X Y Z U")
    distinguished_parent = f_mod(q0, prime)
    if distinguished_parent not in parents:
        return None

    norm_rows = {"y": [], "z": [], "u": []}
    l_root_rows = []
    for parent in parents:
        if l_mod(parent, prime) == 0 or parent[2] == 0:
            return None
        outer = poly_mod(
            l_mod(parent, prime) * X**3
            + (4 - 3 * parent[1] * parent[2]) * X
            - 2 * parent[2],
            X,
            prime,
        )
        if outer.degree() != 3 or sp.gcd(outer, outer.diff()).degree() != 0:
            return None
        try:
            relation, qx = cubic_inverse_point(parent, prime)
        except RuntimeError:
            return None

        l_element = c_l_value(qx, relation, prime)
        l_gcd = sp.gcd(
            outer,
            sp.Poly(sum(l_element[degree] * X**degree for degree in range(3)), X, modulus=prime),
        ).monic()
        l_root_rows.append((parent, l_gcd))

        for name, variable in (
            ("y", Y),
            ("z", Z),
            ("u", U),
        ):
            try:
                norm_poly = quotient_norm_poly(parent, name, variable, prime)
            except RuntimeError:
                return None
            norm_rows[name].append(norm_poly)

    nontrivial_l_rows = [(parent, value) for parent, value in l_root_rows if value.degree() > 0]
    if len(nontrivial_l_rows) != 1:
        return None
    parent, l_gcd = nontrivial_l_rows[0]
    expected_l_root = sp.Poly(X - q0[0], X, modulus=prime).monic()
    if parent != distinguished_parent or l_gcd != expected_l_root:
        return None

    answer = {}
    for name, variable, constructor in (
        ("y", Y, y_block),
        ("z", Z, z_block),
        ("u", U, reciprocal_x_block),
    ):
        total = sp.Poly(1, variable, modulus=prime)
        for norm_poly in norm_rows[name]:
            total *= norm_poly
        total = total.monic()
        gcd = sp.gcd(total, total.diff()).monic()
        q0_poly = poly_mod(constructor(q0, variable), variable, prime)
        expected_gcd = sp.gcd(q0_poly, q0_poly.diff()).monic()
        if total.degree() != 81 or expected_gcd.degree() != 1 or gcd != expected_gcd:
            return None
        answer[name] = {
            "degree": total.degree(),
            "derivative_gcd_degree": gcd.degree(),
            "derivative_gcd": str(gcd.as_expr()),
            "nine_parent_norm_degrees": [row.degree() for row in norm_rows[name]],
            "coefficient_digest": digest([int(c) % prime for c in total.all_coeffs()]),
        }
    return answer


def split_tree(reverse, target, depth):
    level = [target]
    for _ in range(depth):
        next_level = []
        for point in level:
            fibre = reverse.get(point, ())
            if len(fibre) != 3:
                return None
            next_level.extend(fibre)
        if len(set(next_level)) != len(next_level):
            return None
        level = next_level
    return tuple(sorted(level))


def finite_field_witnesses(prime, depths):
    reverse = {}
    l_points = []
    for x in range(prime):
        for y in range(prime):
            for z in range(prime):
                point = (x, y, z)
                reverse.setdefault(f_mod(point, prime), []).append(point)
                if l_mod(point, prime) == 0:
                    l_points.append(point)

    pending = set(depths)
    witnesses = {}
    for q0 in l_points:
        if not pending:
            break
        # Generic one-step ramified-block gates for y and reciprocal x.
        if (
            q0[0] == 0
            or q0[2] == 0
            or (4 - 3 * q0[1] * q0[2]) % prime == 0
            or s_mod(q0, prime) == 0
            or z_square_factor_mod(q0, prime) == 0
            or not l_smooth_mod(q0, prime)
        ):
            continue
        for depth in sorted(tuple(pending)):
            eta = q0
            for _ in range(depth):
                eta = f_mod(eta, prime)
            if depth == 3:
                parents = split_tree(reverse, eta, 2)
                if parents is None:
                    continue
                collision = quotient_depth_three_collision_row(parents, q0, prime)
                if collision is None:
                    continue
                witnesses[depth] = {
                    "level": 4,
                    "q0": q0,
                    "eta": eta,
                    "split_parent_blocks": len(parents),
                    "algebraic_preceding_blocks": 27,
                    "unique_L_block": True,
                    "method": "nine split parents plus exact cubic quotient norms",
                    "collision": collision,
                }
                pending.remove(depth)
                continue
            blocks = split_tree(reverse, eta, depth)
            if blocks is None or q0 not in blocks:
                continue
            l_blocks = [point for point in blocks if l_mod(point, prime) == 0]
            if l_blocks != [q0]:
                continue
            if any(point[2] == 0 for point in blocks):
                continue
            collision = coordinate_collision_row(blocks, prime)
            quotient_crosscheck = split_quotient_crosscheck(blocks, prime)
            if not quotient_crosscheck:
                continue
            if not all(
                row["degree"] == 3 ** (depth + 1)
                and row["derivative_gcd_degree"] == 1
                and row["pairwise_coprime"]
                for row in collision.values()
            ):
                continue
            witnesses[depth] = {
                "level": depth + 1,
                "q0": q0,
                "eta": eta,
                "preceding_blocks": len(blocks),
                "unique_L_block": True,
                "split_quotient_crosscheck": True,
                "collision": collision,
            }
            pending.remove(depth)
    return {
        "prime": prime,
        "L_point_count": len(l_points),
        "witnesses": witnesses,
        "missing_depths": sorted(pending),
    }


def raw_reciprocal_ledger():
    rows = []
    for level in range(1, 9):
        degree = 3**level
        monic_x_baseline = 3 - 2 * degree
        scalar_clearing_shift = 2 * degree - 2
        raw_baseline = monic_x_baseline + scalar_clearing_shift
        require(raw_baseline == 1, (level, degree, raw_baseline))
        for index_length in (0, 1, 2, 7):
            reciprocal_exponent = 1 + 2 * index_length
            monic_x_exponent = monic_x_baseline + 2 * index_length
            raw_x_exponent = monic_x_exponent + scalar_clearing_shift
            require(raw_x_exponent == reciprocal_exponent, (level, index_length))
        rows.append(
            {
                "level": level,
                "degree": degree,
                "v_norm_x": -1,
                "monic_x_discriminant_baseline": monic_x_baseline,
                "raw_scalar_shift": scalar_clearing_shift,
                "reciprocal_and_raw_baseline": raw_baseline,
            }
        )
    return rows


def main() -> None:
    rational = rational_level_two_witness()
    finite_rows = []
    missing = {1, 2, 3}
    for prime in (31, 41, 43, 47, 53, 59, 61):
        if not missing:
            break
        row = finite_field_witnesses(prime, sorted(missing))
        finite_rows.append(row)
        missing -= set(row["witnesses"])

    record = {
        "verdict": "exact cross-block noncollision witnesses; all-level criterion remains separate",
        "rational_level_two": rational,
        "finite_field_rows": finite_rows,
        "reciprocal_raw_clearing_ledger": raw_reciprocal_ledger(),
        "unwitnessed_depths": sorted(missing),
    }
    record["semantic_sha256"] = digest(record)
    require(record["semantic_sha256"] == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", record["semantic_sha256"]))
    print(json.dumps(record, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
