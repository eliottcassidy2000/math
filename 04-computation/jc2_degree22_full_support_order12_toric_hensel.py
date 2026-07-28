#!/usr/bin/env python3
"""Exact full-support BCDEW order-12 Hensel/toric obstruction.

This is a downstream companion to THM-2636 and THM-2683.  In the remaining
degree-22 coefficient chart put

    B=t^2, C=c t^3, D=d t^4, E=e t^5, W=w t^6.

The order-11 Hensel equations have the nine-monomial support

    E, EW, DE, C, CW, CD, CD^2, C^2E, C^3.

On C != 0 set a=E/C and h=C^2 and divide the equations by C.  Their monomial
vector becomes

    (a, aW, aD, 1, W, D, D^2, ah, h).

Four independent terminal equations are linear in W and h.  Two of them
solve W,h; the other two give quadratic compatibility polynomials F(a,D)
and G(a,D).  Their D-resultant contains the square of the chosen 2 by 2
pivot determinant.  Exact saturation leaves a squarefree degree-seven
polynomial H7(a).  A linear subresultant recovers D in K[a]/(H7), after which
W and h are recovered as well.  The five order-12 Hensel rows have no common
zero in this seven-dimensional algebra.  The root field is checked directly.
For the larger pair field, the carrier and its equality across all six pivot
charts are constructed exactly first; a good simple degree-one residue then
certifies the quotient-algebra unit identity.  Degree preservation,
squarefreeness, every denominator, and every torus unit are audited, so this
is a characteristic-zero good-reduction proof rather than a numerical probe.

The degree seven is intrinsic.  Elimination gives the eight-binomial toric
ideal, and the faithful exponent lattice (a,D,W,h) has a regular unimodular
triangulation with seven maximal simplices.  We audit both the degree-five
root field and the degree-ten unordered-root-pair field.
"""

from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import io
import math
from itertools import combinations
from pathlib import Path

from gmpy2 import mpq as Fraction
import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "jc2_degree22_bcd_weighted_hensel_kummer_thm2636.py"
BASE_SHA256 = "0866a29f665aedc6d2c226f35943852e56907ff821e705a0dbca2651e71fa15c"
require(
    hashlib.sha256(BASE_PATH.read_bytes()).hexdigest() == BASE_SHA256,
    "audited THM-2636 dependency changed",
)

spec = importlib.util.spec_from_file_location("thm2636_full_support_base", BASE_PATH)
require(spec is not None and spec.loader is not None, "cannot load THM-2636")
base = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(base)

t, v = base.t, base.v
c, d, e, w = base.c, base.d, base.e, base.w


class GenericHensel(base.HenselSystem):
    """THM-2636's exact Hensel engine with arbitrary exponent-tuple length."""

    def __init__(self, name, field, q0, s0, parameters):
        self.parameters = tuple(parameters)
        self.zero_exponent = (0,) * len(self.parameters)
        super().__init__(name, field, q0, s0)

    def pc(self, terms=None):
        result = {}
        for exponent, coefficient in (terms or {}).items():
            require(
                len(exponent) == len(self.parameters),
                "parameter exponent has the wrong length",
            )
            reduced = self.field.reduce(coefficient)
            if not self.field.zero(reduced):
                result[tuple(exponent)] = reduced
        return result

    def pc_const(self, coefficient):
        return self.pc({self.zero_exponent: s.sympify(coefficient)})

    def pc_mul(self, left, right):
        result = {}
        for exponent1, coefficient1 in left.items():
            for exponent2, coefficient2 in right.items():
                exponent = tuple(a + b for a, b in zip(exponent1, exponent2))
                result[exponent] = self.field.reduce(
                    result.get(exponent, 0) + coefficient1 * coefficient2
                )
        return self.pc(result)

    def v_divmod_monic(self, dividend, divisor):
        divisor = self.vtrim(divisor)
        require(
            divisor[-1] == {self.zero_exponent: s.Integer(1)},
            "Hensel divisor is not monic",
        )
        remainder = self.vtrim(dividend[:])
        quotient = [
            {} for _ in range(max(1, len(remainder) - len(divisor) + 1))
        ]
        while len(remainder) >= len(divisor) and not (
            len(remainder) == 1 and not remainder[0]
        ):
            shift = len(remainder) - len(divisor)
            leading = remainder[-1]
            quotient[shift] = self.pc_add(quotient[shift], leading)
            subtraction = [{} for _ in range(shift)] + [
                self.pc_mul(coefficient, leading) for coefficient in divisor
            ]
            remainder = self.v_add(remainder, self.v_neg(subtraction))
        return self.vtrim(quotient), self.vtrim(remainder)

    def inverse_s0_mod_q0(self, q0_field, s0_field):
        q0 = self.v_from_field(q0_field)
        s0 = self.v_from_field(s0_field)
        columns = []
        for exponent in range(self.k):
            monomial = [{} for _ in range(exponent)] + [self.pc_const(1)]
            remainder = self.v_remainder(self.v_mul(s0, monomial), q0)
            entries = []
            for index in range(self.k):
                coefficient = remainder[index] if index < len(remainder) else {}
                require(
                    set(coefficient).issubset({self.zero_exponent}),
                    "fixed Hensel inverse acquired parameters",
                )
                entries.append(
                    coefficient.get(self.zero_exponent, s.Integer(0))
                )
            columns.append(entries)
        if self.k == 1:
            inverse_coefficients = [self.field.inverse(columns[0][0])]
        else:
            a00, a10 = columns[0]
            a01, a11 = columns[1]
            determinant = self.field.reduce(a00 * a11 - a01 * a10)
            require(not self.field.zero(determinant), "Hensel inverse is singular")
            inverse_coefficients = [
                self.field.reduce(a11 * self.field.inverse(determinant)),
                self.field.reduce(-a10 * self.field.inverse(determinant)),
            ]
        inverse = self.v_from_field(inverse_coefficients)
        require(
            self.v_equal(
                self.v_remainder(self.v_mul(s0, inverse), q0),
                self.v_from_field([1]),
            ),
            "fixed Hensel inverse failed",
        )
        return inverse

    def expression_to_vpoly(self, expression):
        poly = s.Poly(expression, v, *self.parameters, domain=s.QQ)
        result = [{} for _ in range(max(0, poly.degree(v)) + 1)]
        for monomial, coefficient in poly.terms():
            iv, *parameter_exponents = monomial
            result[iv] = self.pc_add(
                result[iv],
                self.pc({tuple(parameter_exponents): coefficient}),
            )
        return self.vtrim(result)

    def lift(self, coefficients, through):
        fixed = self.expression_to_vpoly(coefficients[0])
        require(
            self.v_equal(self.v_mul(self.q0, self.s0), fixed),
            f"{self.name}: fixed factorization failed",
        )
        qs, ss = [self.q0], [self.s0]
        zero = self.v_from_field([0])
        for n in range(1, through + 1):
            rn = (
                self.expression_to_vpoly(coefficients[n])
                if n < len(coefficients)
                else zero
            )
            convolution = zero
            for index in range(1, n):
                convolution = self.v_add(
                    convolution, self.v_mul(qs[index], ss[n - index])
                )
            residual = self.v_add(rn, self.v_neg(convolution))
            qn = self.v_remainder(
                self.v_mul(residual, self.inv_s0_mod_q0), self.q0
            )
            numerator = self.v_add(
                residual, self.v_neg(self.v_mul(qn, self.s0))
            )
            sn, remainder = self.v_divmod_monic(numerator, self.q0)
            require(
                all(not coefficient for coefficient in remainder),
                f"{self.name}: Hensel division failed at order {n}",
            )
            qs.append(qn)
            ss.append(sn)
            reconstructed = zero
            for index in range(n + 1):
                reconstructed = self.v_add(
                    reconstructed, self.v_mul(qs[index], ss[n - index])
                )
            require(
                self.v_equal(reconstructed, rn),
                f"{self.name}: product control failed at order {n}",
            )
        return qs, ss

    def matrix_rref(self, matrix):
        reduced = [[self.field.reduce(value) for value in row] for row in matrix]
        rows = len(reduced)
        columns = len(reduced[0])
        pivots = []
        pivot_row = 0
        for column in range(columns):
            pivot = next(
                (
                    row
                    for row in range(pivot_row, rows)
                    if not self.field.zero(reduced[row][column])
                ),
                None,
            )
            if pivot is None:
                continue
            reduced[pivot_row], reduced[pivot] = (
                reduced[pivot],
                reduced[pivot_row],
            )
            inverse = self.field.inverse(reduced[pivot_row][column])
            reduced[pivot_row] = [
                self.field.reduce(value * inverse)
                for value in reduced[pivot_row]
            ]
            for row in range(rows):
                if row == pivot_row or self.field.zero(reduced[row][column]):
                    continue
                multiplier = reduced[row][column]
                reduced[row] = [
                    self.field.reduce(
                        reduced[row][index]
                        - multiplier * reduced[pivot_row][index]
                    )
                    for index in range(columns)
                ]
            pivots.append(column)
            pivot_row += 1
            if pivot_row == rows:
                break
        return tuple(pivots), reduced


class FastField:
    """Power-basis arithmetic in the audited number field using Fractions."""

    def __init__(self, slow_field):
        self.generator = slow_field.generator
        self.degree = slow_field.modulus.degree()
        self.modulus_tail = tuple(
            self._fraction(slow_field.modulus.nth(index))
            for index in range(self.degree)
        )
        self.zero = (Fraction(0),) * self.degree
        self.one = (Fraction(1),) + (Fraction(0),) * (self.degree - 1)
        self.inverse_cache = {}

    @staticmethod
    def _fraction(value):
        value = s.Rational(value)
        return Fraction(int(value.p), int(value.q))

    def from_sympy(self, value):
        poly = s.Poly(s.expand(value), self.generator, domain=s.QQ)
        require(poly.degree() < self.degree, "field element was not reduced")
        return tuple(
            self._fraction(poly.nth(index)) for index in range(self.degree)
        )

    def iszero(self, value):
        return value == self.zero

    def add(self, left, right):
        return tuple(a + b for a, b in zip(left, right))

    def neg(self, value):
        return tuple(-a for a in value)

    def sub(self, left, right):
        return tuple(a - b for a, b in zip(left, right))

    def scale_int(self, value, scalar):
        return tuple(Fraction(scalar) * entry for entry in value)

    def mul(self, left, right):
        work = [Fraction(0)] * (2 * self.degree - 1)
        for i, left_coefficient in enumerate(left):
            if left_coefficient:
                for j, right_coefficient in enumerate(right):
                    work[i + j] += left_coefficient * right_coefficient
        for exponent in range(2 * self.degree - 2, self.degree - 1, -1):
            coefficient = work[exponent]
            if coefficient:
                shift = exponent - self.degree
                for index, modulus_coefficient in enumerate(self.modulus_tail):
                    work[shift + index] -= coefficient * modulus_coefficient
        return tuple(work[: self.degree])

    def inv(self, value):
        require(not self.iszero(value), "attempted inversion of zero")
        if value in self.inverse_cache:
            return self.inverse_cache[value]
        columns = []
        for column in range(self.degree):
            basis = tuple(
                Fraction(1) if index == column else Fraction(0)
                for index in range(self.degree)
            )
            columns.append(self.mul(value, basis))
        matrix = [
            [columns[column][row] for column in range(self.degree)]
            + [Fraction(1) if row == 0 else Fraction(0)]
            for row in range(self.degree)
        ]
        for column in range(self.degree):
            pivot = next(
                row
                for row in range(column, self.degree)
                if matrix[row][column]
            )
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            pivot_value = matrix[column][column]
            matrix[column] = [value / pivot_value for value in matrix[column]]
            for row in range(self.degree):
                if row == column or not matrix[row][column]:
                    continue
                multiplier = matrix[row][column]
                matrix[row] = [
                    value - multiplier * pivot_entry
                    for value, pivot_entry in zip(matrix[row], matrix[column])
                ]
        result = tuple(matrix[index][-1] for index in range(self.degree))
        require(self.mul(value, result) == self.one, "field inverse failed")
        self.inverse_cache[value] = result
        return result


class PrimeField:
    """The tiny interface used by the sparse polynomial code over F_p."""

    def __init__(self, prime):
        self.prime = prime
        self.zero = 0
        self.one = 1

    def iszero(self, value):
        return value % self.prime == 0

    def add(self, left, right):
        return (left + right) % self.prime

    def neg(self, value):
        return (-value) % self.prime

    def sub(self, left, right):
        return (left - right) % self.prime

    def scale_int(self, value, scalar):
        return scalar * value % self.prime

    def mul(self, left, right):
        return (left * right) % self.prime

    def inv(self, value):
        require(value % self.prime, "prime-field inverse of zero")
        return pow(int(value), -1, self.prime)


def exact_toric_certificate():
    """Eliminate the torus parameters and certify the degree-seven lattice."""

    A, D, W, H, lam = s.symbols("A D W H lam")
    xs = s.symbols("x0:9")
    monomials = (
        A,
        A * W,
        A * D,
        s.Integer(1),
        W,
        D,
        D**2,
        A * H,
        H,
    )
    graph = [xs[index] - lam * monomials[index] for index in range(9)]
    groebner = s.groebner(graph, A, D, W, H, lam, *xs, order="lex")
    eliminated = tuple(
        s.expand(poly.as_expr())
        for poly in groebner.polys
        if not (poly.as_expr().free_symbols & {A, D, W, H, lam})
    )
    expected = {
        s.expand(xs[0] * xs[4] - xs[1] * xs[3]),
        s.expand(xs[0] * xs[5] - xs[2] * xs[3]),
        s.expand(xs[0] * xs[6] - xs[2] * xs[5]),
        s.expand(xs[0] * xs[8] - xs[3] * xs[7]),
        s.expand(xs[1] * xs[5] - xs[2] * xs[4]),
        s.expand(xs[1] * xs[8] - xs[4] * xs[7]),
        s.expand(xs[2] * xs[8] - xs[5] * xs[7]),
        s.expand(xs[3] * xs[6] - xs[5] ** 2),
    }
    require(set(eliminated) == expected, "projective toric ideal changed")

    # These are the faithful exponents after quotienting the index-two
    # (C,E)->(-C,-E) redundancy via a=E/C and h=C^2.
    points = (
        (1, 0, 0, 0),
        (1, 0, 1, 0),
        (1, 1, 0, 0),
        (0, 0, 0, 0),
        (0, 0, 1, 0),
        (0, 1, 0, 0),
        (0, 2, 0, 0),
        (1, 0, 0, 1),
        (0, 0, 0, 1),
    )
    original_points = (
        (0, 0, 1, 0),  # E
        (0, 0, 1, 1),  # EW
        (0, 1, 1, 0),  # DE
        (1, 0, 0, 0),  # C
        (1, 0, 0, 1),  # CW
        (1, 1, 0, 0),  # CD
        (1, 2, 0, 0),  # CD^2
        (2, 0, 1, 0),  # C^2E
        (3, 0, 0, 0),  # C^3
    )
    weights = (-1, 1, 1, 0, 1, -1, 2, 1, 0)
    lower = []
    original_affine_minor_gcd = 0
    for indices in combinations(range(9), 5):
        matrix = s.Matrix([[*points[index], 1] for index in indices])
        determinant = matrix.det()
        original_determinant = s.Matrix(
            [[*original_points[index], 1] for index in indices]
        ).det()
        if original_determinant:
            original_affine_minor_gcd = math.gcd(
                original_affine_minor_gcd, abs(int(original_determinant))
            )
        if determinant == 0:
            continue
        affine = matrix.LUsolve(s.Matrix([weights[index] for index in indices]))
        gaps = []
        for index, point in enumerate(points):
            predicted = sum(affine[j] * point[j] for j in range(4)) + affine[4]
            gaps.append(s.Rational(weights[index]) - predicted)
        if all(gap >= 0 for gap in gaps):
            require(
                all(gaps[index] > 0 for index in range(9) if index not in indices),
                "lower triangulation acquired a nonsimplicial cell",
            )
            lower.append((indices, abs(determinant), abs(original_determinant)))
    expected_simplices = (
        (0, 1, 2, 5, 7),
        (0, 1, 4, 5, 8),
        (0, 1, 5, 7, 8),
        (0, 3, 4, 5, 8),
        (1, 2, 5, 6, 7),
        (1, 4, 5, 6, 8),
        (1, 5, 6, 7, 8),
    )
    require(
        tuple(indices for indices, _, _ in lower) == expected_simplices,
        "regular lower triangulation changed",
    )
    require(
        all(determinant == 1 for _, determinant, _ in lower),
        "lower triangulation is no longer unimodular",
    )
    require(
        all(determinant == 2 for _, _, determinant in lower)
        and original_affine_minor_gcd == 2,
        "original CDEW affine lattice no longer has index two",
    )
    return {
        "ideal_generators": len(eliminated),
        "simplices": expected_simplices,
        "degree": sum(determinant for _, determinant, _ in lower),
        "original_index": original_affine_minor_gcd,
    }


TORIC = exact_toric_certificate()
require(TORIC["degree"] == 7, "terminal toric degree changed")


# Bivariate polynomials in (a,D) are sparse dictionaries (ia,iD)->K.
def bclean(field, poly):
    return {monomial: value for monomial, value in poly.items() if not field.iszero(value)}


def badd(field, left, right):
    result = dict(left)
    for monomial, value in right.items():
        result[monomial] = field.add(result.get(monomial, field.zero), value)
    return bclean(field, result)


def bneg(field, poly):
    return bclean(field, {monomial: field.neg(value) for monomial, value in poly.items()})


def bsub(field, left, right):
    return badd(field, left, bneg(field, right))


def bmul(field, left, right):
    result = {}
    for (ia, iD), left_value in left.items():
        for (ja, jD), right_value in right.items():
            monomial = (ia + ja, iD + jD)
            result[monomial] = field.add(
                result.get(monomial, field.zero),
                field.mul(left_value, right_value),
            )
    return bclean(field, result)


def bdet(field, matrix):
    size = len(matrix)
    if size == 0:
        return {(0, 0): field.one}
    total = {}
    for column in range(size):
        minor = [row[:column] + row[column + 1 :] for row in matrix[1:]]
        term = bmul(field, matrix[0][column], bdet(field, minor))
        total = badd(field, total, term if column % 2 == 0 else bneg(field, term))
    return total


def bdegree_D(poly):
    return max((iD for _, iD in poly), default=-1)


def bcoefficient_D(poly, exponent):
    return {(ia, 0): value for (ia, iD), value in poly.items() if iD == exponent}


def bresultant_D(field, left, right):
    left_degree = bdegree_D(left)
    right_degree = bdegree_D(right)
    require(left_degree >= 0 and right_degree >= 0, "zero resultant input")
    left_coefficients = [
        bcoefficient_D(left, exponent)
        for exponent in range(left_degree, -1, -1)
    ]
    right_coefficients = [
        bcoefficient_D(right, exponent)
        for exponent in range(right_degree, -1, -1)
    ]
    size = left_degree + right_degree
    sylvester = []
    for shift in range(right_degree):
        sylvester.append(
            [{}] * shift
            + left_coefficients
            + [{}] * (size - shift - len(left_coefficients))
        )
    for shift in range(left_degree):
        sylvester.append(
            [{}] * shift
            + right_coefficients
            + [{}] * (size - shift - len(right_coefficients))
        )
    result = bdet(field, sylvester)
    require(all(iD == 0 for _, iD in result), "D-resultant retained D")
    return result


# Univariate polynomials in a are dictionaries ia->K.
def b_to_u(poly):
    require(all(iD == 0 for _, iD in poly), "bivariate polynomial retained D")
    return {ia: value for (ia, _), value in poly.items()}


def u_to_b(poly):
    return {(ia, 0): value for ia, value in poly.items()}


def uclean(field, poly):
    return {degree: value for degree, value in poly.items() if not field.iszero(value)}


def udegree(poly):
    return max(poly, default=-1)


def uadd(field, left, right):
    result = dict(left)
    for degree, value in right.items():
        result[degree] = field.add(result.get(degree, field.zero), value)
    return uclean(field, result)


def uneg(field, poly):
    return uclean(field, {degree: field.neg(value) for degree, value in poly.items()})


def usub(field, left, right):
    return uadd(field, left, uneg(field, right))


def umul(field, left, right):
    result = {}
    for left_degree, left_value in left.items():
        for right_degree, right_value in right.items():
            degree = left_degree + right_degree
            result[degree] = field.add(
                result.get(degree, field.zero),
                field.mul(left_value, right_value),
            )
    return uclean(field, result)


def uscale(field, poly, scalar):
    return uclean(
        field,
        {degree: field.mul(value, scalar) for degree, value in poly.items()},
    )


def udivmod(field, dividend, divisor):
    require(divisor, "univariate division by zero")
    quotient = {}
    remainder = dict(dividend)
    divisor_degree = udegree(divisor)
    inverse_leading = field.inv(divisor[divisor_degree])
    while remainder and udegree(remainder) >= divisor_degree:
        shift = udegree(remainder) - divisor_degree
        scalar = field.mul(remainder[udegree(remainder)], inverse_leading)
        quotient[shift] = field.add(quotient.get(shift, field.zero), scalar)
        subtraction = {
            degree + shift: field.mul(value, scalar)
            for degree, value in divisor.items()
        }
        remainder = usub(field, remainder, subtraction)
    return uclean(field, quotient), uclean(field, remainder)


def umonic(field, poly):
    require(poly, "cannot normalize zero polynomial")
    inverse = field.inv(poly[udegree(poly)])
    return uscale(field, poly, inverse)


def ugcd(field, left, right):
    left, right = uclean(field, left), uclean(field, right)
    while right:
        _, remainder = udivmod(field, left, right)
        left, right = right, remainder
    return umonic(field, left) if left else {}


def uderivative(field, poly):
    return uclean(
        field,
        {
            degree - 1: field.scale_int(value, degree)
            for degree, value in poly.items()
            if degree
        },
    )


def usaturate(field, poly, divisor):
    result = umonic(field, poly)
    removed = 0
    while True:
        common = ugcd(field, result, divisor)
        if udegree(common) <= 0:
            return umonic(field, result), removed
        quotient, remainder = udivmod(field, result, common)
        require(not remainder, "exact saturation division failed")
        removed += udegree(common)
        result = quotient


def urem(field, dividend, modulus):
    return udivmod(field, dividend, modulus)[1]


def uxgcd(field, left, right):
    old_r, current_r = dict(left), dict(right)
    old_s, current_s = {0: field.one}, {}
    old_t, current_t = {}, {0: field.one}
    while current_r:
        quotient, remainder = udivmod(field, old_r, current_r)
        old_r, current_r = current_r, remainder
        old_s, current_s = current_s, usub(
            field, old_s, umul(field, quotient, current_s)
        )
        old_t, current_t = current_t, usub(
            field, old_t, umul(field, quotient, current_t)
        )
    require(old_r, "extended gcd of two zero polynomials")
    inverse = field.inv(old_r[udegree(old_r)])
    return (
        uscale(field, old_r, inverse),
        uscale(field, old_s, inverse),
        uscale(field, old_t, inverse),
    )


def uinverse_mod(field, value, modulus):
    gcd, coefficient, _ = uxgcd(field, value, modulus)
    require(udegree(gcd) == 0, "nonunit in terminal quotient algebra")
    return urem(field, coefficient, modulus)


def udigest(poly):
    payload = repr(tuple(sorted(poly.items()))).encode()
    return hashlib.sha256(payload).hexdigest()


def terminal_field_certificate(
    name, slow_field, q0, s0, coefficients, quotient_mode
):
    require(quotient_mode in {"exact", "good_residue"}, "unknown quotient mode")
    system = GenericHensel(name, slow_field, q0, s0, (c, d, e, w))
    qs, ss = system.lift(coefficients, 12)
    equations11 = qs[11] + ss[11]
    equations12 = qs[12] + ss[12]
    require(len(equations11) == len(equations12) == 5, "Hensel row count changed")

    support11 = tuple(sorted(set().union(*(set(row) for row in equations11))))
    support12 = tuple(sorted(set().union(*(set(row) for row in equations12))))
    expected11 = (
        (0, 0, 1, 0),
        (0, 0, 1, 1),
        (0, 1, 1, 0),
        (1, 0, 0, 0),
        (1, 0, 0, 1),
        (1, 1, 0, 0),
        (1, 2, 0, 0),
        (2, 0, 1, 0),
        (3, 0, 0, 0),
    )
    expected12 = (
        (0, 0, 0, 0),
        (0, 0, 0, 1),
        (0, 0, 0, 2),
        (0, 0, 2, 0),
        (0, 1, 0, 0),
        (0, 1, 0, 1),
        (0, 2, 0, 0),
        (0, 3, 0, 0),
        (1, 0, 1, 0),
        (1, 1, 1, 0),
        (2, 0, 0, 0),
        (2, 0, 0, 1),
        (2, 1, 0, 0),
        (4, 0, 0, 0),
    )
    require(support11 == expected11, f"{name}: order-11 support changed")
    require(support12 == expected12, f"{name}: order-12 support changed")

    matrix11 = [
        [row.get(monomial, s.Integer(0)) for monomial in support11]
        for row in equations11
    ]
    pivots11, rref11 = system.matrix_rref(matrix11)
    require(pivots11 == (0, 1, 2, 3), f"{name}: terminal rank changed")
    require(
        all(slow_field.zero(value) for value in rref11[4]),
        f"{name}: terminal rank exceeds four",
    )
    selected_rows = next(
        rows
        for rows in combinations(range(5), 4)
        if len(system.matrix_rref([matrix11[index] for index in rows])[0]) == 4
    )
    selected_raw = [matrix11[index] for index in selected_rows]

    matrix12 = [
        [row.get(monomial, s.Integer(0)) for monomial in support12]
        for row in equations12
    ]
    pivots12, rref12 = system.matrix_rref(matrix12)
    require(pivots12 == (0, 1, 2, 3), f"{name}: order-12 rank changed")
    require(
        all(slow_field.zero(value) for value in rref12[4]),
        f"{name}: order-12 rank exceeds four",
    )

    field = FastField(slow_field)
    selected = [
        [field.from_sympy(slow_field.reduce(value)) for value in row]
        for row in selected_raw
    ]

    def linear(constant, coefficient):
        return bclean(field, {(0, 0): constant, (1, 0): coefficient})

    components = []
    for row in selected:
        coefficient_W = linear(row[4], row[1])
        coefficient_h = linear(row[8], row[7])
        remainder = bclean(
            field,
            {
                (0, 0): row[3],
                (1, 0): row[0],
                (0, 1): row[5],
                (1, 1): row[2],
                (0, 2): row[6],
            },
        )
        components.append((coefficient_W, coefficient_h, remainder))

    deltas = {}
    for left, right in combinations(range(4), 2):
        deltas[(left, right)] = bsub(
            field,
            bmul(field, components[left][0], components[right][1]),
            bmul(field, components[right][0], components[left][1]),
        )
    delta_gcd = None
    for delta_poly in deltas.values():
        delta_u = b_to_u(delta_poly)
        delta_gcd = (
            umonic(field, delta_u)
            if delta_gcd is None
            else ugcd(field, delta_gcd, delta_u)
        )
    require(
        delta_gcd is not None and udegree(delta_gcd) == 0,
        f"{name}: all W,h pivot minors vanish at a common a",
    )

    pivot_results = []
    exact_carriers = []
    residue_payload = None
    for pair, delta_bivariate in deltas.items():
        remaining = tuple(index for index in range(4) if index not in pair)
        compatibilities = [
            bdet(
                field,
                [
                    list(components[pair[0]]),
                    list(components[pair[1]]),
                    list(components[index]),
                ],
            )
            for index in remaining
        ]
        F, G = compatibilities
        require(
            bdegree_D(F) == bdegree_D(G) == 2,
            f"{name}/{pair}: compatibility degree changed",
        )
        raw_resultant = b_to_u(bresultant_D(field, F, G))
        require(raw_resultant, f"{name}/{pair}: pivot chart is not finite")
        delta = b_to_u(delta_bivariate)
        H7, removed_degree = usaturate(field, raw_resultant, delta)
        artifact, artifact_remainder = udivmod(field, raw_resultant, H7)
        require(
            udegree(raw_resultant) == 11,
            f"{name}/{pair}: raw resultant degree changed",
        )
        require(
            udegree(delta) == 2
            and removed_degree == 4
            and not artifact_remainder
            and umonic(field, artifact)
            == umonic(field, umul(field, delta, delta)),
            f"{name}/{pair}: pivot artifact is not exactly delta squared",
        )
        require(
            udegree(H7) == TORIC["degree"],
            f"{name}/{pair}: terminal degree is not seven",
        )
        require(
            udegree(ugcd(field, H7, delta)) == 0,
            f"{name}/{pair}: pivot vanishes on its saturated carrier",
        )
        exact_carriers.append(H7)
        if quotient_mode == "good_residue":
            if residue_payload is None:
                residue_payload = (pair, H7, delta, F, G)
            pivot_results.append(
                {
                    "pair": pair,
                    "raw_resultant_degree": udegree(raw_resultant),
                    "removed_degree": removed_degree,
                    "H7_degree": udegree(H7),
                    "H7_digest": udigest(H7),
                    "subresultant_U_degree": None,
                    "order12_gcd_degrees": (),
                    "order12_residue_degrees": (),
                    "order12_residue_digests": (),
                }
            )
            continue

        require(
            udegree(ugcd(field, H7, uderivative(field, H7))) == 0,
            f"{name}/{pair}: terminal degree-seven algebra is nonreduced",
        )

        leading_F = b_to_u(bcoefficient_D(F, 2))
        leading_G = b_to_u(bcoefficient_D(G, 2))
        infinity = ugcd(field, H7, ugcd(field, leading_F, leading_G))
        require(
            udegree(infinity) == 0,
            f"{name}/{pair}: terminal point escaped to D=infinity",
        )

        # Cancel D^2 to obtain U(a)D+V(a) on the common zero scheme.
        linear_subresultant = bsub(
            field,
            bmul(field, u_to_b(leading_G), F),
            bmul(field, u_to_b(leading_F), G),
        )
        require(
            bdegree_D(linear_subresultant) == 1,
            f"{name}/{pair}: linear subresultant degree changed",
        )
        U = b_to_u(bcoefficient_D(linear_subresultant, 1))
        V = b_to_u(bcoefficient_D(linear_subresultant, 0))
        require(
            udegree(ugcd(field, H7, U)) == 0,
            f"{name}/{pair}: D reconstruction denominator vanished",
        )

        def Aadd(left, right):
            return urem(field, uadd(field, left, right), H7)

        def Aneg(value):
            return urem(field, uneg(field, value), H7)

        def Amul(left, right):
            return urem(field, umul(field, left, right), H7)

        def Apow(value, exponent):
            result = {0: field.one}
            factor = urem(field, value, H7)
            power = exponent
            while power:
                if power & 1:
                    result = Amul(result, factor)
                factor = Amul(factor, factor)
                power //= 2
            return result

        def Ainv(value):
            return uinverse_mod(field, urem(field, value, H7), H7)

        def Aconstant(value):
            return {} if field.iszero(value) else {0: value}

        a_element = {1: field.one}
        D_element = Aneg(Amul(V, Ainv(U)))

        def evaluate_bivariate(poly, D_value):
            result = {}
            max_a = max((ia for ia, _ in poly), default=0)
            max_D = max((iD for _, iD in poly), default=0)
            a_powers = [
                Apow(a_element, exponent) for exponent in range(max_a + 1)
            ]
            D_powers = [
                Apow(D_value, exponent) for exponent in range(max_D + 1)
            ]
            for (ia, iD), coefficient in poly.items():
                term = Amul(a_powers[ia], D_powers[iD])
                term = Amul(Aconstant(coefficient), term)
                result = Aadd(result, term)
            return result

        # Verify the subresultant reconstruction before using it.
        require(
            not evaluate_bivariate(F, D_element)
            and not evaluate_bivariate(G, D_element),
            f"{name}/{pair}: compatibility reconstruction failed",
        )

        left, right = pair
        W_numerator = bsub(
            field,
            bmul(field, components[left][1], components[right][2]),
            bmul(field, components[right][1], components[left][2]),
        )
        h_numerator = bsub(
            field,
            bmul(field, components[right][0], components[left][2]),
            bmul(field, components[left][0], components[right][2]),
        )
        delta_inverse = Ainv(delta)
        W_element = Amul(
            evaluate_bivariate(W_numerator, D_element), delta_inverse
        )
        h_element = Amul(
            evaluate_bivariate(h_numerator, D_element), delta_inverse
        )

        # The quotient reconstruction must satisfy all independent terminal rows.
        for coefficient_W, coefficient_h, remainder in components:
            residual = evaluate_bivariate(remainder, D_element)
            residual = Aadd(
                residual,
                Amul(
                    evaluate_bivariate(coefficient_W, D_element), W_element
                ),
            )
            residual = Aadd(
                residual,
                Amul(
                    evaluate_bivariate(coefficient_h, D_element), h_element
                ),
            )
            require(
                not residual,
                f"{name}/{pair}: terminal quotient reconstruction failed",
            )

        units = {
            "a": a_element,
            "D": D_element,
            "W": W_element,
            "h": h_element,
        }
        for unit_name, value in units.items():
            require(
                udegree(ugcd(field, H7, value)) == 0,
                f"{name}/{pair}: {unit_name} vanished on a terminal point",
            )

        order12_residues = []
        gcd_degrees = []
        aggregate = H7
        for equation in equations12:
            result = {}
            for (ic, iD, ie, iW), coefficient in equation.items():
                require(
                    (ic + ie) % 2 == 0,
                    f"{name}: order-12 row retained the C-sign coordinate",
                )
                term = Aconstant(
                    field.from_sympy(slow_field.reduce(coefficient))
                )
                term = Amul(term, Apow(a_element, ie))
                term = Amul(term, Apow(D_element, iD))
                term = Amul(term, Apow(W_element, iW))
                term = Amul(term, Apow(h_element, (ic + ie) // 2))
                result = Aadd(result, term)
            order12_residues.append(result)
            aggregate = ugcd(field, aggregate, result)
            gcd_degrees.append(udegree(aggregate))
        require(
            udegree(aggregate) == 0,
            f"{name}/{pair}: terminal points survive every order-12 row",
        )
        pivot_results.append(
            {
                "pair": pair,
                "raw_resultant_degree": udegree(raw_resultant),
                "removed_degree": removed_degree,
                "H7_degree": udegree(H7),
                "H7_digest": udigest(H7),
                "subresultant_U_degree": udegree(U),
                "order12_gcd_degrees": tuple(gcd_degrees),
                "order12_residue_degrees": tuple(
                    udegree(value) for value in order12_residues
                ),
                "order12_residue_digests": tuple(
                    udigest(value) for value in order12_residues
                ),
            }
        )

    require(
        exact_carriers
        and all(carrier == exact_carriers[0] for carrier in exact_carriers[1:]),
        f"{name}: the six saturated terminal carriers differ",
    )

    # The exact gcd of all six pivot minors is one, so the six charts cover
    # every terminal point.  Their exact monic carriers coincide, and each
    # carrier is coprime to its own pivot minor.  Consequently the first chart
    # already contains the entire terminal scheme; one good-reduction quotient
    # computation there closes every pair-field pivot at once.

    residue_result = None
    if quotient_mode == "good_residue":
        require(residue_payload is not None, f"{name}: missing residue payload")
        pair, H7, delta, F, G = residue_payload

        # Choose a degree-one prime of the exact number field at which every
        # input coefficient is integral.  The exact monic H7 is formed before
        # reduction; preserving its degree and squarefreeness makes this a
        # good-reduction certificate, not a heuristic specialization.  The
        # exact H7 is monic and p-integral.  Any nonconstant common factor in
        # characteristic zero may therefore be chosen monic over the local
        # DVR; degree preservation forces its reduction to retain positive
        # degree.  The residue gcd one below rules that out, uniformly in all
        # complex embeddings of the abstract number field.
        elements = []
        for poly in [*exact_carriers, *deltas.values(), F, G]:
            elements.extend(poly.values())
        for component in components:
            for poly in component:
                elements.extend(poly.values())
        equations12_fast = []
        for equation in equations12:
            converted = {
                monomial: field.from_sympy(slow_field.reduce(value))
                for monomial, value in equation.items()
            }
            equations12_fast.append(converted)
            elements.extend(converted.values())

        good_prime = None
        generator_residue = None
        for prime in s.primerange(101, 20000):
            if any(
                int(coefficient.denominator) % prime == 0
                for element in elements + [field.modulus_tail]
                for coefficient in element
            ):
                continue
            modulus = [
                int(coefficient.numerator)
                * pow(int(coefficient.denominator), -1, prime)
                % prime
                for coefficient in field.modulus_tail
            ]
            for residue in range(prime):
                value = (
                    pow(residue, field.degree, prime)
                    + sum(
                        modulus[index] * pow(residue, index, prime)
                        for index in range(field.degree)
                    )
                ) % prime
                derivative = (
                    field.degree * pow(residue, field.degree - 1, prime)
                    + sum(
                        index
                        * modulus[index]
                        * pow(residue, index - 1, prime)
                        for index in range(1, field.degree)
                    )
                ) % prime
                if value == 0 and derivative != 0:
                    good_prime, generator_residue = prime, residue
                    break
            if good_prime is not None:
                break
        require(good_prime is not None, f"{name}: no good degree-one prime")

        prime_field = PrimeField(good_prime)

        def scalar_mod(element):
            total = 0
            power = 1
            for coefficient in element:
                total = (
                    total
                    + int(coefficient.numerator)
                    * pow(int(coefficient.denominator), -1, good_prime)
                    * power
                ) % good_prime
                power = power * generator_residue % good_prime
            return total

        def mod_b(poly):
            return bclean(
                prime_field,
                {monomial: scalar_mod(value) for monomial, value in poly.items()},
            )

        def mod_u(poly):
            return uclean(
                prime_field,
                {degree: scalar_mod(value) for degree, value in poly.items()},
            )

        Hp = mod_u(H7)
        deltap = mod_u(delta)
        Fp, Gp = mod_b(F), mod_b(G)
        require(
            udegree(Hp) == TORIC["degree"],
            f"{name}: exact H7 lost degree modulo the good prime",
        )
        require(
            udegree(ugcd(prime_field, Hp, uderivative(prime_field, Hp))) == 0,
            f"{name}: exact H7 lost squarefreeness modulo the good prime",
        )
        require(
            udegree(ugcd(prime_field, Hp, deltap)) == 0,
            f"{name}: pivot delta is not a residue unit",
        )

        leading_F = b_to_u(bcoefficient_D(Fp, 2))
        leading_G = b_to_u(bcoefficient_D(Gp, 2))
        linear_subresultant = bsub(
            prime_field,
            bmul(prime_field, u_to_b(leading_G), Fp),
            bmul(prime_field, u_to_b(leading_F), Gp),
        )
        require(
            bdegree_D(linear_subresultant) == 1,
            f"{name}: residue subresultant is not linear",
        )
        U = b_to_u(bcoefficient_D(linear_subresultant, 1))
        V = b_to_u(bcoefficient_D(linear_subresultant, 0))

        def Aadd(left, right):
            return urem(prime_field, uadd(prime_field, left, right), Hp)

        def Aneg(value):
            return urem(prime_field, uneg(prime_field, value), Hp)

        def Amul(left, right):
            return urem(prime_field, umul(prime_field, left, right), Hp)

        def Apow(value, exponent):
            result = {0: 1}
            factor = urem(prime_field, value, Hp)
            power = exponent
            while power:
                if power & 1:
                    result = Amul(result, factor)
                factor = Amul(factor, factor)
                power //= 2
            return result

        def Ainv(value):
            return uinverse_mod(
                prime_field, urem(prime_field, value, Hp), Hp
            )

        a_element = {1: 1}
        D_element = Aneg(Amul(V, Ainv(U)))

        def evaluate_bivariate(poly, D_value):
            result = {}
            for (ia, iD), coefficient in poly.items():
                term = {ia: coefficient}
                term = Amul(term, Apow(D_value, iD))
                result = Aadd(result, term)
            return result

        require(
            not evaluate_bivariate(Fp, D_element)
            and not evaluate_bivariate(Gp, D_element),
            f"{name}: residue compatibility reconstruction failed",
        )

        components_p = [tuple(mod_b(poly) for poly in component)
                        for component in components]
        left, right = pair
        W_numerator = bsub(
            prime_field,
            bmul(prime_field, components_p[left][1], components_p[right][2]),
            bmul(prime_field, components_p[right][1], components_p[left][2]),
        )
        h_numerator = bsub(
            prime_field,
            bmul(prime_field, components_p[right][0], components_p[left][2]),
            bmul(prime_field, components_p[left][0], components_p[right][2]),
        )
        delta_inverse = Ainv(deltap)
        W_element = Amul(
            evaluate_bivariate(W_numerator, D_element), delta_inverse
        )
        h_element = Amul(
            evaluate_bivariate(h_numerator, D_element), delta_inverse
        )

        for coefficient_W, coefficient_h, remainder in components_p:
            residual = evaluate_bivariate(remainder, D_element)
            residual = Aadd(
                residual,
                Amul(evaluate_bivariate(coefficient_W, D_element), W_element),
            )
            residual = Aadd(
                residual,
                Amul(evaluate_bivariate(coefficient_h, D_element), h_element),
            )
            require(not residual, f"{name}: residue terminal row survived")

        units = {
            "a": a_element,
            "D": D_element,
            "W": W_element,
            "h": h_element,
            "delta": deltap,
            "subresultant": U,
        }
        for unit_name, value in units.items():
            inverse = Ainv(value)
            require(
                Amul(value, inverse) == {0: 1},
                f"{name}: residue {unit_name} inverse failed",
            )

        aggregate = Hp
        gcd_degrees = []
        residue_degrees = []
        residue_digests = []
        for equation in equations12_fast:
            result = {}
            for (ic, iD, ie, iW), coefficient in equation.items():
                require(
                    (ic + ie) % 2 == 0,
                    f"{name}: order-12 row retained the C-sign coordinate",
                )
                term = {0: scalar_mod(coefficient)}
                term = Amul(term, Apow(a_element, ie))
                term = Amul(term, Apow(D_element, iD))
                term = Amul(term, Apow(W_element, iW))
                term = Amul(term, Apow(h_element, (ic + ie) // 2))
                result = Aadd(result, term)
            residue_degrees.append(udegree(result))
            residue_digests.append(udigest(result))
            aggregate = ugcd(prime_field, aggregate, result)
            gcd_degrees.append(udegree(aggregate))
        require(
            udegree(aggregate) == 0,
            f"{name}: good reduction retains an order-12 terminal point",
        )
        residue_result = {
            "pair": pair,
            "prime": good_prime,
            "generator_residue": generator_residue,
            "H7_degree": udegree(Hp),
            "order12_gcd_degrees": tuple(gcd_degrees),
            "order12_residue_degrees": tuple(residue_degrees),
            "order12_residue_digests": tuple(residue_digests),
            "units": tuple(units),
        }

    return {
        "field_degree": slow_field.modulus.degree(),
        "selected_rows": selected_rows,
        "support11": support11,
        "support12": support12,
        "rank11": len(pivots11),
        "rank12": len(pivots12),
        "delta_cover_gcd_degree": udegree(delta_gcd),
        "pivot_results": tuple(pivot_results),
        "quotient_units": ("a", "D", "W", "h"),
        "quotient_mode": quotient_mode,
        "residue_result": residue_result,
    }


leading_v = s.Poly(base.R_universal, v).coeff_monomial(v**5)
require(
    len(s.Poly(base.R_universal, t, v, c, d, e, w, domain=s.QQ).terms()) == 60,
    "universal eliminant term count changed",
)
require(leading_v == -88239118492602, "universal v-leading coefficient changed")
P = s.cancel(base.R_universal / leading_v)
require(s.Poly(P, v).LC() == 1, "full-support eliminant is not v-monic")
T_DEGREE = s.Poly(P, t).degree()
require(T_DEGREE == 10, "full-support t-degree changed")
COEFFICIENTS = [
    s.Poly(P, t).coeff_monomial(t**order)
    for order in range(T_DEGREE + 1)
]
require(s.expand(COEFFICIENTS[0] - base.P5_expr) == 0,
        "fixed quintic changed")
require(
    base.P5.degree() == 5
    and base.P5.is_irreducible
    and s.gcd(base.P5, base.P5.diff()).degree() == 0,
    "fixed-section root-subset taxonomy changed",
)

CERTIFICATES = []
for field_name, field, q0, s0, quotient_mode in (
    ("root", base.root_field, base.root_q0, base.root_s0, "exact"),
    ("pair", base.pair_field, base.pair_q0, base.pair_s0, "good_residue"),
):
    CERTIFICATES.append(
        (
            field_name,
            terminal_field_certificate(
                f"BCDEW_{field_name}",
                field,
                q0,
                s0,
                COEFFICIENTS,
                quotient_mode,
            ),
        )
    )


print("degree-22 full-support BCDEW order-12 Hensel/toric obstruction")
print(f"base_thm2636_sha256={BASE_SHA256}")
print("universal_terms=60")
print(f"universal_degrees=t:{T_DEGREE},v:5")
print(f"constant_v_leading={leading_v}")
print("fixed_quintic=degree5,irreducible,squarefree")
print("proper_factor_subset_sizes_up_to_complement=1,2")
print("terminal_support=E,EW,DE,C,CW,CD,CD^2,C^2E,C^3")
print("terminal_torus_coordinates=a=E/C,D,W,h=C^2")
print("terminal_monomials=a,aW,aD,1,W,D,D^2,ah,h")
print(f"toric_ideal_quadrics={TORIC['ideal_generators']}")
print(
    "toric_regular_lower_simplices="
    + ";".join("".join(map(str, simplex)) for simplex in TORIC["simplices"])
)
print("toric_lower_simplex_determinants=1,1,1,1,1,1,1")
print("original_CDEW_lower_simplex_determinants=2,2,2,2,2,2,2")
print(f"original_CDEW_affine_lattice_index={TORIC['original_index']}")
print(f"faithful_toric_degree={TORIC['degree']}")
for field_name, result in CERTIFICATES:
    print(f"{field_name}_field_degree={result['field_degree']}")
    print(
        f"{field_name}_order11_shape=5x{len(result['support11'])}:"
        f"rank={result['rank11']}:selected_rows={result['selected_rows']}"
    )
    print(
        f"{field_name}_order12_shape=5x{len(result['support12'])}:"
        f"rank={result['rank12']}"
    )
    print(
        f"{field_name}_delta_cover_pairs=01,02,03,12,13,23:"
        f"gcd_degree={result['delta_cover_gcd_degree']}"
    )
    print(f"{field_name}_six_exact_H7_carriers_identical=True")
    print(f"{field_name}_quotient_mode={result['quotient_mode']}")
    print(f"{field_name}_quotient_units={','.join(result['quotient_units'])}")
    for pivot in result["pivot_results"]:
        pair_label = "".join(map(str, pivot["pair"]))
        print(
            f"{field_name}_pivot{pair_label}:"
            f"raw_resultant_degree={pivot['raw_resultant_degree']}:"
            f"pivot_artifact_removed_degree={pivot['removed_degree']}:"
            f"H7_degree={pivot['H7_degree']}:H7_squarefree=True:"
            f"H7_digest={pivot['H7_digest']}"
        )
        if result["quotient_mode"] == "exact":
            print(
                f"{field_name}_pivot{pair_label}:"
                f"linear_subresultant_U_degree="
                f"{pivot['subresultant_U_degree']}:"
                f"order12_residue_degrees="
                + ",".join(map(str, pivot["order12_residue_degrees"]))
                + ":aggregate_gcd_degrees="
                + ",".join(map(str, pivot["order12_gcd_degrees"]))
            )
            print(
                f"{field_name}_pivot{pair_label}:order12_residue_digests="
                + ",".join(pivot["order12_residue_digests"])
            )
    if result["residue_result"] is not None:
        residue = result["residue_result"]
        pair_label = "".join(map(str, residue["pair"]))
        print(
            f"{field_name}_good_residue_pivot={pair_label}:"
            f"prime={residue['prime']}:"
            f"generator={residue['generator_residue']}:simple=True:"
            f"exact_H7_degree_preserved={residue['H7_degree']}"
        )
        print(
            f"{field_name}_good_residue_units="
            + ",".join(residue["units"])
        )
        print(
            f"{field_name}_good_residue_order12_degrees="
            + ",".join(map(str, residue["order12_residue_degrees"]))
            + ":aggregate_gcd_degrees="
            + ",".join(map(str, residue["order12_gcd_degrees"]))
        )
        print(
            f"{field_name}_good_residue_order12_digests="
            + ",".join(residue["order12_residue_digests"])
        )
    print(f"{field_name}_all_six_pivot_charts_closed=True")
    print(f"{field_name}_full_support_factor_trajectory=EMPTY")
print("root_and_pair_factor_trajectories=EMPTY")
print("full_support_BCDEW_eliminant=ABSOLUTELY_IRREDUCIBLE")
print("ALL CHECKS PASSED")
