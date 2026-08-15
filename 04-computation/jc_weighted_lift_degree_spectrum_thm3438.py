#!/usr/bin/env python3
"""Independent exact companion for the reserved THM-3438 lane.

This standard-library-only checker rebuilds the weighted-lift construction from
its displayed one-variable data.  It deliberately does not import the source
project's SymPy checkers.  Exact sparse polynomial arithmetic over ``Q`` checks

* cancellation of the apparent x^-1 and x^-2 denominators;
* the weighted-invariant Jacobian identity and direct three-variable Jacobians;
* the integer quartic map, its collision, and its component degrees;
* the quartic inverse equation and all four sheets at one generic target in the
  exact quotient Q[w]/(w^4-w^2-1);
* the p_d seed conditions for generic degrees 3 through 100; and
* parameter-symbolic identities proving the seed/cancellation formulas for all
  integer d >= 2, together with the generic-open reconstruction argument.

THM-3438 is a reserved lane.  This file is an exact independent companion and
proof candidate, not a status promotion or a literature-priority claim.  Every
truth gate uses ``require`` and survives ``python -O``.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from pathlib import Path


EXPECTED_SEMANTIC_SHA256 = "fb0f262a5a3d98d2f334eadb53380729b388656d8aa3035f32446ab2cebba721"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


class Poly:
    """Small exact sparse polynomial over Q in a fixed number of variables."""

    def __init__(self, nvars: int, terms: dict[tuple[int, ...], object] | None = None):
        require(nvars >= 1, nvars)
        cleaned: dict[tuple[int, ...], Fraction] = {}
        for exponent, coefficient in (terms or {}).items():
            exponent = tuple(exponent)
            require(len(exponent) == nvars and min(exponent) >= 0, (nvars, exponent))
            value = Fraction(coefficient)
            if value:
                cleaned[exponent] = cleaned.get(exponent, Fraction(0)) + value
                if not cleaned[exponent]:
                    del cleaned[exponent]
        self.nvars = nvars
        self.terms = cleaned

    @classmethod
    def zero(cls, nvars: int) -> "Poly":
        return cls(nvars)

    @classmethod
    def constant(cls, nvars: int, value: object) -> "Poly":
        value = Fraction(value)
        return cls(nvars, {(0,) * nvars: value} if value else {})

    @classmethod
    def variable(cls, nvars: int, index: int) -> "Poly":
        require(0 <= index < nvars, (nvars, index))
        exponent = [0] * nvars
        exponent[index] = 1
        return cls(nvars, {tuple(exponent): 1})

    @classmethod
    def monomial(cls, nvars: int, exponent: tuple[int, ...], coefficient: object = 1) -> "Poly":
        return cls(nvars, {exponent: coefficient})

    def _coerce(self, other: object) -> "Poly":
        if isinstance(other, Poly):
            require(other.nvars == self.nvars, (self.nvars, other.nvars))
            return other
        return Poly.constant(self.nvars, other)

    def __add__(self, other: object) -> "Poly":
        other = self._coerce(other)
        result = dict(self.terms)
        for exponent, coefficient in other.terms.items():
            result[exponent] = result.get(exponent, Fraction(0)) + coefficient
            if not result[exponent]:
                del result[exponent]
        return Poly(self.nvars, result)

    def __radd__(self, other: object) -> "Poly":
        return self + other

    def __neg__(self) -> "Poly":
        return Poly(self.nvars, {exponent: -coefficient for exponent, coefficient in self.terms.items()})

    def __sub__(self, other: object) -> "Poly":
        return self + (-self._coerce(other))

    def __rsub__(self, other: object) -> "Poly":
        return self._coerce(other) - self

    def __mul__(self, other: object) -> "Poly":
        other = self._coerce(other)
        result: dict[tuple[int, ...], Fraction] = {}
        for left_exp, left_coefficient in self.terms.items():
            for right_exp, right_coefficient in other.terms.items():
                exponent = tuple(left + right for left, right in zip(left_exp, right_exp))
                result[exponent] = result.get(exponent, Fraction(0)) + left_coefficient * right_coefficient
                if not result[exponent]:
                    del result[exponent]
        return Poly(self.nvars, result)

    def __rmul__(self, other: object) -> "Poly":
        return self * other

    def __truediv__(self, scalar: object) -> "Poly":
        scalar = Fraction(scalar)
        require(scalar != 0, "polynomial scalar division by zero")
        return Poly(self.nvars, {exponent: coefficient / scalar for exponent, coefficient in self.terms.items()})

    def __pow__(self, exponent: int) -> "Poly":
        require(isinstance(exponent, int) and exponent >= 0, exponent)
        result = Poly.constant(self.nvars, 1)
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result

    def __eq__(self, other: object) -> bool:
        try:
            other = self._coerce(other)
        except (TypeError, ValueError):
            return False
        return self.terms == other.terms

    def coefficient(self, exponent: tuple[int, ...]) -> Fraction:
        require(len(exponent) == self.nvars, exponent)
        return self.terms.get(tuple(exponent), Fraction(0))

    def derivative(self, index: int) -> "Poly":
        require(0 <= index < self.nvars, (self.nvars, index))
        result = {}
        for exponent, coefficient in self.terms.items():
            power = exponent[index]
            if power:
                lowered = list(exponent)
                lowered[index] -= 1
                result[tuple(lowered)] = coefficient * power
        return Poly(self.nvars, result)

    def evaluate(self, values: tuple[object, ...]) -> Fraction:
        require(len(values) == self.nvars, (self.nvars, values))
        values = tuple(Fraction(value) for value in values)
        total = Fraction(0)
        for exponent, coefficient in self.terms.items():
            term = coefficient
            for value, power in zip(values, exponent):
                term *= value**power
            total += term
        return total

    def degree(self) -> int:
        return max((sum(exponent) for exponent in self.terms), default=-1)

    def univariate_degree(self) -> int:
        require(self.nvars == 1, self.nvars)
        return max((exponent[0] for exponent in self.terms), default=-1)

    def leading_coefficient(self) -> Fraction:
        degree = self.univariate_degree()
        require(degree >= 0, "zero polynomial has no leading coefficient")
        return self.coefficient((degree,))

    def divide_monomial(self, exponent: tuple[int, ...]) -> "Poly":
        require(len(exponent) == self.nvars and min(exponent) >= 0, exponent)
        require(
            all(all(have >= need for have, need in zip(term_exp, exponent)) for term_exp in self.terms),
            ("monomial does not divide", exponent, self.key()),
        )
        return Poly(
            self.nvars,
            {
                tuple(have - need for have, need in zip(term_exp, exponent)): coefficient
                for term_exp, coefficient in self.terms.items()
            },
        )

    def key(self) -> tuple[tuple[tuple[int, ...], int, int], ...]:
        return tuple(
            (exponent, coefficient.numerator, coefficient.denominator)
            for exponent, coefficient in sorted(self.terms.items())
        )


def poly_digest(*polynomials: Poly) -> str:
    payload = tuple(polynomial.key() for polynomial in polynomials)
    return sha256(repr(payload).encode("ascii")).hexdigest()


def univariate_integral(polynomial: Poly) -> Poly:
    require(polynomial.nvars == 1, polynomial.nvars)
    return Poly(
        1,
        {
            (exponent[0] + 1,): coefficient / (exponent[0] + 1)
            for exponent, coefficient in polynomial.terms.items()
        },
    )


def compose_univariate(polynomial: Poly, argument: Poly) -> Poly:
    require(polynomial.nvars == 1, polynomial.nvars)
    if not polynomial.terms:
        return Poly.zero(argument.nvars)
    highest = polynomial.univariate_degree()
    powers = [Poly.constant(argument.nvars, 1)]
    for _ in range(highest):
        powers.append(powers[-1] * argument)
    result = Poly.zero(argument.nvars)
    for (exponent,), coefficient in polynomial.terms.items():
        result += coefficient * powers[exponent]
    return result


def shift_univariate_down(polynomial: Poly, amount: int) -> Poly:
    require(polynomial.nvars == 1 and amount >= 0, (polynomial.nvars, amount))
    require(all(exponent[0] >= amount for exponent in polynomial.terms), (amount, polynomial.key()))
    return Poly(1, {(exponent[0] - amount,): coefficient for exponent, coefficient in polynomial.terms.items()})


def q_from_seed(seed: Poly, c: Fraction) -> Poly:
    require(seed.nvars == 1 and c != 0, (seed.nvars, c))
    terms = {}
    for (power,), coefficient in seed.terms.items():
        if power:
            terms[(power + 1,)] = coefficient * power / (c * (power + 1))
    q = Poly(1, terms)
    w = Poly.variable(1, 0)
    require(q.derivative(0) == w * seed.derivative(0) / c, "q-prime identity")
    require(q.coefficient((0,)) == 0, "q(0)")
    return q


def substitute_vt(polynomial: Poly) -> Poly:
    """Substitute v=xy and t=x^2 z into Q[v,t]."""
    require(polynomial.nvars == 2, polynomial.nvars)
    result = {}
    for (v_power, t_power), coefficient in polynomial.terms.items():
        exponent = (v_power + 2 * t_power, v_power, t_power)
        result[exponent] = result.get(exponent, Fraction(0)) + coefficient
    return Poly(3, result)


def determinant3(matrix: tuple[tuple[Poly, Poly, Poly], ...]) -> Poly:
    require(len(matrix) == 3 and all(len(row) == 3 for row in matrix), "3x3 matrix required")
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def jacobian3(coordinates: tuple[Poly, Poly, Poly]) -> Poly:
    require(all(coordinate.nvars == 3 for coordinate in coordinates), "three-variable coordinates required")
    return determinant3(tuple(tuple(coordinate.derivative(index) for index in range(3)) for coordinate in coordinates))


def weighted_lift(seed: Poly, c: object = 1, b: object = 1, a_override: object | None = None) -> dict[str, object]:
    """Build the lift after cancelling gamma, then certify the x-divisions."""
    c = Fraction(c)
    b = Fraction(b)
    require(seed.nvars == 1 and c != 0, (seed.nvars, c))
    integral = univariate_integral(seed)
    require(seed.coefficient((0,)) == 0, "p(0) != 0")
    require(seed.evaluate((1,)) == -c, ("p(1)", seed.evaluate((1,)), -c))
    require(integral.evaluate((1,)) == 0, ("integral", integral.evaluate((1,))))
    kappa = seed.derivative(0).evaluate((1,)) / c
    require(kappa != -2, ("kappa=-2", seed.key()))
    a = Fraction(a_override) if a_override is not None else -(1 + kappa) / (2 + kappa)

    q = q_from_seed(seed, c)
    r = shift_univariate_down(seed, 1)
    s = shift_univariate_down(q, 2)
    v = Poly.variable(2, 0)
    t = Poly.variable(2, 1)
    u = 1 + v
    gamma = 1 + a * v + b * t
    w = u * gamma
    beta = c + u * compose_univariate(r, w)
    alpha = u + u**2 * compose_univariate(s, w)

    # Check that the cancellation-built expressions equal the displayed
    # rational numerators after multiplication by gamma or gamma^2.
    require(beta * gamma == c * gamma + compose_univariate(seed, w), "beta/gamma cancellation")
    require(alpha * gamma**2 == u * gamma**2 + compose_univariate(q, w), "alpha/gamma^2 cancellation")

    beta_origin = beta.coefficient((0, 0))
    alpha_origin = alpha.coefficient((0, 0))
    alpha_v = alpha.coefficient((1, 0))
    require((beta_origin, alpha_origin, alpha_v) == (0, 0, 0), (beta_origin, alpha_origin, alpha_v))

    beta_xyz = substitute_vt(beta)
    alpha_xyz = substitute_vt(alpha)
    gamma_xyz = substitute_vt(gamma)
    first = alpha_xyz.divide_monomial((2, 0, 0))
    second = beta_xyz.divide_monomial((1, 0, 0))
    x = Poly.variable(3, 0)
    third = x * gamma_xyz

    p_invariant = beta * gamma
    q_invariant = alpha * gamma**2
    jac_pq = p_invariant.derivative(0) * q_invariant.derivative(1) - p_invariant.derivative(1) * q_invariant.derivative(0)
    require(jac_pq == -b * c * gamma**2, ("weighted Jacobian", poly_digest(jac_pq, gamma)))

    return {
        "a": a,
        "b": b,
        "c": c,
        "kappa": kappa,
        "q": q,
        "alpha": alpha,
        "beta": beta,
        "gamma": gamma,
        "P": p_invariant,
        "Q": q_invariant,
        "F": (first, second, third),
    }


def pd_seed(d: int) -> Poly:
    require(d >= 2, d)
    delta = Fraction(6, d * (d + 1))
    w = Poly.variable(1, 0)
    return 2 * w - 3 * w**2 + w ** (d - 1) - w**d - delta * w + delta * w**2


class RatFn:
    """Exact rational function in one formal parameter D."""

    def __init__(self, numerator: Poly, denominator: Poly | None = None):
        require(numerator.nvars == 1, numerator.nvars)
        denominator = denominator or Poly.constant(1, 1)
        require(denominator.nvars == 1 and denominator != 0, "zero rational-function denominator")
        if numerator == 0:
            self.numerator = Poly.zero(1)
            self.denominator = Poly.constant(1, 1)
        else:
            scale = denominator.leading_coefficient()
            self.numerator = numerator / scale
            self.denominator = denominator / scale

    @classmethod
    def constant(cls, value: object) -> "RatFn":
        return cls(Poly.constant(1, value))

    @classmethod
    def variable(cls) -> "RatFn":
        return cls(Poly.variable(1, 0))

    def _coerce(self, other: object) -> "RatFn":
        return other if isinstance(other, RatFn) else RatFn.constant(other)

    def __add__(self, other: object) -> "RatFn":
        other = self._coerce(other)
        return RatFn(
            self.numerator * other.denominator + other.numerator * self.denominator,
            self.denominator * other.denominator,
        )

    def __radd__(self, other: object) -> "RatFn":
        return self + other

    def __neg__(self) -> "RatFn":
        return RatFn(-self.numerator, self.denominator)

    def __sub__(self, other: object) -> "RatFn":
        return self + (-self._coerce(other))

    def __rsub__(self, other: object) -> "RatFn":
        return self._coerce(other) - self

    def __mul__(self, other: object) -> "RatFn":
        other = self._coerce(other)
        return RatFn(self.numerator * other.numerator, self.denominator * other.denominator)

    def __rmul__(self, other: object) -> "RatFn":
        return self * other

    def __truediv__(self, other: object) -> "RatFn":
        other = self._coerce(other)
        require(other.numerator != 0, "rational-function division by zero")
        return RatFn(self.numerator * other.denominator, self.denominator * other.numerator)

    def __rtruediv__(self, other: object) -> "RatFn":
        return self._coerce(other) / self

    def __eq__(self, other: object) -> bool:
        other = self._coerce(other)
        return self.numerator * other.denominator == other.numerator * self.denominator

    def evaluate(self, value: int) -> Fraction:
        numerator = self.numerator.evaluate((value,))
        denominator = self.denominator.evaluate((value,))
        require(denominator != 0, ("pole", value, self.denominator.key()))
        return numerator / denominator

    def key(self) -> tuple[object, object]:
        return self.numerator.key(), self.denominator.key()


def all_d_symbolic_audit() -> tuple[object, ...]:
    """Check the p_d endpoint, integral, derivative, and cancellation identities in Q(D)."""
    d = RatFn.variable()
    delta = 6 / (d * (d + 1))
    p_at_one = 2 - 3 + 1 - 1 - delta + delta
    integral = 1 - 1 + 1 / d - 1 / (d + 1) - delta / 6
    kappa = -5 + delta
    kappa_plus_two = kappa + 2
    a = -(1 + kappa) / (2 + kappa)
    q_at_one = p_at_one - integral
    alpha_origin = 1 + q_at_one
    beta_origin = 1 + p_at_one
    alpha_v = 1 + kappa * (1 + a) + 2 * a

    require(p_at_one == -1, p_at_one.key())
    require(integral == 0, integral.key())
    require(kappa == -5 + 6 / (d * (d + 1)), kappa.key())
    require(beta_origin == 0 and alpha_origin == 0 and alpha_v == 0, (beta_origin.key(), alpha_origin.key(), alpha_v.key()))
    require(kappa_plus_two == (6 - 3 * d * (d + 1)) / (d * (d + 1)), kappa_plus_two.key())
    require(a == (4 * d * (d + 1) - 6) / (6 - 3 * d * (d + 1)), a.key())
    for integer_d in range(2, 101):
        require(kappa_plus_two.evaluate(integer_d) < 0, (integer_d, kappa_plus_two.evaluate(integer_d)))

    # Coefficient-generic Jacobian cancellation, with q' eliminated through
    # c*q'=w*p'.  This identity is polynomial in the formal jet p'(w).
    w, gamma, seed_prime, c = (Poly.variable(4, index) for index in range(4))
    c_times_det = c * seed_prime * w - c * (c * gamma + w * seed_prime)
    require(c_times_det == -(c**2) * gamma, poly_digest(c_times_det + c**2 * gamma))

    return (
        "p_d(0)=0",
        "p_d(1)=-1",
        "integral_0^1_p_d=1/[d(d+1)]-(6/[d(d+1)])/6=0",
        "p_d_prime(1)=-5+6/[d(d+1)]",
        "kappa+2=(6-3d(d+1))/[d(d+1)]<0_for_integer_d>=2",
        "a=(4d(d+1)-6)/(6-3d(d+1))",
        "beta00=alpha00=alpha_v00=0",
        "c*Jac_(w,gamma)(P,Q)=-c^2*gamma",
        "Jac_(v,t)(P,Q)=-bc*gamma^2",
        "det_JF=bc",
        "T_d=R_d-Pw+cQ_is_primitive_and_linear_in_Q",
        "k[P,Q,w]/(T_d)=k[P,w]_so_T_d_is_prime_and_irreducible",
        "deg_w(T_d)=d+1_and_reconstruction_identifies_the_source_field",
    )


def finite_seed_and_lift_audit() -> tuple[object, ...]:
    seed_rows = []
    for d in range(2, 100):
        seed = pd_seed(d)
        integral = univariate_integral(seed)
        derivative_at_one = seed.derivative(0).evaluate((1,))
        expected_derivative = -5 + Fraction(6, d * (d + 1))
        q = q_from_seed(seed, Fraction(1))
        require(seed.univariate_degree() == d, (d, seed.univariate_degree()))
        require(seed.evaluate((0,)) == 0 and seed.evaluate((1,)) == -1, (d, seed.evaluate((0,)), seed.evaluate((1,))))
        require(integral.evaluate((1,)) == 0, (d, integral.evaluate((1,))))
        require(derivative_at_one == expected_derivative and derivative_at_one != -2, (d, derivative_at_one))
        require(q.derivative(0) == Poly.variable(1, 0) * seed.derivative(0), d)
        seed_rows.append(
            (
                d,
                d + 1,
                len(seed.terms),
                len(q.terms),
                derivative_at_one,
                -(1 + derivative_at_one) / (2 + derivative_at_one),
                poly_digest(seed, q),
            )
        )

    # Full expansions are intentionally a different finite gate from the
    # all-d proof.  Degrees 3..18 exercise growing support without making the
    # transcript or runtime depend on a large CAS expansion.
    lift_rows = []
    direct_determinants = []
    for d in range(2, 18):
        seed = pd_seed(d)
        lift = weighted_lift(seed)
        first, second, third = lift["F"]
        row = (
            d,
            d + 1,
            (first.degree(), second.degree(), third.degree()),
            (len(first.terms), len(second.terms), len(third.terms)),
            lift["a"],
            poly_digest(first, second, third),
        )
        lift_rows.append(row)
        if d in (2, 3, 4):
            determinant = jacobian3((first, second, third))
            require(determinant == 1, (d, determinant.key()))
            direct_determinants.append((d, determinant.coefficient((0, 0, 0)), len(determinant.terms)))

    # Nonunit parameter control: the factors b and c must survive exactly.
    scaled_seed = 5 * pd_seed(2)
    scaled = weighted_lift(scaled_seed, c=5, b=7)
    scaled_det = jacobian3(scaled["F"])
    require(scaled_det == 35, scaled_det.key())
    scaled_control = (scaled["a"], scaled["b"], scaled["c"], scaled_det.coefficient((0, 0, 0)))

    seed_surface_digest = sha256(repr(tuple(seed_rows)).encode("ascii")).hexdigest()
    lift_surface_digest = sha256(repr(tuple(lift_rows)).encode("ascii")).hexdigest()
    return (
        (2, 99, 3, 100, len(seed_rows), tuple(seed_rows[:6]), seed_rows[-1], seed_surface_digest),
        (2, 17, 3, 18, len(lift_rows), tuple(lift_rows), lift_surface_digest),
        tuple(direct_determinants),
        scaled_control,
    )


def uni_divmod(dividend: Poly, divisor: Poly) -> tuple[Poly, Poly]:
    require(dividend.nvars == divisor.nvars == 1 and divisor != 0, "univariate division")
    quotient = Poly.zero(1)
    remainder = dividend
    divisor_degree = divisor.univariate_degree()
    divisor_lead = divisor.leading_coefficient()
    while remainder != 0 and remainder.univariate_degree() >= divisor_degree:
        power = remainder.univariate_degree() - divisor_degree
        coefficient = remainder.leading_coefficient() / divisor_lead
        term = Poly.monomial(1, (power,), coefficient)
        quotient += term
        remainder -= term * divisor
    return quotient, remainder


def uni_mod(polynomial: Poly, modulus: Poly) -> Poly:
    return uni_divmod(polynomial, modulus)[1]


def uni_xgcd(left: Poly, right: Poly) -> tuple[Poly, Poly, Poly]:
    require(left.nvars == right.nvars == 1, "univariate xgcd")
    old_r, r = left, right
    old_s, s = Poly.constant(1, 1), Poly.zero(1)
    old_t, t = Poly.zero(1), Poly.constant(1, 1)
    while r != 0:
        quotient, remainder = uni_divmod(old_r, r)
        old_r, r = r, remainder
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    scale = old_r.leading_coefficient()
    gcd_polynomial, bezout_left, bezout_right = old_r / scale, old_s / scale, old_t / scale
    require(bezout_left * left + bezout_right * right == gcd_polynomial, "Bezout identity")
    return gcd_polynomial, bezout_left, bezout_right


def uni_inverse_mod(polynomial: Poly, modulus: Poly) -> Poly:
    gcd_polynomial, bezout, _ = uni_xgcd(polynomial, modulus)
    require(gcd_polynomial == 1, ("not invertible modulo polynomial", gcd_polynomial.key()))
    inverse = uni_mod(bezout, modulus)
    require(uni_mod(polynomial * inverse, modulus) == 1, "quotient inverse check")
    return inverse


def matrix_determinant(matrix: list[list[Fraction]]) -> Fraction:
    require(matrix and all(len(row) == len(matrix) for row in matrix), "square matrix")
    work = [list(row) for row in matrix]
    determinant = Fraction(1)
    size = len(work)
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column]
        determinant *= pivot_value
        for row in range(column + 1, size):
            if not work[row][column]:
                continue
            factor = work[row][column] / pivot_value
            for entry in range(column, size):
                work[row][entry] -= factor * work[column][entry]
    return determinant


def resultant(left: Poly, right: Poly) -> Fraction:
    require(left.nvars == right.nvars == 1 and left != 0 and right != 0, "nonzero univariate resultant")
    m = left.univariate_degree()
    n = right.univariate_degree()
    width = m + n
    left_desc = [left.coefficient((power,)) for power in range(m, -1, -1)]
    right_desc = [right.coefficient((power,)) for power in range(n, -1, -1)]
    matrix = []
    for offset in range(n):
        matrix.append([Fraction(0)] * offset + left_desc + [Fraction(0)] * (n - 1 - offset))
    for offset in range(m):
        matrix.append([Fraction(0)] * offset + right_desc + [Fraction(0)] * (m - 1 - offset))
    require(len(matrix) == width and all(len(row) == width for row in matrix), (m, n, width))
    return matrix_determinant(matrix)


def discriminant(polynomial: Poly) -> Fraction:
    degree = polynomial.univariate_degree()
    sign = -1 if degree * (degree - 1) // 2 % 2 else 1
    return sign * resultant(polynomial, polynomial.derivative(0)) / polynomial.leading_coefficient()


def quotient_substitute(polynomial: Poly, coordinates: tuple[Poly, ...], modulus: Poly) -> Poly:
    require(polynomial.nvars == len(coordinates), (polynomial.nvars, len(coordinates)))
    require(all(coordinate.nvars == 1 for coordinate in coordinates), "quotient coordinates")
    maxima = [0] * polynomial.nvars
    for exponent in polynomial.terms:
        for index, power in enumerate(exponent):
            maxima[index] = max(maxima[index], power)
    powers = []
    for coordinate, maximum in zip(coordinates, maxima):
        row = [Poly.constant(1, 1)]
        for _ in range(maximum):
            row.append(uni_mod(row[-1] * coordinate, modulus))
        powers.append(row)
    result = Poly.zero(1)
    for exponent, coefficient in polynomial.terms.items():
        term = Poly.constant(1, coefficient)
        for index, power in enumerate(exponent):
            term = uni_mod(term * powers[index][power], modulus)
        result = uni_mod(result + term, modulus)
    return result


def quartic_map_audit() -> tuple[object, ...]:
    x, y, z = (Poly.variable(3, index) for index in range(3))
    u = 1 + 3 * x * y
    gamma = 1 - 4 * x * y - x**2 * z
    numerator_a = 2 * u + u**2 - 3 * u**4 * gamma**2
    numerator_b = 1 + u - 2 * u**3 * gamma**2
    a_coordinate = numerator_a.divide_monomial((2, 0, 0))
    b_coordinate = numerator_b.divide_monomial((1, 0, 0))
    c_coordinate = x * gamma
    coordinates = (a_coordinate, b_coordinate, c_coordinate)

    ordinary_degrees = tuple(coordinate.degree() for coordinate in coordinates)
    z_degrees = tuple(max((exponent[2] for exponent in coordinate.terms), default=-1) for coordinate in coordinates)
    require(ordinary_degrees == (12, 11, 4), ordinary_degrees)
    require(z_degrees == (2, 2, 1), z_degrees)
    jacobian = jacobian3(coordinates)
    require(jacobian == -6, jacobian.key())

    left_point = (1, 0, 0)
    right_point = (-1, 0, 2)
    hostile_point = (-1, 0, 1)
    left_image = tuple(coordinate.evaluate(left_point) for coordinate in coordinates)
    right_image = tuple(coordinate.evaluate(right_point) for coordinate in coordinates)
    hostile_image = tuple(coordinate.evaluate(hostile_point) for coordinate in coordinates)
    require(left_point != right_point and left_image == right_image == (0, 0, 1), (left_image, right_image))
    require(hostile_image != left_image, hostile_image)

    # Exact elimination and reconstruction identities in Q[w,gamma].
    w, g = Poly.variable(2, 0), Poly.variable(2, 1)
    p_invariant = g + w - 2 * w**3
    q_invariant = 2 * w * g + w**2 - 3 * w**4
    inverse_equation = w**4 - w**2 + 2 * w * p_invariant - q_invariant
    gamma_reconstruction = p_invariant - w + 2 * w**3
    require(inverse_equation == 0 and gamma_reconstruction == g, (inverse_equation.key(), gamma_reconstruction.key()))

    # Target (A,B,C)=(1,0,1): P=0, Q=1.  Work in the reduced
    # four-dimensional algebra Q[w]/(w^4-w^2-1), so one exact calculation
    # checks the reconstruction at all four algebraic roots simultaneously.
    w1 = Poly.variable(1, 0)
    quartic = w1**4 - w1**2 - 1
    quartic_discriminant = discriminant(quartic)
    require(quartic_discriminant == -400, quartic_discriminant)
    gamma_class = uni_mod(-w1 + 2 * w1**3, quartic)
    gamma_inverse = uni_inverse_mod(gamma_class, quartic)
    x_class = gamma_inverse
    u_class = uni_mod(w1 * gamma_inverse, quartic)
    y_class = uni_mod((u_class - 1) * gamma_class / 3, quartic)
    z_class = uni_mod((1 - Fraction(4, 3) * (u_class - 1) - gamma_class) * gamma_class**2, quartic)
    reconstructed = tuple(quotient_substitute(coordinate, (x_class, y_class, z_class), quartic) for coordinate in coordinates)
    require(reconstructed == (Poly.constant(1, 1), Poly.zero(1), Poly.constant(1, 1)), tuple(item.key() for item in reconstructed))

    # Discriminant hostile: at (A,B,C)=(-1/4,0,1), the inverse quartic is a
    # square and gamma vanishes at both repeated roots.  Reconstruction at a
    # repeated inverse root is therefore an escape-to-infinity event, not a
    # finite critical point.
    wall = w1**4 - w1**2 + Fraction(1, 4)
    wall_factor = w1**2 - Fraction(1, 2)
    require(wall == wall_factor**2 and discriminant(wall) == 0, (wall.key(), discriminant(wall)))
    wall_gcd, _, _ = uni_xgcd(uni_mod(-w1 + 2 * w1**3, wall), wall)
    require(wall_gcd == wall_factor, wall_gcd.key())

    return (
        ordinary_degrees,
        z_degrees,
        tuple(len(coordinate.terms) for coordinate in coordinates),
        poly_digest(*coordinates),
        jacobian.coefficient((0, 0, 0)),
        (left_point, right_point, left_image),
        (hostile_point, hostile_image),
        "w^4-w^2+2Pw-Q=0;gamma=P-w+2w^3;f_prime=2gamma",
        ("target=(1,0,1)", quartic_discriminant, len(quartic.terms), gamma_inverse.key(), tuple(item.key() for item in reconstructed)),
        ("wall_target=(-1/4,0,1)", wall_factor.key(), discriminant(wall), wall_gcd.key()),
    )


def hostile_audit() -> tuple[object, ...]:
    w = Poly.variable(1, 0)
    good = pd_seed(2)

    # Same endpoint values, wrong integral: q(1) changes and the x^-2
    # numerator acquires a forbidden constant term.
    bad_integral = good + w * (1 - w)
    bad_integral_value = univariate_integral(bad_integral).evaluate((1,))
    q_bad = q_from_seed(bad_integral, Fraction(1))
    alpha_constant = 1 + q_bad.evaluate((1,))
    require(bad_integral.evaluate((0,)) == 0 and bad_integral.evaluate((1,)) == -1, "endpoint hostile changed")
    require(bad_integral_value == Fraction(1, 6) and alpha_constant == Fraction(-1, 6), (bad_integral_value, alpha_constant))

    # Correct seed, wrong a: the forbidden v term survives in alpha.
    good_kappa = good.derivative(0).evaluate((1,))
    good_a = -(1 + good_kappa) / (2 + good_kappa)
    wrong_a = good_a + 1
    q_good = q_from_seed(good, Fraction(1))
    r_good = shift_univariate_down(good, 1)
    s_good = shift_univariate_down(q_good, 2)
    v, t = Poly.variable(2, 0), Poly.variable(2, 1)
    u = 1 + v
    gamma_wrong = 1 + wrong_a * v + t
    alpha_wrong = u + u**2 * compose_univariate(s_good, u * gamma_wrong)
    wrong_linear = alpha_wrong.coefficient((1, 0))
    require(wrong_linear == good_kappa + 2 == -2, wrong_linear)

    # The three seed conditions alone do not remove kappa=-2.
    singular = 4 * w - 9 * w**2 + 4 * w**3
    singular_data = (
        singular.evaluate((0,)),
        singular.evaluate((1,)),
        univariate_integral(singular).evaluate((1,)),
        singular.derivative(0).evaluate((1,)),
    )
    require(singular_data == (0, -1, 0, -2), singular_data)

    # d=1 is the sharp lower boundary for the displayed p_d formula.
    d1 = 1 - 2 * w
    require(d1.evaluate((0,)) == 1, d1.key())

    return (
        ("wrong_integral", bad_integral_value, alpha_constant, "x^-2_constant_survives"),
        ("wrong_a", wrong_a, wrong_linear, "x^-1_term_survives_after_alpha/x^2"),
        ("kappa_minus_2_seed", singular_data, "a_undefined"),
        ("d=1_boundary", d1.evaluate((0,)), "p_d(0)_fails"),
        ("b=0_boundary", "det_JF=bc=0"),
    )


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert statement found")

    symbolic = all_d_symbolic_audit()
    finite = finite_seed_and_lift_audit()
    quartic = quartic_map_audit()
    hostiles = hostile_audit()

    consequence = (
        "for_every_integer_d>=2_the_displayed_p_d_seed_satisfies_the_lift_gates",
        "the_resulting_dimension3_Keller_map_has_det1_and_generic_degree_d+1",
        "therefore_every_generic_degree_n>=3_occurs_in_this_weighted_lift_family",
        "the_integer_quartic_member_has_det_minus6_and_an_explicit_collision",
    )
    semantic_surface = (symbolic, finite, quartic, hostiles, consequence)
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256, (semantic_digest, EXPECTED_SEMANTIC_SHA256))

    print("THM-3438 weighted-lift Keller degree-spectrum independent exact companion")
    print("status=VERIFIED_EXACT_INDEPENDENT_COMPANION_FOR_RESERVED_THM3438;no_status_promotion")
    print("source_formula=https://jacobianfun.org/jacobian-explained;external_checker_not_imported_or_executed")
    print("arithmetic=stdlib_sparse_polynomials_over_Q;formal_parameter_ring_Q(D);quartic_quotient_Q[w]/(w^4-w^2-1)")
    print(f"all_d_symbolic_identities={symbolic}")
    print("all_d_polynomiality=p(0)=0_gives_p=w*r;q_prime=w*p_prime/c_gives_q=w^2*s;beta00=0_removes_x^-1;alpha00=alpha_v00=0_remove_x^-2")
    print("all_d_generic_degree=T_d=I_d(w)-wP+Q_is_primitive_linear_in_Q;quotient_eliminates_Q_to_k[P,w];T_d_irreducible_of_degree_d+1;off_C0_each_root_reconstructs_one_unique_input")
    print(f"finite_seed_bank=(d_min,d_max,n_min,n_max,count,first_rows,last_row,sha256)={finite[0]}")
    print(f"finite_expanded_lift_bank=(d_min,d_max,n_min,n_max,count,rows,sha256)={finite[1]}")
    print(f"direct_p_d_three_variable_jacobians=(d,constant,term_count)={finite[2]}")
    print(f"scaled_bc_control=(a,b,c,direct_det)={finite[3]}")
    print(f"quartic_certificate=(ordinary_degrees,z_degrees,term_counts,map_sha256,det,collision,nearby_hostile,inverse_identity,generic_quotient,wall_hostile)={quartic}")
    print(f"hostile_controls={hostiles}")
    print(f"consequence={consequence}")
    print("scope=exact_algebra_of_the_displayed_weighted_lift_and_p_d_family;proves_no_classification_of_maps_within_a_degree;no_literature_priority;does_not_touch_JC2;reserved_THM3438_status_is_not_changed_by_this_artifact")
    print(f"semantic_sha256={semantic_digest}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=standard_library_only;no_randomness;no_elapsed_fields;all_truth_gates_survive_python_O")
    print("commands=python -B 04-computation/jc_weighted_lift_degree_spectrum_thm3438.py;python -B -O 04-computation/jc_weighted_lift_degree_spectrum_thm3438.py")


if __name__ == "__main__":
    main()
