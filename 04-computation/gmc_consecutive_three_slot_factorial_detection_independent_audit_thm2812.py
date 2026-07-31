from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb, factorial


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


class Poly:
    """A tiny exact Z[n,t] implementation, independent of a CAS."""

    def __init__(self, terms: dict[tuple[int, int], int | Fraction] | None = None):
        clean: dict[tuple[int, int], Fraction] = {}
        for exponent, coefficient in (terms or {}).items():
            value = Fraction(coefficient)
            if value:
                clean[exponent] = clean.get(exponent, Fraction(0)) + value
        self.terms = {exponent: value for exponent, value in clean.items() if value}

    @staticmethod
    def coerce(value: Poly | int | Fraction) -> Poly:
        return value if isinstance(value, Poly) else Poly({(0, 0): value})

    def __add__(self, other: Poly | int | Fraction) -> Poly:
        other = Poly.coerce(other)
        answer = dict(self.terms)
        for exponent, coefficient in other.terms.items():
            answer[exponent] = answer.get(exponent, Fraction(0)) + coefficient
        return Poly(answer)

    __radd__ = __add__

    def __neg__(self) -> Poly:
        return Poly({exponent: -coefficient for exponent, coefficient in self.terms.items()})

    def __sub__(self, other: Poly | int | Fraction) -> Poly:
        return self + (-Poly.coerce(other))

    def __rsub__(self, other: Poly | int | Fraction) -> Poly:
        return Poly.coerce(other) - self

    def __mul__(self, other: Poly | int | Fraction) -> Poly:
        other = Poly.coerce(other)
        answer: dict[tuple[int, int], Fraction] = {}
        for (an, at), ac in self.terms.items():
            for (bn, bt), bc in other.terms.items():
                exponent = (an + bn, at + bt)
                answer[exponent] = answer.get(exponent, Fraction(0)) + ac * bc
        return Poly(answer)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> Poly:
        require(exponent >= 0, "negative polynomial exponent")
        answer = Poly.coerce(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                answer = answer * base
            base = base * base
            power //= 2
        return answer

    def coefficient_in_t(self, exponent: int) -> Poly:
        return Poly(
            {
                (n_degree, 0): coefficient
                for (n_degree, t_degree), coefficient in self.terms.items()
                if t_degree == exponent
            }
        )

    def max_t_degree(self) -> int:
        return max((t_degree for _, t_degree in self.terms), default=-1)

    def is_zero(self) -> bool:
        return not self.terms


ONE = Poly.coerce(1)
N = Poly({(1, 0): 1})
T = Poly({(0, 1): 1})


def rising_linear(multiplier: int, length: int) -> Poly:
    answer = ONE
    for offset in range(1, length + 1):
        answer *= multiplier * N + offset
    return answer


def cleared_reduced_moment(power: int) -> Poly:
    """D^power F_power after F_1=0, with no rational functions."""
    denominator = (N + 1) * (N + 2)
    z_numerator = -(1 + (N + 1) * T)
    answer = Poly()
    for b_count in range(power + 1):
        for c_count in range(power - b_count + 1):
            a_count = power - b_count - c_count
            multinomial = (
                factorial(power)
                // (factorial(a_count) * factorial(b_count) * factorial(c_count))
            )
            shift = b_count + 2 * c_count
            answer += (
                multinomial
                * T**b_count
                * z_numerator**c_count
                * denominator ** (power - c_count)
                * rising_linear(power, shift)
            )
    return answer


D = (N + 1) * (N + 2)
G1 = cleared_reduced_moment(1)
G2 = cleared_reduced_moment(2)
G3 = cleared_reduced_moment(3)
require(G1.is_zero(), "cleared first moment did not vanish")

q2 = 2 * (N + 1) ** 2 * (2 * N + 1)
q1 = 6 * (N + 1) * (2 * N + 1)
q0 = 9 * N + 5
Q = q2 * T**2 + q1 * T + q0
require(((N + 1) * G2 - D**2 * Q).is_zero(), "second conic identity failed")

A = 36 * N**4 + 57 * N**3 + 15 * N**2 - 9 * N - 3
B = 52 * N**3 + 49 * N**2 + 8 * N - 1
K = (N + 1) ** 2 * (2 * N + 1)

# The claimed congruence is F_3 == -(A t+B)/K modulo Q.  Clear D^3 K,
# then test divisibility of the resulting cubic by Q without polynomial
# division.  If C=sum c_i t^i and Q=q_2 t^2+q_1 t+q_0, divisibility is
# equivalent to the two cross-multiplied identities below.
cubic = K * G3 + D**3 * (A * T + B)
require(cubic.max_t_degree() <= 3, "unexpected degree in cleared cubic")
c0, c1, c2, c3 = (cubic.coefficient_in_t(index) for index in range(4))
linear_gate = c1 * q2**2 - c3 * q0 * q2 - (c2 * q2 - c3 * q1) * q1
constant_gate = c0 * q2**2 - (c2 * q2 - c3 * q1) * q0
require(linear_gate.is_zero(), "linear pseudo-remainder gate failed")
require(constant_gate.is_zero(), "constant pseudo-remainder gate failed")

# For Q=q_2 t^2+q_1 t+q_0 and ell=A t+B,
# Res_t(Q,ell)=q_2 B^2-q_1 A B+q_0 A^2.
resultant = q2 * B**2 - q1 * A * B + q0 * A**2
positive_factor = (
    16 * N**5
    + 344 * N**4
    + 865 * N**3
    + 731 * N**2
    + 247 * N
    + 29
)
require(
    (resultant - (N + 1) ** 4 * positive_factor).is_zero(),
    "resultant factorization failed",
)


@dataclass(frozen=True)
class Quadratic:
    """a+b*r in Q[r]/(r^2+d), so r^2=-d."""

    a: Fraction
    b: Fraction
    d: int

    @staticmethod
    def rational(value: int | Fraction, d: int) -> Quadratic:
        return Quadratic(Fraction(value), Fraction(0), d)

    def __add__(self, other: Quadratic | int | Fraction) -> Quadratic:
        other = (
            other
            if isinstance(other, Quadratic)
            else Quadratic.rational(other, self.d)
        )
        require(self.d == other.d, "quadratic fields differ")
        return Quadratic(self.a + other.a, self.b + other.b, self.d)

    __radd__ = __add__

    def __neg__(self) -> Quadratic:
        return Quadratic(-self.a, -self.b, self.d)

    def __sub__(self, other: Quadratic | int | Fraction) -> Quadratic:
        return self + (-other if isinstance(other, Quadratic) else -Fraction(other))

    def __mul__(self, other: Quadratic | int | Fraction) -> Quadratic:
        other = (
            other
            if isinstance(other, Quadratic)
            else Quadratic.rational(other, self.d)
        )
        require(self.d == other.d, "quadratic fields differ")
        return Quadratic(
            self.a * other.a - self.d * self.b * other.b,
            self.a * other.b + self.b * other.a,
            self.d,
        )

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> Quadratic:
        require(exponent >= 0, "negative quadratic exponent")
        answer = Quadratic.rational(1, self.d)
        base = self
        power = exponent
        while power:
            if power & 1:
                answer = answer * base
            base = base * base
            power //= 2
        return answer

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0


def rising_integer(base: int, length: int) -> int:
    answer = 1
    for offset in range(1, length + 1):
        answer *= base + offset
    return answer


def normalized_moment_at_root(
    power: int, n_value: int, t_value: Quadratic, z_value: Quadratic
) -> Quadratic:
    answer = Quadratic.rational(0, t_value.d)
    for b_count in range(power + 1):
        for c_count in range(power - b_count + 1):
            a_count = power - b_count - c_count
            multinomial = (
                factorial(power)
                // (factorial(a_count) * factorial(b_count) * factorial(c_count))
            )
            answer += (
                multinomial
                * t_value**b_count
                * z_value**c_count
                * rising_integer(power * n_value, b_count + 2 * c_count)
            )
    return answer


# Independent exact hostile evaluations of both conjugate sharp lines,
# including n=0 and a very large translate.  No large factorials are formed.
tested_translates = (0, 1, 2, 3, 7, 31, 257, 10**6)
for n_value in tested_translates:
    d_value = 2 * n_value + 1
    for sign in (-1, 1):
        t_value = Quadratic(
            Fraction(-3, 2 * (n_value + 1)),
            Fraction(sign, 2 * (n_value + 1) * d_value),
            d_value,
        )
        z_value = Quadratic(
            Fraction(1, 2 * (n_value + 1) * (n_value + 2)),
            Fraction(-sign, 2 * d_value * (n_value + 1) * (n_value + 2)),
            d_value,
        )
        moments = tuple(
            normalized_moment_at_root(power, n_value, t_value, z_value)
            for power in (1, 2, 3)
        )
        require(moments[0].is_zero(), f"F1 root failure at n={n_value}")
        require(moments[1].is_zero(), f"F2 root failure at n={n_value}")
        require(not moments[2].is_zero(), f"F3 hostile failure at n={n_value}")


def bivariate_multiply(
    left: dict[tuple[int, int], Quadratic],
    right: dict[tuple[int, int], Quadratic],
) -> dict[tuple[int, int], Quadratic]:
    d_value = next(iter(left.values())).d
    zero = Quadratic.rational(0, d_value)
    answer: dict[tuple[int, int], Quadratic] = {}
    for (lz, lw), lc in left.items():
        for (rz, rw), rc in right.items():
            exponent = (lz + rz, lw + rw)
            answer[exponent] = answer.get(exponent, zero) + lc * rc
    return {exponent: coefficient for exponent, coefficient in answer.items() if not coefficient.is_zero()}


def bivariate_power(
    polynomial: dict[tuple[int, int], Quadratic], exponent: int
) -> dict[tuple[int, int], Quadratic]:
    d_value = next(iter(polynomial.values())).d
    answer = {(0, 0): Quadratic.rational(1, d_value)}
    for _ in range(exponent):
        answer = bivariate_multiply(answer, polynomial)
    return answer


def gaussian_contraction(
    polynomial: dict[tuple[int, int], Quadratic],
) -> Quadratic:
    d_value = next(iter(polynomial.values())).d
    answer = Quadratic.rational(0, d_value)
    for (z_degree, w_degree), coefficient in polynomial.items():
        if z_degree == w_degree:
            answer += factorial(z_degree) * coefficient
    return answer


# Directly expand the two-charge polynomial in independent Z,W variables
# and apply E[Z^a W^b]=delta_(a,b) a!.  This audits the convention and the
# binomial channel formula independently of the one-variable derivation.
gaussian_cases = 0
for n_value in (1, 2, 7):
    d_value = 2 * n_value + 1
    t_value = Quadratic(
        Fraction(-3, 2 * (n_value + 1)),
        Fraction(1, 2 * (n_value + 1) * d_value),
        d_value,
    )
    z_value = Quadratic(
        Fraction(1, 2 * (n_value + 1) * (n_value + 2)),
        Fraction(-1, 2 * d_value * (n_value + 1) * (n_value + 2)),
        d_value,
    )
    a_value = 2
    lift = {
        (0, 1): Quadratic.rational(a_value, d_value),
        (n_value, n_value - 1): Quadratic.rational(1, d_value),
        (n_value + 1, n_value): t_value,
        (n_value + 2, n_value + 1): z_value,
    }
    radial_moments = {
        power: normalized_moment_at_root(power, n_value, t_value, z_value)
        * factorial(power * n_value)
        for power in (1, 2, 3)
    }
    for moment_order in range(1, 7):
        observed = gaussian_contraction(bivariate_power(lift, moment_order))
        if moment_order % 2:
            expected = Quadratic.rational(0, d_value)
        else:
            half = moment_order // 2
            expected = (
                comb(moment_order, half)
                * a_value**half
                * radial_moments[half]
            )
        require(
            observed == expected,
            f"direct Gaussian contraction failed at n={n_value}, m={moment_order}",
        )
        gaussian_cases += 1

# Boundary hostile controls: after L(A s^u+B s^v)=0 and A=1,
# the second moment is the exact positive integral of
# (s^u-(u!/v!)s^v)^2.
boundary_pairs = 0
for u in range(0, 41):
    for v in range(u + 1, 42):
        ratio = Fraction(factorial(u), factorial(v))
        second = (
            factorial(2 * u)
            - 2 * ratio * factorial(u + v)
            + ratio**2 * factorial(2 * v)
        )
        require(second > 0, f"boundary positivity failed at {(u, v)}")
        boundary_pairs += 1

print("THM-2812 INDEPENDENT HOSTILE AUDIT (NO CAS)")
print("cleared_first_moment=PASS")
print("cleared_second_conic=PASS")
print("division_free_cubic_remainder=PASS")
print(
    "resultant=(n+1)^4*(16*n^5+344*n^4+865*n^3+"
    "731*n^2+247*n+29)"
)
print("sharp_projective_lines=both conjugates exact")
print(f"hostile_translates={tested_translates}")
print(f"boundary_two_slot_pairs={boundary_pairs}")
print(f"direct_gaussian_contraction_cases={gaussian_cases}")
print("scope_guard=genuinely two-charge lift only (H != 0)")
print("ALL INDEPENDENT EXACT CHECKS PASSED")
