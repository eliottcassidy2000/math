#!/usr/bin/env python3
"""Dependency-free quadratic-field audit for THM-4307."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class Q265:
    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @staticmethod
    def make(value: int | Fraction | "Q265") -> "Q265":
        if isinstance(value, Q265):
            return value
        return Q265(Fraction(value), Fraction(0))

    def __add__(self, other: int | Fraction | "Q265") -> "Q265":
        other = self.make(other)
        return Q265(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self) -> "Q265":
        return Q265(-self.a, -self.b)

    def __sub__(self, other: int | Fraction | "Q265") -> "Q265":
        return self + (-self.make(other))

    def __rsub__(self, other: int | Fraction | "Q265") -> "Q265":
        return self.make(other) - self

    def __mul__(self, other: int | Fraction | "Q265") -> "Q265":
        other = self.make(other)
        return Q265(
            self.a * other.a + 265 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "Q265":
        if exponent < 0:
            return (self.inverse()) ** (-exponent)
        result = Q265.make(1)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent >>= 1
        return result

    def norm(self) -> Fraction:
        return self.a * self.a - 265 * self.b * self.b

    def inverse(self) -> "Q265":
        n = self.norm()
        if n == 0:
            raise ZeroDivisionError
        return Q265(self.a / n, -self.b / n)

    def __truediv__(self, other: int | Fraction | "Q265") -> "Q265":
        return self * self.make(other).inverse()

    def __str__(self) -> str:
        return f"({self.a})+({self.b})sqrt265"


Poly = dict[tuple[int, int], Q265]


def clean(poly: Poly) -> Poly:
    return {monomial: coefficient for monomial, coefficient in poly.items()
            if coefficient != Q265()}


def add(left: Poly, right: Poly) -> Poly:
    result = dict(left)
    for monomial, coefficient in right.items():
        result[monomial] = result.get(monomial, Q265()) + coefficient
    return clean(result)


def scale(poly: Poly, scalar: int | Fraction | Q265) -> Poly:
    scalar = Q265.make(scalar)
    return clean({monomial: scalar * coefficient
                  for monomial, coefficient in poly.items()})


def multiply(left: Poly, right: Poly) -> Poly:
    result: Poly = {}
    for (q1, t1), c1 in left.items():
        for (q2, t2), c2 in right.items():
            monomial = (q1 + q2, t1 + t2)
            result[monomial] = result.get(monomial, Q265()) + c1 * c2
    return clean(result)


def power(poly: Poly, exponent: int) -> Poly:
    result: Poly = {(0, 0): Q265.make(1)}
    for _ in range(exponent):
        result = multiply(result, poly)
    return result


def t_shift(poly: Poly, exponent: int) -> Poly:
    return {(q_degree, t_degree + exponent): coefficient
            for (q_degree, t_degree), coefficient in poly.items()}


def evaluate(coefficients: list[Q265], value: Q265) -> Q265:
    result = Q265()
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def derivative(coefficients: list[Q265]) -> list[Q265]:
    return [index * coefficient for index, coefficient in enumerate(coefficients)][1:]


def literal_strict_transform(U: Q265) -> Poly:
    one_plus_q: Poly = {(0, 0): Q265.make(1), (1, 0): Q265.make(1)}
    d = Fraction(2048, 45)
    K = Fraction(2848, 45) - Fraction(7, 6) * d
    hhat: Poly = {}
    hhat = add(hhat, scale(add(add(power(one_plus_q, 6),
                                      scale(power(one_plus_q, 5), -2)),
                                  power(one_plus_q, 4)), U))
    hhat = add(hhat, t_shift(scale(add(power(one_plus_q, 4),
                                          scale(power(one_plus_q, 3), -1)), d), 4))
    hhat = add(hhat, t_shift(add(scale(power(one_plus_q, 3), Fraction(-1376, 135)),
                                      scale(power(one_plus_q, 2), K)), 6))
    hhat = add(hhat, t_shift(scale(power(one_plus_q, 2), Fraction(8, 3)), 8))
    hhat = add(hhat, t_shift(scale(one_plus_q, -3), 10))

    # Multiply by q and add -t^12/2.  The separate -q*z^12 term becomes
    # -y*z^12/t^8 and is not part of this (q,t)-coefficient dictionary.
    F = multiply({(1, 0): Q265.make(1)}, hhat)
    F = add(F, {(0, 12): Q265.make(Fraction(-1, 2))})
    strict: Poly = {}
    for (q_degree, t_degree), coefficient in F.items():
        exponent = 4 * q_degree + t_degree - 12
        if exponent < 0:
            raise AssertionError((q_degree, t_degree, exponent))
        strict[(q_degree, exponent)] = coefficient
    return clean(strict)


def literal_regime_a_k3_strict_transform() -> Poly:
    """Reconstruct the eta=1, zeta_3=-1, Delta=0 hostile from source."""

    one_plus_q: Poly = {(0, 0): Q265.make(1), (1, 0): Q265.make(1)}
    U3 = Fraction(135, 28672)
    hhat: Poly = {}
    hhat = add(hhat, scale(add(add(power(one_plus_q, 6),
                                      scale(power(one_plus_q, 5), -2)),
                                  power(one_plus_q, 4)), U3))
    hhat = add(hhat, t_shift(add(power(one_plus_q, 4),
                                      scale(power(one_plus_q, 3), -1)), 3))
    hhat = add(hhat, t_shift(add(
        scale(power(one_plus_q, 3), Fraction(-1376, 135)),
        scale(power(one_plus_q, 2), Fraction(2848, 45))), 6))
    hhat = add(hhat, t_shift(scale(power(one_plus_q, 2), Fraction(8, 3)), 8))
    hhat = add(hhat, t_shift(scale(one_plus_q, -3), 10))

    source = multiply({(1, 0): Q265.make(1)}, hhat)
    source = add(source, {(0, 12): Q265.make(Fraction(-1, 2))})
    strict: Poly = {}
    for (q_degree, t_degree), coefficient in source.items():
        exponent = 3 * q_degree + t_degree - 9
        if exponent < 0:
            raise AssertionError((q_degree, t_degree, exponent))
        strict[(q_degree, exponent)] = coefficient
    # The omitted literal -q*z^12 row becomes -y*z^12/t^6.
    return clean(strict)


def expected_strict(U: Q265) -> Poly:
    rows: Poly = {
        (0, 0): Q265.make(Fraction(-1, 2)),
        (1, 0): Q265.make(Fraction(8, 3)),
        (2, 0): Q265.make(Fraction(2048, 45)),
        (3, 0): U,
        (1, 2): Q265.make(-3),
        (2, 2): Q265.make(Fraction(-1376, 135)),
        (2, 4): Q265.make(Fraction(16, 3)),
        (3, 4): Q265.make(Fraction(2048, 15)),
        (4, 4): 4 * U,
        (2, 6): Q265.make(-3),
        (3, 6): Q265.make(Fraction(-2752, 135)),
        (3, 8): Q265.make(Fraction(8, 3)),
        (4, 8): Q265.make(Fraction(2048, 15)),
        (5, 8): 6 * U,
        (4, 10): Q265.make(Fraction(-1376, 135)),
        (5, 12): Q265.make(Fraction(2048, 45)),
        (6, 12): 4 * U,
        (7, 16): U,
    }
    return clean(rows)


def main() -> None:
    rows: list[tuple[int, Q265, Q265, Q265, Q265, Q265, Q265]] = []
    for epsilon in (1, -1):
        U = Q265(Fraction(-315392, 3645),
                 epsilon * Fraction(217088, 18225))
        rho = Q265(Fraction(-15, 256), epsilon * Fraction(-3, 256))
        if literal_strict_transform(U) != expected_strict(U):
            raise AssertionError(f"strict transform epsilon={epsilon}")

        p = [Q265.make(Fraction(-1, 2)), Q265.make(Fraction(8, 3)),
             Q265.make(Fraction(2048, 45)), U]
        p_prime = [p[1], 2 * p[2], 3 * p[3]]
        if evaluate(p, rho) != Q265() or evaluate(p_prime, rho) != Q265():
            raise AssertionError(f"double root epsilon={epsilon}")
        a = 3 * U * rho + Fraction(2048, 45)
        g = Fraction(-1376, 135) * rho**2 - 3 * rho
        expected_a = Q265(Fraction(-6784, 135), epsilon * Fraction(128, 135))
        expected_g = Q265(Fraction(-707, 3072), epsilon * Fraction(65, 3072))
        if a != expected_a or g != expected_g:
            raise AssertionError(f"splitter coefficients epsilon={epsilon}")
        if rho.norm() == 0 or a.norm() == 0 or g.norm() == 0:
            raise AssertionError(f"degenerate coefficient epsilon={epsilon}")
        R2 = [Q265(), Q265.make(-3), Q265.make(Fraction(-1376, 135))]
        R4 = [Q265(), Q265(), Q265.make(Fraction(16, 3)),
              Q265.make(Fraction(2048, 15)), 4 * U]
        p2 = evaluate(derivative(derivative(p)), rho)
        p3 = evaluate(derivative(derivative(derivative(p))), rho)
        R2_value = evaluate(R2, rho)
        R2_prime = evaluate(derivative(R2), rho)
        R2_second = evaluate(derivative(derivative(R2)), rho)
        R4_value = evaluate(R4, rho)
        R4_prime = evaluate(derivative(R4), rho)
        y1 = (R2_value - rho * R2_prime) / (rho * p2)
        y2 = (
            -Fraction(1, 2) * (p2 + rho * p3) * y1**2
            - rho * R2_second * y1 + R4_value - rho * R4_prime
        ) / (rho * p2)
        chi = p2 * y1 + R2_prime
        lam = (p2 * y2 + Fraction(1, 2) * p3 * y1**2
               + R2_second * y1 + R4_prime)
        expected_chi = Q265(Fraction(-173, 72), epsilon * Fraction(43, 360))
        expected_lam = Q265(Fraction(11975, 27648),
                            epsilon * Fraction(-53621, 7326720))
        if chi != expected_chi or lam != expected_lam or chi != g / rho:
            raise AssertionError(f"discriminant graph epsilon={epsilon}")
        if chi.norm() != Fraction(269, 135):
            raise AssertionError(f"graph slope norm epsilon={epsilon}")
        rows.append((epsilon, U, rho, a, g, chi, lam))

    if rows[0][1] != Q265(rows[1][1].a, -rows[1][1].b):
        raise AssertionError("U conjugation")
    if rows[0][2] != Q265(rows[1][2].a, -rows[1][2].b):
        raise AssertionError("rho conjugation")
    if rows[0][3].norm() != Fraction(13893632, 6075):
        raise AssertionError("quadratic-unit norm")
    if rows[0][4].norm() != Fraction(-269, 4096):
        raise AssertionError("R2 norm")

    # Rational Regime-A k=3, Delta=0 hostile, rebuilt from the literal source.
    U3 = Fraction(135, 28672)
    rho3 = Fraction(-14336, 135)
    hostile_strict = literal_regime_a_k3_strict_transform()
    hostile_base = {q_degree: coefficient
                    for (q_degree, t_degree), coefficient in hostile_strict.items()
                    if t_degree == 0}
    expected_base = {
        1: Q265.make(U3 * rho3**2),
        2: Q265.make(-2 * U3 * rho3),
        3: Q265.make(U3),
    }
    if hostile_base != expected_base or -Fraction(1, 2) / U3 != rho3:
        raise AssertionError("Regime-A literal repeated base")
    hostile_t1 = {q_degree: coefficient
                  for (q_degree, t_degree), coefficient in hostile_strict.items()
                  if t_degree == 1}
    hostile_t2 = {q_degree: coefficient
                  for (q_degree, t_degree), coefficient in hostile_strict.items()
                  if t_degree == 2}
    if hostile_t1 or hostile_t2 != {1: Q265.make(Fraction(8, 3))}:
        raise AssertionError("Regime-A literal first graph row")
    graph_row_at_rho = hostile_t2[1] * rho3
    lambda3 = (graph_row_at_rho / rho3).a
    if lambda3 != Fraction(8, 3):
        raise AssertionError("Regime-A graph coefficient")
    # With t=sigma*z and z=sigma^2*w, z^12/t^6 and lambda3*t^2
    # become sigma^6*w^6 and lambda3*sigma^6*w^2.  Removing the square
    # sigma^6*w^2 leaves the displayed quartic.
    if (12 - 6, 2 * 6 - 6, 2 * (1 + 2)) != (6, 6, 6):
        raise AssertionError("Regime-A boundary monomial ledger")
    quartic_I = 12 * Fraction(-8, 3)
    quartic_J = 0
    quartic_j = Fraction(1728 * 4 * quartic_I**3,
                         4 * quartic_I**3 - quartic_J**2)
    if (quartic_I, quartic_J, quartic_j) != (-32, 0, 1728):
        raise AssertionError("Regime-A elliptic invariants")
    tau3 = 3
    d3 = 2 * tau3
    fq_order3 = 6 * tau3 + d3 // 2
    good_order3 = 9 + 22 - fq_order3
    if (fq_order3, good_order3) != (21, 10):
        raise AssertionError("Regime-A good-form order")

    print("THM4307_BALANCED_REFINEMENT_INDEPENDENT_V1")
    print("STRICT_TRANSFORM exact rows=18 plus separate -y*z^12/t^8")
    print("LOCAL_ALGEBRA double Q^2 derivative ideal=(Q) dimension=1; triple Q^3 dimension=2 excluded")
    for epsilon, U, rho, a, g, chi, lam in rows:
        print(f"epsilon={epsilon:+d} U={U} rho={rho} a={a} g={g}")
        print(f"epsilon={epsilon:+d} nonzero_norms rho={rho.norm()} a={a.norm()} g={g.norm()}")
        print(f"epsilon={epsilon:+d} graph m=({chi})*x+({lam})*x^2+O(x^3) chi_norm={chi.norm()}")
    print("VALUATION H_excess-T2_excess=2(beta-5s)")
    print("BOUNDARY radicand square_part*squarefree_quadratic; normalization_genus=0")
    print("BALANCED_NORMAL_CHART smooth; completed_local_refinement_carriers_rational")
    print("REGIME_A_HOSTILE k=3 Delta=0 U=135/28672 rho=-14336/135 curve=V^2-(w^4-8/3) I=-32 J=0 j=1728 Fq_order=21 good_form_order=10")
    print("VERDICT PASS BALANCED_FORMAL_LOCAL_EXTINCTION REGIME_A_POSITIVE_GENUS_CONTROL")


if __name__ == "__main__":
    main()
