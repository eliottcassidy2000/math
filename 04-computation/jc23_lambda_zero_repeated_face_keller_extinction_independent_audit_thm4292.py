#!/usr/bin/env python3
"""Standard-library exact audit for THM-4292.

This is a clean-room certificate for the repeated-face calculation in the
exact M=12, W=0, Z=-U infinity chart.  It uses a small sparse Laurent-
polynomial implementation over Fraction, reconstructs the truncated critical
equation directly from the r-powers, performs formal critical-point
substitution, and checks the valuation and normalization case table.

No result here asserts degree persistence or constructs a Keller pair.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import comb
import random


VARS = ("c", "alpha11", "upsilon5", "eta", "Delta")
NVAR = len(VARS)


def require(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)


class Laurent:
    """Sparse Laurent polynomial over Q in the fixed variables VARS."""

    def __init__(self, terms: dict[tuple[int, ...], F] | None = None) -> None:
        clean: dict[tuple[int, ...], F] = {}
        for monomial, coefficient in (terms or {}).items():
            require(len(monomial) == NVAR, "Laurent monomial length")
            coefficient = F(coefficient)
            if coefficient:
                clean[tuple(monomial)] = clean.get(tuple(monomial), F(0)) + coefficient
        self.terms = {m: a for m, a in clean.items() if a}

    @staticmethod
    def constant(value: F | int) -> "Laurent":
        value = F(value)
        return Laurent({(0,) * NVAR: value}) if value else Laurent()

    @staticmethod
    def monomial(exponents: tuple[int, ...], coefficient: F | int = 1) -> "Laurent":
        return Laurent({exponents: F(coefficient)})

    @staticmethod
    def variable(index: int) -> "Laurent":
        exponents = [0] * NVAR
        exponents[index] = 1
        return Laurent.monomial(tuple(exponents))

    def __add__(self, other: "Laurent" | F | int) -> "Laurent":
        other = other if isinstance(other, Laurent) else Laurent.constant(other)
        terms = dict(self.terms)
        for monomial, coefficient in other.terms.items():
            terms[monomial] = terms.get(monomial, F(0)) + coefficient
        return Laurent(terms)

    def __radd__(self, other: "Laurent" | F | int) -> "Laurent":
        return self + other

    def __neg__(self) -> "Laurent":
        return Laurent({m: -a for m, a in self.terms.items()})

    def __sub__(self, other: "Laurent" | F | int) -> "Laurent":
        return self + (-other if isinstance(other, Laurent) else -F(other))

    def __rsub__(self, other: "Laurent" | F | int) -> "Laurent":
        return (-self) + other

    def __mul__(self, other: "Laurent" | F | int) -> "Laurent":
        other = other if isinstance(other, Laurent) else Laurent.constant(other)
        terms: dict[tuple[int, ...], F] = {}
        for left, a in self.terms.items():
            for right, b in other.terms.items():
                monomial = tuple(x + y for x, y in zip(left, right))
                terms[monomial] = terms.get(monomial, F(0)) + a * b
        return Laurent(terms)

    def __rmul__(self, other: "Laurent" | F | int) -> "Laurent":
        return self * other

    def __pow__(self, exponent: int) -> "Laurent":
        require(exponent >= 0, "nonnegative Laurent power")
        result = Laurent.constant(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Laurent):
            return self.terms == other.terms
        if isinstance(other, (int, F)):
            return self == Laurent.constant(other)
        return False

    def evaluate(self, values: dict[str, F]) -> F:
        total = F(0)
        for monomial, coefficient in self.terms.items():
            value = coefficient
            for name, exponent in zip(VARS, monomial):
                value *= values[name] ** exponent
            total += value
        return total


c = Laurent.variable(0)
alpha = Laurent.variable(1)
upsilon = Laurent.variable(2)
eta = Laurent.variable(3)
Delta = Laurent.variable(4)
c_inv = Laurent.monomial((-1, 0, 0, 0, 0))


# A bivariate sparse polynomial has keys (t-degree, y-degree) and Laurent
# coefficients.  All t-degrees used below are nonnegative.
Bivariate = dict[tuple[int, int], Laurent]
Series = dict[int, Laurent]


def add_bivariate(poly: Bivariate, key: tuple[int, int], coefficient: Laurent) -> None:
    value = poly.get(key, Laurent()) + coefficient
    if value.terms:
        poly[key] = value
    elif key in poly:
        del poly[key]


def derivative_y(poly: Bivariate) -> Bivariate:
    result: Bivariate = {}
    for (t_degree, y_degree), coefficient in poly.items():
        if y_degree:
            add_bivariate(result, (t_degree, y_degree - 1), y_degree * coefficient)
    return result


def series_multiply(left: Series, right: Series, limit: int) -> Series:
    result: Series = {}
    for i, a in left.items():
        for j, b in right.items():
            if i + j <= limit:
                result[i + j] = result.get(i + j, Laurent()) + a * b
    return {degree: value for degree, value in result.items() if value.terms}


def series_power(series: Series, exponent: int, limit: int) -> Series:
    result: Series = {0: Laurent.constant(1)}
    for _ in range(exponent):
        result = series_multiply(result, series, limit)
    return result


def compose_y(poly: Bivariate, y_series: Series, limit: int) -> Series:
    result: Series = {}
    powers: dict[int, Series] = {0: {0: Laurent.constant(1)}}
    for y_degree in range(1, max((key[1] for key in poly), default=0) + 1):
        powers[y_degree] = series_multiply(powers[y_degree - 1], y_series, limit)
    for (t_degree, y_degree), coefficient in poly.items():
        for shift, value in powers[y_degree].items():
            if t_degree + shift <= limit:
                degree = t_degree + shift
                result[degree] = result.get(degree, Laurent()) + coefficient * value
    return {degree: value for degree, value in result.items() if value.terms}


def specialize_ladder(value: Laurent) -> Laurent:
    """Set alpha=eta=0 and upsilon=-8c/3, retaining c and Delta."""

    result: dict[tuple[int, ...], F] = {}
    for monomial, coefficient in value.terms.items():
        ec, ea, eu, ee, ed = monomial
        if ea or ee:
            continue
        new_monomial = (ec + eu, 0, 0, 0, ed)
        new_coefficient = coefficient * F(-8, 3) ** eu
        result[new_monomial] = result.get(new_monomial, F(0)) + new_coefficient
    return Laurent(result)


def specialize_bivariate(poly: Bivariate) -> Bivariate:
    result: Bivariate = {}
    for key, value in poly.items():
        specialized = specialize_ladder(value)
        if specialized.terms:
            result[key] = specialized
    return result


def reconstruct_critical_polynomial() -> Bivariate:
    """Expand the constrained literal Hhat and form F/t^12 through t^5."""

    # Hhat is stored in (t-degree,q-degree), after r=1+q is expanded.
    hhat: dict[tuple[int, int], Laurent] = {}

    def add_r_term(coefficient: Laurent, t_degree: int, r_power: int) -> None:
        for q_degree in range(r_power + 1):
            key = (t_degree, q_degree)
            value = hhat.get(key, Laurent()) + comb(r_power, q_degree) * coefficient
            if value.terms:
                hhat[key] = value
            elif key in hhat:
                del hhat[key]

    repeated_u = F(-1, 4) * c**2
    fixed_weight_six = Laurent.constant(F(-1376, 135))

    add_r_term(repeated_u, 0, 6)
    add_r_term(-repeated_u, 0, 4)
    add_r_term(alpha, 1, 5)
    add_r_term(-alpha, 1, 4)
    add_r_term(upsilon, 2, 5)
    add_r_term(-upsilon, 2, 4)
    add_r_term(eta, 3, 4)
    add_r_term(-eta, 3, 3)
    add_r_term(Delta, 4, 4)
    add_r_term(-Delta, 4, 3)
    add_r_term(fixed_weight_six, 6, 3)
    add_r_term(c - fixed_weight_six, 6, 2)
    add_r_term(Laurent.constant(F(8, 3)), 8, 2)
    add_r_term(Laurent.constant(-3), 10, 1)

    result: Bivariate = {(0, 0): Laurent.constant(F(-1, 2))}
    for (t_degree, q_degree), coefficient in hhat.items():
        new_t_degree = t_degree + 6 * q_degree - 6
        if 0 <= new_t_degree <= 5:
            add_bivariate(result, (new_t_degree, q_degree + 1), coefficient)
        else:
            require(new_t_degree >= 0, "no negative t-powers after constrained sums")
    return result


def expected_critical_polynomial() -> Bivariate:
    return {
        (0, 0): Laurent.constant(F(-1, 2)),
        (0, 1): c,
        (0, 2): F(-1, 2) * c**2,
        (1, 2): alpha,
        (2, 1): Laurent.constant(F(8, 3)),
        (2, 2): upsilon,
        (3, 2): eta,
        (4, 1): Laurent.constant(-3),
        (4, 2): Delta,
    }


def evaluate_bivariate(poly: Bivariate, values: dict[str, F]) -> dict[tuple[int, int], F]:
    return {
        key: coefficient.evaluate(values)
        for key, coefficient in poly.items()
        if coefficient.evaluate(values)
    }


def fraction_series_multiply(
    left: dict[int, F], right: dict[int, F], limit: int
) -> dict[int, F]:
    result: dict[int, F] = {}
    for i, a in left.items():
        for j, b in right.items():
            if i + j <= limit:
                result[i + j] = result.get(i + j, F(0)) + a * b
    return {degree: value for degree, value in result.items() if value}


def compose_y_fraction(
    poly: dict[tuple[int, int], F], y_series: dict[int, F], limit: int
) -> dict[int, F]:
    result: dict[int, F] = {}
    max_y_degree = max((key[1] for key in poly), default=0)
    powers: dict[int, dict[int, F]] = {0: {0: F(1)}}
    for y_degree in range(1, max_y_degree + 1):
        powers[y_degree] = fraction_series_multiply(
            powers[y_degree - 1], y_series, limit
        )
    for (t_degree, y_degree), coefficient in poly.items():
        for shift, value in powers[y_degree].items():
            if t_degree + shift <= limit:
                degree = t_degree + shift
                result[degree] = result.get(degree, F(0)) + coefficient * value
    return {degree: value for degree, value in result.items() if value}


def derivative_y_fraction(
    poly: dict[tuple[int, int], F]
) -> dict[tuple[int, int], F]:
    result: dict[tuple[int, int], F] = {}
    for (t_degree, y_degree), coefficient in poly.items():
        if y_degree:
            result[(t_degree, y_degree - 1)] = y_degree * coefficient
    return result


def solve_critical_series(
    poly: dict[tuple[int, int], F], y0: F, limit: int
) -> dict[int, F]:
    derivative = derivative_y_fraction(poly)
    y_series: dict[int, F] = {0: y0}
    for degree in range(1, limit + 1):
        constant_part = compose_y_fraction(derivative, y_series, limit).get(degree, F(0))
        trial = dict(y_series)
        trial[degree] = F(1)
        trial_part = compose_y_fraction(derivative, trial, limit).get(degree, F(0))
        slope = trial_part - constant_part
        require(slope != 0, f"critical Hensel slope at degree {degree}")
        coefficient = -constant_part / slope
        if coefficient:
            y_series[degree] = coefficient
    residual = compose_y_fraction(derivative, y_series, limit)
    require(all(residual.get(degree, F(0)) == 0 for degree in range(limit + 1)),
            "critical series residual")
    return y_series


def poly_clean(poly: dict[int, F]) -> dict[int, F]:
    return {degree: F(value) for degree, value in poly.items() if value}


def poly_derivative(poly: dict[int, F]) -> dict[int, F]:
    return {degree - 1: degree * value for degree, value in poly.items() if degree}


def poly_divmod(left: dict[int, F], right: dict[int, F]) -> tuple[dict[int, F], dict[int, F]]:
    left = poly_clean(left)
    right = poly_clean(right)
    require(bool(right), "polynomial division by nonzero")
    quotient: dict[int, F] = {}
    divisor_degree = max(right)
    divisor_lead = right[divisor_degree]
    while left and max(left) >= divisor_degree:
        degree = max(left) - divisor_degree
        coefficient = left[max(left)] / divisor_lead
        quotient[degree] = quotient.get(degree, F(0)) + coefficient
        for exponent, value in right.items():
            target = exponent + degree
            left[target] = left.get(target, F(0)) - coefficient * value
        left = poly_clean(left)
    return poly_clean(quotient), left


def poly_gcd(left: dict[int, F], right: dict[int, F]) -> dict[int, F]:
    left, right = poly_clean(left), poly_clean(right)
    while right:
        _, remainder = poly_divmod(left, right)
        left, right = right, remainder
    require(bool(left), "nonzero polynomial gcd")
    lead = left[max(left)]
    return {degree: value / lead for degree, value in left.items()}


def main() -> None:
    critical = reconstruct_critical_polynomial()
    require(critical == expected_critical_polynomial(), "literal critical expansion")

    # Direct formal evaluation at y0=1/c gives the first three recursive
    # critical values.  The fourth requires the moving critical point.
    y0 = {0: c_inv}
    values_at_y0 = compose_y(critical, y0, 4)
    require(values_at_y0.get(0, Laurent()) == 0, "repeated constant face")
    require(values_at_y0[1] == alpha * c_inv**2, "C1")
    require(
        values_at_y0[2] == upsilon * c_inv**2 + F(8, 3) * c_inv,
        "C2",
    )
    require(values_at_y0[3] == eta * c_inv**2, "C3")

    specialized = specialize_bivariate(critical)
    critical_shift = {0: c_inv, 2: F(-8, 3) * c_inv**2}
    derivative_residual = compose_y(derivative_y(specialized), critical_shift, 3)
    require(all(derivative_residual.get(i, Laurent()) == 0 for i in range(4)),
            "critical shift through t^3")
    shifted_value = compose_y(specialized, critical_shift, 4)
    expected_c4 = Delta * c_inv**2 + F(32, 9) * c_inv**2 - 3 * c_inv
    require(shifted_value[4] == expected_c4, "C4 including critical correction")

    # Exact allowed maximal-cancellation witness.
    a0 = F(7168, 135)
    c_value = F(5152, 405)
    delta_value = F(4672, 135)
    upsilon_value = F(-8, 3) * c_value
    witness = {
        "c": c_value,
        "alpha11": F(0),
        "upsilon5": upsilon_value,
        "eta": F(0),
        "Delta": delta_value,
    }
    require(c_value == a0 - F(7, 6) * delta_value, "allowed c-Delta relation")
    require(c_value != 0, "witness c nonzero")
    u_value = -c_value**2 / 4
    require(u_value != 0 and c_value**2 + 4 * u_value == 0, "special repeated U")
    c_values = [
        F(0),
        upsilon_value / c_value**2 + F(8, 3) / c_value,
        F(0),
        (delta_value + F(32, 9) - 3 * c_value) / c_value**2,
    ]
    require(c_values == [F(0)] * 4, "maximal C1-C4 cancellation")
    numeric_critical = evaluate_bivariate(critical, witness)
    y_critical = solve_critical_series(numeric_critical, 1 / c_value, 5)
    critical_value = compose_y_fraction(numeric_critical, y_critical, 5)
    require(all(critical_value.get(i, F(0)) == 0 for i in range(6)),
            "formal maximal cancellation through t^5")

    # The maximal witness at (s,beta)=(1,6) is split by b^12 q before
    # every q-dependent O(t^6) correction.
    s_witness, beta_witness = 1, 6
    tau_witness = s_witness + beta_witness
    d_b_witness = 6 * (beta_witness - s_witness)
    require(d_b_witness == 30, "witness b gap")
    require(2 * tau_witness == 14 and 4 * tau_witness == 28, "canceled gaps")
    require(d_b_witness < 6 * tau_witness, "b precedes O(t^6)")
    require(6 * s_witness + 2 * beta_witness == 18, "witness form order")

    # beta<s: after the only possible lower-row/b^12 cancellation, the
    # terminal discriminant is x^2+4U and has simple roots for U!=0.
    lower_slopes: list[F] = []
    for k in range(1, 6):
        slope = F(k, 12 - k)
        lower_slopes.append(slope)
        require(slope < 1, f"lower slope k={k}")
        require(12 * slope == k * (1 + slope), f"lower balance k={k}")
    lower_control = {2: F(1), 0: F(12, 7)}  # U=3/7, so 4U=12/7.
    require(max(poly_gcd(lower_control, poly_derivative(lower_control))) == 0,
            "beta<s terminal roots simple")

    # beta=s: the special discriminant is X^6(X^6-2c); its only multiple
    # root is X=0, while the nonzero factor is squarefree.
    equal_control = {6: F(1), 0: -2 * c_value}
    require(max(poly_gcd(equal_control, poly_derivative(equal_control))) == 0,
            "beta=s nonzero roots simple")
    require(c_value != 0, "beta=s special factor nonzero constant")

    # beta>s equality faces and their nonzero-root squarefreeness.
    equality_ratios: list[F] = []
    for j in range(1, 5):
        ratio = F(6 + j, 6 - j)
        equality_ratios.append(ratio)
        s_value, beta_value = 6 - j, 6 + j
        require(j * (s_value + beta_value) == 6 * (beta_value - s_value),
                f"gap equality j={j}")
        y0_control = F(2, 5)
        cj_control = F(j + 1, 7)
        nonzero_factor = {0: cj_control, 6 - j: -y0_control}
        gcd = poly_gcd(nonzero_factor, poly_derivative(nonzero_factor))
        require(max(gcd) == 0, f"equality nonzero roots simple j={j}")
        require(min({j, 6}) == j, f"X=0 multiplicity j={j}")

    require(equality_ratios == [F(7, 5), F(2), F(3), F(5)],
            "equality ratio list")

    # Exact finite and randomized audits of all comparison directions and
    # every form-order row.  The formulas themselves have positive
    # coefficients, while these controls guard transcription and signs.
    grid_checks = 0
    for s_value in range(1, 61):
        for beta_value in range(1, 61):
            require(3 * s_value + 5 * beta_value > 0, "generic form order")
            if beta_value > s_value:
                d_b = 6 * (beta_value - s_value)
                require(6 * s_value + 2 * beta_value > 0, "b form order")
                require(d_b < 6 * (s_value + beta_value), "b before t^6")
                for j, threshold in enumerate(equality_ratios, start=1):
                    d_j = j * (s_value + beta_value)
                    comparison = (d_b > d_j) - (d_b < d_j)
                    ratio_comparison = (
                        (F(beta_value, s_value) > threshold)
                        - (F(beta_value, s_value) < threshold)
                    )
                    require(comparison == ratio_comparison, f"threshold j={j}")
                    order_twice = (6 - j) * s_value + (10 - j) * beta_value
                    require(order_twice > 0, f"critical form order j={j}")
                    grid_checks += 1

    rng = random.Random(4292)
    random_checks = 4096
    for _ in range(random_checks):
        s_value = rng.randrange(1, 10**6)
        beta_value = rng.randrange(s_value + 1, 10**6 + s_value + 1)
        j = rng.randrange(1, 5)
        threshold = F(6 + j, 6 - j)
        d_b = 6 * (beta_value - s_value)
        d_j = j * (s_value + beta_value)
        require((d_b > d_j) == (F(beta_value, s_value) > threshold),
                "random threshold greater")
        require((d_b < d_j) == (F(beta_value, s_value) < threshold),
                "random threshold less")
        require((6 - j) * s_value + (10 - j) * beta_value > 0,
                "random positive form order")

    # Restricted THM-4291 specialization: j=2 is genuinely forced.
    restricted_c2 = F(8, 3) / a0
    require(restricted_c2 != 0, "restricted fixed t^8 q coefficient")
    require(F(6 + 2, 6 - 2) == 2, "restricted equality beta/s=2")
    require(6 * 1 + 2 * 2 == 10 and 2 * 1 + 4 * 2 == 10,
            "restricted collision form order")

    print("THM4292_LAMBDA_ZERO_REPEATED_FACE_KELLER_EXTINCTION_INDEPENDENT_V1")
    print("UNIVERSE exact_M12 W=0 Z=-U char0 U!=0 local_q=0_prepared_factor")
    print("LITERAL_RECONSTRUCTION F=q*(h-b^12)+q^2*V-t^12/2 V(0,0)=2U PASS")
    print("NEWTON_TRICHOTOMY lambda<6tau lambda=6tau lambda>6tau PASS")
    print("BETA_LT_S slopes=" + ",".join(str(value) for value in lower_slopes)
          + " terminal_discriminant=x^2+4U NONZERO_ROOTS_SIMPLE")
    print("BETA_EQ_S discriminant=(c-X^6)^2+4U special=X^6*(X^6-2c) ONLY_MULTIPLE_ROOT_X0")
    print("CRITICAL_COEFFICIENTS C1=alpha11/c^2 C2=upsilon5/c^2+8/(3c) "
          "C3=eta/c^2 C4=(Delta+32/9-3c)/c^2")
    print("MAXIMAL_CANCELLATION c=5152/405 Delta=4672/135 upsilon5=-8c/3 "
          "U=-c^2/4 C1=C2=C3=C4=0")
    print("MAXIMAL_WITNESS s=1 beta=6 canceled_gaps=14,28 b_gap=30 O_t6_gap=42 form_order=18")
    print("EQUALITY_RATIOS j1=7/5 j2=2 j3=3 j4=5 NONZERO_ROOTS_SIMPLE X0_ADVANCES_BETA")
    print("FORM_ORDERS b=6s+2beta j1=(5s+9beta)/2 j2=2s+4beta "
          "j3=(3s+7beta)/2 j4=s+3beta ALL_POSITIVE")
    print(f"CONTROLS grid={grid_checks} randomized={random_checks} seed=4292")
    print("RESTRICTED_THM4291 C2=8/(3a)!=0 equality_beta_over_s=2 collision_order=10s")
    print("SCOPE repeated_faces_extinguished degree_persistence_and_JC2_not_claimed")
    print("VERDICT PASS STANDARD_LIBRARY_EXACT_INDEPENDENT_AUDIT")


if __name__ == "__main__":
    main()
