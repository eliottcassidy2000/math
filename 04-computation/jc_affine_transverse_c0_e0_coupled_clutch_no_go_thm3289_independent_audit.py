#!/usr/bin/env python3
"""Independent exact audit of the provisional THM-3289 scout.

This audit deliberately does not import or execute the primary THM-3289
script.  It reconstructs the two cubic response pairs from their coefficient
definitions, asks SymPy's resultant engine for the localized elimination,
derives the response and clutch jets with truncated coefficient lists, and
uses a small independent dense Euclidean algorithm for the two final wall
polynomials.  In particular, no pre-factored fifth- or sixth-jet formula is
used as input.

The audited statement remains a critical-point obstruction for one explicit
first-coordinate family.  It is not an inverse-cover construction and does
not prove the planar Jacobian conjecture.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import comb
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    "04-computation/jc_centered_heptic_source_morse_obstruction_thm3212.py":
        "03d121e57dd2edaece7cd8693d792349f03757c6e781eb5d9d0c897fcc889448",
    "05-knowledge/results/jc_centered_heptic_source_morse_obstruction_thm3212.out":
        "729e0c7b9fa51fa5c4ac5f18f50dc4413c8a8bb7bf5f0ebf2a7709650304bc85",
    "04-computation/jc_heptic_constant_parameter_discriminant_audit_thm3212.py":
        "adeb3c548d5fe3966eefc7ec4eeadfd1410a62356eca8a6c387e39cbe8fc6122",
    "05-knowledge/results/jc_heptic_constant_parameter_discriminant_audit_thm3212.out":
        "d170cf2212848ef76722579a40b65506bedf6e65a031012ca06c27dcd1ef77bb",
    "04-computation/jc_affine_transverse_c0_clutch_no_go_thm3279.py":
        "06820b2476fc0f2cefe3982d054a7db09bb88b4892503580550a1b154564508a",
    "05-knowledge/results/jc_affine_transverse_c0_clutch_no_go_thm3279.out":
        "4a88b5ab31eed4c9a5f90f814a6a24a614db0afecc4cf1cab7fa32dae7c991c4",
    "04-computation/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289.py":
        "f63ff06e3f5ed30f3f6bc5be99756c347d6af5f8e9b220ce8336abff2cd2ca31",
    "05-knowledge/results/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289.out":
        "1aef4341650cdfaf1043a8699e3a1725a0100af6d9848d99dfa924b6f054dba1",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def series_product(first, second, order):
    """Multiply coefficient lists modulo t**(order+1)."""

    zero = first[0] * 0 + second[0] * 0
    result = [zero for _ in range(order + 1)]
    for left_degree, left_coefficient in enumerate(first):
        if left_coefficient == 0:
            continue
        ceiling = order + 1 - left_degree
        for right_degree, right_coefficient in enumerate(second[:ceiling]):
            if right_coefficient != 0:
                result[left_degree + right_degree] += (
                    left_coefficient * right_coefficient
                )
    return result


def series_power(coefficients, exponent, order, one):
    result = [one] + [one * 0 for _ in range(order)]
    for _ in range(exponent):
        result = series_product(result, coefficients, order)
    return result


def evaluate_factor_series(K_poly, base_series, order, one):
    """Evaluate a sparse abstract polynomial in truncated coefficient lists."""

    caches = {
        (base_index, exponent): series_power(series, exponent, order, one)
        for base_index, series in enumerate(base_series)
        for exponent in range(7)
    }
    result = [one * 0 for _ in range(order + 1)]
    for monomial, coefficient in K_poly.terms():
        term = [one * coefficient] + [one * 0 for _ in range(order)]
        for base_index, exponent in enumerate(monomial):
            term = series_product(
                term, caches[(base_index, exponent)], order
            )
        for degree in range(order + 1):
            result[degree] += term[degree]
    return result


def derive_localized_factor():
    """Use the builtin resultant route, not a literal Sylvester determinant."""

    y, A, C, V, derivative_V, d, k = sp.symbols("y A C V Vp d k")
    derivative_A = sp.symbols("Ap")
    L = y**2 + y + C * V
    first = 2 * L * (2 * y + 1) + V * A
    second = V**3 * k + V**2 * y + L * (-derivative_V * y + 2 * V**2 * d)
    require(sp.Poly(first, y).LC() == 4, "R1 has nonzero projective leading term")

    resultant = sp.expand(sp.resultant(first, second, y))
    quotient, remainder = sp.div(
        sp.Poly(resultant, V), sp.Poly(V**3, V)
    )
    require(remainder.is_zero, "builtin resultant has exact V^3 factor")
    K_expression = sp.expand(quotient.as_expr())
    K_poly = sp.Poly(K_expression, A, C, V, derivative_V, d, k)
    require(len(K_poly.terms()) == 40, "localized factor has 40 terms")

    raw_x_equation = (
        2 * L * (derivative_V * y**2 + d * V**2)
        + derivative_A * y * V**2
        + k * V**3
    )
    response_derivative = (2 * V + A * derivative_V) / (2 * V)
    reduction = sp.cancel(
        raw_x_equation.subs(derivative_A, response_derivative)
        - second
        - derivative_V * y * first / 2
    )
    require(reduction == 0, "localized gradient reduction")

    degree_weights = (8, 1, 16, 15, 0, 0)
    weighted_rows = [
        (
            sum(power * weight for power, weight in zip(monomial, degree_weights)),
            monomial,
            coefficient,
        )
        for monomial, coefficient in K_poly.terms()
    ]
    maximum_weight = max(row[0] for row in weighted_rows)
    maximum_rows = [row for row in weighted_rows if row[0] == maximum_weight]
    require(maximum_weight == 97, "weighted degree is 97")
    require(
        maximum_rows == [(97, (0, 1, 6, 0, 1, 2), sp.Integer(128))],
        "unique degree-97 term is 128*C*V^6*d*k^2",
    )
    print(
        "resultant_route=sympy_builtin;resultant=V^3*K40;"
        "gradient_reduction=raw-R2-(Vp*y/2)*R1=0;"
        "R1_y_lead=4=>no_projective_y_infinity"
    )
    print(
        "unique_weighted_degree97=128*C*V^6*d*k^2;"
        "affine_lead=128*d^2*k^2*V_lead^6"
    )
    return K_poly, (A, C, V, derivative_V, d, k)


def finite_gate_audit() -> None:
    """Check the finite owner calculation before any localized elimination."""

    x, z, A, C, V, derivative_A, derivative_V, c0, d, e0, k = sp.symbols(
        "x z A C V Ap Vp c0 d e0 k"
    )
    inner = V * z**2 + z + C
    derivative_z = 2 * inner * (2 * V * z + 1) + A
    derivative_x = (
        2 * inner * (derivative_V * z**2 + d) + derivative_A * z + k
    )
    owner = {V: 0, A: 0, z: -C}
    require(sp.expand(derivative_z.subs(owner)) == 0, "owner z derivative")
    require(
        sp.expand(derivative_x.subs(owner)) == k - derivative_A * C,
        "owner finite clutch is k-A'*C",
    )
    require(sp.diff(e0 + k * x, x) == k, "e0 is gradient-inert")
    print(
        "coefficient_scope=char0_field_K0;c0,d,e0,k_in_K0;"
        "C=c0+d*x;E=e0+k*x;e0_gradient_inert=1"
    )
    print(
        "finite_owner_gate=Delta=k-Aprime*C;"
        "Delta=0_gives_explicit_critical_point_z=-C;"
        "parameter_lanes=d=0|d!=0,k=0|d*k!=0"
    )


def t_owner_audit(K_poly, symbols) -> None:
    A, C, V, derivative_V, d_symbol, k_symbol = symbols
    t, v, c, d, k = sp.symbols("t v c d k", nonzero=True)
    rows = []
    for multiplicity in (3, 4, 5, 6):
        response_slope = sp.Rational(2, 2 - multiplicity)
        local_expression = sp.expand(
            K_poly.as_expr().subs(
                {
                    A: response_slope * t,
                    C: c + d * t,
                    V: v * t**multiplicity,
                    derivative_V: multiplicity * v * t ** (multiplicity - 1),
                    d_symbol: d,
                    k_symbol: k,
                }
            )
        )
        local = sp.Poly(local_expression, t)
        first_possible = 3 * multiplicity - 1
        expected = (
            sp.Rational(
                16 * multiplicity * (multiplicity - 1), multiplicity - 2
            )
            * v**3
            * (k - response_slope * c)
        )
        require(
            all(local.nth(index) == 0 for index in range(first_possible)),
            f"T lower coefficients vanish for multiplicity {multiplicity}",
        )
        require(
            sp.factor(local.nth(first_possible) - expected) == 0,
            f"T leading coefficient for multiplicity {multiplicity}",
        )
        rows.append((multiplicity, first_possible))
    print(
        f"T_rows={tuple(rows)};"
        "lead=16*m*(m-1)/(m-2)*v^3*(k-Aprime*c);"
        "no_T_escape_under_finite_gate"
    )


def response_jet(v_jet, order, zero, two):
    """Solve 2*V*A'-A*V'=2*V recursively at a simple V zero."""

    coefficients = [zero for _ in range(order + 1)]
    coefficients[1] = two
    for degree in range(2, order + 1):
        subtotal = zero
        for v_index in range(1, degree + 1):
            a_index = degree + 1 - v_index
            if 1 <= a_index < degree:
                subtotal += (
                    (2 * a_index - v_index)
                    * v_jet[v_index]
                    * coefficients[a_index]
                )
        coefficients[degree] = (
            2 * v_jet[degree] - subtotal
        ) / ((2 * degree - 1) * v_jet[1])
    return coefficients


def derive_generic_s_walls(K_poly, symbols):
    """Derive q3 and q4 from lists; leave q5 and q6 entirely unevaluated here."""

    A, C, V, derivative_V, d_symbol, k_symbol = symbols
    order = 6
    c, d, k = sp.symbols("c d k")
    v_symbols = [sp.S.Zero] + list(sp.symbols("v1:8"))
    v1, v2, v3 = v_symbols[1:4]
    a_jet = response_jet(v_symbols, order, sp.S.Zero, sp.Integer(2))
    base_series = [
        a_jet,
        [c, d] + [sp.S.Zero for _ in range(order - 1)],
        v_symbols[:order + 1],
        [index * v_symbols[index] for index in range(1, order + 2)],
        [d] + [sp.S.Zero for _ in range(order)],
        [k] + [sp.S.Zero for _ in range(order)],
    ]
    q = evaluate_factor_series(K_poly, base_series, order, sp.S.One)
    require(
        all(sp.factor(q[index]) == 0 for index in range(3)),
        "universal simple-S order is at least three",
    )

    expected_q3 = (
        sp.Rational(8, 3)
        * v1**2
        * (2 * c - k)
        * (6 * c * v1**2 + 3 * k * v1**2 + 4 * v2)
    )
    require(sp.factor(q[3] - expected_q3) == 0, "q3 factorization")
    q3_as_k = sp.Poly(q[3], k)
    quotient, remainder = q3_as_k.div(sp.Poly(k - 2 * c, k))
    require(remainder.is_zero and quotient.degree() == 1, "remove finite q3 wall")
    live_k = sp.factor(-quotient.nth(0) / quotient.nth(1))
    require(
        sp.factor(live_k + 2 * c + sp.Rational(4, 3) * v2 / v1**2) == 0,
        "derived live k wall",
    )

    q4_live = sp.factor(q[4].subs(k, live_k))
    q4_as_d = sp.Poly(q4_live, d)
    require(q4_as_d.degree() == 1, "live q4 is linear in d")
    live_d = sp.factor(-q4_as_d.nth(0) / q4_as_d.nth(1))
    expected_live_d = (
        15 * v1**3
        - 18 * v1 * v3
        + 16 * v2**2
        - 30 * c**2 * v1**4
        - 10 * c * v1**2 * v2
    ) / (30 * v1**3)
    require(sp.factor(live_d - expected_live_d) == 0, "derived live d wall")
    slope = 3 * c * v1**2 + v2
    require(
        sp.factor(sp.diff(q4_live, d) - sp.Rational(128, 3) * v1**2 * slope)
        == 0,
        "q4 d-slope",
    )
    exceptional_c = -v2 / (3 * v1**2)
    require(
        sp.factor(live_k.subs(c, exceptional_c) - 2 * exceptional_c) == 0,
        "q4 slope zero returns to finite failure k=2c",
    )
    require(
        all(sp.factor(q[index].subs({k: live_k, d: live_d})) == 0
            for index in range(5)),
        "derived walls kill q0 through q4",
    )
    print(
        "S_q3=(8/3)*v1^2*(2c-k)*(6c*v1^2+3k*v1^2+4v2);"
        "finite_hostile=k=2c;live_k=-2c-4v2/(3v1^2)"
    )
    print(
        "S_q4_d_slope=(128/3)*v1^2*(3c*v1^2+v2);"
        "slope_zero=>live_k=2c=>finite_gate_failure;live_d=derived_unique"
    )
    return {
        "c": c,
        "v_symbols": v_symbols,
        "live_k": live_k,
        "live_d": live_d,
        "slope": slope,
    }


def taylor_coefficients(poly, root, field, count):
    """Expand at x=root directly with the binomial theorem."""

    coefficients = [field.zero for _ in range(count)]
    for (power,), coefficient in poly.to_dict().items():
        ceiling = min(power, count - 1)
        for degree in range(ceiling + 1):
            coefficients[degree] += (
                coefficient * comb(power, degree) * root ** (power - degree)
            )
    return coefficients


def evaluate_positive_expression(expression, values, ring):
    """Evaluate an integer-power SymPy expression in a polynomial ring."""

    if expression in values:
        return values[expression]
    if expression.is_Integer or expression.is_Rational:
        return ring.convert(expression)
    if expression.is_Add:
        total = ring.zero
        for argument in expression.args:
            total += evaluate_positive_expression(argument, values, ring)
        return total
    if expression.is_Mul:
        product = ring.one
        for argument in expression.args:
            product *= evaluate_positive_expression(argument, values, ring)
        return product
    if expression.is_Pow:
        require(expression.exp.is_Integer and expression.exp >= 0, "positive power")
        return evaluate_positive_expression(expression.base, values, ring) ** int(
            expression.exp
        )
    raise RuntimeError(f"unsupported expression node: {expression}")


def evaluate_rational_expression(expression, values, ring):
    numerator, denominator = sp.fraction(sp.cancel(expression))
    numerator_value = evaluate_positive_expression(numerator, values, ring)
    denominator_value = evaluate_positive_expression(denominator, values, ring)
    require(denominator_value.degree() == 0, "wall denominator is field-constant")
    inverse = ring.domain.one / denominator_value.LC
    return numerator_value * (ring.one * inverse)


class CubicArithmetic:
    """Tiny rational-triple implementation of one exact cubic field."""

    def __init__(self, field):
        modulus = field.mod.to_sympy_dict()
        leading = self._fraction(modulus[(3,)])
        self.reduction = tuple(
            self._fraction(modulus[(degree,)]) / leading for degree in range(3)
        )
        self.zero = (Fraction(0), Fraction(0), Fraction(0))
        self.one = (Fraction(1), Fraction(0), Fraction(0))

    @staticmethod
    def _fraction(value):
        return Fraction(int(value.numerator), int(value.denominator))

    def convert(self, value):
        table = value.to_sympy_dict()
        return tuple(self._fraction(table.get((degree,), 0)) for degree in range(3))

    def add(self, first, second):
        return tuple(left + right for left, right in zip(first, second))

    def subtract(self, first, second):
        return tuple(left - right for left, right in zip(first, second))

    def multiply(self, first, second):
        raw = [Fraction(0) for _ in range(5)]
        for left_degree, left in enumerate(first):
            for right_degree, right in enumerate(second):
                raw[left_degree + right_degree] += left * right
        for degree in (4, 3):
            coefficient = raw[degree]
            if coefficient:
                offset = degree - 3
                for lower_degree, modulus_coefficient in enumerate(self.reduction):
                    raw[offset + lower_degree] -= coefficient * modulus_coefficient
        return tuple(raw[:3])

    @lru_cache(maxsize=None)
    def inverse(self, value):
        require(value != self.zero, "cubic-field division by zero")
        basis = (
            self.one,
            (Fraction(0), Fraction(1), Fraction(0)),
            (Fraction(0), Fraction(0), Fraction(1)),
        )
        columns = [self.multiply(value, basis_vector) for basis_vector in basis]
        matrix = [
            [columns[column][row] for column in range(3)]
            + [Fraction(1 if row == 0 else 0)]
            for row in range(3)
        ]
        for column in range(3):
            pivot = next(
                (row for row in range(column, 3) if matrix[row][column]), None
            )
            require(pivot is not None, "singular cubic multiplication matrix")
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            pivot_value = matrix[column][column]
            matrix[column] = [entry / pivot_value for entry in matrix[column]]
            for row in range(3):
                if row == column:
                    continue
                factor = matrix[row][column]
                if factor:
                    matrix[row] = [
                        entry - factor * pivot_entry
                        for entry, pivot_entry in zip(matrix[row], matrix[column])
                    ]
        return tuple(matrix[row][3] for row in range(3))

    def divide(self, numerator, denominator):
        return self.multiply(numerator, self.inverse(denominator))


def dense_trim(coefficients, arithmetic):
    result = list(coefficients)
    while result and result[-1] == arithmetic.zero:
        result.pop()
    return result


def dense_from_poly(poly, arithmetic):
    if not poly:
        return []
    table = poly.to_dict()
    return [
        arithmetic.convert(table.get((degree,), poly.ring.domain.zero))
        for degree in range(poly.degree() + 1)
    ]


def dense_monic(coefficients, arithmetic):
    coefficients = dense_trim(coefficients, arithmetic)
    require(bool(coefficients), "cannot normalize zero polynomial")
    inverse = arithmetic.inverse(coefficients[-1])
    return [arithmetic.multiply(coefficient, inverse) for coefficient in coefficients]


def dense_divmod(dividend, divisor, arithmetic):
    """Independent low-to-high dense long division over an exact field."""

    remainder = dense_trim(dividend, arithmetic)
    divisor = dense_trim(divisor, arithmetic)
    require(bool(divisor), "dense division by zero")
    if len(remainder) < len(divisor):
        return [], remainder
    quotient = [
        arithmetic.zero for _ in range(len(remainder) - len(divisor) + 1)
    ]
    while remainder and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        factor = arithmetic.divide(remainder[-1], divisor[-1])
        quotient[shift] = arithmetic.add(quotient[shift], factor)
        for index, coefficient in enumerate(divisor):
            remainder[index + shift] = arithmetic.subtract(
                remainder[index + shift], arithmetic.multiply(factor, coefficient)
            )
        remainder = dense_trim(remainder, arithmetic)
    return dense_trim(quotient, arithmetic), remainder


def dense_prs(first, second, arithmetic):
    first = dense_trim(first, arithmetic)
    second = dense_trim(second, arithmetic)
    if len(first) < len(second):
        first, second = second, first
    require(bool(second), "PRS needs a nonzero second polynomial")
    profile = [(len(first) - 1, len(second) - 1)]
    while second:
        _, remainder = dense_divmod(first, second, arithmetic)
        if not remainder:
            return tuple(profile), dense_monic(second, arithmetic)
        first, second = second, dense_monic(remainder, arithmetic)
        profile.append((len(first) - 1, len(second) - 1))
    raise RuntimeError("unreachable PRS state")


def dense_exact_quotient(dividend, divisor, arithmetic):
    quotient, remainder = dense_divmod(dividend, divisor, arithmetic)
    require(not remainder, "dense exact division has a remainder")
    return quotient


def serialize_dense(coefficients) -> str:
    return "|".join(f"{index}:{coefficient}" for index, coefficient in enumerate(coefficients))


def pair_digest(first, second, arithmetic) -> str:
    payload = (
        serialize_dense(dense_monic(first, arithmetic))
        + "\n"
        + serialize_dense(dense_monic(second, arithmetic))
    )
    return sha256(payload.encode("utf-8")).hexdigest()


def specialize_factor(K_poly, base_values, polynomial_ring):
    result = polynomial_ring.zero
    for monomial, coefficient in K_poly.terms():
        term = polynomial_ring.convert(coefficient)
        for value, exponent in zip(base_values, monomial):
            term *= value**exponent
        result += term
    return result


def canonical_poly_text(poly) -> str:
    return "|".join(
        f"{degree}:{coefficient}"
        for (degree,), coefficient in sorted(poly.monic().to_dict().items())
    )


def build_response_case(name: str):
    """Rebuild a THM-3212 response pair from the stated cubic data."""

    u, x = sp.symbols("u x")
    if name == "4111":
        modulus = sp.Poly(100 * u**3 + 244 * u**2 + 237 * u + 44, u, domain=QQ)
        exponent_a, exponent_b = 4, 1
    elif name == "3211":
        modulus = sp.Poly(75 * u**3 - 89 * u**2 - 31 * u + 61, u, domain=QQ)
        exponent_a, exponent_b = 3, 2
    else:
        raise RuntimeError(f"unknown response case: {name}")

    field = QQ.alg_field_from_poly(modulus, alias="u")
    alpha = field.ext
    x_ring = field.poly_ring(x)
    X = x_ring.gens[0]
    if name == "4111":
        accessory_v = (8 * alpha**2 + 9 * alpha + 8) / 7
        shift = 5 * (alpha + 1) / 7
        A0 = 80 * accessory_v**2 * (alpha + 1) / 343
    else:
        accessory_v = (24 * alpha**2 - 16 * alpha - 16) / 21
        shift = (5 * alpha - 4) / 7
        A0 = 9 * accessory_v**2 * (5 * alpha - 4) / 343

    gamma = -7 * A0
    quadratic = X**2 - alpha * X + accessory_v
    D = X**exponent_a * (X - 1) ** exponent_b * quadratic
    T = X * (X - 1) * quadratic
    E = (
        exponent_a * (X - 1) * quadratic
        + exponent_b * X * quadratic
        + X * (X - 1) * (2 * X - alpha)
    ) / 7
    S = X + shift
    V = 4 * S * D * T**2 / gamma**2
    A = 2 * S * E * T / gamma
    g = S * T
    require((V.degree(), A.degree()) == (16, 8), f"{name} response degrees")
    require(
        2 * V * A.diff(X) - A * V.diff(X) == 2 * V,
        f"{name} response identity",
    )
    radical = V.exquo(V.gcd(V.diff(X))).monic()
    require(radical == g.monic(), f"{name} radical(V)=S*T")
    return {
        "name": name,
        "field": field,
        "x_ring": x_ring,
        "X": X,
        "a": exponent_a,
        "b": exponent_b,
        "shift": field.convert(shift),
        "S": S,
        "T": T,
        "V": V,
        "A": A,
        "g": g,
    }


def audit_case(case, K_poly, wall_data):
    name = case["name"]
    field = case["field"]
    x_ring = case["x_ring"]
    X = case["X"]
    V = case["V"]
    A = case["A"]
    S = case["S"]
    T = case["T"]
    g = case["g"]

    v_jet = taylor_coefficients(V, -case["shift"], field, 8)
    a_global_jet = taylor_coefficients(A, -case["shift"], field, 7)
    require(v_jet[0] == field.zero and v_jet[1] != field.zero, f"{name} simple S")
    a_response_jet = response_jet(v_jet, 6, field.zero, field.convert(2))
    require(
        a_response_jet == a_global_jet,
        f"{name} recursively derived A jet matches rebuilt response",
    )

    c_symbol = wall_data["c"]
    c_ring = field.poly_ring(c_symbol)
    local_c = c_ring.gens[0]
    value_map = {c_symbol: local_c}
    for index in range(1, 8):
        value_map[wall_data["v_symbols"][index]] = c_ring.one * v_jet[index]
    live_k = evaluate_rational_expression(wall_data["live_k"], value_map, c_ring)
    live_d = evaluate_rational_expression(wall_data["live_d"], value_map, c_ring)
    slope = evaluate_rational_expression(wall_data["slope"], value_map, c_ring)
    order = 6
    base_series = [
        [c_ring.one * value for value in a_response_jet],
        [local_c, live_d] + [c_ring.zero for _ in range(order - 1)],
        [c_ring.one * value for value in v_jet[:order + 1]],
        [c_ring.one * (index * v_jet[index]) for index in range(1, order + 2)],
        [live_d] + [c_ring.zero for _ in range(order)],
        [live_k] + [c_ring.zero for _ in range(order)],
    ]
    local_q = evaluate_factor_series(K_poly, base_series, order, c_ring.one)
    require(all(not local_q[index] for index in range(5)), f"{name} q0..q4 walls")
    q5 = local_q[5]
    q6 = local_q[6]
    require((q5.degree(), q6.degree(), slope.degree()) == (3, 4, 1), f"{name} jet degrees")

    arithmetic = CubicArithmetic(field)
    q5_dense = dense_from_poly(q5, arithmetic)
    q6_dense = dense_from_poly(q6, arithmetic)
    slope_dense = dense_from_poly(slope, arithmetic)
    q5_core_dense = dense_exact_quotient(q5_dense, slope_dense, arithmetic)
    require(len(q5_core_dense) - 1 == 2, f"{name} q5 slope quotient degree")
    full_profile, full_gcd = dense_prs(q6_dense, q5_dense, arithmetic)
    live_profile, live_gcd = dense_prs(q6_dense, q5_core_dense, arithmetic)
    slope_profile, slope_gcd = dense_prs(q6_dense, slope_dense, arithmetic)
    require(
        full_profile == ((4, 3), (3, 2), (2, 1), (1, 0)),
        f"{name} full independent PRS profile",
    )
    require(
        live_profile == ((4, 2), (2, 1), (1, 0)),
        f"{name} live independent PRS profile",
    )
    require(
        len(full_gcd) == len(live_gcd) == len(slope_gcd) == 1,
        f"{name} independent PRS terminal units",
    )
    wall_digest = pair_digest(q5_dense, q6_dense, arithmetic)

    control_C = X + 1
    control_d = field.one
    control_k = field.one
    control_K = specialize_factor(
        K_poly,
        (A, control_C, V, V.diff(X), control_d, control_k),
        x_ring,
    )
    finite_gate = control_k - A.diff(X) * control_C
    require(finite_gate.gcd(g).degree() == 0, f"{name} control finite gate")
    require(control_K.degree() == 97, f"{name} control K degree")
    require(control_K.LC == 128 * V.LC**6, f"{name} control leading coefficient")
    boundary = (
        S**3
        * T**8
        * X ** (3 * case["a"] - 3)
        * (X - 1) ** (3 * case["b"] - 3)
    )
    quadratic = T.exquo(X * (X - 1))
    owner_factors = (S, X, X - 1, quadratic)
    require(
        all(
            owner_factors[left].gcd(owner_factors[right]).degree() == 0
            for left in range(len(owner_factors))
            for right in range(left + 1, len(owner_factors))
        ),
        f"{name} four owner factors are pairwise disjoint",
    )
    require(
        3 * (case["a"] + 2) - 1 == 8 + 3 * case["a"] - 3
        and 3 * (case["b"] + 2) - 1 == 8 + 3 * case["b"] - 3,
        f"{name} T local orders produce boundary exponents",
    )
    require(boundary.degree() == 44, f"{name} boundary degree")
    residual = control_K.exquo(boundary)
    require(residual.degree() == 53, f"{name} residual degree")
    require(residual.gcd(g).degree() == 0, f"{name} residual owner-disjoint")
    require(
        residual.gcd(residual.diff(X)).degree() == 0,
        f"{name} residual squarefree",
    )
    residual_digest = sha256(canonical_poly_text(residual).encode("utf-8")).hexdigest()
    boundary_shape = (
        f"S^3*T^8*x^{3 * case['a'] - 3}*(x-1)^{3 * case['b'] - 3}"
    )
    print(
        f"case={name};fresh_jets=(deg_q5,deg_q6)=(3,4);"
        f"PRS_full={full_profile};PRS_live={live_profile};"
        f"PRS_slope={slope_profile};terminal_gcd_degrees=(0,0,0);"
        "unit_invariant_wall_ideal=(q5,q6)=1;profiles_are_route_specific;"
        "MISTAKE360_unit_rule_respected=1;"
        f"wall_digest={wall_digest}"
    )
    print(
        f"case={name};control=(C,k)=(1+x,1);finite_gate=unit_mod_ST;"
        f"rad(V)=S*T;boundary_from_disjoint_local_orders={boundary_shape};"
        "degrees=(K,boundary,H)=(97,44,53);"
        f"H_owner_disjoint_squarefree=1;H_digest={residual_digest}"
    )
    return wall_digest, residual_digest


def source_audit() -> None:
    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    assert_count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    float_count = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    require(assert_count == 0 and float_count == 0, "source AST gate")
    imported_roots = {
        alias.name.split(".")[0]
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    }
    imported_roots.update(
        node.module.split(".")[0]
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module is not None
    )
    require(
        not imported_roots.intersection({"importlib", "runpy", "subprocess"}),
        "primary scout cannot be imported or executed indirectly",
    )
    primary_name = "jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm" + "3289.py"
    require(
        source.count(primary_name) == 1,
        "primary filename occurs only in the hash pin",
    )
    print(
        f"source_ast=(assert_nodes={assert_count},float_literals={float_count},"
        "primary_import_or_execution_paths=0)"
    )


def main() -> None:
    print("THM-3289 independent affine C0-E0 coupled-clutch audit")
    for relative_path, expected_hash in DEPENDENCIES.items():
        require(
            lf_hash(ROOT / relative_path) == expected_hash,
            f"dependency drift: {relative_path}",
        )
    print(f"dependency_hash_checks={len(DEPENDENCIES)}")
    K_poly, symbols = derive_localized_factor()
    finite_gate_audit()
    t_owner_audit(K_poly, symbols)
    wall_data = derive_generic_s_walls(K_poly, symbols)
    case_digests = tuple(
        audit_case(build_response_case(name), K_poly, wall_data)
        for name in ("4111", "3211")
    )
    print(f"independent_case_digests={case_digests}")
    print(
        "audited_consequence=finite_gate_and_d*k!=0=>"
        "ord_S(K)<=6,ord_S(H)<=3,no_T_escape,"
        "off_owner_resultant_root_multiplicity>=53-3=50;not_distinct_point_count"
    )
    print(
        "scope=independent_audit_of_provisional_THM3289_and_inherited_dependencies;"
        "lanes_d0_THM3212|dne0_k0_THM3279|dkne0_THM3289;"
        "no_theorem_or_promotion_pin;not_inverse_cover_not_JC2"
    )
    source_audit()


if __name__ == "__main__":
    main()
