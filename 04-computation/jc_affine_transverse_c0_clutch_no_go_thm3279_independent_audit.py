#!/usr/bin/env python3
"""Independent hostile audit of provisional THM-3279.

This companion does not import THM-3279's formulas or algebraic-number-field
objects.  It derives the universal eliminant as a 6-by-6 Sylvester
determinant, solves the response jet recursively, and evaluates the resulting
factor in a custom exact implementation of Q[u]/(q).  The custom quotient
arithmetic independently reproduces both characteristic-zero residual
digests used by the candidate.
"""

import ast
from fractions import Fraction
from functools import reduce
from hashlib import sha256
from pathlib import Path

import sympy as sp


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
PRIMARY_SCRIPT = ROOT / "04-computation/jc_affine_transverse_c0_clutch_no_go_thm3279.py"
PRIMARY_OUTPUT = ROOT / "05-knowledge/results/jc_affine_transverse_c0_clutch_no_go_thm3279.out"
DEPENDENCIES = {
    PRIMARY_SCRIPT:
        "06820b2476fc0f2cefe3982d054a7db09bb88b4892503580550a1b154564508a",
    PRIMARY_OUTPUT:
        "4a88b5ab31eed4c9a5f90f814a6a24a614db0afecc4cf1cab7fa32dae7c991c4",
    ROOT / "04-computation/jc_centered_heptic_source_morse_obstruction_thm3212.py":
        "03d121e57dd2edaece7cd8693d792349f03757c6e781eb5d9d0c897fcc889448",
    ROOT / "05-knowledge/results/jc_centered_heptic_source_morse_obstruction_thm3212.out":
        "729e0c7b9fa51fa5c4ac5f18f50dc4413c8a8bb7bf5f0ebf2a7709650304bc85",
    ROOT / "04-computation/jc_heptic_constant_parameter_discriminant_audit_thm3212.py":
        "adeb3c548d5fe3966eefc7ec4eeadfd1410a62356eca8a6c387e39cbe8fc6122",
    ROOT / "05-knowledge/results/jc_heptic_constant_parameter_discriminant_audit_thm3212.out":
        "d170cf2212848ef76722579a40b65506bedf6e65a031012ca06c27dcd1ef77bb",
}
EXPECTED_DIGESTS = {
    "4111": "c07d66a389368e57ff32cdcc1c10f134951a0f2008ed65b707033d8ea8844e8e",
    "3211": "41228f21741350a603d05efd59a11ec6157d5d0fb17b6b77f03f360cd81e78ad",
}


def lf_hash(path):
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(data).hexdigest()


for dependency, expected in DEPENDENCIES.items():
    require(lf_hash(dependency) == expected, ("dependency drift", dependency.name))


def source_ast(path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return (
        sum(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        sum(isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)),
    )


require(source_ast(PRIMARY_SCRIPT) == (0, 0), "primary AST gate")
require(source_ast(Path(__file__)) == (0, 0), "audit AST gate")


def sylvester_determinant(first, second, variable):
    """Resultant by a literal Sylvester determinant, not sp.resultant."""
    first_poly = sp.Poly(first, variable)
    second_poly = sp.Poly(second, variable)
    first_degree = first_poly.degree()
    second_degree = second_poly.degree()
    first_coefficients = first_poly.all_coeffs()
    second_coefficients = second_poly.all_coeffs()
    rows = []
    for shift in range(second_degree):
        rows.append(
            [0] * shift + first_coefficients
            + [0] * (second_degree - 1 - shift)
        )
    for shift in range(first_degree):
        rows.append(
            [0] * shift + second_coefficients
            + [0] * (first_degree - 1 - shift)
        )
    require(len(rows) == len(rows[0]) == first_degree + second_degree,
            "Sylvester shape")
    return sp.expand(sp.Matrix(rows).det(method="domain-ge"))


def derive_universal_factor():
    y = sp.symbols("y")
    A, B, C, J, K, V, e, Vp, Cp = sp.symbols(
        "A B C J K V e Vp Cp"
    )
    L = y**2 + B * y + C * V
    R1 = 2 * L * (2 * y + B) + V * A
    R2 = V**3 * e + V**2 * y + L * (J * y + K)
    resultant = sylvester_determinant(R1, R2, y)
    universal_factor = sp.cancel(resultant / V**2)
    require(sp.expand(resultant - V**2 * universal_factor) == 0,
            "universal V^2 factor")
    require(len(sp.Poly(universal_factor, A, B, C, J, K, V, e).terms()) == 40,
            "universal factor term census")
    specialized = sp.expand(universal_factor.subs({
        B: 1, J: -Vp, K: 2 * V**2 * Cp, e: 0,
    }))
    critical_factor = sp.cancel(specialized / V)
    require(sp.expand(specialized - V * critical_factor) == 0,
            "specialized extra V factor")
    critical_poly = sp.Poly(critical_factor, A, C, V, Vp, Cp)
    require(len(critical_poly.terms()) == 20,
            "specialized critical term census")
    weights = (8, 1, 16, 15, 0)
    weighted_terms = tuple(
        (sum(exponent * weight for exponent, weight in zip(powers, weights)),
         powers, coefficient)
        for powers, coefficient in critical_poly.terms()
    )
    top_degree = max(record[0] for record in weighted_terms)
    top_terms = tuple(record for record in weighted_terms
                      if record[0] == top_degree)
    require(top_degree == 96 and top_terms == ((96, (2, 0, 5, 0, 3), 32),),
            "affine infinity ledger")
    return critical_factor, (A, C, V, Vp, Cp)


CRITICAL_FACTOR, CRITICAL_SYMBOLS = derive_universal_factor()


def localization_and_jet_audit():
    y, z = sp.symbols("y z")
    A, Ap, C, Cp, V, Vp, ep = sp.symbols("A Ap C Cp V Vp ep")
    H = V * z**2 + z + C
    pz = 2 * H * (2 * V * z + 1) + A
    px = 2 * H * (Vp * z**2 + Cp) + Ap * z + ep
    R1 = sp.cancel(V * pz.subs(z, y / V))
    raw = sp.cancel(V**3 * px.subs(z, y / V))
    L = y**2 + y + C * V
    expected_R1 = 2 * L * (2 * y + 1) + V * A
    R2 = V**3 * ep + V**2 * y + L * (-Vp * y + 2 * V**2 * Cp)
    response_Ap = (2 * V + A * Vp) / (2 * V)
    require(sp.cancel(R1 - expected_R1) == 0, "localized z gradient")
    require(sp.cancel(raw.subs(Ap, response_Ap) - R2
                      - Vp * y * R1 / 2) == 0,
            "localized x-gradient sign")

    t = sp.symbols("t")
    c, d, v = sp.symbols("c d v", nonzero=True)
    leading_rows = []
    A_symbol, C_symbol, V_symbol, Vp_symbol, Cp_symbol = CRITICAL_SYMBOLS
    for multiplicity in (3, 4, 5, 6):
        local_V = v * t**multiplicity
        local_A = sp.Rational(2, 2 - multiplicity) * t
        local_C = c + d * t
        local_K = sp.Poly(sp.expand(CRITICAL_FACTOR.subs({
            A_symbol: local_A,
            C_symbol: local_C,
            V_symbol: local_V,
            Vp_symbol: sp.diff(local_V, t),
            Cp_symbol: d,
        })), t)
        target = 3 * multiplicity - 1
        expected = (sp.Rational(32 * multiplicity * (multiplicity - 1),
                                (multiplicity - 2)**2) * c * v**3)
        require(all(local_K.nth(index) == 0 for index in range(target)),
                ("T order", multiplicity))
        require(local_K.nth(target) == expected,
                ("T coefficient", multiplicity))
        leading_rows.append((multiplicity, target))

    v_symbols = sp.symbols("v1:5")
    a_symbols = sp.symbols("a1:5")
    local_V = sum(v_symbols[index - 1] * t**index for index in range(1, 5))
    local_A = sum(a_symbols[index - 1] * t**index for index in range(1, 5))
    response = sp.expand(
        2 * local_V * sp.diff(local_A, t)
        - local_A * sp.diff(local_V, t) - 2 * local_V
    )
    response_solution = {}
    for order, coefficient_symbol in enumerate(a_symbols, start=1):
        equation = sp.expand(response.subs(response_solution)).coeff(t, order)
        require(equation.has(coefficient_symbol),
                ("response jet not triangular", order))
        solutions = sp.solve(equation, coefficient_symbol)
        require(len(solutions) == 1, ("response jet wall", order))
        response_solution[coefficient_symbol] = solutions[0]
    local_A = sp.expand(local_A.subs(response_solution))
    local_C = c + d * t
    local_K = sp.expand(CRITICAL_FACTOR.subs({
        A_symbol: local_A,
        C_symbol: local_C,
        V_symbol: local_V,
        Vp_symbol: sp.diff(local_V, t),
        Cp_symbol: d,
    }))
    require(all(sp.expand(local_K).coeff(t, index) == 0
                for index in range(3)),
            "universal S^3 boundary")
    v1, v2, v3, v4 = v_symbols[:4]
    q3 = sp.factor(local_K.coeff(t, 3))
    expected_q3 = sp.Rational(32, 3) * c * v1**2 * (3 * c * v1**2 + 2 * v2)
    require(q3 == expected_q3, "independent first S wall")
    wall_c_solutions = sp.solve(q3 / c, c)
    require(wall_c_solutions == [-2 * v2 / (3 * v1**2)],
            "first S wall uniqueness")
    wall_c = wall_c_solutions[0]
    q4 = sp.factor(local_K.coeff(t, 4).subs(c, wall_c))
    expected_q4 = (-sp.Rational(64, 135) * v2 / v1
                   * (90 * d * v1**3 - 45 * v1**3
                      + 54 * v1 * v3 - 28 * v2**2))
    require(q4 == expected_q4, "independent second S wall")
    wall_d_solutions = sp.solve(q4, d)
    expected_wall_d = ((45 * v1**3 - 54 * v1 * v3 + 28 * v2**2)
                       / (90 * v1**3))
    require(wall_d_solutions == [expected_wall_d],
            "second S wall uniqueness")
    wall_d = wall_d_solutions[0]
    obstruction = (315 * v1**3 * v2 + 360 * v1**2 * v4
                   - 414 * v1 * v2 * v3 + 148 * v2**3)
    q5 = sp.factor(local_K.coeff(t, 5).subs({c: wall_c, d: wall_d}))
    require(q5 == -32 * v2 * obstruction / (315 * v1**2),
            "independent untunable fifth jet")
    return tuple(leading_rows)


T_ROWS = localization_and_jet_audit()


class CubicField:
    """Three-coordinate exact arithmetic in Q[u]/(q_3 u^3+...+q_0)."""

    def __init__(self, coefficients):
        self.q = tuple(Fraction(value) for value in coefficients)
        require(len(self.q) == 4 and self.q[3] != 0,
                "cubic modulus")
        self.relation = tuple(-self.q[index] / self.q[3]
                              for index in range(3))
        self.zero = (Fraction(0), Fraction(0), Fraction(0))
        self.one = (Fraction(1), Fraction(0), Fraction(0))
        self.gen = (Fraction(0), Fraction(1), Fraction(0))

    def scalar(self, value):
        return (Fraction(value), Fraction(0), Fraction(0))

    def add(self, left, right):
        return tuple(a + b for a, b in zip(left, right))

    def neg(self, value):
        return tuple(-entry for entry in value)

    def sub(self, left, right):
        return self.add(left, self.neg(right))

    def mul(self, left, right):
        work = [Fraction(0)] * 5
        for left_degree, left_value in enumerate(left):
            for right_degree, right_value in enumerate(right):
                work[left_degree + right_degree] += left_value * right_value
        for degree in (4, 3):
            coefficient = work[degree]
            work[degree] = Fraction(0)
            for offset, replacement in enumerate(self.relation):
                work[degree - 3 + offset] += coefficient * replacement
        return tuple(work[:3])

    def pow(self, value, exponent):
        require(exponent >= 0, "negative field power")
        answer = self.one
        base = value
        remaining = exponent
        while remaining:
            if remaining % 2:
                answer = self.mul(answer, base)
            base = self.mul(base, base)
            remaining //= 2
        return answer

    def inv(self, value):
        require(value != self.zero, "field division by zero")
        basis = (self.one, self.gen, self.mul(self.gen, self.gen))
        matrix = [
            [self.mul(value, basis[column])[row] for column in range(3)]
            + [Fraction(1 if row == 0 else 0)]
            for row in range(3)
        ]
        for column in range(3):
            pivot_row = next((row for row in range(column, 3)
                              if matrix[row][column]), None)
            require(pivot_row is not None, "nonfield cubic quotient")
            matrix[column], matrix[pivot_row] = matrix[pivot_row], matrix[column]
            pivot = matrix[column][column]
            matrix[column] = [entry / pivot for entry in matrix[column]]
            for row in range(3):
                if row == column:
                    continue
                multiplier = matrix[row][column]
                matrix[row] = [entry - multiplier * pivot_entry
                               for entry, pivot_entry
                               in zip(matrix[row], matrix[column])]
        inverse = tuple(matrix[row][3] for row in range(3))
        require(self.mul(value, inverse) == self.one,
                "cubic inverse audit")
        return inverse

    def div(self, left, right):
        return self.mul(left, self.inv(right))


class PolynomialRing:
    """Dense exact univariate polynomials over CubicField."""

    def __init__(self, field):
        self.field = field
        self.zero = ()
        self.one = (field.one,)
        self.x = (field.zero, field.one)

    def trim(self, coefficients):
        answer = list(coefficients)
        while answer and answer[-1] == self.field.zero:
            answer.pop()
        return tuple(answer)

    def const(self, value):
        return () if value == self.field.zero else (value,)

    def add(self, left, right):
        size = max(len(left), len(right))
        return self.trim(
            self.field.add(
                left[index] if index < len(left) else self.field.zero,
                right[index] if index < len(right) else self.field.zero,
            )
            for index in range(size)
        )

    def neg(self, value):
        return tuple(self.field.neg(coefficient) for coefficient in value)

    def sub(self, left, right):
        return self.add(left, self.neg(right))

    def scale(self, value, scalar):
        if scalar == self.field.zero:
            return ()
        return self.trim(self.field.mul(coefficient, scalar)
                         for coefficient in value)

    def mul(self, left, right):
        if not left or not right:
            return ()
        answer = [self.field.zero] * (len(left) + len(right) - 1)
        for left_degree, left_value in enumerate(left):
            for right_degree, right_value in enumerate(right):
                answer[left_degree + right_degree] = self.field.add(
                    answer[left_degree + right_degree],
                    self.field.mul(left_value, right_value),
                )
        return self.trim(answer)

    def pow(self, value, exponent):
        require(exponent >= 0, "negative polynomial power")
        answer = self.one
        base = value
        remaining = exponent
        while remaining:
            if remaining % 2:
                answer = self.mul(answer, base)
            base = self.mul(base, base)
            remaining //= 2
        return answer

    def derivative(self, value):
        return self.trim(
            self.field.mul(self.field.scalar(index), value[index])
            for index in range(1, len(value))
        )

    def degree(self, value):
        return len(value) - 1

    def leading(self, value):
        require(bool(value), "zero polynomial has no leading coefficient")
        return value[-1]

    def divmod(self, numerator, denominator):
        require(bool(denominator), "polynomial division by zero")
        quotient = [self.field.zero] * max(
            0, self.degree(numerator) - self.degree(denominator) + 1
        )
        remainder = numerator
        denominator_degree = self.degree(denominator)
        denominator_lead = self.leading(denominator)
        while remainder and self.degree(remainder) >= denominator_degree:
            shift = self.degree(remainder) - denominator_degree
            coefficient = self.field.div(self.leading(remainder), denominator_lead)
            quotient[shift] = self.field.add(quotient[shift], coefficient)
            subtraction = (self.field.zero,) * shift + self.scale(
                denominator, coefficient
            )
            remainder = self.sub(remainder, subtraction)
        return self.trim(quotient), remainder

    def exquo(self, numerator, denominator):
        quotient, remainder = self.divmod(numerator, denominator)
        require(not remainder, "inexact polynomial quotient")
        return quotient

    def monic(self, value):
        return self.scale(value, self.field.inv(self.leading(value)))

    def gcd(self, left, right):
        first, second = left, right
        while second:
            _, remainder = self.divmod(first, second)
            first, second = second, remainder
        return self.monic(first) if first else ()

    def shift(self, value, amount):
        linear = (amount, self.field.one)
        answer = ()
        for coefficient in reversed(value):
            answer = self.add(self.mul(answer, linear), self.const(coefficient))
        return answer

    def coefficient(self, value, degree):
        return value[degree] if degree < len(value) else self.field.zero


def field_sum(field, *values):
    return reduce(field.add, values, field.zero)


def field_product(field, *values):
    return reduce(field.mul, values, field.one)


def fraction_text(value):
    if value.denominator == 1:
        return str(value.numerator)
    return "%d/%d" % (value.numerator, value.denominator)


def anp_text(field, value):
    coefficients = list(value)
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    numerator = ", ".join(fraction_text(entry)
                           for entry in reversed(coefficients))
    modulus = ", ".join(fraction_text(entry)
                        for entry in reversed(field.q))
    return "ANP([%s], [%s], QQ)" % (numerator, modulus)


def primary_poly_digest(field, polynomial, ring):
    monic = ring.monic(polynomial)
    text = "|".join(
        "%d:%s" % (degree, anp_text(field, coefficient))
        for degree, coefficient in enumerate(monic)
        if coefficient != field.zero
    )
    return sha256(text.encode("ascii")).hexdigest()


CRITICAL_TERMS = tuple(
    (powers, Fraction(int(coefficient)))
    for powers, coefficient in sp.Poly(
        CRITICAL_FACTOR, *CRITICAL_SYMBOLS
    ).terms()
)


def evaluate_critical_factor(ring, arguments):
    field = ring.field
    power_cache = {}
    for argument_index, argument in enumerate(arguments):
        exponents = {powers[argument_index] for powers, _ in CRITICAL_TERMS}
        for exponent in exponents:
            power_cache[(argument_index, exponent)] = ring.pow(argument, exponent)
    answer = ring.zero
    for powers, coefficient in CRITICAL_TERMS:
        term = ring.const(field.scalar(coefficient))
        for argument_index, exponent in enumerate(powers):
            term = ring.mul(term, power_cache[(argument_index, exponent)])
        answer = ring.add(answer, term)
    return answer


def exact_case(name):
    if name == "4111":
        field = CubicField((44, 237, 244, 100))
        exponent_a, exponent_b = 4, 1
        boundary_extras = (9, 0)
    else:
        field = CubicField((61, -31, -89, 75))
        exponent_a, exponent_b = 3, 2
        boundary_extras = (6, 3)
    ring = PolynomialRing(field)
    scalar = field.scalar
    u = field.gen
    u2 = field.mul(u, u)
    if name == "4111":
        accessory_v = field.div(
            field_sum(field, field.mul(scalar(8), u2),
                      field.mul(scalar(9), u), scalar(8)),
            scalar(7),
        )
        shift = field.div(field.mul(scalar(5), field.add(u, field.one)),
                          scalar(7))
        A0 = field.div(field_product(
            field, scalar(80), field.pow(accessory_v, 2), field.add(u, field.one)
        ), scalar(343))
    else:
        accessory_v = field.div(
            field_sum(field, field.mul(scalar(24), u2),
                      field.mul(scalar(-16), u), scalar(-16)),
            scalar(21),
        )
        shift = field.div(field.add(field.mul(scalar(5), u), scalar(-4)),
                          scalar(7))
        A0 = field.div(field_product(
            field, scalar(9), field.pow(accessory_v, 2),
            field.add(field.mul(scalar(5), u), scalar(-4))
        ), scalar(343))
    gamma = field.mul(scalar(-7), A0)
    x = ring.x
    x_minus_one = ring.sub(x, ring.one)
    quadratic = ring.add(
        ring.sub(ring.pow(x, 2), ring.scale(x, u)),
        ring.const(accessory_v),
    )
    D = field_product_polynomial(
        ring, ring.pow(x, exponent_a), ring.pow(x_minus_one, exponent_b), quadratic
    )
    T = field_product_polynomial(ring, x, x_minus_one, quadratic)
    S = ring.add(x, ring.const(shift))
    E_numerator = ring.add(
        ring.add(
            ring.scale(ring.mul(x_minus_one, quadratic), scalar(exponent_a)),
            ring.scale(ring.mul(x, quadratic), scalar(exponent_b)),
        ),
        ring.mul(ring.mul(x, x_minus_one),
                 ring.sub(ring.scale(x, scalar(2)), ring.const(u))),
    )
    E = ring.scale(E_numerator, field.inv(scalar(7)))
    require(ring.add(D, ring.const(A0)) == ring.mul(S, ring.pow(E, 2)),
            ("accessory identity", name))
    V = ring.scale(
        field_product_polynomial(ring, S, D, ring.pow(T, 2)),
        field.div(scalar(4), field.pow(gamma, 2)),
    )
    A = ring.scale(field_product_polynomial(ring, S, E, T),
                   field.div(scalar(2), gamma))
    require((ring.degree(V), ring.degree(A)) == (16, 8),
            ("response degrees", name))
    owner = ring.mul(S, T)
    require(ring.degree(owner) == 5
            and ring.degree(ring.gcd(owner, ring.derivative(owner))) == 0,
            ("owner is not squarefree degree five", name))
    # The displayed construction is V=unit*S*D*T^2, while D and T use the
    # same three factors x, x-1 and the quadratic.  Once ST is squarefree,
    # this exact factorization gives rad(V)=ST without a costly degree-16 gcd.
    response = ring.sub(
        ring.sub(ring.scale(ring.mul(V, ring.derivative(A)), scalar(2)),
                 ring.mul(A, ring.derivative(V))),
        ring.scale(V, scalar(2)),
    )
    require(not response, ("response identity", name))
    shifted_V = ring.shift(V, field.neg(shift))
    v1, v2, v3, v4 = tuple(ring.coefficient(shifted_V, index)
                           for index in range(1, 5))
    require(v1 != field.zero and v2 != field.zero,
            ("nonzero first jets", name))
    wall_c = field.div(field.mul(scalar(-2), v2),
                       field.mul(scalar(3), field.pow(v1, 2)))
    wall_d = field.div(
        field_sum(
            field,
            field.mul(scalar(45), field.pow(v1, 3)),
            field.mul(scalar(-54), field.mul(v1, v3)),
            field.mul(scalar(28), field.pow(v2, 2)),
        ),
        field.mul(scalar(90), field.pow(v1, 3)),
    )
    obstruction = field_sum(
        field,
        field.mul(scalar(315), field.mul(field.pow(v1, 3), v2)),
        field.mul(scalar(360), field.mul(field.pow(v1, 2), v4)),
        field.mul(scalar(-414), field_product(field, v1, v2, v3)),
        field.mul(scalar(148), field.pow(v2, 3)),
    )
    require(wall_c != field.zero and wall_d != field.zero
            and obstruction != field.zero,
            ("wall or obstruction vanished", name))
    C = ring.add(ring.const(wall_c), ring.scale(S, wall_d))
    require(ring.degree(ring.gcd(C, ring.mul(S, T))) == 0,
            ("wall control not unit on ST", name))
    K = evaluate_critical_factor(
        ring, (A, C, V, ring.derivative(V), ring.const(wall_d))
    )
    require(ring.degree(K) == 96, ("critical degree", name))
    expected_lead = field_product(
        field, scalar(32), field.pow(ring.leading(A), 2),
        field.pow(ring.leading(V), 5), field.pow(wall_d, 3)
    )
    require(ring.leading(K) == expected_lead,
            ("critical leading coefficient", name))
    boundary = field_product_polynomial(
        ring, ring.pow(S, 3), ring.pow(T, 8),
        ring.pow(x, boundary_extras[0]),
        ring.pow(x_minus_one, boundary_extras[1]),
    )
    require(ring.degree(boundary) == 44, ("boundary degree", name))
    H = ring.exquo(K, boundary)
    require(ring.degree(H) == 52, ("residual degree", name))
    require(ring.degree(ring.gcd(H, T)) == 0,
            ("residual T escape", name))
    shifted_H = ring.shift(H, field.neg(shift))
    require(ring.coefficient(shifted_H, 0) == field.zero
            and ring.coefficient(shifted_H, 1) == field.zero
            and ring.coefficient(shifted_H, 2) != field.zero,
            ("sharp S order", name))
    shifted_K = ring.shift(K, field.neg(shift))
    expected_q5 = field.div(
        field.mul(scalar(-32), field.mul(v2, obstruction)),
        field.mul(scalar(315), field.pow(v1, 2)),
    )
    require(ring.coefficient(shifted_K, 5) == expected_q5,
            ("fifth jet in custom field", name))
    digest = primary_poly_digest(field, H, ring)
    require(digest == EXPECTED_DIGESTS[name],
            ("independent residual digest", name, digest))
    certificate = sha256(repr((v1, v2, v3, v4, wall_c, wall_d,
                               obstruction)).encode("ascii")).hexdigest()
    return digest, certificate


def field_product_polynomial(ring, *values):
    return reduce(ring.mul, values, ring.one)


CASE_RESULTS = tuple((name,) + exact_case(name)
                     for name in ("4111", "3211"))


print("THM-3279 AFFINE TRANSVERSE INDEPENDENT HOSTILE AUDIT")
print("dependency_hash_checks=%d;primary_script_output_hashes=PASS;theorem_status_not_pinned" %
      len(DEPENDENCIES))
print("source_ast=(primary_assert0_float0,audit_assert0_float0)")
print("localization=raw_minus_R2=+(Vprime*y/2)*R1;gradient_ideal_on_V_nonzero=PASS")
print("universal=Sylvester6x6;Res=V^2*F40;specialized=V^3*K20")
print("affine_degree=96;unique_top=32*A_lead^2*V_lead^5*d^3")
print("T_rows=%s;lead=32*m*(m-1)/(m-2)^2*c*v^3" % (T_ROWS,))
print("S_walls=(c=-2*v2/(3*v1^2),d=(45*v1^3-54*v1*v3+28*v2^2)/(90*v1^3))")
print("S_fifth=-32*v2*(315*v1^3*v2+360*v1^2*v4-414*v1*v2*v3+148*v2^3)/(315*v1^2)")
for name, digest, certificate in CASE_RESULTS:
    print("case=%s custom_Q[u]/q=(v1,v2,R,c,d_nonzero);unit_C_mod_ST;deg_K_H=96_52;ordS_H=2;gcd_H_T=1" % name)
    print("case=%s residual_digest=%s;field_certificate_sha256=%s" %
          (name, digest, certificate))
print("consequence=at_most_2_of_52_residual_units_on_S_and_none_on_T;off_owner_multiplicity_at_least_50")
print("constant_C_lane=inherited_THM3212;nonconstant_C_meeting_g_has_explicit_(alpha,0)_critical_point")
print("scope=B1_constant_E0_affine_C0_only;resultant_multiplicity_not_distinct_points;no_inverse_cover_no_JC2")
print("all_independent_checks=PASS")
