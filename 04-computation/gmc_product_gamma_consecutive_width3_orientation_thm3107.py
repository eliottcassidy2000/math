#!/usr/bin/env python3
"""Exact coefficient certificate for THM-3107.

For the normalized support {0,1,2}, derive THM-2824's two binary
quadratic/cubic divisibility numerators.  Substitute the finite product-Gamma
moment ratio q_n=prod_j(n+theta_j), remove the manifest positive factors, and
check coefficientwise strict negativity for one through six Gamma layers.
"""

import sympy as sp
from flint import fmpz_mpoly_ctx


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


x = sp.symbols("x")
moment_symbols = sp.symbols("w1:7")
w = (sp.Integer(1),) + moment_symbols


def moment(poly):
    return sp.expand(
        sum(
            coefficient * w[degree]
            for (degree,), coefficient in sp.Poly(sp.expand(poly), x).terms()
        )
    )


U = x / w[1] - 1
V = x**2 / w[2] - x / w[1]
g11 = moment(U**2)
g12 = moment(U * V)
g22 = moment(V**2)
t111 = moment(U**3)
t112 = moment(U**2 * V)
t122 = moment(U * V**2)
t222 = moment(V**3)

I1 = 3 * t112 * g11 * g22 - t222 * g11**2 - 2 * t111 * g12 * g22
I2 = 3 * t122 * g11 * g22 - 2 * t222 * g12 * g11 - t111 * g22**2
N1, denominator1 = sp.cancel(I1).as_numer_denom()
N2, denominator2 = sp.cancel(I2).as_numer_denom()
require(denominator1 == w[1] ** 5 * w[2] ** 3, "I1 denominator")
require(denominator2 == w[1] ** 5 * w[2] ** 4, "I2 denominator")


def terms(poly):
    answer = []
    for exponents, coefficient in sp.Poly(poly, *moment_symbols).terms():
        require(coefficient.is_Integer, "integral universal numerator")
        answer.append((int(coefficient), tuple(int(value) for value in exponents)))
    return answer


term_data = (
    ("P", terms(N1), 4, 1, 16),
    ("Q", terms(N2), 5, 2, 12),
)
require(len(term_data[0][1]) == 20, "N1 term census")
require(len(term_data[1][1]) == 23, "N2 term census")

expected_terms = (0, 2, 34, 308, 2331, 16681, 117439)
records = []

for layer_count in range(1, 7):
    names = tuple(f"theta_{index}" for index in range(1, layer_count + 1))
    ring = fmpz_mpoly_ctx.get(names)
    variables = ring.gens()

    ratios = []
    for n in range(6):
        ratio = ring.constant(1)
        for variable in variables:
            ratio *= variable + n
        ratios.append(ratio)

    moments = [ring.constant(1)]
    for n in range(6):
        moments.append(moments[-1] * ratios[n])

    layer_record = []
    for label, universal_terms, q0_power, q1_power, expected_minimum in term_data:
        numerator = ring.constant(0)
        for coefficient, exponents in universal_terms:
            monomial = ring.constant(coefficient)
            for index, exponent in enumerate(exponents, start=1):
                if exponent:
                    monomial *= moments[index] ** exponent
            numerator += monomial

        positive_divisor = ratios[0] ** q0_power * ratios[1] ** q1_power
        quotient, remainder = divmod(numerator, positive_divisor)
        require(not remainder, f"A={layer_count} {label} manifest exact division")

        oriented = -quotient
        coefficients = [int(coefficient) for coefficient in oriented.coeffs()]
        require(coefficients, f"A={layer_count} {label} nonempty quotient")
        require(all(coefficient > 0 for coefficient in coefficients),
                f"A={layer_count} {label} coefficientwise orientation")
        require(len(oriented) == expected_terms[layer_count],
                f"A={layer_count} {label} term census")
        require(oriented.total_degree() == 6 * layer_count - 5,
                f"A={layer_count} {label} total degree")
        require(min(coefficients) == expected_minimum,
                f"A={layer_count} {label} sharp coefficient floor")

        # The certificate is symmetric in the labelled Gamma layers.
        if layer_count >= 2:
            swap = list(variables)
            swap[0], swap[-1] = swap[-1], swap[0]
            require(oriented.compose(*swap) == oriented,
                    f"A={layer_count} {label} endpoint permutation symmetry")

        layer_record.append((len(oriented), oriented.total_degree(), min(coefficients)))
    require(layer_record[0][:2] == layer_record[1][:2],
            f"A={layer_count} common support census")
    records.append(tuple(layer_record))

# Closed one-layer controls and exact shape-two boundary from THM-3100.
theta = sp.symbols("theta")
one_layer_P = 16 * theta + 40
one_layer_Q = 12 * theta + 20
require(one_layer_P.subs(theta, sp.Rational(-5, 2)) == 0,
        "P sign-boundary hostile")
require(one_layer_Q.subs(theta, sp.Rational(-5, 3)) == 0,
        "Q sign-boundary hostile")
require(sp.cancel(I1.subs({
    moment_symbols[index - 1]: sp.rf(2, index) for index in range(1, 7)
})) == sp.Rational(-1, 2), "shape-two I1 control")
require(sp.cancel(I2.subs({
    moment_symbols[index - 1]: sp.rf(2, index) for index in range(1, 7)
})) == sp.Rational(-11, 36), "shape-two I2 control")

print("universal_numerators=N1_terms20;N2_terms23")
print("denominators=I1:w1^5*w2^3;I2:w1^5*w2^4")
for layer_count, record in enumerate(records, start=1):
    p_record, q_record = record
    print(
        f"A={layer_count} "
        f"P_terms={p_record[0]} Q_terms={q_record[0]} "
        f"degree={p_record[1]} min_coefficients={p_record[2]}/{q_record[2]}"
    )
print("factorization=N1=-q0^4*q1*P_A;N2=-q0^5*q1^2*Q_A")
print("one_layer=P1=16*theta+40;Q1=12*theta+20")
print("shape2_support012=I1:-1/2;I2:-11/36")
print("sign_hostiles=P1(-5/2)=0;Q1(-5/3)=0")
print("conclusion=positive_shapes_A1_to_A6_support012_detected")
print("boundary=A0_degenerate;A7plus_open;nonconsecutive_open;soft_MID_insufficient")
