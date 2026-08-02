#!/usr/bin/env python3
"""Rigorous symbolic/Arb certificate for THM-3104.

The proof uses exact SymPy algebra to derive the two width-three divisibility
numerators and 256-bit Arb balls to certify a Poincare--Miranda box.  The
printed transcript is deliberately made of fixed rational bounds, so normal
and optimized runs are byte-stable.
"""

from fractions import Fraction

import sympy as sp
from flint import arb, ctx


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


# Derive THM-2824's two division-free invariants from the moment symbols.
x = sp.symbols("x")
w_symbols = sp.symbols("w1:7", positive=True)
w = (sp.Integer(1),) + w_symbols


def moment(poly):
    expansion = sp.Poly(sp.expand(poly), x)
    return sp.expand(
        sum(coefficient * w[degree] for (degree,), coefficient in expansion.terms())
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
N1, D1 = sp.cancel(I1).as_numer_denom()
N2, D2 = sp.cancel(I2).as_numer_denom()
require(sp.expand(D1 - w[1] ** 5 * w[2] ** 3) == 0, "I1 denominator")
require(sp.expand(D2 - w[1] ** 5 * w[2] ** 4) == 0, "I2 denominator")


def terms(poly):
    answer = []
    for exponents, coefficient in sp.Poly(poly, *w_symbols).terms():
        require(coefficient.is_Integer, "integer numerator coefficient")
        answer.append((int(coefficient), tuple(int(value) for value in exponents)))
    return answer


n1_terms = terms(N1)
n2_terms = terms(N2)
require(len(n1_terms) == 20 and len(n2_terms) == 23, "numerator term census")
require(sp.simplify(sp.gamma(sp.Rational(-3, 2)) - 4 * sp.sqrt(sp.pi) / 3) == 0,
        "Levy normalization")


# Each monomial c product_n w_n^e becomes
# c exp(a sum e_n n^(3/2) + b sum e_n n^2).
ctx.prec = 256
stable_exponents = [arb(n) * arb(n).sqrt() for n in range(1, 7)]


def exponential_modes(term_list):
    answer = []
    for coefficient, exponents in term_list:
        stable = sum(
            (exponents[index] * stable_exponents[index] for index in range(6)),
            arb(0),
        )
        gaussian = arb(
            sum(exponents[index] * (index + 1) ** 2 for index in range(6))
        )
        answer.append((coefficient, stable, gaussian))
    return answer


F_modes = exponential_modes(n1_terms)
N2_modes = exponential_modes(n2_terms)
H_modes = (
    [(9 * coefficient, stable, gaussian) for coefficient, stable, gaussian in N2_modes]
    + [(-10 * coefficient, stable, gaussian) for coefficient, stable, gaussian in F_modes]
)


def evaluate(mode_list, a_value, b_value, da=0, db=0):
    return sum(
        (
            coefficient
            * stable**da
            * gaussian**db
            * (a_value * stable + b_value * gaussian).exp()
            for coefficient, stable, gaussian in mode_list
        ),
        arb(0),
    )


a0_fraction = Fraction(101128546502785376754010485709, 10**31)
b0_fraction = Fraction(193641702957876518308716902515, 10**31)
ha_fraction = Fraction(1, 10**12)
hb_fraction = Fraction(1, 10**13)
a0 = arb(f"{a0_fraction.numerator}/{a0_fraction.denominator}")
b0 = arb(f"{b0_fraction.numerator}/{b0_fraction.denominator}")
ha = arb(f"{ha_fraction.numerator}/{ha_fraction.denominator}")
hb = arb(f"{hb_fraction.numerator}/{hb_fraction.denominator}")
a_box = arb(a0, ha)
b_box = arb(b0, hb)


def centered_face_enclosures(mode_list, axis):
    value = evaluate(mode_list, a0, b0)
    derivative_a = evaluate(mode_list, a0, b0, 1, 0)
    derivative_b = evaluate(mode_list, a0, b0, 0, 1)
    bound_aa = evaluate(mode_list, a_box, b_box, 2, 0).abs_upper()
    bound_ab = evaluate(mode_list, a_box, b_box, 1, 1).abs_upper()
    bound_bb = evaluate(mode_list, a_box, b_box, 0, 2).abs_upper()
    remainder = (
        bound_aa * ha**2 + 2 * bound_ab * ha * hb + bound_bb * hb**2
    ) / 2
    error = arb(0, remainder)
    if axis == "a":
        negative_face = value - ha * derivative_a + arb(0, hb) * derivative_b + error
        positive_face = value + ha * derivative_a + arb(0, hb) * derivative_b + error
    else:
        negative_face = value + arb(0, ha) * derivative_a - hb * derivative_b + error
        positive_face = value + arb(0, ha) * derivative_a + hb * derivative_b + error
    return value, derivative_a, derivative_b, remainder, negative_face, positive_face


F_data = centered_face_enclosures(F_modes, "a")
H_data = centered_face_enclosures(H_modes, "b")
F_value, _, _, F_remainder, F_left, F_right = F_data
H_value, _, _, H_remainder, H_bottom, H_top = H_data


def qarb(numerator, denominator):
    return arb(f"{numerator}/{denominator}")


# Fixed conservative rational enclosures for all four faces.
require(F_left.lower() > qarb(-713, 10**18), "N1 left lower bound")
require(F_left.upper() < qarb(-4, 10**16), "N1 left sign")
require(F_right.lower() > qarb(4, 10**16), "N1 right sign")
require(F_right.upper() < qarb(713, 10**18), "N1 right upper bound")
require(H_bottom.lower() > qarb(-324, 10**18), "H bottom lower bound")
require(H_bottom.upper() < qarb(-322, 10**18), "H bottom sign")
require(H_top.lower() > qarb(322, 10**18), "H top sign")
require(H_top.upper() < qarb(324, 10**18), "H top upper bound")
require(F_remainder < qarb(18, 10**27), "N1 Taylor remainder")
require(H_remainder < qarb(12, 10**26), "H Taylor remainder")
require(abs(F_value).upper() < qarb(3, 10**35), "N1 center residual")
require(abs(H_value).upper() < qarb(2, 10**34), "H center residual")


# Strict derivative and Jacobian bounds make the box zero unique.
F_a_box = evaluate(F_modes, a_box, b_box, 1, 0)
F_b_box = evaluate(F_modes, a_box, b_box, 0, 1)
H_a_box = evaluate(H_modes, a_box, b_box, 1, 0)
H_b_box = evaluate(H_modes, a_box, b_box, 0, 1)
jacobian = F_a_box * H_b_box - F_b_box * H_a_box
require(F_a_box.lower() > qarb(1, 2000), "N1 increasing in a")
require(jacobian.lower() > qarb(17, 10**7), "positive Jacobian")


# The normalized U,V Gram form remains strictly positive on the whole box.
def moment_ball(index):
    return (
        a_box * stable_exponents[index - 1] + b_box * index**2
    ).exp()


w1, w2, w3, w4 = (moment_ball(index) for index in range(1, 5))
gram11 = (w2 - w1**2) / w1**2
gram12 = (w1 * w3 - w2**2) / (w1**2 * w2)
gram22 = (w1**2 * w4 - 2 * w1 * w2 * w3 + w2**3) / (w1**2 * w2**2)
gram_det = gram11 * gram22 - gram12**2
require(gram11.lower() > qarb(241, 5000), "g11 lower bound")
require(gram11.upper() < qarb(483, 10000), "g11 upper bound")
require(gram22.lower() > qarb(263, 5000), "g22 lower bound")
require(gram22.upper() < qarb(527, 10000), "g22 upper bound")
require(gram_det.lower() > qarb(295, 10**6), "Gram determinant lower bound")
require(gram_det.upper() < qarb(297, 10**6), "Gram determinant upper bound")


print("THM-3104 MULTIPLICATIVE-ID STIELTJES WIDTH-THREE NULLITY")
print("symbolic_I_numerators=20,23;denominators=w1^5*w2^3,w1^5*w2^4")
print("levy_normalization=Gamma(-3/2)=4*sqrt(pi)/3")
print(f"a_center={a0_fraction};a_radius={ha_fraction}")
print(f"b_center={b0_fraction};b_radius={hb_fraction}")
print("N1_left=(-7.13e-16,-4.00e-16);N1_right=(4.00e-16,7.13e-16)")
print("H_bottom=(-3.24e-16,-3.22e-16);H_top=(3.22e-16,3.24e-16)")
print("taylor_remainders=N1<1.8e-26,H<1.2e-25")
print("poincare_miranda_common_zero=PASS")
print("N1_a>1/2000;jacobian_d(N1,H)>17/10^7;unique_box_zero=PASS")
print("g11=(0.0482,0.0483);g22=(0.0526,0.0527)")
print("gram_determinant=(0.000295,0.000297);strict_positive=PASS")
print("scope=full_support_strict_Stieltjes_and_multiplicative_ID_not_product_Gamma")
print("all_rigorous_symbolic_and_Arb_checks=PASS")
