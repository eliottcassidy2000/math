#!/usr/bin/env python3
"""Exact companion for THM-3470.

Universe: the normalized standard triangle with the three-real-parameter
Hermitian Keller shear

    U = x-y,
    M = x+y-2/3,
    g = a M + b(U^2-1/6) + c(U^4-1/15) + i U.

There is no floating-point arithmetic, random sampling, or assertion-sensitive
truth gate in this file.  The printed decimal boxes are outward-rounded views
of exact rational intervals.
"""

from hashlib import sha256

import sympy as sp


a, b, c, M, U = sp.symbols("a b c M U", real=True)
x, y = sp.symbols("x y", real=True)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def triangle_monomial(i, j):
    """Normalized average of x^i y^j on x,y>=0, x+y<=1."""
    return sp.Rational(2) * sp.factorial(i) * sp.factorial(j) / sp.factorial(i + j + 2)


def direct_MU_average(j, k):
    """Independent binomial expansion in x,y, without the (s,u) integral."""
    answer = sp.S.Zero
    for h in range(j + 1):
        centered_coefficient = sp.binomial(j, h) * (-sp.Rational(2, 3)) ** (j - h)
        for r in range(h + 1):
            sum_coefficient = sp.binomial(h, r)
            for ell in range(k + 1):
                difference_coefficient = sp.binomial(k, ell) * (-1) ** (k - ell)
                answer += (
                    centered_coefficient
                    * sum_coefficient
                    * difference_coefficient
                    * triangle_monomial(r + ell, h - r + k - ell)
                )
    return sp.factor(answer)


def mu_MU(j, k):
    """Normalized <M^j U^k>, from s=x+y and u=(x-y)/(x+y)."""
    if k % 2:
        return sp.S.Zero
    return sp.factor(
        sum(
            sp.binomial(j, h)
            * (-sp.Rational(2, 3)) ** (j - h)
            * sp.Rational(2, (h + k + 2) * (k + 1))
            for h in range(j + 1)
        )
    )


def MU_average(poly):
    poly = sp.Poly(sp.expand(poly), M, U)
    return sp.factor(sum(coeff * mu_MU(j, k) for (j, k), coeff in poly.terms()))


def primitive_numerator(poly):
    denominator, integral = sp.Poly(poly, a, b, c, domain=sp.QQ).clear_denoms()
    content, primitive = integral.primitive()
    return denominator, content, primitive


def interval_multiply(left, right):
    values = (
        left[0] * right[0],
        left[0] * right[1],
        left[1] * right[0],
        left[1] * right[1],
    )
    return min(values), max(values)


def polynomial_interval(poly, interval):
    """Exact Horner interval evaluation for a univariate QQ polynomial."""
    coefficients = poly.all_coeffs()
    result = (coefficients[0], coefficients[0])
    for coefficient in coefficients[1:]:
        result = interval_multiply(result, interval)
        result = (result[0] + coefficient, result[1] + coefficient)
    return result


def outward_decimal_box(interval, digits):
    denominator = 10**digits
    lower = int(sp.floor(interval[0] * denominator))
    upper = int(sp.ceiling(interval[1] * denominator))
    return (lower, upper, denominator)


def polynomial_digest(*polynomials):
    payload = "\n".join(str(sp.Poly(poly).as_expr()) for poly in polynomials)
    return sha256(payload.encode("ascii")).hexdigest()


# Independent audit of the (M,U) table against direct simplex monomials.
table_checks = 0
for total in range(21):
    for j in range(total + 1):
        k = total - j
        direct = direct_MU_average(j, k)
        require(direct == mu_MU(j, k), f"mixed-moment mismatch at {(j, k)}")
        table_checks += 1


g = (
    a * M
    + b * (U**2 - sp.Rational(1, 6))
    + c * (U**4 - sp.Rational(1, 15))
    + sp.I * U
)
moments = {m: sp.factor(MU_average(g**m)) for m in range(1, 6)}
require(moments[1] == 0, "the family was not centered")
require(all(sp.im(moment) == 0 for moment in moments.values()), "a moment was not real")

denominators = {}
primitive_moments = {}
for m in range(2, 6):
    denominator, content, primitive = primitive_numerator(moments[m])
    require(content == 1, f"unexpected numerator content at moment {m}")
    denominators[m] = denominator
    primitive_moments[m] = primitive.as_expr()
    require(
        sp.factor(moments[m] - primitive.as_expr() / denominator) == 0,
        f"primitive normalization failed at moment {m}",
    )


# Direct Jacobian audit in the original x,y coordinates.
U_xy = x - y
M_xy = x + y - sp.Rational(2, 3)
p_xy = (
    a * M_xy
    + b * (U_xy**2 - sp.Rational(1, 6))
    + c * (U_xy**4 - sp.Rational(1, 15))
)
q_xy = U_xy
jacobian = sp.factor(sp.diff(p_xy, x) * sp.diff(q_xy, y) - sp.diff(p_xy, y) * sp.diff(q_xy, x))
require(jacobian == -2 * a, "constant-Jacobian identity failed")


P2, P3, P4, P5 = (primitive_moments[m] for m in range(2, 6))

# First exact elimination: the moment-five gate itself.
gate_basis = sp.groebner((P2, P3, P4, P5), a, b, c, order="grevlex")
require(len(gate_basis.polys) == 1, "moment-five basis was not a singleton")
require(gate_basis.polys[0].as_expr() == 1, "moments two through five did not generate one")

# Second exact path: triangularize the moment-two-through-four locus, then
# prove that the moment-five remainder is coprime to its terminal polynomial.
prefix_basis = sp.groebner((P2, P3, P4), a, b, c, order="lex")
require(len(prefix_basis.polys) == 3, "prefix lex basis did not have triangular length three")
lex_a, lex_b, lex_c = (poly.as_expr() for poly in prefix_basis.polys)
require(sp.degree(lex_a, a) == 1 and not lex_a.has(b), "first lex equation was not a-linear")
require(sp.degree(lex_b, b) == 1 and not lex_b.has(a), "second lex equation was not b-linear")
require(not lex_c.has(a) and not lex_c.has(b), "terminal lex equation was not univariate")
require(not sp.diff(lex_a, a).has(c), "a-linear coefficient depended on c")
require(not sp.diff(lex_b, b).has(c), "b-linear coefficient depended on c")

A_of_c = sp.Poly(-lex_a.subs(a, 0) / sp.diff(lex_a, a), c, domain=sp.QQ)
B_of_c = sp.Poly(-lex_b.subs(b, 0) / sp.diff(lex_b, b), c, domain=sp.QQ)
C_terminal = sp.Poly(lex_c, c, domain=sp.QQ)
require(C_terminal.degree() == 24, "terminal polynomial did not have degree 24")
require(sp.gcd(C_terminal, C_terminal.diff()).degree() == 0, "terminal polynomial was not squarefree")
require(all(exponent[0] % 2 == 0 for exponent, _ in C_terminal.terms()), "terminal polynomial was not even")

P5_remainder = sp.Poly(prefix_basis.reduce(P5)[1], c, domain=sp.QQ)
require(P5_remainder.degree() == 23, "unexpected moment-five remainder degree")
require(sp.gcd(C_terminal, P5_remainder).degree() == 0, "moment-five remainder met the prefix locus")

# All real roots are isolated exactly.  Each rational interval below contains
# exactly one root, and the full exact interval routine finds no others.
root_scale = 10**14
positive_root_boxes = (
    (sp.Rational(216566442254669, root_scale), sp.Rational(216566442254671, root_scale)),
    (sp.Rational(393762058418513, root_scale), sp.Rational(393762058418515, root_scale)),
)
root_boxes = tuple((-upper, -lower) for lower, upper in positive_root_boxes) + positive_root_boxes
isolated_real_roots = sp.intervals(C_terminal, eps=sp.Rational(1, 10**12))
require(sum(multiplicity for _, multiplicity in isolated_real_roots) == 4, "terminal real-root count was not four")
require(len(isolated_real_roots) == 4, "terminal real roots were not distinct")
for interval in root_boxes:
    require(C_terminal.count_roots(interval[0], interval[1]) == 1, f"bad root box {interval}")

# At a prefix root, moment five is its exact lex remainder divided by the
# displayed denominator.  Interval evaluation supplies rigorous point and
# failure boxes without introducing floating point.
M5_on_prefix = sp.Poly(P5_remainder.as_expr() / denominators[5], c, domain=sp.QQ)
real_point_boxes = []
for c_interval in root_boxes:
    a_interval = polynomial_interval(A_of_c, c_interval)
    b_interval = polynomial_interval(B_of_c, c_interval)
    m5_interval = polynomial_interval(M5_on_prefix, c_interval)
    require(a_interval[0] * a_interval[1] > 0, "a-box met the non-Keller wall")
    require(m5_interval[0] * m5_interval[1] > 0, "moment-five box met zero")
    real_point_boxes.append(
        (
            outward_decimal_box(a_interval, 9),
            outward_decimal_box(b_interval, 9),
            outward_decimal_box(c_interval, 14),
            outward_decimal_box(m5_interval, 9),
        )
    )

# Quadratic boundary control: c=0 has real moment-one-through-three survivors,
# but moment four already closes that lower-degree face.
Q2, Q3, Q4 = (sp.expand(poly.subs(c, 0)) for poly in (P2, P3, P4))
quadratic_prefix = sp.groebner((Q2, Q3), a, b, order="lex")
quadratic_gate = sp.groebner((Q2, Q3, Q4), a, b, order="grevlex")
require(len(quadratic_prefix.polys) == 2, "quadratic prefix control lost its finite locus")
require(quadratic_prefix.is_zero_dimensional, "quadratic prefix control was not zero-dimensional")
require(len(quadratic_gate.polys) == 1 and quadratic_gate.polys[0].as_expr() == 1, "quadratic moment-four gate failed")
quadratic_a_equation, quadratic_b_equation = (poly.as_expr() for poly in quadratic_prefix.polys)
require(sp.degree(quadratic_a_equation, a) == 1, "quadratic control was not a-linear")
require(not sp.diff(quadratic_a_equation, a).has(b), "quadratic a-linear coefficient depended on b")
quadratic_A_of_b = sp.Poly(
    -quadratic_a_equation.subs(a, 0) / sp.diff(quadratic_a_equation, a),
    b,
    domain=sp.QQ,
)
quadratic_terminal = sp.Poly(quadratic_b_equation, b, domain=sp.QQ)
require(quadratic_terminal.degree() == 6, "quadratic terminal polynomial did not have degree six")
require(sp.gcd(quadratic_terminal, quadratic_terminal.diff()).degree() == 0, "quadratic terminal was not squarefree")
quadratic_isolations = sp.intervals(quadratic_terminal, eps=sp.Rational(1, 10**12))
require(sum(multiplicity for _, multiplicity in quadratic_isolations) == 2, "quadratic real-root count was not two")
quadratic_scale = 10**14
quadratic_positive_box = (
    sp.Rational(118352913031470, quadratic_scale),
    sp.Rational(118352913031480, quadratic_scale),
)
quadratic_root_boxes = (
    (-quadratic_positive_box[1], -quadratic_positive_box[0]),
    quadratic_positive_box,
)
quadratic_M4_remainder = sp.Poly(
    quadratic_prefix.reduce(Q4)[1] / denominators[4], b, domain=sp.QQ
)
quadratic_real_point_boxes = []
for b_interval in quadratic_root_boxes:
    require(quadratic_terminal.count_roots(*b_interval) == 1, "bad quadratic root box")
    a_interval = polynomial_interval(quadratic_A_of_b, b_interval)
    m4_interval = polynomial_interval(quadratic_M4_remainder, b_interval)
    require(a_interval[0] * a_interval[1] > 0, "quadratic witness met the non-Keller wall")
    require(m4_interval[0] * m4_interval[1] > 0, "quadratic fourth moment met zero")
    quadratic_real_point_boxes.append(
        (
            outward_decimal_box(a_interval, 9),
            outward_decimal_box(b_interval, 14),
            outward_decimal_box(m4_interval, 9),
        )
    )

# The c U^4 top term has nonzero projections to all three C3 characters:
# U is a nonzero combination of the conjugate source characters, and the five
# terms of U^4 have weights 0,1,2 modulo three.
degree_four_character_weights = sorted({(2 * k - 4) % 3 for k in range(5)})
require(degree_four_character_weights == [0, 1, 2], "quartic top term was not mixed-character")


print("THM-3470 MIXED-CHARACTER TRIANGULAR QUARTIC MOMENT-FIVE GATE")
print("universe=normalized_standard_triangle;real_a_b_c;g=a*M+b*(U^2-1/6)+c*(U^4-1/15)+i*U")
print("coordinates=U=x-y;M=x+y-2/3;dagger=coefficient_conjugation")
print(f"mixed_moment_table_checks={table_checks};total_degree<=20;independent_direct_simplex_path=PASS")
print(f"centering=M1_identically_{moments[1]};reflection_all_M1_to_M5_real=PASS")
print(f"real_jacobian={jacobian};keller_iff=a_nonzero")
for m in range(2, 6):
    print(f"M{m}=P{m}/{denominators[m]}")
    print(f"P{m}={primitive_moments[m]}")
print("moment_five_gate_grevlex_basis=[1]")
print(
    "prefix_lex_shape="
    f"a_linear,b_linear,c_degree_{C_terminal.degree()};"
    f"degrees_A_B=({A_of_c.degree()},{B_of_c.degree()});"
    "terminal_squarefree=PASS;terminal_even=PASS"
)
print(f"prefix_lex_basis_sha256={polynomial_digest(lex_a, lex_b, lex_c)}")
print(f"terminal_polynomial_sha256={polynomial_digest(C_terminal.as_expr())}")
print(f"M5_remainder_degree={P5_remainder.degree()};gcd_terminal_remainder=1")
print("M2_M3_M4_complex_locus=24_distinct_reduced_points")
print(f"M2_M3_M4_real_locus_count={len(root_boxes)}")
for index, point_box in enumerate(real_point_boxes, start=1):
    print(f"real_point_{index}_boxes_a_b_c_M5={point_box}")
print("each_real_point_has_nonzero_a=PASS;each_real_point_has_nonzero_M5=PASS")
print("quartic_top_C3_character_weights=(0,1,2);genuinely_mixed_when_c_nonzero")
print("quadratic_boundary=M1_to_M3_two_real_survivors;M1_to_M4_grevlex_basis_[1]")
for index, point_box in enumerate(quadratic_real_point_boxes, start=1):
    print(f"quadratic_real_point_{index}_boxes_a_b_M4={point_box}")
print(f"semantic_sha256={polynomial_digest(P2, P3, P4, P5, lex_c, P5_remainder.as_expr())}")
print("STATUS=PASS")
