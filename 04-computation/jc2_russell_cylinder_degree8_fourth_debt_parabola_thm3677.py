"""Exact symbolic companion for THM-3677."""

from math import factorial

import sympy as sp


CHECKS = 0


def require(label, condition):
    global CHECKS
    if not condition:
        raise RuntimeError(f"FAILED: {label}")
    CHECKS += 1


x, q_source, t, X = sp.symbols("x q_source t X")
p, r = sp.symbols("p r")
points = (-1, 0, 1)
lambda_row = (sp.Rational(5, 18), -1, sp.Rational(13, 18))

Q1 = (
    x**5
    + sp.Rational(9, 2) * x**4
    - 2 * x**3
    - sp.Rational(27, 4) * x**2
    + x
    - sp.Rational(3, 4)
)
P = sp.expand(x**2 * (x**2 - 1) ** 2)
Q6 = sp.expand(Q1 - sp.Rational(259, 36) * P)
R1 = sp.expand(P * (1 - x**2))
R2 = sp.expand(P * (4 - 9 * x))
Qfamily = sp.expand(Q6 + p * R1 + r * R2)

Qstar = (
    -x**7
    - sp.Rational(27, 4) * x**6
    + 3 * x**5
    + 18 * x**4
    - 3 * x**3
    - sp.Rational(27, 2) * x**2
    + x
    - sp.Rational(3, 4)
)

D = 1 + x**2 * q_source
b = sp.expand((D - 1) * (D + 2) ** 2)
c = sp.expand(x * D * (D + 2))
e = sp.expand(q_source * (D + 3))
y_general = c / 3
z_general = e + 3


def multiply(left, right, cutoff=4):
    answer = {}
    for (i, j), left_value in left.items():
        for (u, v), right_value in right.items():
            if i + j + u + v > cutoff:
                continue
            key = (i + u, j + v)
            answer[key] = answer.get(key, 0) + left_value * right_value
    return {key: sp.factor(value) for key, value in answer.items() if value}


def power(value, exponent):
    answer = {(0, 0): sp.S.One}
    for _ in range(exponent):
        answer = multiply(answer, value)
    return answer


def jet(expr, point):
    shifted = sp.Poly(
        sp.expand(expr.subs(x, X + point)),
        X,
        t,
        domain=sp.QQ.frac_field(p, r),
    )
    return {
        (source_degree, stable_degree): sp.factor(coefficient)
        for (source_degree, stable_degree), coefficient in shifted.terms()
        if source_degree + stable_degree <= 4
    }


print("THM-3677 exact companion -- degree-eight fourth-debt parabola")
print("status=PROVED VERIFIED-EXACT PENDING-INDEPENDENT-HOSTILE-AUDIT")

require("compiler surface relation", sp.expand(c**2 * e - b * (b + 4)) == 0)

values = tuple(sp.factor(Qfamily.subs(x, point)) for point in points)
slopes = tuple(sp.factor(sp.diff(Qfamily, x).subs(x, point)) for point in points)
curvatures = tuple(sp.factor(sp.diff(Qfamily, x, 2).subs(x, point)) for point in points)
require("family retained values", values == (-3, sp.Rational(-3, 4), -3))
require("family retained slopes", slopes == (sp.Rational(-9, 2), 1, sp.Rational(9, 2)))
require("family zero second debt", sp.factor(5 * curvatures[0] + 13 * curvatures[2] + 243) == 0)

# Completeness inside degree <=8: six value/slope rows plus the second-debt
# row have rank seven, and R1,R2 are an independent nullspace basis.
monomials = [x**degree for degree in range(9)]
rows = []
for point in points:
    rows.append([monomial.subs(x, point) for monomial in monomials])
    rows.append([sp.diff(monomial, x).subs(x, point) for monomial in monomials])
rows.append(
    [
        5 * sp.diff(monomial, x, 2).subs(x, -1)
        + 13 * sp.diff(monomial, x, 2).subs(x, 1)
        for monomial in monomials
    ]
)
constraint_matrix = sp.Matrix(rows)
require("degree-eight constraint rank", constraint_matrix.rank() == 7)
coefficient = lambda polynomial: sp.Matrix(
    [sp.Poly(polynomial, x).nth(degree) for degree in range(9)]
)
require("R1 in constraint kernel", constraint_matrix * coefficient(R1) == sp.zeros(7, 1))
require("R2 in constraint kernel", constraint_matrix * coefficient(R2) == sp.zeros(7, 1))
require("R1 R2 independent", sp.Matrix.hstack(coefficient(R1), coefficient(R2)).rank() == 2)
print("PASS family_complete=Q6+pR1+rR2_dimension2_zero_second_debt")

pulled_y = sp.expand(y_general.subs(q_source, Qfamily + t**2))
pulled_z = sp.expand(z_general.subs(q_source, Qfamily + t**2))
y_x = sp.diff(pulled_y, x)
y_t = sp.diff(pulled_y, t)
z_x = sp.diff(pulled_z, x)
z_t = sp.diff(pulled_z, t)
area = sp.expand(y_x * z_t - y_t * z_x)

packets = {}
for point in points:
    y_jet = jet(pulled_y, point)
    z_jet = jet(pulled_z, point)
    packets[point] = (
        [power(y_jet, degree) for degree in range(5)],
        [power(z_jet, degree) for degree in range(5)],
        (jet(area, point), jet(y_x, point), jet(z_x, point)),
    )


def monomial_pullback(kind, y_degree, z_degree, w_degree, point):
    y_powers, z_powers, bases = packets[point]
    value = multiply(y_powers[y_degree], z_powers[z_degree])
    value = multiply(value, {(0, w_degree): sp.S.One})
    return multiply(value, bases[kind])


Aminus = 8 * (
    1819584 * p * r
    + 404352 * p
    - 89092224 * r**2
    - 13115349 * r
    + 8006926
) / 3326427
Bminus = (177840 * r - 26159) / 28431
Azero = -8 * (
    1819584 * p * r
    - 2552472 * p
    + 81746496 * r**2
    - 74734101 * r
    - 15181226
) / 3326427
Bzero = 64 * (9 * r - 1) / 81
Bplus = -13 * (144 * r + 65) / 2187
Cminus = -8 * (2340 * r - 503) / 3159
Czero = 8 * (2340 * r - 503) / 3159

# Coefficients multiply Taylor coefficients J_s^(d)/d!, not raw derivatives.
identity = {
    (0, -1, 0): Aminus,
    (0, -1, 1): Bminus,
    (0, -1, 2): -sp.Rational(10, 81),
    (0, 0, 0): Azero,
    (0, 0, 1): Bzero,
    (0, 1, 1): Bplus,
    (0, 1, 2): -sp.Rational(26, 81),
    (2, -1, 0): Cminus,
    (2, -1, 1): sp.Rational(5, 27),
    (2, 0, 0): Czero,
    (2, 1, 1): -sp.Rational(13, 27),
}

monomial_count = 0
mutation_detected = False
for kind in range(3):
    for y_degree in range(5):
        for z_degree in range(5 - y_degree):
            for w_degree in range(5 - y_degree - z_degree):
                branch = {
                    point: monomial_pullback(kind, y_degree, z_degree, w_degree, point)
                    for point in points
                }
                left = sum(
                    lambda_row[index] * branch[point].get((0, 4), 0)
                    for index, point in enumerate(points)
                )
                right = sum(
                    coefficient_value
                    * branch[point].get((source_degree, stable_degree), 0)
                    for (stable_degree, point, source_degree), coefficient_value in identity.items()
                )
                require(
                    f"symbolic two-form monomial {kind,y_degree,z_degree,w_degree}",
                    sp.factor(left - right) == 0,
                )
                mutated = right + branch[-1].get((0, 0), 0)
                if sp.factor(left - mutated) != 0:
                    mutation_detected = True
                monomial_count += 1

require("complete two-form universe", monomial_count == 105)
require("active identity mutation", mutation_detected)
debt = sp.factor(Aminus + Azero)
expected_debt = 64 * (729 * p - 42120 * r**2 + 15192 * r + 5717) / 6561
require("fourth debt formula", sp.factor(debt - expected_debt) == 0)
print(f"PASS symbolic_two_form_monomials={monomial_count}_fourth_debt={expected_debt}")

require("Q6 specialization", sp.factor(debt.subs({p: 0, r: 0}) - sp.Rational(365888, 6561)) == 0)
require("Qstar family coordinate", sp.expand(Qfamily.subs({p: 0, r: sp.Rational(1, 9)}) - Qstar) == 0)
require("Qstar specialization", sp.factor(debt.subs({p: 0, r: sp.Rational(1, 9)}) - sp.Rational(5440, 81)) == 0)

p_parabola = sp.Rational(520, 9) * r**2 - sp.Rational(1688, 81) * r - sp.Rational(5717, 729)
require("zero parabola", sp.factor(debt.subs(p, p_parabola)) == 0)
Qdagger = sp.expand(Qfamily.subs({p: -sp.Rational(5717, 729), r: 0}))
Qdagger_expected = sp.expand(
    (
        22868 * x**8
        - 89583 * x**6
        + 2916 * x**5
        + 123684 * x**4
        - 5832 * x**3
        - 63530 * x**2
        + 2916 * x
        - 2187
    )
    / 2916
)
require("Qdagger exact polynomial", sp.expand(Qdagger - Qdagger_expected) == 0)
require("Qdagger degree", sp.Poly(Qdagger, x).degree() == 8)
print("PASS zero_fourth_debt_parabola=p=520r^2/9-1688r/81-5717/729")
print("PASS Qdagger=(22868x8-89583x6+2916x5+123684x4-5832x3-63530x2+2916x-2187)/2916")

print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
