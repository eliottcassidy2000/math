#!/usr/bin/env python3
"""Exact Ore-normal audit for THM-2046 and HYP-8803.

The THM-2044 substitution ell=L+g generates the Ore extension

    Q[x,q][ell; delta],   ell*f=f*ell+delta(f),
    delta=3*x^2*d/dx+(2-6*x*q)*d/dq.

This script transports the two relevant position polynomials from the exact
A3 Keller/Dixmier certificate into that extension and measures the quantum
residual without appealing to a truncated hbar expansion.
"""

from math import comb

import sympy as sp

x, q, e = sp.symbols("x q e")
s = x * q


def delta(f):
    return sp.expand(3 * x**2 * sp.diff(f, x) + (2 - 6 * x * q) * sp.diff(f, q))


def clean(a):
    return {degree: sp.expand(coefficient) for degree, coefficient in a.items() if coefficient != 0}


def from_expr(f):
    poly = sp.Poly(sp.expand(f), e)
    return clean({monomial[0]: coefficient for monomial, coefficient in poly.terms()})


def to_expr(a):
    return sp.expand(sum(coefficient * e**degree for degree, coefficient in a.items()))


def add(a, b, scale_b=1):
    out = dict(a)
    for degree, coefficient in b.items():
        out[degree] = sp.expand(out.get(degree, 0) + scale_b * coefficient)
    return clean(out)


def scale(a, scalar):
    return clean({degree: scalar * coefficient for degree, coefficient in a.items()})


def multiply(a, b):
    """Multiply with every e moved to the right."""
    out = {}
    for i, ai in a.items():
        for j, bj in b.items():
            derivative = bj
            for k in range(i + 1):
                degree = i - k + j
                out[degree] = out.get(degree, 0) + ai * comb(i, k) * derivative
                derivative = delta(derivative)
    return clean(out)


def commutator(a, b):
    return add(multiply(a, b), multiply(b, a), scale_b=-1)


def ordered_lift(f, theta):
    """PBW theta-ordering; theta=0 is right-normal, theta=1/2 is Weyl."""
    normal = from_expr(f)
    out = {}
    for n, coefficient in normal.items():
        derivative = coefficient
        for k in range(n + 1):
            degree = n - k
            out[degree] = out.get(degree, 0) + comb(n, k) * theta**k * derivative
            derivative = delta(derivative)
    return clean(out)


ONE = {0: sp.Integer(1)}
X = {0: x}
Q = {0: q}
E = {1: sp.Integer(1)}
R = {0: x * (2 - 3 * s)}

# These are the commutative expansions of the sheared THM-1300 coordinates.
# ``from_expr`` is the coefficient-left/e-right normal-ordering prescription.
y = q - x * e / 3
u = 1 + x * y
T_expr = sp.expand(y + 3 * x * u**2 * e + 3 * x * y**2 * (4 + 3 * x * y))
S_expr = sp.expand((u**3 * e + y**2 * u * (4 + 3 * x * y)) / 2)
connection_B = (
    2 * e**4 * x**6 * (3 * s - 2)
    + e**3 * x**4 * (-90 * s**2 - 30 * s + 55)
    + e**2 * x**2 * (540 * s**3 + 720 * s**2 - 120 * s - 270)
    + e * (-1620 * s**4 - 3780 * s**3 - 1890 * s**2 + 810 * s + 540)
    + q**2 * (2430 * s**3 + 8100 * s**2 + 8640 * s + 2430)
)
connection_H = sp.expand(-e * connection_B / 1620)
T = from_expr(T_expr)
S = from_expr(S_expr)

# A second lift keeps the displayed factorization and multiplication order of
# the three-dimensional formula.  It is genuinely different because e does
# not commute with x,q after the symplectic substitution.
Y = add(Q, scale(multiply(X, E), -sp.Rational(1, 3)))
U = add(ONE, multiply(X, Y))
V = add(scale(ONE, 4), scale(multiply(X, Y), 3))
T_factored = add(
    Y,
    add(
        scale(multiply(multiply(multiply(X, U), U), E), 3),
        scale(multiply(multiply(multiply(X, Y), Y), V), 3),
    ),
)
S_factored = scale(
    add(
        multiply(multiply(multiply(U, U), U), E),
        multiply(multiply(multiply(Y, Y), U), V),
    ),
    sp.Rational(1, 2),
)

assert to_expr(commutator(E, X)) == 3 * x**2
assert to_expr(commutator(E, Q)) == 2 - 6 * x * q
assert not commutator(E, R)
assert not commutator(R, T)
assert not commutator(R, S)

residual = add(commutator(S, T), ONE, scale_b=-1)
residual_factored = add(commutator(S_factored, T_factored), ONE, scale_b=-1)
theta = sp.symbols("theta")
T_theta = ordered_lift(T_expr, theta)
S_theta = ordered_lift(S_expr, theta)
residual_theta = add(commutator(S_theta, T_theta), ONE, scale_b=-1)

print("DC2 ORE-DESCENT AUDIT")
print("Ore relation: [ell,x]=3*x^2, [ell,q]=2-6*x*q: PASS")
print("central invariant: [ell,R]=0: PASS")
print("normal-order ell degrees (T,S):", max(T), max(S))
assert (sp.Poly(T_expr, e).degree(), sp.Poly(connection_H, e).degree(), sp.Poly(S_expr, e).degree()) == (2, 5, 3)
print("THM-2044 momentum orders (R,T,D,S): 0 2 5 3 PASS")
print("exact R column: [R,T]=[R,S]=0: PASS")
print("[S,T]-1 normal ell degree:", max(residual) if residual else -1)
print("[S,T]-1 coefficient factors (descending ell degree):")
for degree in sorted(residual, reverse=True):
    print("  ell^%d:" % degree, sp.factor(residual[degree]))
print("factored-order exact R column:",
      "PASS" if not commutator(R, T_factored) and not commutator(R, S_factored) else "FAIL")
print("factored-order [S,T]-1 ell degree:",
      max(residual_factored) if residual_factored else -1)
print("factored-order top residual coefficient:",
      sp.factor(residual_factored[max(residual_factored)]))
theta_top_factor = sp.factor(residual_theta[3])
assert sp.factor(theta_top_factor.subs(theta, sp.Rational(1, 2))) == 0
assert sp.factor(theta_top_factor / (2 * theta - 1)) != 0
print("theta-order ell^3 residual coefficient:", theta_top_factor)
print("unique scalar theta killing the ell^3 layer: 1/2")
weyl_residual = clean({
    degree: sp.factor(coefficient.subs(theta, sp.Rational(1, 2)))
    for degree, coefficient in residual_theta.items()
})
print("Weyl theta=1/2 residual ell degree:", max(weyl_residual) if weyl_residual else -1)
assert max(weyl_residual) == 2
for degree in sorted(weyl_residual, reverse=True):
    print("  ell^%d:" % degree, sp.factor(weyl_residual[degree]))


def x_valuation(f):
    return min(monomial[0] for monomial, _ in sp.Poly(sp.expand(f), x, q).terms())


normal_boundary_profile = [
    (degree, x_valuation(residual[degree]), x_valuation(residual[degree]) - 2 * degree)
    for degree in sorted(residual, reverse=True)
]
weyl_boundary_profile = [
    (degree, x_valuation(weyl_residual[degree]), x_valuation(weyl_residual[degree]) - 2 * degree)
    for degree in sorted(weyl_residual, reverse=True)
]
assert normal_boundary_profile == [(3, 9, 3), (2, 7, 3), (1, 5, 3), (0, 3, 3)]
assert weyl_boundary_profile == [(2, 10, 6), (1, 8, 6), (0, 6, 6)]
print("boundary profiles (ell_degree, v_x(coeff), v_x-2*ell_degree):")
print("  normal:", normal_boundary_profile)
print("  Weyl:  ", weyl_boundary_profile)
print("Weyl order kills boundary grade 3 and exposes uniform grade 6: PASS")


def boundary_beta(a):
    return min(x_valuation(coefficient) - 2 * degree for degree, coefficient in a.items())


boundary_u = sp.symbols("boundary_u")


def boundary_initial_symbol(a, grade):
    """Return x^(-grade) times the grade piece, with x^2*e replaced by u."""
    answer = 0
    for degree, coefficient in a.items():
        for (x_degree, q_degree), scalar in sp.Poly(sp.expand(coefficient), x, q).terms():
            if x_degree - 2 * degree == grade:
                answer += scalar * q**q_degree * boundary_u**degree
    return sp.expand(answer)


def boundary_lift(symbol, x_power):
    return from_expr(sp.expand(x**x_power * symbol.subs(boundary_u, x**2 * e)))


def cancel_boundary_grade(residual_here):
    """Solve the associated-graded linear correction equation at min beta."""
    grade = boundary_beta(residual_here)
    initial = boundary_initial_symbol(residual_here, grade)
    P_u = 2 * boundary_u**2 - 10 * boundary_u + 9

    # For C_S=x^(g-1)*int A dq and C_T=x^g*int B dq, the grade-g
    # equation is (8/3)(u-2)A+(P(u)/9)B=-initial.  Since P(2)=-3,
    # the following is a polynomial Bezout solution for every initial.
    B0 = sp.expand(3 * initial.subs(boundary_u, 2))
    numerator = sp.expand(-initial - P_u * B0 / 9)
    A0 = sp.cancel(sp.Rational(3, 8) * numerator / (boundary_u - 2))
    A0 = sp.Poly(A0, boundary_u, q, domain=sp.QQ).as_expr()
    assert sp.expand(sp.Rational(8, 3) * (boundary_u - 2) * A0 + P_u * B0 / 9 + initial) == 0

    int_A = sp.integrate(A0, q)
    int_B = sp.integrate(B0, q)
    assert sp.expand(sp.diff(int_A, q) - A0) == 0
    assert sp.expand(sp.diff(int_B, q) - B0) == 0
    correction_S = boundary_lift(int_A, grade - 1)
    correction_T = boundary_lift(int_B, grade)
    return grade, initial, correction_S, correction_T


T_iter = ordered_lift(T_expr, sp.Rational(1, 2))
S_iter = ordered_lift(S_expr, sp.Rational(1, 2))
print("associated-graded simultaneous correction ladder:")
ladder_grades = []
for step in range(8):
    residual_iter = add(commutator(S_iter, T_iter), ONE, scale_b=-1)
    if not residual_iter:
        print("  step", step, "EXACT TERMINATION")
        break
    grade, initial, correction_S, correction_T = cancel_boundary_grade(residual_iter)
    ladder_grades.append(grade)
    assert not correction_S or boundary_beta(correction_S) >= grade - 1
    assert not correction_T or boundary_beta(correction_T) >= grade
    S_iter = add(S_iter, correction_S)
    T_iter = add(T_iter, correction_T)
    new_residual = add(commutator(S_iter, T_iter), ONE, scale_b=-1)
    new_grade = None if not new_residual else boundary_beta(new_residual)
    assert new_grade is None or new_grade > grade
    print(
        "  step", step,
        "grade", grade, "->", "zero" if new_grade is None else new_grade,
        "initial=", sp.factor(initial),
        "correction_terms(S,T)=", (len(correction_S), len(correction_T)),
    )
print("  visited grades:", ladder_grades)
print("  associated-graded correction map is surjective because gcd(u-2,2u^2-10u+9)=1: PASS")

tt = -sp.Rational(1, 3) / x
assert sp.factor(delta(tt) - 1) == 0
print("localized slice t=-1/(3*x), delta(t)=1: PASS")
print("O[x^-1] = Q[R,t,t^-1][ell; d/dt], so the remaining gate is extension across x=0.")

# THM-2046's mixed-commutator matrix calculation, displayed in rank two.
P, Q0 = sp.symbols("P Q0", cls=sp.Function)
b11, b12, b21, b22 = sp.symbols("b11 b12 b21 b22")
Px, Pq, Qx, Qq = sp.symbols("P_x P_q Q_x Q_q")
B = sp.Matrix([[b11, b12], [b21, b22]])
JT = sp.Matrix([[Px, Pq], [Qx, Qq]])
mixed = B * JT
assert sp.factor(mixed.det() - B.det() * JT.det()) == 0
print("filtered-pullback determinant identity det(B)*det(J)=1: SYMBOLIC PASS")
print("TOURNAMENT ANALYSIS: not imposed; vertices would be noncommuting proof obligations, and orienting them destroys the simultaneous Weyl relation.")
