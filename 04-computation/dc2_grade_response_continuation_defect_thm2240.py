#!/usr/bin/env python3
"""Exact Ore audit for THM-2240.

This is intentionally standalone: it rebuilds the THM-2049 Weyl pair and the
coefficient-left Ore product instead of importing the older descent script.
All checks use exact SymPy arithmetic and explicit ``require`` calls so that
normal and ``python -O`` executions have identical logical content.
"""

from math import comb

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


x, q, ell = sp.symbols("x q ell")
u = sp.symbols("u")


def delta(coefficient):
    return sp.expand(
        3 * x**2 * sp.diff(coefficient, x)
        + (2 - 6 * x * q) * sp.diff(coefficient, q)
    )


def clean(operator):
    return {
        degree: sp.expand(coefficient)
        for degree, coefficient in operator.items()
        if sp.expand(coefficient) != 0
    }


def from_commutative_expr(expression):
    """Read ell as a commuting placeholder, then put it coefficient-left."""
    polynomial = sp.Poly(sp.expand(expression), ell)
    return clean(
        {monomial[0]: coefficient for monomial, coefficient in polynomial.terms()}
    )


def add(left, right, scale_right=1):
    answer = dict(left)
    for degree, coefficient in right.items():
        answer[degree] = sp.expand(
            answer.get(degree, 0) + scale_right * coefficient
        )
    return clean(answer)


def ore_multiply(left, right):
    """Multiply with ell*f=f*ell+delta(f), returning coefficient-left form."""
    answer = {}
    for i, a_i in left.items():
        for j, b_j in right.items():
            derivative = b_j
            for r in range(i + 1):
                degree = i - r + j
                answer[degree] = sp.expand(
                    answer.get(degree, 0)
                    + a_i * comb(i, r) * derivative
                )
                derivative = delta(derivative)
    return clean(answer)


def commutator(left, right):
    return add(ore_multiply(left, right), ore_multiply(right, left), -1)


def weyl_order(expression):
    """The theta=1/2 PBW ordering used in THM-2049."""
    normal = from_commutative_expr(expression)
    answer = {}
    for degree, coefficient in normal.items():
        derivative = coefficient
        for r in range(degree + 1):
            target_degree = degree - r
            answer[target_degree] = sp.expand(
                answer.get(target_degree, 0)
                + comb(degree, r) * sp.Rational(1, 2) ** r * derivative
            )
            derivative = delta(derivative)
    return clean(answer)


def x_valuation(coefficient):
    return min(
        monomial[0]
        for monomial, _ in sp.Poly(sp.expand(coefficient), x, q).terms()
    )


def beta(operator):
    return min(
        x_valuation(coefficient) - 2 * degree
        for degree, coefficient in operator.items()
    )


def homogeneous_symbol(operator, grade):
    """Return H_g with in_g(operator)=x^g H_g(q,u) in Rees notation."""
    answer = 0
    for degree, coefficient in operator.items():
        for (x_degree, q_degree), scalar in sp.Poly(
            sp.expand(coefficient), x, q
        ).terms():
            if x_degree - 2 * degree == grade:
                answer += scalar * q**q_degree * u**degree
    return sp.expand(answer)


def symbol_lift(symbol, x_power):
    """Map sum f_k(q)u^k to sum f_k(q)x^(x_power+2k)ell^k.

    This is a coefficient-left symbol lift. It is not literal evaluation of
    a polynomial at the noncommutative product x^2*ell.
    """
    return from_commutative_expr(
        sp.expand(x**x_power * symbol.subs(u, x**2 * ell))
    )


ONE = {0: sp.Integer(1)}

# Rebuild the commutative sheared positions and then Weyl-order them.
y = q - x * ell / 3
w = 1 + x * y
T_expr = sp.expand(y + 3 * x * w**2 * ell + 3 * x * y**2 * (4 + 3 * x * y))
S_expr = sp.expand((w**3 * ell + y**2 * w * (4 + 3 * x * y)) / 2)
T_W = weyl_order(T_expr)
S_W = weyl_order(S_expr)


def residual(S_operator, T_operator):
    return add(commutator(S_operator, T_operator), ONE, -1)


D = 2 * u**2 - 10 * u + 9
K = 4 * u**2 - 13 * u + 13

R_W = residual(S_W, T_W)
require(beta(R_W) == 6, "Weyl residual must begin in beta grade six")
require(
    sp.expand(homogeneous_symbol(R_W, 6) - 2 * (u - 1) * (u - 2)) == 0,
    "unexpected Weyl grade-six residual",
)

# The standard grade-six section of THM-2049.
A_standard = -sp.Rational(3, 4) * q * (u - 1)
B_standard = sp.Integer(0)
S_standard = add(S_W, symbol_lift(A_standard, 5))
T_standard = T_W
R_standard = residual(S_standard, T_standard)
H7_standard = -sp.Rational(3, 2) * q * (10 * u**2 - 36 * u + 29)
require(beta(R_standard) == 7, "standard correction must expose grade seven")
require(
    sp.expand(homogeneous_symbol(R_standard, 7) - H7_standard) == 0,
    "unexpected standard grade-seven residual",
)

# The C=1 syzygy gauge, with zero q-integration constants.
eta_S = symbol_lift(q * D, 5)
eta_T = symbol_lift(-24 * q * (u - 2), 6)
S_C1 = add(S_standard, eta_S)
T_C1 = add(T_standard, eta_T)
R_C1 = residual(S_C1, T_C1)
H7_C1 = q * (34 * u**2 - 100 * u + 121) / 2
C1_splitter = 8 * q * K

piece_S = commutator(eta_S, T_W)
piece_T = commutator(S_W, eta_T)
require(
    sp.expand(
        homogeneous_symbol(piece_S, 6) - sp.Rational(8, 3) * (u - 2) * D
    )
    == 0,
    "C=1 S-piece has wrong grade-six response",
)
require(
    sp.expand(
        homogeneous_symbol(piece_T, 6) + sp.Rational(8, 3) * (u - 2) * D
    )
    == 0,
    "C=1 T-piece has wrong grade-six response",
)
require(
    sp.expand(
        homogeneous_symbol(piece_S, 7)
        - 8 * q * (2 * u - 5) * (2 * u**2 - 8 * u + 5)
    )
    == 0,
    "C=1 S-piece has wrong grade-seven response",
)
require(
    sp.expand(
        homogeneous_symbol(piece_T, 7)
        + 8 * q * (u - 2) * (4 * u**2 - 22 * u + 19)
    )
    == 0,
    "C=1 T-piece has wrong grade-seven response",
)
require(beta(R_C1) == 7, "C=1 correction must expose grade seven")
require(homogeneous_symbol(R_C1, 6) == 0, "C=1 failed to kill grade six")
require(
    sp.expand(homogeneous_symbol(R_C1, 7) - H7_C1) == 0,
    "unexpected C=1 grade-seven residual",
)
require(
    sp.expand(H7_C1 - H7_standard - C1_splitter) == 0,
    "C=1 splitter identity failed",
)

# Derive the full first-subleading response using the first two Rees layers.
t = sp.Rational(2, 3) * u * (4 - u)
v = (u - 1) ** 2
s = u * (u - 3) * (2 * u - 9) / 54
r = -u * (u**2 - 3 * u - 3) / 18

F_fun = sp.Function("F")(q, u)
G_fun = sp.Function("G")(q, u)
M_S = (
    (2 * v - 3 * t) * sp.diff(F_fun, u)
    - 15 * sp.diff(t, u) * F_fun
    + q * (6 * sp.diff(t, u) - 2 * sp.diff(v, u)) * sp.diff(F_fun, q)
)
M_T = (
    18 * sp.diff(s, u) * G_fun
    + (6 * s - 2 * r) * sp.diff(G_fun, u)
    + q * (2 * sp.diff(r, u) - 6 * sp.diff(s, u)) * sp.diff(G_fun, q)
)

J_fun = sp.Function("J")(q, u)
general_syzygy_response = sp.expand(
    M_S.subs(F_fun, D * J_fun).doit()
    + M_T.subs(G_fun, -24 * (u - 2) * J_fun).doit()
)
expected_syzygy_response = sp.expand(
    4 * K * (J_fun + q * sp.diff(J_fun, q)) + 18 * sp.diff(J_fun, u)
)
require(
    sp.expand(general_syzygy_response - expected_syzygy_response) == 0,
    "arbitrary-q syzygy response identity failed",
)


def L_a(a):
    return sp.expand(
        20 * (u - 2) * a
        + 2 * (2 * u**2 - 6 * u + 1) * sp.diff(a, u)
    )


def L_b(b):
    return sp.expand(
        D * b + u * (u - 2) * (u - 4) * sp.diff(b, u) / 3
    )


def L_J(J):
    return sp.expand(4 * K * (J + q * sp.diff(J, q)) + 18 * sp.diff(J, u))


# Axiswise leading-term audits, including the constant edge n=0.
for n in range(9):
    image_a = sp.Poly(L_a(u**n), u)
    require(image_a.degree() == n + 1, f"L_a degree failed at n={n}")
    require(image_a.LC() == 4 * (n + 5), f"L_a LC failed at n={n}")

    image_b = sp.Poly(L_b(u**n), u)
    require(image_b.degree() == n + 2, f"L_b degree failed at n={n}")
    require(
        image_b.LC() == sp.Rational(n + 6, 3),
        f"L_b LC failed at n={n}",
    )

for q_degree in range(6):
    for u_degree in range(7):
        c = u**u_degree
        layer = sp.expand(
            (
                4 * (q_degree + 2) * K * c
                + 18 * sp.diff(c, u)
            )
            / (q_degree + 1)
        )
        polynomial = sp.Poly(layer, u)
        require(
            polynomial.degree() == u_degree + 2,
            "q-layer syzygy degree audit failed",
        )
        require(
            polynomial.LC()
            == sp.Rational(16 * (q_degree + 2), q_degree + 1),
            "q-layer syzygy leading coefficient audit failed",
        )

# Referee correction: the sum of the two integration-constant axes is not
# injective, even though each axis separately is.
a_kernel = -146 + 140 * u - 25 * u**2
b_kernel = -680 + 300 * u
require(a_kernel != 0 and b_kernel != 0, "kernel witness must be nonzero")
require(
    sp.expand(L_a(a_kernel) + L_b(b_kernel)) == 0,
    "combined (a,b) kernel witness failed",
)

# Full-Ore hostile cross-check with mixed q/u degrees and mixed signs.
C_test = (
    2
    - 3 * u
    + 5 * u**2
    + q * (-1 + 2 * u**3)
    + q**2 * (7 - 4 * u + u**4)
)
J_test = sp.integrate(C_test, q)
require(sp.expand(J_test.subs(q, 0)) == 0, "J_test must have zero q-constant")
a_test = -1 + 4 * u + 2 * u**3
b_test = 3 - u + 6 * u**2 + u**5
A_test = sp.expand(A_standard + D * J_test + a_test)
B_test = sp.expand(-24 * (u - 2) * J_test + b_test)
S_test = add(S_W, symbol_lift(A_test, 5))
T_test = add(T_W, symbol_lift(B_test, 6))
R_test = residual(S_test, T_test)
require(homogeneous_symbol(R_test, 6) == 0, "hostile gauge failed at grade six")
require(
    sp.expand(
        homogeneous_symbol(R_test, 7)
        - H7_standard
        - L_J(J_test)
        - L_a(a_test)
        - L_b(b_test)
    )
    == 0,
    "full-Ore hostile grade-seven formula failed",
)

# A common next-rung correction shifts both residuals equally at grade seven;
# it cannot erase the splitter without first knowing which representative was
# chosen.
next_S = symbol_lift(q * (u + 1), 6)
next_T = symbol_lift(q**2 + u, 7)
R_standard_next = residual(add(S_standard, next_S), add(T_standard, next_T))
R_C1_next = residual(add(S_C1, next_S), add(T_C1, next_T))
require(
    sp.expand(
        homogeneous_symbol(R_C1_next, 7)
        - homogeneous_symbol(R_standard_next, 7)
        - C1_splitter
    )
    == 0,
    "common next-rung correction erased the splitter",
)

# Pure-rung no-skip: higher q layers are injective and the q^1 layer cannot
# cancel the standard residual.  The u^2 coefficient would force c=15/32,
# after which the u coefficient remains 21/4.
c_from_u2 = sp.Rational(15, 32)
remaining_u_coefficient = sp.expand(54 - 104 * c_from_u2)
require(
    remaining_u_coefficient == sp.Rational(21, 4),
    "constant c incompatibility calculation failed",
)

print("THM-2240 DC2 GRADE-RESPONSE CONTINUATION DEFECT")
print("ambient=Q[x,q][ell;delta], delta=3*x^2*d_x+(2-6*x*q)*d_q")
print("symbol_lift=sum f_k(q)u^k -> sum f_k(q)x^(p+2k)ell^k")
print("weyl_grade6=", sp.factor(homogeneous_symbol(R_W, 6)))
print("standard_grade7=", sp.factor(H7_standard))
print("C1_grade6_piece_S=", sp.factor(homogeneous_symbol(piece_S, 6)))
print("C1_grade6_piece_T=", sp.factor(homogeneous_symbol(piece_T, 6)))
print("C1_grade7_piece_S=", sp.factor(homogeneous_symbol(piece_S, 7)))
print("C1_grade7_piece_T=", sp.factor(homogeneous_symbol(piece_T, 7)))
print("C1_grade7=", sp.factor(H7_C1))
print("C1_splitter=", sp.factor(C1_splitter))
print("arbitrary_q_syzygy_response=4*K*(J+q*J_q)+18*J_u")
print("q_layer_response=[4*(r+2)*K*c_r+18*c_r']/(r+1)")
print("L_a=20*(u-2)*a+2*(2*u^2-6*u+1)*a'")
print("L_b=D*b+u*(u-2)*(u-4)*b'/3")
print("L_a_axis_leading=4*(n+5)*lc(a), degree=n+1")
print("L_b_axis_leading=(n+6)*lc(b)/3, degree=n+2")
print("q_layer_numerator_leading=16*(r+2)*lc(c_r), degree=n+2")
print("axiswise_injectivity_checks=n:0..8; r:0..5; u_degree:0..6 PASS")
print("combined_ab_kernel_a=", sp.factor(a_kernel))
print("combined_ab_kernel_b=", sp.factor(b_kernel))
print("combined_ab_kernel_response=", sp.expand(L_a(a_kernel) + L_b(b_kernel)))
print("hostile_full_ore_arbitrary_q_crosscheck=PASS")
print("common_next_rung_splitter=", sp.factor(C1_splitter))
print("constant_c_from_u2=", c_from_u2)
print("remaining_u_coefficient=", remaining_u_coefficient)
print("pure_homogeneous_grade6_rung_cannot_skip_grade7=PASS")
print("scope=next-rung beta6/beta7 attachments excluded; formal recursion unchanged")
print("all_exact_checks=PASS")
