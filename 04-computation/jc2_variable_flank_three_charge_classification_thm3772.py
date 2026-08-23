#!/usr/bin/env python3
"""Exact companion for THM-3772's variable-flank classification."""

from __future__ import annotations

import ast
import hashlib
from itertools import product
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


X, T, Z, L = sp.symbols("X T Z L")


def jac(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, X) * sp.diff(second, T)
        - sp.diff(first, T) * sp.diff(second, X)
    )


def component(a_profile: sp.Expr, chi: sp.Expr, b_profile: sp.Expr) -> sp.Expr:
    z_source = X * T
    return sp.expand(
        X * a_profile.subs(Z, z_source)
        + chi.subs(Z, z_source)
        + T * b_profile.subs(Z, z_source)
    )


def direct_smooth(q_value: sp.Expr) -> bool:
    critical = sp.groebner(
        (sp.diff(q_value, X), sp.diff(q_value, T)),
        X,
        T,
        order="grevlex",
    )
    return critical.contains(sp.Integer(1))


def pure_profile_smooth(profile: sp.Expr) -> bool:
    if sp.expand(profile) == 0 or sp.expand(profile.subs(Z, 0)) == 0:
        return False
    polynomial = sp.Poly(profile, Z, domain=sp.QQ)
    return sp.gcd(polynomial, polynomial.diff()).degree() == 0


z_source = X * T

# Universal quadratic generic-fibre identities with algebraically independent
# coefficients.  They verify the sign before any special profile is imposed.
a_coefficients = sp.symbols("A_0:5")
b_coefficients = sp.symbols("B_0:5")
chi_coefficients = sp.symbols("chi_0:5")
a_generic = sum(a_coefficients[i] * Z**i for i in range(5))
b_generic = sum(b_coefficients[i] * Z**i for i in range(5))
chi_generic = sum(chi_coefficients[i] * Z**i for i in range(5))
f_generic = sp.expand(Z * a_generic * b_generic)
q_generic = component(a_generic, chi_generic, b_generic)
y_generic = sp.expand(
    X * a_generic.subs(Z, z_source)
    - T * b_generic.subs(Z, z_source)
)
delta_generic = sp.expand((L - chi_generic) ** 2 - 4 * f_generic)
gate(sp.expand(jac(q_generic, z_source) - y_generic) == 0,
     "universal Jacobian coordinate")
gate(
    sp.expand(
        y_generic**2
        - delta_generic.subs({L: q_generic, Z: z_source})
    )
    == 0,
    "universal generic-fibre equation",
)

# The repeated-root eliminant of Delta is the same F'^2-F chi'^2 object as
# in the normalized theorem, despite both flanks now varying.
g_generic = sp.diff(chi_generic, Z)
r_generic = sp.expand(
    sp.diff(f_generic, Z) ** 2 - f_generic * g_generic**2
)
w = sp.symbols("W")
delta_w = sp.expand((w - chi_generic) ** 2 - 4 * f_generic)
delta_w_prime = sp.diff(delta_w, Z)
gate(
    sp.expand(
        sp.resultant(delta_w, delta_w_prime, w)
        - 16 * r_generic
    )
    == 0,
    "generic repeated-root eliminant",
)

# Mixed rational-exact boundary, first orientation: A is constant and B is
# affine.  The second orientation is its literal X<->T dual.
a, b, g, h, p, c = sp.symbols("a b g h p c", nonzero=True)
chi_boundary = h + g * Z
product_boundary = p + g**2 * Z / 4

a_first = sp.sympify(a)
b_first = sp.expand(product_boundary / a)
q_first = component(a_first, chi_boundary, b_first)
y_first = sp.expand(
    X * a_first - T * b_first.subs(Z, z_source)
)
d_first = sp.expand(g * (q_first - h) + 2 * p)
p_first = sp.cancel(c * y_first / d_first)
gate(
    sp.expand(
        q_first
        - h
        - p * T / a
        - a * X * (1 + g * T / (2 * a)) ** 2
    )
    == 0,
    "constant-A boundary form",
)
gate(sp.factor(jac(p_first, q_first) - c) == 0,
     "constant-A rational primitive")

b_second = sp.sympify(b)
a_second = sp.expand(product_boundary / b)
q_second = component(a_second, chi_boundary, b_second)
y_second = sp.expand(
    X * a_second.subs(Z, z_source) - T * b_second
)
d_second = sp.expand(g * (q_second - h) + 2 * p)
p_second = sp.cancel(c * y_second / d_second)
gate(
    sp.expand(
        q_second
        - h
        - p * X / b
        - b * T * (1 + g * X / (2 * b)) ** 2
    )
    == 0,
    "constant-B boundary form",
)
gate(sp.factor(jac(p_second, q_second) - c) == 0,
     "constant-B rational primitive")

# Sharp smoothness and opposite principal parts on the nonlinear mixed wall.
s_first = sp.expand(1 + g * T / (2 * a))
residual_first = sp.expand(2 * p / g + a * X * s_first)
lambda_zero = sp.expand(h - 2 * p / g)
gate(sp.expand(q_first - lambda_zero - s_first * residual_first) == 0,
     "constant-A split vertical fibre")
w_zero = sp.expand(lambda_zero - h)
gate(sp.expand((y_first + w_zero).subs(T, -2 * a / g)) == 0,
     "constant-A first component Y value")
x_on_residual = sp.cancel(-2 * p / (g * a * s_first))
gate(
    sp.factor((y_first - w_zero).subs(X, x_on_residual)) == 0,
    "constant-A residual component Y value",
)
gate(
    sp.expand(sp.diff(q_first, X).subs(T, -2 * a / g)) == 0
    and sp.expand(
        sp.diff(q_first, T).subs({T: -2 * a / g, p: 0})
    )
    == 0,
    "constant-A p-zero singular boundary",
)

s_second = sp.expand(1 + g * X / (2 * b))
residual_second = sp.expand(2 * p / g + b * T * s_second)
gate(sp.expand(q_second - lambda_zero - s_second * residual_second) == 0,
     "constant-B split vertical fibre")
gate(sp.expand((y_second - w_zero).subs(X, -2 * b / g)) == 0,
     "constant-B first component Y value")
t_on_residual = sp.cancel(-2 * p / (g * b * s_second))
gate(
    sp.factor((y_second + w_zero).subs(T, t_on_residual)) == 0,
    "constant-B residual component Y value",
)
gate(
    sp.expand(sp.diff(q_second, T).subs(X, -2 * b / g)) == 0
    and sp.expand(
        sp.diff(q_second, X).subs({X: -2 * b / g, p: 0})
    )
    == 0,
    "constant-B p-zero singular boundary",
)

# At g=0 the product is a nonzero constant and both orientations collapse to
# a linear form with a polynomial mate.
q_linear = h + a * X + b * T
p_linear = c * (a * X - b * T) / (2 * a * b)
gate(sp.expand(jac(p_linear, q_linear) - c) == 0,
     "mixed linear endpoint")

# Pure-flank exactness and signs.  Nonconstant profiles remain rationally
# exact but their vertical poles cannot be removed by a target-only shear.
profile_coefficients = sp.symbols("profile_0:6")
profile = sum(profile_coefficients[i] * Z**i for i in range(6))
q_pure_a = component(profile, sp.sympify(h), sp.Integer(0))
p_pure_a = sp.cancel(-c * T / profile.subs(Z, z_source))
q_pure_b = component(sp.Integer(0), sp.sympify(h), profile)
p_pure_b = sp.cancel(c * X / profile.subs(Z, z_source))
gate(sp.factor(jac(p_pure_a, q_pure_a) - c) == 0,
     "pure positive-flank rational primitive")
gate(sp.factor(jac(p_pure_b, q_pure_b) - c) == 0,
     "pure negative-flank rational primitive")

pure_smooth_controls = 0
pure_singular_controls = 0
for degree in range(1, 8):
    good_profile = 1 + Z**degree
    for q_value in (
        component(good_profile, sp.Integer(0), sp.Integer(0)),
        component(sp.Integer(0), sp.Integer(0), good_profile),
    ):
        gate(direct_smooth(q_value), "pure squarefree smooth control")
        pure_smooth_controls += 1
for bad_profile in (Z, (1 + Z) ** 2, Z * (1 + Z), (1 - Z**2) ** 2):
    for q_value in (
        component(bad_profile, sp.Integer(0), sp.Integer(0)),
        component(sp.Integer(0), sp.Integer(0), bad_profile),
    ):
        gate(not direct_smooth(q_value), "pure singular boundary control")
        pure_singular_controls += 1

# Exhaust all linear A,B,chi profiles over {-1,0,1}.  On every directly
# smooth mixed member the generic Delta is squarefree; exactness occurs
# exactly at AB=p+g^2 z/4.  Pure members agree with the separate classification.
values = (-1, 0, 1)
profiles = [u + v * Z for u, v in product(values, repeat=2)]
grid_total = 0
grid_smooth = 0
grid_rational_exact = 0
mixed_squarefree_controls = 0
lambda_domain = sp.QQ.frac_field(L)
for a_profile in profiles:
    for b_profile in profiles:
        for chi in profiles:
            q_value = component(a_profile, chi, b_profile)
            smooth = direct_smooth(q_value)
            grid_total += 1
            grid_smooth += int(smooth)
            if not smooth:
                continue

            a_zero = sp.expand(a_profile) == 0
            b_zero = sp.expand(b_profile) == 0
            if a_zero or b_zero:
                if a_zero and b_zero:
                    predicted_exact = False
                else:
                    predicted_exact = sp.degree(chi, Z) <= 0
                    active = b_profile if a_zero else a_profile
                    gate(pure_profile_smooth(active),
                         "smooth pure profile criterion")
            else:
                f_value = sp.expand(Z * a_profile * b_profile)
                delta_value = sp.expand((L - chi) ** 2 - 4 * f_value)
                delta_poly = sp.Poly(delta_value, Z, domain=lambda_domain)
                gate(
                    sp.gcd(delta_poly, delta_poly.diff()).degree() == 0,
                    "smooth mixed generic squarefreeness",
                )
                mixed_squarefree_controls += 1
                g_value = sp.diff(chi, Z)
                p_value = sp.expand((a_profile * b_profile).subs(Z, 0))
                predicted_exact = sp.expand(
                    a_profile * b_profile - p_value - g_value**2 * Z / 4
                ) == 0
                gate(
                    predicted_exact == (delta_poly.degree() <= 1),
                    "mixed degree-one exactness boundary",
                )

            grid_rational_exact += int(bool(predicted_exact))

# Finite polynomial-mate hostiles on one representative from each nonlinear
# rational-exact family supplement the all-degree divisor argument.
mixed_sample = sp.expand(q_first.subs({a: 1, g: 2, h: 0, p: 1}))
pure_sample = component(1 + Z, sp.Integer(0), sp.Integer(0))
mate_obstructions = 0
for sample_index, q_value in enumerate((mixed_sample, pure_sample)):
    for bound in range(0, 10):
        basis = [
            X**i * T**j
            for i in range(bound + 1)
            for j in range(bound + 1 - i)
        ]
        coefficients = sp.symbols(
            f"mate_{sample_index}_{bound}_0:{len(basis)}"
        )
        prospective = sum(
            coefficient * monomial
            for coefficient, monomial in zip(coefficients, basis)
        )
        equation = sp.Poly(jac(prospective, q_value) - 1, X, T)
        system = [coefficient for _, coefficient in equation.terms()]
        gate(sp.linsolve(system, coefficients) == sp.EmptySet,
             "finite polynomial-mate obstruction")
        mate_obstructions += 1

semantic_rows = (
    "component:Q=X*A(XT)+chi(XT)+T*B(XT)",
    "generic:Y=X*A-T*B;Y^2=(Lambda-chi)^2-4zAB;J(Q,z)=Y",
    "mixed-rational:chi=h+gz;AB=p+g^2*z/4;p!=0",
    "mixed-polynomial:only-g=0-linear;nonlinear-split-fibre-opposite-values",
    "pure-rational:chi-constant;active-profile-squarefree-origin-unit",
    "pure-polynomial:active-profile-constant-only",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3772-variable-flank-three-charge-rational-exact-classification")
print("scope=algebraically_closed_characteristic_zero;arbitrary_A,B,chi_in_k[z]")
print("mixed_rational_mate=iff_chi=h+gz_and_AB=p+g^2*z/4_with_p_nonzero")
print("pure_rational_mate=iff_chi_constant_and_active_flank_smooth")
print("polynomial_mate=iff_Q_is_linear")
print("mechanisms=quadratic_time_form;factorization;vertical_principal_parts")
print(
    f"linear_profile_grid={grid_total};smooth={grid_smooth};"
    f"rational_exact={grid_rational_exact};"
    f"mixed_squarefree={mixed_squarefree_controls}"
)
print(
    f"pure_smooth_controls={pure_smooth_controls};"
    f"pure_singular_controls={pure_singular_controls};"
    f"mate_obstructions={mate_obstructions}"
)
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
