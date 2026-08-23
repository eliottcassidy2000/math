#!/usr/bin/env python3
"""Exact companion for THM-3765's normalized three-charge classification."""

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


def source_component(chi: sp.Expr, psi: sp.Expr) -> sp.Expr:
    return sp.expand(X + chi.subs(Z, X * T) + T * psi.subs(Z, X * T))


def predicted_smooth(chi: sp.Expr, psi: sp.Expr) -> bool:
    """Exact root-support criterion of THM-3765 over an algebraic closure."""
    chi = sp.expand(chi)
    psi = sp.expand(psi)
    g_profile = sp.diff(chi, Z)
    if psi == 0:
        return sp.expand(g_profile.subs(Z, 0)) == 0

    psi_zero = sp.expand(psi.subs(Z, 0)) == 0
    axis_safe = (
        not psi_zero
        or (
            sp.expand(sp.diff(psi, Z).subs(Z, 0)) == 0
            and sp.expand(g_profile.subs(Z, 0)) == 0
        )
    )
    if not axis_safe:
        return False

    product_profile = sp.expand(Z * psi)
    eliminant = sp.expand(
        sp.diff(product_profile, Z) ** 2
        - product_profile * g_profile**2
    )
    if eliminant == 0:
        return False

    eliminant_poly = sp.Poly(eliminant, Z, domain=sp.QQ)
    radical = eliminant_poly.exquo(
        sp.gcd(eliminant_poly, eliminant_poly.diff())
    )
    if radical.eval(0) == 0:
        radical = radical.exquo(sp.Poly(Z, Z, domain=sp.QQ))
    multiple_root_part = sp.gcd(
        sp.Poly(psi, Z, domain=sp.QQ),
        sp.Poly(sp.diff(psi, Z), Z, domain=sp.QQ),
    )
    return multiple_root_part.rem(radical).is_zero


# Universal identities with algebraically independent profile coefficients.
chi_coefficients = sp.symbols("chi_0:6")
psi_coefficients = sp.symbols("psi_0:6")
chi_generic = sum(chi_coefficients[i] * Z**i for i in range(6))
psi_generic = sum(psi_coefficients[i] * Z**i for i in range(6))
g_generic = sp.diff(chi_generic, Z)
f_generic = sp.expand(Z * psi_generic)
q_generic = source_component(chi_generic, psi_generic)
z_source = X * T
qx_expected = sp.expand(
    1
    + T * g_generic.subs(Z, z_source)
    + T**2 * sp.diff(psi_generic, Z).subs(Z, z_source)
)
qt_expected = sp.expand(
    X * g_generic.subs(Z, z_source)
    + sp.diff(f_generic, Z).subs(Z, z_source)
)
gate(sp.expand(sp.diff(q_generic, X) - qx_expected) == 0,
     "universal X derivative")
gate(sp.expand(sp.diff(q_generic, T) - qt_expected) == 0,
     "universal T derivative")

critical_x = sp.expand(
    X**2 + Z * g_generic * X + Z**2 * sp.diff(psi_generic, Z)
)
critical_t = sp.expand(g_generic * X + sp.diff(f_generic, Z))
eliminant_generic = sp.expand(
    sp.diff(f_generic, Z) ** 2 - f_generic * g_generic**2
)
gate(
    sp.expand(sp.resultant(critical_x, critical_t, X) - eliminant_generic)
    == 0,
    "universal torus resultant",
)

y_source = sp.expand(X - T * psi_generic.subs(Z, z_source))
delta_generic = sp.expand((L - chi_generic) ** 2 - 4 * f_generic)
gate(sp.expand(jac(q_generic, z_source) - y_source) == 0,
     "generic-fibre Jacobian sign")
gate(
    sp.expand(
        y_source**2
        - delta_generic.subs({L: q_generic, Z: z_source})
    )
    == 0,
    "generic-fibre quadratic equation",
)

# Exhaust every degree-at-most-two profile over {-1,0,1}.  Direct affine
# Groebner smoothness agrees with the axis/resultant/absorbed-root theorem.
grid_total = 0
grid_smooth = 0
values = (-1, 0, 1)
for chi_values in product(values, repeat=3):
    chi = sum(chi_values[i] * Z**i for i in range(3))
    for psi_values in product(values, repeat=3):
        psi = sum(psi_values[i] * Z**i for i in range(3))
        component = source_component(chi, psi)
        critical_ideal = sp.groebner(
            (sp.diff(component, X), sp.diff(component, T)),
            X,
            T,
            order="grevlex",
        )
        directly_smooth = critical_ideal.contains(sp.Integer(1))
        predicted = predicted_smooth(chi, psi)
        gate(directly_smooth == predicted,
             "degree-two smoothness classification")
        grid_total += 1
        grid_smooth += int(directly_smooth)

# The complete degree-at-most-one generic-fibre boundary and its rational
# primitive.  The parameters are kept symbolic, so this includes all of its
# linear degenerations as identities before denominators are specialized.
g, h, p, c = sp.symbols("g h p c", nonzero=True)
chi_boundary = h + g * Z
f_boundary = p * Z + g**2 * Z**2 / 4
psi_boundary = sp.expand(f_boundary / Z)
delta_boundary = sp.expand((L - chi_boundary) ** 2 - 4 * f_boundary)
gate(sp.degree(delta_boundary, Z) == 1, "linear generic-fibre boundary")
gate(
    sp.expand(
        delta_boundary
        - (L - h) ** 2
        + (2 * g * (L - h) + 4 * p) * Z
    )
    == 0,
    "degree-one discriminant",
)
q_boundary = source_component(chi_boundary, psi_boundary)
y_boundary = sp.expand(X - T * psi_boundary.subs(Z, z_source))
d_boundary = sp.expand(g * (q_boundary - h) + 2 * p)
p_boundary = sp.cancel(c * y_boundary / d_boundary)
gate(sp.factor(jac(p_boundary, q_boundary) - c) == 0,
     "nonlinear rational primitive")
gate(
    sp.expand(q_boundary - h - p * T - X * (1 + g * T / 2) ** 2)
    == 0,
    "affine-in-X boundary form",
)
gate(
    sp.expand(sp.diff(q_boundary, X) - (1 + g * T / 2) ** 2) == 0,
    "boundary X derivative",
)
gate(
    sp.expand(
        sp.diff(q_boundary, T) - p - g * X * (1 + g * T / 2)
    )
    == 0,
    "boundary T derivative",
)

# Linear endpoints and the singular p=0 nonlinear endpoint.
pl = sp.symbols("pl", nonzero=True)
q_linear = X + h + pl * T
p_linear = c * (X - pl * T) / (2 * pl)
gate(sp.expand(jac(p_linear, q_linear) - c) == 0,
     "linear two-sided mate")
gate(sp.expand(jac(-c * T, X + h) - c) == 0,
     "linear one-sided mate")
q_singular_edge = sp.expand(
    q_boundary.subs(p, 0)
)
gate(
    sp.expand(sp.diff(q_singular_edge, X).subs(T, -2 / g)) == 0
    and sp.expand(
        sp.diff(q_singular_edge, T).subs({T: -2 / g, X: 0})
    )
    == 0,
    "nonlinear p-zero singular edge",
)

# The top-X descent identity is verified at several arbitrary depths.  It is
# an all-degree coefficient identity in the proof, not a bounded ansatz gate.
A = sp.Function("A")(T)
for depth in range(1, 9):
    top = sp.Function(f"p_{depth}")(T)
    top_equation = sp.expand(
        depth * sp.diff(A, T) * top - A * sp.diff(top, T)
    )
    normalized_derivative = sp.diff(top / A**depth, T)
    gate(
        sp.simplify(top_equation + A ** (depth + 1) * normalized_derivative)
        == 0,
        "top-X descent identity",
    )

# Pure-profile boundary: the generic time form has denominator L-chi.  These
# controls check the Jacobian coordinate and separable simple roots in several
# degrees; the nonzero-residue argument in the theorem is all-degree.
pure_controls = 0
for degree in range(1, 8):
    chi_pure = Z**degree + Z + 1
    q_pure = source_component(chi_pure, sp.Integer(0))
    gate(sp.expand(jac(q_pure, z_source) - X) == 0,
         "pure-profile Jacobian coordinate")
    pure_denominator = sp.Poly(L - chi_pure, Z, domain=sp.QQ.frac_field(L))
    gate(
        sp.gcd(pure_denominator, pure_denominator.diff()).degree() == 0,
        "pure-profile generic simple roots",
    )
    pure_controls += 1

# Pell-Chebyshev hostile smooth profiles from THM-3757.  Their eliminant is
# exactly -psi and every root is absorbed with multiplicity at least two.
tau = sp.symbols("tau")
pell_controls = 0
for depth in range(1, 7):
    cheb_t = sp.chebyshevt(depth, tau)
    cheb_u = sp.chebyshevu(depth - 1, tau)
    g_tau = sp.expand(cheb_t + (tau - 1) * cheb_u)
    s_tau = sp.expand(
        sp.chebyshevu(depth, tau) / (2 * (depth + 1))
        - sp.chebyshevu(depth - 1, tau) / (2 * depth)
    )
    substitution = {tau: 2 * Z - 1}
    g_profile = sp.expand(g_tau.subs(substitution))
    s_profile = sp.expand(s_tau.subs(substitution))
    psi_profile = sp.expand((Z - 1) * s_profile**2)
    f_profile = sp.expand(Z * psi_profile)
    r_profile = sp.expand(
        sp.diff(f_profile, Z) ** 2 - f_profile * g_profile**2
    )
    gate(sp.expand(r_profile + psi_profile) == 0,
         "Pell absorbed-root eliminant")
    gate(predicted_smooth(sp.integrate(g_profile, Z), psi_profile),
         "Pell root-support smoothness")
    pell_controls += 1

# A distinct high-multiplicity absorbed-root family realizes R=-16m^2 F^2.
# It is a hostile control that the criterion permits far more than the Pell
# locus even though the generic differential still obstructs a mate.
absorbed_controls = 0
for multiplicity in range(2, 7):
    base = Z * (Z - 1)
    h_profile = sp.expand(base**multiplicity)
    f_profile = sp.expand(h_profile**2)
    psi_profile = sp.expand(f_profile / Z)
    g_profile = sp.expand(
        -2
        * multiplicity
        * base ** (multiplicity - 1)
        * (Z**2 + (Z - 1) ** 2)
    )
    chi_profile = sp.integrate(g_profile, Z)
    r_profile = sp.expand(
        sp.diff(f_profile, Z) ** 2 - f_profile * g_profile**2
    )
    gate(
        sp.expand(r_profile + 16 * multiplicity**2 * f_profile**2) == 0,
        "power-absorbed eliminant",
    )
    gate(predicted_smooth(chi_profile, psi_profile),
         "power-absorbed root-support smoothness")
    absorbed_controls += 1

# Finite hostile polynomial-mate gates on the first nonlinear rational-exact
# member supplement, but do not replace, the all-degree top-X descent.
nonlinear_sample = sp.expand(q_boundary.subs({g: 2, h: 0, p: 1}))
mate_obstructions = 0
for bound in range(0, 10):
    basis = [
        X**i * T**j
        for i in range(bound + 1)
        for j in range(bound + 1 - i)
    ]
    coefficients = sp.symbols(f"mate_{bound}_0:{len(basis)}")
    prospective = sum(
        coefficient * monomial
        for coefficient, monomial in zip(coefficients, basis)
    )
    equation = sp.Poly(jac(prospective, nonlinear_sample) - 1, X, T)
    system = [coefficient for _, coefficient in equation.terms()]
    gate(sp.linsolve(system, coefficients) == sp.EmptySet,
         "finite nonlinear-boundary polynomial no-mate")
    mate_obstructions += 1

semantic_rows = (
    "component:Q=X+chi(XT)+T*psi(XT)",
    "smooth:axis-safe;R=Fprime^2-F*chiprime^2;nonzero-roots-absorbed",
    "generic:Y^2=(Lambda-chi)^2-4F;J(Q,z)=Y;time=-dz/Y",
    "rational-mate:chi=h+gz;psi=p+g^2*z/4;smooth-iff-g=0-or-p!=0",
    "polynomial-mate:iff-g=0;Q=X+h+pT",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3765-normalized-three-consecutive-charge-radial-classification")
print("scope=algebraically_closed_characteristic_zero;arbitrary_chi,psi_in_k[z]")
print("smoothness=axis_condition_plus_resultant_roots_absorbed_by_multiple_psi_roots")
print("rational_mate=iff_chi=h+gz_and_psi=p+g^2*z/4_with_g=0_or_p_nonzero")
print("polynomial_mate=iff_g=0_equivalently_Q_is_linear")
print("mechanisms=torus_resultant;hyperelliptic_time_form;residue;top_X_descent")
print(f"degree_two_grid={grid_total};smooth={grid_smooth}")
print(
    f"pure_controls={pure_controls};pell_controls={pell_controls};"
    f"power_absorbed_controls={absorbed_controls};"
    f"mate_obstructions={mate_obstructions}"
)
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
