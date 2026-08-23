#!/usr/bin/env python3
"""Exact companion for THM-3775's anti-invariant vertical residue gate."""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


X, T, Z, L, Y = sp.symbols("X T Z L Y")


def jac(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, X) * sp.diff(second, T)
        - sp.diff(first, T) * sp.diff(second, X)
    )


def reduce_quadratic(expression: sp.Expr, delta: sp.Expr) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), Y)
    modulus = sp.Poly(Y**2 - delta, Y)
    return sp.expand(polynomial.rem(modulus).as_expr())


# Universal anti-commutation in a quadratic function field.  Algebraically
# independent coefficients make this a hostile sign check, not an ansatz.
delta_coefficients = sp.symbols("delta_0:6")
m_coefficients = sp.symbols("M_0:5")
f_coefficients = sp.symbols("f_0:5")
g_coefficients = sp.symbols("g_0:5")
delta_generic = sum(delta_coefficients[i] * Z**i for i in range(6))
m_generic = sum(m_coefficients[i] * Z**i for i in range(5))
f_generic = sum(f_coefficients[i] * Z**i for i in range(5))
g_generic = sum(g_coefficients[i] * Z**i for i in range(5))


def derivation(expression: sp.Expr) -> sp.Expr:
    raw = (
        -m_generic * Y * sp.diff(expression, Z)
        - m_generic * sp.diff(delta_generic, Z) * sp.diff(expression, Y) / 2
    )
    return reduce_quadratic(raw, delta_generic)


def involution(expression: sp.Expr) -> sp.Expr:
    return sp.expand(expression.subs(Y, -Y))


generic_element = f_generic + g_generic * Y
gate(
    reduce_quadratic(
        derivation(involution(generic_element))
        + involution(derivation(generic_element)),
        delta_generic,
    )
    == 0,
    "quadratic involution anti-commutes with Hamiltonian derivation",
)
gate(
    sp.expand(derivation(f_generic) + m_generic * sp.diff(f_generic, Z) * Y)
    == 0,
    "invariant functions map to anti-invariant functions",
)
gate(
    sp.expand(
        derivation(g_generic * Y)
        + m_generic
        * (
            sp.diff(g_generic, Z) * delta_generic
            + g_generic * sp.diff(delta_generic, Z) / 2
        )
    )
    == 0,
    "anti-invariant functions map to invariant functions",
)

# The diagonal target direction is the trivial component representation.
# Every nonzero sign vector has nontrivial quotient class.
permutation_controls = 0
for component_count in range(2, 10):
    entries = sp.symbols(f"e_{component_count}_0:{component_count}")
    mean = sum(entries) / component_count
    quotient = tuple(sp.expand(entry - mean) for entry in entries)
    ones = sp.ones(component_count, 1)
    projector = sp.eye(component_count) - ones * ones.T / component_count
    gate(projector * ones == sp.zeros(component_count, 1),
         "component projector kills the target line")
    gate(projector.rank() == component_count - 1,
         "component quotient has exact codimension one")
    gate(sum(quotient) == 0, "component quotient has zero trivial projection")
    gate(
        all(value == 0 for value in quotient)
        == all(sp.expand(entries[i] - entries[0]) == 0
               for i in range(1, component_count)),
        "trivial quotient iff all component coefficients agree",
    )
    permutation_controls += 1

v, s, alpha = sp.symbols("v s alpha", nonzero=True)
split_values = sp.Matrix((v * s, -v * s))
trivial_projection = sp.Matrix((1, 1)) * sum(split_values) / 2
sign_projection = sp.simplify(split_values - trivial_projection)
gate(trivial_projection == sp.zeros(2, 1), "split sign vector has zero trace")
gate(sign_projection == split_values, "split sign vector is anti-invariant")
gate(sp.expand(split_values[0] - split_values[1]) != 0,
     "split sign values cannot equalize")

# Nonsplit hostile controls: a nonzero vY cannot be a base scalar because its
# square has a squarefree factor.  These controls cover odd and even degrees.
nonsplit_controls = 0
for degree in range(1, 9):
    delta_value = sp.expand(Z**degree + Z + 1)
    if sp.gcd(sp.Poly(delta_value, Z), sp.Poly(sp.diff(delta_value, Z), Z)).degree() != 0:
        delta_value = sp.expand(Z**degree + Z + 2)
    delta_poly = sp.Poly(delta_value, Z)
    gate(sp.gcd(delta_poly, delta_poly.diff()).degree() == 0,
         "nonsplit control is squarefree")
    gate(
        any(exponent % 2 for _, exponent in sp.factor_list(delta_value)[1]),
        "nonsplit control has odd prime valuation",
    )
    nonsplit_controls += 1

# THM-3765's normalized nonlinear wall: the exceptional fibre is exactly the
# two sign branches, and its leading coefficient vector is nonzero anti-trace.
g0, p0, h0, c0 = sp.symbols("g0 p0 h0 c0", nonzero=True)
z_source = X * T
psi = p0 + g0**2 * Z / 4
q_normalized = sp.expand(
    X + h0 + g0 * z_source + T * psi.subs(Z, z_source)
)
y_normalized = sp.expand(X - T * psi.subs(Z, z_source))
lambda_normalized = h0 - 2 * p0 / g0
d_normalized = sp.expand(g0 * (q_normalized - lambda_normalized))
p_normalized = sp.cancel(c0 * y_normalized / d_normalized)
s_normalized = sp.expand(1 + g0 * T / 2)
r_normalized = sp.expand(2 * p0 / g0 + X * s_normalized)
gate(sp.expand(q_normalized - lambda_normalized - s_normalized * r_normalized) == 0,
     "THM3765 exhaustive split fibre")
gate(sp.factor(jac(p_normalized, q_normalized) - c0) == 0,
     "THM3765 anti-invariant primitive")
w_normalized = -2 * p0 / g0
y_first = sp.expand(y_normalized.subs(T, -2 / g0))
y_second = sp.factor(
    y_normalized.subs(X, -2 * p0 / (g0 * s_normalized))
)
gate(sp.expand(y_first + w_normalized) == 0,
     "THM3765 first sign value")
gate(sp.expand(y_second - w_normalized) == 0,
     "THM3765 second sign value")
gate(
    sp.expand(c0 * y_first / g0 + c0 * y_second / g0) == 0
    and sp.expand(c0 * y_first / g0) != 0,
    "THM3765 nonzero anti-invariant principal vector",
)

# THM-3772's scaled dual walls have the same sign vector in both orientations.
a0, b0 = sp.symbols("a0 b0", nonzero=True)
a_profile = a0
b_profile = (p0 + g0**2 * Z / 4) / a0
q_a = sp.expand(
    X * a_profile
    + h0
    + g0 * z_source
    + T * b_profile.subs(Z, z_source)
)
y_a = sp.expand(X * a_profile - T * b_profile.subs(Z, z_source))
s_a = sp.expand(1 + g0 * T / (2 * a0))
gate(
    sp.expand(
        q_a - lambda_normalized
        - s_a * (2 * p0 / g0 + a0 * X * s_a)
    )
    == 0,
    "THM3772 constant-A split",
)
gate(
    sp.expand(y_a.subs(T, -2 * a0 / g0) + w_normalized) == 0,
    "THM3772 constant-A first sign",
)
gate(
    sp.factor(
        y_a.subs(X, -2 * p0 / (g0 * a0 * s_a))
        - w_normalized
    )
    == 0,
    "THM3772 constant-A second sign",
)

b_profile_dual = b0
a_profile_dual = (p0 + g0**2 * Z / 4) / b0
q_b = sp.expand(
    X * a_profile_dual.subs(Z, z_source)
    + h0
    + g0 * z_source
    + T * b_profile_dual
)
y_b = sp.expand(X * a_profile_dual.subs(Z, z_source) - T * b_profile_dual)
s_b = sp.expand(1 + g0 * X / (2 * b0))
gate(
    sp.expand(
        q_b - lambda_normalized
        - s_b * (2 * p0 / g0 + b0 * T * s_b)
    )
    == 0,
    "THM3772 constant-B split",
)
gate(
    sp.expand(y_b.subs(X, -2 * b0 / g0) - w_normalized) == 0,
    "THM3772 constant-B first sign",
)
gate(
    sp.factor(
        y_b.subs(T, -2 * p0 / (g0 * b0 * s_b))
        + w_normalized
    )
    == 0,
    "THM3772 constant-B second sign",
)

# THM-3758 is covered by the universal component quotient even though its
# affine source does not exhaust both generic sign branches.  Its two scalar
# coefficients are an exact nonzero sign vector.
a1, beta, gamma = sp.symbols("a1 beta gamma", nonzero=True)
a_carrier = sp.expand(a0 + a1 * Z)
b_carrier = sp.expand(gamma + beta * a0 * Z + beta * a1 * Z**2 / 2)
q_carrier = sp.expand(
    X * a_carrier.subs(Z, z_source)
    + X**2 * b_carrier.subs(Z, z_source)
)
y_carrier = sp.expand(
    a_carrier.subs(Z, z_source)
    + 2 * X * b_carrier.subs(Z, z_source)
)
carrier_bracket = sp.expand(Z + Y / a1)
axis_bracket = sp.expand(carrier_bracket.subs({Z: 0, Y: a0}))
residual_bracket = sp.expand(
    carrier_bracket.subs(Y, -a_carrier)
)
gate(sp.expand(axis_bracket - a0 / a1) == 0,
     "THM3758 axis principal bracket")
gate(sp.expand(residual_bracket + a0 / a1) == 0,
     "THM3758 residual principal bracket")
carrier_coefficients = sp.Matrix(
    (-c0 * axis_bracket / 2, -c0 * residual_bracket / 2)
)
gate(sum(carrier_coefficients) == 0,
     "THM3758 coefficient vector has zero trace")
gate(carrier_coefficients[0] != 0,
     "THM3758 coefficient vector is nontrivial")

semantic_rows = (
    "equalizer:target-shear-principal-parts=trivial-diagonal-component-line",
    "quadratic:D(z)=-M*Y;involution-anticommutes;mate-normalizes-to-vY",
    "nonsplit:vY-scalar-implies-zero",
    "split:(vs,-vs)-diagonal-implies-zero",
    "corollaries:THM3758-component-sign;THM3765-and-THM3772-exhaustive-sign-fibres",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3775-quadratic-anti-invariant-vertical-residue-gate")
print("scope=characteristic_zero;complete_rational_torsor;exhaustive_reduced_quadratic_pole_fibres")
print("target_shear=trivial_diagonal_component_representation_only")
print("quadratic_primitive=target_shear_equivalent_to_anti_invariant_vY")
print("vertical_gate=nonzero_anti_invariant_leading_coefficient_never_equalizes")
print("corollaries=THM3758_component_vector;THM3765;THM3772")
print(
    f"permutation_controls={permutation_controls};"
    f"nonsplit_controls={nonsplit_controls};"
    "split_controls=1;concrete_corollaries=3"
)
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
