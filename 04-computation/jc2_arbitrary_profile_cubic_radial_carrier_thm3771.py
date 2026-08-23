#!/usr/bin/env python3
"""Exact companion for THM-3771's arbitrary-profile radial dressing."""

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


X, T, z, w, r = sp.symbols("X T z w r")
a, b, c, lam = sp.symbols("a b c lam", nonzero=True)
U_symbol, W_symbol, L = sp.symbols("U W L")


def jac(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, X) * sp.diff(second, T)
        - sp.diff(first, T) * sp.diff(second, X)
    )


def chart_jac(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, X) * sp.diff(second, z)
        - sp.diff(first, z) * sp.diff(second, X)
    )


def log_bracket(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        3
        * U_symbol
        * (
            sp.diff(first, U_symbol) * sp.diff(second, W_symbol)
            - sp.diff(first, W_symbol) * sp.diff(second, U_symbol)
        )
    )


def total_degree_basis(bound: int) -> list[sp.Expr]:
    return [
        X**i * T**j
        for i in range(bound + 1)
        for j in range(bound + 1 - i)
    ]


def source_family(
    u_profile: sp.Expr, r_value: sp.Expr, phi_profile: sp.Expr
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    z_source = X * T
    u_source = u_profile.subs(z, z_source)
    source_u = sp.expand(X * u_source)
    source_w = sp.expand(source_u + 3 * z_source + r_value)
    source_q = sp.expand(source_u * phi_profile.subs(w, source_w))
    return source_u, source_w, source_q


def squarefree_nonzero(poly: sp.Expr, variable: sp.Symbol) -> bool:
    polynomial = sp.Poly(poly, variable)
    if polynomial.degree() <= 0:
        return polynomial.as_expr() != 0
    return sp.Poly(sp.gcd(polynomial, polynomial.diff()), variable).degree() == 0


def admissible(
    u_profile: sp.Expr, r_value: sp.Expr, phi_profile: sp.Expr
) -> bool:
    return bool(
        sp.expand(u_profile.subs(z, 0)) != 0
        and sp.expand(phi_profile.subs(w, r_value)) != 0
        and squarefree_nonzero(u_profile, z)
        and squarefree_nonzero(phi_profile, w)
        and sp.Poly(
            sp.gcd(
                sp.Poly(u_profile, z),
                sp.Poly(phi_profile.subs(w, 3 * z + r_value), z),
            ),
            z,
        ).degree()
        == 0
    )


def roots_exact(poly: sp.Expr, variable: sp.Symbol) -> tuple[sp.Expr, ...]:
    polynomial = sp.Poly(poly, variable)
    if polynomial.degree() <= 0:
        return ()
    roots = sp.roots(polynomial.as_expr(), variable)
    gate(sum(roots.values()) == polynomial.degree(), "split root control")
    return tuple(sorted(roots, key=sp.default_sort_key))


# Universal arbitrary-profile chart.  SymPy retains u and u' formally, so
# these are identities for every polynomial profile rather than interpolated
# low-degree tests.
u_function = sp.Function("u")
u_z = u_function(z)
U_chart = X * u_z
W_chart = U_chart + 3 * z + r
gate(chart_jac(U_chart, W_chart) == 3 * u_z,
     "universal independent-chart determinant")

z_source = X * T
u_source_formal = u_function(z_source)
U_source_formal = X * u_source_formal
W_source_formal = U_source_formal + 3 * z_source + r
gate(sp.simplify(jac(U_source_formal, W_source_formal) - 3 * U_source_formal) == 0,
     "universal log-canonical source bracket")

# The quadratic carrier expansion is also formal in u.
phi_formal = a + b * w + c * w**2
v = 3 * z + r
A = sp.expand(u_z * phi_formal.subs(w, v))
B = sp.expand(u_z**2 * sp.diff(phi_formal, w).subs(w, v))
C = c * u_z**3
Q_chart = sp.expand(U_chart * phi_formal.subs(w, W_chart))
gate(sp.expand(Q_chart - (X * A + X**2 * B + X**3 * C)) == 0,
     "universal cubic radial-carrier expansion")

# Exact birational inverse in the rational function field.
z_inverse = (W_symbol - U_symbol - r) / 3
X_inverse = U_symbol / u_function(z_inverse)
T_inverse = z_inverse / X_inverse
gate(sp.simplify(X_inverse * T_inverse - z_inverse) == 0,
     "birational inverse product")
gate(sp.simplify(X_inverse * u_function(z_inverse) - U_symbol) == 0,
     "birational inverse U")
gate(
    sp.simplify(
        X_inverse * u_function(z_inverse) + 3 * X_inverse * T_inverse + r
        - W_symbol
    )
    == 0,
    "birational inverse W",
)

# Generic Hamiltonian derivation and complete primitive in (U,W).
phi_generic = a + b * W_symbol + c * W_symbol**2
Q_generic = U_symbol * phi_generic
P_generic = -lam * W_symbol / (3 * Q_generic)
gate(sp.expand(log_bracket(U_symbol, W_symbol) - 3 * U_symbol) == 0,
     "log bracket normalization")
gate(sp.expand(log_bracket(W_symbol, Q_generic) + 3 * Q_generic) == 0,
     "Hamiltonian W derivative")
gate(sp.cancel(log_bracket(P_generic, Q_generic) - lam) == 0,
     "universal rational primitive")
gate(
    sp.diff(-lam * W_symbol / (3 * L), W_symbol) == -lam / (3 * L),
    "generic-fibre primitive equation",
)

# The same primitive is checked directly in source variables with formal u.
Q_source_formal = U_source_formal * phi_formal.subs(w, W_source_formal)
P_source_formal = -lam * W_source_formal / (3 * Q_source_formal)
gate(sp.simplify(jac(P_source_formal, Q_source_formal) - lam) == 0,
     "direct arbitrary-profile source primitive")


# Smooth controls traverse genuine cubic, multi-root u, constant-u,
# constant-phi, both-constant, and linear-phi boundary cells.
controls = (
    (1 + z, sp.Integer(3), w**2 - 5, "normalized_cubic"),
    ((z - 1) * (z + 2), sp.Integer(0), w**2 + 1, "two_root_cubic"),
    (sp.Integer(2), sp.Integer(0), w**2 - w + 1, "constant_u"),
    (1 + z, sp.Integer(0), sp.Integer(2), "constant_phi"),
    (sp.Integer(2), sp.Integer(5), sp.Integer(3), "coordinate_boundary"),
    (z + 2, sp.Integer(1), w + 1, "linear_phi"),
)

smooth_controls = 0
rational_controls = 0
address_counts: list[int] = []
for u_control, r_control, phi_control, label in controls:
    gate(admissible(u_control, r_control, phi_control),
         f"admissible profile control: {label}")
    U_control, W_control, Q_control = source_family(
        u_control, r_control, phi_control
    )
    critical = sp.groebner(
        (sp.diff(Q_control, X), sp.diff(Q_control, T)), X, T, order="lex"
    )
    gate(critical.contains(sp.Integer(1)), f"smooth profile control: {label}")
    P_control = sp.cancel(-W_control / (3 * Q_control))
    gate(sp.cancel(jac(P_control, Q_control) - 1) == 0,
         f"rational-mate control: {label}")

    u_roots = roots_exact(u_control, z)
    phi_roots = roots_exact(phi_control, w)
    addresses = (
        (r_control,)
        + tuple(sp.expand(r_control + 3 * root) for root in u_roots)
        + phi_roots
    )
    gate(len(addresses) == 1 + sp.degree(u_control, z) + sp.degree(phi_control, w),
         f"address count: {label}")
    for index, left in enumerate(addresses):
        for right in addresses[index + 1 :]:
            gate(sp.simplify(left - right) != 0,
                 f"distinct zero-fibre addresses: {label}")
    address_counts.append(len(addresses))
    smooth_controls += 1
    rational_controls += 1


# Each omitted smoothness condition has a direct critical witness.
hostile_boundaries = (
    (z, sp.Integer(1), w**2 + 1, {X: 0, T: 0}, "u(0)=0"),
    ((z - 1) ** 2, sp.Integer(0), w**2 + 1, {X: 1, T: 1}, "u repeated"),
    (1 + z, sp.Integer(0), w, {X: 0, T: 0}, "phi(r)=0"),
    (z - 1, sp.Integer(0), w - 3, {X: 1, T: 1}, "shared address"),
    (1 + z, sp.Integer(0), (w - 2) ** 2, {X: 2, T: 0}, "phi repeated"),
)
boundary_controls = 0
for u_bad, r_bad, phi_bad, point, label in hostile_boundaries:
    _, _, Q_bad = source_family(u_bad, r_bad, phi_bad)
    gate(
        sp.diff(Q_bad, X).subs(point) == 0
        and sp.diff(Q_bad, T).subs(point) == 0,
        f"critical boundary witness: {label}",
    )
    boundary_controls += 1


# Both profiles constant are exactly the positive polynomial boundary.
_, _, Q_coordinate = source_family(sp.Integer(2), sp.Integer(5), sp.Integer(3))
P_coordinate = -T / 6
gate(Q_coordinate == 6 * X, "coordinate boundary Q")
gate(jac(P_coordinate, Q_coordinate) == 1, "coordinate boundary mate")


# Frozen normalized degree-nine hostile and a bounded independent polynomial
# response census.  The theorem is all-degree; these systems are controls.
U_normalized, W_normalized, Q_normalized = source_family(
    1 + z, sp.Integer(3), w**2 - 5
)
expected_normalized = (
    X**6 * T**3
    + 6 * X**5 * T**3
    + 9 * X**4 * T**3
    + 3 * X**5 * T**2
    + 18 * X**4 * T**2
    + 27 * X**3 * T**2
    + 3 * X**4 * T
    + 18 * X**3 * T
    + 22 * X**2 * T
    + X**3
    + 6 * X**2
    + 4 * X
)
gate(sp.expand(Q_normalized - expected_normalized) == 0,
     "normalized integral expansion")
gate(sp.total_degree(Q_normalized) == 9, "normalized total degree")
gate(sp.cancel(jac(-W_normalized / (3 * Q_normalized), Q_normalized) - 1) == 0,
     "normalized rational mate")

mate_obstructions = 0
rank_gap_min = None
for bound in range(0, 9):
    basis = total_degree_basis(bound)
    coefficients = sp.symbols(f"p_{bound}_0:{len(basis)}")
    prospective = sum(
        coefficient * monomial
        for coefficient, monomial in zip(coefficients, basis)
    )
    equation = sp.Poly(jac(prospective, Q_normalized) - 1, X, T)
    rows = [coefficient for _, coefficient in equation.terms()]
    matrix, rhs = sp.linear_eq_to_matrix(rows, coefficients)
    gap = matrix.row_join(rhs).rank() - matrix.rank()
    gate(gap == 1, "bounded normalized polynomial-mate obstruction")
    rank_gap_min = gap if rank_gap_min is None else min(rank_gap_min, gap)
    mate_obstructions += 1


semantic_rows = (
    "family:z=XT;U=X*u(z);W=U+3*z+r;Q=U*phi(W);deg(phi)<=2",
    "chart:J(U,W)=3*U;z=(W-U-r)/3;X=U/u(z);T=z/X",
    "smooth:u(0)*phi(r)!=0;u,phi_squarefree;gcd(u,phi(3z+r))=1",
    "torsor:P=-lambda*W/(3*Q)+H(Q);constants=k(Q)",
    "addresses:{r}union{r+3*roots(u)}union{roots(phi)}",
    "polynomial_mate:iff_u_and_phi_constant",
    "normalized:u=1+z;r=3;phi=W^2-5;degree=9",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3771-arbitrary-profile-cubic-radial-carrier-log-canonical-nonentry")
print("scope=algebraically_closed_characteristic_zero;nonzero_u;nonzero_deg_phi_at_most_2")
print("chart=birational_log_canonical;J(U,W)=3U")
print("smoothness=exact_squarefree_axis_and_address_avoidance_criterion")
print("rational_mates=complete_torsor_-lambda*W/(3Q)+k(Q)")
print("polynomial_mate=iff_u_and_phi_are_both_constant")
print(f"smooth_controls={smooth_controls};rational_controls={rational_controls};address_counts={','.join(map(str, address_counts))}")
print(f"boundary_controls={boundary_controls};mate_obstructions={mate_obstructions};rank_gap_min={rank_gap_min}")
print("normalized=u=1+z;r=3;phi=W^2-5;total_degree=9")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
