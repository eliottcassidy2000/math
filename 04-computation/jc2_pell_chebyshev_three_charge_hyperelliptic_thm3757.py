#!/usr/bin/env python3
"""Exact companion for THM-3757's Pell-Chebyshev three-charge tower."""

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


X, T, Z, L, tau = sp.symbols("X T Z L tau")


def jac(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, X) * sp.diff(second, T)
        - sp.diff(first, T) * sp.diff(second, X)
    )


def tower(depth: int) -> tuple[sp.Expr, ...]:
    """Return G,H,S,chi,psi,F for the theorem's depth."""
    cheb_t = sp.chebyshevt(depth, tau)
    cheb_u = sp.chebyshevu(depth - 1, tau)
    g_tau = sp.expand(cheb_t + (tau - 1) * cheb_u)
    h_tau = sp.expand(cheb_t + (tau + 1) * cheb_u)
    s_tau = sp.expand(
        sp.chebyshevu(depth, tau) / (2 * (depth + 1))
        - sp.chebyshevu(depth - 1, tau) / (2 * depth)
    )
    substitution = {tau: 2 * Z - 1}
    g_profile = sp.expand(g_tau.subs(substitution))
    h_profile = sp.expand(h_tau.subs(substitution))
    s_profile = sp.expand(s_tau.subs(substitution))
    chi = sp.integrate(g_profile, Z)
    chi = sp.expand(chi - chi.subs(Z, 0))
    psi = sp.expand((Z - 1) * s_profile**2)
    product_profile = sp.expand(Z * psi)
    return g_profile, h_profile, s_profile, chi, psi, product_profile


depths = range(1, 9)
degree_rows: list[str] = []
groebner_controls = 0
mate_obstructions = 0
for depth in depths:
    G, H, S, chi, psi, F = tower(depth)
    cheb_t = sp.chebyshevt(depth, tau)
    cheb_u = sp.chebyshevu(depth - 1, tau)

    # Two equivalent Chebyshev/Pell presentations.
    gate(
        sp.expand(
            (tau + 1) * G.subs(Z, (tau + 1) / 2)
            - sp.chebyshevt(depth + 1, tau)
            - cheb_t
        )
        == 0,
        "G Chebyshev quotient",
    )
    gate(
        sp.expand(
            (tau - 1) * H.subs(Z, (tau + 1) / 2)
            - sp.chebyshevt(depth + 1, tau)
            + cheb_t
        )
        == 0,
        "H Chebyshev quotient",
    )
    gate(
        sp.expand(Z * G**2 - (Z - 1) * H**2 - 1) == 0,
        "polynomial Pell identity",
    )

    # The explicit inverse of the transport operator produces S.
    transported = sp.expand(
        (tau**2 - 1)
        * sp.diff(S.subs(Z, (tau + 1) / 2), tau)
        + tau * S.subs(Z, (tau + 1) / 2)
    )
    gate(
        sp.expand(transported - (tau - 1) * H.subs(Z, (tau + 1) / 2) / 2)
        == 0,
        "Chebyshev transport inverse",
    )
    gate(sp.expand(S.subs(Z, 1)) == 0, "S has the absorbed root")
    gate(sp.expand(S.subs(Z, 0) - (-1) ** depth) == 0,
         "S axis value")
    gate(sp.expand(psi.subs(Z, 0) + 1) == 0, "psi axis value")

    # The eliminant is exactly the negative profile.
    gate(
        sp.expand(sp.diff(F, Z) - (Z - 1) * S * H) == 0,
        "F derivative transport",
    )
    eliminant = sp.expand(sp.diff(F, Z) ** 2 - F * G**2)
    gate(sp.expand(eliminant + psi) == 0,
         "three-charge critical eliminant")
    radical_psi = sp.exquo(sp.Poly(psi, Z), sp.gcd(sp.Poly(psi, Z), sp.Poly(sp.diff(psi, Z), Z)))
    gate(
        sp.rem(sp.Poly(sp.diff(psi, Z), Z), radical_psi).is_zero,
        "every eliminant root is an absorbed multiple profile root",
    )

    # Direct resultant of the two torus critical equations.
    critical_x = sp.expand(
        X**2 + Z * G * X + Z**2 * sp.diff(psi, Z)
    )
    critical_t = sp.expand(G * X + sp.diff(F, Z))
    gate(
        sp.expand(sp.resultant(critical_x, critical_t, X) - eliminant) == 0,
        "direct torus resultant",
    )

    # Generic hyperelliptic time-form model.
    delta = sp.expand((L - chi) ** 2 - 4 * F)
    rational_domain = sp.QQ.frac_field(L)
    delta_poly = sp.Poly(delta, Z, domain=rational_domain)
    gate(
        sp.gcd(delta_poly, sp.Poly(sp.diff(delta, Z), Z, domain=rational_domain)).degree()
        == 0,
        "generic hyperelliptic polynomial is squarefree",
    )
    gate(delta_poly.degree() >= depth + 1,
         "generic-fibre degree lower bound")
    gate(sp.degree(G, Z) == depth, "G degree")
    gate(sp.degree(chi, Z) == depth + 1, "chi degree")
    degree_rows.append(
        f"n{depth}:delta{delta_poly.degree()}:genus{(delta_poly.degree()-1)//2}"
    )

    component = sp.expand(
        X + chi.subs(Z, X * T) + T * psi.subs(Z, X * T)
    )
    z_source = X * T
    y_source = sp.expand(X - T * psi.subs(Z, z_source))
    gate(
        sp.expand(jac(component, z_source) - y_source) == 0,
        "hyperelliptic Jacobian coordinate",
    )
    gate(
        sp.expand(
            y_source**2
            - delta.subs({L: component, Z: z_source})
        )
        == 0,
        "hyperelliptic equation on the source",
    )

    # Direct affine smoothness controls at the first four depths.
    if depth <= 4:
        critical_ideal = sp.groebner(
            (sp.diff(component, X), sp.diff(component, T)),
            X,
            T,
            order="grevlex",
        )
        gate(critical_ideal.contains(sp.Integer(1)), "direct smoothness control")
        groebner_controls += 1

    # Bounded polynomial-mate hostiles supplement the all-degree differential
    # proof for the two smallest, structurally different genus branches.
    if depth <= 2:
        max_bound = 12 if depth == 1 else 9
        for bound in range(max_bound + 1):
            basis = [
                X**i * T**j
                for i in range(bound + 1)
                for j in range(bound + 1 - i)
            ]
            coefficients = sp.symbols(
                f"mate_{depth}_{bound}_0:{len(basis)}"
            )
            prospective = sum(
                coefficient * monomial
                for coefficient, monomial in zip(coefficients, basis)
            )
            equation = sp.Poly(jac(prospective, component) - 1, X, T)
            system = [coefficient for _, coefficient in equation.terms()]
            gate(
                sp.linsolve(system, coefficients) == sp.EmptySet,
                "finite Pell-tower mate obstruction",
            )
            mate_obstructions += 1

# The rational n=1 boundary has a reducible fibre and a logarithmic generic
# time form; both identities are exact regression controls.
G1, H1, S1, chi1, psi1, F1 = tower(1)
component1 = sp.expand(X + chi1.subs(Z, X * T) + T * psi1.subs(Z, X * T))
w = X * T - 1
factor_left = sp.expand(1 + T * w)
factor_right = sp.expand(X + w**2)
gate(sp.expand(component1 + 1 - factor_left * factor_right) == 0,
     "n=1 split fibre")
delta1 = sp.Poly(sp.expand((L - chi1) ** 2 - 4 * F1), Z)
gate(delta1.degree() == 2, "n=1 conic boundary")
gate(sp.expand(delta1.LC() + 4 * L + 3) == 0,
     "n=1 nonzero logarithmic infinity coefficient")

semantic_rows = (
    "pell:z*G_n^2-(z-1)*H_n^2=1",
    "transport:S_n=U_n/[2(n+1)]-U_(n-1)/(2n)",
    "component:Q_n=X+chi_n(XT)+T*(XT-1)*S_n(XT)^2",
    "smooth:resultant=-(z-1)*S_n^2;all_roots_absorbed",
    "mate:none_rational;n1_log_residue;n>=2_holomorphic_hyperelliptic_time_form",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3757-Pell-Chebyshev-three-charge-hyperelliptic-obstruction-tower")
print("scope=algebraically_closed_characteristic_zero;all_integer_depths_n>=1")
print("component=Q_n=X+chi_n(XT)+T*(XT-1)*S_n(XT)^2")
print("smoothness=critical_resultant_equals_negative_profile_with_absorbed_roots")
print("jacobian_mate=no_rational_mate_at_every_depth")
print("mechanism=n1_logarithmic_infinity_residue;n>=2_nonzero_holomorphic_time_form")
print(f"depth_controls={len(tuple(depths))};" + ",".join(degree_rows))
print(
    f"groebner_controls={groebner_controls};"
    f"mate_obstructions={mate_obstructions};n1_split_fibre=1"
)
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
