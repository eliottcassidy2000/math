#!/usr/bin/env python3
"""Exact companion for THM-3782's spectral Danielewski completion."""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def jac(first: sp.Expr, second: sp.Expr, left: sp.Symbol, right: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(first, left) * sp.diff(second, right)
        - sp.diff(first, right) * sp.diff(second, left)
    )


def target_bracket(first: sp.Expr, second: sp.Expr, u: sp.Symbol, p: sp.Symbol) -> sp.Expr:
    return sp.factor(
        sp.diff(first, u) * sp.diff(second, p)
        - sp.diff(first, p) * sp.diff(second, u)
    )


u, p, v, e, z = sp.symbols("u p v e z")
x, y = sp.symbols("x y")

spectral_rows = []

# Universal target-field algebra for one through six distinct values.
for arm_count in range(1, 7):
    roots = tuple(sp.Integer(2 * j - arm_count) for j in range(arm_count))
    Phi = sp.expand(sp.prod(z - root for root in roots))
    Phi_v = Phi.subs(z, v)
    V = u * p
    E = sp.cancel(Phi.subs(z, V) / u)

    gate(sp.cancel(u * E - Phi.subs(z, V)) == 0,
         f"spectral relation h={arm_count}")
    gate(target_bracket(u, V, u, p) == u,
         f"log-canonical bracket h={arm_count}")
    gate(
        sp.cancel(target_bracket(u, E, u, p) - sp.diff(Phi, z).subs(z, V)) == 0,
        f"UE bracket h={arm_count}",
    )
    gate(sp.cancel(target_bracket(V, E, u, p) - E) == 0,
         f"VE bracket h={arm_count}")

    bezout_phi, bezout_derivative, gcd_value = sp.gcdex(Phi_v, sp.diff(Phi_v, v), v)
    gate(sp.expand(gcd_value - 1) == 0, f"squarefree gcd h={arm_count}")
    gate(
        sp.expand(
            bezout_phi * Phi_v
            + bezout_derivative * sp.diff(Phi_v, v)
            - 1
        )
        == 0,
        f"Bezout identity h={arm_count}",
    )
    relation = u * e - Phi_v
    unit_expression = sp.expand(
        bezout_phi * u * e
        + bezout_derivative * sp.diff(Phi_v, v)
        - 1
    )
    gate(
        sp.rem(sp.Poly(unit_expression, e), sp.Poly(relation, e)).as_expr() == 0,
        f"unit differential-minor ideal h={arm_count}",
    )
    spectral_rows.append(f"h={arm_count};roots={roots};Phi={sp.sstr(Phi)}")

# The synchronized one-value branch gives a literal polynomial mate.
lam = sp.symbols("lambda")
U_one = x
P_one = y + lam / x
V_one = sp.expand(U_one * P_one)
E_one = sp.cancel((V_one - lam) / U_one)
gate(sp.cancel(jac(U_one, P_one, x, y) - 1) == 0,
     "one-value rational Keller seed")
gate(sp.expand(V_one.subs(x, 0) - lam) == 0,
     "one-value component spectrum")
gate(sp.cancel(E_one - y) == 0, "one-value polynomial mate formula")
gate(jac(U_one, E_one, x, y) == 1, "one-value polynomial Darboux pair")

# A two-value completion whose image misses a divisor.  It verifies why the
# codimension-one coverage hypothesis cannot be deleted: x=U/V is an extra
# polynomial target-field function and (x,y) itself is a Darboux pair.
A_two = 1 + x * y
U_two = sp.expand(x * A_two)
P_two = 1 / x
V_two = sp.expand(U_two * P_two)
E_two = sp.cancel(V_two * (V_two - 1) / U_two)
gate(sp.cancel(jac(U_two, P_two, x, y) - 1) == 0,
     "two-value rational Keller seed")
gate(sp.expand(V_two - A_two) == 0, "two-value blowdown")
gate(sp.cancel(E_two - y) == 0, "two-value spectral completion")
gate(sp.expand(U_two * E_two - V_two * (V_two - 1)) == 0,
     "two-value Danielewski relation")
gate(sp.cancel(U_two / V_two - x) == 0,
     "failed-coverage extra polynomial observable")
gate(jac(x, E_two, x, y) == 1,
     "failed-coverage enlarged-intersection Darboux control")

# Reducedness is load-bearing.  The single spectral factor only pays the
# radical of U=x^2 and leaves E=y/(2x) nonpolynomial.
U_nonreduced = x**2
P_nonreduced = y / (2 * x)
V_nonreduced = sp.cancel(U_nonreduced * P_nonreduced)
E_nonreduced = sp.cancel(V_nonreduced / U_nonreduced)
gate(sp.cancel(jac(U_nonreduced, P_nonreduced, x, y) - 1) == 0,
     "nonreduced hostile Keller seed")
gate(sp.Poly(V_nonreduced, x, y).as_expr() == V_nonreduced,
     "nonreduced hostile blowdown polynomial")
gate(sp.denom(E_nonreduced).has(x),
     "nonreduced hostile completion nonpolynomial")

# A horizontal pole is not cleared by U and lies outside the theorem.
U_horizontal = x
P_horizontal = y + 1 / (x - 1)
V_horizontal = sp.cancel(U_horizontal * P_horizontal)
gate(sp.cancel(jac(U_horizontal, P_horizontal, x, y) - 1) == 0,
     "horizontal-pole hostile Keller seed")
gate(sp.denom(V_horizontal).has(x - 1),
     "horizontal pole survives blowdown")

# THM-3779's m=2 member is a positive two-value, codimension-one-surjective
# control for the unconditional algebra.
m = 2
A_tower = 1 + x * y
B_tower = sp.expand(1 + x ** (2 * m + 1) * A_tower**m)
U_tower = sp.expand(x * A_tower * B_tower)
P_tower = sp.cancel(((m + 1) * B_tower - 1) / (m * x))
V_tower = sp.expand(U_tower * P_tower)
E_tower = sp.cancel(V_tower * (V_tower - 1) / U_tower)
gate(sp.cancel(jac(U_tower, P_tower, x, y) - 1) == 0,
     "tower positive Keller control")
gate(sp.Poly(E_tower, x, y).as_expr() == sp.expand(E_tower),
     "tower positive spectral polynomial")
gate(sp.expand(U_tower * E_tower - V_tower * (V_tower - 1)) == 0,
     "tower positive relation")
gate(sp.expand(V_tower.subs(x, 0) - 1) == 0,
     "tower axis value")
gate(sp.cancel(V_tower.subs(y, -1 / x)) == 0,
     "tower second value")

semantic_rows = (
    "hypotheses:R=k[x,y];U_reduced;J(U,P)=1;V=UP_in_R;no_horizontal_poles",
    "spectrum:V_mod_D_i=lambda_i;Lambda=distinct_values;Phi=product(z-lambda)",
    "completion:E=Phi(V)/U_in_R;relation=UE=Phi(V)",
    "brackets:{U,V}=U;{U,E}=Phi'(V);{V,E}=E;target_etale",
    "conditional:codim1_surjective_implies_R_intersect_k(U,P)=D_Phi",
    "one_value:E=P-lambda/U;{U,E}=1;polynomial_mate",
    "multi_value_over_C:exponent1_nonexact;no_Darboux_in_D_Phi",
    "boundary:failed_codim1_coverage_can_enlarge_intersection_and_restore_Darboux",
    "boundary:nonreduced_U_only_pays_radical;horizontal_poles_not_cleared",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()
packet = hashlib.sha256("\n".join(spectral_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3782-simple-pole-spectral-Danielewski-gate")
print("field=algebraically_closed_characteristic_zero;de_Rham_consequence_over_C")
print("hypotheses=U_reduced;J(U,P)=1;V=UP_polynomial;no_horizontal_poles")
print("component_spectrum=finite_constant_values_Lambda")
print("completion=Phi(z)=product_(lambda_in_Lambda)(z-lambda);E=Phi(V)/U")
print("surface=D_Phi=k[U,V,E]/(UE-Phi(V));smooth_exponent_one")
print("brackets={U,V}=U;{U,E}=Phi'(V);{V,E}=E;map=etale")
print("conditional_intersection=codim1_surjectivity_implies_k[x,y]_intersect_k(U,P)=D_Phi")
print("one_value=polynomial_mate_E=P-lambda/U")
print("two_or_more_values_over_C=no_Darboux_pair_in_D_Phi")
print("failed_surjectivity_boundary=extra_observables_and_Darboux_pair_may_return")
print("nonreduced_boundary=Phi(V)_pays_only_rad(U)")
print("horizontal_boundary=UP_nonpolynomial_if_uncleared_pole")
print("spectral_controls_h=1,2,3,4,5,6")
print(f"packet_sha256={packet}")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("NO_CLAIM=automatic_codim1_surjectivity_or_planar_Jacobian_counterexample")
print("RESULT=PASS")
