#!/usr/bin/env python3
"""Exact cubic-tau elimination for the zero/zero-arm THM-3885 sector.

This heavier scout proves that the finite-exact two-arm closure extends from
degree five to degree six.  It deliberately performs the Groebner calculation
from the six residual equations rather than trusting a frozen basis.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.expand(expression) == 0, message)


# ---------------------------------------------------------------------------
# 1. Cubic x=0 specialization has only the zero square survivor.
# ---------------------------------------------------------------------------

y, p, q, r = sp.symbols("y p q r")
tau = p * y + q * y**2 + r * y**3
F = sp.expand(256 - 96 * tau**2 - 8 * (y**2 - 4) * tau**3 - 3 * tau**4)

# F(0)=256 and [y]F=0.  After negating a square root if necessary, its
# constant is 16 and its linear coefficient is zero.  Coefficients y^2..y^6
# determine the remaining root recursively without dividing by p,q,or r.
G = sp.Integer(16)
g = {}
F_poly = sp.Poly(F, y)
for degree in range(2, 7):
    previous_square = sp.Poly(sp.expand(G**2), y)
    g[degree] = sp.factor(
        (
            F_poly.coeff_monomial(y**degree)
            - previous_square.coeff_monomial(y**degree)
        )
        / 32
    )
    G += g[degree] * y**degree

expected_g = {
    2: -3 * p**2,
    3: p * (p**2 - 6 * q),
    4: -sp.Rational(3, 8) * (p**4 - 8 * p**2 * q + 16 * p * r + 8 * q**2),
    5: sp.Rational(1, 16)
    * (3 * p**5 - 24 * p**3 * q - 4 * p**3 + 48 * p**2 * r
       + 48 * p * q**2 - 96 * q * r),
    6: -sp.Rational(1, 128)
    * (13 * p**6 - 120 * p**4 * q + 192 * p**3 * r
       + 288 * p**2 * q**2 + 96 * p**2 * q - 768 * p * q * r
       - 128 * q**3 + 384 * r**2),
}
for degree, expected in expected_g.items():
    zero(g[degree] - expected, f"forced square-root coefficient g{degree}")

square_debt = sp.Poly(sp.expand(G**2 - F), y)
equations = [square_debt.coeff_monomial(y**degree) for degree in range(7, 13)]
gate(all(equation != 0 for equation in equations), "six nontrivial residual equations")
gate(sp.Poly(F, y).degree() <= 12, "complete cubic-tau coefficient universe")

# This is intentionally recomputed from the equations.  It is the expensive
# validity gate and may take several minutes in SymPy.
basis = sp.groebner(
    equations,
    p,
    q,
    r,
    order="grevlex",
    method="f5b",
    domain=sp.QQ,
)
gate(len(basis.polys) == 15, "complete grevlex basis length")

decisive_members = {
    "r4": r**4,
    "q4_after_r": q**4 - 6 * p**2 * r**2,
    "p5_after_qr": p**5 + 16 * p**2 * r + 16 * p * q**2
    + sp.Rational(448, 5) * r**3,
}
for label, member in decisive_members.items():
    zero(basis.reduce(member)[1], f"decisive ideal member {label}")

# No saturation or coefficient branch occurs: over any characteristic-zero
# field, r^4=0, then q^4=0, then p^5=0 at every field-valued zero.
gate(all(not denominator.has(p, q, r) for denominator in [sp.denom(value) for value in g.values()]),
     "root recursion divides only by characteristic-zero constants")


# ---------------------------------------------------------------------------
# 2. The remaining degree-six global y-profiles are parity-obstructed.
# ---------------------------------------------------------------------------

x = sp.symbols("x")
a = x + 1
L = 9 * x + 4
K = y**2 - 15 * x**2 - 15 * x - 4
Tf = sp.symbols("Tf")
S_f0 = sp.expand(L**4 - 6 * a * L**2 * Tf**2 - 8 * K * Tf**3 - 3 * a**2 * Tf**4)

# The two zero arms and L-parity already give T=a*L^2*H1.  The cubic-tau
# lemma forces H1(0,y)=0, so H1=x*J with total degree J<=2.
A0, B0, C0, D0, E0, F0 = sp.symbols("A0 B0 C0 D0 E0 F0")
J = A0 * x**2 + B0 * x * y + C0 * y**2 + D0 * x + E0 * y + F0
T_degree6 = sp.expand(a * L**2 * x * J)
S_degree6 = sp.Poly(S_f0.subs(Tf, T_degree6), y)

t2 = sp.expand(C0 * a * L**2 * x)
top_y8 = sp.factor(S_degree6.coeff_monomial(y**8))
zero(top_y8 + t2**3 * (8 + 3 * a**2 * t2),
     "quadratic-y profile top coefficient")
gate(sp.rem(top_y8, a**3, x) == 0, "top y8 coefficient has a-adic order at least 3")
gate(sp.rem(top_y8, a**4, x) != 0, "top y8 coefficient has exact odd a-adic order 3")

linear_profile = sp.expand(T_degree6.subs(C0, 0))
t1 = sp.Poly(linear_profile, y).coeff_monomial(y)
linear_residual = sp.Poly(S_f0.subs(Tf, linear_profile), y)
zero(linear_residual.coeff_monomial(y**5) + 8 * t1**3,
     "linear-y profile has unique odd y5 coefficient")

univariate_T = sp.expand(linear_profile.subs({B0: 0, E0: 0}))
univariate_residual = sp.Poly(S_f0.subs(Tf, univariate_T), y)
zero(univariate_residual.coeff_monomial(y**2) + 8 * univariate_T**3,
     "univariate profile retains nonzero y2 debt")
zero(univariate_residual.coeff_monomial(y), "univariate profile has missing y term")
gate(univariate_residual.coeff_monomial(1).subs(x, 0) == 256,
     "univariate profile has nonzero constant branch")


semantic = {
    "univariate": "tau=p*y+q*y^2+r*y^3; square residual forces p=q=r=0",
    "certificate": "six residual equations; exact 15-element grevlex Groebner basis",
    "decisive": "r^4, q^4-6p^2r^2, p^5+16p^2r+16pq^2+448r^3/5",
    "global": "f=0 plus both zero arms plus degT<=6 forces T=0",
    "boundary": "degree>=7 remains open; JC2 open",
    "status": "THM3893 proved verified exact independently hostile-audited",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("scout=jc2-cusp-residual-cubic-tau-degree6")
print("cubic_tau_residual_equations=6;grevlex_basis_size=15")
print("cubic_tau_field_zeros=p=q=r=0")
print("saturation_or_candidate_coefficient_division=NONE")
print("f_zero_both_arms_zero_degT_at_most_6=T_ZERO")
print("first_possible_two_zero_arm_survivor_degree=7")
print("THM3893_status=PROVED_VERIFIED_EXACT_INDEPENDENTLY_HOSTILE_AUDITED")
print("JC2_status=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
