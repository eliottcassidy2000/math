#!/usr/bin/env python3
"""Independent (s,p)-coordinate audit for THM-4140.

This route does not import the primary checker.  It derives the rational
critical pair after s=XT, p=T+s^2, eliminates s rather than X, and restores
the p=0 and T=0 strata separately.
"""

from __future__ import annotations

import hashlib

import sympy as sy


def check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


s, p, r = sy.symbols("s p r")
t = p - s**2
delta = sy.Rational(5696, 15) / (6 * r**2 + 7)
kappa = sy.cancel(delta * r**2)
phi = sy.cancel(2 * delta * r)

H = (
    -3 * p
    + sy.Rational(8, 3) * p**2
    - sy.Rational(1376, 135) * p**3
    + kappa * s**2 * p**2
    + phi * s * p**3
    + delta * p**4
)
Gsp = -s**2 / (2 * t) + H

A = sy.expand(-s + t**2 * p * (2 * kappa * s + phi * p))
B = sy.expand(
    -6
    + sy.Rational(32, 3) * p
    - sy.Rational(2752, 45) * p**2
    + 6 * kappa * s**2 * p
    + 7 * phi * s * p**2
    + 8 * delta * p**3
)

# These identities prove equivalence to the rational gradient for p*t!=0.
Gs = sy.diff(Gsp, s)
Gp = sy.diff(Gsp, p)
check(sy.factor(t**2 * Gs - p * A) == 0, "A is not the s-gradient numerator")
check(sy.factor(2 * t**2 * Gp - (t**2 * B - s * A)) == 0, "B gradient identity failed")
check(sy.Poly(A, s).degree() == 5, "wrong s-degree for A")
check(sy.Poly(B, s).degree() == 2, "wrong s-degree for B")
check(
    sy.factor(sy.Poly(A, s).LC() - 2 * kappa * p) == 0,
    "wrong leading row for A",
)
check(
    sy.factor(sy.Poly(B, s).LC() - 6 * kappa * p) == 0,
    "wrong leading row for B",
)

resultant_sp = sy.factor(sy.resultant(A, B, s))
R14_expr = sy.cancel(resultant_sp / p**2)
R14 = sy.Poly(R14_expr, p)
check(R14.degree() == 14, "independent residual is not degree fourteen")

expected_r14_lead = -(
    sy.Integer(1539884872666062544522946019328)
    * r**4
    / (sy.Integer(512578125) * (6 * r**2 + 7) ** 6)
)
expected_r14_constant = -(
    sy.Integer(112127901696) * r**4
    / (sy.Integer(25) * (6 * r**2 + 7) ** 2)
)
check(sy.factor(R14.LC() - expected_r14_lead) == 0, "R14 leading row changed")
check(sy.factor(R14.TC() - expected_r14_constant) == 0, "R14 constant row changed")
check(sy.factor(resultant_sp - p**2 * R14.as_expr()) == 0, "independent resultant reconstruction failed")

monic_coefficients = [sy.factor(coefficient / R14.LC()) for coefficient in R14.all_coeffs()]
coefficient_payload = "\n".join(map(str, monic_coefficients)) + "\n"
coefficient_sha256 = hashlib.sha256(coefficient_payload.encode()).hexdigest()
check(
    coefficient_sha256 == "e5d793f30b0db3def23bab3f1ed08ab247b98fb01ab91477099f757b4c623848",
    "independent R14 coefficient ledger changed",
)

# Three modular controls use a different specialization and primes from the
# primary route.  Nonzero constant terms make the p-adic resultant valuation
# exactly two, while squarefreeness guards the residual point count.
modular_controls: list[str] = []
R14_r1 = sy.Poly(R14.as_expr().subs(r, 1), p)
clear_denominator = sy.lcm([sy.denom(value) for value in R14_r1.all_coeffs()])
for prime in (103, 107, 109):
    residual_mod_prime = sy.Poly(
        sy.expand(clear_denominator * R14_r1.as_expr()), p, modulus=prime
    )
    check(residual_mod_prime.degree() == 14, f"degree drop modulo {prime}")
    check(residual_mod_prime.TC() != 0, f"p-valuation exceeds two modulo {prime}")
    check(
        sy.gcd(residual_mod_prime, residual_mod_prime.diff()).degree() == 0,
        f"residual is not squarefree modulo {prime}",
    )
    modular_controls.append(f"{prime}:v_p=2,sqfree14")

# p^2 is a projection artefact: the polynomial pair A,B has no common root
# at p=0.  The rational gradient before this division has two omitted points.
check(sy.factor(A.subs(p, 0)) == -s, "wrong A at p=0")
check(sy.factor(B.subs(p, 0)) == -6, "wrong B at p=0")
check(sy.factor(Gs.subs(p, 0)) == 0, "G_s should vanish identically at p=0")
check(
    sy.factor(2 * Gp.subs(p, 0) - (1 / s**2 - 6)) == 0,
    "wrong direct p=0 equation",
)
check(sy.factor(Gsp.subs(p, 0) - sy.Rational(1, 2)) == 0, "wrong p=0 value")

# The birational chart omits t=0.  Reconstructing the polynomial chart gives
# exactly X^2=-6 there, independently restoring two further points.
X, T = sy.symbols("X T")
P = T + X**2 * T**2
Gxt = sy.cancel(Gsp.subs({s: X * T, p: P}))
check(
    not sy.denom(Gxt).has(X, T),
    "reconstructed (X,T) expression is not polynomial over Q(r)",
)
GX = sy.diff(Gxt, X)
GT = sy.diff(Gxt, T)
check(sy.factor(GX.subs(T, 0)) == 0, "wrong X-gradient at T=0")
check(
    sy.factor(GT.subs(T, 0) + (X**2 + 6) / 2) == 0,
    "wrong T=0 stratum",
)
check(sy.factor(Gxt.subs(T, 0)) == 0, "wrong T=0 critical value")

# At t=0 the independent polynomial A=-s forces s=0, after which B=-6;
# hence no R14 solution was silently counted on the omitted divisor.
check(sy.factor(A.subs(p, s**2)) == -s, "t=0 A specialization changed")
check(sy.factor(B.subs({p: s**2, s: 0})) == -6, "t=0 false common root")

semantic_lines = (
    "scope=delta_only_exact_M8_DeltaD_wall",
    "parameter=delta=5696/(15(6r^2+7));kappa=delta*r^2;phi=2delta*r;r!=0;6r^2+7!=0",
    "critical_projection_residual_degree=14",
    "universal_t=-1/6:length=2",
    "omitted_t=0:length=2",
    "affine_critical_length=18",
    "critical_values=target_nodes_only_under_Keller",
    "verdict=FINITE_CRITICAL_LEDGER_ONLY;DeltaD_wall=OPEN;JC2=OPEN",
)
semantic_sha256 = hashlib.sha256(("\n".join(semantic_lines) + "\n").encode()).hexdigest()

print("THM4140_DELTA_D_CRITICAL_INDEPENDENT")
print("coordinates=s=XT;p=T+s^2;t=p-s^2")
print("gradient_pair=A_degree_s_5;B_degree_s_2")
print("critical_resultant=p^2*R14;degree_R14=14")
print("R14_monic_coefficient_sha256=" + coefficient_sha256)
print("positive_controls=" + ";".join(modular_controls))
print("p=0_direct=s^2=1/6,length=2,value=1/2")
print("t=0_restored=X^2=-6,length=2,value=0")
print("semantic_sha256=" + semantic_sha256)
print("verdict=affine_critical_length=18;DeltaD_wall=OPEN;JC2=OPEN")
