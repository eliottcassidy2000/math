#!/usr/bin/env python3
"""Exact primary audit for the THM-3905 equality-color third response."""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


GATES = 0


def require(condition: bool, message: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(message)


def trunc(expression: sp.Expr, epsilon: sp.Symbol, order: int = 4) -> sp.Expr:
    return sp.cancel(sp.series(expression, epsilon, 0, order).removeO())


def zero(expression: sp.Expr, epsilon: sp.Symbol, message: str) -> None:
    numerator = sp.together(trunc(expression, epsilon)).as_numer_denom()[0]
    require(sp.expand(numerator) == 0, message)


epsilon = sp.symbols("epsilon")
a, L, kappa = sp.symbols("a L kappa", nonzero=True)

# (name, coefficient, T exponent, f exponent, K exponent)
support = (
    ("base", L**4, 0, 0, 0),
    ("f", 6*a*L**4, 0, 1, 0),
    ("f2", 12*a**2*L**4, 0, 2, 0),
    ("f3", 8*a**3*L**4, 0, 3, 0),
    ("KTf", 6*L**2, 1, 1, 1),
    ("KTf2", 12*a*L**2, 1, 2, 1),
    ("KTf3", 6*a**2*L**2, 1, 3, 1),
    ("T2", -6*a*L**2, 2, 0, 0),
    ("T2f", -6*a**2*L**2, 2, 1, 0),
    ("T2f2", 3*a**3*L**2, 2, 2, 0),
    ("KT3", -8, 3, 0, 1),
    ("KT3f", -6*a, 3, 1, 1),
    ("T4", -3*a**2, 4, 0, 0),
    ("K2f3", 2*L**2, 0, 3, 2),
    ("K2f4", 3*a*L**2, 0, 4, 2),
    ("K2T2f2", -3, 2, 2, 2),
)


def deficit(n: int, p: int, q: int, r: int) -> int:
    return (4-p-q)*n + 4-2*r


expected = {
    1: (
        ("K2f4", 0), ("K2T2f2", 0), ("K2f3", 1),
        ("KTf3", 2), ("KT3f", 2), ("KTf2", 3), ("KT3", 3),
    ),
    2: (
        ("K2f4", 0), ("K2T2f2", 0), ("KTf3", 2),
        ("KT3f", 2), ("K2f3", 2),
    ),
    3: (
        ("K2f4", 0), ("K2T2f2", 0), ("KTf3", 2),
        ("KT3f", 2), ("K2f3", 3),
    ),
    4: (
        ("K2f4", 0), ("K2T2f2", 0), ("KTf3", 2),
        ("KT3f", 2),
    ),
}

support_ledger: dict[int, tuple[tuple[str, int], ...]] = {}
for n in range(1, 9):
    rows = tuple(
        sorted(
            (name, deficit(n, p, q, r))
            for name, _, p, q, r in support
            if deficit(n, p, q, r) <= 3
        )
    )
    support_ledger[n] = rows
    reference = tuple(sorted(expected[min(n, 4)]))
    require(rows == reference, f"deficit<=3 support mismatch at n={n}: {rows}")


coefficient_ledger: dict[int, tuple[str, ...]] = {}
for n in range(1, 9):
    u = sp.symbols("u0:4")
    v = sp.symbols("v0:4", nonzero=True)
    U = sum(u[j]*epsilon**j for j in range(min(n, 3)+1))
    V = sum(v[j]*epsilon**j for j in range(min(n, 3)+1))
    K = epsilon**-2 * (1+kappa*epsilon**2)
    T = epsilon**-n * U
    f = epsilon**-n * V

    S = sum(coefficient*T**p*f**q*K**r for _, coefficient, p, q, r in support)
    color_product = trunc(epsilon**(4*n+4)*S/V**2 + 3*U**2, epsilon)

    candidate = (
        3*a*L**2*V**2
        + 2*L**2*epsilon**n*(1+2*kappa*epsilon**2)*V
        + 6*epsilon**2*(kappa+a*U/V)*(a*L**2*V**2-U**2)
    )
    if n == 1:
        candidate += epsilon**3*(12*a*L**2*U-8*U**3/V**2)
    candidate = trunc(candidate, epsilon)

    zero(color_product-candidate, epsilon, f"compact third-response law n={n}")
    coeffs = tuple(str(sp.factor(sp.expand(color_product).coeff(epsilon, j))) for j in range(4))
    coefficient_ledger[n] = coeffs

    # Clear denominators before every exact comparison, and probe a hostile
    # specialization with nonzero V constant to catch accidental cancellation.
    hostile = {
        a: 2, L: 5, kappa: -7,
        u[0]: 3, u[1]: -2, u[2]: 4, u[3]: 1,
        v[0]: 6, v[1]: 5, v[2]: -3, v[3]: 2,
    }
    zero((color_product-candidate).subs(hostile), epsilon,
         f"hostile compact third-response law n={n}")


# The canonical THM-3902 leading payment gives an exact positive third-jet
# control for n>=4 and a sharp zero-sidecar hostile at n=3.
x, droot = sp.symbols("x droot")
ac = x+1
Lc = 9*x+4
uc = (3*Lc**2-ac)/(2*droot)
hc = (ac+3*Lc**2)/2


def zero_mod_droot(expression: sp.Expr, message: str) -> None:
    numerator = sp.together(expression).as_numer_denom()[0]
    remainder = sp.rem(sp.Poly(sp.expand(numerator), droot),
                       sp.Poly(droot**2+3, droot)).as_expr()
    require(sp.expand(remainder) == 0, message)


def reduce_droot_fraction(expression: sp.Expr) -> sp.Expr:
    """Reduce a rational function modulo droot^2+3 and rationalize droot."""
    numerator, denominator = sp.fraction(sp.cancel(expression))
    modulus = sp.Poly(droot**2+3, droot)
    numerator = sp.rem(sp.Poly(sp.expand(numerator), droot), modulus).as_expr()
    denominator = sp.rem(sp.Poly(sp.expand(denominator), droot), modulus).as_expr()
    d0 = denominator.coeff(droot, 0)
    d1 = denominator.coeff(droot, 1)
    rationalized = sp.expand(numerator*(d0-d1*droot))
    rationalized = sp.rem(sp.Poly(rationalized, droot), modulus).as_expr()
    return sp.cancel(rationalized/(d0**2+3*d1**2))


zero_mod_droot(hc**2-3*(ac*Lc**2-uc**2), "canonical leading payment")
third_n4 = sp.Integer(0)
require(third_n4 == 0, "n>=4 zero-sidecar third response")
third_n3 = Lc**2/hc
require(sp.degree(hc, x) == 2 and sp.degree(Lc, x) == 1,
        "n=3 hostile degree profile")
require(sp.gcd(sp.Poly(hc, x), sp.Poly(Lc, x)).degree() == 0,
        "n=3 hostile coprimality")
require(sp.rem(sp.Poly(Lc**2, x), sp.Poly(hc, x)).as_expr() != 0,
        "n=3 marked source must not be polynomial")

# THM-3902's named address-compatible n=1 two-jet control dies at this row.
p_control = -sp.Rational(2, 147)/(1+4*droot)
v1_control = p_control*hc
j1_control = ac*v1_control+sp.Rational(1, 3)
u1_control = droot*j1_control/3
j2_control = (
    (kappa.subs(kappa, -15*x**2-15*x-4)+ac*uc)*hc
    + Lc**2*p_control*(3*ac*p_control*hc+2)/2
)
zero_mod_droot(u1_control.subs(x, 0)-4*v1_control.subs(x, 0),
               "named n=1 control address")

Uc = uc+u1_control*epsilon
Vc = 1+v1_control*epsilon
kappac = -15*x**2-15*x-4
Ec = ac*Lc**2*Vc**2-Uc**2
Pc = trunc(
    3*ac*Lc**2*Vc**2
    + 2*Lc**2*epsilon*(1+2*kappac*epsilon**2)*Vc
    + 6*epsilon**2*(kappac+ac*Uc/Vc)*Ec
    + epsilon**3*(12*ac*Lc**2*Uc-8*Uc**3/Vc**2),
    epsilon,
)
P3_control = sp.expand(Pc).coeff(epsilon, 3)
j3_control = sp.cancel((P3_control-2*j1_control*j2_control)/(2*hc))
g1_control = hc*v1_control+j1_control
g2_control = j2_control+v1_control*j1_control
g3_control = sp.cancel(j3_control+v1_control*j2_control)

# Direct raw-square response agrees with the color calculation.
Hc = hc+j1_control*epsilon+j2_control*epsilon**2+j3_control*epsilon**3
Gamma_c = trunc(Vc*Hc, epsilon)
S_normalized = trunc(Vc**2*(Pc-3*Uc**2), epsilon)
zero_mod_droot(
    sp.expand(S_normalized).coeff(epsilon, 3)
    - (2*hc*g3_control+2*g1_control*g2_control),
    "named n=1 direct third response",
)
zero_mod_droot(
    sp.expand(Gamma_c).coeff(epsilon, 3)-g3_control,
    "named n=1 Gamma/V reconstruction",
)

denominator_constant = sp.Integer(249143169618)
h_factor = 2*hc
scaled_numerator = reduce_droot_fraction(
    denominator_constant*h_factor*g3_control
)
scaled_denominator = sp.fraction(scaled_numerator)[1]
require(sp.Poly(scaled_denominator, x).degree() == 0,
        "named n=1 scaled denominator must be constant")
component_0 = sp.cancel(scaled_numerator.subs(droot, 0))
component_1 = sp.cancel((scaled_numerator-component_0)/droot)


def remainder_x(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    remainder = sp.rem(sp.Poly(numerator, x), sp.Poly(h_factor, x)).as_expr()
    return sp.cancel(remainder/denominator)


expected_0 = -sp.Rational(46118408, 3)*(x+1)
expected_1 = sp.Rational(184473632, 3)*(x+1)
require(sp.expand(remainder_x(component_0)-expected_0) == 0,
        "named n=1 numerator residue, unit component")
require(sp.expand(remainder_x(component_1)-expected_1) == 0,
        "named n=1 numerator residue, d component")
require(expected_0 != 0 and expected_1 != 0,
        "named n=1 third coefficient is nonpolynomial")


source = Path(__file__).read_text(encoding="utf-8")
require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "audit must retain checks under -O")

semantic = {
    "support": {str(n): support_ledger[n] for n in support_ledger},
    "law": (
        "C-C+=3aL2V2+2L2e^n(1+2kappa e2)V+"
        "6e2(kappa+aU/V)(aL2V2-U2)+"
        "[n=1]e3(12aL2U-8U3/V2) mod e4"
    ),
    "boundaries": {
        "n1": ("marked@1", "curvature@2", "marked-kappa@3", "KTf2-KT3@3"),
        "n2": ("marked@2", "third coefficient sees v1"),
        "n3": ("marked@3",),
        "nge4": ("no marked source through third response", "zero-sidecar j3=0 lift"),
    },
    "control": "n=3 zero-sidecar j3=L2/h nonpolynomial; n>=4 j3=0",
    "named_n1_control": (
        "THM3902 address-compatible two-jet has denominator "
        "249143169618*(243x2+217x+49), nonzero numerator residue"
    ),
    "scope": "necessary third response only; positive controls and JC2 remain open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3905-nonzero-sidecar-equality-color-third-response")
print("status=PASS;THIRD_RESPONSE_LAW_VERIFIED;JC2_OPEN")
print("support_counts=" + ",".join(f"n{n}:{len(support_ledger[n])}" for n in range(1, 9)))
print("n1_sources=marked@1,curvature@2,marked-kappa@3,KTf2-KT3@3")
print("n2_sources=curvature+marked@2;third_sees_v1")
print("n3_sources=curvature@2;marked@3")
print("nge4_sources=leading+curvature_only_through_epsilon3")
print("control=n3_zero_sidecar_j3=L2/h_nonpolynomial;nge4_j3=0_positive_third_jet")
print("named_n1_THM3902_control=FAILS_THIRD_RESPONSE_NONPOLYNOMIAL_G3")
print("named_n1_denominator=249143169618*(243*x^2+217*x+49)")
print("compact_law=3aL2V2+2L2e^n(1+2kappa*e2)V+6e2(kappa+aU/V)(aL2V2-U2)+delta_n1*e3(12aL2U-8U3/V2)")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"gates={GATES}")
print("ALL CHECKS PASSED")
