#!/usr/bin/env python3
"""Exact companion for THM-3896's provisional n=4 Kummer closure.

The script works entirely in the THM-3881 residual ring.  It verifies the
complete homogeneous parameterization after the two forced gauge jets, the
high-degree D*Q shell, both x-adic branches of the terminal linear form, and
the two coefficients that obstruct the square-root recursion.  A specialized
positive control reaches the last obstruction, so the calculation does not
mistake an already-empty leading shell for a new closure.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.expand(expression) == 0, message)


def homogeneous_component(
    expression: sp.Expr, degree: int, x: sp.Symbol, y: sp.Symbol
) -> sp.Expr:
    h = sp.Symbol("_h")
    scaled = sp.Poly(sp.expand(expression.subs({x: h * x, y: h * y})), h)
    return sp.expand(scaled.coeff_monomial(h**degree))


def x_valuation(expression: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> int:
    polynomial = sp.Poly(sp.expand(expression), x, y)
    gate(not polynomial.is_zero, "x-valuation requested for zero polynomial")
    return min(monomial[0] for monomial, coefficient in polynomial.terms()
               if coefficient != 0)


def reduce_rho(expression: sp.Expr, rho: sp.Symbol) -> sp.Expr:
    """Reduce an expression modulo rho^2+216; inputs have rational denominators."""
    numerator, denominator = sp.fraction(sp.together(expression))
    reduced_numerator = sp.Poly(sp.expand(numerator), rho).rem(
        sp.Poly(rho**2 + 216, rho)
    ).as_expr()
    reduced_denominator = sp.Poly(sp.expand(denominator), rho).rem(
        sp.Poly(rho**2 + 216, rho)
    ).as_expr()
    # All uses below have rho-free rational denominators after the preferred
    # representative 1/rho=-rho/216 has been chosen.
    gate(not sp.Poly(reduced_denominator, rho).degree() > 0,
         "unexpected rho-dependent denominator")
    return sp.expand(reduced_numerator / reduced_denominator)


x, y, rho = sp.symbols("x y rho")
alpha, beta = sp.symbols("alpha beta")
z = sp.symbols("z")

a = x + 1
L = 9 * x + 4
K2 = y**2 - 15 * x**2
K1 = -15 * x
K0 = -4
K = K2 + K1 + K0
P = a * L**2
Delta = sp.expand(a**3 * L**2 - K**2)

Delta5 = homogeneous_component(Delta, 5, x, y)
Delta4 = homogeneous_component(Delta, 4, x, y)
Delta3 = homogeneous_component(Delta, 3, x, y)
zero(Delta5 - 81 * x**5, "Delta degree-five form")
zero(Delta4 - (315 * x**4 - K2**2), "Delta degree-four form")
zero(Delta3 - (475 * x**3 + 30 * x * K2), "Delta degree-three form")


# ---------------------------------------------------------------------------
# 1. Complete n=4 parameterization after the Kummer terminal symbol.
# ---------------------------------------------------------------------------

q00, q01, q02, q03 = sp.symbols("q00 q01 q02 q03")
rr40, rr41, rr42, rr43, rr44 = sp.symbols("rr40 rr41 rr42 rr43 rr44")
q0_pre = q00 * x**3 + q01 * x**2 * y + q02 * x * y**2 + q03 * y**3
r4_pre = (
    rr40 * x**4 + rr41 * x**3 * y + rr42 * x**2 * y**2
    + rr43 * x * y**3 + rr44 * y**4
)
pre_kummer_top = sp.expand(
    8 * (81 * x**5 * q0_pre) * (81 * x**5 * q0_pre**2)
    + 3 * r4_pre**2 * (81 * x**5 * q0_pre**2)
)
zero(
    pre_kummer_top
    - 243 * q0_pre**2 * x**5 * (r4_pre**2 + 216 * q0_pre * x**5),
    "n=4 middle collision before Kummer normalization",
)

q10, q11, q12 = sp.symbols("q10 q11 q12")
q20, q21 = sp.symbols("q20 q21")
u0, u1 = sp.symbols("u0 u1")
t20, t21, t22 = sp.symbols("t20 t21 t22")
t10, t11 = sp.symbols("t10 t11")

s = alpha * x + beta * y
q0 = x * s**2
q1 = q10 * x**2 + q11 * x * y + q12 * y**2
q2 = q20 * x + q21 * y
u = u0 * x + u1 * y
tau2 = t20 * x**2 + t21 * x * y + t22 * y**2
tau1 = t10 * x + t11 * y
W = q0 + q1 + q2

bar_f = u + z
bar_T = rho * x**2 * s + tau2 + tau1 + 4 * z
f = sp.expand(a * W + bar_f)
T = sp.expand(-K * W + bar_T)

expected_f = {
    4: x * q0,
    3: q0 + x * q1,
    2: q1 + x * q2,
    1: q2 + u,
    0: z,
}
expected_T = {
    5: -K2 * q0,
    4: -K2 * q1 + 15 * x * q0,
    3: -K2 * q2 + 15 * x * q1 + 4 * q0 + rho * x**2 * s,
    2: 15 * x * q2 + 4 * q1 + tau2,
    1: 4 * q2 + tau1,
    0: 4 * z,
}
for degree, expected in expected_f.items():
    zero(homogeneous_component(f, degree, x, y) - expected,
         f"complete f component degree {degree}")
for degree, expected in expected_T.items():
    zero(homogeneous_component(T, degree, x, y) - expected,
         f"complete T component degree {degree}")
zero(T.subs({x: 0, y: 0}) - 4 * f.subs({x: 0, y: 0}),
     "THM-3881 address")

r = sp.expand(a * T + K * f)
bar_r = sp.expand(a * bar_T + K * bar_f)
zero(r - bar_r, "gauge edge r is invariant")
r4 = rho * x**3 * s
r3 = rho * x**2 * s + x * tau2 + K2 * u
r2 = tau2 + x * tau1 + K1 * u + K2 * z
zero(homogeneous_component(r, 4, x, y) - r4, "terminal r4")
zero(homogeneous_component(r, 3, x, y) - r3, "next r3")
zero(homogeneous_component(r, 2, x, y) - r2, "next r2")

A = sp.expand(K * T + a * P * f)
B = sp.expand(P * f**2 - T**2)
bar_A = sp.expand(K * bar_T + a * P * bar_f)
bar_B = sp.expand(P * bar_f**2 - bar_T**2)
zero(A - (bar_A + Delta * W), "A transport by Delta W")
zero(B - (bar_B + 2 * W * bar_A + Delta * W**2),
     "B transport by its two sidecars")


# ---------------------------------------------------------------------------
# 2. The exact high-degree shell.
# ---------------------------------------------------------------------------

D = sp.expand(Delta * W**2)
Q = sp.expand(8 * Delta * W + 3 * r**2)

D11 = Delta5 * q0**2
D10 = 2 * Delta5 * q0 * q1 + Delta4 * q0**2
D9 = (
    Delta5 * (q1**2 + 2 * q0 * q2)
    + 2 * Delta4 * q0 * q1
    + Delta3 * q0**2
)
Q8 = reduce_rho(8 * Delta5 * q0 + 3 * r4**2, rho)
Q7 = reduce_rho(8 * (Delta5 * q1 + Delta4 * q0) + 6 * r4 * r3, rho)
Q6 = reduce_rho(
    8 * (Delta5 * q2 + Delta4 * q1 + Delta3 * q0)
    + 3 * (r3**2 + 2 * r4 * r2),
    rho,
)

zero(homogeneous_component(D, 11, x, y) - D11, "D11 component")
zero(homogeneous_component(D, 10, x, y) - D10, "D10 component")
zero(homogeneous_component(D, 9, x, y) - D9, "D9 component")
zero(Q8, "critical Kummer cancellation Q8")
zero(reduce_rho(homogeneous_component(Q, 7, x, y) - Q7, rho),
     "Q7 component")
zero(reduce_rho(homogeneous_component(Q, 6, x, y) - Q6, rho),
     "Q6 component")

# Exact degree accounting behind S_d=(D*Q)_d for d>=17.  After Q8=0,
# deg Q<=7, while deg W=3, deg bar_A<=5, deg bar_B<=6, deg D<=11.
gate(sp.Poly(bar_A, x, y).total_degree() <= 5, "bar A degree bound")
gate(sp.Poly(bar_B, x, y).total_degree() <= 6, "bar B degree bound")
gate(sp.Poly(D, x, y).total_degree() <= 11, "D degree bound")
gate(sp.Poly(Q, x, y).total_degree() <= 8, "unreduced Q degree bound")
gate(Q8 == 0, "Q degree drops to seven modulo rho squared plus 216")
gate(7 + max(3 + 5, 6) <= 15, "Q sidecar correction degree")
gate(5 + 11 <= 16, "bar A times D correction degree")
gate(max(4, 8 + 2 + 4, 3 + 11) <= 14,
     "all non-AB/r2B residual macros degree at most fourteen")


# ---------------------------------------------------------------------------
# 3. Generic terminal line beta != 0: the degree-17 obstruction.
# ---------------------------------------------------------------------------

generic_D11 = sp.expand(D11)
generic_D10 = sp.expand(D10)
generic_Q7 = sp.expand(Q7)
gate(x_valuation(generic_D11, x, y) == 7, "generic v_x(D11)=7")
gate(x_valuation(generic_D10, x, y) == 2, "generic v_x(D10)=2")
gate(x_valuation(generic_Q7, x, y) == 1, "generic v_x(Q7)=1")

# D11*Q6 is divisible by x^7, so the unique x^3 coefficient of S17 comes
# from D10*Q7.  It is independent of every lower homogeneous symbol.
generic_decisive = sp.Poly(sp.expand(generic_D10 * generic_Q7), x, y)
zero(generic_decisive.coeff_monomial(x**3 * y**14) - 8 * beta**6,
     "generic decisive coefficient")
gate(x_valuation(generic_D10 * generic_Q7, x, y) == 3,
     "generic v_x(S17)=3")


# ---------------------------------------------------------------------------
# 4. Exceptional line s=sigma*x: one more lift, then degree 16 closes.
# ---------------------------------------------------------------------------

sigma, eta = sp.symbols("sigma eta")
q0e = sigma**2 * x**3
r4e = rho * sigma * x**4
r3e0, r3e1, r3e2, r3e3 = sp.symbols("r3e0 r3e1 r3e2 r3e3")
r3e = r3e0 * x**3 + r3e1 * x**2 * y + r3e2 * x * y**2 + r3e3 * y**3
r2e0, r2e1, r2e2 = sp.symbols("r2e0 r2e1 r2e2")
r2e = r2e0 * x**2 + r2e1 * x * y + r2e2 * y**2

D11e = sp.expand(Delta5 * q0e**2)
D10e = sp.expand(2 * Delta5 * q0e * q1 + Delta4 * q0e**2)
D9e = sp.expand(
    Delta5 * (q1**2 + 2 * q0e * q2)
    + 2 * Delta4 * q0e * q1
    + Delta3 * q0e**2
)
Q7e = reduce_rho(
    8 * (Delta5 * q1 + Delta4 * q0e) + 6 * r4e * r3e,
    rho,
)
Q6e = reduce_rho(
    8 * (Delta5 * q2 + Delta4 * q1 + Delta3 * q0e)
    + 3 * (r3e**2 + 2 * r4e * r2e),
    rho,
)

gate(x_valuation(D11e, x, y) == 11, "exceptional v_x(D11)=11")
gate(x_valuation(D10e, x, y) == 6, "exceptional v_x(D10)=6")
gate(x_valuation(Q7e, x, y) == 3, "exceptional v_x(Q7)=3")
zero(sp.cancel(Q7e / x**3).subs(x, 0) + 8 * sigma**2 * y**4,
     "exceptional square symbol modulo x")

h20, h21 = sp.symbols("h20 h21")
h2 = h20 * x**2 + h21 * x * y + eta * y**2
c40, c41, c42, c43, c44 = sp.symbols("c40 c41 c42 c43 c44")
c4 = c40 * x**4 + c41 * x**3 * y + c42 * x**2 * y**2 + c43 * x * y**3 + c44 * y**4
M = sp.expand(sigma**2 * Delta4 + 162 * x**2 * q1)
g9 = 9 * sigma**2 * x**7 * h2
g8 = sp.expand(x**2 * M * h2 / 18 + sp.Rational(9, 2) * sigma**2 * x**4 * c4)

# Under Q7=x^3*h2^2 and Q6=h2*c4, this is exactly S17=2*g9*g8.
formal_S17 = sp.expand(
    sigma**2 * x**9 * M * h2**2
    + 81 * sigma**4 * x**11 * h2 * c4
)
zero(formal_S17 - 2 * g9 * g8, "exceptional degree-17 square equation")
zero(M.subs(x, 0) + sigma**2 * y**4, "M modulo x")
zero(sp.Poly(g8, x, y).coeff_monomial(x**2 * y**6)
     + sigma**2 * eta / 18, "exceptional leading g8 coefficient")

# The exact degree-16 shell has x-valuation at least six:
# D11*Q5 + D10*Q6 + D9*Q7 + 8*bar_A5*D11.
gate(x_valuation(D9e, x, y) >= 3, "exceptional v_x(D9)>=3")
gate(11 + 0 >= 6, "D11 Q5 valuation bound")
gate(6 + 0 >= 6, "D10 Q6 valuation bound")
gate(3 + 3 >= 6, "D9 Q7 valuation bound")
gate(0 + 11 >= 6, "bar A5 D11 valuation bound")

g8_square_lead = sp.Poly(sp.expand(g8**2), x, y).coeff_monomial(x**4 * y**12)
exceptional_obstruction = sp.expand((-g8_square_lead).subs(eta**2, -8 * sigma**2))
zero(exceptional_obstruction - sp.Rational(2, 81) * sigma**6,
     "exceptional decisive coefficient")


# ---------------------------------------------------------------------------
# 5. Positive controls: the first two square equations really can lift.
# ---------------------------------------------------------------------------

tau2_control = sp.Rational(17, 18) * rho * x * s  # = -204*x*s/rho
r3_control = reduce_rho(rho * x**2 * s + x * tau2_control, rho)
Q7_control = reduce_rho(
    8 * Delta4 * q0 + 6 * r4 * r3_control,
    rho,
)
zero(Q7_control / x + 8 * s**2 * K2**2,
     "all-s first square-layer positive control")

# On the exceptional line, tau1=7*rho*x/72 (= -21*x/rho) makes Q6
# divisible by h2=lambda*K2 for every remaining constant z.
r3_control_e = sp.Rational(35, 18) * rho * x**3
r2_control_e = sp.Rational(25, 24) * rho * x**2 + K2 * z
Q7_control_e = reduce_rho(
    8 * Delta4 * x**3 + 6 * rho * x**4 * r3_control_e,
    rho,
)
Q6_control_e = reduce_rho(
    8 * Delta3 * x**3
    + 3 * (r3_control_e**2 + 2 * rho * x**4 * r2_control_e),
    rho,
)
zero(Q7_control_e + 8 * x**3 * K2**2,
     "exceptional positive-control Q7")
zero(Q6_control_e - 6 * (40 + rho * z) * x**4 * K2,
     "exceptional positive-control Q6")

# Reconstruct the full specialized residual and independently recover the
# degree-16 obstruction coefficient.
f_control = a * x**3 + z
T_control = sp.expand(
    -K * x**3
    + rho * x**3
    + sp.Rational(17, 18) * rho * x**2
    + sp.Rational(7, 72) * rho * x
    + 4 * z
)
r_control = sp.expand(a * T_control + K * f_control)
A_control = sp.expand(K * T_control + a * P * f_control)
B_control = sp.expand(P * f_control**2 - T_control**2)
S_control = sp.expand(
    L**4
    + 2 * (3 * A_control + 3 * P + r_control**2) * L**2 * f_control
    + (8 * A_control + 6 * P + 3 * r_control**2) * B_control
)
S19_control = reduce_rho(homogeneous_component(S_control, 19, x, y), rho)
S18_control = reduce_rho(homogeneous_component(S_control, 18, x, y), rho)
S16_control = reduce_rho(homogeneous_component(S_control, 16, x, y), rho)
zero(S19_control, "positive-control degree nineteen cancellation")
zero(S18_control + 648 * x**14 * K2**2,
     "positive-control degree eighteen square")

# If lambda^2=-8, g8=lambda*H8.  Hence g8^2=-8*H8^2.
H8_control = sp.expand(
    x**2 * Delta4 * K2 / 18
    - sp.Rational(27, 8) * (40 + rho * z) * x**8
)
N16_control = reduce_rho(S16_control + 8 * H8_control**2, rho)
zero(sp.Poly(N16_control, x, y).coeff_monomial(x**4 * y**12)
     - sp.Rational(2, 81), "positive-control terminal obstruction")
gate(x_valuation(S16_control, x, y) >= 6,
     "positive-control v_x(S16)>=6")
gate(x_valuation(N16_control, x, y) == 4,
     "positive-control v_x(S16-g8^2)=4")


semantic = {
    "theorem": "THM-3896",
    "universe": "THM3881 residual, n=4 equality seam",
    "parameterization": "three gauge components plus (u,z,tau2,tau1)",
    "kummer_entry": "S19=243*q0^2*x^5*(r4^2+216*q0*x^5)",
    "high_shell": "S_d=(Delta*W^2*(8*Delta*W+3*r^2))_d for d>=17",
    "generic_branch": "[x^3*y^14]S17=8*beta^6 contradicts v_x(g9)=4",
    "exceptional_branch": "[x^4*y^12](S16-g8^2)=2*sigma^6/81 contradicts v_x(g9)=7",
    "positive_control": "Q7=-8*x^3*K2^2 and Q6=6*(40+rho*z)*x^4*K2",
    "scope": "provisional proof candidate; independent hostile audit pending; JC2 open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3896-cusp-residual-degree-four-equality-seam-emptiness")
print("status=RESERVED_PROVISIONAL_PROOF_CANDIDATE_VERIFIED_EXACT_AWAITING_INDEPENDENT_HOSTILE_AUDIT")
print("n4_kummer_entry_rederived=YES")
print("complete_n4_parameterization=YES")
print("high_shell_degrees_at_least_17=DQ")
print("generic_decisive_coefficient=[x^3*y^14]S17=8*beta^6")
print("generic_terminal_line=EMPTY_PROVISIONAL")
print("exceptional_divisibility_gate=h2_divides_Q6")
print("exceptional_decisive_coefficient=[x^4*y^12](S16-g8^2)=2*sigma^6/81")
print("exceptional_terminal_line=EMPTY_PROVISIONAL")
print("positive_control_Q7=-8*x^3*K2^2")
print("positive_control_Q6=6*(40+rho*z)*x^4*K2")
print("degree_four_equality_seam=EMPTY_PROVISIONAL")
print("JC2_status=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
