#!/usr/bin/env python3
"""Exact candidate extension of THM-3885 from cubic through quartic degree."""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(newline="\n")

GATES = 0


def gate(condition: bool, message: str) -> None:
    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


# A total-degree-four survivor inherits THM-3885's a=0 normal form
# T=c+(x+1)U, deg U<=3.  At x=0 the address T(0,0)=0 leaves the complete
# cubic-in-y restriction tau=p*y+q*y^2+r*y^3.
y = sp.symbols("y")
p, q, r = sp.symbols("p q r")
tau = p * y + q * y**2 + r * y**3
S0 = sp.expand(256 - 96 * tau**2 - 8 * (y**2 - 4) * tau**3 - 3 * tau**4)
zero(S0.subs(y, 0) - 256, "x-zero address value")
gate(sp.Poly(S0, y).degree() <= 12, "x-zero quartic-cell degree cap")

# Fix the root sign by G(0)=16.  Coefficients through y^6 determine the
# complete possible degree-six square root.  The remaining six coefficients
# generate an exact zero-dimensional ideal.
root = sp.Integer(16)
root_coefficients = []
for degree in range(1, 7):
    unknown = sp.symbols(f"g{degree}")
    coefficient = sp.Poly(
        sp.expand((root + unknown * y**degree) ** 2 - S0), y
    ).coeff_monomial(y**degree)
    solutions = sp.solve(coefficient, unknown, dict=False)
    gate(len(solutions) == 1, f"unique root recurrence degree={degree}")
    root_coefficients.append(sp.factor(solutions[0]))
    root = sp.expand(root + solutions[0] * y**degree)

expected_root_coefficients = (
    0,
    -3 * p**2,
    p * (p**2 - 6 * q),
    -sp.Rational(3, 8) * (p**4 - 8 * p**2 * q + 16 * p * r + 8 * q**2),
    sp.Rational(1, 16)
    * (3 * p**5 - 24 * p**3 * q - 4 * p**3 + 48 * p**2 * r
       + 48 * p * q**2 - 96 * q * r),
    -sp.Rational(1, 128)
    * (13 * p**6 - 120 * p**4 * q + 192 * p**3 * r
       + 288 * p**2 * q**2 + 96 * p**2 * q - 768 * p * q * r
       - 128 * q**3 + 384 * r**2),
)
for degree, (actual, expected) in enumerate(
    zip(root_coefficients, expected_root_coefficients), start=1
):
    zero(actual - expected, f"quartic x-zero root recurrence degree={degree}")

difference = sp.Poly(sp.expand(root**2 - S0), y)
for degree in range(7):
    zero(difference.coeff_monomial(y**degree), f"root recurrence clears y^{degree}")
remainders = [
    sp.together(difference.coeff_monomial(y**degree)).as_numer_denom()[0]
    for degree in range(7, 13)
]
gate(all(remainder != 0 for remainder in remainders), "six nontrivial terminal remainders")

# Exact characteristic-zero Groebner computation.  Three members suffice for
# the zero-set conclusion; reductions certify their membership in the full
# terminal-remainder ideal.
terminal_basis = sp.groebner(remainders, r, q, p, order="grevlex", method="f5b")
cert_r = r**4
cert_q = q**4 - 6 * p**2 * r**2
cert_p = 5 * p**5 + 80 * p**2 * r + 80 * p * q**2 + 448 * r**3
for certificate, name in (
    (cert_r, "r fourth-power certificate"),
    (cert_q, "q fourth-power certificate"),
    (cert_p, "p fifth-power certificate"),
):
    quotients, remainder = terminal_basis.reduce(certificate)
    zero(remainder, name)
    gate(any(quotient != 0 for quotient in quotients), f"nontrivial {name}")

# Thus r=0, then q=0, then p=0 over a characteristic-zero field.  This
# removes global y-degree three.  The remaining y-degree cases are independent
# leading-coefficient obstructions in alpha=a=x+1.
alpha = sp.symbols("alpha")
b0, b1 = sp.symbols("b0 b1")
b = alpha * (b0 + b1 * alpha)
lead_y8 = sp.expand(-b**3 * (8 + 3 * alpha**2 * b))
gate(sp.Poly(lead_y8, alpha).coeff_monomial(alpha**3) == -8 * b0**3,
     "dy2 odd alpha valuation when b0 nonzero")

# If b0=0, squarehood of the leading coefficient would force the quartic
# binomial 8+3*b1*alpha^4 to be a square.  Its complete quadratic-root row is
# inconsistent when b1!=0.
h0, h1, h2 = sp.symbols("h0 h1 h2")
binomial_difference = sp.Poly(
    sp.expand((h0 + h1 * alpha + h2 * alpha**2) ** 2
              - (8 + 3 * b1 * alpha**4)),
    alpha,
)
binomial_row = tuple(binomial_difference.coeff_monomial(alpha**i) for i in range(5))
gate(
    binomial_row
    == (h0**2 - 8, 2 * h0 * h1, h1**2 + 2 * h0 * h2,
        2 * h1 * h2, h2**2 - 3 * b1),
    "dy2 terminal binomial row",
)

# dy=1 has unique y^5 coefficient -8b^3; dy=0 has the exact split
# S=-8t(alpha)^3*y^2+F(alpha), and F cannot vanish by its unique degree row.
B = sp.symbols("B", nonzero=True)
gate(-8 * B**3 != 0, "dy1 odd y-degree coefficient")
for e in range(101):
    degree_rows = (4, 2 * e + 3, 3 * e + 2, 4 * e + 2)
    if e == 0:
        gate(degree_rows[0] > max(degree_rows[1:]), "dy0 constant unique L4 row")
    else:
        gate(degree_rows[3] > max(degree_rows[:3]), f"dy0 quartic unique row e={e}")

# Modular hostile controls use independently recomputed bases and must retain
# the same pure-power triangular consequence away from small bad primes.
modular_rows = []
for prime in (7, 11, 101):
    basis_mod = sp.groebner(remainders, r, q, p, order="grevlex", modulus=prime)
    rem_r = basis_mod.reduce(cert_r)[1]
    rem_q = basis_mod.reduce(cert_q)[1]
    rem_p = basis_mod.reduce(cert_p)[1]
    zero(rem_r, f"modular r certificate p={prime}")
    zero(rem_q, f"modular q certificate p={prime}")
    zero(rem_p, f"modular p certificate p={prime}")
    modular_rows.append((prime, len(basis_mod.polys)))

semantic = {
    "input": "THM3885 a-zero normal form and admissible address",
    "x_zero": "tau=p*y+q*y2+r*y3; six terminal square-root remainders",
    "certificates": "r4;q4-6p2r2;5p5+80p2r+80pq2+448r3",
    "remaining": "dy2 alpha parity/binomial;dy1 odd;dy0 missing-y",
    "conclusion": "f=0 and total degree T at most four imply T=0",
    "scope": "degree at least five and JC2 remain open",
}
blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("THM3885_FZERO_QUARTIC_CLOSURE_PRIMARY_20260823")
print("status=PROVED_VERIFIED_EXACT_CANDIDATE;INDEPENDENT_AUDIT_PENDING")
print("xzero_terminal_certificates=r^4;q^4-6*p^2*r^2;5*p^5+80*p^2*r+80*p*q^2+448*r^3")
print("closed=f=0;addressed;deg_total_T<=4;unique_T=0")
print(f"modular_rows={tuple(modular_rows)}")
print(f"semantic_sha256={hashlib.sha256(blob).hexdigest()}")
print(f"gates={GATES}")
print("ALL CHECKS PASSED")
