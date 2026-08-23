#!/usr/bin/env python3
"""Exact companion for THM-3806's binomial cubic r-repair obstruction."""

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
    gate(sp.cancel(expression) == 0, message)


r, z, e, u = sp.symbols("r z e u")
mu, lam = sp.symbols("mu lambda", nonzero=True)
variables = (r, z, e)
surface = r**2 * e - z**3 + r
poisson = sp.Matrix(
    [
        [0, 3 * r**2, 9 * z**2],
        [-3 * r**2, 0, 3 + 6 * r * e],
        [-9 * z**2, -3 - 6 * r * e, 0],
    ]
)


def bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    dl = sp.Matrix([sp.diff(left, q) for q in variables])
    dr = sp.Matrix([sp.diff(right, q) for q in variables])
    return sp.expand((dl.T * poisson * dr)[0])


# Universal two-equation compression.
G, D = sp.symbols("G D")
K = 1 + 2 * u
P = sp.expand(G * u**2 - K * (2 * e**3 + u * e * D))
Q = sp.expand(e**2 * K**3 - 729 * G**3 * u**2 * (1 + u)**2)
resultant_generic = sp.factor(sp.resultant(P, Q, u))
H_universal = sp.factor(resultant_generic / (G**3 * e**4))
zero(resultant_generic - G**3 * e**4 * H_universal,
     "universal residual resultant")


generic_rows: list[str] = []
exception_rows: list[str] = []
for j in range(3):
    g = mu * e**j + lam * e**3
    gp = sp.diff(g, e)
    H = sp.Poly(sp.expand(H_universal.subs({G: g, D: gp})), e)
    gate(H.degree() == 17, f"generic H degree j={j}")
    zero(H.LC() - 8503056 * lam**3 * (lam - 1)**2,
         f"generic H leading coefficient j={j}")

    quotient, remainder = sp.div(
        sp.Poly(sp.expand(e * g * sp.diff(H.as_expr(), e)), e), H
    )
    # Division is over Q(mu,lambda), localized at lambda-1.  Freeze one
    # quotient identity and the decisive remainder coefficient in each row.
    zero(
        e * g * sp.diff(H.as_expr(), e)
        - quotient.as_expr() * H.as_expr() - remainder.as_expr(),
        f"generic division identity j={j}",
    )
    generic_denominator = sp.Integer(1)
    for coefficient in quotient.all_coeffs() + remainder.all_coeffs():
        generic_denominator = sp.lcm(
            generic_denominator, sp.denom(sp.cancel(coefficient))
        )
    expected_generic_denominators = (lam - 1, 1, 4 * (lam - 1) ** 3)
    zero(
        generic_denominator - expected_generic_denominators[j],
        f"generic division denominator j={j}",
    )
    if j in (0, 1):
        gate(remainder.coeff_monomial(e**6) == 70 * lam**2,
             f"generic e6 witness j={j}")
        generic_rows.append(f"{j}:17:e6:70lambda^2")
    else:
        expected = (
            -1062882 * mu**4 * lam**3 * (lam - 2) / (lam - 1)**2
        )
        zero(remainder.coeff_monomial(e**16) - expected,
             "generic j=2 e16 witness")
        generic_rows.append(
            "2:17:e16:-1062882mu^4lambda^3(lambda-2)/(lambda-1)^2"
        )

    # lambda=1 is a genuine degree drop.  Substitute into H first and then
    # redo polynomial division; never specialize the localized quotient.
    g_one = sp.expand(g.subs(lam, 1))
    H_one = sp.Poly(sp.expand(H.as_expr().subs(lam, 1)), e)
    quotient_one, remainder_one = sp.div(
        sp.Poly(sp.expand(e * g_one * sp.diff(H_one.as_expr(), e)), e),
        H_one,
    )
    zero(
        e * g_one * sp.diff(H_one.as_expr(), e)
        - quotient_one.as_expr() * H_one.as_expr() - remainder_one.as_expr(),
        f"lambda=1 recomputed division j={j}",
    )
    expected_degrees = (11, 10, 15)
    expected_lcs = (2125764 * mu**2, 26244, 2125764 * mu**2)
    gate(H_one.degree() == expected_degrees[j],
         f"lambda=1 degree j={j}")
    zero(H_one.LC() - expected_lcs[j], f"lambda=1 leading coefficient j={j}")
    if j == 0:
        zero(remainder_one.coeff_monomial(e**2) - sp.Rational(361, 81) / mu,
             "lambda=1 j=0 e2 witness")
        exception_rows.append("lambda1:j0:deg11:e2:361/(81mu)")
    elif j == 1:
        gate(remainder_one.coeff_monomial(e**6) == 35,
             "lambda=1 j=1 e6 witness")
        exception_rows.append("lambda1:j1:deg10:e6:35")
    else:
        gate(remainder_one.coeff_monomial(e**4) == 30 * mu**2,
             "lambda=1 j=2 e4 witness")
        exception_rows.append("lambda1:j2:deg15:e4:30mu^2")

    # Only j=2 has a second generic-witness zero at lambda=2.  Recompute
    # there as well; j=0,1 are already handled by 70lambda^2.
    if j == 2:
        g_two = sp.expand(g.subs(lam, 2))
        H_two = sp.Poly(sp.expand(H.as_expr().subs(lam, 2)), e)
        quotient_two, remainder_two = sp.div(
            sp.Poly(sp.expand(e * g_two * sp.diff(H_two.as_expr(), e)), e),
            H_two,
        )
        zero(
            e * g_two * sp.diff(H_two.as_expr(), e)
            - quotient_two.as_expr() * H_two.as_expr()
            - remainder_two.as_expr(),
            "lambda=2 recomputed division j=2",
        )
        gate(H_two.degree() == 17, "lambda=2 j=2 degree")
        gate(H_two.LC() == 68024448, "lambda=2 j=2 leading coefficient")
        gate(remainder_two.coeff_monomial(e**4) == 30 * mu**2,
             "lambda=2 j=2 e4 witness")
        exception_rows.append("lambda2:j2:deg17:e4:30mu^2")


# The Hamiltonian typing and finite-root reconstruction are replayed rather
# than merely cited.  Here G,D stand for g(eta),g'(eta), with e*G!=0.
g_function = sp.Function("g")(e)
A_generic = e**2 - z / 3 + r * g_function
critical = [sp.factor(bracket(A_generic, q)) for q in variables]
zero(critical[0] - (r**2 - 9 * z**2 * (2 * e + r * sp.diff(g_function, e))),
     "generic r Hamiltonian")
zero(critical[1] - (3 * g_function * r**2
     - 3 * (1 + 2 * r * e) * (2 * e + r * sp.diff(g_function, e))),
     "generic z Hamiltonian")
zero(critical[2] - (9 * g_function * z**2 - (1 + 2 * r * e)),
     "generic e Hamiltonian")
zero(bracket(A_generic, surface), "surface Casimir")

r_rec = u / e
z_rec = 9 * G * u * (1 + u) / (e * K)
zero(z_rec**2 - K / (9 * G) + Q / (9 * G * e**2 * K**2),
     "Q reconstructs z square")
zero(surface.subs({r: r_rec, z: z_rec})
     - u * (1 + u) * Q / (e**3 * K**3),
     "Q reconstructs surface")
A_e_symbol = 2 * e + r_rec * D
C_e_rec = sp.factor(9 * G * z_rec**2 - K)
C_z_rec = sp.factor(3 * G * r_rec**2 - 3 * K * A_e_symbol)
C_r_rec = sp.factor(r_rec**2 - 9 * z_rec**2 * A_e_symbol)
zero(C_e_rec + Q / (e**2 * K**2), "Q kills e Hamiltonian")
zero(C_z_rec - 3 * P / e**2, "P kills z Hamiltonian")
zero(
    C_r_rec - P / (G * e**2)
    - Q * (u**2 - P / G) / (e**4 * K**3),
    "P and Q kill r Hamiltonian",
)
zero(Q.subs(u, 0) - e**2, "exclude u=0")
zero(Q.subs(u, -1) + e**2, "exclude u=-1")
zero(Q.subs(u, -sp.Rational(1, 2)) + 729 * G**3 / 16,
     "exclude K=0")
zero(sp.LC(sp.Poly(Q, u)) + 729 * G**3,
     "Q leading coefficient excludes infinity")


semantic = {
    "carrier": "A=e^2-z/3+r*(mu*e^j+lambda*e^3); j=0,1,2; mu*lambda!=0",
    "mechanism": "boundary-only roots imply H divides e*g*Hprime",
    "generic": "lambda!=1; degH=17; j0,j1 use e6=70lambda^2; j2 uses e16 unless lambda=2",
    "lambda1": "recompute division: j0 e2=361/(81mu), j1 e6=35, j2 e4=30mu^2",
    "lambda2": "j2 recompute: e4=30mu^2",
    "recovery": "r=u/e; z=9g*u*(u+1)/(e*(1+2u))",
    "open": "general cubic with at least three monomials; degree>=4; mixed corrections",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
generic_blob = "\n".join(generic_rows).encode()
exception_blob = "\n".join(exception_rows).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3806-binomial-cubic-r-repairs-of-nodal-carriers-have-critical-points")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("carrier=A=e^2-z/3+r*(mu*e^j+lambda*e^3);j=0,1,2;mu*lambda!=0")
print("boundary_criterion=all_H_roots_in_V(e*g)_implies_H_divides_e*g*Hprime")
print("generic=lambda!=1;deg_H=17;LC_H=8503056lambda^3(lambda-1)^2")
print("generic_witness=j0,j1:e6=70lambda^2;j2:e16=-1062882mu^4lambda^3(lambda-2)/(lambda-1)^2")
print("lambda1=j0:deg11:e2=361/(81mu);j1:deg10:e6=35;j2:deg15:e4=30mu^2")
print("lambda2=j2:deg17:e4=30mu^2")
print("recovery=r=u/e;z=9g*u*(u+1)/(e*(1+2u))")
print("open=three_or_more_cubic_monomials;degree_ge_4;mixed_corrections;different_arm")
print(f"generic_sha256={hashlib.sha256(generic_blob).hexdigest()}")
print(f"exception_sha256={hashlib.sha256(exception_blob).hexdigest()}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
