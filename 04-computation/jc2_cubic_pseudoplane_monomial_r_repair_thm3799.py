#!/usr/bin/env python3
"""Exact companion for THM-3799's monomial r-repair obstruction."""

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
lam = sp.symbols("lambda", nonzero=True)
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


# Generic polynomial g(e) and the exact Hamiltonian reduction.
g = sp.Function("g")(e)
gp = sp.diff(g, e)
A = e**2 - z / 3 + r * g
K_source = 1 + 2 * r * e
critical = [sp.factor(bracket(A, q)) for q in variables]
zero(critical[0] - (r**2 - 9 * z**2 * (2 * e + r * gp)),
     "generic r Hamiltonian")
zero(critical[1] - (3 * g * r**2 - 3 * K_source * (2 * e + r * gp)),
     "generic z Hamiltonian")
zero(critical[2] - (9 * g * z**2 - K_source),
     "generic e Hamiltonian")
zero(bracket(A, surface), "surface is a Casimir")


# Replace g(e),g'(e) by independent coefficient symbols.  With u=re, the
# z- and e-Hamiltonians reduce to P and the square/cube compatibility Q.
G, D = sp.symbols("G D")
K = 1 + 2 * u
P = sp.expand(G * u**2 - K * (2 * e**3 + u * e * D))
Q = sp.expand(e**2 * K**3 - 729 * G**3 * u**2 * (1 + u)**2)
gate(sp.degree(P, u) == 2, "P generic degree two")
gate(sp.degree(Q, u) == 4, "Q generic degree four")
zero(sp.LC(sp.Poly(P, u)) - (G - 2 * e * D), "P leading coefficient")
zero(sp.LC(sp.Poly(Q, u)) + 729 * G**3, "Q leading coefficient")
zero(P.subs(u, 0) + 2 * e**3, "P boundary u=0")
zero(Q.subs(u, 0) - e**2, "Q boundary u=0")
zero(Q.subs(u, -1) + e**2, "Q boundary u=-1")
zero(Q.subs(u, -sp.Rational(1, 2)) + 729 * G**3 / 16,
     "Q boundary K=0")


# The exact universal resultant.  Dividing its forced G^3 e^4 boundary
# factor leaves H; monomial specialization factors H once more as g*J_m.
resultant_generic = sp.factor(sp.resultant(P, Q, u))
gate(sp.rem(resultant_generic, G**3 * e**4, G) == 0,
     "generic resultant forced factor")
H = sp.factor(resultant_generic / (G**3 * e**4))
gate(sp.Poly(H, G, D, e).total_degree() > 0, "generic residual resultant")


def J_formula(m: int) -> sp.Expr:
    """Closed residual resultant after g=lambda*e^m and removing g."""
    value = sp.Integer(1 - 2 * m)
    value += 23328 * (2 * m - 1) * e**7
    value -= 17496 * (2 * m - 1) * (m - 1) * lam * e**(m + 4)
    value += (
        2916 * (m - 1) * (3 * m**2 - 3 * m + 1)
        * lam**2 * e**(2 * m + 1)
    )
    value += 8503056 * lam**2 * e**(2 * m + 8)
    if m > 0:
        value -= 729 * m**2 * (m - 1)**2 * lam**3 * e**(3 * m - 2)
    value -= 8503056 * (m - 1) * lam**3 * e**(3 * m + 5)
    value += 2125764 * (m - 1)**2 * lam**4 * e**(4 * m + 2)
    return sp.expand(value)


monomial_rows: list[str] = []
for m in range(257):
    gm = lam * e**m
    specialized_H = sp.expand(H.subs({G: gm, D: sp.diff(gm, e)}))
    Jm = J_formula(m)
    zero(specialized_H - gm * Jm, f"closed J formula m={m}")
    specialized_resultant = sp.expand(
        resultant_generic.subs({G: gm, D: sp.diff(gm, e)})
    )
    zero(specialized_resultant - lam**4 * e**(4 * m + 4) * Jm,
         f"full resultant factor m={m}")
    gate(sp.expand(Jm).subs(e, 0) == 1 - 2 * m,
         f"J nonzero constant m={m}")
    gate(1 - 2 * m != 0, f"P leading scalar m={m}")
    gate(sp.degree(Jm, e) > 0, f"J nonconstant m={m}")

    if m != 3:
        gate(sp.expand(Jm).coeff(e, 7) == 23328 * (2 * m - 1),
             f"unique e7 coefficient m={m}")
        gate(23328 * (2 * m - 1) != 0,
             f"unique e7 is nonzero m={m}")
    else:
        e14 = sp.factor(sp.expand(Jm).coeff(e, 14))
        e7 = sp.factor(sp.expand(Jm).coeff(e, 7))
        zero(e14 - 8503056 * lam**2 * (lam - 1)**2,
             "m=3 top seam")
        gate(sp.expand(e7.subs(lam, 1)) == 26244,
             "m=3 lambda=1 surviving e7 coefficient")

    # At any nonzero root of J_m, both resultant degrees remain their
    # advertised degrees, since lambda and e are nonzero.
    P_lc = sp.factor((G - 2 * e * D).subs({G: gm, D: sp.diff(gm, e)}))
    Q_lc = sp.factor((-729 * G**3).subs(G, gm))
    zero(P_lc - lam * e**m * (1 - 2 * m),
         f"specialized P leading coefficient m={m}")
    zero(Q_lc + 729 * lam**3 * e**(3 * m),
         f"specialized Q leading coefficient m={m}")
    monomial_rows.append(
        f"{m}:{sp.degree(Jm,e)}:{sp.factor(sp.LC(sp.Poly(Jm,e)))}"
    )


# The exceptional collision is checked independently of the loop so that a
# future rewrite cannot silently route m=3 through the generic e^7 claim.
J3 = sp.factor(J_formula(3))
J3_at_one = sp.factor(J3.subs(lam, 1))
zero(J3_at_one - (26244 * e**7 - 5), "m=3 lambda=1 exact polynomial")
gate(J3_at_one.subs(e, 0) == -5, "m=3 lambda=1 roots nonzero")


# Exact point reconstruction from a common P,Q root satisfying the audited
# denominator conditions e*G*u*(u+1)*K != 0.
r_rec = u / e
z_rec = 9 * G * u * (1 + u) / (e * K)
z_square_error = sp.factor(z_rec**2 - K / (9 * G))
zero(z_square_error + Q / (9 * G * e**2 * K**2),
     "Q reconstructs z square")
surface_rec = sp.factor(surface.subs({r: r_rec, z: z_rec}))
zero(surface_rec - u * (1 + u) * Q / (e**3 * K**3),
     "Q reconstructs source surface")

A_e_symbol = 2 * e + r_rec * D
C_e_rec = sp.factor(9 * G * z_rec**2 - K)
zero(C_e_rec + Q / (e**2 * K**2), "Q kills e Hamiltonian")
C_z_rec = sp.factor(3 * G * r_rec**2 - 3 * K * A_e_symbol)
zero(C_z_rec - 3 * P / e**2, "P kills z Hamiltonian")

# The remaining Hamiltonian is forced by the exact Casimir identity and the
# nonzero coefficient F_r=K.  Check its explicit P,Q reduction as well.
C_r_rec = sp.factor(r_rec**2 - 9 * z_rec**2 * A_e_symbol)
casimir_rec = sp.factor(
    K * C_r_rec - 3 * z_rec**2 * C_z_rec + r_rec**2 * C_e_rec
)
zero(casimir_rec, "reconstructed Casimir relation")
zero(
    C_r_rec - P / (G * e**2)
    - Q * (u**2 - P / G) / (e**4 * K**3),
    "P and Q kill the remaining Hamiltonian",
)


semantic = {
    "carrier": "A=e^2-z/3+lambda*r*e^m; m>=0; lambda!=0; c=1",
    "critical_reduction": "P=g*u^2-(1+2u)(2e^3+u*e*g'); Q=e^2(1+2u)^3-729g^3u^2(1+u)^2",
    "resultant": "Res_u(P,Q)=lambda^4*e^(4m+4)*J_m(e)",
    "root_gate": "J_m(0)=1-2m; J_m nonconstant; fixed P/Q leading coefficients",
    "recovery": "r=u/e; z=9g*u*(u+1)/(e*(1+2u))",
    "exception": "m=3;lambda=1 gives J=26244e^7-5",
    "open": "nonmonomial polynomial g; z^2/r-mixed corrections; different arm profile",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
monomial_blob = "\n".join(monomial_rows).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3799-monomial-r-repairs-of-nodal-carriers-have-critical-points")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("carrier=A=e^2-z/3+lambda*r*e^m;m>=0;lambda!=0")
print("critical_variables=e;u=re")
print("P=g*u^2-(1+2u)*(2e^3+u*e*gprime)")
print("Q=e^2*(1+2u)^3-729*g^3*u^2*(1+u)^2")
print("resultant=lambda^4*e^(4m+4)*J_m(e)")
print("root_gate=J_m(0)=1-2m;J_m_nonconstant")
print("boundary=u_notin_0,-1,-1/2;e*g*(1+2u)!=0")
print("recovery=r=u/e;z=9g*u*(u+1)/(e*(1+2u))")
print("exception=m=3;lambda=1;J=26244e^7-5")
print("open=nonmonomial_g;z2_r_mixed;different_arm_profile;Darboux_pair")
print(f"monomial_sha256={hashlib.sha256(monomial_blob).hexdigest()}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
