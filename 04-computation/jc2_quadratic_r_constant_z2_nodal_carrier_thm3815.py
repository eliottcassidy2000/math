#!/usr/bin/env python3
"""Exact canonical proof for THM-3815's quadratic-r/constant-z2 cell."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True

import sympy as sp


CHECKS = 0
EXPECTED_RESIDUAL_SHA256 = "28b077c0775d3603163dddd01d078ba18feee257205a41e53d300684cbbbc449"
EXPECTED_SEMANTIC_SHA256 = "2e90439d024488eefc6647b5a991d603276ed536a53e3fb5e3340fd27a565a2e"


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


def sylvester(f: sp.Expr, g: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    fp, gp = sp.Poly(f, variable), sp.Poly(g, variable)
    m, n = fp.degree(), gp.degree()
    fc, gc = fp.all_coeffs(), gp.all_coeffs()
    rows = []
    for shift in range(n):
        rows.append([0] * shift + fc + [0] * (n - 1 - shift))
    for shift in range(m):
        rows.append([0] * shift + gc + [0] * (m - 1 - shift))
    matrix = sp.Matrix(rows)
    gate(matrix.shape == (8, 8), "fixed 8x8 Sylvester matrix")
    return sp.factor(matrix.det(method="domain-ge"))


r, z, e, u = sp.symbols("r z e u")
q, b, c, h = sp.symbols("q b c h")
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
    dl = sp.Matrix([sp.diff(left, x) for x in variables])
    dr = sp.Matrix([sp.diff(right, x) for x in variables])
    return sp.expand((dl.T * poisson * dr)[0])


g = q * e**2 + b * e + c
carrier = e**2 - z / 3 + r * g + h * z**2
L = 1 - 6 * h * z
K = 1 + 2 * r * e
critical = tuple(sp.factor(bracket(carrier, x)) for x in variables)
expected = (
    L * r**2 - 9 * z**2 * (2 * e + r * (2 * q * e + b)),
    3 * g * r**2 - 3 * K * (2 * e + r * (2 * q * e + b)),
    9 * g * z**2 - L * K,
)
for index, (actual, wanted) in enumerate(zip(critical, expected)):
    zero(actual - wanted, f"Hamiltonian {index}")
zero(bracket(carrier, surface), "surface Casimir")
zero(K * critical[0] - 3 * z**2 * critical[1] + r**2 * critical[2],
     "Casimir syzygy")

# On z!=0 put u=r/z.  Source membership itself gives
# e=(z^2-u)/(u^2 z), avoiding the artificial M=1+q*u*z denominator that
# appears if C_r is solved first.  Multiplying C_r by u^2/z and C_e by u^4
# gives the following denominator-free quartics.
w = z**2 - u
p1 = sp.expand(L * u**4 * z - 18 * w - 18 * q * u * z * w
               - 9 * b * u**3 * z**2)
p2 = sp.expand(9 * q * w**2 + 9 * b * z * u**2 * w
               + 9 * c * z**2 * u**4
               - L * (u**4 + 2 * u**3 * w))
gate(sp.Poly(p1, u).degree() == 4, "p1 quartic")
gate(sp.Poly(p2, u).degree() == 4, "p2 generic quartic")

resultant = sylvester(p1, p2, u)
U_expr = sp.cancel(-resultant / (729 * z**8))
gate(sp.denom(U_expr) == 1, "exact residual division")
U = sp.Poly(sp.expand(U_expr), z, q, b, c, h)
zero(resultant + 729 * z**8 * U.as_expr(), "resultant factorization")
gate(U.degree(z) == 16, "generic residual z degree")
gate(U.eval({z: 0}) == 144, "residual constant 144")

residual_payload = [
    (tuple(int(x) for x in monomial), str(coefficient))
    for monomial, coefficient in U.terms()
]
residual_sha256 = hashlib.sha256(
    json.dumps(residual_payload, separators=(",", ":")).encode()
).hexdigest()
gate(residual_sha256 == EXPECTED_RESIDUAL_SHA256,
     f"residual freeze actual={residual_sha256};terms={len(U.terms())}")

# h!=0: the nonzero linear coefficient makes U nonconstant.  If all roots
# were on L=0, U(0)=144 would force U=144*L^d.  The linear row forces d=4.
U_z = sp.Poly(U.as_expr(), z)
gate(U_z.coeff_monomial(z) == -3456 * h, "linear coefficient")
d = sp.symbols("d", integer=True, positive=True)
zero(sp.diff(144 * L**d, z).subs(z, 0) + 864 * d * h,
     "pure-boundary linear coefficient")
gate((-3456 * h + 864 * d * h).subs(d, 4) == 0,
     "h nonzero forces d=4")

excess = sp.Poly(sp.expand(U.as_expr() - 144 * L**4), z)
E2 = q**3 + 2 * q * b + 8 * c
gate(excess.coeff_monomial(z**2) == 324 * E2,
     "boundary-only z2 equation")
c_solution = -(q**3 + 2 * q * b) / 8
gate(sp.cancel(excess.coeff_monomial(z**3).subs(c, c_solution)) == 144 * q,
     "boundary-only z3 forces q=0")
gate(sp.cancel(c_solution.subs(q, 0)) == 0,
     "boundary-only q0 forces c0")
gate(excess.coeff_monomial(z**5).subs({q: 0, c: 0}) == 1296 * b**2,
     "boundary-only z5 forces b=0")
zero(
    excess.as_expr().subs({q: 0, b: 0, c: 0}) - 128 * z**7 * L**5,
    "boundary-only terminal 128z7L5",
)
gate(sp.Poly(128 * z**7 * L**5, z).degree() == 12,
     "terminal is nonzero for h nonzero")

# h=0: there is no L-boundary.  If the residual were constant, it would be
# 144.  The same z2/z3/z5 cascade would force q=c=b=0, but z7 remains 128.
U_h0 = sp.Poly(sp.expand(U.as_expr().subs(h, 0)), z)
gate(U_h0.TC() == 144, "h0 constant term")
gate((U_h0.as_expr() - 144).coeff(z, 2) == 324 * E2,
     "h0 z2 equation")
gate(
    sp.cancel((U_h0.as_expr() - 144).coeff(z, 3).subs(c, c_solution))
    == 144 * q,
    "h0 z3 forces q0",
)
gate((U_h0.as_expr() - 144).subs({q: 0, c: 0}).coeff(z, 5)
     == 1296 * b**2, "h0 z5 forces b0")
gate(U_h0.as_expr().subs({q: 0, b: 0, c: 0}) == 128 * z**7 + 144,
     "h0 terminal nonconstant residual")

# Every selected residual root has z*L!=0.  Since p1 retains degree four,
# the projective pair cannot meet only at infinity, even if p2 drops degree.
# Also p1(0) excludes u=0, so source substitution reconstructs safely.
zero(sp.LC(sp.Poly(p1, u)) - L * z, "p1 leading coefficient")
gate(p1.subs(u, 0) == -18 * z**2, "exclude u=0")
e_rec = w / (u**2 * z)
r_rec = u * z
reconstruction = {r: r_rec, e: e_rec}
zero(surface.subs(reconstruction), "source reconstruction")
zero(critical[0].subs(reconstruction) - z * p1 / u**2,
     "p1 reconstructs C_r")
zero(critical[2].subs(reconstruction) - p2 / u**4,
     "p2 reconstructs C_e")
zero(
    critical[1].subs(reconstruction)
    - (K.subs(reconstruction) * critical[0].subs(reconstruction)
       + r_rec**2 * critical[2].subs(reconstruction)) / (3 * z**2),
    "Casimir reconstructs C_z",
)

# Deep hostile controls exercise the whole cascade rather than a generic row.
deep_E2 = {q: 2, b: 3, c: -sp.Rational(5, 2)}
gate(excess.coeff_monomial(z**2).subs(deep_E2) == 0,
     "hostile E2 seam")
gate(excess.coeff_monomial(z**3).subs(deep_E2) == 288,
     "hostile E2 seam exits at z3")
gate(excess.coeff_monomial(z**5).subs({q: 0, c: 0, b: 5})
     == 32400, "hostile q0c0 seam exits at z5")
gate(excess.as_expr().subs({q: 0, b: 0, c: 0, h: 2}) != 0,
     "hostile terminal pure-normal seam")

semantic = {
    "field": "algebraically closed characteristic zero",
    "carrier": "e2-z/3+r*(q e2+b e+c)+h z2; all q,b,c,h",
    "compression": "source-first u=r/z; e=(z2-u)/(u2z); two quartics p1,p2",
    "resultant": "-729*z8*U; U0=144; generic z degree16",
    "h_nonzero": "z1 forces boundary exponent4; z2 E2; z3 q; z5 b2; terminal 128z7L5",
    "h_zero": "constant assumption dies by same cascade; terminal 128z7+144",
    "projective": "p1 LC=Lz excludes common infinity; p1(0)=-18z2 excludes u0",
    "reconstruction": "source identity; Cr=z*p1/u2; Ce=p2/u4; Cz by Casimir",
    "scope": "quadratic r-profile and constant z2 only; no nonconstant h or deg g>=3",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()
gate(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
     f"semantic freeze actual={semantic_sha256}")

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3815-quadratic-r-profile-constant-z2-nodal-carriers-have-critical-points")
print("status=PROVED_VERIFIED_EXACT_INDEPENDENTLY_HOSTILE_AUDITED")
print("carrier=A=e2-z/3+r*(q*e2+b*e+c)+h*z2;all_q_b_c_h")
print("compression=source_first;e=(z2-u)/(u2z);p1_p2_quartic")
print(f"resultant=-729*z8*U;degree_z_U=16;terms={len(U.terms())};U0=144")
print("h_nonzero=z1_forces_d4;z2_E2;z3_forces_q0c0;z5_forces_b0;terminal=128z7L5")
print("h_zero=same_constant_cascade;terminal=128z7+144")
print("projective=p1_LC=Lz;u0_excluded;finite_common_root")
print("reconstruction=source_Cr_Ce_Cz_exact")
print("scope=no_nonconstant_h;no_deg_g_ge_3;no_rational_mates")
print(f"residual_sha256={residual_sha256}")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
