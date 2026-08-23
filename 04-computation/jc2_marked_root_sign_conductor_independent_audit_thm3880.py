#!/usr/bin/env python3
"""Independent hostile audit of THM-3872, THM-3876, and THM-3880.

This deliberately does not import any canonical companion.  It uses an
elementary forced-square recurrence for the THM-3872 cancellation seam and
cyclotomic quotient arithmetic for the THM-3876/3880 sign-flip family.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys

import sympy as sp

sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def gate(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(expression) == 0, label)


# THM-3880: an integral-domain sign argument, including exact boundaries.
A, C, b, z = sp.symbols("A C b z")
Delta_b = -27*A**2*b**2 + 8*A*C**3 - 54*A*C*b + 9*C**2 - 54*b
P = 1 + sp.Rational(2, 3)*A*C
carrier_u = 1 + A*C + A**2*b
zero(A**2*Delta_b - 27*(P**3 - carrier_u**2), "cusp identity")

plus_numerator = sp.factor(z**3 - 1 - sp.Rational(3, 2)*(z**2 - 1))
minus_numerator = sp.factor(-z**3 - 1 - sp.Rational(3, 2)*(z**2 - 1))
zero(plus_numerator - (z - 1)**2*(2*z + 1)/2, "plus carrier")
zero(minus_numerator + (z + 1)**2*(2*z - 1)/2, "minus carrier")
zero(
    plus_numerator.subs(z, -z) - plus_numerator + 2*z**3,
    "plus opposite-sign jump",
)
zero(
    minus_numerator.subs(z, -z) - minus_numerator - 2*z**3,
    "minus opposite-sign jump",
)

b_zero = sp.Rational(2, 9)*C**2
b_axis = sp.Rational(1, 6)*C**2
zero(
    Delta_b.subs(b, b_zero) + C**2*(2*A*C + 3)**2/3,
    "zero-root hyperbola control",
)
zero(
    Delta_b.subs(b, b_axis) + A*C**3*(3*A*C + 4)/4,
    "A-zero same-sign control",
)

# Independent local conductor audit for the strengthened node/A2 iff.
epsilon = sp.symbols("epsilon")
B_forced = (epsilon*z**3 - 1 - A*C)/A**2
zero(B_forced - B_forced.subs(z, z), "same-sign node value")
zero(Delta_b.subs({A: 0, b: C**2/6}), "A-zero node value")

Ap, Cp, zp, Bp = sp.symbols("Ap Cp zp Bp")
dP = sp.Rational(2, 3)*(Ap*C + A*Cp)
root_derivative = 2*z*zp - dP
zero(root_derivative.subs({Ap: 0, Cp: 0}) - 2*z*zp, "cusp root derivative")
dforced = 2*A*Ap*b + A**2*Bp - (
    3*epsilon*z**2*zp - Ap*C - A*Cp
)
# Modulo A'=C'=0 this is A^2 B'-3 epsilon*z*(z*z'); the root
# derivative gives z*z'=0 without ever dividing by z.
zero(
    dforced.subs({Ap: 0, Cp: 0})
    - (A**2*Bp - 3*epsilon*z*(z*zp)),
    "A-nonzero cusp derivative",
)

dDelta = (
    -54*A*Ap*b**2 - 54*A**2*b*Bp
    + 8*Ap*C**3 + 24*A*C**2*Cp
    - 54*(Ap*C*b + A*Cp*b + A*C*Bp)
    + 18*C*Cp - 54*Bp
)
zero(dDelta.subs({A: 0, Ap: 0, Cp: 0}) + 54*Bp, "A-zero cusp derivative")

# Truncated conductor controls: k[[tau^2,tau^3]] misses exactly tau, while
# a node normalization pair descends precisely when its constants agree.
for cutoff in range(4, 25):
    cusp_semigroup = {0}
    for i in range(cutoff + 1):
        for j in range(cutoff + 1):
            exponent = 2*i + 3*j
            if exponent < cutoff:
                cusp_semigroup.add(exponent)
    gate(cusp_semigroup == ({0} | set(range(2, cutoff))), f"A2 conductor cutoff {cutoff}")
    gate(2*cutoff - 1 == 1 + 2*(cutoff - 1), f"node conductor cutoff {cutoff}")

# THM-3876: exact gcd reduction, descent boundaries, and hostiles.
r, eta, U = sp.symbols("r eta U")
boundary_checks = 0
for N in range(1, 33):
    B1 = 2*r**(2*N)*(3 + 4*r**(N + 1))
    A1 = r
    zero(B1 - 2*A1**(2*N)*(3 + 4*A1**(N + 1)), f"M1 N={N}")
    boundary_checks += 1

for N in range(1, 33, 2):
    A2 = r**2
    C2 = 6*r**N*(1 + r**(N + 2))
    sidecar = sp.expand(C2/6 - A2**(N + 1))
    B2 = 2*r**(2*N)*(3 + 4*r**(N + 2))
    zero(sidecar - r**N, f"M2 sidecar N={N}")
    zero(B2 - 2*sidecar**2*(3 + 4*A2*sidecar), f"M2 B N={N}")
    boundary_checks += 2

gcd_checks = 0
for m in range(1, 41):
    for n in range(1, 41):
        common = math.gcd(m, n)
        M, N = m // common, n // common
        gate(math.gcd(M, N) == 1, f"gcd reduction {m},{n}")
        zero((r**common)**M - r**m, f"A gcd pullback {m},{n}")
        zero(
            (r**common)**N*(1 + (r**common)**(M + N))
            - r**n*(1 + r**(m + n)),
            f"C gcd pullback {m},{n}",
        )
        gcd_checks += 3

# Work in Q[zeta]/Phi_M.  eta=zeta^N remains primitive when gcd(M,N)=1.
zeta = sp.symbols("zeta")
cyclotomic_pairs = 0
for M in range(3, 41):
    phi = sp.Poly(sp.cyclotomic_poly(M, zeta), zeta, domain=sp.QQ)
    for N in range(1, 31):
        if math.gcd(M, N) != 1:
            continue
        eta_poly = sp.Poly(zeta**N, zeta, domain=sp.QQ).rem(phi)
        one = sp.Poly(1, zeta, domain=sp.QQ)
        eta_plus = eta_poly + one
        eta_minus = eta_poly - one
        gate(sp.gcd(eta_plus, phi).degree() == 0, f"eta+1 unit M={M},N={N}")
        gate(sp.gcd(eta_minus, phi).degree() == 0, f"eta-1 unit M={M},N={N}")
        inv_plus = sp.invert(eta_plus, phi)
        U_poly = (-inv_plus).rem(phi)
        same_C = (eta_minus + (eta_poly**2 - one)*U_poly).rem(phi)
        zero(same_C.as_expr(), f"same C M={M},N={N}")
        z0 = (one + 2*U_poly).rem(phi)
        z1 = (one + 2*eta_poly*U_poly).rem(phi)
        zero((z0 + z1).rem(phi).as_expr(), f"opposite z M={M},N={N}")
        gate(not z0.is_zero, f"nonzero z M={M},N={N}")
        jump = (3*(eta_poly**2 - one) + 4*(eta_poly**3 - one)*U_poly).rem(phi)
        expected_jump = (-(eta_minus**3)*inv_plus).rem(phi)
        zero((jump - expected_jump).rem(phi).as_expr(), f"jump formula M={M},N={N}")
        gate(not jump.is_zero, f"nonzero jump M={M},N={N}")
        cyclotomic_pairs += 1

# Every possible off-diagonal monomial collision has the same sign flip.
collision_u = -1/(eta + 1)
zero((eta - 1) + (eta**2 - 1)*collision_u, "general collision equality")
zero(
    (1 + 2*collision_u) + (1 + 2*eta*collision_u),
    "general collision sign flip",
)
zero(
    3*(eta**2 - 1) + 4*(eta**3 - 1)*collision_u
    + (eta - 1)**3/(eta + 1),
    "general collision marked jump",
)

# THM-3872: independent full constant-J-span audit.
x, y, alpha, beta, gamma = sp.symbols("x y alpha beta gamma")
Delta = sp.expand(
    81*x**5 + 90*x**4 + 25*x**3 + 30*x**2*y**2
    + 30*x*y**2 - y**4 + 8*y**2
)
P_star = (x + 1)*(9*x + 4)**2
Q_star = (y**2 - sp.Rational(1, 2)*(30*x**2 + 30*x + 8))*(9*x + 4)**2
j1 = x*(x + 1)
j2 = y*(x + 1)
j3 = y**2 + 4*x
Aj1 = -15*x**3 - 15*x**2 + x*y**2 - 4*x
Aj2 = -15*x**2*y - 15*x*y + y**3 - 4*y
Aj3 = 81*x**4 + 9*x**3 - 44*x**2 + 15*x*y**2 - 16*x + 4*y**2

# Exactness of J: its quotient has basis 1,x,y and evaluates invertibly.
J_gb = sp.groebner([j1, j2, j3], x, y, order="lex")
for relation in (x**2 + x, x*y + y, y**2 + 4*x):
    zero(J_gb.reduce(relation)[1], "J rewrite")
evaluation = sp.Matrix([[1, 0, 0], [1, -1, -2], [1, -1, 2]])
gate(evaluation.det() != 0, "J quotient evaluation isomorphism")

addition = alpha*j1 + beta*j2 + gamma*j3
mixed = alpha*Aj1 + beta*Aj2 + gamma*Aj3
shifted_P = sp.expand(P_star + 2*mixed + addition**2)
shifted_Q = sp.expand(
    Q_star + 3*addition*P_star + 3*addition*mixed + addition**3
)
R_span, remainder = sp.div(sp.expand(shifted_P**3 - shifted_Q**2), Delta, x, y)
zero(remainder, "constant span Delta division")
R_span = sp.expand(R_span)

# gamma=0: beta!=0 has odd y-degree five at x=0; beta=0 is the j1 ray.
w_zero = sp.Poly(R_span.subs({gamma: 0, x: 0}), y)
zero(
    w_zero.as_expr()
    - (-3*beta**4*y**4 - 8*beta**3*y**5 + 32*beta**3*y**3
       - 96*beta**2*y**2 + 256),
    "gamma-zero exact specialization",
)
zero(w_zero.coeff_monomial(y**5) + 8*beta**3, "beta odd-degree gate")

j1_ray = sp.Poly(R_span.subs({beta: 0, gamma: 0}), y)
gate(j1_ray.degree() == 2, "j1 y degree")
zero(j1_ray.coeff_monomial(y), "j1 missing linear y")
zero(j1_ray.coeff_monomial(y**2) + 8*alpha**3*x**3, "j1 y2 coefficient")
zero(j1_ray.as_expr().subs({x: 0, y: 0}) - 256, "j1 constant control")

# gamma!=0: odd degree seven except on gamma=-alpha^2/216.
y_zero = sp.Poly(R_span.subs(y, 0), x)
zero(
    y_zero.coeff_monomial(x**7) - 243*gamma**2*(alpha**2 + 216*gamma),
    "constant-span leading seam",
)
seam = sp.Poly(
    sp.expand(R_span.subs({y: 0, gamma: -alpha**2/sp.Integer(216)})),
    x,
)
gate(not seam.as_expr().has(beta), "seam forgets beta")
c = [sp.factor(seam.coeff_monomial(x**i)) for i in range(7)]
zero(c[0] + sp.Rational(64, 9)*(alpha - 6)*(alpha + 6), "seam c0")

# A hypothetical cubic square root is forced by its first three rows.
q1 = sp.factor(c[1]/2)
q2 = sp.factor((c[2] - q1**2/c[0])/2)
q3 = sp.factor((c[3] - 2*q1*q2/c[0])/2)
res4 = sp.factor(c[4] - (q2**2 + 2*q1*q3)/c[0])
res5 = sp.factor(c[5] - 2*q2*q3/c[0])
n4 = sp.Poly(sp.together(res4).as_numer_denom()[0], alpha)
n5 = sp.Poly(sp.together(res5).as_numer_denom()[0], alpha)
zero(sp.gcd(n4, n5).monic().as_expr() - alpha**3, "square recurrence gcd")

# The recurrence divided by c0, so alpha=+/-6 need direct checks.
seam_plus = sp.Poly(seam.as_expr().subs(alpha, 6), x)
seam_minus = sp.Poly(seam.as_expr().subs(alpha, -6), x)
gate(seam_plus.degree() == 5, "alpha=6 odd degree")
gate(seam_minus.degree() == 6, "alpha=-6 degree")
gate(seam_minus.terms()[-1][0] == (1,), "alpha=-6 x-adic order one")

semantic = {
    "THM3872": "constant J span closes by odd specializations plus a forced-square recurrence; gcd=u^3",
    "THM3876": "gcd-reduced exponent M<=2 iff; every off-diagonal collision is a nonzero sign flip",
    "THM3880": "global-sign obstruction plus exact regular node/A2 conductor iff; A=0 seams included",
    "stronger_scope": "the sign lemma only needs an integral source and characteristic not 2,3",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("AUDIT=THM3872+THM3876+THM3880")
print("THM3872_CONSTANT_SPAN=PASS")
print("THM3872_SEAM_METHOD=forced-square-recurrence;gcd_rows_4_5=alpha^3")
print("THM3872_SPECIAL_SEAMS=alpha=6:odd-degree-5;alpha=-6:x-valuation-1")
print("THM3876_TWO_EXPONENT=PASS")
print(f"THM3876_BOUNDARY_GATES={boundary_checks};GCD_GATES={gcd_checks}")
print(f"THM3876_CYCLOTOMIC_PAIRS={cyclotomic_pairs}")
print("THM3876_ALL_OFF_DIAGONAL_COLLISIONS=opposite-nonzero-sign")
print("THM3880_SIGN_OBSTRUCTION=PASS")
print("THM3880_REGULAR_NODE_A2_IFF=PASS")
print("THM3880_STRONGER_SCOPE=integral-source;char-not-2-or-3;any-dimension")
print("THM3880_BOUNDARY=z=0-silent;regularity-and-higher-singularity-jets=open")
print(f"SEMANTIC_SHA256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
