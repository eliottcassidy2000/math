#!/usr/bin/env python3
"""Exact verifier for THM-4385's elliptic sign quotient.

This is a light geometric companion to THM-4380.  It does not rebuild the
source-normal bracket/depth compiler: THM-4380 is the proved input identifying
the two row-ten strata and the row-eleven/row-twelve eliminants.  Here we
independently verify the cubic geometry, quotient branch divisor, and the
seven-to-fourteen lift.
"""

from __future__ import annotations

import sympy as sp


CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise RuntimeError(f"FAILED: {label}")
    CHECKS += 1


P, eta, T, ratio = sp.symbols("P eta T ratio")

a = 7231154026500
b = 50541940696500
c = 6793915500000
d = 353642000625
m = 631918028977864704
n = 91584545734393856

D = sp.Poly(
    a * P**3 + b * P**2 * eta + c * P * eta**2 + d * eta**3
    - m * P - n * eta,
    P,
    eta,
    domain=sp.QQ,
)
Dbar = sp.Poly(
    a * P**3 + b * P**2 * eta + c * P * eta**2 + d * eta**3
    - m * P * T**2 - n * eta * T**2,
    P,
    eta,
    T,
    domain=sp.QQ,
)
F0 = sp.Poly(18612736875 * eta**2 - 4820239249178624, eta, domain=sp.QQ)

check(sp.expand(Dbar.as_expr().subs(T, 1) - D.as_expr()) == 0, "affine chart")
check(
    sp.expand(Dbar.as_expr().subs({P: -P, eta: -eta}) + Dbar.as_expr()) == 0,
    "simultaneous sign preserves cubic zero locus",
)
check(
    sp.expand(D.as_expr().subs(P, 0) - 19 * eta * F0.as_expr()) == 0,
    "P-zero boundary factorization",
)
check(sp.gcd(F0, F0.diff()) == sp.Poly(1, eta, domain=sp.QQ), "boundary quadratic reduced")
check(F0.eval(0) != 0, "omitted origin distinct from boundary pair")

# The affine singular ideal is the unit ideal over Q.  At infinity the binary
# cubic is governed by A; its squarefreeness excludes geometric singularities.
affine_singular = sp.groebner(
    [D.as_expr(), sp.diff(D.as_expr(), P), sp.diff(D.as_expr(), eta)],
    P,
    eta,
    order="lex",
    domain=sp.QQ,
)
check(len(affine_singular.polys) == 1 and affine_singular.polys[0].as_expr() == 1, "affine smoothness over Q")

A = sp.Poly(d * ratio**3 + c * ratio**2 + b * ratio + a, ratio, domain=sp.QQ)
B = sp.Poly(n * ratio + m, ratio, domain=sp.QQ)
H = sp.Poly(Dbar.as_expr().subs(T, 0).subs({P: 1, eta: ratio}), ratio, domain=sp.QQ)
check(H == A, "infinity cubic equals A")
check(sp.gcd(A, A.diff()) == sp.Poly(1, ratio, domain=sp.QQ), "three distinct geometric points at infinity")
check(Dbar.eval({P: 0, eta: 1, T: 0}) != 0, "infinity covered by P-nonzero chart")

# A second good-reduction route guards the exact characteristic-zero
# smoothness computation.  Unit ideal on T=1 plus squarefree infinity is a
# geometric, not merely F_7-rational, certificate.
affine_singular_7 = sp.groebner(
    [D.as_expr(), sp.diff(D.as_expr(), P), sp.diff(D.as_expr(), eta)],
    P,
    eta,
    order="lex",
    modulus=7,
)
A7 = sp.Poly(A.as_expr(), ratio, modulus=7)
check(len(affine_singular_7.polys) == 1 and affine_singular_7.polys[0].as_expr() == 1, "affine smoothness modulo 7")
check(A7.degree() == 3 and sp.gcd(A7, A7.diff()) == sp.Poly(1, ratio, modulus=7), "smooth infinity modulo 7")

check(Dbar.eval({P: 0, eta: 0, T: 1}) == 0, "rational origin lies on cubic")
check(
    Dbar.diff(P).eval({P: 0, eta: 0, T: 1}) == -m
    and Dbar.diff(eta).eval({P: 0, eta: 0, T: 1}) == -n,
    "origin tangent and quotient value",
)

ratio_relation = sp.expand(
    D.as_expr().subs(eta, ratio * P)
    - P * (P**2 * A.as_expr() - B.as_expr())
)
check(ratio_relation == 0, "quadratic quotient equation")
AB = sp.Poly(A.as_expr() * B.as_expr(), ratio, domain=sp.QQ)
check(A.degree() == 3 and B.degree() == 1 and sp.gcd(A, B) == sp.Poly(1, ratio, domain=sp.QQ), "four branch values")
check(sp.gcd(AB, AB.diff()) == sp.Poly(1, ratio, domain=sp.QQ), "branch quartic squarefree")
AB7 = sp.Poly(AB.as_expr(), ratio, modulus=7)
check(AB7.degree() == 4 and sp.gcd(AB7, AB7.diff()) == sp.Poly(1, ratio, modulus=7), "branch quartic good reduction")

K = sp.Poly(
    21252176198679866250754006276839556755825 * ratio**7
    + 799311827675117522149997435401077574131600 * ratio**6
    + 14384863896403857958176347858723924433398460 * ratio**5
    + 142730433788981669223548142320603830110956220 * ratio**4
    + 786944231209420107856657052701244375708027892 * ratio**3
    + 1852564642916723756803328543705267257790149632 * ratio**2
    - 147733098192443646925107791876239203619548432 * ratio
    + 3714269896529642422852685702214695613036368,
    ratio,
    domain=sp.QQ,
)
J = sp.Poly(
    -4607940726340893 * ratio**2
    - 2340230825466590 * ratio
    + 113720641620096958,
    ratio,
    domain=sp.QQ,
)
check(K.degree() == 7 and K.LC() != 0 and K.TC() != 0, "seven finite ratio roots")
check(sp.gcd(K, K.diff()) == sp.Poly(1, ratio, domain=sp.QQ), "ratio roots reduced")
check(sp.gcd(K, AB) == sp.Poly(1, ratio, domain=sp.QQ), "row-eleven divisor avoids ramification")
check(sp.gcd(K, J) == sp.Poly(1, ratio, domain=sp.QQ), "row-twelve divisor disjoint")

print("THM-4385 EXACT: SOURCE-NORMAL ROW-TEN ELLIPTIC SIGN QUOTIENT PASS")
print(f"CHECKS={CHECKS}")
print("CUBIC projective_smooth_over_Q=yes genus=1 rational_origin=(0,0,1)")
print("BOUNDARY D(0,eta)=19*eta*F0; source_projection=affine_cubic_minus_origin=elliptic_curve_minus_E[2]")
print("INVOLUTION sigma(P,eta,T)=(-P,-eta,T); fixed_divisor=origin_plus_three_points_at_infinity=E[2]")
print("QUOTIENT ratio=eta/P; P^2*A(ratio)=B(ratio); fixed_field=Q(ratio); branch_polynomial=A*B squarefree_degree=4")
print("ROW11 K_degree=7 squarefree=yes gcd(K,A*B)=1 etale_lifts=14")
print("ROW12 J_degree=2 gcd(K,J)=1 surviving_lifts=0")
print("MOD7 affine_singular_ideal=unit infinity_gcd=1 branch_gcd=1")
print("SCOPE geometric corollary of THM-4380 only; no chart/seam entry, all-weight lift, Keller pair, JC2, or DC2")
print("RESULT=PASS")
