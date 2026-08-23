#!/usr/bin/env python3
"""Independent assertion-free hostile audit of the THM-3822 candidate.

This companion deliberately does not import or execute the primary script.
It uses an explicit global Bezout reconstruction, an original-packet reverse
substitution, and separate elementary-word calculations.
"""

from __future__ import annotations

import ast
import hashlib
import json
from math import comb
from pathlib import Path

import sympy as sp


CHECKS = 0


def require(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)
    CHECKS += 1


def zero(expression: object, label: str) -> None:
    require(sp.cancel(sp.expand(expression)) == 0, label)


def quotient_zero(basis: sp.GroebnerBasis, expression: object, label: str) -> None:
    remainder = basis.reduce(sp.expand(expression))[1]
    zero(remainder, label)


# Intrinsic coordinates use capital letters to keep this derivation separate
# from the primary companion's source-level formulas.
C, L, H, K, T = sp.symbols("C L H K T")
determinant = C*K-H*L-1
Q = 3*K**2+7*H**2
N = 3+2*H*L
R = H*(14*K+9*H)
S = K*L+3*H*C**2
psi = sp.expand(N*R-Q*S)
eq_Q = Q*T-N
eq_R = R*T-S

# Independent compatibility/elimination check.
zero(sp.resultant(eq_Q, eq_R, T)-psi,
     "linear denominator compatibility is psi")

# Construct a global Bezout inverse explicitly.  First obtain H^3 and K^3
# from Q,R, then expand (CK-HL)^5: every degree-five term contains H^3 or K^3.
u_H = sp.Rational(196, 1615)*H
v_H = -sp.Rational(3, 1615)*(14*K-9*H)
u_K = sp.Rational(1, 4845)*(1615*K+882*H)
v_K = sp.Rational(1, 4845)*(-189*K-686*H)
zero(u_H*Q+v_H*R-H**3, "independent H-cube Bezout reduction")
zero(u_K*Q+v_K*R-K**3, "independent K-cube Bezout reduction")

U = sp.S.Zero
V = sp.S.Zero
for i in range(6):
    term = comb(5, i)*(C*K)**i*(-H*L)**(5-i)
    if i <= 2:
        multiplier = sp.cancel(term/H**3)
        U += multiplier*u_H
        V += multiplier*v_H
    else:
        multiplier = sp.cancel(term/K**3)
        U += multiplier*u_K
        V += multiplier*v_K
U = sp.expand(U)
V = sp.expand(V)
zero(U*Q+V*R-(C*K-H*L)**5,
     "global degree-five Bezout identity")

# On determinant one, D0=UN+VS is the unique simultaneous denominator.
D0 = sp.expand(U*N+V*S)
intrinsic_basis = sp.groebner(
    [determinant, psi], T, L, C, K, H, order="grevlex", domain=sp.QQ
)
quotient_zero(intrinsic_basis, Q*D0-N, "global D0 solves QD=N")
quotient_zero(intrinsic_basis, R*D0-S, "global D0 solves RD=S")

denominator_basis = sp.groebner(
    [determinant, eq_Q, eq_R], T, L, C, K, H,
    order="grevlex", domain=sp.QQ
)
quotient_zero(denominator_basis, T-D0,
              "auxiliary denominator equals global Bezout reconstruction")

# Rebuild the nonlinear Delone--Faddeev packet without using the primary
# companion's reduction order.
A_rec = H*T
W_rec = K*T
Theta_rec = (L-14*H*T)/3
packet = (
    W_rec**2+7*A_rec**2-C*W_rec+A_rec*Theta_rec,
    W_rec*Theta_rec-3*A_rec**2+A_rec*C**2,
    Theta_rec**2-3*A_rec*C+C**3-(C**2-3*A_rec)*W_rec+7*A_rec*Theta_rec,
)
for index, relation in enumerate(packet, 1):
    quotient_zero(denominator_basis, 9*relation,
                  f"reconstructed original packet relation {index}")
calculated_T = C*W_rec-3*A_rec*Theta_rec-14*A_rec**2
quotient_zero(denominator_basis, calculated_T-T,
              "reconstructed derivative denominator is T")

# Reverse direction in the original packet's fraction field.  This guards
# against merely finding a parasitic intrinsic component.
A0, C0, W0, Theta0 = sp.symbols("A0 C0 W0 Theta0")
packet_original = (
    W0**2+7*A0**2-C0*W0+A0*Theta0,
    W0*Theta0-3*A0**2+A0*C0**2,
    Theta0**2-3*A0*C0+C0**3-(C0**2-3*A0)*W0+7*A0*Theta0,
)
T0 = C0*W0-3*A0*Theta0-14*A0**2
H0 = A0/T0
K0 = W0/T0
L0 = 3*Theta0+14*A0
reverse_det_num = sp.together(
    C0*K0-H0*L0-1
).as_numer_denom()[0]
zero(reverse_det_num, "original packet gives determinant one")
reverse_psi_num = sp.together(
    psi.subs({C: C0, L: L0, H: H0, K: K0})
).as_numer_denom()[0]
original_basis = sp.groebner(
    packet_original, Theta0, W0, C0, A0, order="grevlex"
)
quotient_zero(original_basis, reverse_psi_num,
              "original packet gives intrinsic psi")

# The determinant cover D(H) union D(K) has no hidden component.  Each chart
# is a quadratic with a nonsquare squarefree discriminant core.
phi_H = sp.together(
    psi.subs(L, (C*K-1)/H)
).as_numer_denom()[0]
disc_H = sp.factor(sp.discriminant(phi_H, C))
disc_core = sp.factor(disc_H/9)
require(sp.Poly(phi_H, C).degree() == 2, "H-chart is quadratic")
require(sp.Poly(disc_core, H, K).total_degree() > 0,
        "chart discriminant is nonconstant")
require(sp.Poly(disc_core, H, K).sqf_part() == sp.Poly(disc_core, H, K),
        "chart discriminant core is squarefree")
phi_K = sp.together(
    psi.subs(C, (1+H*L)/K)
).as_numer_denom()[0]
zero(sp.discriminant(phi_K, L)-9*K**2*disc_core,
     "K-chart has the same discriminant squareclass")
overlap_equation = sp.factor(
    psi.subs({H: 1, K: 1, C: L+1})
)
zero(overlap_equation+3*(10*L**2+8*L-13),
     "two domain charts have a nonempty overlap")
require(sp.discriminant(10*L**2+8*L-13, L) != 0,
        "overlap point exists over the algebraic closure")
unit_ideal = sp.groebner(
    [determinant, Q, R], L, C, K, H, order="grevlex"
)
require(len(unit_ideal.polys) == 1 and unit_ideal.polys[0].as_expr() == 1,
        "Q,R have no simultaneous SL2 stratum")

# Compare directly with the stronger theorem now present on origin/main.  Its
# first-row square sidecar comes from a quadratic in the reconstructed T.
B0 = K**5-7*H**2*K**3-3*H**2*K**2-6*H**3*K**2-7*H**4
T_quadratic = H**2*Q**2*T**2+2*B0*T+H**2-2*K**3
quotient_zero(denominator_basis, T_quadratic,
              "origin theorem intrinsic denominator quadratic")
zero(sp.discriminant(T_quadratic, T)-4*K**2*disc_core,
     "origin theorem denominator discriminant")
pell_L = K**4-7*H**2*K**2+6*H**2*K-6*H**3*K
pell_R = 49*K**3+(27*H-9)*K**2+49*H**2*K+21*H**3
zero(disc_core-pell_L**2-4*H**4*pell_R,
     "origin theorem Pell completion")
ratio_p, ratio_r = sp.symbols("ratio_p ratio_r")
ratio_cubic = (
    49*ratio_p**3+27*ratio_p**2*ratio_r
    +49*ratio_p*ratio_r**2+21*ratio_r**3
)
require(sp.discriminant(ratio_cubic.subs(ratio_r, 1), ratio_p)
        == -27046348,
        "Pell-zero reduced-ratio cubic has distinct factors")
zero(ratio_cubic.subs(ratio_p, 0)-21*ratio_r**3,
     "Pell-zero ratio denominator is coprime to numerator")

# Punctured-arm quotient.  The Groebner basis itself derives L=0, and adding
# the denominator equations derives T=C^2 globally on H=0.
arm_basis = sp.groebner(
    [determinant, psi, H], T, L, C, K, H, order="grevlex"
)
quotient_zero(arm_basis, L, "H=0 forces L=0")
arm_denominator_basis = sp.groebner(
    [determinant, eq_Q, eq_R, H], T, L, C, K, H, order="grevlex"
)
quotient_zero(arm_denominator_basis, T-C**2, "H=0 forces T=C^2")

# The origin/main theorem's hyperbolic argument uses this monic squarefree
# degree-eight fibre, hence a genus-three smooth completion.
u, tau, Z = sp.symbols("u tau Z")
hyperbolic_fibre = sp.Poly(disc_core.subs({H: -1, K: u}), u)
require(hyperbolic_fibre.degree() == 8 and hyperbolic_fibre.LC() == 1,
        "hyperbolic square sidecar is monic degree eight")
require(sp.gcd(hyperbolic_fibre, hyperbolic_fibre.diff()).degree() == 0,
        "hyperbolic square sidecar is squarefree")
generic_hyperbolic = sp.Poly(
    disc_core.subs({H: tau-1, K: Z}), Z
)
require(generic_hyperbolic.degree() == 8
        and generic_hyperbolic.LC() == 1,
        "generic hyperbolic sidecar remains monic degree eight")
zero(generic_hyperbolic.as_expr().subs(tau, 0)
     -hyperbolic_fibre.as_expr().subs(u, Z),
     "generic hyperbolic fibre has squarefree specialization")

# Direct hostile: determinant and c^2*J(h,c)=1 mod h survive, but psi does not.
x, y = sp.symbols("x y")
Hh = x**4*y-1
Ch = x**3*y
Kh = x*(1-Hh)
Lh = -Hh
J_HC = sp.diff(Hh, x)*sp.diff(Ch, y)-sp.diff(Hh, y)*sp.diff(Ch, x)
zero(Ch*Kh-Hh*Lh-1, "hostile determinant")
zero(J_HC-x**6*y, "hostile Jacobian")
hostile_lift = sp.factor(Ch**2*J_HC-1)
zero(hostile_lift-Hh*(x**8*y**2+x**4*y+1),
     "hostile punctured-arm congruence")
hostile_psi = sp.factor(
    psi.subs({C: Ch, L: Lh, H: Hh, K: Kh})
)
require(hostile_psi != 0, "hostile genuinely fails psi")
require(sp.rem(sp.Poly(hostile_psi, y), sp.Poly(Hh, y)) == 0,
        "hostile failure restricts to zero on the arm")

# Named E_+(p)E_-(s) cell, independently normalized.
p, s = sp.symbols("p s")
forward = {C: 1+p*s, L: p, H: s, K: 1}
zero(determinant.subs(forward), "forward two-shear determinant")
forward_equation = sp.factor(-psi.subs(forward)/3)
a2 = s**3*(7*s**2+3)
b2 = 14*s**4-6*s**3-s**2+1
c2 = s*(7*s**2-9*s-11)
zero(forward_equation-(a2*p**2+b2*p+c2),
     "forward two-shear quadratic")
R7 = sp.factor(sp.discriminant(forward_equation, p))
require(sp.Poly(R7, s).degree() == 7, "forward discriminant degree seven")
require(sp.gcd(sp.Poly(R7, s), sp.Poly(sp.diff(R7, s), s)) == 1,
        "forward discriminant has seven distinct roots")
zero((2*a2*p+b2)**2-R7-4*a2*forward_equation,
     "forward completing-square identity")
require(sp.gcd_list([a2, b2, c2]) == 1,
        "constant-s coefficient packet has no common zero")

# Opposite E_-(s)E_+(p) orientation.  Its restriction is nonzero, hence every
# solution has trdeg k[p,s] <= 1; this already excludes a dominant map into
# the two-dimensional intrinsic surface, though it does not force constants.
opposite = {C: 1, L: p, H: s, K: 1+s*p}
zero(determinant.subs(opposite), "opposite two-shear determinant")
opposite_equation = sp.factor(-psi.subs(opposite)/3)
opposite_expected = (
    p**4*s**3+3*p**3*s**2-4*p**2*s**3+3*p**2*s
    -6*p*s**3-15*p*s**2+p+7*s**3-9*s**2-11*s
)
zero(opposite_equation-opposite_expected,
     "opposite two-shear nonzero relation")
require(opposite_equation != 0, "opposite relation is not tautological")
opposite_disc_s = sp.factor(sp.discriminant(opposite_equation, s))
K4 = 1615*p**4+6554*p**2-16116*p-47069
zero(opposite_disc_s+K4, "opposite cubic discriminant")
require(sp.gcd(sp.Poly(K4, p), sp.Poly(sp.diff(K4, p), p)) == 1,
        "opposite discriminant quartic is squarefree")

# Cheapest longer word E_+(p)E_-(s)E_+(t).  The incoming theorem's G_m-arm
# sidecar actually closes it: on H=s=0 its companion K=1 is constant.  If s
# is constant instead, the nonzero intrinsic relation leaves trdeg at most 1.
t = sp.symbols("t")
long_cell = {
    C: 1+p*s,
    L: p+(1+p*s)*t,
    H: s,
    K: 1+s*t,
}
zero(determinant.subs(long_cell), "length-three determinant")
long_equation = sp.factor(-psi.subs(long_cell)/3)
require(sp.Poly(long_equation, p).degree() == 2,
        "length-three equation is quadratic in p")
long_disc = sp.factor(sp.discriminant(long_equation, p))
zero(long_disc-disc_core.subs({H: s, K: 1+s*t}),
     "length-three discriminant is bivariate chart pullback")
require(sp.Poly(long_disc, s, t).sqf_part() == sp.Poly(long_disc, s, t),
        "length-three bivariate discriminant is squarefree")
zero(long_equation.subs(s, 0)-(p+t),
     "length-three cancellation boundary")
zero(long_disc.subs(s, 0)-1,
     "length-three boundary has a nonconstant square lift")
long_constant_s_coefficients = sp.Poly(long_equation, p, t).coeffs()
require(sp.gcd_list(long_constant_s_coefficients) == 1,
        "no constant s makes the length-three relation tautological")

# The opposite three-shear orientation is already arm-capable.  In the
# origin convention this is E_+(p)E_-(s)E_+(t); modulo its arm H3, K3 has
# inverse C3=1+st and need not be constant.  This is a stopping reason, not
# an atlas construction.
C3_opposite = 1+s*t
L3_opposite = s
H3_opposite = p+t+p*s*t
K3_opposite = 1+p*s
zero(K3_opposite*C3_opposite-1-s*H3_opposite,
     "opposite length-three determinant and nonconstant arm unit")
opposite_long_equation = sp.factor(-psi.subs({
    C: C3_opposite,
    L: L3_opposite,
    H: H3_opposite,
    K: K3_opposite,
})/3)
require(sp.Poly(opposite_long_equation, t).degree() == 5,
        "opposite length-three target relation is nontrivial but dimensional")

source = Path(__file__).read_text(encoding="utf-8")
require(not any(isinstance(node, ast.Assert)
                for node in ast.walk(ast.parse(source))),
        "independent script contains no Python asserts")

semantic = {
    "verdict": "PASS third audit of the stronger origin/main THM3822 theorem",
    "reconstruction": "explicit fifth-power Bezout gives global D0; both original packet directions and both determinant charts pass",
    "origin_sidecars": "denominator quadratic, square H, Pell completion, Gm arm, and genus-three hyperbolic fibre pass independently",
    "two_shear": "both orientations are consistency corollaries of the current arm/square theorem, not new results",
    "length_three_extension": "in the transposed convention Eplus-Eminus-Eplus has H=s and K=1+st; nonconstant s violates the nonconstant-unit arm gate, while constant s leaves a nonzero relation in p,t and cannot dominate",
    "opposite_length_three_stop": "in the origin convention Eplus-Eminus-Eplus has K(1+st)=1 mod H, so its arm can carry a nonconstant unit and the new proof does not close it",
    "scope": "opposite three-shear orientation, arbitrary polynomial SL2 atlases, and longer interacting words remain open",
}
semantic_hash = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("object=THM-3822-independent-hostile-audit")
print("verdict=PASS")
print("intrinsic=GLOBAL_BEZOUT_RECONSTRUCTION_AND_TWO_CHART_DOMAIN_PASS")
print("origin_main=DENOMINATOR_SQUARE;PELL;GM_ARM;GENUS3_HYPERBOLIC_GATE_PASS")
print("arm=NECESSARY_CONGRUENCE_PASS;HOSTILE_CONFIRMS_PSI_LOAD_BEARING")
print(f"forward_discriminant={R7}")
print("two_shear_status=CURRENT_THEOREM_COROLLARY_CHECKS_ONLY")
print(f"opposite_relation={opposite_equation}")
print("length_three_extension=NO_DOMINANT_ATLAS;ARM_CONSTANT_IF_s_NONCONSTANT;TRDEG_AT_MOST_ONE_IF_s_CONSTANT")
print("length_three_hostile=s=0,p=-t_gives_identity_with_nonconstant_word_parameters")
print("opposite_length_three_stop=K3*(1+s*t)=1_mod_H3;NONCONSTANT_ARM_UNIT_POSSIBLE")
print("scope=OPPOSITE_LENGTH_THREE;ARBITRARY_SL2;LONGER_ATLASES_REMAIN_OPEN")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={semantic_hash}")
