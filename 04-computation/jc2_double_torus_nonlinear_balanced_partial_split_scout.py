#!/usr/bin/env python3
"""Independent affine-linear rederivation plus nonlinear split scout.

Reproduction:
  python3 04-computation/jc2_double_torus_nonlinear_balanced_partial_split_scout.py
  python3 -O 04-computation/jc2_double_torus_nonlinear_balanced_partial_split_scout.py
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


P, Q, a = sp.symbols("P Q a", nonzero=True)
h, v, g = sp.symbols("h v g")
omega = (-1 + sp.sqrt(-3))/2
roots = (sp.Integer(1), omega, sp.expand(omega**2))
lines = tuple(Q - root*P for root in roots)


# The independent affine-coordinate case has a squarefree three-line product.
gate(sp.simplify(1 + omega + omega**2) == 0,
     "primitive cube root relation")
gate(sp.expand(sp.prod(lines) - (Q**3 - P**3)) == 0,
     "difference of cubes is the three-line product")
for i in range(3):
    for j in range(i + 1, 3):
        gate(sp.simplify(roots[i] - roots[j]) != 0,
             f"cube-root lines {i},{j} are distinct")


# UFD gives exactly eight allocations of the three prime factors.
branches = {}
for mask in range(8):
    left = sp.Integer(1)
    right = sp.Integer(1)
    for index, line in enumerate(lines):
        if mask & (1 << index):
            left *= line
        else:
            right *= line
    D = sp.expand(a*left)
    S = sp.expand(right/a)
    q0 = sp.expand(S - D)
    q1 = sp.expand(S + D)
    H0 = sp.expand(q0**2 - 4*P**3)
    H1 = sp.expand(q1**2 - 4*Q**3)
    gate(sp.simplify(D*S - (Q**3 - P**3)) == 0,
         f"UFD allocation product mask={mask}")
    gate(sp.simplify(H0 - H1) == 0,
         f"common branch identity mask={mask}")
    branches[mask] = H0

gate(sorted(min(mask.bit_count(), 3-mask.bit_count()) for mask in range(8))
     == [0, 0, 1, 1, 1, 1, 1, 1],
     "complement symmetry leaves extreme and singleton types")
for mask in range(8):
    gate(sp.simplify(branches[mask] - branches[7-mask].subs(a, 1/a)) == 0,
         f"complement plus a->1/a changes q0 sign only mask={mask}")


# Singleton 1+2 split.  Cube-root rescaling makes Q-P representative.
D1 = a*(Q-P)
S1 = (Q**2 + P*Q + P**2)/a
q01 = sp.expand(S1-D1)
q11 = sp.expand(S1+D1)
H_singleton = sp.expand(q01**2 - 4*P**3)
gate(sp.expand(D1*S1 - (Q**3-P**3)) == 0,
     "singleton representative split")
gate(sp.expand(H_singleton - (q11**2-4*Q**3)) == 0,
     "singleton common branch")

quadratic_Q = sp.expand(
    Q**2 + (h**2-a**2)*Q + h**2*(h-a)**2
)
disc_Q = sp.factor(sp.discriminant(quadratic_Q, Q))
gate(sp.expand(disc_Q + (h-a)**3*(3*h+a)) == 0,
     "singleton quadratic discriminant")

conic = sp.expand(v**2 + (h-a)*(3*h+a))
Q_map = sp.expand(((h-a)*v-h**2+a**2)/2)
q01_map = sp.expand(q01.subs({P: h**2, Q: Q_map}))
q01_reduced = sp.rem(sp.Poly(q01_map-2*h**3, v), sp.Poly(conic, v)).as_expr()
H_singleton_map = sp.expand(H_singleton.subs({P: h**2, Q: Q_map}))
H_singleton_reduced = sp.rem(
    sp.Poly(H_singleton_map, v), sp.Poly(conic, v)
).as_expr()
gate(sp.simplify(q01_reduced) == 0,
     "conic map gives q0=2h^3")
gate(sp.simplify(H_singleton_reduced) == 0,
     "conic maps to the singleton branch")
gate(sp.cancel((2*h**3)/(2*h**2) - h) == 0,
     "singleton branch recovers h generically")
gate(sp.simplify(a-(-a/3)) != 0,
     "conic finite branch roots are distinct")
z = sp.symbols("z")
gate(sp.discriminant(z**2+3, z) != 0,
     "conic has two distinct infinity directions")


# Extreme 0+3 split.  Its normalization is a smooth Fermat cubic.
D0 = a
S0 = (Q**3-P**3)/a
q00 = sp.expand(S0-D0)
q10 = sp.expand(S0+D0)
H_extreme = sp.expand(q00**2-4*P**3)
gate(sp.expand(D0*S0-(Q**3-P**3)) == 0,
     "extreme representative split")
gate(sp.expand(H_extreme-(q10**2-4*Q**3)) == 0,
     "extreme common branch")
fermat = g**3-h**3-a
q00_map = sp.expand(q00.subs({P: h**2, Q: g**2}))
q00_reduced = sp.rem(sp.Poly(q00_map-2*h**3, g), sp.Poly(fermat, g)).as_expr()
H_extreme_map = sp.expand(H_extreme.subs({P: h**2, Q: g**2}))
H_extreme_reduced = sp.rem(
    sp.Poly(H_extreme_map, g), sp.Poly(fermat, g)
).as_expr()
gate(sp.simplify(q00_reduced) == 0,
     "Fermat map gives q0=2h^3")
gate(sp.simplify(H_extreme_reduced) == 0,
     "Fermat cubic maps to the extreme branch")
intermediate_extreme = sp.expand(((h**3+a)**2-Q**3).subs(Q, g**2))
gate(sp.rem(sp.Poly(intermediate_extreme, g), sp.Poly(fermat, g)).as_expr() == 0,
     "extreme intermediate equation Q^3=(h^3+a)^2")
gate(sp.gcd(z**3-1, sp.diff(z**3-1, z)) == 1,
     "Fermat cubic has three distinct infinity directions")
gate(sp.solve(
    [sp.diff(g**3-h**3-a*z**3, variable) for variable in (g, h, z)],
    (g, h, z), dict=True
) == [{g: 0, h: 0, z: 0}],
     "projective Fermat cubic gradient vanishes only at the forbidden origin")


# First nonlinear balanced partial split: exact hostile, not a universal
# nonlinear classification.  It already fails by support or branch count.
X, t = sp.symbols("X t")
f = t*(t-1)
p0 = X
p1 = X+f
F0 = f
F1 = f+(1-omega)*X
F2 = f+(1-omega**2)*X
Dn = sp.expand(a*t*F1)
Sn = sp.expand((t-1)*F2/a)
q0n = sp.expand(Sn-Dn)
q1n = sp.expand(Sn+Dn)
Hn = sp.expand(q0n**2-4*p0**3)
gate(sp.simplify(F0*F1*F2-(p1**3-p0**3)) == 0,
     "nonlinear prototype cube-factor identity")
gate(sp.simplify(Dn*Sn-(p1**3-p0**3)) == 0,
     "nonlinear prototype balanced partial split")
gate(sp.simplify(Hn-(q1n**2-4*p1**3)) == 0,
     "nonlinear prototype common branch")


def homogeneous_part(poly, degree):
    expanded = sp.Poly(sp.expand(poly), X, t)
    return sp.expand(sum(
        coefficient*X**monomial[0]*t**monomial[1]
        for monomial, coefficient in expanded.terms()
        if sum(monomial) == degree
    ))


for sign in (1, -1):
    H_sign = sp.expand(Hn.subs(a, sign))
    degree = sp.Poly(H_sign, X, t).total_degree()
    lead = sp.factor(homogeneous_part(H_sign, degree))
    gate(degree == 4, f"nonlinear prototype a={sign} has degree four")
    gate(sp.expand(lead-t**2*(t-sp.sqrt(-3)*X)**2) == 0,
         f"nonlinear prototype a={sign} has two distinct infinity directions")
    gate(sp.factor(lead/t**2).as_poly(X, t).total_degree() == 2,
         f"nonlinear prototype a={sign} has two-line infinity form")
    gate(len(sp.factor_list(lead)[1]) >= 2,
         f"nonlinear prototype a={sign} has at least two infinity supports")

gate(sp.expand(Hn.subs(a, 1)-Hn.subs(a, -1)) == 0,
     "the two support-degenerate parameter signs give the same branch")

# Absolute, rather than merely ground-field, irreducibility at a^2=1.
# On the positive sheet put X=h^2 and h=q0/(2X).  The resulting quadratic
# is nonsquare over k(h), and its squarefree normalization is a smooth conic.
delta = sp.expand(omega-omega**2)
b_special = sp.expand(1-omega**2)
G_special = sp.expand(
    t**2-(1+delta*h**2)*t+b_special*h**2+2*h**3
)
q0_special = sp.expand(q0n.subs({a: 1, X: h**2}))
H_special = sp.expand(Hn.subs(a, 1))
gate(sp.expand(q0_special-2*h**3+G_special) == 0,
     "a=1 branch equation in the recovered square root h")
gate(sp.expand(H_special.subs(X, h**2)
               - G_special*G_special.subs(h, -h)) == 0,
     "quartic pullback splits into the two square-root signs")
disc_special = sp.factor(sp.discriminant(G_special, t))
gate(sp.expand(disc_special+(h+1)**3*(3*h-1)) == 0,
     "a=1 quadratic discriminant")
disc_factors = sp.Poly(disc_special, h, extension=sp.sqrt(-3)).factor_list()[1]
gate(sorted(exponent for _, exponent in disc_factors) == [1, 3],
     "odd valuations prove the discriminant nonsquare over k(h)")
gate(sp.resultant(h+1, 3*h-1, h) != 0,
     "the odd discriminant supports are distinct")
gate(H_special.subs(X, 0) != 0,
     "localizing at X loses no polynomial component")

vn = sp.symbols("vn")
conic_special = sp.expand(vn**2+(h+1)*(3*h-1))
t_special_map = sp.expand(((h+1)*vn+1+delta*h**2)/2)
G_special_reduced = sp.rem(
    sp.Poly(sp.expand(G_special.subs(t, t_special_map)), vn),
    sp.Poly(conic_special, vn),
).as_expr()
q0_special_map = sp.expand(q0n.subs({a: 1, X: h**2, t: t_special_map}))
q0_special_reduced = sp.rem(
    sp.Poly(q0_special_map-2*h**3, vn), sp.Poly(conic_special, vn)
).as_expr()
H_special_map = sp.expand(H_special.subs({X: h**2, t: t_special_map}))
H_special_reduced = sp.rem(
    sp.Poly(H_special_map, vn), sp.Poly(conic_special, vn)
).as_expr()
gate(sp.simplify(G_special_reduced) == 0,
     "special conic maps to the recovered quadratic branch")
gate(sp.simplify(q0_special_reduced) == 0,
     "special conic recovers h=q0/(2X) generically")
gate(sp.simplify(H_special_reduced) == 0,
     "special conic maps to the original irreducible quartic")
special_conic_matrix = sp.Matrix(((3, 0, 1), (0, 1, 0), (1, 0, -1)))
gate(special_conic_matrix.det() == -4,
     "special projective conic is smooth")
gate(sp.discriminant(vn**2+3, vn) != 0,
     "special conic has two distinct infinity points")

C_edge = sp.simplify(1/a-a)
D_edge = sp.simplify((1-omega**2)/a-a*(1-omega))
R = sp.symbols("R")
edge = sp.expand(R*(C_edge*R+D_edge)**2-4)

# At the unique generic projective infinity direction [X:t:z]=[1:0:0],
# put eta=t/X and zz=1/X.  The lower Newton edge has weight
# wt(eta)=1, wt(zz)=2.  Dividing its monomials by zz^3 and setting
# R=eta^2/zz gives exactly the displayed cubic.
eta, zz = sp.symbols("eta zz")
local_curve = sp.Poly(
    sp.cancel(zz**6*Hn.subs({X: 1/zz, t: eta/zz})), eta, zz
)
newton_terms = [
    (monomial, coefficient)
    for monomial, coefficient in local_curve.terms()
    if monomial[0]+2*monomial[1] == 6
]
gate(len(newton_terms) == 4,
     "nonlinear prototype has the expected four-term lower Newton edge")
gate(all(monomial[0] % 2 == 0 for monomial, _ in newton_terms),
     "Newton edge descends to R=eta^2/zz")
edge_from_curve = sp.expand(sum(
    coefficient*R**(monomial[0]//2)
    for monomial, coefficient in newton_terms
))
gate(sp.factor(edge_from_curve-edge) == 0,
     "nonlinear prototype Newton edge is R(CR+D)^2-4")
gate(sp.factor(homogeneous_part(Hn, 6)
               - t**6*(a**2-1)**2/a**2) == 0,
     "away from a^2=1 the prototype has one projective infinity direction")
gate(sp.expand(sp.diff(edge, R)
     -(C_edge*R+D_edge)*(3*C_edge*R+D_edge)) == 0,
     "nonlinear Newton-edge derivative factorization")
gate(edge.subs(R, 0) == -4,
     "nonlinear Newton edge has no zero root")
critical_second = sp.simplify(
    sp.diff(edge, R, 2).subs(R, -D_edge/(3*C_edge))
)
gate(sp.simplify(critical_second-2*C_edge*D_edge) == 0,
     "only possible repeated edge root has explicit second derivative")
gate(sp.factor(sp.discriminant(sp.diff(edge, R), R)
               - 4*C_edge**2*D_edge**2) == 0,
     "a triple edge root would force D=0 when C is nonzero")
gate(sp.factor(a*C_edge-(1-a**2)) == 0,
     "C=0 is exactly the separately handled a^2=1 support degeneration")


# Parallel-affine deformation toward reserved THM-3946.  Separating the two
# affine factors by c!=0 is just the proved two-ended c=1 quartic after
# scaling; their collision at c=0 is exactly THM-3944's doubled conductor.
Y, cpar, uu = sp.symbols("Y cpar uu")
fc = t*(t-cpar)
p0c = Y
p1c = Y+fc
F1c = fc+(1-omega)*Y
F2c = fc+(1-omega**2)*Y
Dc = sp.expand(t*F1c)
Sc = sp.expand((t-cpar)*F2c)
q0c = sp.expand(Sc-Dc)
q1c = sp.expand(Sc+Dc)
Hc = sp.expand(q0c**2-4*p0c**3)
gate(sp.simplify(Dc*Sc-(p1c**3-p0c**3)) == 0,
     "parallel-affine internal-split identity")
gate(sp.simplify(Hc-(q1c**2-4*p1c**3)) == 0,
     "parallel-affine family has a common discriminant")
scaled_Hc = sp.simplify(
    Hc.subs({Y: cpar**2*X, t: cpar*uu})/cpar**6
)
gate(sp.simplify(scaled_Hc-Hn.subs({a: 1, t: uu})) == 0,
     "every c!=0 member scales to the irreducible two-ended c=1 quartic")
gate(sp.expand(Hc.subs(cpar, 0)+Y**2*(4*Y+3*t**2)) == 0,
     "c=0 collision is the doubled-conductor THM3944 branch")


summary = {
    "checks": CHECKS,
    "scope": "independent THM3942 rederivation;one nonlinear balanced hostile only",
    "ufd_allocations": 8,
    "independent_types": "singleton conic two-place;extreme Fermat three-place",
    "dependent_type": "line-degenerate when p1^3!=p0^3",
    "conclusion": "affine whole-factor split closed;nonlinear balanced prototype pays branches",
    "nonlinear_hostile": "a^2=1 irreducible conic has two ends;generic edge has >=2 branches",
    "affine_internal_deformation": "c!=0 two ends;c=0 doubled conductor collision",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("Double-torus nonlinear balanced partial-split exact scout")
print(f"CHECKS={CHECKS}")
print("UFD_ALLOCATIONS=8;TYPES=singleton_or_extreme_up_to_complement")
print("SINGLETON=NORMALIZATION_CONIC;INFINITY_PLACES=2")
print("EXTREME=NORMALIZATION_FERMAT_CUBIC;INFINITY_PLACES=3")
print("DEPENDENT_GRADIENTS=LINE_DEGENERATE_UNDER_GENUINE_SPLIT")
print("AFFINE_LINEAR_ONE_PLACE=NONE_NONLINEAR")
print("NONLINEAR_PROTOTYPE=t(t-1)_balanced_split;A2=1_IRREDUCIBLE_CONIC_TWO_ENDS;GENERIC_BRANCHES>=2")
print("AFFINE_INTERNAL_DEFORMATION=c!=0_TWO_ENDS;c=0_DOUBLED_CONDUCTOR")
print(f"SEMANTIC_SHA256={semantic}")
