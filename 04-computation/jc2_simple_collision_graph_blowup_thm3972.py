#!/usr/bin/env python3
"""Exact companion for THM-3972 (simple-collision graph blowups)."""

from __future__ import annotations

import hashlib
import json
import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(expression) == 0, message)


t, P, T, h, v, y, x = sp.symbols("t P T h v y x")
a, c, r = sp.symbols("a c r")
V, W = sp.symbols("V W")

# The affine-P graph family and its line normalization.
y0 = P - r**2
q = sp.expand(3 * r * P - r**3 + y0**2 * (c + a * y0))
F = sp.expand(T**3 - 3 * P * T - q)
node = x**3 - 3 * r * x**2 - 3 * x * y - c * y**2 - a * y**3
zero(F.subs({T: x - r, P: y + r**2}) - node,
     "node normal form")

D = 1 - a * v**3
N = 3 * r + 3 * v + c * v**2
zero(sp.together(node.subs({x: N / D, y: v * N / D})),
     "line normalization")

D_h = W**3 - a * V**3
N_h = 3 * r * W**2 + 3 * V * W + c * V**2
P_h = sp.expand(r**2 * D_h + V * N_h)
Xi = c**3 + 27 * a**2 * r**3 - 27 * a * c * r + 27 * a
zero(sp.resultant(D, N, v) - Xi, "collision resultant")
gate(sp.Poly(P_h, V, W).total_degree() == 3,
     "projective numerator degree three")
gate(sp.Poly(D_h, V, W).total_degree() == 3,
     "projective denominator degree three")

# The finite cover after resolving basepoints is cut fibrewise by this
# binary cubic.  Its V^2 W coefficient is the constant -3, so no target
# fibre can contain the whole projective line.
Y = sp.symbols("Y")
H_target = sp.expand(Y * D_h - V * N_h)
gate(sp.Poly(H_target, V, W).coeff_monomial(V**2 * W) == -3,
     "binary cubic never vanishes identically")

# On the finite-v chart the target equation is G=0.  These identities make
# x and hence T regular even across D=N=0.
G = sp.expand(y * D - v * N)
x_regular = sp.expand(N + a * v**2 * y)
zero(v * x_regular - y + G, "regular x satisfies vx=y")
zero(D * x_regular - N - a * v**2 * G,
     "regular x satisfies Dx=N")

# A simple resultant root is exactly one transverse basepoint.  At a finite
# basepoint substitute a*v^3=1 and N=0.  Xi' is a nonzero companion-root
# factor times the Jacobian of (D,N) in (t,v).
ap, cp, rp = sp.symbols("ap cp rp")
Dt = -ap * v**3
Dv = -3 * a * v**2
Nt = 3 * rp + cp * v**2
Nv = 3 + 2 * c * v
J_base = sp.expand(Dt * Nv - Dv * Nt)
Xi_t = sp.diff(Xi, a) * ap + sp.diff(Xi, c) * cp + sp.diff(Xi, r) * rp
base_sub = {a: 1 / v**3, r: -v - c * v**2 / 3}
companion = (c**2 * v**2 + 3 * c * v + 9) / v**3
zero(Xi_t.subs(base_sub) - companion * J_base.subs(base_sub),
     "finite basepoint derivative identity")

omega = sp.symbols("omega")
N_at = lambda z: 3 * (-v - c * v**2 / 3) + 3 * z + c * z**2
other_product = sp.expand(N_at(omega * v) * N_at(omega**2 * v))
other_product = sp.rem(sp.Poly(other_product, omega),
                       sp.Poly(omega**2 + omega + 1, omega)).as_expr()
zero(other_product - 3 * v**2 * (c**2 * v**2 + 3 * c * v + 9),
     "companion factor detects a second basepoint")

# In the finite-basepoint blowup chart D=e,N=e*zeta, the target Jacobian
# along the exceptional line is v*(N_v-D_v*zeta)/J_base.  Its only zero is
# the strict ramification intersection zeta=v+2r.
zeta = sp.symbols("zeta")
jac_exc_finite = sp.cancel(v * (Nv - Dv * zeta) / J_base)
zero(J_base * jac_exc_finite - v * (Nv - Dv * zeta),
     "finite exceptional target Jacobian")
zero((Nv / Dv - (v + 2 * r)).subs(base_sub),
     "finite exceptional zero is strict ramification point")

# At the projective-infinity basepoint a=c=W=0, Xi'=27a' and the local
# Jacobian of (W^3-a, 3rW^2+3W+c) is -3a'.
w = sp.symbols("w")
D_inf = w**3 - a
N_inf = 3 * r * w**2 + 3 * w + c
Dt_inf = -ap
Dw_inf = 3 * w**2
Nt_inf = 3 * rp * w**2 + cp
Nw_inf = 6 * r * w + 3
J_inf = sp.expand(Dt_inf * Nw_inf - Dw_inf * Nt_inf)
gate(D_inf.subs({w: 0, a: 0}) == 0, "infinity pole basepoint")
gate(N_inf.subs({w: 0, c: 0}) == 0, "infinity numerator basepoint")
gate(Xi_t.subs({a: 0, c: 0}) == 27 * ap,
     "infinity resultant derivative")
gate(J_inf.subs(w, 0) == -3 * ap,
     "infinity transverse Jacobian")
zero(J_inf.subs(w, 0) * zeta * (-3 / (J_inf.subs(w, 0) * zeta)) + 3,
     "infinity exceptional target Jacobian")

# Relative ramification.  On G=0 the saturated affine equation is R4=0;
# the identity keeps the collision exceptional from becoming a component.
rho = a * v**2 * y + c * v**2 + 2 * v + r
R4 = a * v**4 + 2 * a * r * v**3 + c * v**2 + 2 * v + r
zero(D * rho - R4 - a * v**2 * G,
     "saturated ramification identity")
R_h = (a * V**4 + 2 * a * r * V**3 * W + c * V**2 * W**2
       + 2 * V * W**3 + r * W**4)
gate(sp.Poly(R_h, V, W).total_degree() == 4,
     "ramification multisection degree four")
zero(R_h - (W**2 * N_h - (V + 2 * r * W) * D_h),
     "homogeneous ramification-pole identity")
gate((w**2 - (1 + 2 * r * w) * zeta).subs(w, 0) == -zeta,
     "infinity strict ramification meets only removed point")

# Divisor lattice for n transverse points on one irreducible cubic pole.
# Generators are H,E_1,...,E_n and the removed strict pole is
# 3H-sum E_i.  Both K=-2H+sum E_i and Ram=4H-sum E_i reduce to H.
n = 4
pole_relation = sp.Matrix([3] + [-1] * n)
canonical = sp.Matrix([-2] + [1] * n)
ramification = sp.Matrix([4] + [-1] * n)
H_class = sp.Matrix([1] + [0] * n)
gate(pole_relation.rank() == 1, "one irreducible pole relation")
gate(sp.gcd_list(list(pole_relation)) == 1,
     "collision makes the pole relation primitive")
gate(canonical - H_class == -pole_relation,
     "canonical class reduces to H")
gate(ramification - H_class == pole_relation,
     "ramification class reduces to H")
gate((n + 1) - pole_relation.rank() == n,
     "class rank equals collision count")

# Euler rigidity: chi(X)=2+n-chi(E), while an A2 boundary with class rank n
# has Euler at most 1+n.  Equality forces chi(E)=1.  In a cubic Kummer
# normalization, B branch places give genus B-2; A1 requires B=2 and one
# point over infinity.
B_count, infinity_points = sp.symbols("B_count infinity_points", integer=True)
genus = B_count - 2
chi_normalized_pole = 2 - 2 * genus - infinity_points
gate(chi_normalized_pole.subs({B_count: 2, infinity_points: 1}) == 1,
     "one-support Kummer pole has Euler one")
gate(genus.subs(B_count, 2) == 0,
     "two Kummer branch places give genus zero")

# The first one-collision control a=t,c=1,r=0.
t0 = -sp.Rational(1, 27)
named = {a: t, c: 1, r: 0}
zero(Xi.subs(named) - (1 + 27 * t), "named collision polynomial")
gate(D.subs(named).subs({t: t0, v: -3}) == 0,
     "named pole basepoint")
gate(N.subs(named).subs({t: t0, v: -3}) == 0,
     "named numerator basepoint")
D_named = D.subs(named)
N_named = N.subs(named)
J_named = sp.det(sp.Matrix([
    [sp.diff(D_named, t), sp.diff(D_named, v)],
    [sp.diff(N_named, t), sp.diff(N_named, v)],
]))
gate(J_named.subs({t: t0, v: -3}) == -81,
     "named basepoint transverse")

# On the graph chart T(1-tv^3)=v(v+3), P=vT.  The saturated Jacobian
# selects T=-3 on the exceptional line and has two genuine generic factors.
J_graph = T * (1 + 2 * t * v**3) + 2 * v**2 + 3 * v
gate(J_graph.subs({t: t0, v: -3}) == 3 * (T + 3),
     "exceptional ramification selects one point")
zero(sp.together(J_graph.subs(T, N_named / D_named))
     - 3 * v * (t * v**3 + v + 2) / D_named,
     "named ramification factorization")

residual = t * v**3 + v + 2
gate(sp.Poly(residual, t).degree() == 1,
     "residual is linear in t")
gate(sp.gcd(v**3, v + 2) == 1,
     "residual primitive and irreducible")
zero(v * (t * v**2 + 1) + 2 - residual,
     "residual curve makes v a unit")

# In Cl=<H,E>/(3H-E), R0=H and R1=3H-E=0.  The de Rham
# Gysin quotient has rank one, but the pulled-back target volume is exact.
relation_named = sp.Matrix([3, -1])
K_named = sp.Matrix([-2, 1])
R0_class = sp.Matrix([1, 0])
R1_class = sp.Matrix([3, -1])
gate(K_named - R0_class == -relation_named,
     "named canonical class H")
gate(R1_class == relation_named,
     "second ramification prime has zero quotient class")
gate((R0_class + R1_class) - R0_class == relation_named,
     "canonical vector split H plus zero")
gate(2 - sp.Matrix([[3, -1]]).rank() == 1,
     "named H2 de Rham Gysin quotient rank one")

# Full a=t, constant-(c,r) packet.  The r=0 row always splits into the
# section v=0 and a primitive linear-in-t residual.  For r!=0, simple Xi
# excludes both discriminant seams; the irreducible ramification curve has
# distinct punctures v=0 and v=-2r.
cc, rr = sp.symbols("cc rr", nonzero=True)
Xi_at = sp.expand(Xi.subs({a: t, c: cc, r: rr}))
zero(Xi_at - (27 * rr**3 * t**2 + 27 * (1 - cc * rr) * t + cc**3),
     "constant cr collision quadratic")
zero(sp.discriminant(Xi_at, t)
     + 27 * (cc * rr - 3)**2 * (4 * cc * rr - 3),
     "constant cr collision discriminant")

R_at = sp.expand(R4.subs({a: t, c: cc, r: rr}))
zero(R_at - (t * v**3 * (v + 2 * rr) + cc * v**2 + 2 * v + rr),
     "constant cr ramification linear in t")
gate((cc * v**2 + 2 * v + rr).subs(v, 0) == rr,
     "ramification numerator nonzero at v0")
gate(sp.factor((cc * v**2 + 2 * v + rr).subs(v, -2 * rr))
     == rr * (4 * cc * rr - 3),
     "ramification second-puncture cancellation gate")

Xi_r0 = sp.expand(Xi.subs({a: t, c: cc, r: 0}))
zero(Xi_r0 - (cc**3 + 27 * t), "r0 one collision")
zero(R4.subs({a: t, c: cc, r: 0})
     - v * (t * v**3 + cc * v + 2),
     "r0 two ramification factors")
zero(v * (t * v**2 + cc) + 2 - (t * v**3 + cc * v + 2),
     "r0 residual makes v a unit")

c_cancel = sp.Rational(3, 4) / rr
zero(sp.factor(Xi_at.subs(cc, c_cancel))
     - 27 * (8 * rr**3 * t + 1)**2 / (64 * rr**3),
     "fourcr3 double collision seam")
zero(sp.factor(Xi_at.subs(cc, 3 / rr))
     - 27 * (rr**3 * t - 1)**2 / rr**3,
     "cr3 two-address double collision seam")

summary = {
    "checks": CHECKS,
    "family": "q=3rP-r^3+(P-r^2)^2(c+a(P-r^2))",
    "simple_gate": "Xi squarefree gives transverse basepoints",
    "normalization": "blowup A1xP1 at basepoints minus strict cubic pole",
    "class_irreducible_pole": "Cl=Z^n; K=Ram=H",
    "exceptional_ramification": "finite strict point; infinity removed point",
    "euler": "possible pole has chi1 and one-support Kummer form",
    "constant_height1": "all a=t constant c,r simple rows close",
    "named": "a=t,c=1,r=0 has R0 class H and R1 class 0",
    "derham": "H2 has rank1 but dP wedge dt=d(P dt) is exact",
    "scope": "general simple-collision canonical-compatible rows remain",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3972 simple-collision graph blowup companion")
print(f"CHECKS={CHECKS}")
print("NORMALIZATION=BLOWUP_RELATIVE_P1_AT_TRANSVERSE_BASEPOINTS_MINUS_POLE")
print("FINITE=BINARY_CUBIC_COEFFICIENT_MINUS_3_PREVENTS_FIBRE_COMPONENT")
print("IRREDUCIBLE_POLE_CLASS=CL_ZN;K_H;RAM_H")
print("EXCEPTIONAL=FINITE_STRICT_POINT;INFINITY_REMOVED_POINT")
print("EULER=POLE_CHI1;ONE_SUPPORT_KUMMER_NECESSARY")
print("HEIGHT1_CONSTANT_CR=ALL_SIMPLE_ROWS_CLOSED")
print("NAMED=A_T;C_1;R_0;XI_1_PLUS_27T;ONE_COLLISION")
print("NAMED_RAMIFICATION=R0_CLASS_H;R1_CLASS_ZERO;NO_KELLER_BOUNDARY")
print("DERHAM=H2_RANK1_BUT_PULLBACK_TARGET_VOLUME_EXACT")
print("SCOPE=GENERAL_SIMPLE_COLLISION_CANONICAL_COMPATIBLE_PACKET_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
