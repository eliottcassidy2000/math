#!/usr/bin/env python3
"""Exact companion for the affine point-cap classification in THM-2465."""

import sympy as sp


FAILED = []


def require(label, condition):
    ok = bool(condition)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        FAILED.append(label)


x, y, z = sp.symbols("x y z")


def affine(prefix):
    c0, cx, cy = sp.symbols(f"{prefix}0 {prefix}x {prefix}y")
    return c0 + cx * x + cy * y


def jac(u, v):
    return sp.diff(u, x) * sp.diff(v, y) - sp.diff(u, y) * sp.diff(v, x)


f = affine("f")
B1 = affine("b1")
B2 = affine("b2")
C1 = affine("c1")
C2 = affine("c2")
B3 = affine("b3")
C3 = affine("c3")

F = sp.Matrix((B1 * z + C1, B2 * z + C2, f * z**2 + B3 * z + C3))
detF = sp.expand(F.jacobian((x, y, z)).det())

bx = sp.Matrix((sp.diff(B1, x), sp.diff(B2, x)))
by = sp.Matrix((sp.diff(B1, y), sp.diff(B2, y)))
b = sp.Matrix((B1, B2))
j2 = lambda u, v: sp.det(sp.Matrix.hstack(u, v))
L1 = j2(by, b)
M1 = j2(bx, b)
N2 = j2(bx, by)
D3 = sp.expand(sp.diff(f, x) * L1 - sp.diff(f, y) * M1 + 2 * f * N2)

require("D3 is the z^3 coefficient of the point-cap Jacobian",
        sp.expand(sp.Poly(detF, z).coeff_monomial(z**3) - D3) == 0)
require("D3=-B1^4*Jac(f/B1^2,B2/B1) on B1!=0",
        sp.factor(sp.together(D3 + B1**4 * jac(f / B1**2, B2 / B1))) == 0)

# Rank-two normal form B1=x,B2=y.
aa, bb, cc = sp.symbols("aa bb cc")
f2 = aa * x + bb * y + cc
B12, B22 = x, y
L12 = sp.det(sp.Matrix.hstack(sp.Matrix((0, 1)), sp.Matrix((x, y))))
M12 = sp.det(sp.Matrix.hstack(sp.Matrix((1, 0)), sp.Matrix((x, y))))
N22 = jac(B12, B22)
D32 = sp.expand(sp.diff(f2, x) * L12 - sp.diff(f2, y) * M12 + 2 * f2 * N22)
require("rank-two normal form has D3=aa*x+bb*y+2*cc",
        sp.expand(D32 - (aa * x + bb * y + 2 * cc)) == 0)
require("rank-two D3 vanishes only for f=0",
        sp.Poly(D32, x, y).coeffs() == [aa, bb, 2 * cc])

# Rank-one normal form Bi=ai*ell+ci.
u, v = sp.symbols("u v")
a1, a2, c1, c2 = sp.symbols("a1 a2 c1 c2")
g0, gx, gy = sp.symbols("g0 gx gy")
ell = u * x + v * y
RB1 = a1 * ell + c1
RB2 = a2 * ell + c2
rf = g0 + gx * x + gy * y
delta = a1 * c2 - a2 * c1
rL = sp.det(sp.Matrix.hstack(
    sp.Matrix((sp.diff(RB1, y), sp.diff(RB2, y))), sp.Matrix((RB1, RB2))))
rM = sp.det(sp.Matrix.hstack(
    sp.Matrix((sp.diff(RB1, x), sp.diff(RB2, x))), sp.Matrix((RB1, RB2))))
rN = jac(RB1, RB2)
rD3 = sp.expand(sp.diff(rf, x) * rL - sp.diff(rf, y) * rM + 2 * rf * rN)
require("rank-one law D3=delta*(v*f_x-u*f_y)",
        sp.expand(rD3 - delta * (v * gx - u * gy)) == 0)

m, n = sp.symbols("m n")
require("every affine f=m*ell+n satisfies the rank-one varying-direction law",
        sp.expand(rD3.subs({gx: m * u, gy: m * v, g0: n})) == 0)

# If delta!=0 then (B1,B2) is a basis of span(ell,1).  Its symmetric square
# is therefore a basis of span(ell^2,ell,1); the determinant pins the exact
# nondegeneracy and proves existence/uniqueness of homogeneous Q.
sym2_change = sp.Matrix((
    (a1**2, a1 * a2, a2**2),
    (2 * a1 * c1, a1 * c2 + a2 * c1, 2 * a2 * c2),
    (c1**2, c1 * c2, c2**2),
))
require("Sym^2 basis determinant is delta^3",
        sp.expand(sym2_change.det() - delta**3) == 0)

# The target shear is a genuine polynomial automorphism and kills exactly
# the quadratic z coefficient Q(B1,B2).
P1, P2, P3 = sp.symbols("P1 P2 P3")
q20, q11, q02 = sp.symbols("q20 q11 q02")
Q = q20 * P1**2 + q11 * P1 * P2 + q02 * P2**2
shear = (P1, P2, P3 - Q)
inverse_shear = (P1, P2, P3 + Q)
composition = tuple(sp.expand(expr.subs(
    {P1: shear[0], P2: shear[1], P3: shear[2]}, simultaneous=True))
                    for expr in inverse_shear)
require("target quadratic shear has the displayed polynomial inverse",
        composition == (P1, P2, P3))
Q_on_F12 = Q.subs({P1: B1 * z + C1, P2: B2 * z + C2})
require("target shear changes the z^2 coefficient by f-Q(B1,B2)",
        sp.expand(sp.Poly(F[2] - Q_on_F12, z).coeff_monomial(z**2)
                  - (f - Q.subs({P1: B1, P2: B2}))) == 0)

# Complete the varying-direction branch after normalizing its first two
# z-coefficients to (1,x).  These are differential-polynomial identities, so
# C1,C2,C3,a below are arbitrary polynomial functions rather than a degree
# ansatz.
NC1 = sp.Function("NC1")(x, y)
NC2 = sp.Function("NC2")(x, y)
NC3 = sp.Function("NC3")(x, y)
na = sp.Function("na")(x, y)
NG = sp.Matrix((z + NC1, x * z + NC2, na * z + NC3))
NP1 = z + NC1
nh = NC2 - x * NC1
nk = NC3 - na * NC1
normalized_identity = ((NP1 + sp.diff(nh, x))
                       * (sp.diff(na, y) * NP1 + sp.diff(nk, y))
                       - sp.diff(nh, y)
                       * (sp.diff(na, x) * NP1 + sp.diff(nk, x)))
require("normalized (1,x,a) z-affine determinant identity",
        sp.simplify(NG.jacobian((x, y, z)).det() - normalized_identity) == 0)
require("normalized top coefficient forces a_y=0",
        sp.simplify(sp.Poly(normalized_identity, z).coeff_monomial(z**2)
                    - sp.diff(na, y)) == 0)

ah = sp.Function("ah")(x)
hh = sp.Function("hh")(x, y)
bh = sp.Function("bh")(x)
kh = sp.diff(ah, x) * hh + bh
require("integrated normalized determinant is -h_y*(a''*h+b')",
        sp.simplify(jac(hh, kh)
                    + sp.diff(hh, y)
                    * (sp.diff(ah, x, 2) * hh + sp.diff(bh, x))) == 0)

# A nontrivial member of the resulting universal automorphism family checks
# the explicit inverse, including nonlinear C1 and e(x).
alpha, beta, dd, gamma, delta0 = sp.symbols("alpha beta dd gamma delta0", nonzero=True)
test_C1 = x**2 + x * y + y**2 + 1
test_e = x**2 - 2 * x + 3
test_a = alpha * x + beta
test_h = dd * y + test_e
test_k = alpha * test_h + gamma * x + delta0
test_C2 = x * test_C1 + test_h
test_C3 = test_a * test_C1 + test_k
test_G = sp.Matrix((z + test_C1, x * z + test_C2, test_a * z + test_C3))
require("normalized varying-direction family has determinant -d*gamma",
        sp.simplify(test_G.jacobian((x, y, z)).det() + dd * gamma) == 0)

TP1, TP2, TP3 = sp.symbols("TP1 TP2 TP3")
inv_x = (TP3 - beta * TP1 - alpha * TP2 - delta0) / gamma
inv_y = (TP2 - inv_x * TP1 - test_e.subs(x, inv_x)) / dd
inv_z = TP1 - test_C1.subs({x: inv_x, y: inv_y}, simultaneous=True)
normalized_round_trip = tuple(sp.factor(component.subs(
    {x: inv_x, y: inv_y, z: inv_z}, simultaneous=True) - target)
                              for component, target in zip(test_G, (TP1, TP2, TP3)))
require("normalized varying-direction family has the displayed polynomial inverse",
        normalized_round_trip == (0, 0, 0))

# Rank-zero constant-projective direction.  In the triangular source
# coordinate s=P1, the complete determinant is a quadratic pencil of planar
# Jacobian responses.
RC2 = sp.Function("RC2")(x, y)
Rf = sp.Function("Rf")(x, y)
RG = sp.Function("RG")(x, y)
RR = sp.Function("RR")(x, y)
rank_zero_map = sp.Matrix((z, RC2, Rf * z**2 + RG * z + RR))
rank_zero_response = (jac(RC2, Rf) * z**2
                      + jac(RC2, RG) * z + jac(RC2, RR))
require("rank-zero constant-direction determinant is the planar response pencil",
        sp.simplify(rank_zero_map.jacobian((x, y, z)).det()
                    - rank_zero_response) == 0)

# Nonconstant affine f=x forces the explicit triangular family in the proof.
d0, r0, c0 = sp.symbols("d0 r0 c0", nonzero=True)
rank_C1 = x * y + y**2 + 1
rank_G = x**2 + 1
rank_E = x**3 - x + 2
rank_C2 = d0 * x + c0
rank_R = r0 * y + rank_E
rank_B3 = rank_G + 2 * x * rank_C1
rank_C3 = rank_R + rank_G * rank_C1 + x * rank_C1**2
rank_F = sp.Matrix((z + rank_C1, rank_C2,
                    x * z**2 + rank_B3 * z + rank_C3))
require("rank-zero nonconstant-cap family has determinant d*r",
        sp.simplify(rank_F.jacobian((x, y, z)).det() - d0 * r0) == 0)

RP1, RP2, RP3 = sp.symbols("RP1 RP2 RP3")
rank_inv_x = (RP2 - c0) / d0
rank_inv_y = (RP3 - rank_inv_x * RP1**2
              - rank_G.subs(x, rank_inv_x) * RP1
              - rank_E.subs(x, rank_inv_x)) / r0
rank_inv_z = RP1 - rank_C1.subs(
    {x: rank_inv_x, y: rank_inv_y}, simultaneous=True)
rank_round_trip = tuple(sp.factor(component.subs(
    {x: rank_inv_x, y: rank_inv_y, z: rank_inv_z}, simultaneous=True) - target)
                        for component, target in zip(rank_F, (RP1, RP2, RP3)))
require("rank-zero nonconstant-cap family has the displayed polynomial inverse",
        rank_round_trip == (0, 0, 0))

# Constant f and a centralizer term H(C2) are removed by two target shears,
# leaving the direct suspension of a planar Keller pair.
fconst = sp.symbols("fconst", nonzero=True)
plane_U = x + y**2
plane_R = y
plane_H = plane_U**2 + 1
susp_C1 = x * y + 1
susp_B3 = plane_H + 2 * fconst * susp_C1
susp_C3 = plane_R + plane_H * susp_C1 + fconst * susp_C1**2
susp_F = sp.Matrix((z + susp_C1, plane_U,
                    fconst * z**2 + susp_B3 * z + susp_C3))
require("constant-cap suspension control is Keller",
        sp.simplify(susp_F.jacobian((x, y, z)).det() - 1) == 0)
require("constant-cap target shears recover the planar suspension",
        sp.expand(susp_F[2] - fconst * susp_F[0]**2
                  - (susp_F[1]**2 + 1) * susp_F[0] - plane_R) == 0)

# Rank-one constant-projective scale: normalize (B1,B2)=(x,0) and reproduce
# the full coefficient ledger used by the emptiness proof.
SC1 = sp.Function("SC1")(x, y)
SC2 = sp.Function("SC2")(x, y)
Sf = sp.Function("Sf")(x, y)
SB3 = sp.Function("SB3")(x, y)
SC3 = sp.Function("SC3")(x, y)
scale_F = sp.Matrix((x * z + SC1, SC2, Sf * z**2 + SB3 * z + SC3))
scale_expected = ((2 * Sf * z + SB3)
                  * (z * sp.diff(SC2, y) + jac(SC1, SC2))
                  + x * (z**2 * jac(SC2, Sf)
                         + z * jac(SC2, SB3) + jac(SC2, SC3)))
require("rank-one scale determinant coefficient ledger",
        sp.simplify(scale_F.jacobian((x, y, z)).det() - scale_expected) == 0)

sa, sc = sp.symbols("sa sc")
sb = sp.symbols("sb", nonzero=True)
scale_affine_f = sa * x + sb * y + sc
scale_D = lambda expr: (sb * x * sp.diff(expr, x)
                        + (sa * x + 2 * sb * y + 2 * sc) * sp.diff(expr, y))
scale_Y = y + sc / sb + (sa / sb) * x
require("f_y!=0 scale PDE conjugates to positive Euler weights (1,2)",
        sp.simplify(scale_D(x) - sb * x) == 0
        and sp.simplify(scale_D(scale_Y) - 2 * sb * scale_Y) == 0)

scale_a, scale_c = sp.symbols("scale_a scale_c")
scale_C1 = sp.Function("scale_C1")(x, y)
scale_C2x = sp.Function("scale_C2x")(x)
scale_B3 = sp.Function("scale_B3")(x, y)
scale_C3 = sp.Function("scale_C3")(x, y)
scale_b0_F = sp.Matrix((x * z + scale_C1, scale_C2x,
                        (scale_a * x + scale_c) * z**2
                        + scale_B3 * z + scale_C3))
scale_b0_det = sp.Poly(sp.expand(scale_b0_F.jacobian((x, y, z)).det()), z)
require("f_y=0 scale ledger factors through C2'(x)",
        sp.simplify(scale_b0_det.coeff_monomial(z)
                    - sp.diff(scale_C2x, x)
                    * (x * sp.diff(scale_B3, y)
                       - 2 * (scale_a * x + scale_c)
                       * sp.diff(scale_C1, y))) == 0
        and sp.simplify(scale_b0_det.coeff_monomial(1)
                        - sp.diff(scale_C2x, x)
                        * (x * sp.diff(scale_C3, y)
                           - scale_B3 * sp.diff(scale_C1, y))) == 0)

scale_h = sp.Function("scale_h")(x)
scale_k = sp.Function("scale_k")(x)
scale_r = sp.symbols("scale_r", nonzero=True)
scale_integral = (scale_a * scale_C1**2 + scale_h * scale_C1
                  + scale_r * y + scale_k)
require("rank-one scale equations integrate to the x-divisibility obstruction",
        sp.simplify(sp.diff(scale_integral, y)
                    - ((2 * scale_a * scale_C1 + scale_h)
                       * sp.diff(scale_C1, y) + scale_r)) == 0)

# Minimal constant-direction boundary survivor.  It is Keller and a genuine
# point-cap, but no quadratic shear in its first two targets removes y*z^2.
survivor = sp.Matrix((z, y, x + y * z**2))
require("constant-direction survivor (z,y,x+y*z^2) is Keller",
        sp.expand(survivor.jacobian((x, y, z)).det()) == -1)
A_survivor = tuple(sp.Poly(component, z).coeff_monomial(z**2) for component in survivor)
B_survivor = tuple(sp.Poly(component, z).coeff_monomial(z) for component in survivor[:2, :])
require("survivor has point-cap A=(0,0,y) and constant b=(1,0)",
        A_survivor == (0, 0, y) and B_survivor == (1, 0))
require("no homogeneous quadratic Q(B1,B2) equals the survivor coefficient y",
        sp.expand(y - Q.subs({P1: 1, P2: 0})) != 0)

a, btarget, c = sp.symbols("a btarget c")
inverse = sp.Matrix((c - btarget * a**2, btarget, a))
round_trip = tuple(sp.expand(expr.subs(
    {x: inverse[0], y: inverse[1], z: inverse[2]}, simultaneous=True))
                   for expr in survivor)
require("constant-direction survivor is degree one, with explicit inverse",
        round_trip == (a, btarget, c))

print("FAILED CHECKS:", FAILED if FAILED else "NONE")
assert not FAILED
