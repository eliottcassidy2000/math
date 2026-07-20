#!/usr/bin/env python3
"""jacobian_counterexample_verify_klein_S323.py -- klein-2026-07-19-S323.

INDEPENDENT EXACT VERIFICATION of the claimed counterexample to the Jacobian
Conjecture (owner-supplied, reported found recently):

  F(x,y,z) = ( (1+xy)^3 z + y^2 (1+xy)(4+3xy),
               y + 3x (1+xy)^2 z + 3x y^2 (4+3xy),
               2x - 3x^2 y - x^3 z )        : C^3 -> C^3

Claims checked, all in exact arithmetic (sympy):
  (V1) det Jacobian(F) is the CONSTANT -2  (Keller map).
  (V2) F(0,0,-1/4) = F(1,-3/2,13/2) = F(-1,3/2,13/2) = (-1/4,0,0)
       with the three preimages pairwise distinct  => F not injective.
  (V3) The full fiber F^{-1}(-1/4,0,0) via Groebner basis -- how many points?
  (V4) Generic-fiber degree probe: fiber over a random rational point.

(V1)+(V2) together refute the Jacobian Conjecture in dimension 3 (a polynomial
automorphism is injective; a Keller map that is not injective cannot be an
automorphism, and JC asserts every Keller map is one).
"""
import sympy as sp

x, y, z = sp.symbols('x y z')
w = 1 + x*y
F = [
    w**3 * z + y**2 * w * (4 + 3*x*y),
    y + 3*x * w**2 * z + 3*x * y**2 * (4 + 3*x*y),
    2*x - 3*x**2*y - x**3*z,
]

print("== (V1) Jacobian determinant ==")
J = sp.Matrix([[sp.expand(sp.diff(f, v)) for v in (x, y, z)] for f in F])
det = sp.expand(J.det())
print("det J(F) =", det)
assert det == -2, "DET IS NOT CONSTANT -2 -- CLAIM FAILS"
print("CONFIRMED: det == -2 identically (Keller map).")

print("\n== (V2) the three-point collision ==")
pts = [(0, 0, sp.Rational(-1, 4)),
       (1, sp.Rational(-3, 2), sp.Rational(13, 2)),
       (-1, sp.Rational(3, 2), sp.Rational(13, 2))]
imgs = []
for P in pts:
    img = tuple(sp.nsimplify(sp.expand(f.subs({x: P[0], y: P[1], z: P[2]}))) for f in F)
    imgs.append(img)
    print(f"F{P} = {img}")
assert len(set(pts)) == 3, "preimages not distinct"
assert imgs[0] == imgs[1] == imgs[2] == (sp.Rational(-1, 4), 0, 0), "COLLISION FAILS"
print("CONFIRMED: three distinct points map to (-1/4, 0, 0) -- F is NOT injective.")

print("\n== (V3) the full fiber over (-1/4, 0, 0) ==")
eqs = [sp.expand(F[0] + sp.Rational(1, 4)), sp.expand(F[1]), sp.expand(F[2])]
sols = sp.solve(eqs, [x, y, z], dict=True)
print(f"solve returned {len(sols)} solution(s):")
for s in sols:
    print("  ", {k: sp.nsimplify(v) for k, v in s.items()})

print("\n== (V4) generic fiber probe at target (1, 1, 1) ==")
eqs2 = [sp.expand(F[0] - 1), sp.expand(F[1] - 1), sp.expand(F[2] - 1)]
try:
    G = sp.groebner(eqs2, x, y, z, order='lex')
    last = G.exprs[-1]
    uni = sp.Poly(last, z) if last.free_symbols == {z} else None
    if uni is not None:
        print(f"lex Groebner eliminant in z has degree {uni.degree()} "
              f"(counts fiber points with multiplicity over C)")
    else:
        print("eliminant not univariate in z; exprs tail:", last)
except Exception as e:
    print("groebner probe failed:", e)

print("\nVERDICT: Jacobian Conjecture (n=3) REFUTED by this map, "
      "modulo the checks above (all exact).")
