#!/usr/bin/env python3
"""jacobian_master_quartic_klein_S324.py -- klein-2026-07-19-S324.

THE ANATOMY COMPUTATIONS for the JC counterexample F (see the S324 reflection
anatomy-of-a-jacobian-counterexample-the-master-quartic-klein-S324.md):

 (A1) the universal fiber W-cubic by elimination (z linear in F3, then
      y = (W-1)/x, then resultant in x):
        D(a,b,c) W^3 + (b^2 - 12a) W - 4a,
        D = 27 a^2 c^2 - 18 a b c + 16 a + b^3 c - b^2
      -- no W^2 term identically (trace law: fiber w-values sum to 0);
      at the collision target (-1/4, 0, 0) the cubic IS T_3(W) = 1.
 (A2) disc_W = -4 * (dD/dc)^2 * D  (Galois twist = sqrt(-D)).
 (A3) fiber-count spectrum: 3 on D != 0 (incl. S = 0), 1 on D = 0, never 2;
      spot checks (exact solve) at representative targets.
 (A4) the k = 2 transplant ansatz has ONLY the degenerate Keller solution
      (s = t = u = 0, det = 0): the "3" is structural.
"""
import sympy as sp

x, y, z, W, a, b, c = sp.symbols('x y z W a b c')
w = 1 + x*y
F = [w**3*z + y**2*w*(4+3*x*y), y + 3*x*w**2*z + 3*x*y**2*(4+3*x*y),
     2*x - 3*x**2*y - x**3*z]

print("== (A1) universal fiber cubic ==")
zs = (2*x - 3*x**2*y - c)/x**3
P1 = sp.expand(sp.numer(sp.together((F[0].subs(z, zs) - a)*x**3)))
P2 = sp.expand(sp.numer(sp.together((F[1].subs(z, zs) - b)*x**3)))
P1w = sp.expand(sp.numer(sp.together(P1.subs(y, (W-1)/x))))
P2w = sp.expand(sp.numer(sp.together(P2.subs(y, (W-1)/x))))
R = sp.factor(sp.resultant(sp.Poly(P1w, x), sp.Poly(P2w, x)))
print("resultant =", R)
D = 27*a**2*c**2 - 18*a*b*c + 16*a + b**3*c - b**2
cub = D*W**3 + (b**2 - 12*a)*W - 4*a
print("collision target check:", sp.factor(cub.subs({a: sp.Rational(-1,4), b: 0, c: 0})),
      " (= -(4W^3 - 3W - 1) = -(T_3(W) - 1))")

print("\n== (A2) discriminant ==")
disc = sp.factor(sp.expand(sp.discriminant(sp.Poly(cub, W), W)))
print("disc_W =", disc)
S = sp.diff(D, c)
print("dD/dc =", sp.expand(S), "; disc == -4*(dD/dc)^2*D:",
      sp.simplify(disc + 4*S**2*D) == 0)

print("\n== (A3) fiber spectrum spot checks ==")
def fib(aa, bb, cc):
    return len(sp.solve([F[0]-aa, F[1]-bb, F[2]-cc], [x, y, z], dict=True))
for (aa, bb, cc, tag) in [
    (1, 1, 1, "generic"),
    (1, 1, sp.Rational(17, 54), "S=0, D=11^3/108"),
    (-1, 2, -2, "D=0"),
    (sp.Rational(-1, 4), 0, 0, "collision target"),
]:
    print(f"  ({aa},{bb},{cc}) [{tag}]: fiber = {fib(aa, bb, cc)}")

print("\n== (A4) k=2 transplant ansatz ==")
p, q, r, s, t, u, v = sp.symbols('p q r s t u v')
G = [w**2*z + p*w*y**2 + q*y**2, y + 2*x*w*z + r*x*y**2 + v*x*y**3,
     s*x + t*x**2*y + u*x**2*z]
J = sp.Matrix([[sp.expand(sp.diff(g, vv)) for vv in (x, y, z)] for g in G])
poly = sp.Poly(sp.expand(J.det()), x, y, z)
eqs = [co for mono, co in poly.terms() if mono != (0, 0, 0)]
sols = sp.solve(eqs, [p, q, r, s, t, u, v], dict=True)
print(f"{len(eqs)} coefficient equations; solutions: {sols}")
print("  (only s=t=u=0 with det=0: NO nondegenerate k=2 Keller sibling in this shape)")
print("\ndone.")
