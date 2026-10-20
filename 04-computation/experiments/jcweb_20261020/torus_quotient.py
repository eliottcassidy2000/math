#!/usr/bin/env python3
"""Mixed-weight torus structure of the dimension-three counterexamples.
(A) Alpoge's Keller map F (THM-1300): det JF = -2; equivariance F(lambda.x) = lambda'.F(x) with source weights (1,-1,-2),
    target weights (-2,-1,1); invariant rings C[v,t] (v = xy, t = x^2 z) and C[v',t'] (v' = F2 F3, t' = F1 F3^2);
    the induced quotient map Fbar: C^2 -> C^2, its Jacobian, and the collision descended to the quotient.
(B) The Cima-van den Essen-Gasull-Hubbers-Manosas Markus-Yamabe counterexample field X = (-x + z d^2, -y - d^2, -z),
    d = x + yz: JX = -I + N with N nilpotent; equivariant for weights (1,2,-1); its Keller companion H = x + (X + x) =
    (x + z d^2, y - d^2, z) is a polynomial automorphism (d and z are preserved); invariant ring C[xz, yz^2]."""
import sympy as sp
x, y, z, lam = sp.symbols('x y z lambda')
u = 1 + x*y
F = sp.Matrix([u**3*z + y**2*u*(4 + 3*x*y), y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y), 2*x - 3*x**2*y - x**3*z])
J = F.jacobian([x, y, z]); detJ = sp.expand(J.det()); print("(A) det JF =", detJ)
Fl = F.subs({x: lam*x, y: y/lam, z: z/lam**2}, simultaneous=True)
print("    equivariance (1,-1,-2) -> (-2,-1,1):", [sp.simplify(Fl[i]*lam**w - F[i]) == 0 for i, w in enumerate((2, 1, -1))])
pts = [(0, 0, sp.Rational(-1, 4)), (1, sp.Rational(-3, 2), sp.Rational(13, 2)), (-1, sp.Rational(3, 2), sp.Rational(13, 2))]
print("    collision:", [tuple(F.subs({x: a, y: b, z: c})) for a, b, c in pts])
v, t = sp.symbols('v t')
vp = sp.expand(F[1]*F[2]); tp = sp.expand(F[0]*F[2]**2)
# express in v = xy, t = x^2 z: substitute y = v/x, z = t/x^2 and check no x remains
def descend(expr):
    e = sp.expand(expr.subs({y: v/x, z: t/x**2}, simultaneous=True))
    assert sp.simplify(sp.diff(e, x)) == 0, "not an invariant"
    return sp.expand(e)
Fbar = sp.Matrix([descend(vp), descend(tp)])
print("    quotient map Fbar(v,t) = (F2 F3, F1 F3^2):")
print("      v' =", sp.factor(Fbar[0])); print("      t' =", sp.factor(Fbar[1]))
Jb = sp.factor(Fbar.jacobian([v, t]).det()); print("    det J Fbar =", Jb)
print("    descended collision: Fbar at (v,t) of the three points:", [tuple(Fbar.subs({v: a*b, t: a**2*c})) for a, b, c in pts], " (source (v,t):", [(a*b, a**2*c) for a, b, c in pts], ")")
print("    degrees of Fbar components in (v,t):", [sp.Poly(Fbar[i], v, t).total_degree() for i in range(2)])
# fibre of Fbar over (0,0): solve
sols = sp.solve([Fbar[0], Fbar[1]], [v, t], dict=True); print("    Fbar^{-1}(0,0) =", sols)
# (B) Markus-Yamabe field
d = x + y*z
X = sp.Matrix([-x + z*d**2, -y - d**2, -z]); JX = X.jacobian([x, y, z]); N = JX + sp.eye(3)
print("(B) MY field: N = JX + I nilpotent? N^3 =", sp.simplify(N**3) == sp.zeros(3), "; eigenvalues of JX:", JX.eigenvals())
Xl = X.subs({x: lam*x, y: lam**2*y, z: z/lam}, simultaneous=True)
print("    equivariance weights (1,2,-1):", [sp.simplify(Xl[i]/lam**w - X[i]) == 0 for i, w in enumerate((1, 2, -1))])
H = sp.Matrix([x + z*d**2, y - d**2, z]); print("    Keller companion H = x + X + x: det JH =", sp.expand(H.jacobian([x, y, z]).det()))
print("    H preserves d = x + yz:", sp.expand((H[0] + H[2]*H[1]) - d) == 0, "; inverse: x = x' - z d^2, y = y' + d^2 with d = x' + z y'")
# invariants of (1,2,-1): a = xz (weight 0), b = y z^2 (weight 0); induced map of H on (a,b)
a, b = sp.symbols('a b')
Ha = sp.expand(H[0]*H[2]); Hb = sp.expand(H[1]*H[2]**2)
def descend2(expr):
    e = sp.expand(expr.subs({x: a/z, y: b/z**2}, simultaneous=True)); assert sp.simplify(sp.diff(e, z)) == 0; return sp.expand(e)
print("    H on the quotient (a,b) = (xz, yz^2):", (sp.factor(descend2(Ha)), sp.factor(descend2(Hb))))
# planar mixed-weight check (THM-1345 hyperbolic case) for weights (a,-b): all weighted-homogeneous Keller pairs up to degree 6 are linear
print("(C) planar hyperbolic weights: search for Keller pairs F = (x^i y^j f(x^b y^a), x^k y^l g(x^b y^a)) with constant Jacobian")
import itertools
for (wa, wb) in [(1, 1), (1, 2), (2, 3), (1, 3)]:
    found = []
    tt = sp.symbols('tt')
    for i, j, k, l in itertools.product(range(0, 4), repeat=4):
        if (i, j) == (0, 0) or (k, l) == (0, 0): continue
        if wa*i - wb*j == wa*k - wb*l: continue   # distinct weights needed for a nondegenerate pair? keep all anyway
        cs = sp.symbols('c0:4'); ds = sp.symbols('d0:4')
        f = sum(cs[m]*(x**wb*y**wa)**m for m in range(3)); g = sum(ds[m]*(x**wb*y**wa)**m for m in range(3))
        P = x**i*y**j*f; Q = x**k*y**l*g
        det = sp.expand(sp.diff(P, x)*sp.diff(Q, y) - sp.diff(P, y)*sp.diff(Q, x))
        poly = sp.Poly(det - 1, x, y)
        eqs = [c for (mx, my), c in poly.terms()]
        sol = sp.solve(eqs, list(cs) + list(ds), dict=True)
        for s in sol:
            Ps, Qs = sp.expand(P.subs(s)), sp.expand(Q.subs(s))
            if sp.Poly(Ps, x, y).total_degree() > 1 or sp.Poly(Qs, x, y).total_degree() > 1:
                found.append((Ps, Qs))
    print(f"    weights ({wa},-{wb}): nonlinear weighted-homogeneous Keller pairs found: {len(found)}")
