"""Normalized Apollonius: circle 1 at the origin with radius 1 (translation + scaling), circles 2, 3 free: 6 parameters.
Factor the four discriminants and their product; identify the square classes."""
import sympy as sp
a2,b2,r2,a3,b3,r3,x,y,r = sp.symbols('a2 b2 r2 a3 b3 r3 x y r')
C = [(0,0,1),(a2,b2,r2),(a3,b3,r3)]
def disc(e):
    eqs = [sp.expand((x-a)**2 + (y-b)**2 - (r + ei*rr)**2) for (a,b,rr), ei in zip(C, e)]
    lin = [sp.expand(eqs[0]-eqs[1]), sp.expand(eqs[0]-eqs[2])]
    sol = sp.solve(lin, [x, y], dict=True)[0]
    q = sp.together(sp.expand(eqs[0].subs(sol))); num, den = sp.fraction(q)
    P = sp.Poly(sp.expand(num), r); A, B, Cc = P.all_coeffs()
    return sp.expand(B**2 - 4*A*Cc), sp.factor(den)
classes = [(1,1,1),(1,1,-1),(1,-1,1),(-1,1,1)]
D = {}
for e in classes:
    d, den = disc(e); D[e] = d
    fl = sp.factor_list(d)
    print("class", e, "denominator", den); print("   factors:", [(sp.factor(f), m) for f, m in fl[1]], "const", fl[0]); print()
prod = sp.factor_list(sp.expand(D[classes[0]]*D[classes[1]]*D[classes[2]]*D[classes[3]]))
print("product factor multiplicities:", [(str(f)[:80], m) for f, m in prod[1]], "const", prod[0])
