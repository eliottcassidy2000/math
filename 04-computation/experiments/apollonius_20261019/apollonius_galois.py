"""Galois structure of the problem of Apollonius.  Three circles C_i: centre (a_i, b_i), radius r_i.  A circle (x, y, r) is
tangent to C_i iff (x-a_i)^2 + (y-b_i)^2 = (r + e_i r_i)^2, e_i = +-1.  Global sign flip (e, r) -> (-e, -r) identifies
sign patterns in pairs, leaving 4 classes; each class gives two linear equations (differences) and one quadratic in r:
two solutions conjugate over Q(sqrt Delta_e).  We compute the four discriminants symbolically and test multiplicative
relations (is the product of all four a square? any pairwise product?)."""
import sympy as sp
a1,b1,r1,a2,b2,r2,a3,b3,r3,x,y,r = sp.symbols('a1 b1 r1 a2 b2 r2 a3 b3 r3 x y r')
C = [(a1,b1,r1),(a2,b2,r2),(a3,b3,r3)]
def disc(e):
    eqs = [sp.expand((x-a)**2 + (y-b)**2 - (r + ei*rr)**2) for (a,b,rr), ei in zip(C, e)]
    lin = [sp.expand(eqs[0]-eqs[1]), sp.expand(eqs[0]-eqs[2])]      # linear in x, y, r
    sol = sp.solve(lin, [x, y], dict=True)[0]
    q = sp.expand(eqs[0].subs(sol))
    q = sp.together(q); num, den = sp.fraction(q)
    P = sp.Poly(sp.expand(num), r)
    assert P.degree() == 2
    A, B, Cc = P.all_coeffs()
    return sp.factor(sp.expand(B**2 - 4*A*Cc)), sp.factor(A)
classes = [(1,1,1),(1,1,-1),(1,-1,1),(-1,1,1)]
D = {}; L = {}
for e in classes:
    D[e], L[e] = disc(e)
    print("sign class", e, ": leading coefficient", L[e]); print("   discriminant =", D[e]); print()
def is_square(expr):
    expr = sp.factor(expr)
    if expr == 0: return True
    c, fl = sp.factor_list(expr)
    return all(m % 2 == 0 for _, m in fl) and sp.sqrt(c).is_rational
prod_all = sp.factor(sp.Mul(*[D[e] for e in classes]))
print("product of the four discriminants is a perfect square in Q(a,b,r):", is_square(prod_all))
import itertools
for e, f in itertools.combinations(classes, 2):
    print("pairwise product square?", e, f, is_square(D[e]*D[f]))
for trip in itertools.combinations(classes, 3):
    print("triple product square?", trip, is_square(sp.Mul(*[D[e] for e in trip])))
# the common factor structure: factor each discriminant
for e in classes:
    print("factor_list", e, sp.factor_list(D[e])[1])
