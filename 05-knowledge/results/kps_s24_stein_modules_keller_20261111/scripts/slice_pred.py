"""slice_pred.py - derive the PREDICTION for the Stein dimension of a planar slice of the C^3 Keller map (THM-1300).

F = (F1, F2, F3); master identity (THM-1335): L(F) x^3 + (4 - 3 F2 F3) x - 2 F3 = 0 with L a quartic in the target.
Slice F1 = c: X_c = {1 + x y != 0} (z solved), (u, v) = (F2, F3) restricted; Jacobian -2/(1+xy)^3.
Non-properness curve of the slice: S_c = {L(c, u, v) = 0}; generic fibre over S_c has 1 point instead of 3."""
import sympy as sp

x, y, z, a, b, g, c = sp.symbols("x y z a b g c")
F1 = (1 + x*y)**3 * z + y**2 * (1 + x*y) * (4 + 3*x*y)
F2 = y + 3*x*(1 + x*y)**2 * z + 3*x*y**2 * (4 + 3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
# L(F) x^3 = 2 F3 - (4 - 3 F2 F3) x
num = sp.expand(2*F3 - (4 - 3*F2*F3)*x)
Lx = sp.cancel(num / x**3)
assert sp.denom(Lx) == 1
Lx = sp.expand(Lx)
# express Lx as a polynomial in (F1, F2, F3) of degree <= 4: linear system on coefficients
mons = [(i, j, k) for i in range(5) for j in range(5) for k in range(5) if i + j + k <= 4]
cs = sp.symbols("c0:%d" % len(mons))
expr = sum(ci * F1**i * F2**j * F3**k for ci, (i, j, k) in zip(cs, mons)) - Lx
eqs = sp.Poly(sp.expand(expr), x, y, z).coeffs()
sol = sp.solve(eqs, cs, dict=True)[0]
L = sum(sol.get(ci, 0) * a**i * b**j * g**k for ci, (i, j, k) in zip(cs, mons))
L = sp.factor(sp.expand(L))
print("L(a, b, g) =", L)
Lc = sp.expand(L.subs({a: c}))
u, v = sp.symbols("u v")
Lc = Lc.subs({b: u, g: v})
P = sp.Poly(Lc, u, v)
print("slice curve L(c, u, v) =", sp.collect(sp.expand(Lc), [u, v]))
print("total degree", P.total_degree(), "deg_u", sp.degree(Lc, u), "deg_v", sp.degree(Lc, v))
L1 = sp.Poly(Lc.subs(c, 1), u, v)
print("c = 1: irreducible over Q:", len(sp.factor_list(L1.as_expr())[1]) == 1 and sp.factor_list(L1.as_expr())[1][0][1] == 1)
# slice map
w = sp.symbols("w")
zz = sp.solve(sp.Eq(F1, c), z)[0]
uu = sp.simplify(F2.subs(z, zz))
vv = sp.simplify(F3.subs(z, zz))
print("u =", sp.factor(uu), "  v =", sp.factor(vv))
J = sp.simplify(sp.diff(uu, x)*sp.diff(vv, y) - sp.diff(uu, y)*sp.diff(vv, x))
print("Jacobian of (u, v) in (x, y):", sp.factor(J))
# intersections of S_1 with a generic hyperbola (u - 2)(v - 3) = t: eliminate v
t = sp.symbols("t")
vsol = 3 + t/(u - 2)
res = sp.together(L1.as_expr().subs(v, vsol))
numr = sp.Poly(sp.numer(res), u)
print("number of affine intersections of S_1 with (u-2)(v-3) = t (degree in u of the numerator):", numr.degree(), "; squarefree for t = 5:", sp.Poly(sp.numer(sp.together(L1.as_expr().subs(v, 3 + 5/(u - 2)))), u).discriminant() != 0)
