import sympy as sp
x, y, z = sp.symbols("x y z")
F1 = (1 + x*y)**3 * z + y**2 * (1 + x*y) * (4 + 3*x*y)
F2 = y + 3*x*(1 + x*y)**2 * z + 3*x*y**2 * (4 + 3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
F = sp.Matrix([F1, F2, F3])
J = F.jacobian([x, y, z])
print("det J =", sp.factor(sp.expand(J.det())))
pts = [(0, 0, sp.Rational(-1, 4)), (1, sp.Rational(-3, 2), sp.Rational(13, 2)), (-1, sp.Rational(3, 2), sp.Rational(13, 2))]
for p in pts:
    print(p, "->", [sp.nsimplify(c.subs({x: p[0], y: p[1], z: p[2]})) for c in F])
print("degrees:", [sp.Poly(sp.expand(c), x, y, z).total_degree() for c in F])
