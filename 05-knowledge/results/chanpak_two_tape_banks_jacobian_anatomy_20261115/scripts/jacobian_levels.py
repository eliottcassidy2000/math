# The dimension-3 Keller map of THM-1300 as (fold family) o (inversion), its level surfaces, and a planar pole count.
import sympy as sp
x, y, z, s, w, c = sp.symbols('x y z s w c')
u = 1 + x*y; l = 2 - 3*x*y - x**2*z
F = [u**3*z + y**2*u*(4 + 3*x*y), y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y), 2*x - 3*x**2*y - x**3*z]
print("det JF =", sp.factor(sp.Matrix(F).jacobian([x, y, z]).det()), "; F3 = x*l:", sp.expand(F[2] - x*l) == 0)
# (1) F = Psi o tau
Psi = [s**2 + s*w - c*s**3, 4*s + 2*w - 3*c*s**2, c]
tau = {s: u/x, w: 1/x, c: x*l}
print("F = Psi o tau:", all(sp.simplify(Psi[i].subs(tau) - F[i]) == 0 for i in range(3)))
print("det J Psi =", sp.factor(sp.Matrix(Psi).jacobian([s, w, c]).det()), "; det J tau =", sp.factor(sp.simplify(sp.Matrix([u/x, 1/x, x*l]).jacobian([x, y, z]).det())))
# tau is an isomorphism {x != 0} -> {w != 0}: inverse
inv = {x: 1/w, y: s - w}
zi = sp.solve(sp.Eq((x*l).subs(inv), c), z)[0]
print("inverse of tau: x = 1/w, y = s - w, z =", sp.expand(zi))
# (2) Psi is the root cover of the cubic family  c S^3 - 2 S^2 + b S - 2a
a, b = Psi[0], Psi[1]
print("c s^3 - 2 s^2 + b s - 2a on the image:", sp.expand(c*s**3 - 2*s**2 + b*s - 2*a))
# general Legendre family: Psi_r = (s w + r - s r', 2w - 2r', c) has Jacobian 2w for every r(s, c)
r = sp.Function('r')(s, c)
Pr = [s*w + r - s*sp.diff(r, s), 2*w - 2*sp.diff(r, s), c]
print("Legendre family Jacobian:", sp.simplify(sp.Matrix(Pr).jacobian([s, w, c]).det()))
# (3) level surfaces: on F3 = c the pair (F1, F2) in the plane coordinates (s, w), w != 0, is (Psi_1, Psi_2): Jacobian 2w;
#     on x = 0 the map is (z + 4y^2, y, 0); on l = 0 it is ((1+v)(2+v)/x^2, (6+4v)/x, 0) with v = xy
print("on x = 0:", [sp.expand(f.subs(x, 0)) for f in F])
v = sp.symbols('v')
zl = sp.solve(sp.Eq(l, 0), z)[0]
print("on l = 0:", [sp.factor(sp.simplify(f.subs(z, zl).subs(y, v/x))) for f in F])
# (4) planar pole count for the same ansatz with constant coefficients: psi_r o (g x^3 y + phi(x), 1/x), r quadratic
g_, rho, rho0, rho1, p0 = sp.symbols('gamma rho rho0 rho1 p0')
rq = rho*s**2/2 + rho0*s + rho1
psi = [s*w + rq - s*sp.diff(rq, s), 2*w - 2*sp.diff(rq, s)]
sub = {s: g_*x**3*y + 1/(rho*x) + p0, w: 1/x}
G = [sp.expand(sp.simplify(f.subs(sub))) for f in psi]
print("planar ansatz, quadratic r: second component polynomial:", sp.Poly(sp.together(G[1]).as_numer_denom()[1], x).degree() == 0, "; first component =", sp.collect(G[0], x))
