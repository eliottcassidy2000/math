"""C1: local intersection multiplicities (Serre chi = Tor Euler characteristic)
at THM-1300's triple collision, and the fibre-degree bookkeeping.

F = (u^3 z + y^2 u (4+3xy), y + 3x u^2 z + 3x y^2 (4+3xy), 2x - 3x^2 y - x^3 z), u = 1+xy.

For a target point p, the fibre scheme is Z_p = Spec Q[x,y,z]/I_p, I_p = (F - p).
At an isolated point P of Z_p the three equations form a regular sequence in the
3-dimensional regular local ring O_P (height 3 = dim), so the Koszul complex is a
free resolution of O_P/I_p and Tor_i(O_P/(F1-p1), O_P/(F2-p2,F3-p3)) = 0 for i > 0.
Hence Serre's chi at P equals the plain length dim_Q (O_P / I_p).

We compute: a lex/grevlex Groebner basis of I_p, the normal set (standard monomials),
the multiplication matrices M_x, M_y, M_z on Q[x,y,z]/I_p, and the generalized
eigenspace dimensions of a separating linear form -> local multiplicities.
"""
import sympy as sp
from sympy import Rational as R

x, y, z = sp.symbols('x y z')
u = 1 + x*y
F1 = sp.expand(u**3*z + y**2*u*(4 + 3*x*y))
F2 = sp.expand(y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y))
F3 = sp.expand(2*x - 3*x**2*y - x**3*z)
J = sp.Matrix([[sp.diff(f, v) for v in (x, y, z)] for f in (F1, F2, F3)])
detJ = sp.expand(J.det())
print("det JF =", detJ)
assert detJ == -2

def fibre_algebra(p, verbose=True):
    I = [F1 - p[0], F2 - p[1], F3 - p[2]]
    G = sp.groebner(I, x, y, z, order='grevlex')
    # normal set: monomials not divisible by any leading monomial
    lms = [sp.Poly(g, x, y, z).monoms(order='grevlex')[0] for g in G.exprs]
    # enumerate standard monomials up to a degree bound
    std = []
    B = 12
    for a in range(B):
        for b in range(B):
            for c in range(B):
                if not any(a >= l[0] and b >= l[1] and c >= l[2] for l in lms):
                    std.append((a, b, c))
    if verbose:
        print("  target", p, " leading monomials", lms, " dim =", len(std))
    return G, std

def mult_matrix(G, std, var):
    n = len(std)
    M = sp.zeros(n, n)
    for j, (a, b, c) in enumerate(std):
        mon = x**a*y**b*z**c*var
        r = G.reduce(sp.expand(mon))[1]
        P = sp.Poly(r, x, y, z)
        for mono, coeff in zip(P.monoms(), P.coeffs()):
            i = std.index(mono)
            M[i, j] = coeff
    return M

print("\n[1] fibre over (-1/4, 0, 0)")
p = (R(-1, 4), 0, 0)
G, std = fibre_algebra(p)
Mx, My, Mz = (mult_matrix(G, std, v) for v in (x, y, z))
# separating form l = x + 2y + 3z
L = Mx + 2*My + 3*Mz
lam = sp.symbols('lam')
cp = sp.factor(L.charpoly(lam).as_expr())
print("  char poly of x+2y+3z:", cp)
pts = [(0, 0, R(-1, 4)), (1, R(-3, 2), R(13, 2)), (-1, R(3, 2), R(13, 2))]
for P in pts:
    assert all(sp.simplify(f.subs({x: P[0], y: P[1], z: P[2]}) - q) == 0 for f, q in zip((F1, F2, F3), p))
    val = P[0] + 2*P[1] + 3*P[2]
    # algebraic multiplicity of eigenvalue val = local length dim O_P/I
    n = L.shape[0]
    mult = n - ((L - val*sp.eye(n))**n).rank()
    print("  point", P, " Jacobian det", detJ, " local length (= Serre chi) =", mult)

print("\n[2] generic fibre degree: random targets")
for p in [(R(2), R(3), R(5)), (R(-7, 3), R(1, 2), R(11)), (R(1), R(0), R(0))]:
    G, std = fibre_algebra(p)

print("\n[3] the origin and the a-axis near 0: orbit branch escapes to infinity")
for p in [(0, 0, 0), (R(1, 1000), 0, 0)]:
    G, std = fibre_algebra(p)
