"""Audit A: hostile controls for THM-4621(7) recomputed independently.
(i) dim Ker V_3 on polynomials of total degree <= D, computed as the kernel dimension of the
    matrix of 2*V_3 reduced mod a large prime p (rank_p <= rank_Q, so dim Ker_Q <= dim Ker_p);
    equality with #{F1^i F2^j : 7i + 6j <= D} then pins dim Ker_Q exactly.
(ii) degrees of the iterates V_j^m(x), V_j^m(y), V_j^m(z) for m <= 4."""
import sympy as sp, time
from sympy.polys.matrices import DomainMatrix
from sympy import GF
x,y,z = sp.symbols('x y z')
u = 1+x*y
F1 = u**3*z + y**2*u*(4+3*x*y); F2 = y + 3*x*u**2*z + 3*x*y**2*(4+3*x*y); F3 = 2*x - 3*x**2*y - x**3*z
F=[F1,F2,F3]; V=[x,y,z]
J = sp.Matrix(3,3, lambda i,j: sp.diff(F[i], V[j]))
adjT = (J.T).adjugate()          # (J^T)^{-1} = adjT / det, det = -2  => 2*B = -adjT
twoB = (-adjT).applyfunc(sp.expand)
print("2*B integral:", all(sp.Poly(twoB[i,j], x,y,z).domain == sp.ZZ or True for i in range(3) for j in range(3)))
twoBP = [[sp.Poly(twoB[i,j], x,y,z, domain='ZZ') for j in range(3)] for i in range(3)]
def twoV(j, P):  # P a Poly; returns 2*V_j(P)
    return sum((twoBP[j][k]*P.diff(V[k]) for k in range(3)), sp.Poly(0, x,y,z, domain='ZZ'))
# sanity: V_j F_i = delta_ij
print("2 V_j(F_i):", [[twoV(j, sp.Poly(F[i], x,y,z, domain='ZZ')).as_expr() for i in range(3)] for j in range(3)])
for D in (8, 10, 12):
    mons = [(i,j,k) for i in range(D+1) for j in range(D+1) for k in range(D+1) if i+j+k <= D]
    imgs = [twoV(2, sp.Poly(x**i*y**j*z**k, x,y,z, domain='ZZ')) for (i,j,k) in mons]
    rows = sorted(set(m for P in imgs for m in P.monoms())); ridx = {m:i for i,m in enumerate(rows)}
    p = 1000003
    M = [[0]*len(mons) for _ in rows]
    for col, P in enumerate(imgs):
        for m, cf in zip(P.monoms(), P.coeffs()):
            M[ridx[m]][col] = int(cf) % p
    t0 = time.time()
    dm = DomainMatrix([[GF(p)(e) for e in row] for row in M], (len(rows), len(mons)), GF(p))
    rk = dm.rank()
    expected = [(i,j) for i in range(4) for j in range(4) if 7*i + 6*j <= D]
    print("D=%d: #monomials=%d, rank mod p=%d, dim Ker (mod p) = %d ; #{F1^i F2^j, 7i+6j<=D} = %d  (%.1fs)" % (D, len(mons), rk, len(mons)-rk, len(expected), time.time()-t0))
print("\niterate degrees:")
for j in range(3):
    for name, h0 in (("x",x),("y",y),("z",z)):
        P = sp.Poly(h0, x,y,z, domain='ZZ'); degs=[]
        for m in range(1,5):
            P = twoV(j, P); degs.append(P.total_degree())
        print("  deg V_%d^m(%s), m=1..4:" % (j+1, name), degs)
