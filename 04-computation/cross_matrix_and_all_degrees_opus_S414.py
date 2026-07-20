"""opus-2026-07-20-S414 -- HYP-8390 RESOLVED, AND THE ONE-SIDED CONJECTURE AT ALL DEGREES.

THE LEVER.  For charge span {-1,0,+1}, P = z*A(s) + B(s) + zbar*C(s) with s = z*zbar.
The cross term of E[P^2] is 2*E[ z*A * zbar*C ] = 2*E[ s*A(s)*C(s) ], so with
A = sum a_i s^i and C = sum c_j s^j,

      cross = 2 * A^T M C ,      M_{ij} = E[s^{i+j+1}] = (i+j+1)!

M is the HANKEL MOMENT MATRIX of the measure s*e^{-s} ds -- the GAMMA(2) distribution --
because int_0^inf s^k * s e^{-s} ds = (k+1)!.  A Hankel matrix of a positive measure with
infinite support is POSITIVE DEFINITE, hence nonsingular, at every size.  HYP-8390: YES.

BUT NONSINGULARITY ALONE IS NOT ENOUGH: A^T M C = 0 with M nonsingular does NOT force A=0 or
C=0 (take M = I, A = (1,0), C = (0,1)).  The second moment cannot finish the job.  So go to
ALL moments, where the structure is much better:

THE STRUCTURAL REDUCTION.  A product of m factors has total charge 0 iff it uses equally
many z*A and zbar*C factors, say p of each, and r = m-2p copies of B.  The product is then
(z A)^p (zbar C)^p B^r = s^p A^p C^p B^r = h^p B^r  with  h := s*A*C.  Hence

      E[P^m]  =  sum_{2p+r=m}  m!/(p! p! r!)  E[ h^p B^r ]

-- the nullcone condition depends ONLY on the pair (h, B), not on A and C separately.

CONSEQUENCE (all degrees, no elimination).  If B = 0 the sum collapses to the even terms:
E[P^{2p}] = (2p)!/(p!p!) * E[h^p], so the nullcone forces E[h^p] = 0 for every p >= 1.  By
THM-1570 sB (the Laplace layer: the only g with all int g^p e^{-s} ds = 0 is g = 0) we get
h = 0, i.e. s*A*C = 0, and since C[s] is a DOMAIN, A = 0 or C = 0.  ONE-SIDED, AT EVERY
DEGREE.  This replaces the degree<=1 Groebner result of THM-1570 sA with a structural proof.
"""
import sympy as sp
from math import factorial

z, zb, t, s = sp.symbols('z zb t s')

def E1(e):
    e = sp.expand(e)
    if e == 0: return sp.Integer(0)
    p = sp.Poly(e, z, zb); tot = 0
    for (a, b), c in zip(p.monoms(), p.coeffs()):
        if a == b: tot += c*factorial(a)
    return sp.expand(tot)

def Es(g):
    """E[g(s)] with s ~ Exp(1): E[s^k] = k!"""
    g = sp.expand(g)
    if g == 0: return sp.Integer(0)
    p = sp.Poly(g, s)
    return sp.expand(sum(c*factorial(k) for (k,), c in zip(p.monoms(), p.coeffs())))

print("="*78)
print("(1) HYP-8390 RESOLVED:  M_{ij} = (i+j+1)!  is the GAMMA(2) Hankel matrix")
print("="*78)
import numpy as np
for N in range(1, 8):
    M = np.array([[float(factorial(i+j+1)) for j in range(N)] for i in range(N)])
    ev = np.linalg.eigvalsh(M)
    print(f"   N={N}: M = (i+j+1)!   min eigenvalue = {ev.min():.6g}   pos.def: {bool(ev.min()>0)}")
print("   at N=2: M = [[1,2],[2,6]], and 2*M = [[2,4],[4,12]] -- MATCHES THM-1570's")
print("   observed cross term 2a0c0+4a0c1+4a1c0+12a1c1 exactly.")
print("   => M is nonsingular at every degree.  BUT A^T M C = 0 does NOT force A=0 or C=0")
print("      (M=I, A=(1,0), C=(0,1)), so the SECOND moment alone cannot finish.")

print()
print("="*78)
print("(2) THE STRUCTURAL REDUCTION:  E[P^m] = sum_{2p+r=m} m!/(p!p!r!) E[h^p B^r]")
print("="*78)
D = 2
A = sp.symbols(f'a0:{D+1}'); B = sp.symbols(f'b0:{D+1}'); C = sp.symbols(f'c0:{D+1}')
ss = z*zb
Ap = sum(A[k]*ss**k for k in range(D+1))
Bp = sum(B[k]*ss**k for k in range(D+1))
Cp = sum(C[k]*ss**k for k in range(D+1))
P = sp.expand(z*Ap + Bp + zb*Cp)
As_ = sum(A[k]*s**k for k in range(D+1))
Bs_ = sum(B[k]*s**k for k in range(D+1))
Cs_ = sum(C[k]*s**k for k in range(D+1))
h = sp.expand(s*As_*Cs_)
Pm = sp.Integer(1); ok = True
for m in range(1, 7):
    Pm = sp.expand(Pm*P)
    lhs = E1(Pm)
    rhs = 0
    for p in range(0, m//2 + 1):
        r = m - 2*p
        rhs += sp.Rational(factorial(m), factorial(p)*factorial(p)*factorial(r)) * \
               Es(sp.expand(h**p * Bs_**r))
    good = sp.simplify(sp.expand(lhs - rhs)) == 0
    ok = ok and good
    print(f"   m={m}: identity holds: {good}")
print(f"   => structural reduction VERIFIED through m=6 at degree {D}: {ok}")

print()
print("="*78)
print("(3) THE ALL-DEGREE THEOREM (case B = 0)")
print("="*78)
print("   With B = 0 the sum keeps only r=0, so E[P^{2p}] = (2p)!/(p!p!) * E[h^p].")
print("   Nullcone => E[h^p] = 0 for all p >= 1 => (THM-1570 sB) h = 0 => s*A*C = 0")
print("   => A*C = 0 in the DOMAIN C[s] => A = 0 or C = 0.   ONE-SIDED, ALL DEGREES.")
print()
print("   check the collapse symbolically at B=0:")
sub0 = {B[k]: 0 for k in range(D+1)}
Pm = sp.Integer(1)
for m in range(1, 7):
    Pm = sp.expand(Pm*P)
    lhs = sp.expand(E1(Pm).subs(sub0))
    if m % 2 == 1:
        print(f"   m={m} (odd) : E[P^m] = {lhs}   (should be 0)")
    else:
        p = m//2
        rhs = sp.Rational(factorial(m), factorial(p)**2)*Es(sp.expand(h**p))
        print(f"   m={m} (even): E[P^m] - C(m,p)E[h^p] = {sp.simplify(sp.expand(lhs-rhs))}")
