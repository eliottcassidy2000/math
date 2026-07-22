#!/usr/bin/env python3
"""Exact Paley checks plus an explicitly heuristic LRC(14) comparison (S215).

The executable theorem is only the classical ODD-PRIME Paley law for the
relation ``i -> j iff j-i is a nonzero quadratic residue mod p``:

* p == 3 (mod 4): a Paley tournament, with principal adjacency eigenvalue
  (p-1)/2 and nonprincipal eigenvalues (-1 +- i*sqrt(p))/2;
* p == 1 (mod 4): a Paley graph, with principal eigenvalue (p-1)/2 and
  nonprincipal eigenvalues (-1 +- sqrt(p))/2.

Thus the nonprincipal spectrum is a SHIFTED, HALF-SCALED Gauss sum, not the
Gauss sum itself.  The prime 2 lies outside this odd-prime construction.  The
small-prime/LRC labels below are retained as prompts only: they prove no LRC
implication, and coincidences of ranks or fields require an explicit map.
"""
from cmath import exp, pi

def sep(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)
def legendre(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1
def paley_adj(p):   # A[i][j]=1 iff (j-i) is a QR mod p
    return [[1 if legendre(j-i,p)==1 else 0 for j in range(p)] for i in range(p)]
def charpoly_int(A):  # Faddeev-LeVerrier: exact integer char poly, returns coeffs [c_n,...,c_0] of det(xI-A)
    n=len(A); import copy
    M=[[0]*n for _ in range(n)]; I=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    c=[1]+[0]*n; Mk=[row[:] for row in I]
    for k in range(1,n+1):
        AM=[[sum(A[i][t]*Mk[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        ck=-sum(AM[i][i] for i in range(n))//k
        c[k]=ck
        Mk=[[AM[i][j]+(ck if i==j else 0) for j in range(n)] for i in range(n)]
    return c  # det(xI - A) = sum c[k] x^{n-k}

# ==========================================================================
sep("A  exact odd-prime Paley law: tournament (p=3 mod4) vs graph (p=1 mod4)")
for p in (3,5,7,11,13):
    A=paley_adj(p)
    antisym=all(A[i][j]+A[j][i]==1 for i in range(p) for j in range(p) if i!=j)   # tournament
    sym=all(A[i][j]==A[j][i] for i in range(p) for j in range(p))                  # graph
    kind = "TOURNAMENT (p=3 mod4)" if antisym else ("GRAPH sym (p=1 mod4)" if sym else "?")
    scores=[sum(A[i]) for i in range(p)]
    reg = len(set(scores))==1
    print(f"  p={p:2d} ({p%4} mod4): {kind}; regular (all out-degree {(p-1)//2})? {reg and scores[0]==(p-1)//2}")
print("  => 3,7,11 (=3 mod4) are TOURNAMENTS; 5,13 (=1 mod4) are self-complementary GRAPHS.")
print("     p=2 is outside this odd-prime Paley law; no 'i*sqrt(2)' claim is made.")

# ==========================================================================
sep("B  tournament spectrum: a shifted half-scale of the Gauss sum, plus the principal eigenvalue")
for p in (3,7,11):
    A=paley_adj(p); c=charpoly_int(A)  # det(xI-A)
    # expected quadratic factor x^2 + x + (p+1)/4
    q=(p+1)//4
    # verify char poly == (x-(p-1)/2)*(x^2+x+q)^((p-1)/2) by polynomial multiply
    def pmul(a,b):
        r=[0]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b): r[i+j]+=x*y
        return r
    expected=[1, -(p-1)//2]  # (x - (p-1)/2)  as [1, -(p-1)/2]
    quad=[1,1,q]
    for _ in range((p-1)//2): expected=pmul(expected,quad)
    # c is det(xI-A) high-to-low; expected is high-to-low too
    match = c==expected
    disc = 1-4*q  # discriminant of x^2+x+q  = 1-(p+1) = -p
    print(f"  p={p:2d}: char_A = (x-{(p-1)//2})(x^2+x+{q})^{(p-1)//2} ? {match}; quad disc = {disc} = -p ; roots (-1+-i*sqrt{p})/2")
print("  Paley-3: x^2+x+1 = Phi_3 (primitive cube roots); Paley-7: (-1+-i*sqrt7)/2 (THM-1830);")
print("  Paley-11: (-1+-i*sqrt11)/2. Nonprincipal eigenvalues are (-1 +- g_p)/2, not g_p itself.")

# ==========================================================================
sep("C  the quadratic GAUSS SUM g_p = sum legendre(a,p) zeta_p^a : |g_p|^2=p ; i*sqrt p (3mod4) vs sqrt p (1mod4)")
for p in (3,5,7,11,13):
    g=sum(legendre(a,p)*exp(2j*pi*a/p) for a in range(1,p))
    realish = abs(g.imag)<1e-9
    print(f"  p={p:2d}: g_p = {g.real:+.4f}{g.imag:+.4f}i ; |g_p|^2 = {abs(g)**2:.4f} = p? {abs(abs(g)**2-p)<1e-6} ; "
          + ("REAL sqrt p (=1 mod4)" if realish else "IMAG i*sqrt p (=3 mod4)"))
print("  => 3,7,11 have g_p=i*sqrt(p); 5,13 have g_p=sqrt(p). Adjacency eigenvalues are shifted half-scales.")
print("     For Paley graphs the separate principal adjacency eigenvalue is (p-1)/2.")

# ==========================================================================
sep("D  heuristic small-prime comparisons for LRC(14): arithmetic checks, no implication")
# 7: x^14 - 1 = Phi_1 Phi_2 Phi_7 Phi_14 ; mod 7 Frobenius collapse (x^14-1 == (x-1)^7 (x+1)^7 mod 7)
def poly_mod(coeffs,m): return [c%m for c in coeffs]
# (x-1)^7 (x+1)^7 mod 7 vs x^14 - 1
def pmul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
xm1_7=[1];
for _ in range(7): xm1_7=pmul(xm1_7,[1,-1])   # (x-1)^7
xp1_7=[1]
for _ in range(7): xp1_7=pmul(xp1_7,[1,1])    # (x+1)^7
prod=pmul(xm1_7,xp1_7)
x14m1=[1]+[0]*13+[-1]
print("  7 (apex, 14=2*7): x^14-1 == (x-1)^7 (x+1)^7 (mod 7)?", poly_mod(prod,7)==poly_mod(x14m1,7),
      " [THM-2043 F_7[C_14]=F_7[X]/(X-1)^7 x F_7[X]/(X+1)^7]")
# Q(cos 2pi/7): 2cos(2pi/7) has min poly x^3 + x^2 - 2x - 1, discriminant 49 = 7^2
import cmath
c7=2*cmath.cos(2*pi/7).real
val=c7**3+c7**2-2*c7-1
print(f"  7 cap field Q(cos 2pi/7): 2cos(2pi/7)={c7:.5f} root of x^3+x^2-2x-1? {abs(val)<1e-9} (disc 49=7^2, cubic)")
# 3: argmax t* = 14/Phi_6(14), Phi_6(x)=x^2-x+1
phi6_14=14**2-14+1
print(f"  3 comparison: Phi_6(14)=14^2-14+1={phi6_14} => t*=14/{phi6_14}; Paley-3 has Phi_3, with Phi_6(x)=Phi_3(-x)")
# 5: golden field Q(sqrt5), Fibonacci = the LRC foil
phi=(1+5**0.5)/2
print(f"  5 comparison: Q(sqrt5), phi={phi:.5f}, x^2-x-1=0 -> Fibonacci/golden-field analogy (S206)")
# 11: a numerical rank/scarcity comparison only; no Paley-to-LRC map is supplied.
mult11=[m for m in range(1,15) if m%11==0]
print(f"  11 comparison: rank ker_Z(1,...,12)=11; multiples of 11 in [1,14] = {mult11}; neither fact forces an LRC speed")

sep("SUMMARY -- exact Paley kernel and limits of the LRC analogy")
print("""  EXACT: for odd prime p, the quadratic-residue relation is a tournament for p=3 mod4 and a graph
  for p=1 mod4. Its nonprincipal adjacency spectrum is (-1 +- g_p)/2, while (p-1)/2 is principal.
  OUTSIDE SCOPE: p=2 is not covered. Phi_3 (Paley-3) and Phi_6 (the cited LRC parameter) differ by x -> -x.
  HEURISTIC ONLY: the labels 2/chirality, 3/argmax, 5/foil, 7/apex, and 11/rank are comparison prompts.
  No Paley object, field match, scarcity observation, or equal rank by itself implies an LRC(14) fact.""")
