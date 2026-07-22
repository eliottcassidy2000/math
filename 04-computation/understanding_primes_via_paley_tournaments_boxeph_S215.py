#!/usr/bin/env python3
"""understanding_primes_via_paley_tournaments_boxeph_S215.py -- boxeph-2026-07-21-S215

Understand primes 5,7,11 as well as 2,3 -- through the standing lens 'a set of pairwise relations is a
tournament'. EACH PRIME IS ITS PALEY OBJECT (i -> j iff j-i is a quadratic residue mod p):

  p = 3 mod 4  (3, 7, 11)  -> Paley TOURNAMENT : antisymmetric, self-converse, vertex-transitive,
     doubly-regular; eigenvalues (p-1)/2 and (-1 +- i sqrt p)/2 = the QUADRATIC GAUSS SUM i sqrt p.
     This is the 'symmetric-intransitive' HOT pole, opposite the transitive AP nullcone vertex.
  p = 1 mod 4  (5, 13)     -> Paley GRAPH : symmetric, self-COMPLEMENTARY; eigenvalues (-1 +- sqrt p)/2
     = the REAL Gauss sum sqrt p. For p=5, Q(sqrt5) = the GOLDEN field = Fibonacci = the LRC FOIL (S206).
  p = 2                    -> the INVOLUTION itself (Phi_2 = x+1, the antipode / reversal / chirality)
     that acts on all the above (S210-S213).

LRC(14) = 2*7 : its arithmetic is 7 (Paley-7, i sqrt 7, Phi_7 apex), its chirality is 2, its argmax is
3 (Eisenstein Phi_6, t*=14/183), its foil is 5 (golden sqrt5 = Fibonacci), its rank is 11 (S214).
"""
from cmath import exp, pi
from fractions import Fraction as F

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
sep("A  each prime IS its Paley object: tournament (p=3 mod4) vs self-complementary graph (p=1 mod4)")
for p in (3,5,7,11,13):
    A=paley_adj(p)
    antisym=all(A[i][j]+A[j][i]==1 for i in range(p) for j in range(p) if i!=j)   # tournament
    sym=all(A[i][j]==A[j][i] for i in range(p) for j in range(p))                  # graph
    kind = "TOURNAMENT (p=3 mod4)" if antisym else ("GRAPH sym (p=1 mod4)" if sym else "?")
    scores=[sum(A[i]) for i in range(p)]
    reg = len(set(scores))==1
    print(f"  p={p:2d} ({p%4} mod4): {kind}; regular (all out-degree {(p-1)//2})? {reg and scores[0]==(p-1)//2}")
print("  => 3,7,11 (=3 mod4) are TOURNAMENTS; 5,13 (=1 mod4) are self-complementary GRAPHS. 2 is the involution.")

# ==========================================================================
sep("B  the spectrum IS the Gauss sum: Paley char_A = (x-(p-1)/2)(x^2+x+(p+1)/4)^((p-1)/2), roots (-1+-isqrt p)/2")
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
print("  Paley-3: x^2+x+1 = Eisenstein (roots = primitive cube roots, i*sqrt3); Paley-7: (-1+-i*sqrt7)/2 (THM-1830);")
print("  Paley-11: (-1+-i*sqrt11)/2. The imaginary part sqrt(p)/2 IS the quadratic Gauss sum magnitude.")

# ==========================================================================
sep("C  the quadratic GAUSS SUM g_p = sum legendre(a,p) zeta_p^a : |g_p|^2=p ; i*sqrt p (3mod4) vs sqrt p (1mod4)")
for p in (3,5,7,11,13):
    g=sum(legendre(a,p)*exp(2j*pi*a/p) for a in range(1,p))
    realish = abs(g.imag)<1e-9
    print(f"  p={p:2d}: g_p = {g.real:+.4f}{g.imag:+.4f}i ; |g_p|^2 = {abs(g)**2:.4f} = p? {abs(abs(g)**2-p)<1e-6} ; "
          + ("REAL sqrt p (=1 mod4)" if realish else "IMAG i*sqrt p (=3 mod4)"))
print("  => 3,7,11 -> i*sqrt p (the Paley-tournament eigenvalue); 5,13 -> sqrt p (real). p=5: sqrt5 = GOLDEN field.")

# ==========================================================================
sep("D  the LRC(14)=2*7 roles: 7 apex, 2 chirality, 3 argmax (Eisenstein), 5 golden foil, 11 rank")
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
print(f"  7 cap field Q(cos 2pi/7): 2cos(2pi/7)={c7:.5f} root of x^3+x^2-2x-1? {abs(val)<1e-9} (disc 49=7^2, cubic, h=1)")
# 3: argmax t* = 14/Phi_6(14), Phi_6(x)=x^2-x+1
phi6_14=14**2-14+1
print(f"  3 (argmax, Eisenstein): Phi_6(14)=14^2-14+1={phi6_14} => t*=14/{phi6_14} (Paley-3 = the 3-cycle atom)")
# 5: golden field Q(sqrt5), Fibonacci = the LRC foil
phi=(1+5**0.5)/2
print(f"  5 (foil, golden): Q(sqrt5), phi={phi:.5f}, x^2-x-1=0 -> Fibonacci = the LRC FOIL (S206, loosest)")
# 11: rank 11 (S214), and 11 is SCARCE (only multiple <=14 is 11 itself)
mult11=[m for m in range(1,15) if m%11==0]
print(f"  11 (rank): rank-11 AP-core (S214); multiples of 11 in [1,14] = {mult11} (SCARCE => forces the speed 11)")

sep("SUMMARY -- the periodic table of the small primes for LRC(14)/tournaments")
print("""  2  = the INVOLUTION (Phi_2=x+1, the antipode/reversal/chirality) acting on everything (S210-S213).
  3  = the ATOM: Paley-3 = the 3-CYCLE, char x^2+x+1 = Eisenstein (i*sqrt3); the argmax Phi_6 (t*=14/183).
  5  = the GOLDEN/FOIL prime (=1 mod4): Paley GRAPH (self-complementary), Gauss sum REAL sqrt5 = Q(sqrt5)
       = Fibonacci = the LRC foil (loosest, safest -- S206). Also the Bonferroni certificate depth.
  7  = the APEX of LRC(14)=2*7: Paley-7 TOURNAMENT, eigenvalues (-1+-i*sqrt7)/2 = Gauss sum i*sqrt7 (my
       S212 Euler-branch index); Phi_7/Phi_14 carry the hardness; F_7[C_14]=F_7[X]/(X+-1)^7 (THM-2043);
       cap field Q(cos 2pi/7) cubic disc 49.
  11 = the RANK prime (=3 mod4): Paley-11 TOURNAMENT (i*sqrt11); the rank-11 AP-core / relation code (S214);
       SCARCE (only multiple <=14 is 11) => a forced, rigid speed.
  UNIFYING LAW: prime p IS its Paley object -- a TOURNAMENT (self-converse, i*sqrt p) if p=3 mod4 (3,7,11),
  a self-COMPLEMENTARY GRAPH (real sqrt p) if p=1 mod4 (5,13); 2 is the reversal that pairs them. The
  Gauss sum sqrt p / i*sqrt p is literally the Paley spectrum, and it fixes each prime's LRC(14) role.""")
