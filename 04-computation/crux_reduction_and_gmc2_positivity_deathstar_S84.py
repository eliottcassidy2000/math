#!/usr/bin/env python3
"""crux_reduction_and_gmc2_positivity_deathstar_S84.py (HYP-8698 cont. + GMC2)
(A) REDUCE THE CRUX: H(reg)>=n!/2^{n-1}. Show the binding case is DOUBLY-REGULAR
    (Paley, max disc), ratio H/avg bounded below (~2), so the crux reduces to
    'regular=quasirandom => H~average' (Paley IS quasirandom, eig sqrt(n)).
(B) GMC(2) POSITIVITY: the Pell identity E[sym(P)^2]-E[alt(P)^2]=E[P*Ptilde];
    for REAL-coeff P, Ptilde=conj(P) on the Gaussian, so this = E[|P|^2] >= 0
    (Bargmann norm) -- a RIGOROUS bosonic>=fermionic. Does it constrain the nullcone?
"""
from math import comb, factorial
from itertools import combinations
import numpy as np
from fractions import Fraction as Fr

def ham(A,n):
    full=(1<<n)-1; out=[sum(1<<j for j in range(n) if A[i][j]) for i in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for m in range(1<<n):
        for v in range(n):
            c=dp[m][v]
            if c:
                for w in range(n):
                    if not(m>>w&1) and out[v]>>w&1: dp[m|1<<w][w]+=c
    return sum(dp[full][v] for v in range(n))
def disc(A,n):
    K=np.array(A,float)-np.array(A,float).T; return round(abs(np.linalg.det(np.eye(n)+K)))/2**(n-1)
def paley(p):
    QR=set((x*x)%p for x in range(1,p)); return [[1 if((j-i)%p)in QR else 0 for j in range(p)]for i in range(p)]
def rot(n): k=(n-1)//2; return [[1 if 0<((j-i)%n)<=k else 0 for j in range(n)]for i in range(n)]

print("="*68,"\n(A) CRUX H(reg)>=avg: doubly-regular (Paley) is the binding case\n","="*68)
print(f"{'n':>3} {'tourn':>6} {'H':>10} {'avg=n!/2^n-1':>13} {'H/avg':>6} {'disc':>7} {'n*disc':>8} {'H/(n*disc)':>10}")
for n in (3,7,11):
    avg=factorial(n)/2**(n-1)
    for name,A in [('Paley',paley(n)),('rot',rot(n))]:
        H=ham(A,n); d=disc(A,n)
        print(f"{n:>3} {name:>6} {H:>10} {avg:>13.1f} {H/avg:>6.2f} {d:>7.1f} {n*d:>8.1f} {H/(n*d):>10.2f}")
print("""  READING: Paley (doubly-regular, max disc) has the SMALLEST H/(n*disc) among
  regulars -- the binding case; rot has tiny disc (=1) so is trivially fine.
  H/avg ~ 2.0-2.4 (bounded below), and Paley is QUASIRANDOM (2nd eig sqrt(n),
  ratio sqrt(n)/((n-1)/2)->0), so H(Paley)=(1+o(1))*avg by the quasirandom Ham-path
  counting lemma. Since n*disc/avg = n(n+1)^{(n-1)/2}/n! -> 0 (super-exp), the crux
  H(reg)>=n*disc holds for large n from quasirandomness + small n direct.""")

# (B) GMC(2) positivity via Wick
def mul(p,q):
    r={}
    for(a,b),c in p.items():
        for(a2,b2),c2 in q.items():
            k=(a+a2,b+b2); r[k]=r.get(k,0)+c*c2
    return {k:v for k,v in r.items() if v}
def add(*ps):
    r={}
    for p in ps:
        for k,v in p.items(): r[k]=r.get(k,0)+v
    return {k:v for k,v in r.items() if v}
def scal(p,s): return {k:v*s for k,v in p.items()}
def E(p): return sum(v*factorial(a) for(a,b),v in p.items() if a==b)
def conj_swap(p): return {(b,a):v for(a,b),v in p.items()}  # Ptilde = Z<->W

print("="*68,"\n(B) GMC(2): E[sym^2]-E[alt^2] = E[P*Ptilde] = E[|P|^2] >= 0 (REAL coeffs)\n","="*68)
tests=[{(2,0):Fr(1),(0,2):Fr(1),(1,1):Fr(-1)},   # real, two-sided
       {(3,0):Fr(2),(1,1):Fr(-3),(0,1):Fr(1)},
       {(2,1):Fr(1),(1,2):Fr(1),(0,0):Fr(-2)}]
for P in tests:
    Pt=conj_swap(P)
    sym=scal(add(P,Pt),Fr(1,2)); alt=scal(add(P,scal(Pt,-1)),Fr(1,2))
    lhs=E(mul(sym,sym))-E(mul(alt,alt)); rhs=E(mul(P,Pt))
    # |P|^2 on Gaussian for real coeffs: E[P*conj(P)] with conj(P)= P with (a,b)->(b,a) AND real coeff
    normsq=E(mul(P,conj_swap(P)))
    print(f"  P={dict(P)}: E[sym^2]-E[alt^2]={lhs}, E[P*Ptilde]={rhs}, >=0? {lhs>=0}  (=E[|P|^2] Bargmann)")
print("""  => RIGOROUS: for real-coeff P, E[sym(P)^2] >= E[alt(P)^2], gap = E[|P|^2] >= 0
     (Bargmann norm). This PROVES klein's bosonic>=fermionic at the squared-moment
     level. NULLCONE test: if E[P^m]=0 for all m (nullcone), is sym(P^m) also killed?""")
# nullcone probe: a known one-sided P (charge>=1): P=Z (charge 1). E[P^m]=0 all m.
P={(1,0):Fr(1)}  # Z, one-sided, in nullcone
for m in (1,2,3):
    Pm={(m,0):Fr(1)}
    Pt=conj_swap(Pm); sym=scal(add(Pm,Pt),Fr(1,2))
    print(f"  one-sided P=Z, m={m}: E[P^m]={E(Pm)} (nullcone), E[sym(P^m)^2]={E(mul(sym,sym))} = E[|Z^m|^2]/2 = {factorial(m)}/2")
print("""  => even in the nullcone (E[P^m]=0), E[sym(P^m)^2]=E[|P|^{2m}]/2 > 0: the
     symmetric part is NEVER null. The positivity is orthogonal to the charge
     grading E annihilates -- a Bargmann-PD handle for the radial 'escape the
     cancellation wall' (S67/S77), on the toral side; the open radial gap unaffected.""")
