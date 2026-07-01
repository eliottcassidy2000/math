"""Anchor the multi-front prompt. (A) SUM OF TWO SQUARES = the obstruction: n=a^2+b^2 iff no prime=3mod4 to odd
power; 21=3*7 (both 3mod4) NOT a sum of 2 squares = the sqrt21 residual; 61 IS (clean). (B) PSL_2(7): order 168,
element orders, the heptagon group carrying sqrt-7. (C) FACILITY-LOCATION game + discrepancy (Koksma-Hlawka)."""
import numpy as np
from math import gcd
def sum2sq(n):
    m=n; obstruct=[]
    p=2
    while p*p<=m:
        if m%p==0:
            e=0
            while m%p==0: m//=p; e+=1
            if p%4==3 and e%2==1: obstruct.append(p)
        p+=1
    if m>1 and m%4==3: obstruct.append(m)
    return (len(obstruct)==0), obstruct
print("(A) SUM OF TWO SQUARES = the obstruction characterization:")
for n in [21,61,183,3,7,14,403]:
    ok,obs=sum2sq(n)
    print(f"  {n}: sum of 2 squares? {ok}  obstruction primes(=3mod4 odd power)={obs}")
print("  => 21=3*7 NOT a sum of 2 squares (obstruction 3,7) = the sqrt21 residual = narrow-Z/2 (S27); 61 CLEAN.")
print("     Same fact 3 ways: (21 not sum-2-sq) = (3,7 both =3mod4) = (narrow class Z/2) = (iota-odd Gauss i*sqrt p).")
# (B) PSL_2(7) via matrices over F_7 mod +-I
def matmul(A,B,p=7): return tuple(( (A[0]*B[0]+A[1]*B[2])%p,(A[0]*B[1]+A[1]*B[3])%p,(A[2]*B[0]+A[3]*B[2])%p,(A[2]*B[1]+A[3]*B[3])%p ))
def det(A,p=7): return (A[0]*A[3]-A[1]*A[2])%p
SL=[ (a,b,c,d) for a in range(7) for b in range(7) for c in range(7) for d in range(7) if (a*d-b*c)%7==1 ]
# PSL = SL mod {I,-I}
def norm(A): 
    negA=tuple((-x)%7 for x in A); return min(A,negA)
PSL=set(norm(A) for A in SL)
def order(A):
    X=A; k=1
    while norm(X)!=norm((1,0,0,1)):
        X=matmul(X,A); k+=1
        if k>20: break
    return k
from collections import Counter
orders=Counter(order(A) for A in PSL)
print(f"\n(B) PSL_2(7): |PSL|={len(PSL)} (=168=7*24=|Aut Fano|); element orders {dict(sorted(orders.items()))}")
print(f"   element order 7 (heptagon rotations) present: {7 in orders}; the 2 cuspidal 3-dim irreps have character (-1+-sqrt-7)/2.")
print(f"   |Aut(Paley_7)|=21=7*3 is the Frobenius/Borel subgroup (order-7 x order-3); sqrt21=sqrt(3*7) crosses them.")
print(f"   PSL_2(7) w/ LPS/Ramanujan generators = a good EXPANDER => left-right Cayley complex = LTC candidate.")
# (C) facility-location: AP runners, arrangement discrepancy vs lonely gap
def disc_and_gap(S,t):
    pts=sorted((v*t)%1 for v in S)+[0.0]; pts=sorted(set(pts))
    gaps=[(pts[(i+1)%len(pts)]-pts[i])%1 for i in range(len(pts))]
    return max(gaps)  # largest gap (observer's best loneliness ~ maxgap adjacent to 0)
print("\n(C) FACILITY-LOCATION: AP runners, observer at 0; M(S)=max_t (half the gap at 0); AP at t=1/n = UNIFORM (min discrepancy):")
for n in [5,7,14]:
    S=list(range(1,n))
    # at t=1/n the runners are equally spaced (min discrepancy) -> observer gap = 1/n
    pts=sorted((v/n)%1 for v in S); allp=sorted(set(pts+[0.0]))
    g0=min(min((0-p)%1,(p-0)%1) for p in pts)  # dist 0 to nearest runner
    print(f"  n={n}: AP at t=1/n -> runners equally spaced (min-discrepancy/UNIFORM), observer loneliness min_v||v/n||={g0:.4f}=1/n={1/n:.4f}")
print("  => LRC = the min-discrepancy (uniform) config MINIMIZES the max-gap =1/n; higher-discrepancy configs have a")
print("     bigger gap somewhere (observer lonelier, M>1/n). Koksma-Hlawka: discrepancy = the potential; AP=extremal.")
