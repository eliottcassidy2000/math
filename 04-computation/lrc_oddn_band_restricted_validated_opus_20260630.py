"""Validate search idea A (band-restricted enumeration): iterate (q,a,m), band = allowed small speeds, ENUMERATE
covering sets with residues in band, find min M. Should find n=7 covering-min 2/13 fast."""
import math
from fractions import Fraction
from itertools import combinations
def M_exact(S,Qmax):
    best=Fraction(0)
    for q in range(2,Qmax+1):
        bb=0
        for a in range(1,q):
            if math.gcd(a,q)!=1: continue
            m=q
            for s in S:
                r=(s*a)%q; d=r if r<=q-r else q-r
                if d<m:m=d
                if m<=bb:break
            if m>bb:bb=m
        v=Fraction(bb,q)
        if v>best:best=v
    return best
def is_cov(S,n): return all(any(x%q==0 for x in S) for q in range(2,n+1))
def band_search(n,qrange,smax_mult=2):
    Qmax=2*(n*n-n+1)
    best=Fraction(1); bestS=None
    seen=set()
    for q in qrange:
        for a in range(1,q):
            if math.gcd(a,q)!=1: continue
            for m in range(2, q//(n-1)+2):
                # allowed small speeds: residue s*a mod q in [m, q-m]
                allowed=[s for s in range(1, smax_mult*q) if m<=(s*a)%q<=q-m]
                if len(allowed)<n-1: continue
                # restrict to small allowed for enumeration
                small=[s for s in allowed if s<=2*n+2][:14]
                if len(small)<n-1: continue
                for combo in combinations(small,n-1):
                    if combo in seen: continue
                    seen.add(combo)
                    S=list(combo)
                    if not is_cov(S,n): continue
                    M=M_exact(S,Qmax)
                    if M<best: best=M; bestS=S
    return best,bestS
M,S=band_search(7,range(8,20))
print(f"n=7 band-restricted search: covering-min={M}={float(M):.5f} at {S}  (target 2/13={float(Fraction(2,13)):.5f})")
print(f"   found 2/13: {M==Fraction(2,13)}  => search idea A WORKS (band restricts residues, enumeration collapses).")
print(f"   (the band {{s: m<=sa mod q<=q-m}} = the independent set in the danger circulant C(q;+-1..+-(m-1)).)")
