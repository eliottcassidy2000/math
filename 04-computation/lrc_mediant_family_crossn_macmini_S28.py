#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S28 (HYP-4582) -- the MEDIANT construction across N.

opus-S118: the mediant family is S_N = {1,..,N-2} u {N, 3(N-1)} (binding pair
{5,3(N-1)} sums to q=3N+2).  TEST: M(S_N) achieves the mediant 3/(3N+2) (in gap)
only at N=7 (q=23) and N=13 (q=41) -- both PRIME and N==1 mod 6.  At N=12 it gives
M=3/35 (q=35=5*7, a BETTER WITNESS intrudes) => the mediant 3/38 is NOT achieved by
this natural construction.  Concretely shows the N=12 arithmetic obstruction
(opus-S118 non-monotonic emptiness) vs N=13 prime-q success.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

def Mfast(S):
    S=sorted(set(S)); Q=set()
    for v in S: Q.add(2*v)
    for a,b in combinations(S,2): Q.add(a+b); Q.add(abs(a-b))
    Q.discard(0); best=F(0)
    for q in Q:
        for a in range(1,q):
            mn=min(min((v*a)%q,q-((v*a)%q)) for v in S)
            val=F(mn,q)
            if val>best: best=val
    return best

def factor(m):
    f={}; d=2
    while d*d<=m:
        while m%d==0: f[d]=f.get(d,0)+1; m//=d
        d+=1
    if m>1: f[m]=f.get(m,0)+1
    return f

if __name__=="__main__":
    print("S_N={1,..,N-2}u{N,3(N-1)}: M(S_N) vs mediant 3/(3N+2); q=3N+2\n")
    for N in range(4,17):
        S=sorted(set(list(range(1,N-1))+[N,3*(N-1)]))
        g=reduce(gcd,S); Sp=[x//g for x in S]
        M=Mfast(Sp); med=F(3,3*N+2); lo,hi=F(1,N+1),F(2,2*N+1)
        fac="*".join(f"{p}^{e}" if e>1 else str(p) for p,e in sorted(factor(3*N+2).items()))
        tag="MEDIANT (in gap)" if M==med and lo<M<hi else ("loose" if M>=hi else "")
        print(f"  N={N:>2}: M={str(M):>7} med={str(med):>7} q=3N+2={3*N+2}={fac:<8} {tag}")
    print("=> mediant achieved only at N=7(q23),13(q41) [prime, N==1 mod6]; N=12 gives 3/35 (intruding witness q=35=5*7).")
