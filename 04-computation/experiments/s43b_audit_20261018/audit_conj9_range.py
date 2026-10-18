#!/usr/bin/env python3
"""S43b referee: THM-4686 corollary (D) F_{(Mp-1)/2} == F_{(p-1)/2} F_{(M-1)/2} (mod p^2) in the session's stated
range (M = 3..11 odd, 5 <= p <= 150), independent code (valuation-tracked central binomials mod p^2)."""
from math import isqrt
def primes_in(a,b): return [p for p in range(max(a,2),b+1) if all(p%d for d in range(2,isqrt(p)+1))]
def tvals(p, m, K):
    pp=p*p; w=pow(m%pp,-1,pp); v=0; u=1; pw=1; out=[1]
    for k in range(1,K+1):
        for x,s in ((2*k,1),(2*k-1,1),(k,-1),(k,-1)):
            e=0
            while x%p==0: x//=p; e+=1
            v+=s*e; u=u*(x if s>0 else pow(x,-1,pp))%pp
        pw=pw*w%pp
        out.append(0 if v>=1 else pow(u,3,pp)*pw%pp)
    return out
res={}
for m in [1,-8,16,-64,64,256,-512,4096,3,5,7,10]:
    fails={M:0 for M in (3,5,7,9,11)}; tested=0
    for p in primes_in(5,150):
        if m%p==0: continue
        tested+=1; pp=p*p; K=(11*p-1)//2; t=tvals(p,m,K)
        pref=[]; s=0
        for x in t: s=(s+x)%pp; pref.append(s)
        for M in fails:
            if (pref[(M*p-1)//2]-pref[(p-1)//2]*pref[(M-1)//2])%pp: fails[M]+=1
    res[m]=(tested,fails)
    print('m = %5d: primes tested %d, failures by M %s' % (m, tested, fails))
