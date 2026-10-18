#!/usr/bin/env python3
"""S43b referee: L(m) <=> F_m == (-m/p) F_{4096/m} (mod p^2), for every m with p !| m (consequence of (R));
truth table over CM and generic m, 5 <= p <= 300."""
from math import comb, isqrt
def primes_in(a,b): return [p for p in range(max(a,2),b+1) if all(p%d for d in range(2,isqrt(p)+1))]
def leg(a,p):
    a%=p
    return 0 if a==0 else (1 if pow(a,(p-1)//2,p)==1 else -1)
rows=dis=0; cmok=0; cmrows=0
for p in primes_in(5,300):
    pp=p*p; n=(p-1)//2
    H=[0]*p
    for i in range(1,p): H[i]=(H[i-1]+pow(i,-1,p))%p
    C3=[comb(2*k,k)**3%pp for k in range(n+1)]
    for m in [1,-8,16,-64,64,256,-512,4096,3,5,7,10,-3,2,-1,100,12345,-999]:
        if m%p==0: continue
        im=pow(m%pp,-1,pp); m2=4096*im%pp; im2=pow(m2,-1,pp)
        F=sum(C3[k]*pow(im,k,pp) for k in range(n+1))%pp
        F2=sum(C3[k]*pow(im2,k,pp) for k in range(n+1))%pp
        h=sum(C3[k]*pow(im,k,pp)%p*(H[2*k]-H[k]) for k in range(n+1))%p
        q=((pow(m%pp,p-1,pp)-1)//p)%p
        L=(6*h-q*F)%p==0
        R=(F-leg(-m,p)*F2)%pp==0
        rows+=1; dis+=(L!=R)
        if m in (1,-8,16,-64,64,256,-512,4096): cmrows+=1; cmok+=R
print('L(m) <=> F_m == (-m/p) F_(4096/m) mod p^2: %d (m,p) rows, %d disagreements; holds at %d of %d CM rows' % (rows, dis, cmok, cmrows))
