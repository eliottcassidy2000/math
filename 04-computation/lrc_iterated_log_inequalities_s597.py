#!/usr/bin/env python3
"""Iterated-log abstraction in the repo's framework. (A) worry-set prime complexity
omega(2n-1) ~ loglog n (Hardy-Ramanujan). (B) the iterated-omega 'tower height' ~ log*
(inverse tetration). opus-2026-06-03-S597."""
from math import log
from sympy import factorint, primeomega, primenu
def omega(m): return len(factorint(m))  # distinct prime factors
def iter_omega_height(m):  # apply omega until <=2; count steps (a log*-like height)
    h=0
    while m>2:
        m=omega(m); h+=1
        if m<=2: break
    return h
def main():
    print("(A) worry-set obstruction primes omega(2n-1): average vs loglog n, max order vs log/loglog")
    for N in [100,1000,5000,20000]:
        s=0; mx=0
        for n in range(4,N+1):
            w=omega(2*n-1); s+=w; mx=max(mx,w)
        avg=s/(N-3); llN=log(log(N)) if N>3 else 0
        bound=log(2*N)/log(log(2*N))
        print(f"  N={N:5d}: mean omega(2n-1)={avg:.3f}  vs loglog(N)={llN:.3f} (+Mertens const ~ {avg-llN:.2f}); "
              f"max omega={mx} vs log(2N)/loglog(2N)={bound:.1f}")
    print()
    print("(B) iterated-omega 'tower height' h(2n-1) (~ log* (2n-1), inverse tetration):")
    from collections import Counter
    for N in [1000,20000]:
        c=Counter(iter_omega_height(2*n-1) for n in range(4,N+1))
        print(f"  N={N}: height distribution {dict(sorted(c.items()))}; max height={max(c)} (log* scale, grows ULTRA-slow)")
    print()
    print("(C) the doubling order ord_n(2) (THM-404 orbit size) -- the multiplicative depth:")
    def ord2(n):
        if n%2==0: return None
        x=2%n; o=1
        while x!=1 and o<2*n: x=(2*x)%n; o+=1
        return o if x==1 else None
    for n in [7,13,27,127]:
        print(f"  n={n}: ord_n(2)={ord2(n)} (vs log2(n)={log(n)/log(2):.1f}; n-1={n-1})")
if __name__=='__main__': main()
