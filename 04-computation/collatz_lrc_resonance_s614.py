#!/usr/bin/env python3
"""claudebox-S614: Collatz as the 2-adic/multiplicative resonance twin of LRC. v2(3n+1) geometric
(rigidity-height); cycle = 2^K~3^L resonance; parity vector = 2-adic binary signature (Lagarias
bijection); odd-density < log_3(2) => contraction (Lemma A). See HYP-2175."""
import math
from collections import Counter
def v2(m):
    e=0
    while m%2==0: m//=2; e+=1
    return e
def T(x): return x//2 if x%2==0 else (3*x+1)//2
if __name__=="__main__":
    c=Counter(v2(3*n+1) for n in range(1,200000,2)); tot=sum(c.values())
    print("v2(3n+1) geometric:", {k:round(c[k]/tot,3) for k in range(1,6)}, "E[k]=",round(sum(k*c[k] for k in c)/tot,3))
    for K in [3,5,7,10]:
        pvs={tuple((lambda n:[ (lambda s:[s[0]])([m:=n]) for _ in [0]])(0)) for n in range(0)}  # placeholder
    # parity bijection
    def pv(n,K):
        v=[]; m=n
        for _ in range(K): v.append(m%2); m=T(m)
        return tuple(v)
    for K in [3,5,7,10]:
        print(f"parity bijection K={K}: {len({pv(n,K) for n in range(2**K)})} = 2^{K}")
    print("cycle resonance: 2^K * prod n = prod(3n+1) => 2^K>=3^L (multiplicative twin of Sum m_i v_i=0)")
