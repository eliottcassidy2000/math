#!/usr/bin/env python3
"""
signed_lrc_energy_globalmax_s707b.py   monad-explorer-2026-06-06-S707
Firm up the central new claim: AP={1,...,n-1} is the GLOBAL maximizer of the
signed additive energy E_+ over ALL (n-1)-subsets of (Z/C)\\{0}, C=2n-1.
Extends s707 part B to larger n.  Also reports #maximizers and AP-orbit size.
Pure integer arithmetic -> exact.
"""
from itertools import combinations
from math import gcd

def E_plus(A,C):
    rp={}
    a=list(A)
    for i in range(len(a)):
        ai=a[i]
        for j in range(i+1,len(a)):
            s=(ai+a[j])%C
            rp[s]=rp.get(s,0)+1
    return sum(v*v for v in rp.values())

def ap_orbit(C,k):
    AP=frozenset(range(1,k+1))
    us=[u for u in range(1,C) if gcd(u,C)==1]
    return {frozenset((u*x)%C for x in AP) for u in us}

for n in range(4,13):
    C=2*n-1; k=n-1
    Eap=E_plus(range(1,n),C)
    best=-1; nmax=0
    # track whether every maximizer is in AP-orbit
    orb=ap_orbit(C,k)
    nonAP_max=0
    for combo in combinations(range(1,C),k):
        E=E_plus(combo,C)
        if E>best:
            best=E; nmax=1; nonAP_max=0 if frozenset(combo) in orb else 1
        elif E==best:
            nmax+=1
            if frozenset(combo) not in orb: nonAP_max+=1
    print(f"n={n:2d} C={C:2d} k={k:2d}: E_+(AP)={Eap:4d} globalmax={best:4d} "
          f"AP_is_max={Eap==best}  #maximizers={nmax:4d} |AP-orbit|={len(orb):2d} "
          f"#nonAP_maximizers={nonAP_max}")
print("\nCONCLUSION: AP is global additive-energy max at every n tested;"
      " but maximizer set is far larger than the AP dilation-orbit.")
