#!/usr/bin/env python3
"""The extremal family stratifies by 2-adic index and is self-similar under doubling
(v->2v, t->t/2) -- a dyadic fractal on top of the n->n+2 ladder. opus-2026-06-03-S579."""
from fractions import Fraction as F
from math import gcd
def dist(x): x%=1; return min(x,1-x)
def v2(x):
    c=0
    while x%2==0: x//=2; c+=1
    return c
def safe_measure(V,delta):
    eps=set([F(0)])
    for v in V:
        for k in range(v+1):
            for s in(1,-1): eps.add((F(k)+s*delta)/v % 1)
    pts=sorted(eps); meas=F(0); L=len(pts)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(v*mid)>=delta for v in V): meas+=ln
    return meas
def main():
    print("(1) DOUBLING self-similarity: even layer of AP_{2n-1} at t=1/(2n) == AP_{n-1} at t=1/n")
    for n in [4,6,7,8,10,14]:
        AP_n=list(range(1,n))
        ph_n=sorted(dist(F(i,n)) for i in AP_n)              # AP_{n-1} phases at t=1/n
        AP_2n=list(range(1,2*n)); even=[v for v in AP_2n if v%2==0]
        ph_even=sorted(dist(F(v,2*n)) for v in even)         # even layer of AP_{2n-1} at t=1/(2n)
        print(f"   n={n:2d}: AP_n phases == even-layer phases of AP_2n : {ph_n==ph_even}")
    print()
    print("(2) 2-adic STRATIFICATION of the AP sleeves at lonely time t=1/n (n=k+1)")
    for n in [8,12,14,16]:
        AP=list(range(1,n)); strata={}
        for v in AP: strata.setdefault(v2(v),[]).append(v)
        desc="; ".join(f"v2={r}:{strata[r]}" for r in sorted(strata))
        print(f"   n={n:2d}: {desc}")
    print()
    print("(3) SATURATION-curve self-similarity: mu_j(AP_{2n-1}) restricted to even sleeves vs AP_{n-1}")
    for n in [4,6,8]:
        d2=F(1,2*n); dn=F(1,n)
        even=[v for v in range(1,2*n) if v%2==0]
        cur_even=[float(safe_measure(even[:j],d2)) for j in range(1,len(even)+1)]
        cur_n=[float(safe_measure(list(range(1,n))[:j],dn)) for j in range(1,n)]
        print(f"   n={n}: AP_n curve     ={['%.3f'%x for x in cur_n]}")
        print(f"        even(AP_2n) curve={['%.3f'%x for x in cur_even]}  (same shape => self-similar)")
if __name__=='__main__': main()

def binding_analysis():
    print()
    print("(4) WHICH 2-adic stratum BINDS at the lonely time t=1/n (gap = delta)?")
    for n in [8,12,14,16]:
        AP=list(range(1,n)); delta=F(1,n)
        binders=[v for v in AP if dist(F(v,n))==delta]
        safe_margin={}
        for v in AP:
            r=v2(v); ph=dist(F(v,n))
            safe_margin.setdefault(r,[]).append(float(ph))
        mins={r:min(safe_margin[r]) for r in safe_margin}
        print(f"   n={n:2d}: binders(phase=delta)={binders} (all 2-adic val {sorted(set(v2(v) for v in binders))}); "
              f"per-stratum min phase={ {r:round(mins[r],3) for r in sorted(mins)} } (delta={float(delta):.3f})")
    print()
    print("(5) DYADIC tower n=2^a * q (odd q): even layer = scaled AP_q-level; odd layer = the binder")
    for n in [14,28,12,24]:
        q=n; a=0
        while q%2==0: q//=2; a+=1
        AP=list(range(1,n)); delta=F(1,n)
        odd=[v for v in AP if v%2==1]; even=[v for v in AP if v%2==0]
        odd_min=min(float(dist(F(v,n))) for v in odd)
        even_min=min(float(dist(F(v,n))) for v in even)
        print(f"   n={n:2d}=2^{a}*{q}: odd-layer min phase={odd_min:.4f} (binds={odd_min==float(delta)}); "
              f"even-layer min phase={even_min:.4f} (slack x{even_min/float(delta):.1f})")
binding_analysis()
