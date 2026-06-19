#!/usr/bin/env python3
"""
Verify THM-536-B2 N*(k) = largest N with meas(S7(AP_{N+1})) <= cap_k, for k=8,9,10,
where AP_{N+1}={0,1,...,N}.  Subset domination: E subset {0..N} => meas(S7(E))<=meas(S7(AP_{N+1})).
So span<=N* is PROVED-certified.  Confirm N*=7,8,10 (claim) and show the cliff at N*+1.
"""
from fractions import Fraction as F
from itertools import combinations

def J(A,E):
    E=sorted(set(E)); arcs=[(F(j,7),F(j+1,7)) for j in A]; bp=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for (a,b) in arcs:
            for end in (a,b):
                m=0
                while True:
                    xv=(end+m)/e
                    if xv>=1:break
                    if xv>=0:bp.add(xv)
                    m+=1
    bp=sorted(b for b in bp if 0<=b<1); tot=F(0)
    for lo,hi in zip(bp,bp[1:]+[F(1)]):
        if hi<=lo:continue
        mid=(lo+hi)/2
        if all(not any(a<((e*mid)%1)<b for (a,b) in arcs) for e in E): tot+=hi-lo
    return tot
def measS7(E): return sum(((-1)**r)*J(set(A),E) for r in range(7) for A in combinations(range(1,7),r))

CAP={8:F(2243,5880),9:F(1979,4004),10:F(55,91)}
CLAIM_NSTAR={8:7,9:8,10:10}

# meas(S7(AP_{N+1})) is independent of k (it's just the full consec cluster of size N+1).
# But cap_k depends on k.  N*(k) = largest N with meas(S7({0..N})) <= cap_k.
# Precompute meas(S7({0..N})) for N up to 25.
print("meas(S7({0..N})) for N=6..22:")
ms={}
for N in range(6,23):
    ms[N]=measS7(list(range(N+1)))
    print(f"  N={N:2d}  size={N+1:2d}  meas={float(ms[N]):.5f} ({ms[N]})")

print()
for k in (8,9,10):
    cap=CAP[k]
    Nstar=None
    for N in range(6,23):
        if ms[N]<=cap:
            Nstar=N
        else:
            break
    cliff = ms.get(Nstar+1)
    print(f"k={k}: cap={float(cap):.5f} ({cap})  N*={Nstar} (claim {CLAIM_NSTAR[k]})  "
          f"meas(S7(AP_{Nstar+1}))={float(ms[Nstar]):.5f}<=cap ; "
          f"meas(S7(AP_{Nstar+2}))={float(cliff):.5f}{'>' if cliff>cap else '<='}cap")
    print(f"      N*==claim? {Nstar==CLAIM_NSTAR[k]}")
    # also: at k=8, is N*=7=k-1 (only the AP itself certified by THM-536-B2)?
    if k==8:
        print(f"      k=8 special: N*=={k-1}? {Nstar==k-1}  (THM-536-B2 certifies ONLY the AP, span=7)")
