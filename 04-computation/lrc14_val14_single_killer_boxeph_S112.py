#!/usr/bin/env python3
"""
INV val=14 is the single-killer case (essentially done)  (boxeph-2026-07-18-S112)

val=14 <=> M=14/183 <=> single-killer (v_f > 13*max(core)) = the deep-well scale = the
RESOLVED case of INV (THM-724 uniqueness + HYP-4382 + THM-1007), NOT the open compact crux.
"""
from math import gcd
from fractions import Fraction as Fr

def Mstar(V, QMAX=400):
    best=(Fr(0),1,0)
    for q in range(2,QMAX+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((a*v)%q, q-((a*v)%q)) for v in V)
            if Fr(m,q)>best[0]: best=(Fr(m,q),q,a)
    return best

def is_AP(C):
    C=sorted(C); d=C[1]-C[0]; return all(C[i+1]-C[i]==d for i in range(len(C)-1))

print('val=14 (M=14/183) families are SINGLE-KILLER (v_f > 13*max(core)):')
for name,V in [('deep well {1..12,182}',list(range(1,13))+[182]),
               ('dilate x2 {2..24,364}',[2*k for k in range(1,13)]+[364])]:
    M,q,a=Mstar(V); C=sorted(V)[:-1]; vf=max(V)
    print(f'  {name}: M={M}={float(M):.4f} core-AP={is_AP(C)} v_f={vf}>13*max(core)={13*max(C)} => SK={vf>13*max(C)}')
print()
print('HYP-4382 (M(C)=1/13 <=> C dilated {1..12}) -- non-AP cores are NOT tight:')
for C in [[1,2,3,4,5,6,7,8,9,10,11,13],[1,2,3,4,5,6,8,9,10,11,12,24]]:
    M,q,a=Mstar(C)
    print(f'  C={C}: M(C)={M}={float(M):.4f} (=1/13? {M==Fr(1,13)}, AP? {is_AP(C)})')
print()
print('=> INV val=14 = THM-724 deep-well uniqueness: interval-core (Case1) + dilated (Case2, M>=1/13)')
print('   + tight-non-dilated EMPTY (Case3, HYP-4382), residual empirical but M>1/14 unconditional (THM-1007).')
print('   The OPEN part of INV is the COMPACT (rho<13) case, which is NOT val=14.')
