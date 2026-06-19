#!/usr/bin/env python3
"""
lrc14_asm_widehunt_kps-S10-wf.py  (fast engine)
V4/V5 of the finite-check-to-B audit: HUNT for an over-cap shape with span > B(k),
and probe the wide-spread 'stranger sup L_y' constants.

The claim's OWN gaps say Regime B (span>B) is genuinely infinite and NOT closed by proof.
This is the adversarial half: try hard to breach the cap with span>B.
Families:  resonant/AP, 1-stranger (many scales), 2-stranger, no-scale-sep (shifted AP),
random wide primitive, 7-structured (apex prime), dilated+jitter.
All EXACT (Fraction).
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, comb
from functools import reduce
import random

SECTORS=[(F(j,7),F(j+1,7)) for j in range(1,7)]
def fast_profile(E):
    E=sorted(set(E)); bp=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for j in range(1,7):
            for end in (F(j,7),F(j+1,7)):
                m=0
                while True:
                    xv=(end+m)/e
                    if xv>=1:break
                    if xv>=0:bp.add(xv)
                    m+=1
    bp=sorted(b for b in bp if 0<=b<1); meas=F(0); Sacc=[F(0)]*5
    nz=[e for e in E if e!=0]
    for lo,hi in zip(bp,bp[1:]+[F(1)]):
        if hi<=lo:continue
        L=hi-lo; mid=(lo+hi)/2; fr=[(e*mid)%1 for e in nz]; free=0
        for (a,b) in SECTORS:
            if not any(a<v<b for v in fr): free+=1
        if free==0: meas+=L
        for r in range(5):
            if free>=r: Sacc[r]+=L*comb(free,r)
    return meas,Sacc

DUAL={8:([F(1),F(-1),F(1),F(-9,10),F(3,5)],4),
      9:([F(1),F(-13,18),F(4,9),F(-1,6)],3),
      10:([F(1),F(-13,18),F(4,9),F(-1,6)],3)}
def Ly_from_S(S,k):
    y,R=DUAL[k]; return sum(y[r]*S[r] for r in range(R+1))
CAP={8:F(2243,5880),9:F(1979,4004),10:F(55,91)}
B={8:16,9:15,10:13}
def span(E): return max(E)-min(E)
def prim(E): return reduce(gcd,[e for e in E if e>0])==1

worst={k:(F(0),None) for k in (8,9,10)}
worst_ly={k:(F(0),None) for k in (8,9,10)}
breaches=[]
ncheck={8:0,9:0,10:0}

def test(k,E):
    E=sorted(set(E))
    if 0 not in E or len(E)!=k or not prim(E): return
    if span(E)<=B[k]: return
    ncheck[k]+=1
    m,S=fast_profile(E); ly=Ly_from_S(S,k)
    if m>worst[k][0]: worst[k]=(m,E)
    if ly>worst_ly[k][0]: worst_ly[k]=(ly,E)
    if m>CAP[k] or ly>CAP[k]:
        breaches.append((k,E,float(m),float(ly),m>CAP[k],ly>CAP[k]))

print("="*78)
print("V4: HUNT for over-cap shapes with span>B(k) (the OPEN region)")
print("="*78)

import sys
def flush(msg): print(msg); sys.stdout.flush()

# F1: consec_(k-1) + 1 far stranger, many scales (convergence by N~200 per HYP-2610)
for k in (8,9,10):
    base=list(range(k-1))
    for N in list(range(B[k]+1,260))+[300,400,600,800]:
        test(k,base+[N])
flush("  F1 done (1-stranger)")

# F2: 2-stranger: consec_(k-2) + two satellites
for k in (8,9,10):
    base=list(range(k-2))
    for N1 in [17,19,23,31,41,61,103,151,211]:
        for N2 in [N1+1,N1+2,N1+7,N1+13,2*N1,2*N1+1,3*N1-1,N1+50,N1+101,N1+157]:
            test(k,base+[N1,N2])
flush("  F2 done (2-stranger)")

# F3: shifted AP / tight cluster at M (no scale separation): {0}u{M..M+k-2}
for k in (8,9,10):
    for M in range(8,160):
        test(k,[0]+[M+i for i in range(k-1)])
    flush(f"  F3 k={k} done")
flush("  F3 done (shifted-AP/no-scale-sep)")

# F4: resonant - AP-like with one coprime breaker {0,d,2d,..,(k-2)d, (k-1)d+1}
for k in (8,9,10):
    for d in range(2,30):
        test(k,[0]+[i*d for i in range(1,k-1)]+[(k-1)*d+1])
        test(k,[0]+[i*d for i in range(1,k-1)]+[(k-1)*d-1])

# F5: 7-structured (apex prime) strangers
for k in (8,9,10):
    base=list(range(k-1))
    for N in [21,28,35,42,49,56,63,70,77,84,98,7*17,7*23,7*31]:
        if N>B[k]: test(k,base+[N])
    # two 7-structured
    for N1 in [14,21,28]:
        for N2 in [N1+7,N1+14,N1+21,N1+1]:
            test(k,[0]+list(range(1,k-2))+[N1,N2])

flush("  F4,F5 done (resonant/apex-7)")
# F6: random wide primitive, large span
random.seed(20260619)
for k in (8,9,10):
    for _ in range(1200):
        hi=random.choice([25,40,80,160])
        rest=sorted(random.sample(range(1,hi+1),k-1))
        E=[0]+rest
        if span(E)<=B[k]: continue
        test(k,E)
    flush(f"  F6 k={k} done")
flush("  F6 done (random wide)")

# F7: dilated-consec + jitter (short-relation rich), span pushed beyond B
for k in (8,9,10):
    for d in [2,3,4,5,7]:
        for jit in range(1,2*k):
            test(k, sorted(set([0]+[d*i for i in range(1,k-1)]+[d*(k-1)+jit])))
flush("  F7 done (dilated+jitter)")

# F8: near-consec with a moderate gap inserted (span B+1..B+8)
for k in (8,9,10):
    for extra in range(B[k]+1, B[k]+40):
        test(k, list(range(k-1))+[extra])
    # consec_(k-2) + pair straddling B
    for a in range(k-1, B[k]+5):
        for b in range(a+1, B[k]+12):
            E=sorted(set(list(range(k-2))+[a,b]))
            if len(E)==k: test(k,E)
flush("  F8 done (gap-inserted)")

print(f"\n  shapes tested (span>B): k=8:{ncheck[8]}  k=9:{ncheck[9]}  k=10:{ncheck[10]}")
print("\n--- worst meas(S7) found over span>B ---")
for k in (8,9,10):
    m,E=worst[k]
    print(f"  k={k}: max meas(S7)={float(m):.5f} ({m})  cap={float(CAP[k]):.5f}  "
          f"under-cap by {float(CAP[k]-m):.5f}  E={E}")
print("--- worst L_y found over span>B ---")
for k in (8,9,10):
    ly,E=worst_ly[k]
    print(f"  k={k}: max L_y={float(ly):.5f} ({ly})  cap={float(CAP[k]):.5f}  "
          f"under-cap by {float(CAP[k]-ly):.5f}  E={E}")
print(f"\n  >>> BREACHES (span>B, meas>cap OR L_y>cap): {len(breaches)}")
for b in breaches[:30]:
    print("   ",b)

print()
print("="*78)
print("V5: stranger sup L_y (claim 0.25704/0.48729/0.56241), consec_(k-1)+far N")
print("="*78)
for k in (8,9,10):
    base=list(range(k-1))
    vals=[]
    for N in [50,100,200,400,800,1600,3200]:
        _,S=fast_profile(base+[N]); vals.append((N,float(Ly_from_S(S,k))))
    sup=max(v for _,v in vals)
    claim=[0.25704,0.48729,0.56241][k-8]
    print(f"  k={k}: L_y(consec_{k-1}+N) N={[v[0] for v in vals]} = {[round(v[1],5) for v in vals]}")
    print(f"        observed sup~{sup:.5f}  claim {claim}  cap={float(CAP[k]):.5f}  "
          f"sup<cap? {sup<float(CAP[k])}  sup~claim? {abs(sup-claim)<2e-3}")
