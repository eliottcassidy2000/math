from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import numpy as np

def distZ(num,den):
    r=num%den; return F(min(r,den-r),den)
def reach_M_wit(V):  # exact M + a witness t achieving it
    V=[abs(v) for v in V if v!=0]; n=len(V); dens=set()
    for i in range(n):
        dens.add(2*V[i])
        for j in range(i+1,n):
            dens.add(V[i]+V[j]); 
            if V[i]!=V[j]: dens.add(abs(V[i]-V[j]))
    best=F(0); bt=None
    for d in dens:
        if d==0: continue
        for m in range(1,d):
            mn=min(distZ(v*m,d) for v in V)
            if mn>best: best=mn; bt=F(m,d)
    return best,bt
def M_grid(V,G=2000000):  # independent fine-grid cross-check (upper-ish estimate)
    xs=(np.arange(G)+0.5)/G
    ph=np.mod(np.outer(xs,np.array(V,float)),1.0)
    d=np.minimum(ph,1-ph)
    return float(d.min(axis=1).max())

fam=[1,2,3,4,5,6,7,8,9,10,11,13,24]
M,t=reach_M_wit(fam)
print(f"=== structure of the 2nd tight family (opus-S131) ===")
print(f"V={fam}")
print(f"  exact M={M}={float(M):.6f}; witness t={t}; grid-M={M_grid(fam):.6f} (cross-check ~1/14={1/14:.6f})")
print(f"  V mod 14 = {[v%14 for v in fam]}")
print(f"  at t={t}: values v*t mod 1, dist to Z:")
for v in fam:
    print(f"    v={v:>3}: ||v*t||={distZ(v*(t.numerator),t.denominator)} = {float(distZ(v*t.numerator,t.denominator)):.4f}")

print(f"\n=== FULL census of tight (M=1/14) single-scale 13-families, bounded ===")
thr=F(1,14); tight=[]
# {1..12}+{x}, {1..11}+{x,y}, {1..10}+{x,y,z} bounded
cands=[]
for x in range(13,60): cands.append(list(range(1,13))+[x])
for x in range(12,45):
    for y in range(x+1,50): cands.append(list(range(1,12))+[x,y])
seen=set()
for V in cands:
    if len(set(V))<13: continue
    if reduce(gcd,V)!=1: continue
    key=tuple(sorted(V))
    if key in seen: continue
    seen.add(key)
    M,_=reach_M_wit(V)
    if M==thr: tight.append(key)
print(f"tested {len(seen)} bounded families; TIGHT (M=1/14 exactly): {len(tight)}")
for t in tight: 
    r=[v%14 for v in t]
    print(f"   {t}   mod14={sorted(r)}  (repeat mod14? {len(set(r))<13})")
print(f"\n  AP {{1..13}} mod14 = {sorted([v%14 for v in range(1,14)])}")
print("  => structural question: are all tight families characterized by their mod-14 residue pattern?")
