#!/usr/bin/env python3
"""
klein-2026-07-08-S195: R1 -- general j*=O(k). Efficient sampling version.

j*(E,V) = min{ j>=1 : maxgap{ e_i j mod V } > V/7 }.  HARD = j=1 fails.
KEY STRUCTURAL CLAIM to test: HARD (j=1 fails => V/7-dense) forces the co-offsets
to be QUASI-EVENLY-SPREAD (near-AP), so mac-mini's AP clustering (j*<=ceil(7(k-1)/6))
governs.  Tests: (1) max j* over sampled hard clusters vs bound; (2) are the
worst-j* clusters near-AP? (3) the "hard=>quasi-even" gap statistic.
"""
import numpy as np
from math import gcd, ceil
rng=np.random.default_rng(195195)

def jstar(E,V,Jmax=None):
    E=np.asarray(E); Jmax=Jmax or (3*V)
    thr=V/7.0
    for j in range(1,Jmax+1):
        p=np.sort((E*j)%V); g=np.diff(p); g=np.append(g,V-p[-1]+p[0])
        if g.max()>thr+1e-9: return j
    return None
def is_hard(E,V): return jstar(E,V,Jmax=1) is None

def sample_hard(k,V,ntarget,maxtry=200000):
    """random k-subsets with 0, spread>=6V/7, that are HARD."""
    out=[]; t=0
    while len(out)<ntarget and t<maxtry:
        t+=1
        rest=rng.choice(np.arange(1,V),size=k-1,replace=False)
        E=tuple(sorted([0]+list(int(x) for x in rest)))
        if max(E)<6*V/7: continue
        if not is_hard(E,V): continue
        out.append(E)
    return out

def gap_evenness(E,V):
    """max gap / mean gap at j=1 (1.0 = perfectly even AP-like; larger = lumpy)."""
    p=np.sort(np.array(E)%V); g=np.diff(p); g=np.append(g,V-p[-1]+p[0])
    return g.max()/g.mean()

print("(1) max j* over sampled HARD clusters, vs mac-mini AP bound ceil(7(k-1)/6):")
print(f"{'k':>3} {'V':>5} {'#hard':>6} {'maxj*':>6} {'bound':>6} {'#fail(j*>3V)':>12} {'worst evenness':>14}")
worst_overall={}
for k in (8,9,10,11,12,13):
    bound=ceil(7*(k-1)/6)
    for V in (2*k*7+1, 91, 200, 400):   # include V~14k (forces near-AP hard) + larger
        hs=sample_hard(k,V,300)
        if not hs:
            print(f"{k:>3} {V:>5} {0:>6}  (no hard sampled)"); continue
        js=[jstar(E,V) for E in hs]
        mx=max(j for j in js if j is not None)
        nf=sum(1 for j in js if j is None)
        wi=int(np.argmax([j if j else 0 for j in js]))
        ev=gap_evenness(hs[wi],V)
        flag=" <== j*>bound!" if mx>bound else ""
        print(f"{k:>3} {V:>5} {len(hs):>6} {mx:>6} {bound:>6} {nf:>12} {ev:>14.3f}{flag}")
        worst_overall[(k,V)]=(mx,hs[wi])

print("\n(2) are the WORST-j* hard clusters near-AP?  (evenness ~1 <=> AP-like)")
for k in (11,13):
    V=14*k+3
    hs=sample_hard(k,V,2000)
    scored=sorted(((jstar(E,V) or 0, gap_evenness(E,V), E) for E in hs), reverse=True)
    print(f"  k={k}, V={V}:  top hard clusters by j*:")
    for js,ev,E in scored[:4]:
        d=E[1]-E[0]; isAP=all((E[i+1]-E[i])==d for i in range(len(E)-1))
        print(f"     j*={js:2d}  evenness={ev:.3f}  {'AP(step '+str(d)+')' if isAP else 'near-AP' if ev<1.6 else 'lumpy'}  E={E[:7]}...")
    # correlation: does high j* <=> low evenness (near-AP)?
    arr=[(jstar(E,V) or 0, gap_evenness(E,V)) for E in hs]
    jsv=np.array([a for a,b in arr]); evv=np.array([b for a,b in arr])
    print(f"     corr(j*, evenness) = {np.corrcoef(jsv,evv)[0,1]:+.3f}  (negative => high j* = more even = near-AP)")

print("\n=> if maxj* <= bound everywhere (0 fails) and high-j* <=> near-AP (neg corr),")
print("   then general j*=O(k) reduces to: [hard => quasi-even (near-AP)] + [mac-mini AP clustering].")
