#!/usr/bin/env python3
"""
lrc14_arccount_census_klein_S289.py
===================================
klein-2026-07-13-S289 (owner: prove the arc-count bound r < 3sqrt2 v|G'_{~v}|).

The certificate (THM-732) fires iff  r < 3 sqrt2 v |G'_{~v}|,  r=#arcs of G'_{~v}=good set of W=S\{v},
v=far element=max(S). This census over covering 13-sets: (1) finds the provable UPPER BOUND on r;
(2) locates EXACTLY which covering sets FAIL the certificate (expected: only small-max families);
(3) measures the max-value V0 of max(S) among failures (the finite exceptional list lives at max<=V0).
No FFT needed: r and |G'| are read off a fine grid directly.
"""
import numpy as np, itertools
NG=1<<20
THR=1.0/14.0
t=np.arange(NG,dtype=np.float64)/NG
FR={}  # cache w -> good indicator of speed w
def gind_speed(w):
    if w not in FR:
        fr=(w*t)%1.0; d=np.minimum(fr,1.0-fr); FR[w]=(d>=THR)
    return FR[w]
def good(W):
    g=np.ones(NG,dtype=bool)
    for w in W: g&=gind_speed(w)
    return g
def stats(W):
    g=good(W); L=g.mean()
    gi=g.astype(np.int8); r=int(np.sum(np.diff(np.concatenate([gi,gi[:1]]))==1))
    return L,r
def is_cov(S): return all(any(x%q==0 for x in S) for q in range(2,15))

from math import sqrt
C=3*sqrt(2)

def check(S):
    v=max(S); W=[w for w in S if w!=v]; L,r=stats(W)
    rhs=C*v*L; ratio=(r/rhs) if rhs>0 else 9.9
    return v,W,L,r,ratio, (r<rhs)

# (1) EXHAUSTIVE small-max: all covering 13-subsets of {1..N}
print("="*100)
print("EXHAUSTIVE covering 13-subsets of {1..N}: certificate r<3sqrt2 v|G'| pass/fail, grouped by max(S)")
print("="*100)
fails=[]; allr=[]
for N in [14,15,16,17,18]:
    tot=0; fail=0; maxratio=0; worst=None
    for S in itertools.combinations(range(1,N+1),13):
        if max(S)!=N: continue          # only NEW sets (max exactly N) to avoid recount
        if not is_cov(S): continue
        tot+=1
        v,W,L,r,ratio,ok=check(S); allr.append((r,max(W),L,S))
        if not ok: fail+=1; fails.append((S,v,r,L,ratio))
        if ratio>maxratio: maxratio=ratio; worst=(S,v,r,L,ratio)
    if tot:
        print("N=%2d: %4d covering sets(max=N), %3d FAIL cert; worst ratio=%.3f at %s (v=%d,r=%d,|G'|=%.4f)"%(
            N,tot,fail,maxratio,worst[0],worst[1],worst[2],worst[3]))
print("-"*100)
print("TOTAL failures (max<=18): %d"%len(fails))
for S,v,r,L,ratio in sorted(fails,key=lambda x:-x[4])[:20]:
    print("   FAIL ratio=%.3f  S=%s  v=%d r=%d |G'|=%.4f"%(ratio,S,v,r,L))

# (2) large-max known families (should all PASS)
print("\n"+"="*100)
print("LARGE-max covering families (peel far element) -- expected ALL pass:")
print("="*100)
for S in [(1,2,3,4,5,6,7,8,9,10,11,12,182),(1,2,3,4,5,6,7,8,9,10,11,13,84),
          (1,2,3,4,5,6,7,8,9,10,12,13,154),(1,2,3,4,5,6,7,8,9,11,12,13,70),
          (1,2,3,4,5,6,7,8,9,10,11,13,168),(1,2,3,4,5,6,7,8,9,10,11,12,364)]:
    v,W,L,r,ratio,ok=check(S)
    print("   %-40s v=%3d r=%2d |G'|=%.4f ratio=%.3f %s"%(str(S),v,r,L,ratio,"PASS" if ok else "**FAIL**"))

# (3) provable r upper bound: r vs max(W), 2*max(W), sum-of-speeds
print("\n"+"="*100)
print("r UPPER-BOUND candidates over ALL W seen (max<=18 census):  is r <= max(W)? 2max(W)?")
print("="*100)
viol1=[(r,mw,S) for r,mw,L,S in allr if r>mw]
viol2=[(r,mw,S) for r,mw,L,S in allr if r>2*mw]
rmax=max(allr,key=lambda x:x[0])
print("   max r observed = %d (at W of %s, max(W)=%d)"%(rmax[0],rmax[3],rmax[1]))
print("   #(r > max(W))   = %d / %d"%(len(viol1),len(allr)))
print("   #(r > 2*max(W)) = %d / %d"%(len(viol2),len(allr)))
if viol1[:5]: print("   examples r>max(W):",[(v[0],v[1]) for v in viol1[:5]])
print("\ndone.")
