#!/usr/bin/env python3
"""
lrc14_mu17_gapbound_macmini_0618s5.py  (mac-mini-2026-06-18-S5)

Toward PROVING HYP-2602 (consecutive minimizes mu_{1/7}) by shrinking the residual:
TEST: "if E has a sorted-gap >= G then mu_{1/7}(E) >= thr_k" for small explicit G.
If true with small G, the residual = shapes with ALL sorted-gaps < G (bounded local
structure). Also test the SHARPER claim: mu_{1/7}(E) >= mu_{1/7}(consec_k) whenever E has
a gap >= G (so consec, the all-gaps-=1 shape, is approached from above as gaps shrink).
 (A) for each k=8..12, min mu_{1/7}(E) over shapes with max-sorted-gap = G, G=1..6.
     (G=1 is consec.) Does min mu rise as G grows past 1? (=> only G=1..G0 can be near-min)
 (B) the exact thr_k and whether min-mu-at-gap-G >= thr_k for G>=2 (residual = G=1 nbhd).
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
SEV=F(1,7); H=F(1,14)
def good_set_exact(E, thr=SEV):
    E=sorted(set(E)); k=len(E); diffs=set()
    for a in range(k):
        for b in range(a+1,k): diffs.add(E[b]-E[a])
    bps=set([F(0),F(1)])
    for d in diffs:
        for m in range(0,d+1): bps.add(F(m,d))
    bps=sorted(x for x in bps if 0<=x<=1); good=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        pts=sorted(((e*xm)%1,e) for e in E); order=[e for _,e in pts]
        fl=[int((e*xm)//1) for e in order]
        for idx in range(k):
            ec=order[idx]; fc=fl[idx]
            if idx<k-1: en=order[idx+1]; fn=fl[idx+1]; wrap=F(0)
            else: en=order[0]; fn=fl[0]; wrap=F(1)
            A=F(en-ec); Cc=F(fc-fn)+wrap
            if A==0:
                if Cc>thr: good.append((x0,x1))
                continue
            xb=(thr-Cc)/A
            lo,hi=(max(x0,xb),x1) if A>0 else (x0,min(x1,xb))
            if lo<hi: good.append((lo,hi))
    good=sorted(good); out=[]
    for a,b in good:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def mu17F(E): return sum((b-a for a,b in good_set_exact(E)),F(0))

# thr_k via the m_k sequence (min meas G_P), recomputed
def danger_arcs(u,h=H):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def measGP(P):
    if not P: return F(1)
    dz=merge([iv for u in P for iv in danger_arcs(u)]); s=F(0); prev=F(0)
    for a,b in dz:
        if a>prev: s+=a-prev
        prev=max(prev,b)
    if prev<1: s+=1-prev
    return s
import functools
@functools.lru_cache(None)
def thr(k):
    psz=13-k
    if psz==0: return F(0)
    return 1-min(measGP(P) for P in itertools.combinations(range(1,14),psz))

def maxgap_offsets(E):
    E=sorted(E); return max(E[i+1]-E[i] for i in range(len(E)-1))

print("="*82)
print("(A,B) min mu_{1/7}(E) over shapes whose MAX sorted-gap = G, vs thr_k  (k=8..12)")
print("="*82)
for k in range(8,13):
    tk=thr(k); cons=mu17F(list(range(k)))
    print(f"\n k={k}: thr_k={float(tk):.4f}, mu17(consec=gap1)={float(cons):.4f}")
    # enumerate shapes with bounded spread, group by max-sorted-gap
    SP = k+5   # spread cap for the search
    best={}    # G -> (min mu, E)
    for rest in itertools.combinations(range(1,SP), k-1):
        E=[0]+list(rest); G=maxgap_offsets(E); m=mu17F(E)
        if G not in best or m<best[G][0]: best[G]=(m,E)
    for G in sorted(best):
        if G>6: continue
        m,E=best[G]
        tag = "= consec" if G==1 else (">=thr_k OK" if m>=tk else "BELOW thr_k!")
        tag2 = ">=consec" if m>=cons-F(1,10**9) else "BELOW consec"
        print(f"    max-gap G={G}: min mu17={float(m):.4f}  [{tag}; {tag2}]  worst E={E}")
print("\nINTERPRETATION: if min-mu at G>=2 is >= thr_k, the residual is the G=1 (consec)")
print("neighbourhood only; if also >= mu(consec), consecutive-minimizes holds shape-locally.")
print("\nDONE.")
