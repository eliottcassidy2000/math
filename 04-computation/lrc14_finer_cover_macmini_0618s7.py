#!/usr/bin/env python3
"""
lrc14_finer_cover_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

CREATIVE IMPROVEMENT to the seven-sector route (THM-532 honest gap): S7 is the CRUDEST finite
cover of the 1/7-net N. Use a FINER cover S_L: L equally-spaced arcs A_j=[j/L, j/L+1/7),
j=0..L-1 (overlapping for L>7). A net hits every length-1/7 arc, so
   N(E) ⊆ S_L(E) ⊆ S_7(E),  and meas(S_L) DECREASES to meas(N) as L->inf.
Finer cover = TIGHTER upper bound on meas(N) = MORE slack against cap_k = closes the crude
corr<=C*W gap. Compute meas(S_L) for L=7,14,21,28,42,56 vs meas(N) and cap_k, for consec/
perforated/dissociated k=8..10. Find the smallest L whose AP value has comfortable slack,
and whether the main term stays tiny.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
SEV=F(1,7); H=F(1,14)

def measSL(E, L):
    """meas{x: every arc [j/L, j/L+1/7) (j=0..L-1) hit by {frac(e x): e in E}}."""
    E=sorted(set(E))
    # breakpoints: a point frac(e x) crosses an arc endpoint j/L or j/L+1/7 when
    # e x ≡ j/L (mod 1) or e x ≡ j/L+1/7 (mod 1) -> x = (j/L + s)/e + m/e. Collect all.
    bps=set([F(0),F(1)])
    ends=set()
    for j in range(L):
        ends.add(F(j,L)%1); ends.add((F(j,L)+SEV)%1)
    for e in E:
        if e==0: continue
        for t in ends:
            # x with frac(e x)=t  => x=(t+m)/e, m=0..e-1
            for m in range(e):
                bps.add((t+m)/e)
    bps=sorted(z for z in bps if 0<=z<=1)
    total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        pts=[(e*xm)%1 for e in E]
        # check every arc A_j hit
        ok=True
        for j in range(L):
            lo=F(j,L); hi=lo+SEV
            hit=False
            for p in pts:
                # p in [lo,hi) mod 1 ?
                if hi<=1:
                    if lo<=p<hi: hit=True; break
                else:
                    if p>=lo or p<hi-1: hit=True; break
            if not hit: ok=False; break
        if ok: total+=x1-x0
    return total

# meas(N) = 1 - mu_{1/7}
def good_set(E, thr=SEV):
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
def measN(E): return F(1)-sum((b-a for a,b in good_set(E)),F(0))

H1=F(1,14)
def danger(u):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-H1/u)%1; b=(c+H1/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def mgm(iv):
    iv=sorted(iv); o=[]
    for a,b in iv:
        if o and a<=o[-1][1]: o[-1]=(o[-1][0],max(o[-1][1],b))
        else: o.append((a,b))
    return o
def measGP(P):
    dz=mgm([iv for u in P for iv in danger(u)]); s=F(0); prev=F(0)
    for a,b in dz:
        if a>prev: s+=a-prev
        prev=max(prev,b)
    if prev<1: s+=1-prev
    return s
import functools
@functools.lru_cache(None)
def cap(k):
    return min(measGP(list(P)) for P in itertools.combinations(range(1,14),13-k))

shapes={8:[("consec",list(range(8))),("perf",[0,2,3,4,5,6,7,9]),("dissoc",[0,1,3,7,15,31,63,127])],
        9:[("consec",list(range(9)))],
        10:[("consec",list(range(10)))]}
print("meas(S_L) decreasing from S_7 toward meas(N) as L grows; slack against cap_k")
for k in (8,9,10):
    capk=cap(k)
    print(f"\n=== k={k}, cap_k={float(capk):.4f} ===")
    for name,E in shapes[k]:
        nval=measN(E)
        row=[]
        for L in (7,14,21,28,42):
            try: row.append((L,measSL(E,L)))
            except Exception as ex: row.append((L,None))
        cells=" ".join(f"S_{L}={float(v):.4f}" for L,v in row if v is not None)
        print(f"  {name:<8} {cells}  meas(N)={float(nval):.4f}")
        # slack of finest computed vs cap
        finest=[v for L,v in row if v is not None][-1]
        print(f"           finest slack cap-S_L = {float(capk-finest):.4f}  (S_7 slack {float(capk-row[0][1]):.4f})")
print("\nIf S_14 (or S_21) AP value has slack >> S_7's 0.054, the finer cover buys room to")
print("close the crude corr<=C*W certificate. Main term stays tiny (dissoc ~ small).")
print("\nDONE.")
