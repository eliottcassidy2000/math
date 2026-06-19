#!/usr/bin/env python3
"""
lrc14_mu17_k13_k8hard_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5)
Finish k=13 descent and HAMMER the tightest case k=8 (largest thr_k) with heavy descent.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(20262026)
ONE7=F(1,7); H=F(1,14)
def merge(iv):
    iv=sorted(iv);out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def meas(arcs): return sum((b-a for a,b in arcs),F(0))
def complement(arcs):
    arcs=merge(arcs);out=[];prev=F(0)
    for a,b in arcs:
        if a>prev: out.append((prev,a))
        prev=max(prev,b)
    if prev<1: out.append((prev,F(1)))
    return out
def danger_arcs(u,h=H):
    iv=[]
    for j in range(u):
        c=F(j,u);a=(c-h/u)%1;b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1)));iv.append((F(0),b))
    return iv
def safe_set(P,h=H):
    if not P: return [(F(0),F(1))]
    return complement(merge([iv for u in P for iv in danger_arcs(u,h)]))
def mu17_exact(E):
    E=sorted(set(E));k=len(E);diffs=set()
    for a in range(k):
        for b in range(a+1,k): diffs.add(E[b]-E[a])
    bps=set([F(0),F(1)])
    for d in diffs:
        for m in range(0,d+1): bps.add(F(m,d))
    bps=sorted(x for x in bps if 0<=x<=1);good=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        pts=sorted(((E[t]*xm)%1,E[t]) for t in range(k))
        order=[e for _,e in pts];floors=[int((e*xm)//1) for e in order]
        for idx in range(k):
            e_cur=order[idx];f_cur=floors[idx]
            if idx<k-1: e_nx=order[idx+1];f_nx=floors[idx+1];wrap=F(0)
            else: e_nx=order[0];f_nx=floors[0];wrap=F(1)
            A=F(e_nx-e_cur);Cc=F(f_cur-f_nx)+wrap
            if A==0:
                if Cc>ONE7: good.append((x0,x1))
                continue
            xb=(ONE7-Cc)/A
            if A>0: lo=max(x0,xb);hi=x1
            else: lo=x0;hi=min(x1,xb)
            if lo<hi: good.append((lo,hi))
    return meas(merge(good))
NX=20000
XS=[(i+0.5)/NX for i in range(NX)]
def mu17_grid(E):
    Ef=[float(e) for e in E]; th=1.0/7.0; cnt=0
    for x in XS:
        pts=sorted((e*x)%1.0 for e in Ef); mg=0.0; prev=pts[0]
        for t in range(1,len(pts)):
            d=pts[t]-prev
            if d>mg: mg=d
            prev=pts[t]
        d=pts[0]+1.0-pts[-1]
        if d>mg: mg=d
        if mg>th: cnt+=1
    return cnt/NX
thr={}
for k in range(8,14):
    psz=13-k
    thr[k]=1-min(meas(safe_set(list(P))) for P in itertools.combinations(range(1,14),psz))
def descend(k,spread,restarts,iters):
    best=2.0;bestE=None
    for r in range(restarts):
        E=sorted(set([0]+random.sample(range(1,spread+1),k-1)))
        if len(E)<k: continue
        cur=mu17_grid(E)
        for it in range(iters):
            idx=random.randrange(1,len(E)); old=E[idx]
            step=random.choice([-2,-1,1,2,-5,5,-7,7])
            cand=old+step
            if cand<=0 or cand>spread or cand in E: continue
            E2=sorted(E[:idx]+[cand]+E[idx+1:]); v=mu17_grid(E2)
            if v<=cur: E=E2; cur=v
        if cur<best: best=cur; bestE=E[:]
    return best,bestE

# k=13 finish
print("=== k=13 ===")
print(f"thr_13={float(thr[13]):.4f} ({thr[13]})")
for spread in [15,26,39]:
    bv,bE=descend(13,spread,restarts=14,iters=100)
    if bE:
        ex=mu17_exact(bE)
        print(f"  sp<={spread}: grid~{bv:.4f} EXACT={float(ex):.4f}({ex}) E={bE} {'*** BELOW thr ***' if ex<thr[13] else 'ok'}")

# k=8 HARD (largest thr_k=3637/5880~0.6185): push restarts/iters and many spreads
print("\n=== k=8 HARD (tightest thr_k) ===")
print(f"thr_8={float(thr[8]):.4f} ({thr[8]}); consec mu_1/7=691/735={float(F(691,735)):.4f}")
gmin=2.0; gminE=None
for spread in [10,12,16,20,24,32,40]:
    bv,bE=descend(8,spread,restarts=40,iters=250)
    if bE:
        ex=mu17_exact(bE)
        if ex<gmin: gmin=ex; gminE=bE
        print(f"  sp<={spread:3d}: grid~{bv:.4f} EXACT={float(ex):.4f}({ex}) E={bE} {'*** BELOW thr ***' if ex<thr[8] else 'ok'}")
print(f"\n  GLOBAL k=8 descent min exact = {float(gmin):.4f} ({gmin}) at E={gminE}")
print(f"  thr_8={float(thr[8]):.4f}; slack = {float(gmin-thr[8]):.4f}  -> {'SURVIVES' if gmin>=thr[8] else 'BROKEN'}")
