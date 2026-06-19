#!/usr/bin/env python3
"""
lrc14_mu17_adversary2_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5)
Fast float descent to MINIMIZE mu_{1/7}, then exact-confirm the best per k.
Goal: break mu_1/7(E) >= thr_k if possible.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(424242)
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

# precompute a fixed grid of x for float screening (deterministic, fine)
NX=20000
XS=[(i+0.5)/NX for i in range(NX)]
def mu17_grid(E):
    Ef=[float(e) for e in E]
    th=1.0/7.0; cnt=0
    for x in XS:
        pts=sorted((e*x)%1.0 for e in Ef)
        mg=0.0
        prev=pts[0]
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

def descend(k, spread, restarts, iters):
    best=2.0; bestE=None
    for r in range(restarts):
        E=sorted(set([0]+random.sample(range(1,spread+1),k-1)))
        if len(E)<k: continue
        cur=mu17_grid(E)
        for it in range(iters):
            idx=random.randrange(1,len(E))
            old=E[idx]
            step=random.choice([-2,-1,1,2,-5,5,-spread//3 or -1,spread//3 or 1])
            cand=old+step
            if cand<=0 or cand>spread or cand in E: continue
            E2=sorted(E[:idx]+[cand]+E[idx+1:])
            v=mu17_grid(E2)
            if v<=cur:
                E=E2; cur=v
        if cur<best: best=cur; bestE=E[:]
    return best,bestE

print("=== FAST float-grid descent to MINIMIZE mu_1/7 (NX=%d), exact-confirm best ===\n"%NX)
overall_min_ratio=10.0
for k in [8,9,10,11,12,13]:
    cands=[]
    for spread in [k+2, 2*k, 3*k]:
        bv,bE=descend(k,spread,restarts=18,iters=120)
        if bE is not None: cands.append((bv,bE,spread))
    cands.sort()
    print(f"k={k}: thr_k={float(thr[k]):.4f} ({thr[k]})")
    for bv,bE,spread in cands[:3]:
        ex=mu17_exact(bE)
        below = ex<thr[k]
        ratio=float(ex/thr[k]) if thr[k]>0 else float('inf')
        overall_min_ratio=min(overall_min_ratio, ratio if thr[k]>0 else 10.0)
        print(f"   sp<={spread:3d} grid~{bv:.4f} EXACT={float(ex):.4f}({ex}) E={bE} {'*** BELOW thr ***' if below else 'ok'}")
    print()
print(f"OVERALL: min(exact mu_1/7 / thr_k) over all descents = {overall_min_ratio:.3f}  (>1 means floor SURVIVES)")
