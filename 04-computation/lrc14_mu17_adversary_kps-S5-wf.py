#!/usr/bin/env python3
"""
lrc14_mu17_adversary_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5)

DECISIVE TEST: can an aggressive hostile local-descent MINIMIZE mu_{1/7}(E) below thr_k,
the way the analogous descent broke mu_{2/7} (drove it to ~0.0436 << 1/14)?

If mu_{1/7} also collapses under large-spread descent, the 1/7 global-witness route is ALSO
broken and LRC(14) is NOT reduced. If mu_{1/7} stays >= thr_k (with healthy slack) under the
SAME attack that crushed mu_{2/7}, the 1/7 route's spread bound is strongly corroborated.

We use exact Fraction mu for confirmation but float screening for the descent (speed).
thr_k computed exactly. We minimize mu_{1/7} over integer E, 0 in E, |E|=k, large spread,
via multi-restart randomized coordinate descent (the exact attack profile that worked on 2/7).
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(20260618)

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

def mu17_float(E, N=120000):
    Ef=[float(e) for e in E]
    th=1.0/7.0
    cnt=0
    rnd=random.Random(777)
    for _ in range(N):
        x=rnd.random()
        pts=sorted((e*x)%1.0 for e in Ef)
        mg=0.0
        for t in range(len(pts)-1):
            d=pts[t+1]-pts[t]
            if d>mg: mg=d
        d=pts[0]+1.0-pts[-1]
        if d>mg: mg=d
        if mg>th: cnt+=1
    return cnt/N

# thresholds exact
thr={}
for k in range(8,14):
    psz=13-k
    thr[k]=1-min(meas(safe_set(list(P))) for P in itertools.combinations(range(1,14),psz))

def descend(k, spread, restarts=40, iters=300):
    """Randomized coordinate descent minimizing mu17_float over integer E, 0 in E, max<=spread."""
    best_val=2.0; best_E=None
    for r in range(restarts):
        E=[0]+sorted(random.sample(range(1,spread+1), k-1))
        E=sorted(set(E))
        if len(E)<k:
            continue
        cur=mu17_float(E, N=40000)
        for it in range(iters):
            idx=random.randrange(1,len(E))  # never move the 0
            old=E[idx]
            cand=old+random.choice([-3,-2,-1,1,2,3,-spread//4,spread//4])
            if cand<=0 or cand>spread: continue
            if cand in E: continue
            E2=sorted(E[:idx]+[cand]+E[idx+1:])
            v=mu17_float(E2, N=40000)
            if v<cur:
                E=E2; cur=v
        if cur<best_val:
            best_val=cur; best_E=E[:]
    return best_val, best_E

print("=== ADVERSARIAL DESCENT ON mu_{1/7}: try to break mu_1/7(E) >= thr_k ===")
print("(Same attack profile that drove mu_{2/7} to ~0.0436 << 1/14)\n")
for k in [8,9,10,11,12,13]:
    print(f"--- k={k}, thr_k = {float(thr[k]):.4f} ({thr[k]}) ---")
    for spread in [k+2, 2*k, 3*k, 4*k]:
        bv, bE = descend(k, spread, restarts=25, iters=200)
        if bE is None:
            print(f"   spread<= {spread:3d}: (no valid restart)"); continue
        ex = mu17_exact(bE)
        flag = "*** BELOW thr_k ***" if ex < thr[k] else "ok"
        print(f"   spread<= {spread:3d}: float_min~{bv:.4f}  exact={float(ex):.4f} ({ex})  E={bE}  {flag}")
    print()
print("VERDICT: if every exact mu_1/7 stays >= thr_k, the 1/7 spread floor survives the attack")
print("that destroyed the 2/7 floor.")
