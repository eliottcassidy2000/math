"""
Faster adversarial hunt for low mu(E), k=13. Vectorized float screen with numpy,
exact-confirm the lowest finds. Goal: how low can mu go? Is there a positive floor?
kps-S5-wf
"""
from fractions import Fraction as F
from math import gcd
import random, time
import numpy as np

def mu_fast_np(E, N=40000):
    E = np.array(sorted(set(int(e) for e in E)), dtype=np.float64)
    x = (np.arange(N)+0.5)/N
    # pts[i,s] = frac(E[i]*x[s])
    P = np.mod(np.outer(E, x), 1.0)           # shape (k, N)
    P.sort(axis=0)
    gaps = np.diff(P, axis=0)                  # (k-1, N)
    wrap = P[0]+1.0-P[-1]                      # (N,)
    mg = np.maximum(gaps.max(axis=0), wrap)
    return float(np.mean(mg > 2.0/7.0))

def mu_exact(E):
    E = sorted(set(int(e) for e in E)); k=len(E)
    if k<=1: return F(1) if k==1 else F(0)
    TH=F(2,7); bp=set([F(0),F(1)]); diffs=set()
    for i in range(k):
        for j in range(i+1,k): diffs.add(E[j]-E[i])
    for d in diffs:
        for m in range(0,d+1): bp.add(F(m,d))
    obp=sorted(b for b in bp if F(0)<=b<=F(1)); refined=set(obp)
    for a,b in zip(obp,obp[1:]):
        mid=(a+b)/2; floors={e:(e*mid).__floor__() for e in E}
        order=sorted(E,key=lambda e:e*mid-floors[e])
        for t in range(k):
            if t==k-1:
                ef,el=order[0],order[-1]; slope=ef-el; const=F(1)-floors[ef]+floors[el]
            else:
                eh,elo=order[t+1],order[t]; slope=eh-elo; const=-floors[eh]+floors[elo]
            if slope!=0:
                xb=(TH-const)/slope
                if a<xb<b: refined.add(xb)
    refined=sorted(refined); tot=F(0)
    for a,b in zip(refined,refined[1:]):
        mid=(a+b)/2; pts=sorted(set((e*mid)%1 for e in E))
        if len(pts)==1: mg=F(1)
        else:
            g=[pts[t+1]-pts[t] for t in range(len(pts)-1)]; g.append(pts[0]+1-pts[-1]); mg=max(g)
        if mg>TH: tot+=(b-a)
    return tot

def primitive(E):
    g=0
    for e in E: g=gcd(g,e)
    return g==1

if __name__=="__main__":
    random.seed(7)
    k=13
    best=F(7037,59976); bestE=[0,1,2,4,6,9,11,12,13,15,16,17,18]
    print(f"champion spread-18: {float(best):.6f}", flush=True)
    # big random screen, several spread bands
    t0=time.time()
    candidates=[]
    for band in [(15,25),(20,40),(30,60),(40,90),(60,150),(100,300)]:
        for _ in range(1500):
            S=random.randint(*band)
            E=sorted(set([0]+random.sample(range(1,S+1),12)))
            if len(E)!=13: continue
            mf=mu_fast_np(E,20000)
            candidates.append((mf,E))
    candidates.sort(key=lambda z:z[0])
    print(f"screened {len(candidates)} in {time.time()-t0:.0f}s. Exact-confirm top 25:", flush=True)
    for mf,E in candidates[:25]:
        if not primitive(E): continue
        me=mu_exact(E); mark=" <== NEW LOW" if me<best else ""
        print(f"  E={E} float={mf:.5f} exact={float(me):.6f} ({me}){mark}", flush=True)
        if me<best: best,bestE=me,E
    print(f"\nBEST k=13: {best}={float(best):.6f} at {bestE}", flush=True)

    # Also: local optimization around the best — try all single-element moves
    print("\nLocal descent from best:", flush=True)
    improved=True; cur=bestE[:]; curmu=best; rounds=0
    while improved and rounds<6:
        improved=False; rounds+=1
        S=max(cur)
        for idx in range(1,len(cur)):  # don't move the 0
            for newv in range(1, S+6):
                if newv in cur: continue
                cand=sorted(set(cur[:idx]+[newv]+cur[idx+1:]))
                if len(cand)!=13 or not primitive(cand): continue
                mf=mu_fast_np(cand,20000)
                if mf < float(curmu)-0.0005:
                    me=mu_exact(cand)
                    if me<curmu:
                        cur,curmu=cand,me; improved=True
                        print(f"  descent -> {float(me):.6f} ({me}) E={cand}", flush=True)
                        break
            if improved: break
    print(f"\nFINAL k=13 champion after descent: {curmu}={float(curmu):.6f} at {cur}", flush=True)
