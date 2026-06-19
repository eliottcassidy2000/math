"""
Push mu(E) for k=13 as LOW as possible via aggressive multi-restart local descent +
structured constructions. CRITICAL question: can mu(E) drop below 1/14 ~ 0.071429?

If yes, the pure floor lemma B(k) with c0 >= anything near F(13) is FALSE, and moreover
the per-shape window arguments cannot save it. We must find the true behavior.

We also test whether mu->0 is approachable (no positive floor) using structured
"runner-evasion" shapes: pick E so that for most x, all phase points cluster (low maxgap).
kps-S5-wf
"""
from fractions import Fraction as F
from math import gcd
import random, time
import numpy as np

def mu_fast_np(E, N=40000):
    E = np.array(sorted(set(int(e) for e in E)), dtype=np.float64)
    x=(np.arange(N)+0.5)/N
    P=np.mod(np.outer(E,x),1.0); P.sort(axis=0)
    gaps=np.diff(P,axis=0); wrap=P[0]+1.0-P[-1]
    mg=np.maximum(gaps.max(axis=0),wrap)
    return float(np.mean(mg>2.0/7.0))

def mu_exact(E):
    E=sorted(set(int(e) for e in E)); k=len(E)
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

def descend(start, maxspread, rounds=40, grid=30000):
    cur=sorted(set(start)); curmu=mu_exact(cur)
    for r in range(rounds):
        improved=False
        for idx in range(1,len(cur)):
            cand_best=None
            for newv in range(1, maxspread+1):
                if newv in cur: continue
                cand=sorted(set(cur[:idx]+[newv]+cur[idx+1:]))
                if len(cand)!=len(cur) or not primitive(cand): continue
                mf=mu_fast_np(cand,grid)
                if mf < float(curmu)-0.0008:
                    me=mu_exact(cand)
                    if me<curmu:
                        cand_best=(me,cand); break
            if cand_best:
                curmu,cur=cand_best; improved=True; break
        if not improved: break
    return curmu,cur

if __name__=="__main__":
    random.seed(123)
    k=13
    THRESH114=F(1,14)
    print(f"1/14 = {float(THRESH114):.6f}", flush=True)
    global_best=F(3303,52780); global_bestE=[0,7,13,14,15,16,17,19,21,24,27,29,40]
    print(f"prior champion: {float(global_best):.6f}", flush=True)

    # multi-restart descent from random seeds, allowing spread up to 60
    t0=time.time()
    for trial in range(40):
        S=random.randint(25,55)
        start=sorted(set([0]+random.sample(range(1,S+1),12)))
        if len(start)!=13: continue
        mu,E=descend(start, S+8, rounds=25, grid=25000)
        if mu<global_best:
            global_best,global_bestE=mu,E
            print(f"  [{trial}] NEW GLOBAL LOW {float(mu):.6f} ({mu}) E={E}  spread={max(E)}", flush=True)
        if time.time()-t0>300:
            print("  (time budget for restarts reached)", flush=True)
            break

    # also descend the prior champion further with bigger spread allowance
    mu,E=descend(global_bestE, max(global_bestE)+25, rounds=40, grid=30000)
    if mu<global_best:
        global_best,global_bestE=mu,E
    print(f"\n==== LOWEST mu(k=13) found: {global_best}={float(global_best):.6f} at {global_bestE} ====", flush=True)
    print(f"  spread = {max(global_bestE)}", flush=True)
    print(f"  Below 1/14 ({float(THRESH114):.6f})? {global_best < THRESH114}", flush=True)
    print(f"  Below cap-14 (0.1532)? {global_best < F(5367,35035)}", flush=True)

    # Structured: scaled copies / subtori to test approach to 0
    print("\nStructured 'evasion' constructions (approach to 0?):", flush=True)
    cons=[]
    # tight cluster + few outliers at large scale
    for M in [50,100,200,400,800]:
        E=sorted(set([0,1,2,3,4,5,6,7,8,9,10,11,M]))
        cons.append((f"cluster0-11 + {M}", E))
    # arithmetic-ish dense band scaled
    for sc in [10,20,50]:
        E=sorted(set([0,sc,2*sc,3*sc,3*sc+1,3*sc+2,3*sc+3,3*sc+4,3*sc+5,3*sc+6,3*sc+7,3*sc+8,3*sc+9]))
        cons.append((f"band scale {sc}", E))
    for name,E in cons:
        if len(E)!=13 or not primitive(E):
            print(f"  {name}: invalid", flush=True); continue
        me=mu_exact(E)
        print(f"  {name}: E={E} mu={float(me):.6f} ({me})", flush=True)
