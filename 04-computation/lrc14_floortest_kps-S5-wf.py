"""
DECISIVE TEST: does inf_E mu(E) over k=13 admissible shapes have a POSITIVE floor,
or can mu be driven to 0 (or below any c0)?

The claimed advance asserts a "moderate-spread dip strictly below F(13) then rises back",
implying a positive floor exists (just unknown). Our descent already reached ~0.056 << 1/14.

Strategy here:
 1. Confirm exact mu of the lowest champions.
 2. CONSTRUCT a family with mu provably -> 0 to settle whether B(k)>0 at all:
    Use a 2-scale construction. Take E = {0} U A U (Q + base) where points cluster.
    Specifically test E_t = {0,1,...,c} U {large widely spaced}, and "near-subgroup" shapes
    E subset of multiples that make {frac(e x)} cluster for most x.
 3. The KEY structural test: take E = d * {0,1,...,k-1}  -- by L1 mu = mu(consecutive)=const,
    so scaling alone doesn't lower it. The low-mu shapes are NON-scaled, genuinely spread.
    Test whether the descent floor stabilizes by descending with much larger spread budgets
    from the best known shape.
kps-S5-wf
"""
from fractions import Fraction as F
from math import gcd
import numpy as np, random, time

def mu_fast_np(E,N=40000):
    E=np.array(sorted(set(int(e) for e in E)),dtype=np.float64)
    x=(np.arange(N)+0.5)/N
    P=np.mod(np.outer(E,x),1.0); P.sort(axis=0)
    mg=np.maximum(np.diff(P,axis=0).max(axis=0), P[0]+1.0-P[-1])
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

if __name__=="__main__":
    # 1. confirm champions
    print("Confirm champion exact values:", flush=True)
    for E in [[0,7,13,14,15,16,17,19,21,24,27,29,40],
              [0,7,11,20,21,23,28,33,35,36,39,42,45]]:
        print(f"  E={E} mu={mu_exact(E)}={float(mu_exact(E)):.6f}  prim={primitive(E)} spread={max(E)}", flush=True)

    # 2. KEY: construct shapes engineered for low maxgap-exceedance.
    # Idea: if E points are 'balanced' (close to a perfect difference set / near regular),
    # then for most x the phase points are well-spread => maxgap SMALL => rarely > 2/7.
    # A perfect spread (all gaps ~ 1/13 < 2/7) means maxgap rarely exceeds 2/7.
    # Test: E = a near-arithmetic-progression-mod-1 generator. E = {round(13*j*alpha)} type.
    print("\nEngineered near-equidistributing shapes (low mu target):", flush=True)
    best=F(23059,412335); bestE=[0,7,11,20,21,23,28,33,35,36,39,42,45]
    # try E_j = floor(j * p / 13) style and golden-ratio spreads
    import math
    constructions=[]
    phi=(1+5**0.5)/2
    for scale in [13,26,39,52,91,143]:
        E=sorted(set(int(round(j*scale/13.0)) for j in range(13)))
        if len(E)==13: constructions.append((f"j*{scale}/13",E))
    for M in [21,34,55,89]:
        E=sorted(set(int((j*M)%(M)) for j in range(13)))  # degenerate, skip if <13
    # golden: e_j = round(j * M * frac(phi)) won't be integers nicely; use Beatty
    for M in [40,60,90,120]:
        E=sorted(set(int(round(j*M*0.3819660113)) for j in range(13)))
        if len(E)==13: constructions.append((f"Beatty {M}",E))
    for name,E in constructions:
        if not primitive(E): continue
        me=mu_exact(E)
        mark=" <== LOW" if me<best else ""
        print(f"  {name}: E={E} mu={float(me):.6f} ({me}){mark}", flush=True)
        if me<best: best,bestE=me,E

    # 3. Heavy descent from current best, large spread budget
    print("\nHeavy descent from best (spread budget +30):", flush=True)
    cur=sorted(set(bestE)); curmu=mu_exact(cur); t0=time.time()
    rounds=0
    while rounds<60 and time.time()-t0<360:
        rounds+=1; improved=False
        for idx in range(1,len(cur)):
            for newv in range(1, max(cur)+30):
                if newv in cur: continue
                cand=sorted(set(cur[:idx]+[newv]+cur[idx+1:]))
                if len(cand)!=13 or not primitive(cand): continue
                mf=mu_fast_np(cand,25000)
                if mf<float(curmu)-0.001:
                    me=mu_exact(cand)
                    if me<curmu:
                        cur,curmu=cand,me; improved=True
                        print(f"  -> {float(me):.6f} ({me}) spread={max(cand)}", flush=True)
                        break
            if improved: break
        if not improved: break
    print(f"\nLOWEST mu(k=13): {curmu}={float(curmu):.6f} at {cur} spread={max(cur)}", flush=True)
    print(f"  < 1/14? {curmu < F(1,14)}   < 1/20? {curmu < F(1,20)}   < 1/30? {curmu < F(1,30)}", flush=True)
