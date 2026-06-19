#!/usr/bin/env python3
"""
DECISIVE: does the EXCESS-GAP lower bound Psi(E) -> 0 as spread -> infinity?
If yes, Psi (like the single-arc mu_0) is NOT a uniform floor -> this covering-system
functional, though a rigorous lower bound (Psi<=mu), does NOT close B(k) by itself.
If Psi stays bounded below uniformly, it WOULD close B(k).

Psi(E) = int_0^1 sum_j (gap_j(x)-2/7)_+ dx, exact rational (midpoint rule on atomic cells).
"""
from fractions import Fraction as F
from math import gcd
G0 = F(2,7)
def _frac(q): return q - q.__floor__()
def _collision_bps(E):
    E=sorted(set(E)); k=len(E); bp=set([F(0),F(1)])
    for i in range(k):
        for j in range(i+1,k):
            d=E[j]-E[i]
            if d==0: continue
            for m in range(0,d+1): bp.add(F(m,d))
    return bp
def _gap_eq_bps(E):
    E=sorted(set(E)); bp=set(); diffs=set()
    for i in range(len(E)):
        for j in range(len(E)):
            if i!=j: diffs.add(abs(E[j]-E[i]))
    for D in diffs:
        if D==0: continue
        n=0
        while True:
            cand=[(F(n)+F(2,7))/D,(F(n)+F(5,7))/D];
            for x in cand:
                if F(0)<=x<F(1): bp.add(x)
            if min(cand)>=F(1): break
            n+=1
            if n>D+2: break
    return bp
def gaps_at(E,x):
    pts=sorted(set(_frac(e*x) for e in E))
    if len(pts)==1: return [F(1)]
    g=[pts[t+1]-pts[t] for t in range(len(pts)-1)]; g.append(pts[0]+1-pts[-1]); return g
def excess_at(E,x): return sum((g-G0) for g in gaps_at(E,x) if g>G0)
def psi_exact(E):
    bp=sorted(_collision_bps(E)|_gap_eq_bps(E)); bp=[b for b in bp if F(0)<=b<=F(1)]
    if bp[0]!=F(0): bp=[F(0)]+bp
    if bp[-1]!=F(1): bp=bp+[F(1)]
    tot=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        tot+=excess_at(E,(a+b)/2)*(b-a)
    return tot
def mu_exact(E):
    bp=sorted(_collision_bps(E)|_gap_eq_bps(E)); bp=[b for b in bp if F(0)<=b<F(1)]
    tot=F(0); pts=bp+[F(1)]
    for a,b in zip(bp,pts[1:]):
        if b<=a: continue
        if max(gaps_at(E,(a+b)/2))>G0: tot+=(b-a)
    return tot
def normalize(E):
    E=sorted(set(E)); g=0
    for e in E: g=gcd(g,e)
    return [e//g for e in E] if g else E
def header(t): print("\n"+"="*74); print(t); print("="*74)

if __name__=="__main__":
    import random
    header("Psi SPREAD DECAY at fixed k -- the decisive uniformity test")
    print("  k=5 AP scaled E={0,N,2N,3N,4N}: (AP -> Psi const? since gaps scale-replicate)")
    for N in (1,2,5,10,50,200):
        E=normalize([i*N for i in range(5)])
        print(f"    N={N:4d}: Psi={float(psi_exact(E)):.6f}  mu={float(mu_exact(E)):.6f}")
    print("  k=5 perforated-AP scaled {0,N,2N,3N,5N} (introduce a hole), push N:")
    for N in (1,2,5,10,50):
        E=normalize([0,N,2*N,3*N,5*N])
        print(f"    N={N:4d}: Psi={float(psi_exact(E)):.6f}  mu={float(mu_exact(E)):.6f}")

    header("Worst-case Psi over random bounded shapes, increasing spread budget (k=13)")
    random.seed(424242)
    for sm in (1,2,4,8,16):
        worst=None; nt=600
        for _ in range(nt):
            hi=13*sm
            if hi<12: hi=12
            pool=random.sample(range(1,hi+1),12)
            E=normalize([0]+pool)
            if len(E)<13: continue
            ps=psi_exact(E)
            if worst is None or ps<worst[1]: worst=(E,ps)
        print(f"  spread_mult={sm:3d} (hi={max(13*sm,12):4d}): worst Psi={float(worst[1]):.6f}={worst[1]}")
        print(f"      E={worst[0]}")

    header("Worst-case Psi: random bounded k from 3..13, larger spread, find global min")
    random.seed(999)
    gw=None; nt=4000
    for _ in range(nt):
        k=random.randint(3,13)
        hi=random.randint(k, 6*k+20)
        pool=random.sample(range(1,hi+1),k-1)
        E=normalize([0]+pool)
        if len(E)<3: continue
        ps=psi_exact(E)
        if gw is None or ps<gw[1]: gw=(E,ps)
    print(f"  scanned {nt}: global worst Psi={float(gw[1]):.6f}={gw[1]}")
    print(f"     E={gw[0]}  mu={float(mu_exact(gw[0])):.6f}")
