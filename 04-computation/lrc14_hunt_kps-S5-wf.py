"""
ADVERSARIAL HUNT: try to drive mu(E) as LOW as possible for k=13 (and 12),
testing whether there's a positive uniform floor (lemma B(k)).

Strategies:
  (A) exhaustive spread sweep cap=15..22 for k=13 (sampled/partial where huge)
  (B) random large-spread shapes (up to spread ~120), k=13
  (C) structured: perforated APs, common-factor subtori, two-cluster, dilations
  (D) the iid-ceiling argument: F(13)=0.2263 is a CEILING -> mu cannot exceed it in the
      equidistributed limit, but can dip below. We seek the dip floor.
  (E) does mu EVER hit 0 or get arbitrarily small? L2 says mu>=5/(7 maxE) per shape, so for
      fixed maxE bounded it cannot be 0; but as spread->inf the per-shape bound ->0. We test
      whether the ACTUAL min stabilizes (positive floor) or keeps dropping.

We use a FLOAT fast mu for the search, then EXACT-confirm the best finds.
kps-S5-wf
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import random, time

# fast float mu via deterministic fine grid (lower bound on resolution but fine enough to rank)
def mu_fast(E, N=200000):
    E = sorted(set(int(e) for e in E))
    th = 2.0/7.0
    cnt = 0
    for s in range(N):
        x = (s+0.5)/N
        pts = sorted((e*x) % 1.0 for e in E)
        mg = 0.0; prev = pts[0]
        for t in range(1, len(pts)):
            g = pts[t]-prev
            if g > mg: mg = g
            prev = pts[t]
        w = pts[0]+1.0-pts[-1]
        if w > mg: mg = w
        if mg > th: cnt += 1
    return cnt/N

def mu_exact(E):
    E = sorted(set(int(e) for e in E))
    k = len(E)
    if k <= 1: return F(1) if k==1 else F(0)
    TH = F(2,7)
    bp = set([F(0),F(1)]); diffs=set()
    for i in range(k):
        for j in range(i+1,k): diffs.add(E[j]-E[i])
    for d in diffs:
        for m in range(0,d+1): bp.add(F(m,d))
    obp = sorted(b for b in bp if F(0)<=b<=F(1))
    refined=set(obp)
    for a,b in zip(obp,obp[1:]):
        mid=(a+b)/2
        floors={e:(e*mid).__floor__() for e in E}
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
        mid=(a+b)/2
        pts=sorted(set((e*mid)%1 for e in E))
        if len(pts)==1: mg=F(1)
        else:
            gaps=[pts[t+1]-pts[t] for t in range(len(pts)-1)]; gaps.append(pts[0]+1-pts[-1]); mg=max(gaps)
        if mg>TH: tot+=(b-a)
    return tot

def primitive(E):
    g=0
    for e in E: g=gcd(g,e)
    return g==1

CAP14_13 = F(5367,35035)
SPREAD18 = F(7037,59976)

if __name__ == "__main__":
    random.seed(1)
    k = 13
    best_mu = SPREAD18  # known champion
    best_E = [0,1,2,4,6,9,11,12,13,15,16,17,18]
    print(f"Starting champion (spread-18): mu={best_mu}={float(best_mu):.6f}", flush=True)

    # (C1) DILATION sweep of the cap-14 minimizer by p/q for many ratios, rounding
    E13min = [0,1,2,3,4,5,6,7,8,9,12,13,14]
    print("\n(C1) rational-dilation sweep of cap-14 minimizer (round):", flush=True)
    dil_best = None
    for q in range(2,9):
        for p in range(q+1, 4*q+1):
            lam = F(p,q)
            Ed = sorted(set(round(lam*e) for e in E13min))
            if len(Ed)!=13 or not primitive(Ed): continue
            m = mu_exact(Ed)
            if dil_best is None or m < dil_best[0]:
                dil_best = (m, Ed, lam)
    if dil_best:
        print(f"   best dilation: lam={dil_best[2]} mu={dil_best[0]}={float(dil_best[0]):.6f} E={dil_best[1]}", flush=True)
        if dil_best[0] < best_mu:
            best_mu, best_E = dil_best[0], dil_best[1]

    # (C2) random large-spread search (float screen, exact confirm top hits)
    print("\n(C2) random large-spread screen (k=13, spread up to 130):", flush=True)
    t0=time.time(); screened=[]
    TRIALS=4000
    for _ in range(TRIALS):
        S = random.randint(15, 130)
        elems = random.sample(range(1, S+1), 12)
        E = sorted(set([0]+elems))
        if len(E)!=13: continue
        if max(E)!=S:  # ensure spread
            pass
        mf = mu_fast(E, 60000)
        screened.append((mf, E))
    screened.sort(key=lambda z:z[0])
    print(f"   screened {len(screened)} in {time.time()-t0:.0f}s; exact-confirming top 12:", flush=True)
    for mf, E in screened[:12]:
        if not primitive(E): continue
        me = mu_exact(E)
        mark = " <== NEW LOW" if me < best_mu else ""
        print(f"     E={E} float={mf:.5f} exact={me}={float(me):.6f}{mark}", flush=True)
        if me < best_mu:
            best_mu, best_E = me, E

    # (C3) structured perforated-AP / two-cluster families
    print("\n(C3) structured families:", flush=True)
    fam = []
    # perforated near-AP: 0..L with a few holes, various L
    for L in range(14, 30):
        full = list(range(L+1))
        if len(full) < 13: continue
        # remove (len-13) interior points to minimize -- try a few removal patterns
        nrem = len(full)-13
        if nrem < 0 or nrem > 8: continue
        # spread out removals
        for shift in range(0, max(1,nrem+2)):
            interior = full[1:-1]
            if nrem > len(interior): continue
            step = max(1, len(interior)//max(1,nrem))
            rem = set(interior[shift::step][:nrem]) if nrem>0 else set()
            E = sorted(set(full)-rem)
            if len(E)!=13 or not primitive(E): continue
            fam.append(E)
    seen=set();
    for E in fam:
        tE=tuple(E)
        if tE in seen: continue
        seen.add(tE)
        me = mu_exact(E)
        if me < best_mu:
            print(f"   perforated NEW LOW E={E} mu={me}={float(me):.6f}", flush=True)
            best_mu, best_E = me, E

    print(f"\n==== BEST k=13 mu found: {best_mu}={float(best_mu):.6f} at {best_E} ====", flush=True)
    print(f"  cap-14 value was {CAP14_13}={float(CAP14_13):.6f}", flush=True)
    print(f"  F(13) ceiling = {F(3132376013,13841287201)}={float(F(3132376013,13841287201)):.6f}", flush=True)
