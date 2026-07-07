"""
mac-mini-2026-07-07-S41 (HYP-4837, part 4) -- seeding the FREIMAN-DICHOTOMY completion
of the k=13 tail floor, on top of kps-S59's diameter floor (HYP-4797).

kps-S59 proved: mu_{1/7}(E) >= mu_{1/7}(AP_{D+1}) for primitive diameter D, >= m_P for
D <= 75 (superset monotonicity: E subset {0..D}, finer orbit => smaller maxgap).
Residual: D > 75.  PROPOSED COMPLETION (this session):
  EITHER E is covered by a small 2-dim generalized AP  G = {i*d1 + j*d2 : 0<=i<n1, 0<=j<n2}
     -> mu_{1/7}(E) >= mu_{1/7}(G) by the SAME superset monotonicity; G's orbit is a sumset
        of two AP orbits (2-dim three-distance, bounded complexity) -> a finite mu-ledger
        over (n1, n2, d2/d1-class);
  OR E has small additive structure (not GAP-coverable; Freiman)
     -> sparse balanced lattice -> E[U]/PZ near the iid values (HYP-4837 frame) -> floor.
QUESTION HERE: is the GAP side plausible? i.e. do 2-dim GAPs (the coarse covers of the
structured large-diameter families) keep mu_{1/7} >= m_P robustly, including the shapes
that cover the known record families?
NB: mu(G) depends on (d1,d2) only through the primitive pair class; scan d2/d1.
"""
import numpy as np
from math import gcd

GRID=120_000
xs=(np.arange(GRID)+0.5)/GRID
THETA=1/7; MP=14249/252252

def mu_EU(pts):
    pts=sorted(set(pts))
    ph=np.mod(np.outer(xs,np.array(pts,float)),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    mg=g.max(axis=1)
    U=np.clip(g-THETA,0,None).sum(axis=1)
    return (mg>THETA).mean(), U.mean()

print("=== 2-dim GAP mu-ledger probe: G={i*d1+j*d2}, sizes covering >=13 points ===")
print(f"m_P = {MP:.5f}")
print(f"{'(n1,n2)':>8s} {'d1':>4s} {'d2':>5s} {'#pts':>5s} {'mu_1/7':>7s} {'E[U]':>7s} {'mu>=m_P?':>8s}")
CASES=[]
for (n1,n2) in [(13,1),(7,2),(5,3),(4,4),(9,2),(11,2),(3,5),(2,7)]:
    for (d1,d2) in [(1,50),(1,76),(1,100),(1,137),(2,77),(3,100),(1,29),(5,77)]:
        if gcd(d1,d2)!=1: continue
        pts=sorted(set(i*d1+j*d2 for i in range(n1) for j in range(n2)))
        if len(pts)<13: continue
        mu,EU=mu_EU(pts)
        CASES.append(((n1,n2),d1,d2,len(pts),mu,EU))
        print(f"({n1:2d},{n2:2d}) {d1:4d} {d2:5d} {len(pts):5d} {mu:7.4f} {EU:7.4f} {str(mu>=MP):>8s}")

worst=min(CASES,key=lambda c:c[4])
print(f"\nworst GAP mu = {worst[4]:.4f} at (n1,n2)={worst[0]}, d=(({worst[1]},{worst[2]}))  "
      f"({worst[4]/MP:.1f}x over m_P)")

# GAP covers of the known record families (do the covers stay above the bar?)
print("\n=== GAP covers of the record families ===")
FAMS={
 'monad record': [2,4,6,8,10,11,12,13,14,16,18,20,22],
 'EU-min 3-adic': [0,30,36,45,50,54,60,63,70,72,81,90,108],
 'death-star': [2*i for i in range(1,13)]+[13],
}
for name,E in FAMS.items():
    # crude 2-dim cover: base step d1 = gcd of most differences; here use the known
    # structure: evens+odds -> {i*1 + j*? }; use d1 = smallest positive difference class
    E0=[e-min(E) for e in E]
    mu,EU=mu_EU(E0)
    # cover by {i*g + j*h}: try g = most common small difference, h = the 'patch' offset
    diffs=sorted({b-a for a in E0 for b in E0 if b>a})
    g0=diffs[0]
    # minimal GAP cover heuristic: G = {i*g0 : i <= D/g0} u shifted copy
    D=max(E0)
    cover=sorted(set(list(range(0,D+1,g0)) + [e for e in E0]))
    muC,EUC=mu_EU(cover)
    print(f"{name:16s} family mu={mu:.4f}; heuristic cover ({len(cover)} pts) mu={muC:.4f} "
          f">= m_P? {muC>=MP}")
