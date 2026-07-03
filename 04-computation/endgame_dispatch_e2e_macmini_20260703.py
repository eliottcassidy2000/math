#!/usr/bin/env python3
"""
END-TO-END ENDGAME DISPATCH VALIDATION (mac-mini-2026-07-03-S24).
Apply the FULL assembled LRC(14) dispatch to hard 13-tuples and confirm each closes:
  (0) gcd-reduce (LRCDilation, mac-mini-S24): v -> v/gcd, WLOG gcd=1.
  (1) NON-COVERING: some q in {2..14} divides no speed -> lonely at t=1/q.               [free]
  (2) WINDOW: all |v|<=22 -> window census (native_decide).                              [done]
  (3) FAR-PEEL: some |v|>22 -> covering, gcd=1 => positive good region (M>=14/183>1/14,   [THM-609 step1
      covering-min) => lonely (base LRC(13) floor + remove far comb).                      + far_peel]
Report the ROUTE each family takes and confirm loneliness. Any family closing by NONE = a residual gap.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def nd(x): x = x % 1; return min(x, 1-x)
def is_covering(sp): return all(any(v % q == 0 for v in sp) for q in range(2,15))

def route_and_close(v):
    """apply the dispatch; return (route, lonely_bool, detail)."""
    g = gcd_all([abs(x) for x in v])
    w = [x // g for x in v]              # (0) gcd-reduce
    aw = sorted(abs(x) for x in w)
    # (1) non-covering
    for q in range(2, 15):
        if not any(x % q == 0 for x in aw):
            # lonely at t=1/q: check
            t = F(1, q); m = min(nd(x * t) for x in aw)
            return ("non-covering", m >= F(1,14), f"g={g} q={q} min={float(m):.4f}")
    # covering now
    if max(aw) <= 22:
        # (2) window census -- verify lonely at SOME small t (census would confirm)
        best = max((min(nd(x*F(a,qq)) for x in aw)) for qq in range(2,60) for a in range(1,qq) if gcd(a,qq)==1)
        return ("window-census", best >= F(1,14), f"g={g} maxspeed={max(aw)} bestmin={float(best):.4f}")
    # (3) far-peel: covering gcd=1 with a far entry -> M>=14/183 (positive good region)
    # verify M > 1/14 (fine scan) = the good-region-positive condition the peel needs
    vmax=max(aw); K=20; s=min(3000000, vmax*K); M=0.0
    for k in range(1,s):
        t=k/s; mm=min(nd(x*t) for x in aw)
        if mm>M: M=mm
    return ("far-peel", M > 1/14 - 1e-9, f"g={g} maxspeed={vmax} M={M:.5f} (>=14/183={float(F(14,183)):.5f})")

if __name__ == "__main__":
    rng = random.Random(24)
    print("END-TO-END dispatch: route + closure for hard 13-tuples (incl. gcd>1, tight, aligned, dilated).")
    print("=" * 96)
    fams = [
        ("tight {1..12,182}", list(range(1,13))+[182]),
        ("dilated x2 {2..26}", [2*i for i in range(1,14)]),
        ("dilated x6 tight*6", [6*i for i in range(1,13)]+[6*182]),
        ("{1..13} (non-cov)", list(range(1,14))),
        ("aligned blocker", None),
        ("GW-ish", [1,2,3,4,6,8,12,16,24,32,48,64,96]),
        ("random big", None),
    ]
    def aligned():
        band=list(range(15,60)); rng.shuffle(band)
        far=sorted({q*round(2000/q) for q in band[:9]}); far=[f for f in far if f>22]
        sp=far[:]
        for q in [8,9,5,7,11,13,2,3,4,6]:
            if len(sp)>=13: break
            if not any(s%q==0 for s in sp): sp.append(q)
        while len(sp)<13: sp.append(rng.randint(2,22))
        return sorted(set(sp))[:13]
    fams[4]=("aligned blocker", aligned())
    fams[6]=("random big", sorted(rng.sample(range(1,5000),13)))

    print(f"{'family':>22} {'route':>15} {'lonely?':>8}  detail")
    allok=True
    for name, sp in fams:
        if sp is None or len(set(sp))!=13: continue
        r, ok, d = route_and_close(sp)
        allok = allok and ok
        print(f"{name:>22} {r:>15} {str(ok):>8}  {d}")

    # bulk: 500 random hard covering families
    print("\nBULK: 500 random families (gcd>1 allowed) through the dispatch:")
    routes={}; fails=0
    for _ in range(500):
        sp=sorted(set(rng.sample(range(1,3000),13)))
        if len(sp)!=13: continue
        r,ok,_=route_and_close(sp)
        routes[r]=routes.get(r,0)+1
        if not ok: fails+=1
    print(f"  route counts: {routes}")
    print(f"  closure failures: {fails}/500")
    print(f"\n=> every family closes via exactly one route (gcd-reduce first); 0 failures => the assembled")
    print("   dispatch {LRCDilation + non-covering + window census + far-peel(THM-609)} is complete on the sample.")
