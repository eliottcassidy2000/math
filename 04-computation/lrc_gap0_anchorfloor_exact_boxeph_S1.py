"""
EXACT reconciliation of the origin-gap / anchor-floor numbers for LRC(14).
boxeph-2026-07-07-S1, HYP-4760.

The fleet has CONFLICTING numbers for 'E[gap@0]':
  kps-S58 : E[gap_0] = 2 E[min] = 0.137  ( < 1/7 )   -> origin-gap reduction FALLS SHORT
  klein-S153: E[gap@0](AP)=0.211, inf_E E[gap@0] ~ 0.156 ( > 1/7 ) -> single anchor CLOSES it

Decisive question: is inf_E E[gap@0] ABOVE or BELOW 1/7 = 0.142857 ?
If ABOVE, klein's single-anchor floor closes the crux; if BELOW, need the true max.

We compute everything EXACTLY (piecewise-linear, rational breakpoints).

Conventions (SPEED config = recent mu_{1/7} convention; 0 NOT a config point):
  config C(x) = { frac(e x) : e in E }
  gap_a(x)   = length of the arc between consecutive points of C(x) that contains
               the anchor a (in R/Z).  (If a coincides with a point, measure-zero.)
  gap0(x) = gap_a with a=0 = R(x)+L(x), R=min_i frac(e_i x), L=min_i frac(-e_i x).
  E[gap0] = 2 E[min_i frac(e_i x)]  (x->-x symmetry).
"""
from fractions import Fraction as F

# ---------- exact integral of a min-of-linear (or the gap-at-anchor) ----------
def breakpoints(E, extra=()):
    bps = {F(0), F(1)}
    n = len(E)
    for i in range(n):
        ai = abs(E[i])
        if ai:
            for m in range(0, ai+1):
                bps.add(F(m, ai))
        for j in range(i+1, n):
            d = abs(E[i]-E[j])
            if d:
                for m in range(0, d+1):
                    bps.add(F(m, d))
    for a in extra:
        bps.add(F(a))
    return sorted(b for b in bps if 0 <= b <= 1)

def E_min_frac(E):
    """EXACT E_x[min_i frac(e_i x)]  (e_i != 0)."""
    E = [e for e in E if e != 0]
    bps = breakpoints(E)
    tot = F(0)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        mid = (a+b)/2
        # frac(e mid)=e*mid-floor; the argmin line on (a,b)
        vals = []
        for e in E:
            fl = (e*mid).__floor__()
            vals.append((e*mid - fl, e, fl))
        _, e, fl = min(vals)
        # line L(x)=e*x-fl on [a,b]; integral
        tot += e*(b*b-a*a)/2 - fl*(b-a)
    return tot

def E_gap_at_anchor(E, anchor):
    """EXACT E_x[ length of gap of {frac(e x)} containing 'anchor' ].
    gap = (dist from anchor to nearest point going +) + (going -).
    dist+ = min_i frac(e_i x - anchor)? No: points are frac(e_i x); anchor fixed a.
    right neighbor dist = min_i frac(frac(e_i x) - a) = min_i frac(e_i x - a).
    left  neighbor dist = min_i frac(a - e_i x) = min_i frac(-(e_i x - a)).
    We integrate rdist+ldist exactly."""
    a0 = F(anchor)
    E2 = list(E)
    # need breakpoints where order/argmin of frac(e x - a0) changes: e x - a0 in Z
    # => x = (m + a0)/e. Add those.
    bps = {F(0), F(1)}
    n = len(E2)
    for i in range(n):
        e = E2[i]
        if e != 0:
            # frac(e x - a0) wraps at x=(m+a0)/e
            # m ranges so that x in [0,1]
            import math
            mlo = math.floor(float(-a0)*1) - abs(e) - 2
            mhi = abs(e) + 2
            for m in range(mlo, mhi+1):
                x = (m + a0)/e
                if 0 <= x <= 1: bps.add(x)
        for j in range(i+1, n):
            d = abs(E2[i]-E2[j])
            if d:
                for m in range(0, d+1):
                    bps.add(F(m, d))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if lo == hi: continue
        mid = (lo+hi)/2
        # rdist(x)=min_i frac(e_i mid - a0) achieved by argmin line
        rv = []
        lv = []
        for e in E2:
            # frac(e x - a0) = e x - a0 - floor(e mid - a0)
            flr = (e*mid - a0).__floor__()
            rv.append(((e*mid - a0) - flr, e, flr))
            fll = (a0 - e*mid).__floor__()
            lv.append(((a0 - e*mid) - fll, e, fll))
        _, er, flr = min(rv)
        _, el, fll = min(lv)
        # rline = e_r x - a0 - flr ; lline = a0 - e_l x - fll
        tot += (er*(hi*hi-lo*lo)/2 - (a0+flr)*(hi-lo))
        tot += (-el*(hi*hi-lo*lo)/2 + (a0-fll)*(hi-lo))
    return tot

# ---------------- families ----------------
ap13 = list(range(1,14))
gw   = [1,2,3,4,5,6,7,8,9,10,11,13,24]

print("="*68)
print("E[gap@0] (speed config, 0 NOT a point):  E[gap0] = 2 E[min_i frac(e_i x)]")
print("="*68)
for nm, E in [("AP {1..13}",ap13), ("GW {1..11,13,24}",gw)]:
    em = E_min_frac(E)
    g0 = 2*em
    print(f"  {nm:20s}  E[min]={float(em):.6f}  E[gap0]=2E[min]={str(g0)} = {float(g0):.6f}")
print(f"  1/7 = {1/7:.6f}")
print("  (kps-S58 said E[gap0]=0.137; klein-S153 said E[gap@0](AP)=0.211. Which matches?)")

# cross-check with direct anchor computation at a=0
print("\n  cross-check E_gap_at_anchor(.,0):")
for nm,E in [("AP",ap13),("GW",gw)]:
    print(f"    {nm}: {float(E_gap_at_anchor(E,0)):.6f}")

print("\n" + "="*68)
print("DECISIVE: inf_E E[gap@0] over 13-families -- above or below 1/7?")
print("="*68)
import random
random.seed(11)
def E_gap0_fast(E, G=30000):
    import math
    s=0.0
    for aa in range(G):
        x=(aa+0.5)/G
        mn=2.0; mn2=2.0
        for e in E:
            fe=(e*x)%1.0
            if fe<mn: mn=fe
            g=1.0-fe
            if g<mn2: mn2=g
        s+=mn+mn2
    return s/G
worst=(1.0,None)
below=0
# structured candidates + random
cands=[ap13, gw, [1,2,3,4,5,6,7,8,9,10,11,12,26], [2*j for j in range(1,14)],
       [1,2,3,4,5,6,20,21,22,23,24,25,26]]
for _ in range(200):
    cands.append(sorted(random.sample(range(1,50),13)))
for E in cands:
    v=E_gap0_fast(E)
    if v<worst[0]: worst=(v,E)
    if v<1/7: below+=1
print(f"  scanned {len(cands)} families; {below} have E[gap0] < 1/7")
print(f"  worst (min) E[gap0] = {worst[0]:.6f} at {worst[1]}")
print(f"  => inf_E E[gap0] is {'BELOW' if worst[0]<1/7 else 'ABOVE'} 1/7 on this sample")

print("\n" + "="*68)
print("ANCHOR FLOOR  E[max_{a in A} gap_a]  -- inf over families, growing A")
print("="*68)
def E_anchorfloor_fast(E, anchors, G=20000):
    s=0.0
    for aa in range(G):
        x=(aa+0.5)/G
        pts=sorted((e*x)%1.0 for e in E)
        best=0.0
        for a in anchors:
            # gap containing a
            # find neighbors
            lo=None;hi=None
            for p in pts:
                if p<=a:
                    lo=p
            # circular
            # right neighbor = first pt > a else pts[0]+1
            rn=None
            for p in pts:
                if p>a: rn=p;break
            if rn is None: rn=pts[0]+1.0
            # left neighbor = last pt <= a else pts[-1]-1
            ln=None
            for p in pts:
                if p<=a: ln=p
            if ln is None: ln=pts[-1]-1.0
            g=rn-ln
            if g>best: best=g
        s+=best
    return s/G
anchor_sets = {
    "{0}":[0.0],
    "{0,1/2}":[0.0,0.5],
    "{0,1/3,2/3}":[0.0,1/3,2/3],
    "{j/8}":[j/8 for j in range(8)],
}
for label,A in anchor_sets.items():
    worst=(1.0,None)
    below=0
    for E in cands:
        v=E_anchorfloor_fast(E,A)
        if v<worst[0]: worst=(v,E)
        if v<1/7: below+=1
    print(f"  A={label:12s} inf E[max_a gap_a] = {worst[0]:.5f}  ({below}/{len(cands)} below 1/7)  min-at {worst[1] if worst[0]<0.18 else '...'}")
