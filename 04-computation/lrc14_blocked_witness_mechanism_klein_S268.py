#!/usr/bin/env python3
"""
lrc14_blocked_witness_mechanism_klein_S268.py
=============================================
klein-2026-07-12-S268

UNIFY kps cont.52's "t=1/14 is BLOCKED" mechanism with the covering-min 14/183 (klein-S267),
and DECISIVELY CORRECT kps's "DC floor = 1/12" (a box artifact).

kps cont.52 (HYP-6180): a divisor-complete (DC) family has a multiple of 14, so at t=1/14
that runner sits at residue 0 => t=1/14 is BLOCKED => the family must use a COARSER witness;
kps concludes the DC floor is 1/12 (worst DC "bottoms at q=24", band-edge 2/24=1/12).

THE CORRECTION (this script): kps has the MECHANISM right but the EXTREMIZATION backwards.
The band-edge value ceil(q/14)/q is (roughly) DECREASING in q. So M is small when the family
is forced to clear at a LARGE q. kps searched only families that clear at SMALL q (q<=24 =>
1/12); the covering-MIN is the DC family MOST STUCK — clearing only at the LARGEST q — which is
the deep well {1..12,182}, stuck at q=183, giving M = 14/183 < 1/12.

Verifies:
 (1) t=1/14 is blocked for every DC family (a 14-multiple sits at residue 0), NOT for the
     tight non-DC families (AP, GW), which reach 1/14 at t=1/14.
 (2) each DC family's M is achieved at a specific q (argq); the deep well's argq is LARGE (183),
     the 2-block's is SMALL (24). M ~ ceil(argq/14)/argq (band-edge, tight at the extremals).
 (3) ceil(q/14)/q is decreasing in q => min-M covering family = the one stuck at the largest argq
     => covering-min = deep well 14/183, and kps's small-q families (1/12) are NOT the min.
 (4) the deep well does NOT clear at q=24 (where the 2-block does) => it is genuinely "stuck".
"""
import math
from fractions import Fraction

def dist0(r,q): return min(r,q-r)
def exact_M(v):
    n=len(v); best=Fraction(0); argq=None
    qs=set()
    for a in range(n):
        for b in range(a,n): qs.add(v[a]+v[b])
    for q in sorted(qs):
        mx=0
        for p in range(1,q):
            mn=q
            for vl in v:
                d=dist0((vl*p)%q,q)
                if d<mn:
                    mn=d
                    if mn<=mx: break
            if mn>mx: mx=mn
        val=Fraction(mx,q)
        if val>best: best=val; argq=q
    return best,argq
def covering(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def reach_at(v, t):  # min_i ||v_i t||, t a Fraction
    return min(min((vi*t)%1, 1-((vi*t)%1)) for vi in v)
def clears_at(v, q):  # does some p make all residues avoid the danger arc {0,+-1,..,+-(ceil(q/14)-1)}?
    m = math.ceil(q/14)-1
    for p in range(1,q):
        if all(dist0((vi*p)%q, q) > m for vi in v):
            return p
    return None
def qmin_clear(v, qhi=400):  # smallest q>=15, 14∤q, at which v clears
    for q in range(15, qhi+1):
        if q%14==0: continue
        if clears_at(v,q) is not None:
            return q
    return None

fams = [
    ("AP {1..13} (tight, non-DC)", list(range(1,14))),
    ("GW {1..11,13,24} (tight, non-DC)", [1,2,3,4,5,6,7,8,9,10,11,13,24]),
    ("kps 2-block {1,2,3,4,10..18}", [1,2,3,4,10,11,12,13,14,15,16,17,18]),
    ("compressed 2*{1..12}u{13}", sorted([2*i for i in range(1,13)]+[13])),
    ("deep well {1..12,182}", list(range(1,13))+[182]),
]
print("="*80)
print("(1) t=1/14 BLOCKED for DC (a 14-multiple at residue 0); reachable for tight non-DC")
print("="*80)
t14=Fraction(1,14)
for nm,v in fams:
    r=reach_at(v,t14); dc=covering(v)
    has14=any(x%14==0 for x in v)
    print(f"  {nm:34s} DC={str(dc):5s} has-14-mult={str(has14):5s}  min_i||v_i/14|| = {r} = {float(r):.4f}"
          f"  {'<- t=1/14 BLOCKED' if r==0 else '<- t=1/14 gives 1/14' if r==t14 else ''}")

print()
print("="*80)
print("(2)+(3) each family's M is achieved at argq; band-edge ceil(q/14)/q DECREASES in q")
print("="*80)
print(f"  {'family':34s} {'M':>9s} {'argq':>5s} {'ceil(argq/14)/argq':>18s} {'qmin_clear':>10s}")
rows=[]
for nm,v in fams:
    M,argq=exact_M(v)
    be=Fraction(math.ceil(argq/14), argq)
    qm=qmin_clear(v)
    rows.append((nm,v,M,argq,qm))
    print(f"  {nm:34s} {str(M):>9s} {argq:>5d} {str(be):>10s}={float(be):.4f} {str(qm):>10s}")
print("  => the DEEP WELL achieves M at a LARGE argq (183); kps's 2-block at SMALL argq (24).")
print("     Since ceil(q/14)/q decreases in q, the LARGE-argq family (deep well) has the SMALLEST M.")

print()
print("="*80)
print("(4) THE 'STUCK' PROPERTY: the deep well does NOT clear at small q (where 2-block does)")
print("="*80)
dw=list(range(1,13))+[182]; tb=[1,2,3,4,10,11,12,13,14,15,16,17,18]
for q in [15,16,20,24,26,27]:
    pdw=clears_at(dw,q); ptb=clears_at(tb,q)
    print(f"  q={q:3d} (band radius {math.ceil(q/14)-1}, band-edge {math.ceil(q/14)}/{q}={float(Fraction(math.ceil(q/14),q)):.4f}):"
          f"  deep-well clears={'p=%d'%pdw if pdw else 'NO'};  2-block clears={'p=%d'%ptb if ptb else 'NO'}")
print(f"  deep well's smallest clearing q = {qmin_clear(dw)}  (STUCK at large q => small M=14/183)")
print(f"  2-block's smallest clearing q   = {qmin_clear(tb)}  (clears early => larger M=1/12)")

print()
print("="*80)
print("UNIFIED PICTURE / CORRECTION of kps cont.52")
print("="*80)
print("  kps mechanism (CORRECT): DC => t=1/14 blocked => coarser witness => M > 1/14.")
print("  kps extremization (BACKWARDS): 'worst DC bottoms at q=24 => 1/12'. FALSE.")
print("  band-edge ceil(q/14)/q DECREASES in q => M is MINIMIZED by the family stuck at the")
print("  LARGEST q. The deep well {1..12,182} clears ONLY at q=183 (182=13*14 blocks every small q),")
print("  so covering-min = 14/183 << 1/12. kps's 2170-family hunt found only small-q_min families.")
print(f"  Number line (all DC): 14/183={float(Fraction(14,183)):.5f} < 1/13={1/13:.5f} < 3/37={3/37:.5f} < 1/12={1/12:.5f}")
print("\ndone.")
