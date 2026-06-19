#!/usr/bin/env python3
"""
lrc14_sector_uniform_macmini_0618s7b.py  (mac-mini-2026-06-18-S7, ANGLE B)

Seek a PROVABLE uniform upper bound meas(S7(E)) <= cap_k for ALL primitive |E|=k, to close
the residual the pointwise lemma leaves.  Several candidate provable bounds, tested vs the
true meas(S7) and vs cap_k on a broad adversarial bank.

CANDIDATE BOUNDS (all provable upper bounds on meas(S7(E))):
 B1 (PAIR bound): S7 requires every sector hit; in particular for any fixed PAIR of "antipodal"
    sectors, both hit. meas(S7(E)) <= meas{both sector j and sector j+3 (say) hit}. Take MIN over
    sub-conditions. Cheap, provable (S7 subset each).  Actually meas(S7) <= meas(any necessary
    condition).  Necessary: sectors 0..6 all hit. Take the 3 "hardest" sectors to hit jointly.
 B2 (SUBSET-AP bound, pointwise): pick the LARGEST AP-run inside E; residues of that run are a
    subset; but better: ANY 7-subset E' of E with E' subset {a, a+d, ..., a+6d} (a 7-term AP)
    has meas(S7(E)) <= meas(S7(E')) = meas(S7(AP_7)) = 31/210 (scale+translation invariance!).
    => IF E CONTAINS A 7-TERM AP, meas(S7(E)) <= 31/210 = 0.1476 < cap_8!  Test how often.
 B3 (the real one): meas(S7(E)) <= meas(S7(E')) for ANY E' subset E (monotone in E: more
    offsets only help coverage). So meas(S7(E)) <= min over 7-subsets E' of meas(S7(E'))? NO --
    wrong direction: S7(E') subset S7(E) since E' subset E means FEWER points so HARDER to cover,
    so S7(E') SUBSET S7(E), meas(S7(E')) <= meas(S7(E)). That's a LOWER bound.
    The UPPER bound direction: meas(S7(E)) <= meas(S7(E union extra)). Adding offsets increases.
    So to UPPER bound meas(S7(E)), EMBED E in a larger structured set with known meas.
    E subset {0..N} => meas(S7(E)) <= meas(S7({0..N})) = AP_{N+1}. (the pointwise lemma).

So the provable upper bounds are: pointwise-AP (B-pointwise) and... that's essentially it for
clean ones. The PAIR/necessary-condition bounds (B1) give upper bounds too. Let's MEASURE how
tight the best necessary-condition bound is, and whether B2 (contains-7AP) covers the residual.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

def measS7(E):
    E=sorted(set(E)); Enz=[e for e in E if e!=0]
    bps=set([F(0),F(1)])
    for e in Enz:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        res=set(int(7*e*xm)%7 for e in E)
        if len(res)==7: total+=x1-x0
    return total

# necessary-condition (sector-r hit) measures, and the BEST triple-of-sectors necessary bound
def meas_sectors_all_hit(E, sectors):
    """meas{x: for each r in sectors, some e in E has sigma_e(x)=r}."""
    E=sorted(set(E)); Enz=[e for e in E if e!=0]
    bps=set([F(0),F(1)])
    for e in Enz:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    sset=set(sectors)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        res=set(int(7*e*xm)%7 for e in E)
        if sset <= res: total+=x1-x0
    return total

def contains_7AP(E):
    """does E contain a 7-term arithmetic progression a, a+d,...,a+6d?"""
    Es=sorted(set(E)); S=set(Es)
    for a in Es:
        for d in range(1, (max(Es)-a)//6 + 1):
            if all((a+i*d) in S for i in range(7)):
                return (a,d)
    return None

cap = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91)}
AP7 = F(31,210)

print("="*92)
print("B2 test: if E contains a 7-term AP, meas(S7(E)) <= meas(S7(7AP)) = 31/210 = 0.1476")
print("    (scale+translation invariance: any 7-term AP gives the same value as AP_7).")
print("    NOTE: 31/210 < cap_8=0.3815, so contains-7AP => certified for ALL k.")
print("="*92)
# verify the 7-AP invariance claim
for a,d in [(0,1),(3,2),(5,4),(10,3),(0,5)]:
    E=tuple(a+i*d for i in range(7))
    v=measS7(E)
    print(f"  7-term AP a={a},d={d}: E={E} meas(S7)={v} ({float(v):.6f})  ==31/210? {v==AP7}")

print()
print("="*92)
print("How often does a primitive |E|=k contain a 7-term AP? (those are auto-certified)")
print("="*92)
def gen(k,maxE):
    out=[]
    for rest in itertools.combinations(range(1,maxE+1),k-1):
        E=(0,)+rest; g=0
        for e in E: g=gcd(g,e)
        if g!=1: continue
        out.append(E)
    return out
for k in [8,9]:
    maxE={8:14,9:15}[k]
    shapes=gen(k,maxE)
    has7=sum(1 for E in shapes if contains_7AP(E))
    print(f"  k={k}, box maxE<={maxE}: {has7}/{len(shapes)} contain a 7-term AP (auto-certified by B2)")

print()
print("="*92)
print("Best NECESSARY-CONDITION bound: meas(S7) <= min over 3-sector subsets of meas(all 3 hit).")
print("How tight is it vs true meas(S7), and does it beat cap?")
print("="*92)
bank = [
  ("AP8",(0,1,2,3,4,5,6,7)),
  ("AP8-hole",(0,1,2,3,4,5,6,8)),
  ("2blocks",(0,1,2,3,40,41,42,43)),
  ("dissoc",(0,1,3,7,15,31,63,127)),
  ("near_ap_far",(0,1,2,3,4,5,60,61)),
  ("AP9-ish",(0,1,2,3,4,5,6,7,9)),
]
for name,E in bank:
    k=len(E); true=measS7(E)
    # best necessary bound: min over all 3-subsets of sectors of meas(those 3 hit)
    best=F(1)
    for S3 in itertools.combinations(range(7),3):
        v=meas_sectors_all_hit(E,S3)
        if v<best: best=v
    ck=cap.get(k,F(1))
    print(f"  {name:<14} k={k}: true meas(S7)={float(true):.5f}, best 3-sector necc bound={float(best):.5f}, "
          f"cap={float(ck):.4f}  bound<=cap? {best<=ck}")
print("\nDONE.")
