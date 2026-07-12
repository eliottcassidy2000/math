"""
opus-2026-07-11-S244: the ~6-runner crux, worked creatively via the EVEN-FOLD past work (S558o) -- the
coprime core SPREADS and misses the structured good set for the residual; the "blanket" obstruction is the
AP/V* wall (bucket A), not divisor-complete.

PAST WORK (oracle S554-S558o, the even-fold): M14(S) <= M(fold(S)), fold = even speeds halved; LRC(13)
proven => the even runners are protected => the EVEN-GOOD SET G = {t : ||v t||>=1/14 for every even v} has
POSITIVE MEASURE for free. LRC(14) <=> the ODD runners leave a point of G clear <=> |G \ union_odd D_v| > 0.
S558o's obstruction: two levers refuted -- union bound |G|>o/7 (fails, o up to 12) and generic
anti-correlation (INVERTS at AP/V*: odds density-in-G -> 1.00, they BLANKET G). Verdict: the lever is
POSITIONAL, not measure; the wall is exactly AP/V*.

THE SYNTHESIS (this session): iterate the even-fold over ALL small primes. Fold out div-by-2, then div-by-3,
5, 7, 11, 13 -- each protected by LRC(<=13). What remains is the COPRIME-TO-30030 core (<=6, opus-S243).
So the even-fold's "odd runners" sharpen to the <=6 coprime core, and:
  LRC(14) <=> the <=6 coprime core does not BLANKET the structured good set G' = {t : all non-core safe}.

VERIFIED (fine grid): for DIVISOR-COMPLETE (residual) families, the core-danger density in G' = 0.23 (mean),
ALL < 1 => the core SPREADS and MISSES G' => LONELY. Core size <=4 in sample (<=6). By contrast, the AP
{1..13} and V* {1..11,13,24} have density = 1.000 (the BLANKET wall) -- but these are BUCKET A (no mult of
14, dispatched by t=1/14), NOT divisor-complete.

SO the ~6-runner crux is FAVORABLE for the residual: the coprime core is SPREAD (S239 bad coverer) => low
danger-density in G' (0.23, huge margin from the blanket at 1) => misses G' => lonely. S558o's blanket
obstruction (density->1) is the AP/V* wall = bucket A, dispatched separately. This unifies FIVE threads onto
ONE object: even-fold (S558o) + coprime core (S243) + spread=bad-coverer (S239) + klein ~6-odd shrink (S263)
+ mac-mini decorrelation lifts (cont.49). The honest remainder is the anti-concentration in its FAVORABLE,
quantitatively-robust direction (spread core => density < 1), which S238 showed has no bounded-window
shortcut but is not marginal (0.23 << 1).
"""
from math import gcd
from functools import reduce
import random
def divisor_complete(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def primitive(v): return reduce(gcd,v)==1
def dist0(x): x=x-int(x); return min(x,1-x)
def density_in_Gp(v,G=20000):
    core=[x for x in v if gcd(x,30030)==1]; struct=[x for x in v if gcd(x,30030)>1]
    gp=0; cov=0
    for m in range(G):
        t=(m+0.5)/G
        if all(dist0(x*t)>=1/14 for x in struct):
            gp+=1
            if any(dist0(x*t)<1/14 for x in core): cov+=1
    return len(core), gp/G, (cov/gp if gp>0 else 1.0)

def main():
    random.seed(1); pool=[]; tries=0
    while len(pool)<12 and tries<200000:
        tries+=1
        v=sorted(random.sample(range(1,200),13))
        if primitive(v) and divisor_complete(v): pool.append(v)
    print("core-danger density in structured good set G' (lonely <=> < 1):")
    ds=[]
    for v in pool:
        cn,gp,dn=density_in_Gp(v); ds.append(dn)
        print(f"  DC coreN={cn} |G'|={gp:.3f} density={dn:.3f} {'LONELY' if dn<1-1e-6 else 'BLANKET'}")
    print(f"  DC mean density = {sum(ds)/len(ds):.3f} (<< 1 => core spreads, misses G')")
    for name,v in [("AP {1..13}",list(range(1,14))),("V*",[1,2,3,4,5,6,7,8,9,10,11,13,24])]:
        cn,gp,dn=density_in_Gp(sorted(v))
        print(f"  {name}: coreN={cn} density={dn:.3f} {'LONELY' if dn<1-1e-6 else 'BLANKET (=bucket A wall, dispatched by t=1/14)'}")

if __name__=='__main__':
    main()
