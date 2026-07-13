"""
opus-2026-07-11-S258: attempting to prove the LOOSE-stratum anti-concentration bound (measure{W=0}>0 for loose
covering families, S257's second half). Outcome: it is NOT easy -- corrects S257. The loose stratum reduces
(via auto-safe) to the <=6 coprime-to-30030 core anti-concentration, the project's known hard core.

W(t) = #dangerous runners at level c=14/183; measure{W=0}>0 <=> M(v)>=14/183.

(1) FAR-PEEL LEMMA (rigorous but too weak). If sub-family v has a safe interval of width >= 1/V at level c,
    then M(v u {V}) >= c (the interval holds a full period of V's comb => a safe window). BUT a 12-runner
    sub-family's safe interval at 14/183 has width w ~ 2(1/13 - c)/max(v) = 2/(2379*max(v)); far-peeling V
    needs w >= 1/V, i.e. V >= ~1190*max(rest) -- a GIANT runner. ~0% of covering families qualify (the deep
    well V=182 is correctly non-peelable). Handles only S243 Case B (giant runner).

(2) SECOND MOMENT gives the WRONG direction. Coprime family: pairwise-independent danger arcs, E[W]=2ck=1.989,
    Var(W)=2ck(1-2c)=1.685. Chebyshev: measure{W=0} <= Var/E[W]^2 = 0.43 (UPPER bound). Paley-Zygmund lower-
    bounds measure{W>=1}. Neither proves measure{W=0}>0. LLL/Janson fail (p=2c~0.15 not small, arcs heavily
    dependent). The (1-2c)^13~0.116 heuristic is only PAIRWISE independence (1-dim orbit), not a proof.

(3) WHAT REMAINS: covering => <=6 effective core (99%, verified 20664 covering families speeds<300; core-size
    dist peaks at 1-2, <=6 for 99%, only 5 have 7-8 none far-peelable). Auto-safe (S241/S243) collapses to the
    <=6 coprime-to-30030 core. BUT reducing that to M>=14/183 IS the <=6-core anti-concentration (S558o: core
    must not blanket the even-fold good set G) -- the KNOWN hard core (S242-S245), not an easy bound.

NET: S257's 'loose favorable' CORRECTED. Loose stratum = the <=6-core anti-concentration. Elementary tools all
ruled out across S256-S258 (balance-as-bound, single dual certificate, far-peel, second moment). Honest state
of the covering-min: tight/deep-well PROVED (S255); loose = <=6-core anti-concentration, OPEN. The arc maps the
difficulty precisely to where the fleet's anti-concentration thread already sits.

-> opus-S257 (split, loose half corrected), opus-S255 (deep-well proved), opus-S241/S243 (auto-safe, <=6 core),
s558o (even-fold), opus-S242-S245 (anti-concentration core), mac-mini S40.
"""
from math import gcd
from functools import reduce
from fractions import Fraction
import random
def divcomplete(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def primitive(v): return reduce(gcd,v)==1
def eff_core(v): return [x for x in v if gcd(x,30030)==1]
def main():
    c=Fraction(14,183); k=13
    EW=2*c*k; Var=2*c*k*(1-2*c)
    print(f"(2) second moment: E[W]={float(EW):.4f}, Var={float(Var):.4f}, Chebyshev measure{{W=0}}<=Var/E[W]^2={float(Var/EW**2):.4f} (UPPER bound, wrong direction)")
    print(f"    far-peel needs V >= ~1190*max(rest) (giant); deep well V=182 non-peelable (tight).")
    random.seed(21); dist={}; ntot=0
    for _ in range(120000):
        v=sorted(random.sample(range(1,300),13))
        if not primitive(v) or not divcomplete(v): continue
        ntot+=1; ec=len(eff_core(v)); dist[ec]=dist.get(ec,0)+1
        if ntot>=12000: break
    le6=sum(x for e,x in dist.items() if e<=6)
    print(f"(3) covering => <=6 effective core: {le6}/{ntot} ({100*le6//ntot}%); core-size dist {dict(sorted(dist.items()))}")
    print("    => loose stratum reduces to the <=6-core anti-concentration (S558o) = the known hard core. S257 corrected.")
if __name__=='__main__':
    main()
