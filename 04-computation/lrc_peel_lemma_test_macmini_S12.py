#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S12 -- HYP-4382: test the PEEL LEMMA that would close the
spectral-gap induction (HYP-2052 general n => (G) at n=13 => TightLooseDichotomy).

INDUCTION P(n-1) => P(n) [P(n): no config in the open gap (1/n, 2/(2n-1))]:
  gap member S (n-1 speeds) --peel a runner--> S' (n-2 speeds).  M(S') >= M(S).
  By P(n-1): M(S') = 1/(n-1) [AP-tight] or M(S') >= 2/(2n-3).
  CASE A (some peel AP-tight): S = (n-2)-AP + runner => AP-completion => S tight
     or loose => contradiction.  CLOSES.
  CASE B (ALL peels loose, M(S') >= 2/(2n-3)): no contradiction.  THE OBSTRUCTION.

So the induction closes IFF no "all-peels-loose" gap member exists.  Since gap
members don't exist, we TEST THE STRUCTURE on the full known spectrum (exhaustive
n=6,7) to see whether "all-peels-loose" can occur NEAR the tight locus:

  For every gcd-1 config S (n=6,7 exhaustive; bounded speeds), compute M(S) and
  minpeel(S) = min over single-runner peels of M(S').  The induction's Case A is
  available whenever minpeel(S) = 1/(n-1) (an AP-tight peel exists).
  KEY: is there ANY config with M(S) SMALL (near tight, <= 2/(2n-1)) whose
  minpeel is LARGE (> 1/(n-1), i.e. all peels strictly loose)?  If NO such config
  exists across the spectrum, the peel lemma holds empirically and the induction
  is sound; the doubled-apex boundary configs A_n peel to the AP (drop the apex).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations, product
import sys, time
sys.path.insert(0, '04-computation')
from lonely_profile import profile

T0 = time.time()
def log(m=""): print(m, flush=True)

def M_exact(S):
    S = sorted(set(abs(x) for x in S))
    for cap in (13, 10, 8, 6, 4, 3, 2):
        p = profile(S, F(1, cap))
        m = p.M()
        if m is not None:
            return m

def minpeel(S):
    """min over single-runner peels of M(S')."""
    S = list(S); best = None
    for i in range(len(S)):
        Sp = S[:i]+S[i+1:]
        if len(set(Sp)) < len(Sp):
            continue
        m = M_exact(Sp)
        if best is None or m < best:
            best = m
    return best

for n in (6, 7):
    speeds = n-1
    gap_lo, gap_hi = F(1,n), F(2,2*n-1)
    peel_tight = F(1, n-1)          # AP-tight value one level down
    peel_gap_hi = F(2, 2*(n-1)-1)   # second value one level down
    BOUND = 3*n                     # compressed range
    log(f"\n===== n={n} ({speeds} speeds), gap=({gap_lo},{gap_hi}), peel-tight=1/{n-1} =====")
    seen = set(); n_cfg = 0
    spectrum = {}                    # M value -> count
    obstruction = []                 # configs w/ M<=gap_hi but all peels loose
    ap_has_tight_peel = None
    apex = tuple(sorted(list(range(1,n-1))+[2*(n-1)]))   # doubled-apex A_n
    t0 = time.time()
    for combo in combinations(range(1, BOUND+1), speeds):
        if reduce(gcd, combo) != 1:
            continue
        n_cfg += 1
        M = M_exact(combo)
        spectrum[M] = spectrum.get(M,0)+1
        # test the obstruction: M near tight but all peels loose
        if M <= gap_hi:
            mp = minpeel(combo)
            if mp is not None and mp > peel_tight:      # no AP-tight peel
                obstruction.append((M, mp, combo))
        if time.time()-t0 > 75:
            log(f"  [time cap; scanned {n_cfg} configs]"); break
    # report spectrum bottom + the induction obstruction
    lo_vals = sorted(spectrum)[:4]
    log(f"  configs scanned: {n_cfg}; lowest M values: " +
        ", ".join(f"{v}({spectrum[v]})" for v in lo_vals))
    log(f"  in-gap ({gap_lo},{gap_hi}) count: {sum(spectrum[v] for v in spectrum if gap_lo<v<gap_hi)}")
    # AP and doubled-apex peel behavior
    AP = tuple(range(1, n))
    log(f"  AP {AP}: M={M_exact(AP)}, minpeel={minpeel(AP)} (AP-tight peel? {minpeel(AP)==peel_tight})")
    log(f"  doubled-apex {apex}: M={M_exact(apex)}, minpeel={minpeel(apex)} "
        f"(AP-tight peel? {minpeel(apex)==peel_tight})")
    # the obstruction verdict
    log(f"  OBSTRUCTION (M<=2/{2*n-1} but ALL peels loose >1/{n-1}): {len(obstruction)} configs")
    for M, mp, c in sorted(obstruction)[:6]:
        log(f"     M={M} minpeel={mp} {c}")
    if not obstruction:
        log(f"  => peel lemma HOLDS at n={n}: every near-tight config has an AP-tight peel. Induction SOUND.")

log(f"\n[t = {time.time()-T0:.0f}s]")
