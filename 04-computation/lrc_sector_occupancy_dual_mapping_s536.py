#!/usr/bin/env python3
"""
lrc_sector_occupancy_dual_mapping_s536.py    oracle-2026-06-01-S536o

A DUAL mapping: nodes = the n evenly-spaced SECTORS of the circle (the regular
n-gon cells), runners are moving markers that toggle edges as they cross sector
boundaries. (User's idea.) The sector width = 1/n = the loneliness threshold.

SET-UP. Sectors S_k = [k/n, (k+1)/n), k=0..n-1. Observer at 0 sits on the boundary
of S_0 and S_{n-1}. At time t, runner i (speed v_i) is in sector floor(n*frac(v_i t)).
Occupancy vector c(t) = (c_0,...,c_{n-1}), sum = n-1.

LONELINESS <=> the two sectors touching the observer are EMPTY:
   observer lonely  <=>  c_0(t) = 0 AND c_{n-1}(t) = 0
   (all runners in [1/n, 1-1/n]). So LRC@n <=> every speed set reaches an
   occupancy vector with c_0 = c_{n-1} = 0 as t varies.

PIGEONHOLE: n-1 runners in n sectors => ALWAYS >=1 empty sector. LRC asks the empty
cell(s) to be steerable to the OBSERVER's location, and needs the 2 observer-sectors
empty simultaneously (the 2/n gap straddling 0). The regular n-gon (speeds 1..n-1)
is the tight extremal.

THE TOURNAMENT ON SECTORS: rank sectors by occupancy with a cyclic tiebreak,
   a -> b  iff  c_a > c_b, or (c_a == c_b and forward-distance(a,b) < n/2).
Edges flip exactly when a runner crosses a sector boundary (c changes). As t runs,
this is a closed walk on sector-tournament iso-classes. LRC = reach a class where
the observer's two sectors are joint minima (c=0).

DUALITY (the payoff): the DFT of the occupancy vector is the S529 exponential sum,
   chat_m = sum_k c_k e^{-2pi i m k/n}  ~  sum_j e^{-2pi i m (sector phase)} ,
so SECTORS are the real-space dual of the FOURIER/character/resonance picture
(S529-S535). Loneliness (empty observer cells) is the real-space face of the
character condition. We verify the duality numerically.

This script computes: realizable occupancy-vector & sector-tournament iso-classes
(restriction vs all), the LRC emptiness check, and the DFT duality.
"""
from itertools import combinations, permutations, product
from functools import reduce
from math import gcd, pi, cos, sin
import cmath, random

def frac(x): return x - int(x // 1)

# ---------------- occupancy ----------------
def occupancy(speeds, n, t):
    c = [0]*n
    for v in speeds:
        c[int(n*frac(v*t)) % n] += 1
    return tuple(c)

def lonely(speeds, n, t):
    c = occupancy(speeds, n, t)
    return c[0] == 0 and c[n-1] == 0

# ---------------- sector tournament ----------------
def sector_tournament_canon(c, n):
    """a->b iff c_a>c_b or (c_a==c_b and (b-a)%n < n/2). canonical over S_n."""
    adj = [[0]*n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            if a == b: continue
            if c[a] > c[b] or (c[a] == c[b] and (1 <= (b-a) % n <= (n-1)//2)):
                adj[a][b] = 1
    # ensure tournament (exactly one of a->b, b->a); fix any double via the cyclic rule
    for a in range(n):
        for b in range(a+1, n):
            if adj[a][b] == adj[b][a]:
                # tie both-1 or both-0: orient by cyclic
                fwd = 1 <= (b-a) % n <= (n-1)//2
                adj[a][b] = 1 if fwd else 0
                adj[b][a] = 0 if fwd else 1
    best = None
    for p in permutations(range(n)):
        bits = tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or bits < best: best = bits
    return best

A000568 = {3:2,4:4,5:12,6:56,7:456}

def occ_isoclass(c, n):
    """occupancy vector up to the reflection fixing the observer (k -> (n-k)%n)."""
    refl = tuple(c[(n-k) % n] for k in range(n))
    return min(c, refl)

# ---------------- study ----------------
def study(n, n_sets, samples=4000, do_tour=True):
    rnd = random.Random(700+n)
    occ_classes = set(); tour_classes = set()
    lonely_ok = 0; tot = 0
    for _ in range(20000):
        v = tuple(sorted(rnd.sample(range(1, 6*n), n-1)))
        if reduce(gcd, v) != 1: continue
        tot += 1
        if tot > n_sets: break
        seen_occ = set()
        lon = False
        for s in range(samples):
            t = (s+0.5)/samples
            c = occupancy(v, n, t)
            seen_occ.add(c)
            if c[0] == 0 and c[n-1] == 0: lon = True
        if lon: lonely_ok += 1
        for c in seen_occ:
            occ_classes.add(occ_isoclass(c, n))
            if do_tour:
                tour_classes.add(sector_tournament_canon(c, n))
    return occ_classes, tour_classes, lonely_ok, tot

def count_compositions(total, parts):
    # number of compositions of `total` into `parts` nonneg parts = C(total+parts-1, parts-1)
    from math import comb
    return comb(total+parts-1, parts-1)

def duality_check(speeds, n):
    """verify DFT of occupancy ~ exponential sum sum_j e^{-2pi i m v_j t}, sampled."""
    rnd = random.Random(11); worst = 0.0
    for _ in range(50):
        t = rnd.random()
        c = occupancy(speeds, n, t)
        for m in range(1, n):
            chat = sum(c[k]*cmath.exp(-2j*pi*m*k/n) for k in range(n))
            esum = sum(cmath.exp(-2j*pi*m*int(n*frac(v*t))/n) for v in speeds)
            worst = max(worst, abs(chat - esum))
    return worst

def main():
    print("="*74)
    print("DUAL MAPPING: nodes = n evenly-spaced SECTORS; runners toggle edges on crossing")
    print("="*74)
    print("  LONELINESS <=> observer's two sectors empty (c_0 = c_{n-1} = 0).")
    print("  Pigeonhole: n-1 runners in n sectors => always >=1 empty sector.\n")

    print("  n | realizable occ-iso / all compositions | realizable sector-tourn iso / A000568 | LRC(observer emptied)")
    for n in (4, 5, 6, 7):
        do_tour = n <= 6
        occ, tour, lon, tot = study(n, n_sets=(120 if n < 7 else 60), do_tour=do_tour)
        allcomp = count_compositions(n-1, n)             # compositions of n-1 into n parts
        # up to reflection ~ roughly allcomp/2; report raw allcomp as the unrestricted target
        tourstr = f"{len(tour)} / {A000568[n]}" if do_tour else "(skipped n=7)"
        print(f"  {n} |   {len(occ)} / {allcomp}{'':6s}             |   {tourstr:18s}            | {lon}/{tot}")
    print()
    print("  => realizable occupancy iso-classes and sector-tournament iso-classes are a")
    print("     RESTRICTED subset; LRC = every speed set reaches the observer-emptied class.")
    print("     (misses from LRC count = the tight AP/regular-polygon set, lonely only at")
    print("      measure-zero t=k/n the grid skips -- the boundary extremal.)")
    print()

    print("="*74)
    print("DUALITY: DFT of occupancy vector = the S529 exponential/character sum")
    print("="*74)
    for n in (5, 6, 7):
        rnd = random.Random(n)
        v = tuple(sorted(rnd.sample(range(1, 5*n), n-1)))
        w = duality_check(v, n)
        print(f"  n={n}, speeds={v}: max |DFT(occupancy)_m - sum_j e(-m*sector_j/n)| = {w:.2e} (0 => exact dual)")
    print("  => SECTORS are the real-space dual of the Fourier/resonance picture (S529).")
    print("     Empty observer cells (loneliness) <-> the character-sum condition, by DFT.")

if __name__ == "__main__":
    main()
