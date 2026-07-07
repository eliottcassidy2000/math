"""
lrc_exact_mu_ordercells_opus_S136.py

THE ORDER-CELL (MOVIE) FRAME for the density floor — the owner's tiling-model directive
applied to lonely runners — plus the first EXACT-RATIONAL mu_theta engine for GENERAL E.

FRAME (the tournament-tiling dictionary, see reflection):
  * Sorted speeds v_1<...<v_k = the base Hamiltonian path; gap vector d_j = v_j - v_{j-1}
    = the base-path arcs (13 edges on 14 events incl the observer 0).
  * Pairwise differences v_j - v_i = interval sums of d = the TILES.
  * As x varies, the config {frac(v_i x)} changes circular order exactly at the WALLS
    x = m/(v_j - v_i): one adjacent transposition = one tile flip.  The x-circle is a
    CLOSED WALK in the pair-flip hypercube (THM-373's phase clock, already canon).
  * ORDER CELLS = maximal wall-free intervals = the movie's states.  #cells (counted with
    wall multiplicity) = sum over pairs of (v_j - v_i) = the total collision count.
    The AP minimizes this (consecutive integers minimize all pairwise differences):
    #cells(AP_13) = sum_{d=1}^{12} d(13-d) = 364.  The AP is the SIMPLEST MOVIE —
    the analog of the transitive tournament = the all-zeros tiling.
  * Within a cell every gap is affine in x; maxgap = max of k affines; mu_theta(E) =
    sum over cells of the superlevel measure of that max — EXACT RATIONAL ARITHMETIC.
  * For the AP the cells are Farey-k cells and the max collapses to the roof (THM-637);
    for general E this is the "general roof".

ENGINE: adapts death-star-S1's corrected exact E[maxgap] integrator (order cells over all
Farey breakpoints m/d, d <= max difference; per-cell affine gap list; sub-breakpoints at
pairwise crossings) — replacing the integral with the exact superlevel measure.
New capability: all previous general-E mu values in the repo were grid-numeric (~4 digits).

RUNS:
 (1) VALIDATION: exact mu_{1/7}(AP_k) == roof superlevel values (691/735 ... 477/1078), and
     grid cross-checks for GW / records.
 (2) THE BOARD, EXACT: mu_{1/7} for every family on the fleet's minimizer board, as exact
     rationals (hardening the numeric census), + cell counts, distinct-difference counts,
     additive energy.
 (3) EXHAUSTIVE EXACT (A') at k=8, height <= 14: all normalized 8-subsets — verify
     mu_{1/7}(E) >= 691/735 with equality analysis (exact, so equality cases are certified,
     which no grid can do).  Extends opus-S131's numeric exhaustive to EXACT.
 (4) FRAME DIAGNOSTICS: mu vs #cells / #distinct differences / energy — the Freiman
     (difference-multiplicity) route to (A') on exact values.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import sys, time

THETA = F(1, 7)

def order_cells(E, kdenom=None):
    """Yield (a, b) consecutive Farey breakpoints covering all order changes."""
    E = list(E)
    if kdenom is None:
        kdenom = max(max(abs(j) for j in E),
                     max(abs(x - y) for x in E for y in E))
    bps = set([F(0), F(1)])
    for d in range(1, kdenom + 1):
        for m in range(1, d):
            bps.add(F(m, d))
    bps = sorted(bps)
    return list(zip(bps, bps[1:]))

def cell_gap_affines(E, a, b):
    """Affine gap functions (c, b0) valid on the open cell (a, b): gap = c*x + b0."""
    mid = (a + b) / 2
    fl = {j: (j * mid).__floor__() for j in E}
    order = sorted(E, key=lambda j: (j * mid - fl[j]))
    gaps = []
    n = len(order)
    for s in range(n):
        j1 = order[s]; j2 = order[(s + 1) % n]
        if s < n - 1:
            c = F(j2 - j1); b0 = F(-(fl[j2] - fl[j1]))
        else:
            c = F(order[0] - order[-1]); b0 = F(-(fl[order[0]] - fl[order[-1]]) + 1)
        gaps.append((c, b0))
    return gaps

def mu_exact(E, theta=THETA, kdenom=None):
    """EXACT Leb{x in [0,1]: maxgap of {frac(e x): e in E} > theta}."""
    total = F(0)
    for a, b in order_cells(E, kdenom):
        gaps = cell_gap_affines(E, a, b)
        # sub-breakpoints where the argmax can change
        subbp = set([a, b])
        for i in range(len(gaps)):
            ci, bi = gaps[i]
            for j in range(i + 1, len(gaps)):
                cj, bj = gaps[j]
                if ci != cj:
                    xc = (bj - bi) / (ci - cj)
                    if a < xc < b:
                        subbp.add(xc)
        subbp = sorted(subbp)
        for u, v in zip(subbp, subbp[1:]):
            m2 = (u + v) / 2
            cbest, bbest = max(gaps, key=lambda cb: cb[0] * m2 + cb[1])
            # superlevel of the affine cbest*x+bbest > theta on (u, v)
            fu = cbest * u + bbest; fv = cbest * v + bbest
            if fu > theta and fv > theta:
                total += v - u
            elif fu <= theta and fv <= theta:
                pass
            else:
                xs = (theta - bbest) / cbest
                if fu > theta:
                    total += max(F(0), min(v, xs) - u)
                else:
                    total += max(F(0), v - max(u, xs))
    return total

def mu_numeric(E, theta=float(THETA), res=200_000):
    tot = 0
    for j in range(res):
        x = (j + 0.6180339887) / res
        ph = sorted((e * x) % 1.0 for e in E)
        g = max(max(ph[i+1] - ph[i] for i in range(len(ph)-1)), ph[0] + 1 - ph[-1])
        if g > theta: tot += 1
    return tot / res

def normalize(E):
    E = sorted(set(E)); d = [e - E[0] for e in E[1:]]
    g = 0
    for z in d: g = gcd(g, z)
    if g == 0: g = 1
    return tuple([1] + [1 + z // g for z in d])

def cellcount(E):
    E = sorted(E)
    return sum(E[j] - E[i] for i in range(len(E)) for j in range(i+1, len(E)))

def distinct_diffs(E):
    E = sorted(E)
    return len({E[j] - E[i] for i in range(len(E)) for j in range(i+1, len(E))})

def energy(E):
    from collections import Counter
    c = Counter(a + b for a in E for b in E)
    return sum(v * v for v in c.values())

def main():
    t0 = time.time()
    print("=" * 100)
    print("(1) VALIDATION: exact mu_1/7(AP_k) vs roof values; general-E vs grid")
    print("=" * 100)
    roof_vals = {8: F(691,735), 9: F(247,294), 10: F(38,49),
                 11: F(1381,2205), 12: F(13823,24255), 13: F(477,1078)}
    for k in range(8, 14):
        m = mu_exact(list(range(1, k+1)))
        ok = "MATCH" if m == roof_vals[k] else f"*** MISMATCH vs {roof_vals[k]} ***"
        print(f"   AP_{k}: exact mu = {m} = {float(m):.6f}   {ok}")
    for name, Efam in (("GW", list(range(1,12))+[13,24]),
                       ("monad rec", [2,4,6,8,10,12,14,16,18,20,22,11,13])):
        m = mu_exact(Efam); n = mu_numeric(Efam)
        print(f"   {name}: exact mu = {m} = {float(m):.6f}   grid {n:.4f}  "
              f"({'consistent' if abs(float(m)-n) < 2e-3 else '*** CHECK ***'})")

    print()
    print("=" * 100)
    print("(2) THE BOARD, EXACT (mu_1/7 as exact rationals; frame stats)")
    print("=" * 100)
    board = {
        "AP_13":            list(range(1, 14)),
        "GW":               list(range(1, 12)) + [13, 24],
        "2AP+13":           [2,4,6,8,10,12,14,16,18,20,22,24,13],
        "monad rec":        [2,4,6,8,10,12,14,16,18,20,22,11,13],
        "kps adv42":        [2,6,8,10,11,12,14,16,18,20,22,26,42],
        "stretch":          [1,3,4,5,6,7,8,9,10,11,13,18,29],
        "primes":           [2,3,5,7,11,13,17,19,23,29,31,37,41],
        "S1min CE":         [6,8,10,11,12,13,14,16,18,25,36,43,61],
    }
    print(f"{'family':<12} {'exact mu_1/7':>28} {'~':>9} {'cells':>7} {'#diffs':>7} {'E+':>7} {'diam':>5}")
    for name, Efam in board.items():
        m = mu_exact(Efam)
        print(f"{name:<12} {str(m):>28} {float(m):>9.6f} {cellcount(Efam):>7} "
              f"{distinct_diffs(Efam):>7} {energy(Efam):>7} {max(Efam)-min(Efam):>5}")
    print(f"   [AP cell count 364 = sum_(d=1..12) d(13-d): minimal possible over 13-sets]")

    print()
    print("=" * 100)
    print("(3) EXHAUSTIVE EXACT (A') at k=8: all normalized 8-subsets of {1..14}")
    print("=" * 100)
    bar8 = F(691, 735)
    seen = set(); below = []; eq = []
    cnt = 0
    for comb in combinations(range(1, 15), 8):
        En = normalize(comb)
        if En in seen: continue
        seen.add(En)
        cnt += 1
        m = mu_exact(list(En))
        if m < bar8: below.append((En, m))
        elif m == bar8: eq.append(En)
    print(f"   normalized classes checked: {cnt}")
    print(f"   BELOW bar 691/735: {len(below)}")
    for En, m in below[:5]:
        print(f"      {En}: {m}")
    print(f"   EQUAL to bar (exact ties): {len(eq)}")
    for En in eq[:8]:
        print(f"      {En}")
    print("   => (A') k=8 verified EXACT on this box iff below-count = 0;")
    print("      equality classes certify the tight locus (grids cannot).")

    print()
    print(f"[total time {time.time()-t0:.1f}s]")

if __name__ == "__main__":
    main()
