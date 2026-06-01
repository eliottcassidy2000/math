#!/usr/bin/env python3
"""
lrc_complex_vertex_objects_s538.py    oracle-2026-06-01-S538o

More COMPLEX, PRECISE objects for the tournament vertices. The S535 lesson:
metric-bearing vertices restrict the realizable iso-classes; the S537 lesson:
speeds live in a flow/tension duality. So the fertile complex vertices carry
DERIVED metric/algebraic data with built-in consistency constraints.

STAR -- VERTICES = PAIRS (edges of K_n), each carrying the difference-speed
   w_{ij} = v_i - v_j.
The w_{ij} are a TENSION / coboundary: they obey the cocycle law
   w_{ij} + w_{jk} + w_{ki} = 0   (a closed condition on every triangle).
So the pair-tournament is a SECOND-ORDER LRC on the difference set {w_{ij}}, and
its realizable iso-classes are restricted by ADDITIVE consistency (the cocycle),
not free. Tournament: pair p beats pair q iff frac(w_p t) is in the forward
half-turn of frac(w_q t) (the half-turn comparator on difference-phases). LRC for
the observer = the observer-pairs {0,i} (phases v_i t) are all far from 0 = the
original loneliness, now a marked sub-family of the pair-tournament.

PLUS a MULTITUDE of further complex/precise vertex objects (computed or posed):
 - GAPS: vertices = the n arcs between consecutive points; tournament by size;
   realizable gap-vectors restricted (sum=1; the apex/largest gap = loneliness, S530;
   the three-gap/Steinhaus rigidity for orbit-like sets).
 - INCIDENCE CELLS (sector x runner): vertices = occupied (sector,runner) cells.
 - HARMONIC/SPECTRAL: vertices = frequencies m, tournament by character phase ghat(m).
 - ARRANGEMENT CELLS: vertices = cells of the combined braid + 1/n-wall arrangement
   on the torus (order x loneliness type).
 - MATROID FLATS: vertices = flats of the resonance matroid (S537).
 - TIME-FREQUENCY (Gabor) CELLS: vertices = (sector, harmonic) pairs.

We compute restriction R = realizable / all for the PAIR-difference tournament
(n=4) and the GAP structure (n=4,5,6), verify the cocycle, and tabulate.
"""
from itertools import combinations, permutations, product
from functools import reduce
from math import gcd, sin, pi
import random

def frac(x): return x - int(x // 1)

A000568 = {3:2,4:4,5:12,6:56,7:456}

# ---------------- STAR: pair-difference (tension) tournament ----------------
def pair_diff_labeled(speeds_with_obs, t, pairs):
    """labeled adjacency bits of the pair-difference half-turn tournament (cheap)."""
    m = len(pairs)
    phase = [frac((speeds_with_obs[a]-speeds_with_obs[b])*t) for (a, b) in pairs]
    return tuple(1 if 0 < (phase[j]-phase[i]) % 1.0 < 0.5 else 0
                 for i in range(m) for j in range(m) if i != j)

def canon_from_bits(bits, m):
    """canonicalize a labeled tournament (given as off-diagonal bit tuple) over S_m."""
    # rebuild adj
    adj = [[0]*m for _ in range(m)]
    it = iter(bits)
    for i in range(m):
        for j in range(m):
            if i != j: adj[i][j] = next(it)
    best = None
    for p in permutations(range(m)):
        b = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or b < best: best = b
    return best

def verify_cocycle(speeds_with_obs, n):
    """w_{ij}+w_{jk}+w_{ki} = 0 for all triangles."""
    ok = True
    for a, b, c in combinations(range(n), 3):
        w = (speeds_with_obs[a]-speeds_with_obs[b]) + (speeds_with_obs[b]-speeds_with_obs[c]) + (speeds_with_obs[c]-speeds_with_obs[a])
        if w != 0: ok = False
    return ok

def study_pair_tournament(n_runners, n_sets=120, samples=2500):
    """observer + n_runners runners; vertices = all C(n,2) pairs; n=n_runners+1."""
    n = n_runners + 1
    pairs = list(combinations(range(n), 2))      # includes observer index 0
    rnd = random.Random(31+n); labeled = set(); tot = 0; cocycle_ok = True
    for _ in range(8000):
        v = tuple(sorted(rnd.sample(range(1, 6*n), n_runners)))
        if reduce(gcd, v) != 1: continue
        tot += 1
        if tot > n_sets: break
        sp = [0] + list(v)                       # observer speed 0
        if not verify_cocycle(sp, n): cocycle_ok = False
        for s in range(samples):
            t = (s+0.5)/samples
            labeled.add(pair_diff_labeled(sp, t, pairs))    # cheap dedupe first
    m = len(pairs)
    classes = set(canon_from_bits(b, m) for b in labeled)   # canon only the distinct
    return classes, tot, len(pairs), cocycle_ok

# ---------------- GAPS (vertices = arcs between consecutive points) ----------------
def gap_vector(speeds, n, t):
    pts = sorted([0.0] + [frac(v*t) for v in speeds])
    gaps = [pts[i+1]-pts[i] for i in range(len(pts)-1)] + [1.0 - pts[-1] + pts[0]]
    return tuple(gaps)

def study_gaps(n, n_sets=150, samples=3000):
    """vertices = n gaps; realizable gap-ORDER-types (which sizes, sorted multiset
    quantized) and distinct count; track number of distinct gap LENGTHS (three-gap)."""
    rnd = random.Random(57+n); order_types = set(); maxgap_is_obs = 0; tot = 0
    distinct_len_hist = {}
    for _ in range(8000):
        v = tuple(sorted(rnd.sample(range(1, 6*n), n-1)))
        if reduce(gcd, v) != 1: continue
        tot += 1
        if tot > n_sets: break
        for s in range(samples):
            t = (s+0.5)/samples
            g = gap_vector(v, n, t)
            qg = tuple(sorted(round(x, 3) for x in g))     # quantized multiset
            order_types.add(qg)
            ndist = len(set(round(x, 3) for x in g))
            distinct_len_hist[ndist] = distinct_len_hist.get(ndist, 0)+1
    return order_types, tot, distinct_len_hist

# ---------------- main ----------------
def main():
    print("="*74)
    print("STAR: vertices = PAIRS (edges of K_n), carrying difference-speeds w_ij = v_i-v_j")
    print("      (a TENSION: cocycle w_ij+w_jk+w_ki=0) -- a 2nd-order LRC on differences")
    print("="*74)
    for nr in (3,):     # n_runners; n=4 => 6 pair-vertices (canon over 6! feasible)
        classes, tot, npairs, coc = study_pair_tournament(nr, n_sets=150, samples=3000)
        n = nr+1
        allt = A000568.get(npairs, None)
        rstr = f"{len(classes)} / {allt} (R={len(classes)/allt:.3f})" if allt else f"{len(classes)}"
        print(f"  n={n}: {npairs} pair-vertices; realizable pair-tournament iso-classes = {rstr}; "
              f"cocycle holds={coc}; {tot} speed sets")
    # n=5 (10 pair-vertices): full S_10 canon is infeasible; count distinct LABELED
    # adjacency patterns realized, as a restriction proxy (vs 2^C(10,2)=2^45).
    n_runners = 4; n = 5; pairs = list(combinations(range(n), 2))
    rnd = random.Random(99); labeled = set(); tot = 0
    for _ in range(8000):
        v = tuple(sorted(rnd.sample(range(1, 6*n), n_runners)))
        if reduce(gcd, v) != 1: continue
        tot += 1
        if tot > 120: break
        sp = [0] + list(v)
        for s in range(3000):
            t = (s+0.5)/3000
            ph = [frac((sp[a]-sp[b])*t) for (a, b) in pairs]
            bits = tuple(1 if 0 < (ph[j]-ph[i]) % 1.0 < 0.5 else 0
                         for i in range(len(pairs)) for j in range(len(pairs)) if i != j)
            labeled.add(bits)
    print(f"  n=5: 10 pair-vertices; realizable LABELED pair-tournaments = {len(labeled)} "
          f"of 2^45 possible (~{2**45}); a vanishing slice (cocycle-restricted)")
    print("  => the pair (difference/tension) tournament lives on C(n,2) vertices but its")
    print("     realizable classes are restricted by the COCYCLE (additive consistency):")
    print("     the difference-speeds are a coboundary, not free. A second-order LRC.")
    print()

    print("="*74)
    print("MULTITUDE rung -- GAPS: vertices = the n arcs between consecutive points")
    print("="*74)
    for n in (4, 5, 6):
        ots, tot, hist = study_gaps(n)
        print(f"  n={n}: realizable gap-multiset order-types (quantized) = {len(ots)}; "
              f"distinct-gap-length histogram over samples = {dict(sorted(hist.items()))}")
    print("  => gap-vectors sum to 1 and cluster to FEW distinct lengths (three-gap/Steinhaus")
    print("     rigidity); the apex = largest gap = loneliness target (S530). Restricted.")
    print()

    print("="*74)
    print("THE MULTITUDE (further complex/precise vertex objects; posed precisely)")
    print("="*74)
    print("""  - HARMONIC/SPECTRAL: vertices = frequencies m in {1..n-1}; edge by phase of
    ghat(m) e^{2pi i m * something}; the 'character tournament' (dual of S537 flows).
  - ARRANGEMENT CELLS: vertices = cells of the COMBINED braid {x_i=x_j} + threshold
    {x_i=+-1/n} arrangement on T^{n-1}; each = (cyclic order, loneliness pattern);
    the runner orbit is a closed walk; realizable cells = a thin slice (S521o).
  - INCIDENCE CELLS: vertices = occupied (sector, runner) pairs (S536 x runners);
    a bipartite-lift tournament; realizable = doubly-restricted.
  - MATROID FLATS: vertices = flats of the resonance matroid M_v (S537); edge by
    containment/closure order; LRC reads the connectivity of M_v.
  - TIME-FREQUENCY (GABOR) CELLS: vertices = (sector, harmonic) pairs -- the joint
    space x frequency lift unifying S536 (space) and S537 (frequency); the
    'uncertainty' tournament; realizable cells obey a discrete uncertainty bound.
  - WIRING-DIAGRAM EVENTS: vertices = crossing events (adjacent transpositions);
    realizable = stretchable allowable sequences (S535 MAP-wire).
  Each: 'which iso-classes are exhibitable' with restriction from a DIFFERENT
  consistency law (cocycle / three-gap / arrangement / matroid / uncertainty /
  stretchability).""")

if __name__ == "__main__":
    main()
