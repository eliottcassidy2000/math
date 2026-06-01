#!/usr/bin/env python3
"""
recursive_apex_hierarchy_s531.py    oracle-2026-06-01-S531

HOW THE APEX (S530) APPLIES RECURSIVELY.

S530: the polygon OUTSIDE = base path (n-1) + the apex/source-sink arc (n,1);
the apex is the unique boundary tile, the cut/cycle hinge (its fundamental cycle
is the whole n-gon), and the LRC observer's loneliness gap.

THE RECURSIVE SEED: every tile (x,y) (range x-y >= 2) is the apex of the
SUB-POLYGON on the contiguous rank-block [y, x]. So the outside/inside split
recurs at every scale. This script develops that recursion across THREE
quantities and looks for self-similar laws:

  A. COMBINATORICS. tiles  <->  sub-intervals [y,x] of length >= 3  <->
     sub-polygon apexes. The nesting is laminar; the recursive sub-ranking tree
     (S520o) IS the polygon's recursive triangulation. Self-similar counts.

  B. THE HAMILTONIAN-PATH COUNT H (loneliness meter, S26). Flipping a single
     tile of range r on the transitive tournament gives H = 1 + 2^{r-1}
     (apex r=n-1 -> 1+2^{n-2}, parallel-S530). We test whether flipping a
     LAMINAR (nested/disjoint) family of apexes obeys a RECURSIVE law, building
     H of a block from H of its sub-blocks.

  C. LRC LONELINESS, recursively. The apex = observer's gap = the surviving
     feasible arc. Processing runners (the cascade, S527) carves the feasible
     set; the surviving arc is a NESTED sequence of sub-apexes. We track the
     recursive carve tree and its shrink factor (~ (n-2)/n, S527) and tie the
     carve tree to the combinatorial sub-polygon tree.
"""
from itertools import combinations, permutations
from functools import lru_cache
from math import gcd

# ======================================================================
# A. COMBINATORICS: tiles = sub-polygon apexes; the recursive hierarchy
# ======================================================================
def tiles(n):
    return [(x, y) for y in range(1, n+1) for x in range(y+2, n+1)]

def subpolygon_of_tile(x, y):
    """tile (x,y) is the apex of the sub-polygon on the contiguous block [y..x]."""
    return tuple(range(y, x+1))

def laminar_check(intervals):
    """are the intervals pairwise nested or disjoint (laminar)?"""
    for (a, b), (c, d) in combinations(intervals, 2):
        # intervals [a,b],[c,d]; overlap-but-not-nested = crossing
        nested = (a <= c and d <= b) or (c <= a and b <= d)
        disjoint = b < c or d < a
        if not (nested or disjoint):
            return False
    return True

def combinatorics(n):
    tl = tiles(n)
    # self-similar census: a size-s sub-polygon has apex(1)+subpath(s-1)+inside(s(s-3)/2)
    by_size = {}
    for (x, y) in tl:
        s = x - y + 1                       # #vertices in the sub-polygon
        by_size[s] = by_size.get(s, 0) + 1  # = n - s + 1 blocks of length s
    return tl, by_size

# ======================================================================
# B. H = directed-Hamiltonian-path count, and the recursive apex law
# ======================================================================
def transitive_adj(n):
    """vertices 1..n, transitive i->j iff i>j. base path n->...->1."""
    beats = [[False]*(n+1) for _ in range(n+1)]
    for i in range(1, n+1):
        for j in range(1, n+1):
            if i > j:
                beats[i][j] = True
    return beats

def flip_tiles(beats, flips, n):
    """flip the orientation of each tile (x,y) in flips (x>y): make y->x instead."""
    b = [row[:] for row in beats]
    for (x, y) in flips:
        b[x][y] = False
        b[y][x] = True
    return b

def count_ham_paths(beats, n):
    """directed Hamiltonian path count via subset DP. O(2^n n^2)."""
    full = (1 << n) - 1
    # dp[mask][v]: # paths covering exactly the vertices in mask, ending at v (0-indexed v=vertex-1)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            cur = dp[mask][v]
            if not cur:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if beats[v+1][u+1]:        # arc v -> u
                    dp[mask | (1 << u)][u] += cur
    return sum(dp[full][v] for v in range(n))

def H_after_flips(n, flips):
    return count_ham_paths(flip_tiles(transitive_adj(n), flips, n), n)

# ======================================================================
# C. LRC: recursive carving of the feasible (loneliness) arc set
# ======================================================================
def feasible_intervals_after_runners(speeds, n, samples=4000):
    """Approximate, by sampling t in [0,1): the maximal sub-arcs of t where ALL
    processed runners are 'safe' (||v t|| >= 1/n). Returns, per prefix of runners,
    the number of surviving feasible components and their total measure -- the
    recursive carve."""
    thr = 1.0/n
    out = []
    safe_prev = [True]*samples
    for k, v in enumerate(speeds):
        safe = []
        for i in range(samples):
            t = (i + 0.5)/samples
            x = (v*t) % 1.0
            ok = thr <= min(x, 1-x)
            safe.append(safe_prev[i] and ok)
        # count components + measure
        comps = 0; meas = 0; inrun = False
        for s in safe:
            if s:
                meas += 1
                if not inrun: comps += 1; inrun = True
            else:
                inrun = False
        out.append((k+1, comps, meas/samples))
        safe_prev = safe
    return out

# ======================================================================
def main():
    print("="*74)
    print("A. RECURSIVE HIERARCHY: tiles = sub-polygon apexes (self-similar census)")
    print("="*74)
    for n in range(3, 9):
        tl, by_size = combinatorics(n)
        census = "  ".join(f"size{ s}:{c}(=n-s+1={n-s+1})" for s, c in sorted(by_size.items()))
        total = len(tl)
        print(f"  n={n}: {total} tiles = sub-polygons by size:  {census}")
    print("  Each size-s sub-polygon recurses: apex(1) + sub-base-path(s-1) + inside(s(s-3)/2).")
    print("  Laminar nesting of these intervals = the recursive sub-ranking tree (S520o)")
    print("  = a recursive triangulation of the n-gon. Example laminar family n=7:")
    fam = [(7,1),(5,1),(3,1)]   # nested apexes
    ivs = [(y,x) for (x,y) in fam]
    print(f"     flips {fam} -> intervals {[subpolygon_of_tile(x,y) for (x,y) in fam]}, "
          f"laminar={laminar_check(ivs)}")
    print()

    print("="*74)
    print("B. RECURSIVE H: single-tile law  H = 1 + 2^(r-1)  (r = range), then nesting")
    print("="*74)
    print("  Single tile of range r flipped on transitive (verify 1+2^(r-1)):")
    for n in [8]:
        for r in range(2, n):
            x, y = (1 + r, 1)      # tile of range r
            H = H_after_flips(n, [(x, y)])
            print(f"     n={n} flip tile range r={r} (={(x,y)}): H={H:5d}  1+2^(r-1)={1+2**(r-1)}  "
                  f"match={H==1+2**(r-1)}")
    print()
    print("  NESTED apexes (laminar chain), seeking a recursive product/sum law.")
    print("  Notation: flips given as intervals [y..x]; H reported.")
    experiments = [
        (5, [(5,1)]),
        (5, [(5,1),(3,1)]),
        (5, [(5,1),(4,2)]),
        (5, [(5,1),(3,1),(4,2)]),
        (6, [(6,1)]),
        (6, [(6,1),(4,1)]),
        (6, [(6,1),(4,1),(3,1)]),
        (6, [(6,1),(5,2)]),
        (6, [(6,1),(5,2),(4,3)]),
        (7, [(7,1)]),
        (7, [(7,1),(5,1)]),
        (7, [(7,1),(5,1),(3,1)]),
        (7, [(7,1),(6,2),(5,3)]),
    ]
    for n, flips in experiments:
        ivs = [(y, x) for (x, y) in flips]
        H = H_after_flips(n, flips)
        ranges = sorted((x-y) for (x, y) in flips)
        prod_single = 1
        for r in ranges:
            prod_single *= (1 + 2**(r-1))
        sum_single = 1 + sum(2**(r-1) for r in ranges)
        print(f"  n={n} flips(ivs)={[ (y,x) for (x,y) in flips]:} lam={laminar_check(ivs)}  "
              f"H={H}   [Prod(1+2^(r-1))={prod_single}, 1+Sum2^(r-1)={sum_single}]")
    print()

    print("="*74)
    print("C. LRC RECURSIVE CARVING: feasible (loneliness) arcs vs runner prefix")
    print("="*74)
    cases = {5: [(1,2,3,4),(1,2,3,5)], 7: [(1,2,3,4,5,6),(1,2,3,4,5,7)]}
    for n, sets in cases.items():
        for v in sets:
            ap = " [AP=regular polygon]" if list(v)==list(range(1,n)) else ""
            carve = feasible_intervals_after_runners(v, n)
            chain = "  ".join(f"k{k}:{c}c/{m:.3f}" for (k,c,m) in carve)
            print(f"  n={n} v={v}{ap}")
            print(f"     carve (after each runner: #components / measure): {chain}")
            # recursive shrink ratios
            meas = [m for (_,_,m) in carve]
            ratios = [meas[i+1]/meas[i] if meas[i]>0 else float('nan') for i in range(len(meas)-1)]
            print(f"     shrink ratios: {[round(r,3) for r in ratios]}   (S527 predicts ~ (n-2)/n = {(n-2)/n:.3f})")
    print()
    print("  Each runner SPLITS the surviving feasible arc into sub-arcs: the carve")
    print("  tree mirrors the recursive sub-polygon tree of part A. The observer's")
    print("  apex (largest surviving arc) is a NESTED sequence of sub-apexes; LRC =")
    print("  'the recursive carving never empties the apex' (AP = tight: empties to 0).")

if __name__ == "__main__":
    main()
