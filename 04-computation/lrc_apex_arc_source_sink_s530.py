#!/usr/bin/env python3
"""
lrc_apex_arc_source_sink_s530.py    oracle-2026-06-01-S530

REFINEMENT of S529: the polygon OUTSIDE = base path + the ONE source-sink arc.

S529 said "outside (polygon sides) = base path." That is incomplete. The polygon
boundary is a Hamiltonian CYCLE (n arcs); the tiling-model base path is only n-1
arcs (n -> n-1 -> ... -> 1). The missing arc -- the one that CLOSES the path into
the cycle -- is the arc between the SOURCE (vertex n, top of ranking) and the SINK
(vertex 1, bottom). This script proves that this closing arc:

  (i)  is, in tiling combinatorics, a TILE (range x-y = n-1), NOT a base-path arc
       -- in fact the MAXIMAL-range tile = the APEX of the staircase triangle;
  (ii) is the UNIQUE tile that lies on the polygon boundary (cyclic skip 1); every
       OTHER tile is an interior diagonal. So:
           OUTSIDE (n sides) = base path (n-1) + apex tile (1);
           INSIDE (diagonals) = all OTHER tiles = C(n-1,2) - 1 = n(n-3)/2.
  (iii) has fundamental cycle = the WHOLE polygon (the n-cycle); every other tile
        (x,y) has fundamental cycle = the contiguous (x-y+1)-sub-polygon. So the
        apex is the cut/cycle HINGE: a cycle-space generator that is geometrically
        a boundary (tree-like) edge.
  (iv) reconciles S529's "inside born at n=4": for n=3 the ONLY tile IS the apex
        (C(n-1,2)=1), which is on the boundary -> no interior -> no inside debt.

LRC APPLICATION
---------------
Root the cyclic order of the n points (observer at 0 + the n-1 runners) at the
observer. The source = leading point, sink = trailing point; the SOURCE-SINK ARC
is the largest "wrap-around" gap -- the arc you cut to linearize the cycle into the
ranking path. The OBSERVER lives inside this apex arc. Observer-loneliness
(observer >= 1/n from every runner) <=> the apex arc that brackets the observer is
wide enough (each side >= 1/n). So:

   The apex (source-sink) arc IS the observer's loneliness gap.
   LRC = "the dynamically-changing apex arc opens to clearance 1/n on both sides."

The n arcs around the circle (n-1 inter-runner gaps = the runner base-path, plus
the 1 observer-bracketing apex) sum to 1. Loneliness is decided entirely by the
apex. At the regular polygon every gap = 1/n, so the apex sits EXACTLY at the
critical clearance -- the tight case, again.

Verifications below.
"""
from math import gcd
from fractions import Fraction
from itertools import combinations

# ----------------------------------------------------------------------
# (i)-(iv) The combinatorial dictionary
# ----------------------------------------------------------------------
def tiles(n):
    """tiles = arcs (x,y), x>y, x-y>=2, vertices 1..n. range = x-y."""
    return [(x, y) for y in range(1, n+1) for x in range(y+2, n+1)]

def base_path(n):
    """consecutive-rank arcs (k,k-1), k=2..n. n-1 of them."""
    return [(k, k-1) for k in range(2, n+1)]

def cyclic_skip(x, y, n):
    """vertices placed on circle in ranking order (vertex k at position k-1);
    chord skip = min(rank-gap, n - rank-gap)."""
    d = abs(x - y)
    return min(d, n - d)

def fundamental_cycle_len(x, y):
    """adding non-tree edge (x,y) to the base path makes the contiguous cycle
    x -> x-1 -> ... -> y -> x, of length (x-y+1)."""
    return (x - y) + 1

def check_dictionary(n):
    bp = base_path(n)
    tl = tiles(n)
    apex = (n, 1)
    # boundary sides = arcs with cyclic skip 1
    bp_skips = {cyclic_skip(x, y, n) for (x, y) in bp}
    apex_skip = cyclic_skip(*apex, n)
    inside = [(x, y) for (x, y) in tl if cyclic_skip(x, y, n) >= 2]
    boundary_tiles = [(x, y) for (x, y) in tl if cyclic_skip(x, y, n) == 1]
    return {
        "n": n,
        "n_basepath": len(bp),
        "basepath_all_skip1": bp_skips == {1},
        "apex": apex,
        "apex_range": n - 1,
        "apex_is_tile": apex in tl,
        "apex_skip": apex_skip,
        "boundary_tiles": boundary_tiles,           # should be exactly [apex]
        "n_tiles": len(tl),
        "n_inside": len(inside),
        "predicted_inside_n(n-3)/2": n*(n-3)//2,
        "boundary_total(n-1)+1": len(bp) + 1,
        "apex_fund_cycle_len": fundamental_cycle_len(*apex),  # = n (whole polygon)
    }

# ----------------------------------------------------------------------
# LRC: the apex arc = observer's bracketing gap; gaps sum to 1; loneliness
# ----------------------------------------------------------------------
def frac(x):
    return x - int(x // 1)

def observer_apex_gap(speeds, t):
    """The observer is at 0. Among runners, nearest ahead a and nearest behind b.
    Apex arc (bracketing observer) = [-b, a], width a+b. Returns (a, b, width)."""
    xs = [frac(s * t) for s in speeds]
    ahead = [x for x in xs if x > 0]
    a = min(ahead) if ahead else 1.0
    b = min((1.0 - x) for x in xs) if xs else 1.0
    return a, b, a + b

def is_lonely(speeds, t, n):
    a, b, _ = observer_apex_gap(speeds, t)
    return a >= 1.0/n - 1e-12 and b >= 1.0/n - 1e-12

def gaps_sum_check(speeds, t):
    """all n arcs (observer + runners) around the circle sum to 1."""
    pts = sorted([0.0] + [frac(s*t) for s in speeds])
    g = [pts[i+1]-pts[i] for i in range(len(pts)-1)] + [1.0 - pts[-1] + pts[0]]
    return sum(g), len(g)

def scan_lonely(speeds, n, G=200000):
    """find a lonely time and report the apex gap there; also max apex clearance."""
    best_clear = -1.0; best_t = None
    found = None
    for i in range(G):
        t = (i + 0.5)/G
        a, b, w = observer_apex_gap(speeds, t)
        clear = min(a, b)
        if clear > best_clear:
            best_clear = clear; best_t = t
        if found is None and clear >= 1.0/n - 1e-9:
            found = (t, a, b, w)
    return found, best_clear, best_t

# ----------------------------------------------------------------------
def main():
    print("="*74)
    print("(i)-(iv) THE DICTIONARY: outside = base path + apex(source-sink) tile")
    print("="*74)
    for n in range(3, 9):
        d = check_dictionary(n)
        ok_inside = d["n_inside"] == d["predicted_inside_n(n-3)/2"]
        ok_boundary = d["boundary_total(n-1)+1"] == n
        print(f"  n={n}:")
        print(f"     base path: {d['n_basepath']} arcs, all cyclic-skip 1: {d['basepath_all_skip1']}")
        print(f"     apex {d['apex']}: range x-y = {d['apex_range']} (=n-1, MAX tile), "
              f"is a tile: {d['apex_is_tile']}, cyclic skip = {d['apex_skip']} (boundary side)")
        print(f"     tiles on boundary (skip 1): {d['boundary_tiles']}  <- should be just the apex")
        print(f"     #tiles={d['n_tiles']}, #inside-diagonals={d['n_inside']} "
              f"(predicted n(n-3)/2={d['predicted_inside_n(n-3)/2']}: {ok_inside})")
        print(f"     boundary sides = (n-1)+1 = {d['boundary_total(n-1)+1']} = n: {ok_boundary}; "
              f"apex fundamental cycle length = {d['apex_fund_cycle_len']} (= whole n-gon)")
    print()
    print("  => n=3: the ONLY tile is the apex (3,1); it is on the boundary, so the")
    print("     triangle has NO interior diagonal -> no inside debt (matches S529).")
    print("     n=4: 3 tiles = apex (4,1) + 2 interior diagonals (3,1),(4,2) -> inside born.")
    print()

    print("="*74)
    print("LRC: the apex arc = observer's bracketing gap; gaps sum to 1")
    print("="*74)
    cases = {3: [(1,2),(1,3)], 4: [(1,2,3),(1,3,5)], 5: [(1,2,3,4),(1,2,3,5)]}
    for n, sets in cases.items():
        for v in sets:
            # gap-sum sanity at a generic time
            s, k = gaps_sum_check(v, 0.1234567)
            found, best_clear, best_t = scan_lonely(v, n)
            ap = " [AP=regular polygon]" if list(v)==list(range(1,n)) else ""
            print(f"  n={n} v={v}{ap}: {k} circle-arcs sum to {s:.6f}; "
                  f"max apex clearance min(a,b) = {best_clear:.5f} (threshold 1/n={1.0/n:.5f})")
            if found:
                t,a,b,w = found
                print(f"        LONELY at t={t:.5f}: nearest ahead a={a:.5f}, behind b={b:.5f}, "
                      f"apex width a+b={w:.5f} (>=2/n={2.0/n:.5f})")
            else:
                print(f"        max clearance {best_clear:.5f} at t={best_t:.5f} "
                      f"-> tight/boundary (regular polygon reaches exactly 1/n)")
    print()
    print("  => Loneliness is decided ENTIRELY by the apex (source-sink) arc: the")
    print("     observer needs its bracketing gap to clear 1/n on both sides. The base")
    print("     path = the n-1 inter-runner gaps; the apex = the observer's gap; they")
    print("     sum to 1. At the regular polygon every gap = 1/n, so the apex sits")
    print("     exactly at the critical clearance -- the tight case once more.")

if __name__ == "__main__":
    main()
