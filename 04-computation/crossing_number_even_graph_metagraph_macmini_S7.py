#!/usr/bin/env python3
"""
CROSSING NUMBER of K_n  vs  the EVEN-GRAPH metagraph E_n and OCF/H tournament data.
mac-mini-2026-06-20-S7   THREAD: Guy's conjecture <-> even graphs / tournaments.

Guy's conjecture (= cylindrical = 2-page crossing number, proved as UPPER bound;
Abrego et al proved nu_2(K_n) = Z(n)):
    Z(n) = (1/4) floor(n/2) floor((n-1)/2) floor((n-2)/2) floor((n-3)/2).

Repo objects probed:
  - Convex (1-page) drawing of K_n: cr_convex = C(n,4)  (one crossing per 4-set).
  - 2-page book drawing: edges (= CHORDS of the n-gon) get a GF(2) page label.
    Two SAME-PAGE edges cross iff their endpoints interleave around the spine.
    This is EXACTLY a GF(2) 2-coloring of the chord set -> the tiling-hypercube /
    cut+cycle structure that carries tournaments (orientation of K_n) and even
    graphs E_n in this repo.
  - OCF / c3:  Kendall-Babington-Smith  c3 = C(n,3) - sum_v C(s_v,2).
  - THM-555  E_tiling[c3] = (C(n,3)+(n-2))/4.

HYPOTHESES TESTED (all FALSIFIABLE):

H1 (structural identity, the real claim):
    The optimal 2-page crossing number nu_2(K_n) = Guy's Z(n), and the OPTIMAL
    page-bipartition of the chords is a SPECIFIC GF(2) vector in the chord space.
    Claim: Z(n) is reproduced by the cut/cycle GF(2) optimization over chord
    2-colorings, and the minimizing 2-coloring is the "balanced alternating"
    bipartition. TEST: brute-force optimal 2-page over all 2-colorings for small n,
    confirm = Z(n), and identify the optimizer.

H2 (even-graph correspondence):
    A 2-coloring of the edges of K_n that is "even on each page" (every vertex has
    even page-degree on each page) <-> a pair (even graph, complementary even graph).
    The page-degree-even colorings are exactly the EVEN GRAPHS used to build E_n.
    Claim: among 2-page drawings, the EVEN ones (both pages even graphs) have a
    crossing count whose MINIMUM equals Z(n) (i.e. the optimum can be taken even).
    TEST: restrict the brute force to even-page colorings; is min still Z(n)?

H3 (c3 / OCF numeric link -- the "surface analogy" guard):
    Is Z(n) any clean function of E_tiling[c3]=(C(n,3)+(n-2))/4, of C(n,4), or of H?
    TEST: ratios Z(n)/C(n,4), Z(n)/E[c3], differences; look for exact identity.
"""

from fractions import Fraction
from itertools import combinations
from math import comb
import sys

def guy(n):
    return (n//2)*((n-1)//2)*((n-2)//2)*((n-3)//2)//4

# ---------------------------------------------------------------------------
# 2-PAGE (book) crossing number.
# Vertices 0..n-1 on a spine (a line). Each edge {i,j} is drawn in page 0 (above)
# or page 1 (below). Two edges on the SAME page cross iff their endpoints
# interleave: i<k<j<l (strict alternation around the spine).
# nu_2(K_n) = min over page assignments of (# same-page interleaving pairs).
# ---------------------------------------------------------------------------

def interleave(e, f):
    (a, b) = e; (c, d) = f
    # assume a<b, c<d. They cross on a line iff a<c<b<d or c<a<d<b.
    return (a < c < b < d) or (c < a < d < b)

def crossing_pairs(n):
    """All pairs of edges that interleave (independent of page)."""
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    cps = []
    for x in range(len(edges)):
        for y in range(x+1, len(edges)):
            e, f = edges[x], edges[y]
            # only count edge pairs that share no vertex AND interleave
            if len(set(e) | set(f)) == 4 and interleave(e, f):
                cps.append((x, y))
    return edges, cps

def two_page_crossings(n, page):
    """page: dict edge_index -> 0/1. Count same-page interleaving pairs."""
    edges, cps = crossing_pairs(n)
    c = 0
    for (x, y) in cps:
        if page[x] == page[y]:
            c += 1
    return c

def brute_two_page(n, even_only=False):
    """Brute force optimal 2-page crossing number.
       even_only: restrict to assignments where each page is an EVEN graph
       (every vertex even degree on each page)."""
    edges, cps = crossing_pairs(n)
    m = len(edges)
    best = None
    best_assign = None
    # vertex incidence
    inc = [[] for _ in range(n)]
    for idx, (i, j) in enumerate(edges):
        inc[i].append(idx); inc[j].append(idx)
    total = 1 << m
    if total > 2**22:
        return None, None  # too big to brute force
    for mask in range(total):
        if even_only:
            ok = True
            for v in range(n):
                deg0 = sum(1 for idx in inc[v] if not (mask >> idx) & 1)
                if deg0 % 2 != 0:  # page-0 degree even => page-1 degree also even (total deg n-1)
                    # only consistent when n-1 even, i.e. n odd; for n even degrees can't both be even
                    ok = False; break
            if not ok:
                continue
        c = 0
        for (x, y) in cps:
            if ((mask >> x) & 1) == ((mask >> y) & 1):
                c += 1
        if best is None or c < best:
            best = c; best_assign = mask
            if best == guy(n):
                # keep searching only to confirm; we can early accept the value
                pass
    return best, best_assign

# ---------------------------------------------------------------------------
# OCF / c3 / convex data
# ---------------------------------------------------------------------------

def report():
    print("="*78)
    print("Crossing-number vs even-graph / OCF data")
    print("="*78)
    print(f"{'n':>2} {'Guy Z(n)':>9} {'C(n,4)':>8} {'E[c3]':>10} {'Z/C(n,4)':>10} "
          f"{'2page brute':>11} {'2page=Z?':>9} {'evenpage min':>12}")
    for n in range(3, 11):
        Z = guy(n)
        c4 = comb(n, 4)
        Ec3 = Fraction(comb(n,3) + (n-2), 4)
        ratio = Fraction(Z, c4) if c4 else Fraction(0)
        bp = ""
        eq = ""
        evp = ""
        edges = [(i,j) for i in range(n) for j in range(i+1,n)]
        if len(edges) <= 18:
            best, _ = brute_two_page(n, even_only=False)
            bp = str(best)
            eq = "YES" if best == Z else f"NO({best} vs {Z})"
            be, _ = brute_two_page(n, even_only=True)
            evp = str(be) if be is not None else "(none/NA)"
        print(f"{n:>2} {Z:>9} {c4:>8} {str(Ec3):>10} {str(ratio):>10} "
              f"{bp:>11} {eq:>9} {evp:>12}")

    print()
    print("-"*78)
    print("H1 detail: is the optimal 2-page bipartition the 'balanced alternating' one?")
    print("-"*78)
    for n in range(4, 9):
        edges = [(i,j) for i in range(n) for j in range(i+1,n)]
        if len(edges) > 18:
            print(f"  n={n}: too large to brute force ({len(edges)} edges)")
            continue
        best, assign = brute_two_page(n, even_only=False)
        # describe the optimizer: page label by edge "length" j-i (mod n) parity?
        page_by_len = {}
        consistent_len = True
        for idx,(i,j) in enumerate(edges):
            L = j - i
            p = (assign >> idx) & 1
            if L in page_by_len and page_by_len[L] != p:
                consistent_len = False
            page_by_len[L] = p
        print(f"  n={n}: nu_2={best} (Z={guy(n)})  optimizer page-determined-by-edge-length? {consistent_len}")
        if consistent_len:
            print(f"         page(L) for L=1..{n-1}: " +
                  " ".join(f"{L}->{page_by_len[L]}" for L in range(1, n)))

if __name__ == "__main__":
    report()
