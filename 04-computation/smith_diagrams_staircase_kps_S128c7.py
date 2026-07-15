#!/usr/bin/env python3
"""
smith_diagrams_staircase_kps_S128c7.py
======================================
kind-pasteur-2026-07-15-S128 (cont.7).  SMITH DIAGRAMS x THE STAIRCASE (owner: squared squares).

Brooks-Smith-Stone-Tutte read a squared rectangle as a planar electrical network: nodes = maximal
horizontal segments, edges = squares, currents = sizes.  The project's staircase Delta_{n-2} IS a
dissection domain (unit cells = tiles = tournament arc slots), so it carries a canonical Smith
network N_n: nodes = the maximal horizontal segments of the staircase region, edges = cells.
This script computes, exactly:

 (1) N_n for n=4..10: |V|, |E|=C(n-1,2), cycle rank; kappa(N_n) = spanning trees (integer
     Kirchhoff); the POLE RESISTANCE R_n (top segment -> bottom segment) as an exact Fraction --
     the staircase's "electrical aspect ratio" -- and the BSST current vector at the poles
     (would-be square sizes; the dissection's Smith currents).
 (2) SELF-DUALITY test: the anti-diagonal reflection swaps horizontal/vertical segments; compare
     N_n against the vertical-segment network (the BSST dual): isomorphic as two-pole networks?
     (Numerically: same kappa, same R, degree sequences.)
 (3) THE BATTERY TABLE (n=5,6): per merged class K, flux Phi(K) = Sum_{fiber t} (e1(t) - e_n(t))
     (opus THM-790 leg law: Delta E4 = 8(e1-e_n) per line).  Checks: global sum = 0; the BLUE
     identity e1 + e_n = n-2 on every grid-sym tiling (one-line derivation from rho o op); Phi by
     type (PB/MX/PBk) -- where do the metagraph's batteries live?
 (4) Bouwkamp-style fiber code: canonical class rep = lexicographically least tiling word in the
     fiber (staircase-native; explorer-compatible); how it interacts with the cont.5 ordering.
"""
import sys
from fractions import Fraction as F
from math import comb
from itertools import permutations, combinations
from collections import Counter, defaultdict

sys.stdout.reconfigure(line_buffering=True)

def staircase_cells(n):
    """cells (x,y) of Delta_{n-2}: row y (1..n-2), x from y+2..n. Geometric placement:
    put cell (x,y) at column c = x-y-1 (1..n-1-y), row r = y (rows stack upward: row y at height y).
    This is the explorer's staircase: row y has n-1-y cells, right angle at bottom-left."""
    return [(x, y) for y in range(1, n - 1) for x in range(y + 2, n + 1)]

def smith_network(n, horizontal=True):
    """Nodes = maximal horizontal segments of the staircase region (or vertical if not horizontal).
    Cell (x,y) at grid position col=x-y-1, row=y. Horizontal segment at height h spans the
    union of cell-tops (height y+... define: cell (col,row) occupies [col-1,col]x[row-1,row].
    Bottom boundary of row y cells at height y-1; top at height y.
    Returns (adjacency edges list mapping cell->(nodeA,nodeB), #nodes, pole_top, pole_bottom)."""
    cells = staircase_cells(n)
    pos = {}
    for (x, y) in cells:
        col = x - y - 1
        row = y
        pos[(x, y)] = (col, row)
    # horizontal lines: at integer heights h=0..n-2; at height h, segments = maximal runs of
    # columns covered by (cell tops at h) union (cell bottoms at h)
    if horizontal:
        # line at height h: columns where a cell has top=h (row h) or bottom=h (row h+1)
        seg_id = {}
        nid = 0
        for h in range(0, n - 1):
            cols = set()
            for (x, y), (c, r) in pos.items():
                if r == h or r == h + 1:      # top of row h cells is height h ... careful:
                    pass
            # top of cell in row r is height r; bottom is r-1
            cols = sorted({c for (x, y), (c, r) in pos.items() if r == h + 1} |
                          {c for (x, y), (c, r) in pos.items() if r == h})
            # actually: line at height h touches cells with top=h (r=h) or bottom=h (r=h+1)
            cols = sorted({c for (x, y), (c, r) in pos.items() if r == h or r == h + 1})
            if not cols:
                continue
            runs = []
            start = prev = cols[0]
            for c in cols[1:]:
                if c == prev + 1:
                    prev = c
                else:
                    runs.append((start, prev)); start = prev = c
            runs.append((start, prev))
            for (a, b) in runs:
                seg_id[(h, a, b)] = nid
                nid += 1
        def find_seg(h, c):
            for (hh, a, b), i in seg_id.items():
                if hh == h and a <= c <= b:
                    return i
            raise KeyError((h, c))
        edges = {}
        for (x, y), (c, r) in pos.items():
            top = find_seg(r, c)
            bot = find_seg(r - 1, c)
            edges[(x, y)] = (top, bot)
        # poles: the top-most segment (height n-2) and the bottom-most (height 0)
        top_pole = find_seg(n - 2, pos[(n, n - 2)][0]) if (n, n - 2) in pos else find_seg(n - 2, 1)
        bot_pole = find_seg(0, 1)
        return edges, nid, top_pole, bot_pole
    else:
        # vertical: mirror the construction on columns
        seg_id = {}
        nid = 0
        maxcol = n - 2
        for vx in range(0, n - 1):
            rows = sorted({r for (x, y), (c, r) in pos.items() if c == vx or c == vx + 1})
            if not rows:
                continue
            runs = []
            start = prev = rows[0]
            for c in rows[1:]:
                if c == prev + 1:
                    prev = c
                else:
                    runs.append((start, prev)); start = prev = c
            runs.append((start, prev))
            for (a, b) in runs:
                seg_id[(vx, a, b)] = nid
                nid += 1
        def find_segv(vx, r):
            for (xx, a, b), i in seg_id.items():
                if xx == vx and a <= r <= b:
                    return i
            raise KeyError((vx, r))
        edges = {}
        for (x, y), (c, r) in pos.items():
            left = find_segv(c - 1, r)
            right = find_segv(c, r)
            edges[(x, y)] = (left, right)
        left_pole = find_segv(0, 1)
        right_pole = find_segv(n - 2, pos[(n, 1)][1]) if (n, 1) in pos else find_segv(n - 2, 1)
        return edges, nid, left_pole, right_pole

def kirchhoff(edges, nv, a=None, b=None):
    """spanning tree count (a=b=None) or effective resistance a-b (unit conductances), exact."""
    L = [[F(0)] * nv for _ in range(nv)]
    for (u, v) in edges.values():
        L[u][u] += 1; L[v][v] += 1
        L[u][v] -= 1; L[v][u] -= 1
    def det_minor(drop):
        idx = [i for i in range(nv) if i not in drop]
        M = [[L[i][j] for j in idx] for i in idx]
        d = F(1)
        k = len(M)
        for c in range(k):
            piv = None
            for r in range(c, k):
                if M[r][c] != 0:
                    piv = r; break
            if piv is None:
                return F(0)
            if piv != c:
                M[c], M[piv] = M[piv], M[c]
                d = -d
            d *= M[c][c]
            inv = 1 / M[c][c]
            for r in range(c + 1, k):
                f = M[r][c] * inv
                if f:
                    for cc in range(c, k):
                        M[r][cc] -= f * M[c][cc]
        return d
    kappa = det_minor({0})
    if a is None:
        return kappa, None
    # R_ab = kappa(G/ab-identified... standard: R = det L^{(a,b)} / det L^{(a)}
    num = det_minor({a, b})
    return kappa, num / kappa

print("=" * 96)
print("PART 1+2: the Smith network of the staircase, n=4..10; self-duality")
print("=" * 96)
for n in range(4, 11):
    eh, nvh, tp, bp = smith_network(n, True)
    ev, nvv, lp, rp = smith_network(n, False)
    kh, Rh = kirchhoff(eh, nvh, tp, bp)
    kv, Rv = kirchhoff(ev, nvv, lp, rp)
    degh = sorted(Counter([u for (u, v) in eh.values()] + [v for (u, v) in eh.values()]).values())
    degv = sorted(Counter([u for (u, v) in ev.values()] + [v for (u, v) in ev.values()]).values())
    m = comb(n - 1, 2)
    print("n=%2d m=%2d | H-net: V=%2d rank=%2d kappa=%-10d R=%-8s | V-net: V=%2d kappa=%-10d R=%-8s | dual-iso? deg:%s kappa:%s R:%s"
          % (n, m, nvh, m - nvh + 1, kh, Rh, nvv, kv, Rv,
             degh == degv, kh == kv, Rh == Rv))

print()
print("=" * 96)
print("PART 3: the battery table (leg-law KCL bookkeeping), n=5,6")
print("=" * 96)
for n in [5, 6]:
    cells = staircase_cells(n)
    m = len(cells)
    tidx = {t: i for i, t in enumerate(cells)}
    gmap = [tidx[(n - y + 1, n - x + 1)] for (x, y) in cells]
    pairs = [(i, j) for i in range(n) for j in range(n) if i < j]
    pidx = {p: k for k, p in enumerate(pairs)}
    perms = list(permutations(range(n)))

    def beats_of(tv):
        B = [[False] * n for _ in range(n)]
        for k in range(2, n + 1):
            B[k - 1][k - 2] = True
        for (x, y), i in tidx.items():
            if tv[i] == 1:
                B[x - 1][y - 1] = True
            else:
                B[y - 1][x - 1] = True
        return B

    def canon_key(B):
        best = None
        for pm in perms:
            x = 0
            for k, (i, j) in enumerate(pairs):
                pi, pj = pm[i], pm[j]
                lo, hi = (pi, pj) if pi < pj else (pj, pi)
                bit = B[i][j] if pm[i] == lo else (not B[i][j])
                x |= (1 if bit else 0) << pidx[(lo, hi)]
            if best is None or x < best:
                best = x
        return best

    blue_id_ok = True
    flux = defaultdict(int)
    gsn = defaultdict(int)
    tot = defaultdict(int)
    least_word = {}
    for t in range(1 << m):
        tv = [(t >> i) & 1 for i in range(m)]
        gs = all(tv[i] == tv[gmap[i]] for i in range(m))
        B = beats_of(tv)
        key = canon_key(B)
        okey = canon_key([[B[v][u] for u in range(n)] for v in range(n)])
        mk = min(key, okey)
        # legs: e1 = out-tile-degree of vertex 1 = #tiles (x,1) with bit 0 (1 beats x? convention:
        # bit 1 = x beats y; vertex 1 wins tile (x,1) iff bit=0). e_n = #tiles (n,y) with bit 1.
        e1 = sum(1 for (x, y), i in tidx.items() if y == 1 and tv[i] == 0)
        en = sum(1 for (x, y), i in tidx.items() if x == n and tv[i] == 1)
        if gs and e1 + en != n - 2:
            blue_id_ok = False
            print("  BLUE IDENTITY FAILS at tiling %d: e1=%d en=%d" % (t, e1, en))
        flux[mk] += e1 - en
        gsn[mk] += gs
        tot[mk] += 1
        if mk not in least_word or t < least_word[mk]:
            least_word[mk] = t
    total_flux = sum(flux.values())
    bytype = defaultdict(list)
    for mk in flux:
        ty = "PB" if gsn[mk] == tot[mk] else ("PBk" if gsn[mk] == 0 else "MX")
        bytype[ty].append(flux[mk])
    print("n=%d: blue identity e1+e_n=n-2 on all grid-sym tilings: %s ; GLOBAL flux sum = %d (expect 0)"
          % (n, blue_id_ok, total_flux))
    for ty in ["PB", "MX", "PBk"]:
        fl = sorted(bytype[ty])
        nz = sum(1 for f in fl if f != 0)
        print("   %-3s classes=%2d  zero-flux=%2d  nonzero fluxes=%s"
              % (ty, len(fl), len(fl) - nz, [f for f in fl if f != 0]))
    # PART 4 note: Bouwkamp fiber code = least_word; report count of distinct (trivially all) +
    # whether least-word order refines type order
    lw = sorted(least_word.items(), key=lambda kv: kv[1])
    seq = []
    for mk, w in lw:
        ty = "PB" if gsn[mk] == tot[mk] else ("PBk" if gsn[mk] == 0 else "MX")
        seq.append(ty)
    print("   Bouwkamp-order type word: %s" % "".join({"PB": "B", "MX": "M", "PBk": "K"}[t] for t in seq))
print("\nDONE")
