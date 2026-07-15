#!/usr/bin/env python3
"""THM-787 study: the flow of transitivity on the (merged) metagraph (opus-S304).

Owner directive: organize iso-class nodes on the spectrum from the TRANSITIVE
node to the maximally-distributed score class; trace the flow of transitivity
along BLUE lines (grid-symmetric tiling <-> its tile-complement) out of pure
blue classes into mixed classes, then along BLACK lines into pure black
classes; quantify the left-right symmetry (blue) and imbalance (black).

Conventions (CLAUDE.md, strict): base path n -> n-1 -> ... -> 1 (arcs k -> k-1
fixed); tiles = pairs (x,y), x - y >= 2, bit 1 = arc x -> y, bit 0 = arc y -> x;
a LINE joins tiling t to flip-ALL-tiles(t) (the d = m waggly layer); the line is
BLUE iff isGridSym(t): t(x,y) = t(n-y+1, n-x+1) for every tile; PURE BLUE class:
all fiber tilings grid-sym; PURE BLACK: none; MIXED: both.

AXIS: x(class) = sum_i (2 s_i - (n-1))^2 (integer; transitive = max; regular /
near-regular = min; invariant under reversal => well-defined on the merged
metagraph G_n/Z_2).

Iso classes via a strong invariant, VALIDATED by the known counts A000568 =
1, 2, 4, 12, 56, 456 (n = 2..7) and by the fiber x Aut = H identity per class.
"""
import sys, itertools
from collections import defaultdict

def study(n, out):
    V = list(range(1, n + 1))
    tiles = [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2]
    m = len(tiles)
    tile_idx = {t: i for i, t in enumerate(tiles)}
    # grid reflection on tiles: (x,y) -> (n-y+1, n-x+1)
    refl = [tile_idx[(n - y + 1, n - x + 1)] for (x, y) in tiles]

    def tourn(bits):
        """adjacency: adj[u] = bitmask of out-neighbours (vertices 1..n -> bits 0..n-1)."""
        adj = [0] * (n + 1)
        for k in range(n, 1, -1): adj[k] |= 1 << (k - 2)          # path arc k -> k-1
        for i, (x, y) in enumerate(tiles):
            if bits >> i & 1: adj[x] |= 1 << (y - 1)
            else: adj[y] |= 1 << (x - 1)
        return adj

    def scores(adj): return sorted(bin(adj[v]).count('1') for v in V)

    def invariant(adj):
        s = {v: bin(adj[v]).count('1') for v in V}
        # 3-cycle count through each vertex + out/in score multisets, two refinement rounds
        prof = {}
        for v in V:
            outs = sorted(s[u] for u in V if adj[v] >> (u - 1) & 1)
            ins = sorted(s[u] for u in V if adj[u] >> (v - 1) & 1)
            c3 = 0
            for u in V:
                if adj[v] >> (u - 1) & 1:
                    for w in V:
                        if adj[u] >> (w - 1) & 1 and adj[w] >> (v - 1) & 1: c3 += 1
            prof[v] = (s[v], tuple(outs), tuple(ins), c3)
        # round 2: multiset of neighbours' round-1 profiles
        prof2 = {}
        for v in V:
            po = sorted(prof[u] for u in V if adj[v] >> (u - 1) & 1)
            pi = sorted(prof[u] for u in V if adj[u] >> (v - 1) & 1)
            prof2[v] = (prof[v], tuple(po), tuple(pi))
        # arc-type multiset: for each arc u->v, common out/in neighbourhood sizes
        arcs = []
        for u in V:
            for v in V:
                if adj[u] >> (v - 1) & 1:
                    cuo = bin(adj[u] & adj[v]).count('1')
                    a_in_u = 0
                    for w in V:
                        if adj[w] >> (u - 1) & 1 and adj[w] >> (v - 1) & 1: a_in_u += 1
                    arcs.append((s[u], s[v], cuo, a_in_u))
        return (tuple(sorted(prof2.values())), tuple(sorted(arcs)))

    def ham_paths(adj):
        full = (1 << n) - 1
        dp = defaultdict(int)
        for v in V: dp[(1 << (v - 1), v)] = 1
        for size in range(1, n):
            items = [(k, c) for k, c in dp.items() if bin(k[0]).count('1') == size]
            for (mask, v), c in items:
                for u in V:
                    b = 1 << (u - 1)
                    if not mask & b and adj[v] & b:
                        dp[(mask | b, u)] += c
        return sum(c for (mask, v), c in dp.items() if mask == full)

    def reversed_adj(adj):
        radj = [0] * (n + 1)
        for u in V:
            for v in V:
                if adj[u] >> (v - 1) & 1: radj[v] |= 1 << (u - 1)
        return radj

    # sweep all tilings (invariant strengthened with H: buckets are unions of true
    # classes, so bucket-count == A000568(n) certifies EXACT classification)
    classes = {}                      # invariant -> class id
    cls_of = {}                       # tiling -> class id
    fiber = defaultdict(int)
    gs_count = defaultdict(int)
    rep = {}
    for bits in range(1 << m):
        adj = tourn(bits)
        inv = invariant(adj) + (ham_paths(adj),)
        if inv not in classes:
            classes[inv] = len(classes); rep[classes[inv]] = bits
        c = classes[inv]
        cls_of[bits] = c
        fiber[c] += 1
        if all((bits >> i & 1) == (bits >> refl[i] & 1) for i in range(m)):
            gs_count[c] += 1
    C = len(classes)
    counts = {2:1, 3:2, 4:4, 5:12, 6:56, 7:456}
    ok_count = (C == counts.get(n, C))
    assert ok_count, f"n={n}: {C} buckets != A000568 -- invariant incomplete"
    # class data
    data = {}
    for c in range(C):
        adj = tourn(rep[c])
        sc = scores(adj); H = ham_paths(adj)
        assert H % fiber[c] == 0, (n, c, H, fiber[c])   # fiber * Aut = H
        aut = H // fiber[c]
        x = sum((2 * si - (n - 1)) ** 2 for si in sc)
        typ = ('pureblue' if gs_count[c] == fiber[c] else
               'pureblack' if gs_count[c] == 0 else 'mixed')
        radj = reversed_adj(adj)
        data[c] = dict(sc=tuple(sc), H=H, fib=fiber[c], aut=aut, x=x, typ=typ,
                       rcls=None, rinv=invariant(radj) + (ham_paths(radj),))
    inv_to_c = {inv: c for inv, c in classes.items()}
    for c in range(C): data[c]['rcls'] = inv_to_c[data[c]['rinv']]
    SC = sum(1 for c in range(C) if data[c]['rcls'] == c)
    # lines: t <-> all-flip(t)
    fullmask = (1 << m) - 1
    lines = []                        # (colour, cA, cB) unordered per tiling-pair
    seen = set()
    for bits in range(1 << m):
        if bits in seen: continue
        pb = bits ^ fullmask
        seen.add(bits); seen.add(pb)
        blue = all((bits >> i & 1) == (bits >> refl[i] & 1) for i in range(m))
        lines.append((blue, cls_of[bits], cls_of[pb]))
    # assemble the report
    P = lambda *a: print(*a, file=out)
    P(f"===== n = {n}: {C} classes (A000568 check {'OK' if ok_count else 'FAIL'}), "
      f"SC = {SC}, merged nodes = {(C + SC) // 2}, tiles m = {m}, tilings = {1 << m}")
    xs = sorted(set(d['x'] for d in data.values()), reverse=True)
    P(f"  axis levels (x = sum(2s-(n-1))^2, transitive -> distributed): {xs}")
    for typ in ('pureblue', 'mixed', 'pureblack'):
        cs = [c for c in range(C) if data[c]['typ'] == typ]
        hist = defaultdict(int)
        for c in cs: hist[data[c]['x']] += 1
        P(f"  {typ:<9}: {len(cs):4d} classes; axis histogram {dict(sorted(hist.items(), reverse=True))}")
    # line statistics
    stats = defaultdict(int); dx_blue = defaultdict(int); dx_black = defaultdict(int)
    pair_blue = defaultdict(int); pair_black = defaultdict(int)
    for blue, a, b in lines:
        ta, tb = data[a]['typ'], data[b]['typ']
        key = (('BLUE' if blue else 'BLACK'),) + tuple(sorted([ta, tb]))
        stats[key] += 1
        dxa = data[b]['x'] - data[a]['x']
        pair = tuple(sorted([data[a]['x'], data[b]['x']]))
        if blue: dx_blue[abs(dxa)] += 1; pair_blue[pair] += 1
        else: dx_black[abs(dxa)] += 1; pair_black[pair] += 1
    P("  lines by (colour, endpoint types):")
    for k in sorted(stats): P(f"    {k}: {stats[k]}")
    P(f"  BLUE |dx| histogram: {dict(sorted(dx_blue.items()))}")
    P(f"  BLACK |dx| histogram: {dict(sorted(dx_black.items()))}")
    # level vs skew lines; net black flux per axis position
    lvl_b = sum(v for k, v in dx_blue.items() if k == 0); tot_b = sum(dx_blue.values())
    lvl_k = sum(v for k, v in dx_black.items() if k == 0); tot_k = sum(dx_black.values())
    P(f"  BLUE lines level (dx=0): {lvl_b}/{tot_b};  BLACK level: {lvl_k}/{tot_k}")
    # score-multiset comparison across lines: does B's score multiset = A's or its reversal?
    same_sc_b = same_sc_k = 0
    for blue, a, b in lines:
        eq = data[a]['sc'] == data[b]['sc']
        if blue: same_sc_b += eq
        else: same_sc_k += eq
    P(f"  lines with EQUAL score multisets: blue {same_sc_b}/{tot_b}, black {same_sc_k}/{tot_k}")
    # black net flux: orient each black line from higher x to lower x; flux across levels
    flux = defaultdict(int)
    for blue, a, b in lines:
        if blue: continue
        xa, xb = data[a]['x'], data[b]['x']
        if xa != xb:
            hi, lo = max(xa, xb), min(xa, xb)
            flux[(hi, lo)] += 1
    P(f"  BLACK downhill line counts by (x_hi -> x_lo): {dict(sorted(flux.items(), reverse=True))}")
    # the transitive node's lines
    tc = [c for c in range(C) if data[c]['H'] == 1][0]
    tl = [(('BLUE' if blue else 'BLACK'), data[b]['x'], data[b]['typ'])
          for blue, a, b in lines if a == tc] + \
         [(('BLUE' if blue else 'BLACK'), data[a]['x'], data[a]['typ'])
          for blue, a, b in lines if b == tc and a != tc]
    P(f"  transitive node (x = {data[tc]['x']}, type {data[tc]['typ']}): partner profile {sorted(set(tl))}")
    return dict(n=n, C=C, SC=SC, data=data, lines=lines)

if __name__ == '__main__':
    out = sys.stdout
    for n in range(3, 8):
        study(n, out)
