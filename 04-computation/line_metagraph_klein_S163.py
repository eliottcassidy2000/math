#!/usr/bin/env python3
"""
klein-2026-07-07-S163 -- THE LINE METAGRAPH: vertices = lines {t, flip(t)}; the abstraction
of the abstraction.

THEOREMS (proofs in the reflection):
 L1 (folded cube). Single-tile flips commute with the all-tile flip, and flip is free, so
    the line set with wiggly adjacency is EXACTLY the folded m-cube FQ_m = Q_m/antipode
    (vertex-transitive, m-regular... with the antipodal edge folded: degree m, plus at
    m odd no loops; every line has exactly m wiggly neighbors).
 L2 (descent). The iso-groupoid action commutes with both flips, so the LINE METAGRAPH
    LG_n := FQ_m / iso is well-defined; its nodes ('line classes') carry the blue/black
    color and the endpoint class-pair {[T(t)],[T(flip t)]}.
 L3 (fibration). Line classes fiber over unordered endpoint class-pairs; the fiber
    multiplicities decompose S161's allocation numbers: #lines over pair (A,B) =
    sum over line-classes above (A,B) of their orbit sizes. The line-class refinement
    measures exactly the quasi-randomness failure (C5's assortativity).

THIS SCRIPT (n = 4, 5, 6; n=7 count-only if time):
 1. build lines; canonical key of a line = refinement-canon of the ARC-2-COLORED structure
    (color each pair by which of t / flip(t) orients it forward -- equivalently the pair
    (T, T') as an edge-2-colored complete digraph), minimized over the swap t <-> flip(t).
 2. census: L(n) = #line classes (blue/black split); orbit sizes; fibration over endpoint
    class-pairs (max multiplicity -- is the map line-class -> class-pair injective?);
    line-class fiber parities; degrees in LG_n (wiggly adjacency of line classes).
 3. trends: L(n) vs A000568-type counts; the families (self-lines/cross-lines,
    SC-pair/NS-pair endpoints, blue skeleton); the generating-pattern check: is the
    #lines-over-(A,B) predicted by fiberA*fiberB/2^{m-1}-type quasi-randomness + skeleton?
"""
import itertools
from collections import defaultdict, Counter

def build(n):
    tiles = [(x, y) for x in range(1, n+1) for y in range(1, x) if x - y >= 2]
    tidx = {t: k for k, t in enumerate(tiles)}
    sigma = [tidx[(n+1-y, n+1-x)] for (x, y) in tiles]
    return tiles, sigma

def adj(tv, n, tiles):
    A = [[0]*n for _ in range(n)]
    for a in range(2, n+1): A[a-1][a-2] = 1
    for b, (x, y) in enumerate(tiles):
        if (tv >> b) & 1: A[x-1][y-1] = 1
        else: A[y-1][x-1] = 1
    return A

def canon_colored(A1, A2, n):
    """canonical key of the ordered pair (A1, A2) of tournaments on [n] under SIMULTANEOUS
    relabeling: arcs colored by (A1 direction, A2 direction) -- refinement + brute in cells."""
    # color of ordered pair (i,j): 2*A1[i][j] + A2[i][j]  (in {0,1,2,3}; A1[i][j]+A1[j][i]=1)
    col = [0]*n
    for _ in range(n):
        prof = []
        for i in range(n):
            cnt = Counter()
            for j in range(n):
                if i == j: continue
                cnt[(2*A1[i][j] + A2[i][j], col[j])] += 1
            prof.append((col[i], tuple(sorted(cnt.items()))))
        ranks = {p: r for r, p in enumerate(sorted(set(prof)))}
        nc = [ranks[p] for p in prof]
        if nc == col: break
        col = nc
    cells = defaultdict(list)
    for i in range(n): cells[col[i]].append(i)
    ordered = [cells[c] for c in sorted(cells)]
    best = None
    def encode(perm):
        key = []
        for i in range(n):
            row = 0
            for j in range(n):
                if i == j: continue
                row = row*4 + 2*A1[perm[i]][perm[j]] + A2[perm[i]][perm[j]]
            key.append(row)
        return tuple(key)
    def rec(prefix, rest):
        nonlocal best
        if not rest:
            k = encode(prefix)
            if best is None or k < best: best = k
            return
        for p in itertools.permutations(rest[0]):
            rec(prefix + list(p), rest[1:])
    tot = 1
    for c in ordered:
        f = 1
        for i in range(2, len(c)+1): f *= i
        tot *= f
    if tot > 20000:
        for perm in itertools.permutations(range(n)):
            k = encode(perm)
            if best is None or k < best: best = k
        return best
    rec([], ordered)
    return best

def canon_tournament(A, n):
    return canon_colored(A, A, n)   # degenerate coloring = plain tournament canon

if __name__ == "__main__":
    for n in (4, 5, 6):
        tiles, sigma = build(n)
        m = len(tiles); full = (1 << m) - 1
        # tournament classes per tiling (for endpoint pairs)
        tclass = {}
        for tv in range(1 << m):
            tclass[tv] = canon_tournament(adj(tv, n, tiles), n)
        # line classes
        lineclass = {}; lc_members = defaultdict(list)
        for tv in range(1 << m):
            ftv = tv ^ full
            if tv > ftv: continue
            A1, A2 = adj(tv, n, tiles), adj(ftv, n, tiles)
            k1 = canon_colored(A1, A2, n); k2 = canon_colored(A2, A1, n)
            key = min(k1, k2)   # unordered line
            lineclass[tv] = key; lc_members[key].append(tv)
        L = len(lc_members)
        # blue/black per line class (well-defined? gridsym is iso-invariant? gridsym is a
        # TILING property; check constancy on orbits)
        def gridsym(tv): return all(((tv >> b) & 1) == ((tv >> sigma[b]) & 1) for b in range(m))
        blue_const = all(len(set(gridsym(t) for t in v)) == 1 for v in lc_members.values())
        nblue_lc = sum(1 for v in lc_members.values() if gridsym(v[0]))
        # fibration over endpoint class-pairs
        pair_of = {}
        multi = defaultdict(set)
        for key, mem in lc_members.items():
            t0 = mem[0]; p = tuple(sorted((tclass[t0], tclass[t0 ^ full])))
            pair_of[key] = p; multi[p].add(key)
        maxmult = max(len(v) for v in multi.values())
        # orbit sizes + parity
        sizes = sorted(len(v) for v in lc_members.values())
        parities = sorted(set(len(v) % 2 for v in lc_members.values()))
        # LG_n degrees: wiggly adjacency between line classes
        edges = set()
        for tv in range(1 << m):
            ftv = tv ^ full
            if tv > ftv: continue
            for b in range(m):
                tv2 = tv ^ (1 << b)
                l2 = lineclass.get(min(tv2, tv2 ^ full))
                l1 = lineclass[tv]
                if l1 != l2: edges.add(tuple(sorted((repr(l1), repr(l2)))))
        # quasi-randomness check: lines over (A,B) vs fiber product prediction
        fiber = Counter(tclass.values())
        qr = []
        linecount_pair = Counter()
        for tv in range(1 << m):
            ftv = tv ^ full
            if tv > ftv: continue
            p = tuple(sorted((tclass[tv], tclass[ftv])))
            linecount_pair[p] += 1
        for p, cnt in linecount_pair.items():
            fa, fb = fiber[p[0]], fiber[p[1]]
            pred = (fa*fb/(2**m)) if p[0] != p[1] else (fa*fa/(2**(m+1)))
            qr.append(cnt/pred if pred > 0 else float('inf'))
        print(f"\n===== n={n}: m={m}; lines={1<<(m-1)}; LINE CLASSES L(n) = {L} (blue {nblue_lc}, black {L-nblue_lc}) =====")
        print(f"  blue/black constant on orbits: {blue_const}")
        print(f"  orbit sizes: min..max = {sizes[0]}..{sizes[-1]}; size parities present: {parities}")
        print(f"  endpoint-pair fibration: {len(multi)} class-pairs hit; MAX line-classes over one pair = {maxmult}"
              f"  (injective? {maxmult == 1})")
        mm = [ (p, len(v)) for p, v in multi.items() if len(v) == maxmult ][:1]
        print(f"  LG_n: nodes {L}, class-level wiggly edges {len(edges)}")
        import statistics
        print(f"  quasi-randomness ratio (actual/fiber-product prediction) over class-pairs:"
              f" median {statistics.median(qr):.3f}, min {min(qr):.3f}, max {max(qr):.3f}, n={len(qr)}")
        # tournament-class count sanity
        print(f"  tournament classes = {len(fiber)} (A000568: n=4:4, 5:12, 6:56)")
