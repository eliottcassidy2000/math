#!/usr/bin/env python3
"""
The star (point-adjacency) quotient of the tiling hypercube   (mac-mini-2026-07-20-S128)
========================================================================================
Owner's question, three parts:
  (1) single-tile flip = EDGE adjacency: generates everything, so no invariants.
      The map-graph move is POINT adjacency: a point of K_n is a vertex, its clique
      is a STAR, stars are INCIDENCE ROWS -- whence duality.  Does the clique-at-a-
      point feature produce the invariants?
  (2) Is the hypercube quotient by the star group the merged metagraph, a refinement,
      or unrelated?
  (3) Do the cycle invariants descend to iso classes (tournament functions), or are
      they only tiling functions?

Setup.  Vertices 0..n-1, base path n-1 -> n-2 -> ... -> 0.  TILES = pairs (i,j), i<j,
with j-i >= 2 (i.e. K_n MINUS the base path); m = C(n-1,2).  A tiling is a mask in
F_2^m.  STAR(v) = the tiles incident to v = the incidence row of the TILE GRAPH at v.

Claim under test: <stars> = the CUT space of the tile graph (dim n-1 when connected),
so the quotient is dual to the CYCLE space and #orbits = 2^(m-n+1).  The orbit
invariant is then the F_2 holonomy ("Wilson loop") around each cycle of the tile graph.
"""
import numpy as np
from itertools import permutations

def build(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    idx = {p: k for k, p in enumerate(pairs)}
    tiles = [(i, j) for (i, j) in pairs if j - i >= 2]      # K_n minus the base path
    return pairs, idx, tiles, len(pairs), len(tiles)

def tile_graph_dims(n):
    """cut dim = n - components; cycle dim = m - cut dim."""
    _, _, tiles, _, m = build(n)
    adj = {v: set() for v in range(n)}
    for (i, j) in tiles:
        adj[i].add(j); adj[j].add(i)
    seen, comps = set(), 0
    for v in range(n):
        if v in seen: continue
        comps += 1; stack = [v]; seen.add(v)
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if w not in seen: seen.add(w); stack.append(w)
    cut = n - comps
    return m, comps, cut, m - cut

def rref_f2(vecs, m):
    basis, piv = [], []
    for v in vecs:
        for b, p in zip(basis, piv):
            if v >> p & 1: v ^= b
        if v:
            p = v.bit_length() - 1
            basis.append(v); piv.append(p)
    return basis, piv

def reduce_f2(x, basis, piv):
    for b, p in zip(basis, piv):
        if x >> p & 1: x ^= b
    return x

def analyse(n):
    pairs, idx, tiles, E, m = build(n)
    tpos = [idx[t] for t in tiles]                       # tile -> edge position
    N = 1 << m
    # orientation array: bit=1 means i->j (low beats high); base edges fixed at 0
    A = np.zeros((N, E), dtype=np.uint8)
    masks = np.arange(N, dtype=np.int64)
    for k, e in enumerate(tpos):
        A[:, e] = ((masks >> k) & 1).astype(np.uint8)
    pow2 = (1 << np.arange(E, dtype=np.int64))

    # --- iso classes: canonical code = min over all n! relabelings -----------------
    canon = None
    for p in permutations(range(n)):
        src = np.empty(E, dtype=np.int64); fl = np.zeros(E, dtype=np.uint8)
        for e, (i, j) in enumerate(pairs):
            a, b = p[i], p[j]
            t = idx[(min(a, b), max(a, b))]
            src[t] = e; fl[t] = 1 if a > b else 0
        code = (A[:, src] ^ fl) @ pow2
        canon = code if canon is None else np.minimum(canon, code)
    cls = {}
    for x, c in enumerate(canon.tolist()):
        cls.setdefault(c, []).append(x)

    # --- star group (point adjacency) ---------------------------------------------
    stars = []
    for v in range(n):
        b = 0
        for k, (i, j) in enumerate(tiles):
            if i == v or j == v: b |= 1 << k
        stars.append(b)
    basis, piv = rref_f2(stars, m)
    orb = {}
    for x in range(N):
        orb.setdefault(reduce_f2(x, basis, piv), []).append(x)

    # --- the tournament 3-cycle count: C(n,3) - sum_v C(s_v,2) ---------------------
    scores = np.zeros((N, n), dtype=np.int64)
    for e, (i, j) in enumerate(pairs):
        bit = A[:, e].astype(np.int64)
        scores[:, i] += bit; scores[:, j] += 1 - bit
    c3 = n * (n - 1) * (n - 2) // 6 - (scores * (scores - 1) // 2).sum(axis=1)

    def const_on(parts, f):
        return all(len(set(f[x] for x in P)) == 1 for P in parts.values())
    def refines(fine, coarse_of):                      # every block of `fine` inside one block
        return all(len(set(coarse_of[x] for x in P)) == 1 for P in fine.values())

    cls_of = {x: c for c, P in cls.items() for x in P}
    orb_of = {x: o for o, P in orb.items() for x in P}
    m2, comps, cut, cyc = tile_graph_dims(n)
    return dict(n=n, m=m, comps=comps, cut=cut, cyc=cyc,
                dim=len(basis), orbits=len(orb), classes=len(cls),
                pred=1 << cyc,
                star_ref_iso=refines(orb, cls_of), iso_ref_star=refines(cls, orb_of),
                c3_on_iso=const_on(cls, c3), c3_on_star=const_on(orb, c3))

print("=== the tile graph = K_n minus the base path: cut/cycle dimensions ===")
print(f"{'n':>3} {'m':>4} {'comps':>6} {'cut=n-1':>8} {'cycle=m-n+1':>12} {'2^cycle':>10}")
for n in range(4, 10):
    m, comps, cut, cyc = tile_graph_dims(n)
    print(f"{n:>3} {m:>4} {comps:>6} {cut:>8} {cyc:>12} {1<<cyc:>10}")

print()
print("=== star group vs isomorphism, by exhaustive enumeration ===")
for n in range(4, 8):
    r = analyse(n)
    print(f"n={r['n']}: m={r['m']}  star-group dim={r['dim']} (cut dim n-1={r['cut']})"
          f"  orbits={r['orbits']} (predicted 2^(m-n+1)={r['pred']})  iso classes={r['classes']}")
    print(f"      star-orbits refine iso-classes? {r['star_ref_iso']}"
          f"   iso-classes refine star-orbits? {r['iso_ref_star']}")
    print(f"      3-cycle count constant on ISO classes? {r['c3_on_iso']}"
          f"   on STAR orbits? {r['c3_on_star']}")
    assert r['dim'] == r['cut'] and r['orbits'] == r['pred'], "cut/cycle identification FAILED"

print()
print("VERDICT")
print("  (1) <stars> = the CUT space of the tile graph (dim n-1, tile graph connected).")
print("      Quotient = dual of the CYCLE space, #orbits = 2^(m-n+1) -- exact match.")
print("      So point adjacency DOES produce invariants: the F_2 holonomies (cycle parities).")
print("  (2) NOT the merged metagraph and NOT a refinement: for n>=5 neither partition")
print("      refines the other. The star quotient is TRANSVERSE to isomorphism.")
print("  (3) Tournament cycle invariants (3-cycle count) descend to ISO classes but NOT")
print("      to star orbits; the star-orbit class is a TILING function that does not")
print("      descend to iso classes. The two invariant systems are transverse.")
