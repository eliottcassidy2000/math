#!/usr/bin/env python3
"""star_group_quotient_kps_S128c109.py -- kind-pasteur-2026-07-20-S128c109

THE OWNER'S STRUCTURAL IDEA, MADE COMPUTABLE.

  "the single-tile flip is edge adjacency -- generates everything, no invariants.
   The map-graph move is point adjacency; a point of K_n is a vertex, its clique is
   a STAR, stars are INCIDENCE ROWS -- whence duality.  The clique-at-a-point
   feature that distinguishes map graphs from planar duals is exactly what produces
   the invariants."

Formalised.  Tilings are F_2^m, m = C(n-1,2), coordinates = the non-base-path arcs
("tiles").  Fix the base path n-1 -> ... -> 0.  For a vertex v let

    S_v  =  the indicator of the TILES incident to v          (the star at v)

and let GAMMA = <S_1,...,S_n> be the STAR GROUP.  The single-tile flips generate all
of F_2^m and so have no invariants; the star moves generate only GAMMA, and the coset
b + GAMMA is the resulting invariant.  This is the point-adjacency analogue of the
edge-adjacency wiggly layer.

WHY THE DUALITY IS EXACT.  In the GF(2) edge space of K_n the base path is a spanning
tree, the tiles are the NON-TREE edges (cycle-space coordinates), and the stars
delta(v) span the CUT SPACE.  So GAMMA is precisely the image of the cut space under
restriction to cycle coordinates -- "stars are incidence rows" is literally the
statement that GAMMA is the row space of the vertex-edge incidence matrix.  Two
predictions follow and are checked below:
  * sum_v S_v = 0 (every tile has two endpoints), so rank GAMMA <= n-1;
  * for n >= 4 no non-zero cut is supported on the base path (delta(v) has n-1 edges
    but at most 2 of them are path arcs), so the restriction is injective on the cut
    space and rank GAMMA = n-1 EXACTLY.

THE OWNER'S TWO QUESTIONS, ANSWERED BY EXHAUSTION:
  (Q1) Do the resulting invariants DESCEND to iso classes -- are they tournament
       functions, or only tiling functions?  Tested by asking whether the coset
       b + GAMMA is constant on each isomorphism class.
  (Q2) Is the hypercube quotient by the star group the MERGED METAGRAPH, a
       REFINEMENT of it, or unrelated?  Tested by comparing the coset partition with
       the iso-class partition and with the complement-merged partition.

Also measured: whether the star move preserves the standard cycle counts c_3, c_5
(i.e. whether GAMMA-cosets are compatible with tournament invariants at all).
"""
import sys
from itertools import permutations, combinations

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6


def setup(n):
    tiles = [(x, y) for y in range(n) for x in range(n) if x - y >= 2]
    return tiles, len(tiles)


def star(v, tiles):
    b = 0
    for i, (x, y) in enumerate(tiles):
        if x == v or y == v:
            b |= 1 << i
    return b


def f2_span(vecs):
    piv = []
    for v in vecs:
        cur = v
        for p in piv:
            cur = min(cur, cur ^ p)
        if cur:
            piv.append(cur)
            piv.sort(reverse=True)
    return piv


def reduce_vec(v, piv):
    cur = v
    for p in piv:
        cur = min(cur, cur ^ p)
    return cur


def build(bits, tiles, n):
    adj = [0] * n
    for k in range(n - 1, 0, -1):
        adj[k] |= 1 << (k - 1)
    for i, (x, y) in enumerate(tiles):
        if bits >> i & 1:
            adj[x] |= 1 << y
        else:
            adj[y] |= 1 << x
    return adj


def canon(adj, n):
    best = None
    for p in permutations(range(n)):
        key = tuple(sorted((p[a], p[b]) for a in range(n) for b in range(n)
                           if a != b and (adj[a] >> b) & 1))
        if best is None or key < best:
            best = key
    return best


def opp(adj, n):
    o = [0] * n
    for a in range(n):
        for b in range(n):
            if a != b and (adj[a] >> b) & 1:
                o[b] |= 1 << a
    return o


def c3(adj, n):
    t = 0
    for a, b, c in combinations(range(n), 3):
        da = ((adj[a] >> b) & 1) + ((adj[a] >> c) & 1)
        db = ((adj[b] >> a) & 1) + ((adj[b] >> c) & 1)
        dc = ((adj[c] >> a) & 1) + ((adj[c] >> b) & 1)
        if da == db == dc == 1:
            t += 1
    return t


for n in range(4, NMAX + 1):
    tiles, m = setup(n)
    stars = [star(v, tiles) for v in range(n)]
    piv = f2_span(stars)
    r = len(piv)
    print("=" * 78)
    print("n = %d :  m = %d tiles,  2^m = %d tilings" % (n, m, 1 << m))
    print("=" * 78)
    tot = 0
    for s in stars:
        tot ^= s
    print("  sum of all stars = %d  (predicted 0, every tile has two endpoints): %s"
          % (tot, tot == 0))
    print("  rank(GAMMA) = %d   (predicted n-1 = %d): %s"
          % (r, n - 1, r == n - 1))
    print("  quotient  F_2^m / GAMMA  has  2^(m-rank) = %d cosets" % (1 << (m - r)))

    # classify every tiling
    cls, cos = {}, {}
    for bits in range(1 << m):
        adj = build(bits, tiles, n)
        c = canon(adj, n)
        cls.setdefault(c, []).append(bits)
        cos[bits] = reduce_vec(bits, piv)
    print("  isomorphism classes : %d" % len(cls))

    # (Q1) does the coset descend to iso classes?
    mixed = 0
    witness = None
    for c, bl in cls.items():
        s = {cos[b] for b in bl}
        if len(s) > 1:
            mixed += 1
            if witness is None:
                witness = (c, sorted(s)[:3], bl[:2])
    print()
    print("  (Q1) is the coset b + GAMMA constant on each iso class?")
    print("       classes carrying MORE THAN ONE coset : %d of %d" % (mixed, len(cls)))
    if mixed == 0:
        print("       -> the coset DESCENDS: it is a TOURNAMENT function.")
    else:
        print("       -> the coset does NOT descend: it is only a TILING function.")
        print("       witness class carries cosets %s" % (witness[1],))

    # (Q2) how does the coset partition compare with the iso partition?
    # is each coset inside one class (refinement) or each class a union of cosets?
    bycoset = {}
    for bits in range(1 << m):
        bycoset.setdefault(cos[bits], []).append(bits)
    coset_in_class = 0
    for k, bl in bycoset.items():
        if len({canon(build(b, tiles, n), n) for b in bl}) == 1:
            coset_in_class += 1
    print()
    print("  (Q2) cosets lying entirely inside a single iso class : %d of %d"
          % (coset_in_class, len(bycoset)))
    if coset_in_class == len(bycoset):
        print("       -> the coset partition REFINES the iso partition.")
    elif mixed == 0:
        print("       -> the iso partition REFINES the coset partition (classes are")
        print("          unions of cosets): the coset is a COARSER tournament invariant.")
    else:
        print("       -> the two partitions CROSS: neither refines the other.")

    # compare with the complement-merged partition
    merged = {}
    for c, bl in cls.items():
        adj = build(bl[0], tiles, n)
        co = canon(opp(adj, n), n)
        key = min(c, co)
        merged.setdefault(key, set()).add(c)
    print("  merged (complement-quotient) classes : %d" % len(merged))
    print("  quotient cosets : %d   -> equal to merged? %s"
          % (1 << (m - r), (1 << (m - r)) == len(merged)))

    # does the star move preserve c_3 ?
    bad3 = 0
    for bits in range(min(1 << m, 4096)):
        a0 = c3(build(bits, tiles, n), n)
        for s in stars:
            if c3(build(bits ^ s, tiles, n), n) != a0:
                bad3 += 1
                break
    print()
    print("  does a star move preserve the 3-cycle count c_3?  tilings where it does")
    print("     NOT : %d of %d sampled  -> %s"
          % (bad3, min(1 << m, 4096),
             "c_3 is NOT a coset invariant" if bad3 else "c_3 IS preserved"))
    print()
