#!/usr/bin/env python3
"""
overlap_moment_matrix_n6.py  (mac-mini 2026-06-15)

FRESH, self-contained investigation. No project imports.

Question (CREATIVE prompt):
  Can a flag-algebra-style moment matrix built from OVERLAP/conflict densities
  (NOT the skew spectrum) CUT a baby-Hodge hole?

Plan:
  (0) Enumerate all 56 unlabeled tournaments on n=6 via gentourng (upper-tri ascii).
  (1) For each tournament compute:
        - score sequence (out-degrees), sorted
        - c3 = number of 3-cycles (cyclic triangles)
        - c5 = number of 5-cycles (directed/cyclic 5-cycles, as subgraphs)
        - p33 = number of UNORDERED pairs of cyclic triangles that SHARE >=1 vertex
                (intersecting triangle pairs = the OCF conflict / overlap carrier)
        - alpha_2 = number of UNORDERED pairs of cyclic triangles that are
                    VERTEX-DISJOINT (the OCF independent-set / disjoint odd-cycle
                    pairs structure). [n=6 => 3+3 partition possible]
        - H(T) = number of independent sets in the odd-cycle conflict graph,
                 i.e. the Redei-Berge / OCF count; computed directly via the
                 product-over-Hamiltonian-path identity is overkill -- instead we
                 use the canonical OCF count H = #{ collections of vertex-disjoint
                 odd cycles } weighted appropriately. We compute H as the value
                 of the OCF: H(T) = sum over sets S of vertex-disjoint odd cycles
                 of 2^{|S|} ... BUT to stay honest & self-contained we instead
                 compute H via the determinant-free combinatorial definition that
                 we can verify on small n: H(T) = #Hamiltonian-path-consistent ...
                 We use the cleanest verifiable surrogate available without project
                 code: see compute_H below (documented).
  (2) Restrict to the c3 == 8 fiber.  At n=6 the regular score (2,2,2,3,3,3)?? --
      NOTE: n=6 has even score sum; "regular" tournament needs odd n. The prompt
      says regular score (2,2,2,3,3,3) which is the near-regular score for n=6.
      We tabulate (c5, p33, alpha_2, H) for that fiber and show co-variation.
  (3) Build 2x2 and 3x3 moment (Gram) matrices from centered overlap densities and
      test PSD-cutting of any c5 value.
"""

import subprocess, itertools, sys
from fractions import Fraction

N = 6
GENTOURNG = "/opt/homebrew/bin/gentourng"

# ---------------------------------------------------------------------------
# (0) Enumerate tournaments. Upper-tri ascii: pairs (i,j), i<j, row by row.
#     Order produced by gentourng default ascii: for i in 0..n-1: for j in i+1..n-1.
#     bit '1' at (i,j) => arc i->j (i beats j). '0' => arc j->i.
# ---------------------------------------------------------------------------
def load_tournaments(n):
    out = subprocess.run([GENTOURNG, "-q", str(n)],
                         capture_output=True, text=True)
    tours = []
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    for line in out.stdout.split():
        line = line.strip()
        if len(line) != len(pairs):
            continue
        # adjacency: beats[i][j] = True if i->j
        beats = [[False]*n for _ in range(n)]
        for bit, (i, j) in zip(line, pairs):
            if bit == '1':
                beats[i][j] = True
            else:
                beats[j][i] = True
        tours.append(beats)
    return tours

def out_degrees(beats, n):
    return [sum(1 for j in range(n) if beats[i][j]) for i in range(n)]

def is_cyclic_triple(beats, a, b, c):
    # a triple of 3 vertices induces either a transitive triple or a 3-cycle.
    # 3-cycle iff each vertex has out-degree exactly 1 within the triple.
    verts = (a, b, c)
    deg = {}
    for x in verts:
        deg[x] = sum(1 for y in verts if y != x and beats[x][y])
    return all(deg[x] == 1 for x in verts)

def cyclic_triangles(beats, n):
    tris = []
    for a, b, c in itertools.combinations(range(n), 3):
        if is_cyclic_triple(beats, a, b, c):
            tris.append(frozenset((a, b, c)))
    return tris

def count_c3(tris):
    return len(tris)

def is_cyclic_5cycle_on(beats, verts):
    """
    Count directed Hamiltonian cycles (cyclic 5-cycles) on the induced
    sub-tournament of these 5 vertices. A 'cyclic 5-cycle' as a SUBGRAPH:
    number of directed cycles of length 5 using exactly these 5 vertices.
    Returns the count of such 5-cycles.
    """
    vs = list(verts)
    # count directed 5-cycles: cyclic orderings v0->v1->...->v4->v0 with all arcs present.
    # Each undirected cyclic sequence counted once. Fix v0 = smallest to avoid rotation,
    # divide reflection by checking direction. We'll count directed cycles and not divide
    # by 2 -- but a directed 5-cycle and its reverse are different arc-sets; in a tournament
    # only one direction of any given pair exists, so a given cyclic vertex-order is a valid
    # cycle iff all 5 consecutive arcs exist. Count distinct cyclic vertex sequences
    # (up to rotation, NOT reflection) that are fully present.
    v0 = vs[0]
    rest = vs[1:]
    cnt = 0
    for perm in itertools.permutations(rest):
        seq = [v0] + list(perm)
        ok = True
        for k in range(5):
            x = seq[k]; y = seq[(k+1) % 5]
            if not beats[x][y]:
                ok = False
                break
        if ok:
            cnt += 1
    # each cyclic 5-cycle (as an arc set) counted once per rotation start fixed (v0 fixed),
    # so rotations already removed; reflections give reverse cycle which can't coexist in a
    # tournament. So cnt = number of directed 5-cycles on these 5 vertices.
    return cnt

def count_c5(beats, n):
    total = 0
    for verts in itertools.combinations(range(n), 5):
        total += is_cyclic_5cycle_on(beats, verts)
    return total

def count_p33_alpha2(tris):
    """
    p33     = # unordered pairs of cyclic triangles sharing >=1 vertex (overlap)
    alpha_2 = # unordered pairs of cyclic triangles that are vertex-disjoint
    """
    p33 = 0
    alpha2 = 0
    for t1, t2 in itertools.combinations(tris, 2):
        if t1 & t2:
            p33 += 1
        else:
            alpha2 += 1
    return p33, alpha2

# ---------------------------------------------------------------------------
# H(T): the OCF / Redei-Berge count.
# Canonical OCF (Odd-Cycle Collection Formula) value:
#   H(T) = sum over collections C of pairwise vertex-disjoint odd cycles
#          (including the empty collection) of  2^{|C|}.
# This is the standard "number of independent sets in the odd-cycle conflict
# graph, weighted by 2^{size}" object central to this project. We compute it
# directly and self-contained: enumerate all odd cycles (length 3 and 5 at n=6),
# build the conflict graph (two odd cycles conflict iff they SHARE a vertex),
# and sum 2^{|independent set|} over independent sets of that conflict graph.
# (Each odd cycle as a distinct arc-set node.)
# ---------------------------------------------------------------------------
def all_odd_cycles(beats, n):
    """Return list of odd cycles as (vertex_frozenset) ; for the conflict graph
    only the vertex support matters for disjointness. But two distinct 5-cycles
    on the SAME 5 vertices are different cycles -> different nodes, but they share
    all vertices so they always conflict. We list each directed odd cycle as a
    node, keyed by its vertex support (support is what determines conflict)."""
    nodes = []  # each node: (support_frozenset)
    # triangles
    for a, b, c in itertools.combinations(range(n), 3):
        if is_cyclic_triple(beats, a, b, c):
            nodes.append(frozenset((a, b, c)))
    # 5-cycles
    for verts in itertools.combinations(range(n), 5):
        k = is_cyclic_5cycle_on(beats, verts)
        for _ in range(k):
            nodes.append(frozenset(verts))
    return nodes

def compute_H(beats, n):
    nodes = all_odd_cycles(beats, n)
    m = len(nodes)
    # sum over independent sets (pairwise vertex-disjoint) of 2^{size}
    # DP via recursion with pruning. m small (<= a few dozen).
    total = 0
    def rec(idx, used_mask, size):
        nonlocal total
        total += (1 << size)  # count current independent set (incl empty)
        for k in range(idx, m):
            sup = nodes[k]
            bm = 0
            for v in sup:
                bm |= (1 << v)
            if bm & used_mask:
                continue
            rec(k+1, used_mask | bm, size+1)
    # The above double-counts: standard pattern is to add 2^size at entry for the
    # set built so far, then extend. Empty set contributes 2^0 = 1 once. Correct.
    rec(0, 0, 0)
    return total

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    tours = load_tournaments(N)
    rows = []
    for idx, beats in enumerate(tours):
        scores = sorted(out_degrees(beats, N))
        tris = cyclic_triangles(beats, N)
        c3 = count_c3(tris)
        c5 = count_c5(beats, N)
        p33, alpha2 = count_p33_alpha2(tris)
        H = compute_H(beats, N)
        rows.append(dict(idx=idx, scores=tuple(scores), c3=c3, c5=c5,
                         p33=p33, alpha2=alpha2, H=H))

    # global sanity
    print("=== ALL 56 TOURNAMENTS (n=6) ===")
    print(f"{'idx':>3} {'scores':>18} {'c3':>3} {'c5':>4} {'p33':>4} {'a2':>3} {'H':>6}")
    for r in rows:
        print(f"{r['idx']:>3} {str(r['scores']):>18} {r['c3']:>3} {r['c5']:>4} "
              f"{r['p33']:>4} {r['alpha2']:>3} {r['H']:>6}")

    # c3 distribution
    from collections import Counter
    c3dist = Counter(r['c3'] for r in rows)
    print("\n=== c3 distribution ===")
    for k in sorted(c3dist):
        print(f"  c3={k}: {c3dist[k]} tournaments")

    # the c3 == 8 fiber
    fiber = [r for r in rows if r['c3'] == 8]
    print(f"\n=== c3 == 8 FIBER  ({len(fiber)} tournaments) ===")
    print(f"{'idx':>3} {'scores':>18} {'c5':>4} {'p33':>4} {'alpha2':>6} {'H':>6}")
    for r in sorted(fiber, key=lambda r: (r['c5'], r['p33'])):
        print(f"{r['idx']:>3} {str(r['scores']):>18} {r['c5']:>4} {r['p33']:>4} "
              f"{r['alpha2']:>6} {r['H']:>6}")

    # c5 distribution within the fiber
    c5dist = Counter(r['c5'] for r in fiber)
    print("\n=== c5 distribution within c3==8 fiber ===")
    for k in sorted(c5dist):
        ps = sorted(set((rr['p33'], rr['alpha2']) for rr in fiber if rr['c5'] == k))
        print(f"  c5={k}: {c5dist[k]} tournaments ; (p33,alpha2) realized = {ps}")

    # which c5 values are MISSING (the holes) within the fiber's c5 range
    realized = sorted(set(r['c5'] for r in fiber))
    if realized:
        full = list(range(realized[0], realized[-1]+1))
        holes = [v for v in full if v not in realized]
        print(f"\n  c5 realized in fiber: {realized}")
        print(f"  c5 HOLES in fiber range: {holes}")

    return rows, fiber

if __name__ == "__main__":
    main()
