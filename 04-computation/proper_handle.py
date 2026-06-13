"""
proper_handle.py -- kind-pasteur-2026-03-14-S110e
GET A PROPER HANDLE ON EVERYTHING

The complete picture:
H_hat(S) = (1/2^{n-1}) * N(S) where N(S) = sum_{P:S sub arcs(P)} (-1)^{des_S(P)}
E_{2k} = sum_{|S|=2k} H_hat(S)^2 = (1/2^{2(n-1)}) * sum_{|S|=2k} N(S)^2

N(S) != 0 iff S can be the edge set of 2k consecutive pairs in some permutation.
Equivalently: S is a union of DISJOINT PATH FRAGMENTS that can be embedded
into a Hamiltonian path.

A "path fragment" = a simple path (sequence of edges forming a trail).
2k edges forming path fragments = a DISJOINT UNION OF PATHS covering
some subset of vertices.

For each such S: |N(S)| = 2 * (number of Ham paths containing S, mod reversal,
with consistent signs).

THE FORMULA: E_{2k}/E_0 = 2*(n-2k)^k / P(n,2k).

PLAN:
1. Characterize which 2k-edge subsets S have N(S) != 0
2. Compute |N(S)| for each type
3. Sum N(S)^2 over all nonzero S
4. Show the sum equals 2*(n-2k)^k * (n!)^2 / (P(n,2k) * 2^{2(n-1)}) ... wait, simpler:
   sum N(S)^2 = 2*(n-2k)^k * (n!)^2 / P(n,2k)
"""

import sys, math
from itertools import permutations, combinations
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding='utf-8')

def main():
    print("=" * 70)
    print("PROPER HANDLE ON EVERYTHING")
    print("kind-pasteur-2026-03-14-S110e")
    print("=" * 70)

    # ============================================================
    # STEP 1: What makes N(S) nonzero?
    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 1: CHARACTERIZE NONZERO N(S)")
    print(f"{'='*70}")

    print(f"""
  S is a set of 2k UNDIRECTED edges from K_n.
  N(S) != 0 iff there exists a Hamiltonian path P of K_n
  such that all edges of S are consecutive pairs in P.

  Equivalently: the subgraph S must be a UNION OF DISJOINT PATHS
  (a "linear forest"), because consecutive pairs in a permutation
  form paths. And these paths must be embeddable into a single
  Hamiltonian path (which they always are, by inserting free vertices).

  So: N(S) != 0 iff S forms a LINEAR FOREST (disjoint union of paths).

  A linear forest on 2k edges and v vertices has:
  - v = 2k + (number of components): each path of length L has L+1 vertices.
  - Components: paths of lengths L_1, ..., L_r with sum L_i = 2k.
  - Number of components r: between 1 (single path of length 2k) and 2k (matching).
  - Vertices: v = 2k + r.
  - Free vertices: n - v = n - 2k - r.

  For a MATCHING (r = 2k, all paths of length 1):
  v = 2k + 2k = 4k. Free = n - 4k. Need n >= 4k.

  For a SINGLE PATH (r = 1):
  v = 2k + 1. Free = n - 2k - 1.
    """)

    # ============================================================
    # STEP 2: Compute |N(S)| for each type at level 4 (2k=4)
    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 2: |N(S)| FOR EACH LINEAR FOREST TYPE AT LEVEL 4")
    print(f"{'='*70}")

    # Level 4: 2k=4 edges. Possible linear forests:
    # - Single path of length 4: r=1, v=5, free=n-5.
    # - Path of length 3 + isolated edge: r=2, v=3+2=5+... wait
    #   Path L=3 has 4 vertices, edge L=1 has 2 vertices.
    #   Total v = 4+2 = 6 if disjoint. free = n-6.
    # - Path of length 2 + path of length 2: r=2, v=3+3=6 if disjoint. free=n-6.
    # - Path of length 2 + two isolated edges: r=3, v=3+2+2=7. Need n>=7.
    # - Matching of 4 edges: r=4, v=8. Need n>=8.
    # etc.

    # At n=5: only single path (r=1, v=5, free=0). 60 subsets. |N|=2.
    # At n=6: single path (r=1, v=5, free=1) + two length-2 paths (r=2, v=6, free=0).
    #         360 + 90 = 450 subsets.
    # At n=7: single path (r=1, v=5, free=2) + {path3+edge1 (v=6,free=1), 2xpath2 (v=6,free=1)}
    #         + ... more types.

    # Let me compute at n=7 to see the full picture.
    for n in [5, 6, 7]:
        arcs = [(i,j) for i in range(n) for j in range(i+1,n)]
        m = len(arcs)
        all_perms = list(permutations(range(n)))

        # Classify all nonzero 4-edge subsets
        type_counts = Counter()
        type_N2_sum = Counter()

        for s_idx in combinations(range(m), 4):
            s_set = set(s_idx)
            arc_list = [arcs[e] for e in s_idx]

            # Check if linear forest
            degs = [0]*n
            adj_graph = defaultdict(set)
            for a, b in arc_list:
                degs[a] += 1; degs[b] += 1
                adj_graph[a].add(b); adj_graph[b].add(a)

            # Linear forest: max degree <= 2, no cycles
            if max(degs) > 2:
                continue  # Not a linear forest (has vertex of degree 3+)

            # Check for cycles: connected components with |edges| = |vertices| means cycle
            vertices = set()
            for a, b in arc_list:
                vertices.update([a, b])

            visited = set()
            components = []
            for v in vertices:
                if v in visited: continue
                comp_v = set(); comp_e = 0
                stack = [v]
                while stack:
                    u = stack.pop()
                    if u in visited: continue
                    visited.add(u); comp_v.add(u)
                    for w in adj_graph[u]:
                        if w not in visited:
                            stack.append(w)
                for a, b in arc_list:
                    if a in comp_v and b in comp_v:
                        comp_e += 1
                if comp_e >= len(comp_v):  # cycle!
                    break
                components.append((len(comp_v), comp_e))
            else:
                # All components are paths. Classify.
                comp_lengths = sorted([e for v, e in components], reverse=True)
                forest_type = tuple(comp_lengths)
                n_components = len(components)
                n_vertices = len(vertices)
                free = n - n_vertices

                # Compute N(S)
                signed = 0
                for perm in all_perms:
                    arc_sgn = {}
                    for i in range(n-1):
                        a, b = perm[i], perm[i+1]
                        if a < b: arc_sgn[arcs.index((a,b))] = +1
                        else: arc_sgn[arcs.index((b,a))] = -1
                    if s_set.issubset(arc_sgn.keys()):
                        sign = 1
                        for e in s_idx: sign *= arc_sgn[e]
                        signed += sign

                if signed != 0:
                    key = (forest_type, free)
                    type_counts[key] += 1
                    type_N2_sum[key] += signed**2
                continue

        print(f"\n  n={n}: Nonzero level-4 subsets by (forest_type, free):")
        total_N2 = 0
        for key in sorted(type_counts.keys()):
            forest, free = key
            count = type_counts[key]
            N2_sum = type_N2_sum[key]
            avg_N2 = N2_sum // count
            total_N2 += N2_sum
            print(f"    forest={forest}, free={free}: {count} subsets, "
                  f"|N|^2={avg_N2} each, total N^2={N2_sum}")

        # E_4/E_0 from this
        E0 = (math.factorial(n) / 2**(n-1))**2
        E4 = total_N2 / 2**(2*(n-1))
        ratio = E4 / E0
        pred = 2*(n-4)**2 / math.perm(n, 4) if n > 4 else 0
        print(f"    Total sum N^2 = {total_N2}")
        print(f"    E_4/E_0 = {ratio:.10f}, formula = {pred:.10f}, "
              f"match = {abs(ratio - pred) < 1e-8}")

    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 3: THE PATTERN IN |N(S)|^2")
    print(f"{'='*70}")

    # From the data:
    # n=5: (4,), free=0: 60 subsets, |N|^2=4 each. Total = 240.
    # n=6: (4,), free=1: 360, |N|^2=16. Total = 5760.
    #       (2,2), free=0: 90, |N|^2=64. Total = 5760.
    # n=7: (4,), free=2: ?, |N|^2=? ...

    print(f"""
  |N(S)|^2 for a LINEAR FOREST of type (L_1, ..., L_r) with f free vertices:

  Path reversal gives factor 2: |N(S)| = 2 * M(S).
  So |N(S)|^2 = 4 * M(S)^2.

  M(S) = number of INDEPENDENT Ham paths (mod reversal) containing S.

  For a single path of length 4 (type (4,)) with f free vertices:
  The path uses 5 vertices. The f free vertices must be inserted.
  Each free vertex can be inserted at any of the existing positions
  in the Hamiltonian path (but must avoid creating new path fragments).
  Actually: a Hamiltonian path containing the length-4 subpath
  just needs to INSERT the f free vertices into the path.
  Each insertion is at one of (5 + 0, 5+1, ...) positions.
  Wait: the original path has 5 vertices in order. To extend to n=5+f:
  insert f vertices. Each can go at any of 6 positions (before, between, after).
  But the positions are DISTINGUISHABLE: insert at position 0, 1, ..., 5.
  That's 6^f? No, each vertex is distinct.

  Actually: the Ham path must visit all n vertices.
  It already has the subpath a0-a1-a2-a3-a4 as consecutive.
  The remaining f vertices must be placed SOMEWHERE ELSE in the path.
  They go either before a0, after a4, or... no, the 4 edges must be
  CONSECUTIVE. So the remaining vertices are NOT between the subpath vertices.
  They're before the subpath, after it, or... wait, the subpath occupies
  5 consecutive positions in the Ham path. The other f vertices fill
  the remaining positions (before or after the subpath block).

  Ham path of length n-1: _ _ ... _ [a0 a1 a2 a3 a4] _ _ ... _
  The subpath block is somewhere in the middle (or at the ends).
  The f free vertices fill the non-block positions.

  Number of ways: Choose WHICH f positions (out of n-1 total path positions)
  are NOT in the block. The block has 5 consecutive positions.
  So the block starts at position p, 0 <= p <= n-5 = f.
  For each p: the f free positions are {0,...,p-1, p+5,...,n-1}.
  The f free vertices can go in these positions in f! ways.
  Total: (f+1) * f! = (f+1)!.

  But wait: each such Ham path uses the SAME 4 edges (the subpath).
  M(S) = number of such paths / (accounting for the subpath being fixed).
  Actually M(S) = number of permutations (mod reversal) with the block fixed.
  = (f+1)! / 2? No. Let me just compute.

  At n=5 (f=0): (0+1)! = 1. M = 1. |N| = 2*1 = 2. |N|^2 = 4. CHECK!
  At n=6 (f=1): (1+1)! = 2. M = 2. |N| = 2*2 = 4. |N|^2 = 16. CHECK!
  At n=7 (f=2): (2+1)! = 6. M = 6. |N| = 2*6 = 12? Or...

  Wait: M should be the number of INDEPENDENT paths (mod reversal).
  A reversal of the whole path also reverses the subpath block.
  For even |S|: reversal preserves the sign. So P and P^rev give same N.
  But they give the SAME Ham path (undirected). So we shouldn't double-count.
  M(S) = total directed paths / 2 = ?

  At n=5: 2 directed paths (path and reverse). M = 1. |N| = 2. CHECK.
  At n=6: (f+1)! = 2 positions for block * 1! for free vertex = 2 arrangements.
  Each arrangement gives 2 directed paths. Total directed = 4. M = 2. |N| = 4. CHECK.
  At n=7: (f+1)! = 6 arrangements. Total directed = 12. M = 6. |N| = 12. |N|^2 = 144.

  Hmm: but actually each arrangement might give DIFFERENT directed paths
  (not related by global reversal). The reversal flips the WHOLE path,
  not just the block. So some arrangements are related by reversal.

  Let me just check at n=7 what |N| is for type (4,) with free=2.
    """)

    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 4: COMPUTE |N(S)| FOR TYPE (4,) AT n=7")
    print(f"{'='*70}")

    n = 7
    arcs = [(i,j) for i in range(n) for j in range(i+1,n)]
    all_perms = list(permutations(range(n)))

    # Find ONE type-(4,) subset and compute N
    # A single path of length 4: say 0-1-2-3-4 (edges 01,12,23,34).
    test_edges = [(0,1), (1,2), (2,3), (3,4)]
    test_S = set(arcs.index(e) for e in test_edges)

    signed = 0
    path_count = 0
    for perm in all_perms:
        arc_sgn = {}
        for i in range(n-1):
            a, b = perm[i], perm[i+1]
            if a < b: arc_sgn[arcs.index((a,b))] = +1
            else: arc_sgn[arcs.index((b,a))] = -1
        if test_S.issubset(arc_sgn.keys()):
            sign = 1
            for e in test_S: sign *= arc_sgn[e]
            signed += sign
            path_count += 1

    print(f"  S = edges of path 0-1-2-3-4 at n=7:")
    print(f"  N(S) = {signed}, paths containing S = {path_count}")
    print(f"  |N(S)|^2 = {signed**2}")

    # At n=7: f=2 free vertices (5 and 6).
    # Arrangements: the block 0-1-2-3-4 in 3 positions (start at 0, 1, 2 of 7-slot path).
    # For each position, 2 free vertices fill the remaining 2 slots in 2! = 2 ways.
    # Total: 3 * 2 = 6 directed paths. But reversal pairs: 6/2 = 3 independent.
    # M = 3? But then |N| = 2*3 = 6? Or |N| = signed count which depends on signs.

    # Let me look at paths containing S and their signs.
    print(f"\n  Paths containing S = {{(0,1),(1,2),(2,3),(3,4)}} at n=7:")
    for perm in all_perms:
        arc_sgn = {}
        for i in range(n-1):
            a, b = perm[i], perm[i+1]
            if a < b: arc_sgn[arcs.index((a,b))] = +1
            else: arc_sgn[arcs.index((b,a))] = -1
        if test_S.issubset(arc_sgn.keys()):
            sign = 1
            for e in test_S: sign *= arc_sgn[e]
            print(f"    P={perm}, sign_product={sign}")

    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 5: FORMULA FOR |N(S)| FROM PATH COUNTING")
    print(f"{'='*70}")

    # From the paths listed above, I can see the structure.
    # The paths must contain 0-1-2-3-4 as consecutive elements.
    # The other elements (5, 6) go before or after the block.
    # Each gives a specific sign_product.

    # The sign_product for the block: it depends on whether the block
    # is traversed in ascending (0->1->2->3->4) or descending (4->3->2->1->0) order.
    # Ascending: all signs +1. Product = +1.
    # Descending: all signs -1. Product = (-1)^4 = +1.
    # So the block always contributes +1 regardless of direction.
    # The free vertices don't affect the sign (they're not in S).
    # Therefore: EVERY path containing S contributes +1 to N(S).
    # N(S) = (number of paths containing S) = (total directed paths through S).

    # For a single path block of length L (L+1 vertices):
    # The block occupies L+1 consecutive positions in the Ham path.
    # The remaining f = n - L - 1 vertices fill the other positions.
    # The block can start at positions 0, 1, ..., f.
    # For each starting position: the f free vertices fill the remaining slots.
    # There are f! arrangements of the free vertices.
    # Total directed paths through S: (f+1) * f! = (f+1)!
    # And N(S) = (f+1)! (all contribute +1).
    # |N(S)|^2 = ((f+1)!)^2.

    # VERIFY at n=5: f=0. |N|=(0+1)!=1. But we found |N|=2.
    # PROBLEM: I said all signs are +1, but the REVERSE path also uses S.
    # 0-1-2-3-4 (ascending) -> sign product = +1^4 = +1.
    # 4-3-2-1-0 (descending) -> sign product = (-1)^4 = +1.
    # Both contribute +1. Total = 2. And (f+1)! = 1! = 1. But N = 2.
    # So the formula should be (f+1)! * 2? No: both the ascending and
    # descending traversals are DIFFERENT directed paths.
    # (f+1) positions * f! arrangements * 2 (for reversal) = 2 * (f+1)!.
    # N(S) = 2 * (f+1)! (since sign is always +1).
    # |N|^2 = 4 * ((f+1)!)^2.

    # n=5: |N|^2 = 4*1 = 4. CHECK!
    # n=6: |N|^2 = 4*4 = 16. CHECK!
    # n=7: |N|^2 = 4*36 = 144. Let me verify.

    print(f"  Single-path-block formula: |N(S)|^2 = 4 * ((f+1)!)^2")
    print(f"    n=5, f=0: 4*1 = 4. Matches |N|=2.")
    print(f"    n=6, f=1: 4*4 = 16. Matches |N|=4.")
    print(f"    n=7, f=2: 4*36 = 144. Need to verify: |N| should be 12.")
    print(f"    Computed |N(S)| = {abs(signed)}, |N|^2 = {signed**2}")
    print(f"    Match: {signed**2 == 4 * math.factorial(3)**2}")

    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 6: COUNT NONZERO SUBSETS OF EACH TYPE")
    print(f"{'='*70}")

    # Type (4,): single path of length 4 (using 5 specific vertices).
    # Choose 5 vertices from n: C(n, 5).
    # Number of Hamiltonian paths on those 5 vertices: 5!/2 = 60.
    # But as UNDIRECTED EDGE SETS: each set of 4 edges = 1 Ham path edge set.
    # So: C(n, 5) * 60 undirected edge sets... wait.
    # Actually: number of Ham path edge sets on K_5 = 60 (we computed this).
    # For each 5-vertex subset: 60 edge sets.
    # Total type-(4,) subsets: C(n, 5) * 60.
    # But wait: at n=5, C(5,5)*60 = 60. CHECK.
    # At n=6: C(6,5)*60 = 6*60 = 360. CHECK!
    # At n=7: C(7,5)*60 = 21*60 = 1260.

    # Type (2,2): two disjoint paths of length 2 (each using 3 vertices).
    # Choose 6 vertices (3+3, disjoint): C(n, 3) * C(n-3, 3) / 2
    #   (divide by 2 since the two paths are unordered).
    # For each pair of 3-vertex sets: paths of length 2 on 3 vertices.
    #   Number of Ham path edge sets on K_3 = 3 (= 3!/2).
    # So: total type-(2,2) subsets = C(n,3)*C(n-3,3)/2 * 3 * 3.
    # Wait: each 3-vertex set has 3 path edge sets. Two independent sets.
    # Total = [C(n,3)*C(n-3,3)/2] * 3 * 3 = C(n,3)*C(n-3,3)/2 * 9.
    # At n=6: C(6,3)*C(3,3)/2 * 9 = 20*1/2 * 9 = 90. CHECK!!!

    print(f"  Type (4,) count: C(n,5) * 60")
    for n in [5, 6, 7, 8]:
        count = math.comb(n, 5) * 60
        print(f"    n={n}: {count}")

    print(f"\n  Type (2,2) count: C(n,3)*C(n-3,3)/2 * 9")
    for n in [6, 7, 8]:
        count = math.comb(n, 3) * math.comb(n-3, 3) // 2 * 9
        print(f"    n={n}: {count}")

    # Type (3,1): path of length 3 + isolated edge.
    # Choose 4 vertices for path: C(n,4). Path edge sets on K_4: 4!/2 = 12.
    # Wait: Hamiltonian paths on 4 vertices = 4!/2 = 12 directed, 12/2 = 6 undirected?
    # Hmm: 4! = 24 perms. Each gives a directed Ham path.
    # Undirected edge sets: 24/2 = 12? No.
    # A path on 4 vertices has 3 edges. 4 vertices, 3 edges.
    # Number of SPANNING paths (edge sets) on K_4:
    # Each perm gives an edge set. {a,b,c,d} -> {ab,bc,cd}, {a,c,b,d} -> {ac,cb,bd}=different.
    # Total undirected = 4!/2 = 12.
    # Then choose 2 MORE vertices for the isolated edge: C(n-4, 2).
    # The isolated edge: 1 edge on the 2 vertices (just 1 way).
    # Total: C(n,4) * 12 * C(n-4, 2) * 1.
    # But: does the isolated edge need to be DISJOINT from the path? YES (disjoint paths).
    # And: the path+edge must use 4+2 = 6 vertices total.
    # At n=7: C(7,4)*12*C(3,2) = 35*12*3 = 1260.

    print(f"\n  Type (3,1) count: C(n,4)*12*C(n-4,2)")
    for n in [7, 8]:
        count = math.comb(n, 4) * 12 * math.comb(n-4, 2)
        print(f"    n={n}: {count}")

    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 7: TOTAL SUM N^2 AT LEVEL 4")
    print(f"{'='*70}")

    # At n=7, level 4:
    # Type (4,): 1260 subsets, |N|^2 = 4*(f+1)!^2 = 4*6 = ... wait.
    # f = n - 5 = 2. |N|^2 = 4*(2+1)!^2 = 4*36 = 144 (if my formula is right).
    # BUT: does the sign product depend on the specific edge set?
    # For type (4,): single path block. Sign = +1 always (since |S|=4 even).
    # So |N| = number of directed paths through S.
    # For 5 specific vertices forming a path, with f=2 free vertices:
    # Directed paths = 2*(f+1)! = 2*6 = 12. |N| = 12. |N|^2 = 144.

    # Type (2,2): two disjoint length-2 paths on 6 vertices.
    # f = n - 6 = 1 free vertex at n=7.
    # How many directed paths go through BOTH 2-edge subpaths?
    # Each 2-edge subpath: a-b-c (3 vertices, 2 edges).
    # The two subpaths use 6 vertices, leaving 1 free vertex.
    # A Ham path must contain both a1-b1-c1 and a2-b2-c2 as subpaths.
    # These are 2 blocks of 3 consecutive vertices each, plus 1 free vertex.
    # Arrangements: interleave the 2 blocks and 1 free vertex.
    # Possible orderings: [block1, block2, free], [block1, free, block2],
    #   [block2, block1, free], [block2, free, block1],
    #   [free, block1, block2], [free, block2, block1].
    #   But each block can be traversed forward or backward.
    #   Forward: a-b-c. Backward: c-b-a.
    #   4 combinations * 6 orderings... hmm.
    # Wait: the UNDIRECTED edge set is fixed. We count DIRECTED paths.
    # Each block has 2 directions. 2 blocks = 4 direction combos.
    # For each direction combo: (2+1) positions for blocks and free = 3! = 6 arrangements.
    # But some arrangements may repeat.
    # Total directed: 4 * arrangements (accounting for the free vertex position).
    # Actually: 3 items to arrange (block1, block2, free). 3! = 6 orderings.
    # Each block has 2 directions: 2^2 = 4.
    # Total directed: 6 * 4 = 24.
    # But global reversal: each directed path paired with its reverse.
    # 24/2 = 12 independent. M(S) = 12.
    # |N(S)| = 2*12 = 24... hmm that seems too high.

    # Let me check: at n=6, type (2,2) with f=0: |N|=8.
    # f=0: no free vertices. 2 blocks to arrange. 2! = 2 orderings.
    # 2 direction combos per block: 2^2 = 4.
    # Total directed: 2 * 4 = 8. |N| = 8. CHECK!
    # All contribute +1 (sign product = (-1)^0 * ... = +1 since each block has even length).
    # Wait: each block has length 2 (even), so (-1)^{des} for 2 descents = (-1)^2 = +1.
    # And (-1)^0 = +1 for 0 descents. So sign = +1 always. N = 8. CHECK!

    print(f"  Type (2,2) at n=6, f=0:")
    print(f"    Directed paths = 2! * 2^2 = 8. N = 8. |N|^2 = 64. CHECK!")

    print(f"\n  Type (2,2) at n=7, f=1:")
    print(f"    3 items (block1, block2, free). 3! = 6 arrangements.")
    print(f"    2 directions per block: 2^2 = 4.")
    print(f"    Total directed = 6 * 4 = 24. |N| = 24. |N|^2 = 576.")

    # Let me verify at n=7 with a specific type (2,2) subset.
    # S = {(0,1),(1,2),(3,4),(4,5)}: two paths 0-1-2 and 3-4-5.
    test_S2 = set(arcs.index(e) for e in [(0,1),(1,2),(3,4),(4,5)])
    signed2 = 0
    for perm in all_perms:
        arc_sgn = {}
        for i in range(n-1):
            a, b = perm[i], perm[i+1]
            if a < b: arc_sgn[arcs.index((a,b))] = +1
            else: arc_sgn[arcs.index((b,a))] = -1
        if test_S2.issubset(arc_sgn.keys()):
            sign = 1
            for e in test_S2: sign *= arc_sgn[e]
            signed2 += sign
    print(f"    Computed |N| = {abs(signed2)}, |N|^2 = {signed2**2}")
    print(f"    Matches 24? {abs(signed2) == 24}")

    # Now: total sum N^2 at n=7, level 4:
    # Type (4,): 1260 * 144 = 181440
    # Type (2,2): 315 * 576 = 181440 (if count=315 and |N|^2=576)
    # Type (3,1): 1260 * ??? ... need to compute.

    # Wait: I need to check type (3,1) too.
    # Path of length 3 on 4 vertices + edge on 2 vertices. Total 6 vertices.
    # f = n-6 = 1 at n=7.
    # Items: block_3 (4 vertices), block_1 (2 vertices), free (1 vertex). 3 items.
    # Orderings: 3! = 6.
    # Block_3 has 2 directions. Block_1 has 2 directions. Total = 2^2 = 4.
    # Directed paths: 6 * 4 = 24. |N| = 24. |N|^2 = 576.

    test_S3 = set(arcs.index(e) for e in [(0,1),(1,2),(2,3),(4,5)])
    signed3 = 0
    for perm in all_perms:
        arc_sgn = {}
        for i in range(n-1):
            a, b = perm[i], perm[i+1]
            if a < b: arc_sgn[arcs.index((a,b))] = +1
            else: arc_sgn[arcs.index((b,a))] = -1
        if test_S3.issubset(arc_sgn.keys()):
            sign = 1
            for e in test_S3: sign *= arc_sgn[e]
            signed3 += sign
    print(f"\n  Type (3,1) at n=7: |N| = {abs(signed3)}, |N|^2 = {signed3**2}")

    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 8: THE GENERAL FORMULA")
    print(f"{'='*70}")

    # For a linear forest with r components of lengths L_1,...,L_r:
    # Total vertices: v = sum(L_i + 1) = sum(L_i) + r = 2k + r.
    # Free vertices: f = n - v = n - 2k - r.
    # Number of directed paths through S: (r + f)! * 2^r / ???
    # Wait: items to arrange = r blocks + f free vertices = r + f items.
    # Orderings: (r+f)!.
    # Each of r blocks has 2 directions: 2^r.
    # Total directed: (r+f)! * 2^r.
    # Sign: each block of length L_i contributes (-1)^{0 or L_i} = +1
    #   (since L_i is always positive and the sign product for even |S|
    #    is +1 for both directions of each block... wait, each block's
    #    sign product = (-1)^{descents in block} = (-1)^0 = +1 (ascending)
    #    or (-1)^{L_i} (all descending). For the block edges: the UNDIRECTED
    #    edges are fixed. Ascending: all +1. Descending: all -1.
    #    Product for block of length L_i: (+1)^{L_i} = +1 or (-1)^{L_i}.
    #    Since we're summing over BOTH directions of each block:
    #    contribution per block = (+1)^{L_i} + (-1)^{L_i} = 1 + (-1)^{L_i}.
    #    For EVEN L_i: 1 + 1 = 2.
    #    For ODD L_i: 1 + (-1) = 0!

    # WAIT: if L_i is ODD, the two directions CANCEL.
    # A block of ODD length has sign_ascending = +1 and sign_descending = -1.
    # Their sum = 0. So the TOTAL signed count = 0!

    # But at level 4: L_1+...+L_r = 2k = 4. All L_i are positive.
    # Possible partitions of 4: (4), (3,1), (2,2), (2,1,1), (1,1,1,1).
    # For (3,1): L_1=3 (odd), L_2=1 (odd). Both odd -> both contribute 0.
    # Product = 0. N(S) = 0!

    # FOR (4): L_1=4 (even). Contributes 2. N = (f+1)! * 2. |N|^2 = 4*(f+1)!^2.
    # FOR (2,2): L_1=L_2=2 (both even). Each contributes 2. N = (f+r)! * 4.
    #   r=2, so N = (f+2)! * 4 / ??? Wait: directed paths = (r+f)! * 2^r = (2+f)! * 4.
    #   But the sign product = +1 for EACH direction combo (even blocks).
    #   So N = (f+2)! * 4 (ALL contribute +1).

    # FOR (3,1): L_1=3 (odd), L_2=1 (odd). The sign product depends on directions.
    #   Block 1 ascending: +1. Descending: (-1)^3 = -1.
    #   Block 2 ascending: +1. Descending: (-1)^1 = -1.
    #   Combos: (asc,asc)=+1*+1=+1. (asc,desc)=+1*(-1)=-1.
    #           (desc,asc)=(-1)*+1=-1. (desc,desc)=(-1)*(-1)=+1.
    #   Sum: +1-1-1+1 = 0. N = 0!

    # FOR (2,1,1): L_1=2 (even), L_2=L_3=1 (odd).
    #   Block 1: 2 directions, both contribute same sign (+1). Factor 2.
    #   Block 2: 2 directions, signs +1 and -1. Cancel. Factor 0.
    #   Product = 0. N = 0.

    # FOR (1,1,1,1): all L=1 (odd). Each block cancels. N = 0.

    print(f"""
  *** THE SIGN CANCELLATION RULE ***

  For a linear forest with component lengths L_1, ..., L_r:
  The signed count N(S) is ZERO unless ALL L_i are EVEN.

  This is because:
  - Each block of EVEN length L_i: ascending and descending both give
    sign product = +1. Both directions contribute. Factor = 2.
  - Each block of ODD length L_i: ascending gives +1, descending gives -1.
    They CANCEL. Factor = 0.

  Since N(S) = 0 whenever ANY L_i is odd:
  ONLY forests with ALL EVEN component lengths contribute.

  Partitions of 2k into even parts:
  2k = 2*a_1 + 2*a_2 + ... + 2*a_r (each a_i >= 1).
  Equivalently: a_1 + a_2 + ... + a_r = k.
  This is a PARTITION OF k!

  At level 4 (k=2): partitions of 2: (2) and (1,1).
    (2) -> one block of length 4. r=1.
    (1,1) -> two blocks of length 2. r=2.
    These are the ONLY nonzero types! CONFIRMED by computation!

  At level 6 (k=3): partitions of 3: (3), (2,1), (1,1,1).
    (3) -> one block of length 6. r=1.
    (2,1) -> blocks of lengths 4 and 2. r=2.
    (1,1,1) -> three blocks of length 2. r=3.

  For partition (a_1, ..., a_r) of k (each a_i >= 1):
  - Component lengths: 2*a_i (each even).
  - Total edges: 2k.
  - Number of vertices: 2k + r.
  - Free vertices: f = n - 2k - r.
  - Items to arrange: r blocks + f free = r + f items.
  - Directed paths: (r+f)! * 2^r (all signs +1).
  - |N(S)|^2 = ((r+f)! * 2^r)^2 = 4^r * ((r+f)!)^2.

  Number of such subsets: multinomial counting of vertex assignments.

  TOTAL sum N^2 at level 2k = sum over partitions (a_1,...,a_r) of k:
    (count of subsets of this type) * 4^r * ((n-2k-r+r)!)^2
    = (count) * 4^r * ((n-2k)!)^2 ... wait, items = r+f = r + n-2k-r = n-2k.
    So (r+f)! = (n-2k)!. INDEPENDENT of r!

  THEREFORE: |N(S)|^2 = 4^r * ((n-2k)!)^2 for ALL nonzero types.
  BUT (n-2k)! IS the same for all r! So the sum is:
  sum N^2 = ((n-2k)!)^2 * sum over partitions of k: (count of type) * 4^r.

  AND: E_2k/E_0 = sum N^2 / (2^(2(n-1)) * E_0)
     = ((n-2k)!)^2 * [sum (count)*4^r] / ((n!)^2 / 2^(2(n-1)) * 2^(2(n-1)))
     = ((n-2k)!)^2 * [sum (count)*4^r] / (n!)^2.

  The formula says this = 2*(n-2k)^k / P(n,2k) = 2*(n-2k)^k * (n-2k)! / n!.
  So: (n-2k)! * [sum (count)*4^r] / n! = 2*(n-2k)^k / P(n,2k).
  => [sum (count)*4^r] = 2*(n-2k)^k * n! / (P(n,2k) * (n-2k)!)
     = 2*(n-2k)^k * n! / (n! * (n-2k)! / (n-2k)!)... hmm let me simplify.
  P(n,2k) = n!/(n-2k)!.
  sum (count)*4^r = 2*(n-2k)^k * n! / (n!/(n-2k)!) / (n-2k)!... no.

  Wait: E_2k/E_0 = sum N^2 / (n!)^2 (from the formula I derived).
  sum N^2 = sum (count * |N|^2) = sum (count * 4^r * ((n-2k)!)^2).
  (n-2k)!^2 is constant across types!
  So: sum N^2 = ((n-2k)!)^2 * sum_types count * 4^r.
  E_2k/E_0 = ((n-2k)!)^2 * sum count*4^r / (n!)^2.
  = ((n-2k)!/n!)^2 * sum count*4^r.
  = (1/P(n,2k))^2 * sum count*4^r.
  Wait no: (n-2k)!/n! = 1/P(n, n-(n-2k)) = ... this is getting circular.
  Let me just verify numerically.
  """)

    # Numerical verification
    for n_val in [5, 6, 7]:
        f = n_val - 5  # for type (4,)
        nk = (n_val - 4)  # = n - 2k where k=2
        items = nk  # r + f = r + n-2k-r = n-2k

        # |N|^2 should be 4^r * (n-2k)!^2 for each type.
        # Type (4,): r=1. |N|^2 = 4 * (nk)!^2? But nk = n-4.
        # At n=5: 4 * 1!^2 = 4. CHECK.
        # At n=6: 4 * 2!^2 = 16. CHECK.
        # At n=7: 4 * 3!^2 = 144. Need to verify.

        # Type (2,2): r=2. |N|^2 = 16 * (nk)!^2?
        # At n=6: nk=2. 16 * 2!^2 = 64. CHECK!
        # At n=7: nk=3. 16 * 3!^2 = 576. Need to verify.

        print(f"  n={n_val}, k=2: (n-2k) = {nk}")
        print(f"    Type (4,): r=1. Predicted |N|^2 = 4*{nk}!^2 = {4*math.factorial(nk)**2}")
        print(f"    Type (2,2): r=2. Predicted |N|^2 = 16*{nk}!^2 = {16*math.factorial(nk)**2}")

    # Verify at n=7
    print(f"\n  At n=7: type (4,) |N|^2 = {signed**2}, predicted = {4*math.factorial(3)**2}")
    print(f"  At n=7: type (2,2) |N|^2 = {signed2**2}, predicted = {16*math.factorial(3)**2}")
    print(f"  At n=7: type (3,1) |N|^2 = {signed3**2}, predicted = 0 (odd lengths cancel)")

    print(f"\n{'='*70}")
    print("*** THE FORMULA IS: |N(S)|^2 = 4^r * ((n-2k)!)^2 ***")
    print("*** WHERE r = NUMBER OF EVEN-LENGTH COMPONENTS ***")
    print("*** AND ONLY ALL-EVEN FORESTS CONTRIBUTE ***")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
