#!/usr/bin/env python3
"""
Refined proof that β₁(T) ≤ 1 for tournaments.

Key finding from previous script: the edge-sharing graph of 3-cycles
is NOT always connected. So we need a DIFFERENT approach to show all
3-cycles are homologous.

NEW APPROACH: Work directly with the quotient Z₁/B₁.
Show that B₁ + span(any single 3-cycle) = B₁ + span(all 3-cycles).

The KEY LEMMA (sharing an edge a→b) is proved algebraically.
The question is: can we CHAIN these lemma applications to connect
any two 3-cycles?

Two 3-cycles sharing a vertex v always share a DIRECTED edge:
  C₁ = (...→v→... or ...→?→v→...) uses two edges through v.
  C₂ = (...→v→... or ...→?→v→...) uses two edges through v.
  v has one outgoing and one incoming edge in each 3-cycle.

  C₁ uses edge (v→a) and (b→v).
  C₂ uses edge (v→c) and (d→v).

  If a=c: shared outgoing edge v→a=v→c. Apply KEY LEMMA.
  If b=d: shared incoming edge b→v=d→v. Apply KEY LEMMA.
  If a=d: C₁ has v→a, C₂ has a→v. But tournament has exactly one direction!
    So v→a and a→v can't both exist. CONTRADICTION (since they share vertex v
    AND vertex a=d).
  Wait: if a=d, then C₁ has edge v→a and C₂ has edge d→v = a→v.
  But v→a and a→v: tournament has exactly one. If v→a (from C₁), then
  a→v is impossible. So C₂ can't have a→v. Contradiction.

  Actually: v has out-degree 1 and in-degree 1 within each 3-cycle.
  C₁: v→a₁→a₂→v. Outgoing: v→a₁, Incoming: a₂→v.
  C₂: v→b₁→b₂→v. Outgoing: v→b₁, Incoming: b₂→v.

  Shared vertex from {a₁,a₂} ∩ {b₁,b₂}:
    - a₁=b₁: shared outgoing edge v→a₁. KEY LEMMA applies directly.
    - a₂=b₂: shared incoming edge a₂→v. KEY LEMMA applies directly.
    - a₁=b₂: C₂ has b₂→v = a₁→v, but C₁ has v→a₁. Contradiction.
    - a₂=b₁: C₂ has v→b₁ = v→a₂, but C₁ has a₂→v. Contradiction.

  So: two 3-cycles sharing a vertex ALWAYS share a DIRECTED EDGE!
  And the KEY LEMMA then applies.

  Wait, but the script showed FAILURES for vertex-sharing pairs?!
  Let me check what went wrong.

opus-2026-03-08
"""
import numpy as np
from collections import defaultdict, Counter
from math import comb

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def edges_list(A, n):
    return [(i,j) for i in range(n) for j in range(n) if i != j and A[i][j]]

def build_boundary_2(A, n):
    tt = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[a][c]:
                    tt.append((a,b,c))
    el = edges_list(A, n)
    edge_idx = {e: i for i, e in enumerate(el)}
    M = np.zeros((len(el), len(tt)), dtype=float)
    for j, (a,b,c) in enumerate(tt):
        M[edge_idx[(b,c)], j] += 1
        M[edge_idx[(a,c)], j] -= 1
        M[edge_idx[(a,b)], j] += 1
    return M, tt, el, edge_idx

def cycle_vec(edge_idx, cycle_verts, ne):
    v = np.zeros(ne)
    for i in range(len(cycle_verts)):
        a = cycle_verts[i]
        b = cycle_verts[(i+1) % len(cycle_verts)]
        v[edge_idx[(a,b)]] = 1
    return v

def is_in_image(vec, M):
    if M.shape[1] == 0:
        return np.allclose(vec, 0)
    aug = np.column_stack([M, vec.reshape(-1,1)])
    r1 = np.linalg.matrix_rank(M, tol=1e-10)
    r2 = np.linalg.matrix_rank(aug, tol=1e-10)
    return r1 == r2

print("="*72)
print("REFINED β₁ ≤ 1 PROOF")
print("="*72)

# =================================================================
# VERIFY: Two 3-cycles sharing a vertex ALWAYS share a directed edge
# =================================================================
print("\n--- Do vertex-sharing 3-cycles always share a directed edge? ---")

for n in [3, 4, 5]:
    fail = 0
    total = 0
    for A in all_tournaments(n):
        cycles3 = []
        for i in range(n):
            for j in range(n):
                if j == i or not A[i][j]: continue
                for k in range(n):
                    if k == i or k == j: continue
                    if A[j][k] and A[k][i]:
                        cycles3.append((i,j,k))

        for ci in range(len(cycles3)):
            for cj in range(ci+1, len(cycles3)):
                c1, c2 = cycles3[ci], cycles3[cj]
                if not (set(c1) & set(c2)):
                    continue
                total += 1

                # Check if they share a directed edge
                edges1 = set()
                for idx in range(3):
                    edges1.add((c1[idx], c1[(idx+1)%3]))
                edges2 = set()
                for idx in range(3):
                    edges2.add((c2[idx], c2[(idx+1)%3]))

                shared_edges = edges1 & edges2
                if not shared_edges:
                    fail += 1
                    if fail <= 3:
                        print(f"  n={n}: NO shared edge: C1={c1}, C2={c2}")
                        print(f"    edges1={edges1}")
                        print(f"    edges2={edges2}")
                        print(f"    shared vertices: {set(c1) & set(c2)}")

    print(f"n={n}: {total} vertex-sharing pairs, {fail} without shared directed edge")


# =================================================================
# WHY the earlier test showed "failures"
# =================================================================
print(f"\n\n{'='*72}")
print("DEBUGGING: Why did the earlier vertex-sharing test fail?")
print("="*72)

# The issue: two ORDERED 3-cycles (i,j,k) and (i',j',k') can share
# a vertex but traverse it in opposite directions. As ORIENTED cycles,
# (a→b→c→a) has edges (a,b), (b,c), (c,a). Another cycle using the
# same vertex set but opposite orientation would be (a→c→b→a) with
# edges (a,c), (c,b), (b,a). These share vertices {a,b,c} but NO
# directed edges (all reversed).
#
# But for TOURNAMENTS: if (a→b→c→a) is a 3-cycle, then (a→c→b→a)
# needs a→c. Since c→a in the first cycle, we have c→a, so a→c is
# impossible. So the REVERSE cycle doesn't exist.
#
# What about the EARLIER test: It checked C1-C2 ∈ B_1 where C1 and C2
# share a vertex but may not share an edge. My analysis above shows
# they MUST share an edge if both are directed 3-cycles in the same
# tournament. So what went wrong?

print("\nLet me check the specific failure case from earlier:")
print("C1=(0,4,1), C2=(1,3,2), shared={1}")

# Construct a tournament where both exist
n = 5
for A in all_tournaments(n):
    # Check if (0,4,1) is a 3-cycle: 0→4, 4→1, 1→0
    if not (A[0][4] and A[4][1] and A[1][0]):
        continue
    # Check if (1,3,2) is a 3-cycle: 1→3, 3→2, 2→1
    if not (A[1][3] and A[3][2] and A[2][1]):
        continue

    print(f"\nFound tournament with C1=(0,4,1) and C2=(1,3,2)")
    # C1 edges: (0,4), (4,1), (1,0)
    # C2 edges: (1,3), (3,2), (2,1)
    print(f"  C1 edges: (0→4), (4→1), (1→0)")
    print(f"  C2 edges: (1→3), (3→2), (2→1)")
    print(f"  Shared vertex: 1")
    print(f"  Shared directed edges: NONE")
    print(f"  C1 has edges through 1: (4,1) incoming, (1,0) outgoing")
    print(f"  C2 has edges through 1: (1,3) outgoing, (2,1) incoming")
    print(f"  C1 outgoing from 1: 1→0")
    print(f"  C2 outgoing from 1: 1→3")
    print(f"  DIFFERENT outgoing! C1 incoming to 1: 4→1, C2 incoming: 2→1. DIFFERENT!")
    print(f"  So they share vertex 1 but NO directed edge.")

    # But my analysis said this is impossible! Let me recheck.
    # C1 = (0→4→1→0): v=1, out-neighbor=0, in-neighbor=4
    # C2 = (1→3→2→1): v=1, out-neighbor=3, in-neighbor=2
    # shared vertex = 1
    # From C1: v=1, a₁=0 (out), a₂=4 (in): C1 is (v→a₁→...→v) NO!
    # C1 = (0,4,1) means 0→4→1→0. Vertex 1 has incoming 4→1 and outgoing 1→0.
    # C2 = (1,3,2) means 1→3→2→1. Vertex 1 has outgoing 1→3 and incoming 2→1.
    #
    # My analysis: C₁ = v→a₁→a₂→v. But C1 has 0→4→1→0.
    # The shared vertex is 1, and in C1, vertex 1 is the "target" (a₂=1).
    # In C2, vertex 1 is the "start" (v=1).
    #
    # So: C1 = (...→4→1→0→...) cycle through {0,4,1}
    #     C2 = (1→3→2→1) cycle through {1,3,2}
    #
    # In C1, vertex 1's role: in=4→1, out=1→0
    # In C2, vertex 1's role: in=2→1, out=1→3
    # Neither edge through 1 is shared!
    #
    # This violates my earlier analysis. Let me recheck.
    # The analysis said a₁=b₂ leads to contradiction.
    # C1 written as (1→0→4→1): v=1, a₁=0, a₂=4
    # C2 written as (1→3→2→1): v=1, b₁=3, b₂=2
    # Shared vertex from {a₁,a₂}∩{b₁,b₂} = {0,4}∩{3,2} = ∅!
    #
    # My earlier analysis assumed the shared vertex is v (the first vertex).
    # But they share vertex 1, which is the "first" in one but "last" in the other!
    # When I normalize both to start at the shared vertex:
    # C1 starting at 1: (1→0→4→1), so v=1, a₁=0, a₂=4
    # C2 starting at 1: (1→3→2→1), so v=1, b₁=3, b₂=2
    # {a₁,a₂}∩{b₁,b₂} = {0,4}∩{3,2} = EMPTY!
    #
    # So they share vertex v but {a₁,a₂} and {b₁,b₂} are disjoint.
    # This means the two 3-cycles share ONLY one vertex.
    # My analysis only covered the case where {a₁,a₂}∩{b₁,b₂} ≠ ∅
    # (i.e., they share 2 vertices, hence an edge).
    #
    # CONCLUSION: Two 3-cycles CAN share exactly 1 vertex without
    # sharing any directed edge!

    print(f"\n  INSIGHT: They share EXACTLY 1 vertex (vertex 1)")
    print(f"  and the other 2 vertices of each cycle are disjoint.")
    print(f"  So n≥5 is needed for this pattern.")

    # Check if they're still homologous
    bd2, tt, el, edge_idx = build_boundary_2(A, n)
    ne = len(el)
    v1 = cycle_vec(edge_idx, (0,4,1), ne)
    v2 = cycle_vec(edge_idx, (1,3,2), ne)
    diff = v1 - v2
    hom = is_in_image(diff, bd2)
    print(f"  C1-C2 ∈ B₁? {hom}")

    if not hom:
        # Can we find a BRIDGE cycle sharing an edge with both?
        cycles3 = []
        for i in range(n):
            for j in range(n):
                if j == i or not A[i][j]: continue
                for k in range(n):
                    if k == i or k == j: continue
                    if A[j][k] and A[k][i]:
                        cycles3.append((i,j,k))

        c1_edges = {(0,4), (4,1), (1,0)}
        c2_edges = {(1,3), (3,2), (2,1)}

        for c3 in cycles3:
            c3_edges = set()
            for idx in range(3):
                c3_edges.add((c3[idx], c3[(idx+1)%3]))

            shared_with_c1 = c1_edges & c3_edges
            shared_with_c2 = c2_edges & c3_edges

            if shared_with_c1 and shared_with_c2:
                print(f"  BRIDGE: {c3} shares {shared_with_c1} with C1 "
                      f"and {shared_with_c2} with C2")

            if shared_with_c1:
                v3 = cycle_vec(edge_idx, c3, ne)
                d1 = v1 - v3
                if is_in_image(d1, bd2):
                    d2 = v3 - v2
                    if is_in_image(d2, bd2):
                        print(f"  CHAIN: C1 ~ {c3} ~ C2 via edge-sharing ✓")
                        break

    break


# =================================================================
# THE CORRECT APPROACH: Two 3-cycles sharing exactly 1 vertex
# =================================================================
print(f"\n\n{'='*72}")
print("CASE: Two 3-cycles sharing exactly 1 vertex")
print("="*72)

print("""
C₁ = (v→a→b→v) and C₂ = (v→c→d→v) where {a,b}∩{c,d} = ∅.

We need to find a CHAIN of 3-cycles connecting C₁ and C₂ where
consecutive cycles share a directed edge.

The vertices are v, a, b, c, d (5 vertices).
Consider 3-cycles on subsets of these 5 vertices.

Key observation: v→a (from C₁) and v→c (from C₂).
Consider vertices a and c: either a→c or c→a.

Case a→c: v→a→c is a path. Is (v,a,c) transitive? Need v→c: YES (from C₂).
So ∂(v,a,c) = (a,c)-(v,c)+(v,a) ∈ B₁.
This says v→a + a→c ≡ v→c mod B₁. But this gives us a transitive
triangle, not a 3-cycle bridge.

For a BRIDGE 3-cycle, we need a 3-cycle sharing edges with both C₁ and C₂.
Consider: does the 5-vertex tournament on {v,a,b,c,d} have a 3-cycle
that shares an edge with C₁ and an edge with C₂?

Let's enumerate. C₁ edges: (v,a), (a,b), (b,v).
C₂ edges: (v,c), (c,d), (d,v).

A bridge cycle must use one edge from each set.
For example: (v,a,?) using edge (v,a) from C₁ and edge (?,v) from C₂ = (d,v)?
Bridge = (v→a→?→d→v): this is a 4-cycle, not a 3-cycle.

3-cycle bridge using (v,a) from C₁ and one edge from C₂:
  - (v→a→?→v) using (?,v) which is (d,v) from C₂. Then ? = d: (v,a,d) cycle.
    Need v→a, a→d, d→v. Have v→a (from C₁), d→v (from C₂). Need a→d.
    If a→d: (v→a→d→v) is a 3-cycle sharing (v,a) with C₁ and (d,v) with C₂.

  - Using (v,c) from C₂: (v→a→c→v)? Need a→c, c→v. But c→v contradicts v→c.
    Using (c,d): (a→c→d→a)? Need a→c, c→d, d→a. From C₂: c→d. Need a→c, d→a.

So the simplest bridge: (v→a→d→v) exists iff a→d.

If a→d: BRIDGE via (v→a→d→v), sharing edge (v,a) with C₁ and (d,v) with C₂.
  KEY LEMMA: C₁ ≡ (v→a→d→v) mod B₁ [shared edge (v,a)]
  KEY LEMMA: (v→a→d→v) ≡ C₂ mod B₁ [shared edge (d,v)]
  Therefore C₁ ≡ C₂ mod B₁.

If d→a: try bridge (v→c→b→v). Need c→b, b→v. b→v from C₁. Need c→b.
  If c→b: (v→c→b→v) shares (v,c) with C₂ and (b,v) with C₁. Done.
  If b→c: try bridge (v→a→?→v) with ? ∈ {c,d}:
    (v→a→c→v): need a→c, c→v. c→v contradicts v→c. NO.
    Already tried (v→a→d→v): needs a→d, but d→a. NO.
    Try (v→c→?→v) with ? ∈ {a,b}:
    (v→c→a→v): need c→a, a→v. a→v contradicts v→a (from C₁). NO.
    (v→c→b→v): need c→b, already assumed b→c. NO.
    Try (a→?→b→v→a) cycle: 4-cycle, not 3.

  So if d→a AND b→c: NO direct bridge via (v,a)-(d,v) or (v,c)-(b,v).

  BUT: we can try OTHER edges.
  Try bridge sharing (a,b) with C₁:
    (a→b→?→a) with ? ∈ {c,d}:
    (a→b→c→a): need b→c (we have it!), c→a. From {a,c}: if c→a, YES.
    (a→b→d→a): need b→d, d→a (we have d→a). Need b→d.

  Try bridge sharing (c,d) with C₂:
    (c→d→?→c): need d→?, ?→c:
    (c→d→a→c): need d→a (YES), a→c. Need a→c.
    (c→d→b→c): need d→b, b→c (we have b→c). Need d→b.

  Case d→a, b→c:
    Subcase a→c: bridge (c→d→a→c) shares (c,d) with C₂. Then (c→d→a→c) and C₁
      need a common edge. C₁ = (v→a→b→v). Bridge edges: (c,d),(d,a),(a,c).
      C₁ edges: (v,a),(a,b),(b,v). Shared: NONE (a appears in both but no common edge).
      So: bridge (c→d→a→c) shares (c,d) with C₂ but nothing with C₁ directly.
      But (c→d→a→c) also shares vertex a with C₁.
      By KEY LEMMA (edge sharing): C₂ ≡ (c→d→a→c) [share (c,d)].
      Now connect (c→d→a→c) to C₁ = (v→a→b→v). Both use vertex a.

      From (c→d→a→c): edges through a: (d,a) incoming, (a,c) outgoing.
      From C₁: edges through a: (v,a) incoming, (a,b) outgoing.
      No shared edge through a!

      But: bridge (a→b→c→a)? Need a→b (from C₁), b→c (have it), c→a.
      If c→a: contradiction with a→c (assumed in subcase). NO.
      If a→c: bridge (a→c→d→a)? Need a→c (YES), c→d (from C₂), d→a (YES).
      (a→c→d→a) is a 3-cycle! Edges: (a,c),(c,d),(d,a).
      Shares (c,d) and (d,a) with (c→d→a→c)? Wait, (c→d→a→c) has edges (c,d),(d,a),(a,c).
      And (a→c→d→a) has edges (a,c),(c,d),(d,a). SAME edges, same cycle!

      So: connect C₁ to (a,c,d,a).
      C₁ = (v→a→b→v), (a→c→d→a). Shared vertex a but different edges.
      Need another bridge.

      Try: (v→a→c→v)? Need a→c (YES), c→v. But v→c. NO.
      Try: (d→a→b→d)? Need a→b (C₁), b→d. Need b→d.

      ...This case analysis is getting complex. Let me just check computationally
      whether the chain always exists for 5-vertex configurations.

""")

# Complete case analysis for 5-vertex tournaments
print("--- Complete check: 5-vertex cycle chain connectivity ---")

n = 5
chain_lengths = Counter()
no_chain = 0
total_pairs = 0

for A in all_tournaments(n):
    cycles3 = []
    for i in range(n):
        for j in range(n):
            if j == i or not A[i][j]: continue
            for k in range(n):
                if k == i or k == j: continue
                if A[j][k] and A[k][i]:
                    cycles3.append((i,j,k))

    if len(cycles3) < 2:
        continue

    el = edges_list(A, n)
    edge_idx = {e: i for i, e in enumerate(el)}
    bd2, tt, _, _ = build_boundary_2(A, n)
    ne = len(el)

    def get_edges(c):
        return {(c[i], c[(i+1)%3]) for i in range(3)}

    # Build edge-sharing adjacency graph on 3-cycles
    nc = len(cycles3)
    adj = defaultdict(set)
    for i in range(nc):
        for j in range(i+1, nc):
            if get_edges(cycles3[i]) & get_edges(cycles3[j]):
                adj[i].add(j)
                adj[j].add(i)

    # BFS for shortest path between all pairs
    for i in range(nc):
        for j in range(i+1, nc):
            total_pairs += 1
            # BFS from i to j
            visited = {i}
            queue = [(i, 0)]
            found = False
            while queue:
                node, dist = queue.pop(0)
                if node == j:
                    chain_lengths[dist] += 1
                    found = True
                    break
                for nbr in adj[node]:
                    if nbr not in visited:
                        visited.add(nbr)
                        queue.append((nbr, dist+1))

            if not found:
                # Not edge-sharing connected. Are they still homologous?
                v1 = cycle_vec(edge_idx, cycles3[i], ne)
                v2 = cycle_vec(edge_idx, cycles3[j], ne)
                if is_in_image(v1 - v2, bd2):
                    chain_lengths['inf_but_homologous'] += 1
                else:
                    no_chain += 1

print(f"\nn=5: {total_pairs} total 3-cycle pairs")
print(f"  Chain length distribution: {dict(chain_lengths)}")
print(f"  Not homologous: {no_chain}")


# =================================================================
# THE RIGHT PROOF STRATEGY
# =================================================================
print(f"\n\n{'='*72}")
print("CORRECT PROOF APPROACH")
print("="*72)

print("""
The computation shows that edge-sharing chains don't always connect
all 3-cycles. But all 3-cycles are STILL homologous. This means the
boundary space B₁ creates additional relations beyond the edge-sharing
lemma.

BETTER APPROACH: Don't try to chain 3-cycles. Instead, show directly
that dim(H₁) ≤ 1 by counting dimensions.

OBSERVATION from the data:
  β₁ = C(n-1,2) - rank(∂₂)
  rank(∂₂) ∈ {C(n-1,2)-1, C(n-1,2)} for all tournaments.

So: rank(∂₂) always takes one of just TWO values! This is remarkable.

ALTERNATIVE PROOF APPROACH:
Show that for ANY tournament T with t₃ ≥ 1:
  rank(∂₂) ≥ C(n-1,2) - 1.

This is equivalent to: the image of ∂₂ has codimension ≤ 1 in Z₁.

Since ∂₂ maps FROM Ω₂ (which has dimension C(n,3) - t₃) into Z₁
(which has dimension C(n-1,2)):

  rank(∂₂) ≤ min(dim(Ω₂), dim(Z₁)) = min(C(n,3)-t₃, C(n-1,2))

For n ≥ 4: C(n,3)-t₃ ≥ C(n,3)-C(n,3) = 0 but typically ≫ C(n-1,2).
In fact, C(n,3) = n(n-1)(n-2)/6 while C(n-1,2) = (n-1)(n-2)/2.
So dim(Ω₂) ≥ n/3 * dim(Z₁) for small t₃. Plenty of room.

Even for maximal t₃ (regular tournament at n=5: t₃=5, dim(Ω₂)=5,
dim(Z₁)=6): rank(∂₂)=5, so codimension = 1. The 5 transitive triples
are "almost enough" to span Z₁.
""")

# What IS the missing direction?
print("--- What direction does B₁ miss? ---")

for n in [5]:
    for A in all_tournaments(n):
        bd2, tt, el, edge_idx = build_boundary_2(A, n)
        ne = len(el)
        bd1_mat = np.zeros((n, ne))
        for j, (a,b) in enumerate(el):
            bd1_mat[b, j] += 1
            bd1_mat[a, j] -= 1

        U1, S1, Vt1 = np.linalg.svd(bd1_mat, full_matrices=True)
        rank1 = sum(s > 1e-10 for s in S1)
        Z1_basis = Vt1[rank1:].T

        # B_1 in Z_1 coordinates
        B1_Z1 = Z1_basis.T @ bd2

        if B1_Z1.shape[1] == 0:
            continue
        r = np.linalg.matrix_rank(B1_Z1, tol=1e-10)
        dim_Z1 = Z1_basis.shape[1]

        if r < dim_Z1:
            # Find the missing direction
            U, S, Vt = np.linalg.svd(B1_Z1, full_matrices=True)
            missing = U[:, r:]  # columns in Z_1 coords not hit by B_1
            # Convert to edge space
            missing_edge = Z1_basis @ missing

            for k in range(missing.shape[1]):
                vec = missing_edge[:, k]
                vec[np.abs(vec) < 1e-10] = 0
                nz = vec[vec != 0]
                if len(nz) > 0:
                    scale = min(abs(nz))
                    vec = vec / scale
                    vec = np.round(vec).astype(int)

                terms = []
                for i, c in enumerate(vec):
                    if c != 0:
                        if c == 1:
                            terms.append(f"({el[i][0]},{el[i][1]})")
                        elif c == -1:
                            terms.append(f"-({el[i][0]},{el[i][1]})")
                        else:
                            terms.append(f"{c}({el[i][0]},{el[i][1]})")
                # Check if proportional to sum of all 3-cycles
                print(f"  Missing direction: {' + '.join(terms)}")

            # Check: is the missing direction proportional to
            # the "total 3-cycle" = sum of all oriented 3-cycles?
            total_3cycle = np.zeros(ne)
            seen = set()
            for i in range(n):
                for j in range(n):
                    if j == i or not A[i][j]: continue
                    for kk in range(n):
                        if kk == i or kk == j: continue
                        if A[j][kk] and A[kk][i]:
                            triple = frozenset([i,j,kk])
                            if triple not in seen:
                                seen.add(triple)
                                total_3cycle[edge_idx[(i,j)]] += 1
                                total_3cycle[edge_idx[(j,kk)]] += 1
                                total_3cycle[edge_idx[(kk,i)]] += 1

            # Project into Z_1
            total_Z1 = Z1_basis.T @ total_3cycle
            missing_Z1 = missing[:, 0]

            # Check proportionality
            ratio = None
            prop = True
            for i in range(len(total_Z1)):
                if abs(missing_Z1[i]) > 1e-10:
                    r = total_Z1[i] / missing_Z1[i]
                    if ratio is None:
                        ratio = r
                    elif abs(r - ratio) > 1e-8:
                        prop = False
                        break
                elif abs(total_Z1[i]) > 1e-10:
                    prop = False
                    break

            if prop and ratio is not None:
                print(f"  → Proportional to total 3-cycle vector (ratio={ratio:.2f})")
            else:
                print(f"  → NOT proportional to total 3-cycle vector")

            break  # first β₁=1 tournament only


# =================================================================
# KEY COMPUTATION: Universal relation in Z₁
# =================================================================
print(f"\n\n{'='*72}")
print("UNIVERSAL Z₁ RELATION")
print("="*72)

print("""
Define φ ∈ Z₁* (a linear functional on Z₁) by:
  φ(e) = 1 for every directed edge e in T.

Wait — that's not well-defined on Z₁ since Z₁ elements are formal
sums of edges. What I mean is:

Define the vector w ∈ R^{edges} by w_e = 1 for all e.
Then w is NOT in Z₁ (it has nonzero boundary: ∂w = sum of in-degree
minus out-degree at each vertex, which is n-1 - out_deg_v at vertex v).

Actually, in a TOURNAMENT: each vertex has in-deg + out-deg = n-1.
So w ∈ Z₁ iff at every vertex, in-deg = out-deg, i.e., tournament
is regular. For odd n, regular tournaments exist and w ∈ Z₁.

For general tournaments, w ∉ Z₁. So this doesn't directly help.

Let me try a different approach: the "score functional".
Define s_v = out-degree of v. Then ∂_1(sum of edges) at vertex v
gives in-degree(v) - out-degree(v) = (n-1-s_v) - s_v = n-1-2s_v.

A 1-chain z is in Z₁ iff at each vertex v, the sum of coefficients
of incoming edges minus outgoing edges is 0.
""")

# Let me investigate what EXACTLY the missing direction is for each
# tournament type. Is it always the same thing?

print("--- Missing direction for each score sequence at n=5 ---")

n = 5
missing_by_score = defaultdict(list)

for A in all_tournaments(n):
    bd2, tt, el, edge_idx = build_boundary_2(A, n)
    ne = len(el)

    bd1_mat = np.zeros((n, ne))
    for j, (a,b) in enumerate(el):
        bd1_mat[b, j] += 1
        bd1_mat[a, j] -= 1

    U1, S1, Vt1 = np.linalg.svd(bd1_mat, full_matrices=True)
    rank1 = sum(s > 1e-10 for s in S1)
    Z1_basis = Vt1[rank1:].T

    B1_Z1 = Z1_basis.T @ bd2
    if B1_Z1.shape[1] == 0:
        continue
    r = np.linalg.matrix_rank(B1_Z1, tol=1e-10)
    dim_Z1 = Z1_basis.shape[1]

    if r < dim_Z1:
        score = tuple(sorted(sum(A[i]) for i in range(n)))
        t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                 if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]))

        # The missing direction in edge space
        U, S, Vt = np.linalg.svd(B1_Z1, full_matrices=True)
        missing_Z1 = U[:, r:]
        missing_edge = Z1_basis @ missing_Z1

        vec = missing_edge[:, 0].copy()
        vec[np.abs(vec) < 1e-10] = 0
        nz = vec[vec != 0]
        if len(nz) > 0:
            scale = min(abs(nz))
            vec = vec / scale
            vec = np.round(vec, 6)

        # What are the coefficients grouped by vertex degree?
        # For each edge (u,v), what is score[u] and score[v]?
        coeff_by_scores = defaultdict(list)
        for i in range(ne):
            if abs(vec[i]) > 1e-10:
                u, v = el[i]
                su = sum(A[u])
                sv = sum(A[v])
                coeff_by_scores[(su, sv)].append(vec[i])

        missing_by_score[(score, t3)].append(dict(coeff_by_scores))

# Print summary
for (score, t3) in sorted(missing_by_score.keys()):
    examples = missing_by_score[(score, t3)]
    print(f"\n  score={score}, t₃={t3}: {len(examples)} tournaments")
    # Show first example
    ex = examples[0]
    for (su, sv), coeffs in sorted(ex.items()):
        print(f"    edges with out-deg ({su},{sv}): coefficients {coeffs[:3]}...")


# =================================================================
# THE VERTEX FLOW APPROACH
# =================================================================
print(f"\n\n{'='*72}")
print("THE FLOW INTERPRETATION")
print("="*72)

print("""
A 1-cycle z assigns a flow value z_e to each edge e, with flow
conservation at each vertex. The boundary B₁ consists of flows
that can be decomposed into transitive triangle circulations.

∂(a,b,c) = (a,b) + (b,c) - (a,c) sends 1 unit around the path
a→b→c and -1 unit on the shortcut a→c. Net effect: a "triangular
circulation" that transfers flow from the diagonal to the two-step
path.

KEY INSIGHT: B₁ consists of all flows expressible as sums of such
circulations. The quotient H₁ = Z₁/B₁ measures flows that CANNOT
be decomposed into transitive circulations.

A directed 3-cycle (a→b→c→a) sends 1 unit around all three edges.
Since there's no transitive triple in {a,b,c}, this circulation
can't be decomposed into transitive triangles. So it represents
a nontrivial element of H₁ (when β₁=1).

But CAN two independent 3-cycle flows be simultaneously indecomposable?
The answer is NO: any two 3-cycles generate the same 1-dimensional
quotient. This is because:

Consider the "score-weighted Euler class" (will verify):
Define ψ: Z₁ → R by ψ(z) = Σ_e z_e * weight(e)
where weight(a→b) = some function of scores of a and b.

If ψ vanishes on all transitive boundaries and is nonzero on 3-cycles,
then ker(ψ) ⊃ B₁ and dim(Z₁/ker(ψ)) = 1, so dim(H₁) ≤ 1.

Let's check: ∂(a,b,c) = (a,b)+(b,c)-(a,c) with a→b→c, a→c.
ψ(∂(a,b,c)) = w(a,b)+w(b,c)-w(a,c).

For ψ to vanish on all boundaries: w(a,b)+w(b,c) = w(a,c) whenever
a→b→c and a→c (transitive triple).

This is a COCYCLE condition! w: edges → R is a 1-cochain, and
ψ vanishes on B₁ iff w is a 1-cocycle (δw = 0 on Ω₂ = transitive triples).

So: dim(H₁) = dim of space of 1-cocycles vanishing on Ω₂ that are
nonzero on Z₁. This is the first cohomology H¹!

By GLMY duality: H₁ ≅ H¹ (path homology = path cohomology).
So β₁ = dim(H¹) too.

Can we construct the universal cocycle?

For tournaments: try w(a→b) = s_a - s_b (score difference).
Then w(a,b)+w(b,c)-w(a,c) = (s_a-s_b)+(s_b-s_c)-(s_a-s_c) = 0.
So score difference IS a cocycle! But it vanishes on ALL of Z₁
(since it's exact: w = δf where f(v) = s_v).

Try w(a→b) = s_a (just the source score).
w(a,b)+w(b,c)-w(a,c) = s_a+s_b-s_a = s_b. This is NOT zero in general.
So this is NOT a cocycle.

Try w(a→b) = 1 (constant).
w(a,b)+w(b,c)-w(a,c) = 1+1-1 = 1 ≠ 0. NOT a cocycle.

Interesting — so the "constant weight" cocycle doesn't work. But
it gives ψ(3-cycle) = 3 ≠ 0. And ψ(∂) = 1 for every boundary.
So ψ: Z₁ → R takes value 1 on each boundary, value 3 on each 3-cycle,
value k on each k-cycle.

Actually: ψ(z) = Σ z_e (sum of all coefficients) is NOT zero on
boundaries, so it doesn't factor through H₁.

The CORRECT approach: find a weight function that IS zero on boundaries.
""")

# Let me just directly compute H¹ and find the cocycle
print("\n--- Computing the cocycle for β₁=1 tournaments ---")

for n in [5]:
    for A in all_tournaments(n):
        t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                 if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]))
        if t3 != 5:  # regular tournament
            continue

        el = edges_list(A, n)
        edge_idx = {e: i for i, e in enumerate(el)}
        ne = len(el)
        bd2, tt, _, _ = build_boundary_2(A, n)

        # Find cocycle: w ∈ R^{edges} such that for all transitive (a,b,c):
        #   w(a,b) + w(b,c) - w(a,c) = 0
        # This means w^T ∂_2 = 0, i.e., w ∈ ker(∂_2^T) = left null space of ∂_2.

        # Equivalently: w is in the left null space of bd2
        U, S, Vt = np.linalg.svd(bd2, full_matrices=True)
        rank2 = sum(s > 1e-10 for s in S)
        # Left null space: rows of U beyond rank
        left_null = U[:, rank2:]  # columns are left null vectors

        print(f"\nRegular tournament n=5 (t₃=5):")
        print(f"  dim(left null of ∂₂) = {left_null.shape[1]}")
        print(f"  dim(edge space) = {ne}")
        print(f"  rank(∂₂) = {rank2}")

        # Each left null vector is a cocycle (zero on all boundaries)
        # How many are also zero on Z₁? (These are coboundaries.)
        bd1_mat = np.zeros((n, ne))
        for j, (a,b) in enumerate(el):
            bd1_mat[b, j] += 1
            bd1_mat[a, j] -= 1

        # Coboundaries = image of δ_0: R^vertices → R^edges
        # δ_0(f)(a→b) = f(b) - f(a)
        delta0 = bd1_mat.T  # ne × n matrix

        # H¹ = cocycles / coboundaries
        # cocycles = left_null columns
        # Check: how many cocycles are NOT coboundaries?
        cocycles = left_null
        # Project coboundaries into cocycle space
        cobdy_in_cocycle = cocycles.T @ delta0  # cocycle_dim × n

        rank_cobdy = np.linalg.matrix_rank(cobdy_in_cocycle, tol=1e-10)
        dim_H1 = cocycles.shape[1] - rank_cobdy
        print(f"  dim(cocycles) = {cocycles.shape[1]}")
        print(f"  dim(coboundaries in cocycle space) = {rank_cobdy}")
        print(f"  dim(H¹) = {dim_H1}")

        # Find the non-coboundary cocycle
        if dim_H1 > 0:
            U_c, S_c, Vt_c = np.linalg.svd(cobdy_in_cocycle, full_matrices=True)
            noncobdy = U_c[:, rank_cobdy:]  # in cocycle coordinates
            # Convert to edge space
            cocycle_edge = cocycles @ noncobdy

            for k in range(cocycle_edge.shape[1]):
                w = cocycle_edge[:, k]
                w[np.abs(w) < 1e-10] = 0
                nz = w[w != 0]
                if len(nz) > 0:
                    scale = min(abs(nz))
                    w = w / scale
                    w = np.round(w).astype(int)

                print(f"\n  H¹ cocycle (weight function):")
                for i in range(ne):
                    if w[i] != 0:
                        print(f"    w({el[i][0]}→{el[i][1]}) = {w[i]}")

                # Evaluate on a 3-cycle
                for c in [(0,1,2),(0,2,3),(0,3,4),(0,4,1)]:
                    try:
                        val = sum(w[edge_idx[(c[j], c[(j+1)%3])]] for j in range(3))
                        print(f"  w(3-cycle {c}) = {val}")
                    except KeyError:
                        pass

        break

print(f"\n\nDone.")
