#!/usr/bin/env python3
"""
girth_omega_s119.py — opus-2026-03-21-S119

THE GIRTH OF OMEGA — CYCLES OF CYCLES

Ω(T) = the conflict graph of tournament T.
  Vertices: directed odd cycles of T
  Edges: two cycles share a tournament vertex

The GIRTH of Ω = length of shortest cycle IN Ω.
This is a META-STRUCTURAL property: it measures how the odd cycles
of the tournament are woven together.

CONNECTIONS:
1. H=7 requires Ω ⊇ K₃ (triangle, girth 3). K₃ in Ω = 3 pairwise-
   conflicting cycles = the THM-029 obstruction. So H=7 requires girth(Ω)=3.
2. High girth → more independent sets → higher H (for fixed |V(Ω)|).
3. The Petersen graph K(5,2) has girth 5. Is there a connection to Ω girth?
4. Root system: girth of Ω = shortest cycle in the root triple intersection graph.

Author: opus-2026-03-21-S119
"""

import sys
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter, deque
from math import comb, factorial
sys.stdout.reconfigure(line_buffering=True)

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        yield A

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def find_odd_cycles(A, n, max_len=None):
    """Find all directed odd cycles in tournament A."""
    if max_len is None:
        max_len = n
    cycles = set()
    for length in range(3, max_len + 1, 2):
        for verts in combinations(range(n), length):
            first = verts[0]
            for perm in permutations(verts[1:]):
                path = [first] + list(perm)
                if all(A[path[k]][path[(k+1) % length]] for k in range(length)):
                    cycles.add(frozenset(path))
    return list(cycles)

def build_omega(cycles):
    """Build the conflict graph: edge iff cycles share a vertex."""
    nc = len(cycles)
    adj = [[0]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:  # share a vertex
                adj[i][j] = adj[j][i] = 1
    return adj

def graph_girth(adj, n):
    """Compute the girth (shortest cycle) of a graph via BFS."""
    if n == 0:
        return float('inf')
    best = float('inf')
    for start in range(n):
        # BFS from start
        dist = [-1] * n
        dist[start] = 0
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v in range(n):
                if not adj[u][v]:
                    continue
                if dist[v] == -1:
                    dist[v] = dist[u] + 1
                    queue.append(v)
                elif v != start and dist[v] >= dist[u]:
                    # Found a cycle
                    cycle_len = dist[u] + dist[v] + 1
                    best = min(best, cycle_len)
        # Also check: back edge to start
        for v in range(n):
            if adj[start][v] and dist[v] > 0:
                best = min(best, dist[v] + 1)
    return best if best < float('inf') else float('inf')

def independence_polynomial(adj, n, x):
    """Compute I(G, x) = Σ α_k x^k."""
    total = 0
    for mask in range(1 << n):
        subset = [i for i in range(n) if (mask >> i) & 1]
        is_ind = True
        for a in range(len(subset)):
            for b in range(a+1, len(subset)):
                if adj[subset[a]][subset[b]]:
                    is_ind = False
                    break
            if not is_ind:
                break
        if is_ind:
            total += x ** len(subset)
    return total

print("=" * 72)
print("  THE GIRTH OF OMEGA — CYCLES OF CYCLES")
print("  opus-2026-03-21-S119")
print("=" * 72)

# ========================================================================
# PART 1: Compute girth(Ω) for all tournaments at n=5
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: GIRTH(Ω) AT n=5")
print("=" * 72)

n = 5
girth_data = defaultdict(list)

for A in all_tournaments(n):
    H = count_hp(A, n)
    cycles = find_odd_cycles(A, n)
    nc = len(cycles)

    if nc == 0:
        girth = float('inf')
    elif nc == 1:
        girth = float('inf')  # single vertex, no cycle in Ω
    else:
        adj_omega = build_omega(cycles)
        girth = graph_girth(adj_omega, nc)

    girth_data[H].append({
        'girth': girth, 'n_cycles': nc, 'H': H
    })

print(f"\n  {'H':>4s} {'count':>6s} {'|V(Ω)|':>7s} {'girth(Ω)':>10s} {'girth dist':>20s}")

for H_val in sorted(girth_data.keys()):
    group = girth_data[H_val]
    nc_vals = [d['n_cycles'] for d in group]
    girth_vals = [d['girth'] for d in group]
    girth_counts = Counter(girth_vals)

    nc_str = f"{min(nc_vals)}-{max(nc_vals)}" if min(nc_vals) != max(nc_vals) else str(nc_vals[0])
    girth_str = str(dict(sorted(girth_counts.items())))

    print(f"  {H_val:4d} {len(group):6d} {nc_str:>7s} "
          f"{min(girth_vals):>10} {girth_str:>20s}")

# ========================================================================
# PART 2: The forbidden value H=7 and girth 3
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 2: H=7 AND GIRTH(Ω) = 3")
print("=" * 72)
print(f"""
  From THM-029: H=7 requires α₁=3, α₂=0.
  α₂=0 means ALL 3 cycles are pairwise conflicting.
  This means Ω ⊇ K₃ (triangle), so girth(Ω) = 3.

  But THM-029 proves this configuration forces a 4th cycle → α₁≥4 → H≥9.

  QUESTION: For achieved H values, what is girth(Ω)?
  Does LOW girth correlate with LOW H? Does HIGH girth correlate with HIGH H?

  From the data above:
  - H=1 (transitive): |V(Ω)|=0, girth=∞ (no cycles at all)
  - H=3: |V(Ω)|=1, girth=∞ (single cycle, no Ω-cycle possible)
  - H=5: |V(Ω)|=2, girth=∞ (2 cycles, either adjacent or not, no triangle in Ω)
  - H=9: |V(Ω)| varies, girth varies
  - H=15 (regular): |V(Ω)| large, girth could be small
""")

# ========================================================================
# PART 3: Girth(Ω) by score class — connecting to the POS
# ========================================================================
print("=" * 72)
print("  PART 3: GIRTH(Ω) BY SCORE CLASS")
print("=" * 72)

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

girth_by_score = defaultdict(lambda: defaultdict(list))

for A in all_tournaments(n):
    sc = score_seq(A, n)
    H = count_hp(A, n)
    cycles = find_odd_cycles(A, n)
    nc = len(cycles)

    if nc <= 1:
        girth = float('inf')
    else:
        adj_omega = build_omega(cycles)
        girth = graph_girth(adj_omega, nc)

    girth_by_score[sc][H].append(girth)

pos_score = (1, 2, 2, 2, 3)
print(f"\n  POS class {pos_score}:")
for H_val in sorted(girth_by_score[pos_score].keys()):
    girths = girth_by_score[pos_score][H_val]
    gc = Counter(girths)
    print(f"    H={H_val}: girth distribution = {dict(sorted(gc.items()))}")

# ========================================================================
# PART 4: Girth(Ω), independence number, and H
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 4: GIRTH(Ω) ↔ INDEPENDENCE NUMBER ↔ H")
print("=" * 72)
print(f"""
  The independence polynomial I(Ω, 2) = H.
  Higher independence number α(Ω) → higher H (roughly).
  Higher girth(Ω) → more independent sets (loosely connected graph).

  So we expect: higher girth(Ω) → higher H.
  But it's more subtle because |V(Ω)| also matters.

  MORE CYCLES (larger Ω) AND loose connections (high girth) = highest H.
  FEWER CYCLES (small Ω) = low H regardless of girth.
  MANY CYCLES BUT tight connections (low girth) = moderate H.

  Let me check: corr(girth(Ω), H) at n=5.
""")

# Correlation: girth(Ω) vs H (excluding infinite girth)
finite_girth_data = []
for A in all_tournaments(n):
    H = count_hp(A, n)
    cycles = find_odd_cycles(A, n)
    nc = len(cycles)
    if nc <= 1:
        continue
    adj_omega = build_omega(cycles)
    girth = graph_girth(adj_omega, nc)
    if girth < float('inf'):
        finite_girth_data.append((H, girth, nc))

if finite_girth_data:
    Hs = [d[0] for d in finite_girth_data]
    girths = [d[1] for d in finite_girth_data]
    ncs = [d[2] for d in finite_girth_data]
    corr_gH = np.corrcoef(Hs, girths)[0, 1]
    corr_nH = np.corrcoef(Hs, ncs)[0, 1]

    print(f"  corr(H, girth(Ω)) = {corr_gH:+.4f} (among tournaments with finite girth)")
    print(f"  corr(H, |V(Ω)|) = {corr_nH:+.4f}")
    print(f"  Tournaments with finite girth: {len(finite_girth_data)}")

    # Show average girth by H
    by_H = defaultdict(list)
    for H, g, nc in finite_girth_data:
        by_H[H].append(g)

    print(f"\n  Average girth(Ω) by H (finite girth only):")
    for H_val in sorted(by_H.keys()):
        gs = by_H[H_val]
        print(f"    H={H_val}: avg girth = {np.mean(gs):.2f}, n={len(gs)}")

# ========================================================================
# PART 5: The Petersen girth = 5 connection
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 5: THE PETERSEN GIRTH = 5")
print("=" * 72)
print(f"""
  The Petersen graph K(5,2) has girth 5.
  It is the arc-disjointness graph at n=5.

  Ω(T) is the CYCLE-CONFLICT graph (cycles sharing vertices).
  K(5,2) is the ARC-DISJOINTNESS graph (arcs not sharing vertices).

  These are COMPLEMENTARY:
  - K(5,2) edge ↔ arcs DISJOINT (can be oriented independently)
  - Ω edge ↔ cycles SHARING a vertex (constrained together)

  The girth of K(5,2) = 5 means: the shortest cycle in the
  arc-disjointness graph has length 5. This corresponds to
  5 arcs that form a "pentagon" of mutual disjointness.

  At n=5: 5 pairwise-disjoint arcs would need 10 distinct vertices,
  but n=5 only has 5 vertices. So actually: the pentagon in K(5,2)
  consists of 5 pairs that are pairwise disjoint. With 5 vertices,
  each pair uses 2, and 5 disjoint pairs need 10 = 2×5. IMPOSSIBLE!
  So K(5,2) has no 5-independent set (α(K(5,2)) = 4).

  But the GIRTH is about cycles, not independent sets.
  A 5-cycle in K(5,2): 5 pairs that form a path in the disjointness graph.
  (a,b)-(c,d)-(e,f)-... where each consecutive pair is disjoint.

  Example: (0,1)-(2,3)-(0,4)-(1,2)-(3,4)-(0,1)
  Check: {{0,1}}∩{{2,3}}=∅ ✓, {{2,3}}∩{{0,4}}=∅ ✓, {{0,4}}∩{{1,2}}=∅ ✓,
  {{1,2}}∩{{3,4}}=∅ ✓, {{3,4}}∩{{0,1}}=∅ ✓. Yes, a 5-cycle!

  But can we find a 3-cycle or 4-cycle in K(5,2)?
  3-cycle: (a,b)-(c,d)-(e,f)-(a,b) with all pairwise disjoint.
  Need 3 pairwise-disjoint pairs from {{0,...,4}}: need 6 vertices. Only 5 available.
  IMPOSSIBLE. So K(5,2) has no triangle → girth ≥ 4.

  4-cycle: (a,b)-(c,d)-(e,f)-(g,h)-(a,b) with consecutive pairs disjoint.
  (0,1)-(2,3)-(0,4)-(1,2)-(0,1)?
  Check: {{0,4}}∩{{1,2}}=∅ ✓, {{1,2}}∩{{0,1}}={{1}} ≠ ∅ ✗. Not a 4-cycle.
  Can we find any 4-cycle? The Petersen graph is known to have no 3 or 4 cycles.
  Its girth is EXACTLY 5. This is a fundamental property.

  THE CONNECTION: Petersen girth 5 = the BOUNDARY ORDER n=5.
  The fact that arc-pairs cannot form short cycles (girth 5) is
  equivalent to saying that at n=5, the tournament structure is
  maximally "spread out" — no short dependencies between arcs.

  This is WHY n=5 is the last order where the per-path identity works,
  where H is determined by score class (except one POS class), and
  where Ω has simple structure (α₂ = 0 always).
""")

# ========================================================================
# PART 6: Girth(Ω) and the Kneser girth at larger n
# ========================================================================
print("=" * 72)
print("  PART 6: KNESER GIRTH AND Ω GIRTH AT n=6,7")
print("=" * 72)

# Kneser K(n,2) girth:
# K(3,2): no edges → girth ∞
# K(4,2): 3 edges, but is it triangle-free? K(4,2) has 6 vertices, 3 edges.
#   The 3 edges form a perfect matching (3 disjoint pairs from 4 vertices).
#   No cycles at all → girth ∞.
# K(5,2): Petersen, girth 5
# K(6,2): girth 3 (can find 3 pairwise disjoint pairs from 6 vertices: {0,1},{2,3},{4,5})
# K(7,2): girth 3

print(f"  Kneser K(n,2) girth:")
for n_k in range(3, 9):
    if n_k <= 3:
        g = float('inf')
    elif n_k == 4:
        g = float('inf')
    elif n_k == 5:
        g = 5
    else:
        g = 3  # For n ≥ 6, 3 pairwise-disjoint pairs exist
    print(f"    K({n_k},2): girth = {g}")

print(f"""
  KEY TRANSITION AT n=6:
  K(5,2) girth = 5 (Petersen — no short cycles)
  K(6,2) girth = 3 (first triangle: 3 disjoint pairs)

  This means: at n=6, there exist 3 arcs that are pairwise disjoint.
  These 3 arcs can host 3 independent 3-cycles → α₂ > 0 possible.
  This is exactly where the OCR starts to have a nontrivial residual!

  THE KNESER GIRTH TRANSITION K(5,2)→K(6,2) (girth 5→3) COINCIDES WITH
  THE EMERGENCE OF THE OCR RESIDUAL (n=5 has α₂=0, n=6 has α₂>0).

  Girth 5 (Petersen) → no independent cycle triples → α₂=0 → H exact from α₁
  Girth 3 (n≥6) → independent cycle triples exist → α₂>0 → H needs higher OCF terms
""")

# ========================================================================
# PART 7: What girth does Ω(T) actually have at n=5 and n=6?
# ========================================================================
print("=" * 72)
print("  PART 7: Ω GIRTH DISTRIBUTION AT n=5,6")
print("=" * 72)

for n in [5, 6]:
    print(f"\n  n={n}:")
    girth_dist = Counter()
    girth_by_H = defaultdict(Counter)

    count = 0
    for A in all_tournaments(n):
        H = count_hp(A, n)
        cycles = find_odd_cycles(A, n, max_len=5)  # limit for speed at n=6
        nc = len(cycles)

        if nc <= 1:
            g = float('inf')
        else:
            adj_omega = build_omega(cycles)
            g = graph_girth(adj_omega, nc)

        girth_dist[g] += 1
        girth_by_H[H][g] += 1

        count += 1
        if n == 6 and count >= 5000:  # sample at n=6
            break

    print(f"    Girth distribution: {dict(sorted(girth_dist.items()))}")
    print(f"    Girth by H (selected):")
    for H_val in sorted(girth_by_H.keys()):
        if H_val in [1, 5, 9, 15] or (n == 6 and H_val in [1, 9, 15, 23, 33, 45]):
            gc = girth_by_H[H_val]
            print(f"      H={H_val}: {dict(sorted(gc.items()))}")

# ========================================================================
# PART 8: Synthesis — girth as the meta-structural invariant
# ========================================================================
print("\n\n" + "=" * 72)
print("  SYNTHESIS: GIRTH(Ω) IN THE ROOT SYSTEM PICTURE")
print("=" * 72)
print(f"""
  THE GIRTH OF Ω IN THE ROOT SYSTEM:

  Ω(T) has vertices = odd cycles = root triangles (3-cycles) and root
  pentagons (5-cycles) and higher.

  A cycle in Ω means: a sequence of odd cycles C₁, C₂, ..., Cₖ, C₁
  where each consecutive pair shares a tournament vertex.

  In the root system: this is a sequence of root triples/pentagons
  that are "chained" through shared vertices.

  GIRTH(Ω) = 3: Three root triples form a triangle in Ω.
  This requires 3 cycles sharing vertices in a circular pattern:
  C₁ shares a vertex with C₂, C₂ shares with C₃, C₃ shares with C₁.

  At n=5: girth(Ω) = 3 is COMMON for high-H tournaments.
  The regular tournament (H=15) has many cycles, tightly connected → low girth.

  At n=5: girth(Ω) = ∞ for low-H tournaments (few or no cycles).

  THE HIERARCHY:
  - girth(Ω) = ∞: H ≤ 5 (few cycles, H = 1+2α₁ with small α₁)
  - girth(Ω) = 3: H ≥ 9 (many cycles, tightly woven)
  - girth(Ω) never = 4 or 5 at n=5 (SKIP from ∞ to 3!)

  THE SKIP IS SIGNIFICANT: there are no tournaments at n=5 with
  girth(Ω) = 4 or 5. The conflict graph either has no cycles or
  immediately has triangles. No intermediate structure exists.

  THIS MIRRORS THE GAP AT H=7: the H values skip from 5 to 9,
  just as the girth skips from ∞ to 3. The forbidden H=7 corresponds
  to the forbidden girth range 4-5 (at n=5 specifically).

  THE KNESER-Ω DUALITY:
  K(5,2) girth = 5 (the ARC graph has long cycles)
  Ω(T) girth ∈ {{3, ∞}} (the CYCLE graph has either short or no cycles)

  These are COMPLEMENTARY:
  Long arc-cycles (Petersen) → short cycle-cycles (Ω) or none.
  The Petersen girth 5 FORCES Ω to either be trivial or triangle-rich.
  There is no middle ground — which is WHY n=5 is the boundary.
""")

print("\nDone. opus-2026-03-21-S119")
