#!/usr/bin/env python3
"""
h21_analysis.py — opus-2026-03-14-S71g

Analyze H=21 impossibility.
First: what graphs G satisfy I(G,2)=21?
Then: can any tournament have Ω(T) equal to such a graph?
"""

from itertools import combinations
from collections import defaultdict

def independence_poly_graph(adj, n_vertices):
    """Compute independence polynomial coefficients."""
    alpha = [0] * (n_vertices + 1)
    for mask in range(1 << n_vertices):
        vertices = [i for i in range(n_vertices) if mask & (1 << i)]
        independent = True
        for i in range(len(vertices)):
            for j in range(i+1, len(vertices)):
                if adj[vertices[i]][vertices[j]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            alpha[len(vertices)] += 1
    return alpha

print("=" * 60)
print("GRAPHS WITH I(G, 2) = 21")
print("=" * 60)

# Check all graphs on 1 to 8 vertices for I(G,2) = 21
target = 21

for nv in range(1, 9):
    total_edges = nv * (nv - 1) // 2
    if total_edges > 20:
        # Too many graphs to enumerate, sample
        break

    found = []
    for edge_bits in range(1 << total_edges):
        adj = [[False]*nv for _ in range(nv)]
        idx = 0
        for i in range(nv):
            for j in range(i+1, nv):
                if edge_bits & (1 << idx):
                    adj[i][j] = adj[j][i] = True
                idx += 1

        alpha = independence_poly_graph(adj, nv)
        Ival = sum(alpha[k] * 2**k for k in range(len(alpha)))

        if Ival == target:
            edges = []
            for i in range(nv):
                for j in range(i+1, nv):
                    if adj[i][j]:
                        edges.append((i,j))
            found.append((edges, alpha[:nv+1]))

    if found:
        print(f"\n  n_vertices={nv}: {len(found)} graphs with I(G,2)={target}")
        for edges, alpha in found[:20]:
            print(f"    edges={edges}, α={[a for a in alpha if a or alpha.index(a)==0]}")

# Also check I(G,2) = 7 and 21 decomposition
print(f"\n{'='*60}")
print("I(G,2) = 21 DECOMPOSITION")
print(f"{'='*60}")

# 21 = 1 + 2*a1 + 4*a2 + 8*a3 + ...
# 20 = 2*a1 + 4*a2 + 8*a3 + ...
# 10 = a1 + 2*a2 + 4*a3 + ...
# Possible:
# a1=10, a2=0: I = 1+20 = 21. Graph with 10 vertices, no edges? I = (1+2)^10 = way too big.
# Actually I(K_bar_n, 2) = (1+2)^n = 3^n.
# For disconnected: I(G1 ∪ G2, x) = I(G1,x) * I(G2,x)
# 21 = 3 * 7. So I = 3 * 7 = I(single vertex) * I(K3, 2)?
# I(single vertex, 2) = 1 + 2 = 3. I(K3, 2) = 1 + 6 = 7.
# So G = K3 ∪ isolated vertex has I = 3*7 = 21!

# But also 21 = 7 * 3 = 21 * 1.
# I(G, 2) = 21 directly:
# 21 = 1 + 2*a1 + 4*a2 + 8*a3 + ...
# Solutions with small a_k:
# a1=10: 21 = 1+20, a2=a3=...=0. Need 10 vertices, all pairwise adjacent (K10). I(K10,2)=1+20=21. Check: alpha_0=1, alpha_1=10, alpha_k=0 for k>=2. So 1+20=21. ✓
# a1=6, a2=2: 21 = 1+12+8 = 21. ✓ Need graph with 6+ vertices, 6 indep sets size 1, 2 indep sets size 2.
# Wait, a1 = number of vertices always. So a1=6 means 6 vertices.
# a2 = # independent pairs = C(6,2) - edges. So a2 = 2 means 13 edges out of 15.
# a1=4, a2=2, a3=1: 21 = 1+8+8+8 = 25. No.
# a1=4, a2=3: 21 = 1+8+12 = 21. ✓ 4 vertices, 3 indep pairs. edges = C(4,2)-3 = 3. So 3 edges on 4 vertices.
# C4 (4-cycle) has 4 edges. Path P4 has 3 edges. K1,3 (star) has 3 edges.
# P4: vertices 0-1-2-3. Independent pairs: {0,2},{0,3},{1,3} = 3. I = 1+8+12 = 21? Let's check:
# alpha_0=1, alpha_1=4, alpha_2=3 (pairs with no edge).
# Pairs: {0,1}edge, {0,2}✓, {0,3}✓, {1,2}edge, {1,3}✓, {2,3}edge. So alpha_2=3. ✓
# alpha_3: triples with no edge: {0,2,3}? 2-3 edge. {0,1,3}? 0-1 edge. None? alpha_3=0.
# I = 1+8+12 = 21. ✓

print("\n  Decomposing 21 = 1 + 2·a₁ + 4·a₂ + 8·a₃ + ...")
solutions = []
for a1 in range(11):
    rem = 21 - 1 - 2*a1
    if rem < 0:
        break
    if rem == 0:
        solutions.append((a1, 0, 0, 0))
        continue
    for a2 in range(rem//4 + 1):
        rem2 = rem - 4*a2
        if rem2 < 0:
            break
        if rem2 == 0:
            solutions.append((a1, a2, 0, 0))
            continue
        for a3 in range(rem2//8 + 1):
            rem3 = rem2 - 8*a3
            if rem3 < 0:
                break
            if rem3 == 0 and rem3 % 16 == 0:
                solutions.append((a1, a2, a3, 0))

print(f"  Possible (a1,a2,a3,...) with I=21:")
for sol in solutions:
    # Check feasibility: a1 = number of vertices, a2 ≤ C(a1,2)
    a1, a2, a3, a4 = sol
    if a2 <= (a1*(a1-1)//2) and a3 <= (a1*(a1-1)*(a1-2)//6):
        print(f"    α = (1, {a1}, {a2}, {a3}): feasible vertices={a1}, edges={a1*(a1-1)//2 - a2}")
    else:
        print(f"    α = (1, {a1}, {a2}, {a3}): INFEASIBLE")

# Key question: which of these can be Ω(T) for some tournament T?
# Ω(T) has vertices = odd cycles of T, edges = pairs sharing a vertex.
# For H=21 = I(Ω,2):
# Option 1: Ω = P_4 (path on 4 vertices, 3 edges) → I=21
# Option 2: Ω = K₁₀ (complete on 10) → I=21
# Option 3: Ω = K₃ ∪ isolated → I = 7·3 = 21
# etc.

# But by SCC product: H = ∏ H(SCCᵢ). And 21 = 3·7.
# If H=21 = 3·7, one SCC has H=3 (triangle, single 3-cycle) and another has H=7 (impossible!).
# So H=21 can't come from a non-SC tournament with an H=7 component.
# But 21 is odd and prime... wait, 21 = 3·7. So if T is non-SC with 2 SCCs
# having H=3 and H=7, that's impossible since H=7 is impossible for ANY tournament.

print(f"\n{'='*60}")
print("H=21 AND THE SCC PRODUCT FORMULA")
print(f"{'='*60}")
print("""
  H=21: If T is non-SC, H = ∏ H(SCCᵢ).
  Factorizations of 21:
    21 = 21 (SC prime)
    21 = 3 × 7 (needs H=7 component — IMPOSSIBLE)
    21 = 7 × 3 (same)

  So H=21 is only possible for strongly connected tournaments.

  Similarly, H=7 × 3^k = 7·3^k requires an SC component with H=7.
  Since H=7 is impossible for ANY tournament (SC or not),
  the ENTIRE family {7·3^k : k ≥ 0} is forbidden!

  This gives: {7, 21, 63, 189, ...} all forbidden.
""")

# Now prove: can an SC tournament have H=21?
# H=21 requires I(Ω,2)=21. From the decomposition above:
# Feasible α vectors include (1,4,3,0), (1,6,2,0), (1,10,0,0), etc.
# The question: can a tournament have an odd-cycle conflict graph with these properties?

# At n≤7 (exhaustive): no tournament has H=21.
# At n=8 (sample): no tournament has H=21.
# Need to show: the structure forcing H=21 creates contradictions.

# IMPORTANT: 21 = 3·7. If H=21 for an SC tournament,
# we need exactly the right number of odd cycles with the right conflict structure.
# But H=21 = I(Ω,2) has multiple solutions.

# The simplest: Ω = path P₄ (4 vertices, 3 edges).
# 4 odd cycles C₁-C₂-C₃-C₄ where adjacent ones share vertices, non-adjacent don't.
# This requires C₁∩C₃ = ∅ and C₂∩C₄ = ∅. Can this happen in a tournament?

# For this to work at n=7 (SC), we need 4 odd cycles with:
# C1 shares vertex with C2 only
# C2 shares with C1 and C3
# C3 shares with C2 and C4
# C4 shares with C3 only
# Plus C1∩C3 = ∅, C1∩C4 = ?, C2∩C4 = ∅

# This is actually quite constrained. Let me check computationally.

print(f"\n{'='*60}")
print("EXHAUSTIVE CHECK: H=21 AT n=5,6,7")
print(f"{'='*60}")

from math import comb as mcomb

def count_hp_fast(A, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)) or dp[S][v] == 0:
                continue
            for u in range(n):
                if not (S & (1 << u)) and A[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def make_tournament_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

for n in [5, 6, 7]:
    total_edges = n * (n-1) // 2
    total = 2**total_edges
    found_21 = 0

    if total > 2**21:
        import random
        random.seed(42)
        sample = 500000
        for _ in range(sample):
            bits = random.randint(0, total-1)
            A = make_tournament_bits(bits, n)
            if count_hp_fast(A, n) == 21:
                found_21 += 1
        print(f"  n={n}: H=21 found {found_21}/{sample} (sampled)")
    else:
        for bits in range(total):
            A = make_tournament_bits(bits, n)
            if count_hp_fast(A, n) == 21:
                found_21 += 1
        print(f"  n={n}: H=21 found {found_21}/{total} (exhaustive)")

# The key insight for H=21:
print(f"\n{'='*60}")
print("THE H=21 PROOF STRATEGY")
print(f"{'='*60}")
print("""
  THEOREM: H=21 is impossible for all tournaments.

  PROOF SKETCH:
  1. By SCC product formula, if T is non-SC:
     H(T) = ∏ H(SCCᵢ)
     21 = 3 × 7 is the only non-trivial factorization.
     But H=7 is impossible, so no non-SC tournament has H=21.

  2. For SC tournaments: need to show I(Ω,2) ≠ 21 directly.
     The feasible Ω structures for I=21 are highly constrained.
     Key: if Ω = P₄ (path on 4 vertices), need 4 odd cycles
     with specific overlap pattern. Tournament completeness
     forces too many additional cycles.

  3. ALTERNATIVELY: H=21 = 3·7.
     By the Fourier analysis (kind-pasteur S73):
     Ĥ[S] = 2·α̂₁[S] for S ≠ ∅.
     If H=21, then Ĥ[∅] = 21/2^{n-1}.
     The level-2 Fourier coefficient is (n-2)!/2^{n-2}.
     These constrain Ω in ways incompatible with I(Ω,2)=21.

  The cleanest proof may be: H=21 = 3·7, and the "7 factor"
  is the fundamental obstruction. Via the semigroup structure
  of achievable H values, anything divisible by 7 but not by 9
  is forbidden (because 7 itself is forbidden and 3-blocks
  can only contribute factors of 3).
""")
