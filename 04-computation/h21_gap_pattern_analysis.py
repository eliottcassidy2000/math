"""
h21_gap_pattern_analysis.py — Analyze the pattern in permanent H-spectrum gaps.

Known permanent gaps: H=7 (w=3), H=21 (w=10).
n=7 gaps that MAY be permanent: w=3, 10 (known), plus others.

Analyze: what's special about w=3 and w=10 from the OCF perspective?
w = alpha_1 + 2*alpha_2 + 4*alpha_3 + ...

For w=3: possible decompositions with alpha_3=0 (Part C):
  (alpha_1, alpha_2) = (1,1) or (3,0)
For w=10:
  (alpha_1, alpha_2) = (2,4), (4,3), (6,2), (8,1), (10,0)

Key insight from Part C: if alpha_3 >= 1, then H >= 1 + 2*(3 + 3*2) + 8 = 27.
Actually Part C says 3 disjoint cycles give H >= 27. So any w that requires
alpha_3=0 and can't be achieved with alpha_3=0 is blocked.

For w=3 with alpha_3=0: need alpha_1+2*alpha_2=3.
  - (1,1): 1 odd cycle, 1 disjoint pair of cycles. But 1 disjoint pair needs
    2 vertex-disjoint cycles, so alpha_1 >= 2. Contradiction.
  - (3,0): 3 pairwise-conflicting cycles. K_3 in Omega.
    THM-029: 3 pairwise-conflicting cycles force alpha_1 >= 4. Contradiction.
So w=3 is IMPOSSIBLE regardless of tournament size!

For w=10 with alpha_3=0: need alpha_1+2*alpha_2=10.
  - Part D eliminates (4,3)
  - Part F eliminates (6,2)
  - Part N eliminates (8,1) at n=8
  - Part M eliminates (10,0) at n=7

  But at larger n, we need the dichotomy or cycle-rich min-H bound.

Let's see if there's a general number-theoretic pattern.

Author: opus-2026-03-07-S43
"""

# The key observation: H = 1 + 2*w where w = sum_{k>=1} alpha_k * 2^{k-1}
# This is the BINARY representation of (H-1)/2 in a special sense.
# But alpha_k can be > 1, so it's not literally binary.

# For H=7: w=3. Binary: 11. Decompositions: 3=1+2, 3=3.
# For H=21: w=10. Binary: 1010. Decompositions: 10=2+8, 10=2+4+4, etc.
# But with 2^0*alpha_1 + 2^1*alpha_2 + 2^2*alpha_3 + ...

# Actually, w = alpha_1 + 2*alpha_2 + 4*alpha_3 + 8*alpha_4 + ...
# This IS a "mixed base" representation with coefficients alpha_k >= 0.
# The key constraint is: alpha_k <= C(alpha_1, k) (at most that many k-subsets
# of the alpha_1 cycles can be independent).

# For alpha_3 = 0 (Part C at H=21): w = alpha_1 + 2*alpha_2.

# The obstruction is: alpha_2 is bounded by the graph structure of Omega(T).
# alpha_2 = number of independent edges in Omega (pairs of disjoint cycles).

# For a graph G on alpha_1 vertices with alpha_2 independent edges:
# Need alpha_1 + 2*alpha_2 = 10 with alpha_3 = 0 (no independent triangle).
# The constraint alpha_3 = 0 means the graph has no 3 pairwise non-adjacent vertices.
# Equivalently, the complement of Omega has no triangle (complement is K_3-free).
# By Ramsey theory: R(3,3)=6, so alpha_1 >= 6 vertices must have an independent
# triple in the complement... wait, that's the wrong direction.

# alpha_3 = 0 means every set of 3 cycles has at least one adjacent pair.
# Equivalently: the independence number of Omega is <= 2.
# By complement: the clique cover number of complement(Omega) is <= 2.
# So complement(Omega) is a bipartite graph (2-colorable).

# Hmm, wait. alpha_3 = 0 means no independent set of size 3.
# So the independence number of Omega is exactly 2 (or 1 or 0).
# By Ramsey: any graph on >= R(3, k+1) = R(3, 3) = 6 vertices has
# either a K_3 or an independent set of size 3.
# So if alpha_1 >= 6 and independence number <= 2:
# Omega must have a triangle (3 pairwise adjacent cycles = 3 pairwise conflicting).
# But actually, Omega can have triangles — we just need no independent set of size 3.
# Ramsey says: independence_number(G) >= 2 needs alpha_1 <= R(3,3)-1 = 5...
# No. Ramsey says: G on n vertices has either K_r or independent set of size s
# when n >= R(r,s). R(3,3) = 6. So any graph on >= 6 vertices has either a
# triangle or an independent set of size 3.

# BUT: Omega CAN have triangles! It's the INDEPENDENCE number we care about.
# If Omega has independence number <= 2, then by Ramsey, it must have a triangle
# (i.e., alpha_1 < 6, OR Omega has both triangles and independence number <= 2).
# Wait: Ramsey says n >= R(3,3) = 6 implies K_3 OR independence_3.
# If independence <= 2 (no independence_3), then K_3 exists (triangle).
# So for alpha_1 >= 6 and alpha_3 = 0: Omega has a triangle (3 mutually
# conflicting cycles). This is fine — it just means 3 cycles share pairwise
# at least one vertex.

# Key: the independence number of Omega is at most 2 when alpha_3 = 0.
# Turan's theorem: max edges of triangle-free graph on n vertices is n^2/4.
# But we need: what's the max alpha_2 for a graph with independence number <= 2?
# alpha_2 = number of edges in the COMPLEMENT of Omega.
# If independence number of Omega <= 2, then clique cover number of complement <= 2.
# So complement is bipartite. Max edges of a bipartite graph on alpha_1 vertices
# is floor(alpha_1^2 / 4) (complete bipartite K_{a,b} with a+b=alpha_1).
# So alpha_2 <= floor(alpha_1^2 / 4).

# But ALSO: alpha_2 = number of independent PAIRS in Omega, which equals
# number of NON-EDGES of Omega. For a graph G:
# alpha_2 = C(alpha_1, 2) - |E(Omega)|.

# For alpha_1 + 2*alpha_2 = 10 with alpha_3 = 0:
# Enumerate:
print("=== (alpha_1, alpha_2) decompositions for w=10 with alpha_3=0 ===")
for a1 in range(0, 11):
    a2_needed = (10 - a1) / 2
    if a2_needed != int(a2_needed) or a2_needed < 0:
        continue
    a2 = int(a2_needed)
    # Feasibility: need a2 <= C(a1, 2) (at most C(a1,2) pairs)
    from math import comb
    if a2 > comb(a1, 2):
        print(f"  ({a1}, {a2}): INFEASIBLE (a2 > C(a1,2)={comb(a1,2)})")
        continue
    # Need independence number <= 2: complement is bipartite
    # alpha_2 = non-edges of Omega = edges of complement
    # complement is bipartite on a1 vertices with a2 edges
    # Max edges of bipartite graph on a1 vertices: floor(a1^2/4)
    max_bip = a1 * a1 // 4
    if a2 > max_bip:
        print(f"  ({a1}, {a2}): INFEASIBLE by Turan (a2={a2} > floor(a1^2/4)={max_bip})")
        continue
    # Also need: this graph can be realized as Omega(T) for some tournament
    print(f"  ({a1}, {a2}): feasible graph-theoretically, a2/C(a1,2)={a2}/{comb(a1,2)}")

# Now do the same for w=3 (H=7)
print("\n=== (alpha_1, alpha_2) decompositions for w=3 with alpha_3=0 ===")
for a1 in range(0, 4):
    a2_needed = (3 - a1) / 2
    if a2_needed != int(a2_needed) or a2_needed < 0:
        continue
    a2 = int(a2_needed)
    from math import comb
    if a2 > comb(a1, 2):
        print(f"  ({a1}, {a2}): INFEASIBLE (a2 > C(a1,2)={comb(a1,2)})")
        continue
    max_bip = a1 * a1 // 4
    if a2 > max_bip:
        print(f"  ({a1}, {a2}): INFEASIBLE by independence (a2={a2} > {max_bip})")
        continue
    print(f"  ({a1}, {a2}): feasible graph-theoretically")

# KEY INSIGHT: For w=3, decompositions are (1,1) and (3,0).
# (1,1) needs alpha_1=1 and alpha_2=1: 1 cycle with 1 disjoint pair.
# But alpha_2 = 1 means 1 pair of vertex-disjoint cycles.
# That requires alpha_1 >= 2. Contradiction!
# (3,0) needs 3 cycles, all pairwise conflicting (sharing vertex).
# THM-029 shows this forces a 4th cycle. So alpha_1 >= 4. Contradiction.

# For w=10: the analysis is more complex but follows the same pattern.
# Each decomposition is blocked by tournament-specific graph obstructions.

print("\n=== Can we find a UNIFIED obstruction? ===")
print("w=3: blocked by (1,1) infeasible + (3,0) K3 forcing")
print("w=10: blocked by tournament 5-cycle forcing + capacity bounds")
print()
print("Observation: w=3 and w=10 are both values where EVERY")
print("decomposition is simultaneously blocked. For other values,")
print("at least one decomposition is achievable.")

# Let's check: for w=4 (H=9), which decompositions are achievable?
print("\n=== w=4 (H=9) decompositions ===")
for a1 in range(0, 5):
    a2 = (4 - a1) / 2
    if a2 != int(a2) or a2 < 0:
        continue
    a2 = int(a2)
    print(f"  ({a1}, {int(a2)}): a1={a1}, a2={a2}")
    # (0,2): 0 cycles, 2 disjoint pairs — impossible (need cycles first)
    # (2,1): 2 cycles, 1 disjoint pair — possible! 2 disjoint 3-cycles.
    # (4,0): 4 mutually conflicting cycles — possible (star configuration)
# So w=4 is achievable via (2,1) or (4,0). Indeed H=9 exists at n=5.
