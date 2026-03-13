#!/usr/bin/env python3
"""
ALGEBRAIC PROOF that tr(A^5) = 5*c5_dir for tournaments.

Claim: In a tournament, every closed walk of length 5 visits 5 distinct vertices.

Proof: Let w = (i0, i1, i2, i3, i4) be a closed walk (i4→i0, and each i_t→i_{t+1}).
We need A_{i0,i1} = A_{i1,i2} = A_{i2,i3} = A_{i3,i4} = A_{i4,i0} = 1.

Since consecutive vertices must differ (A_{ii}=0), we need each i_t ≠ i_{t+1 mod 5}.

Suppose some non-consecutive pair coincides: i_a = i_b with b - a ≥ 2 (mod 5).

The possible non-consecutive pairs (gap ≥ 2):
  (i0, i2): gap 2. Walk includes i0→i1→i2=i0. Need A_{i1,i0}=1. But also A_{i0,i1}=1.
             So A_{i0,i1}=1 AND A_{i1,i0}=1, contradicting tournament axiom.

  (i1, i3): gap 2. Walk includes i1→i2→i3=i1. Need A_{i2,i1}=1 AND A_{i1,i2}=1. ✗

  (i2, i4): gap 2. Walk includes i2→i3→i4=i2. Need A_{i3,i2}=1 AND A_{i2,i3}=1. ✗

  (i3, i0): gap 2 (since (i3→i4→i0) and i3=i0 means i0→i4→i0). Wait:
             i3 = i0 means the walk is i0→i1→i2→i0→i4→i0.
             Need A_{i4,i0}=1 (last step) AND A_{i0,i4}=1 (step i3→i4 = i0→i4).
             So A_{i0,i4}=1 AND A_{i4,i0}=1. ✗

  (i4, i1): gap 2 (mod 5, since 1-4 = -3 ≡ 2 mod 5). Walk: ...i4=i1→i0→i1...
             Actually the walk is i0→i1→i2→i3→i1→i0 (since i4=i1).
             Step i3→i4=i1: A_{i3,i1}=1. Step i4=i1→i0: A_{i1,i0}=1.
             But also step i0→i1: A_{i0,i1}=1. And step i1→i0: A_{i1,i0}=1.
             So A_{i0,i1}=1 AND A_{i1,i0}=1. ✗

So ALL gap-2 pairs lead to A_{x,y}=A_{y,x}=1, which is impossible in a tournament.

For gap ≥ 3 in a 5-cycle, the only remaining case is gap=3 (since 5-3=2 means the
"other direction" has gap 2, which is already covered).

Wait, let me reconsider. In a 5-step closed walk (i0,...,i4,i0), the possible
coincidences of non-adjacent vertices are:
  (i0,i2), (i1,i3), (i2,i4), (i3,i0), (i4,i1): all gap-2
  (i0,i3), (i1,i4): gap-3

Let me check gap-3:
  (i0, i3): Walk is i0→i1→i2→i0→i4→i0. Steps:
    i2→i3=i0: A_{i2,i0}=1
    i3=i0→i4: A_{i0,i4}=1
    i4→i0: A_{i4,i0}=1
    So A_{i0,i4}=1 AND A_{i4,i0}=1. ✗

  (i1, i4): Walk is i0→i1→i2→i3→i1→i0. Steps:
    i3→i4=i1: A_{i3,i1}=1
    i4=i1→i0: A_{i1,i0}=1
    i0→i1: A_{i0,i1}=1
    So A_{i0,i1}=1 AND A_{i1,i0}=1. ✗

AMAZING! Every possible vertex coincidence in a closed 5-walk leads to a contradiction!

The argument generalizes: in a closed walk of length k on a tournament,
if i_a = i_b with b > a, then the walk contains the sub-walk
i_a → i_{a+1} → ... → i_{b-1} → i_b = i_a (length b-a)
AND the sub-walk i_b → i_{b+1} → ... → i_{k-1} → i_0 → ... → i_a (length k-b+a)

If b-a = 1 or k-(b-a) = 1: consecutive pair, excluded by A_{ii}=0.
If b-a = 2: gap-2, forces A_{i_{a+1},i_a}=1 AND A_{i_a,i_{a+1}}=1. ✗

For k=5: b-a ∈ {2,3} (excluding 1 and 4=k-1). If b-a=3: k-(b-a)=2, same issue.

For k=7: b-a ∈ {2,3,4,5}. If b-a=2 or k-(b-a)=2: contradiction.
  b-a=2: ✗. b-a=5: k-(b-a)=2 ✗.
  b-a=3: k-(b-a)=4. Neither gap is 1 or 2, so NO contradiction!
  b-a=4: k-(b-a)=3. Same: neither gap is 1 or 2. NO contradiction.

So for k=7, vertex coincidences with gap 3 or 4 are possible.
This explains why tr(A^7) ≠ 7*c7_dir but tr(A^5) = 5*c5_dir!

THEOREM: For any tournament T on n vertices:
  tr(A^k) = k * c_{k,dir}  if and only if  min(j, k-j) ≤ 2 for all 2 ≤ j ≤ k-2
  i.e., k ≤ 5 (since for k=6: j=3 gives min(3,3)=3 > 2).

Wait, k=6 is EVEN so doesn't count for odd cycles. Let me reconsider.

For odd k: we need min(j, k-j) ≤ 2 for all j with 2 ≤ j ≤ k-2.
  k=3: j ∈ {}: trivially satisfied. (Actually j=2: min(2,1)=1≤2 ✓)
  k=5: j ∈ {2,3}: min(2,3)=2 ✓, min(3,2)=2 ✓.
  k=7: j ∈ {2,3,4,5}: min(3,4)=3 > 2. FAILS.

So: tr(A^k) = k*c_{k,dir} for k ∈ {3, 5} but NOT for k ≥ 7.
This is a COMPLETE characterization!

opus-2026-03-13-S71c
"""

# Let's also verify for k=4 and k=6 (even cycles)
import numpy as np
from itertools import combinations, permutations

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

print("=" * 60)
print("THEOREM: tr(A^k) = k * c_{k,dir} iff k <= 5")
print("=" * 60)

# Note: for tournaments, c_{even,dir} = 0 trivially (no even directed cycles in tournaments).
# Wait — that's not true! Tournaments CAN have directed 4-cycles.
# E.g.: 0→1→2→3→0 is a 4-cycle if A_{01}=A_{12}=A_{23}=A_{30}=1.

# Actually, do tournaments have even directed cycles? Yes! A tournament on 4 vertices
# with cycle 0→1→2→3→0 is a valid tournament (other edges: 0→2 or 2→0, 1→3 or 3→1).

# Check tr(A^4) at n=4:
n = 4
tb = n*(n-1)//2
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    tr4 = int(np.trace(np.linalg.matrix_power(A, 4)))
    c4 = 0
    for perm in permutations(range(1, n)):
        path = (0,) + perm
        valid = all(A[path[i]][path[(i+1) % n]] for i in range(n))
        if valid: c4 += 1
    if tr4 != 4 * c4:
        print(f"  n=4 bits={bits}: tr4={tr4}, 4*c4={4*c4}")
        break
else:
    print(f"  n=4: tr(A^4) = 4*c4 for ALL 2^6=64 tournaments")

# Check k=6 at n=6:
n = 6
tb = n*(n-1)//2
np.random.seed(42)
fail6 = 0
for trial in range(1000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    tr6 = int(np.trace(np.linalg.matrix_power(A, 6)))
    c6 = 0
    for perm in permutations(range(1, n)):
        path = (0,) + perm
        valid = all(A[path[i]][path[(i+1) % n]] for i in range(n))
        if valid: c6 += 1
    if tr6 != 6 * c6:
        fail6 += 1
        if fail6 == 1:
            print(f"  n=6 first fail: tr6={tr6}, 6*c6={6*c6}, diff={tr6-6*c6}")

print(f"  n=6: tr(A^6) = 6*c6 failures: {fail6}/1000")

# For even k, the analysis is different because tournaments CAN have even directed cycles.
# The gap argument: closed k-walk, i_a = i_b with gap j = b-a.
# Sub-walk of length j from i_a to itself, sub-walk of length k-j from i_a to itself.
# If j=1 or k-j=1: A_{ii}=0. ✗
# If j=2 or k-j=2: A_{xy}=A_{yx}=1. ✗
# So for k=4: j ∈ {2}: min(2,2)=2 ≤ 2. ✓! No degenerate walks.
# For k=6: j ∈ {2,3,4}: j=3 gives min(3,3)=3 > 2. FAILS.

# So the complete answer for ALL k (even and odd):
# tr(A^k) = k * c_{k,dir}  iff  k ≤ 5
# Because for k ≤ 5, every gap j with 2 ≤ j ≤ k-2 has min(j,k-j) ≤ 2.
# For k = 6: j=3 fails. For k ≥ 7: j=3 fails.

print(f"\n  THEOREM: tr(A^k) = k * c_{{k,dir}} for all tournaments iff k <= 5")
print(f"  Proof: vertex repetition in closed k-walk requires gap j with")
print(f"  min(j, k-j) >= 3, which needs k >= 6.")
print(f"  For k=3: no gap j with 2 <= j <= 1 (empty). ✓")
print(f"  For k=4: gap j=2 → min(2,2)=2. ✓")
print(f"  For k=5: gaps j=2,3 → min(2,3)=2, min(3,2)=2. ✓")
print(f"  For k=6: gap j=3 → min(3,3)=3 > 2. ✗")
print(f"  For k>=7: gap j=3 → min(3,k-3) >= min(3,4) = 3 > 2. ✗")

# Now: COROLLARY for c5 being lambda-determined.
# Since tr(A^5) = 5*c5_dir, and tr(A^5) = tr(A · A^4) = tr(A² · A · A²),
# we need to show tr(A^5) is a function of the labeled lambda graph.

# Key: (A²)_{ij} for i≠j is NOT determined by lambda alone.
# But tr(A^5) = Σ_i Σ_j A_{ij} (A^4)_{ji} might still be lambda-determined
# due to cancellation.

# Alternative approach: tr(A^5) = Σ_{i,j,k,l,m distinct} A_{ij}A_{jk}A_{kl}A_{lm}A_{mi}
# (since all 5 vertices must be distinct!)
# This is a polynomial in the adjacency entries.
# Since A_{ij} = (1 + B_{ij})/2 where B_{ij} ∈ {±1} is skew-symmetric:
# The sum can be rewritten in terms of B entries, which are ±1.

# Now: λ(i,j) = #{3-cycles through i,j} = #{w: A_{ij}A_{jw}A_{wi} + A_{ji}A_{iw}A_{wj} > 0}
# More precisely: λ(i,j) = Σ_w (A_{ij}A_{jw}A_{wi} + A_{ji}A_{iw}A_{wj})

# Can tr(A^5) be expressed as a polynomial in λ values?
# tr(A^5) = Σ_{cyc(i,j,k,l,m)} A_{ij}A_{jk}A_{kl}A_{lm}A_{mi}
# Each term is a product of 5 adjacency entries on a directed cycle.

# This is the permanent-like sum over directed Hamilton cycles on 5-vertex subsets.
# It should be expressible via the lambda graph because:
# 1. Lambda determines the LOCAL structure (3-cycles) through each pair
# 2. A 5-cycle can be decomposed into overlapping 3-cycles in a specific way

# Actually, think about it this way: at n=5, the tournament is determined
# (up to isomorphism) by its lambda graph, because lambda captures enough
# information about the adjacency structure.

# Wait — is that true? Two non-isomorphic n=5 tournaments can have the same
# labeled lambda graph? Let me check.

n = 5
tb = n*(n-1)//2
lam_to_bits = {}
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1; L[v][u] += 1
    key = tuple(L[i][j] for i in range(n) for j in range(i+1, n))
    if key in lam_to_bits:
        # Found two tournaments with same labeled lambda — check if they differ
        A_prev = bits_to_adj(lam_to_bits[key], n)
        if not np.array_equal(A, A_prev):
            print(f"  TWO n=5 tournaments with SAME labeled lambda:")
            print(f"    bits1={lam_to_bits[key]}: A={A_prev.tolist()}")
            print(f"    bits2={bits}: A={A.tolist()}")
            print(f"    lambda={key}")
            break
    else:
        lam_to_bits[key] = bits
else:
    print(f"\n  At n=5: labeled lambda UNIQUELY determines the tournament!")
    print(f"  ({len(lam_to_bits)} distinct labeled lambda graphs = {1<<tb} tournaments)")
    print(f"  This trivially proves c5 is lambda-determined at n=5.")

# At n=6: does labeled lambda uniquely determine the tournament?
n = 6
tb = n*(n-1)//2
lam_to_bits_6 = {}
non_unique = 0
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1; L[v][u] += 1
    key = tuple(L[i][j] for i in range(n) for j in range(i+1, n))
    if key in lam_to_bits_6:
        A_prev = bits_to_adj(lam_to_bits_6[key], n)
        if not np.array_equal(A, A_prev):
            non_unique += 1
    else:
        lam_to_bits_6[key] = bits

print(f"\n  At n=6: non-unique labeled lambda pairs: {non_unique}")
if non_unique == 0:
    print(f"  Labeled lambda UNIQUELY determines the tournament at n=6 too!")
    print(f"  ({len(lam_to_bits_6)} distinct labeled lambda graphs = {1<<tb} tournaments)")
    print(f"  This trivially proves c5 (and c7, and EVERYTHING) is lambda-determined at n=6.")
else:
    print(f"  Labeled lambda does NOT uniquely determine tournaments at n=6.")

# At n=7: same check (sampled)
n = 7
tb = n*(n-1)//2
np.random.seed(42)
lam_to_adj_7 = {}
non_unique_7 = 0
for trial in range(100000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1; L[v][u] += 1
    key = tuple(L[i][j] for i in range(n) for j in range(i+1, n))
    if key in lam_to_adj_7:
        if not np.array_equal(A, lam_to_adj_7[key]):
            non_unique_7 += 1
            if non_unique_7 <= 3:
                A_prev = lam_to_adj_7[key]
                # Check if they have same c5 and c7
                tr5_1 = int(np.trace(np.linalg.matrix_power(A, 5)))
                tr5_2 = int(np.trace(np.linalg.matrix_power(A_prev, 5)))
                tr7_1 = int(np.trace(np.linalg.matrix_power(A, 7)))
                tr7_2 = int(np.trace(np.linalg.matrix_power(A_prev, 7)))
                print(f"\n  n=7 non-unique pair #{non_unique_7}:")
                print(f"    tr5: {tr5_1} vs {tr5_2} (same: {tr5_1==tr5_2})")
                print(f"    tr7: {tr7_1} vs {tr7_2} (same: {tr7_1==tr7_2})")
    else:
        lam_to_adj_7[key] = A.copy()

print(f"\n  At n=7 (100k samples): non-unique labeled lambda: {non_unique_7}")

print(f"\nDone.")
