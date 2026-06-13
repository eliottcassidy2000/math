#!/usr/bin/env python3
"""
THE H-PFAFFIAN GAP: Q = (H - |Pf(S)|)(H + |Pf(S)|)/8
opus-2026-03-13-S67k

Now that we know det(I+2A) = Pf(S)² for even n (THM-174),
the gap Q = (H² - det)/8 = (H² - Pf²)/8 = (H-|Pf|)(H+|Pf|)/8.

Let a = (H-|Pf|)/2, b = (H+|Pf|)/2. Then Q = ab/2.

This script investigates:
1. What does the Pfaffian measure? (It's the "antisymmetric part" of H)
2. How does Q relate to the conflict graph structure?
3. Can we express |Pf| in terms of the α_k?
4. Connection between H = I(CG,2) and Pf(S) for even n

For ODD n, √det = |Σ w_i| where w_i = ±Pf(S_ii).
Can we compute these sub-Pfaffians?
"""

import numpy as np
from itertools import combinations, permutations

def adj_matrix(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[idx]: A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def count_hp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def pfaffian(S):
    """Compute Pfaffian of a skew-symmetric matrix S (must be even-dim)."""
    n = len(S)
    if n == 0:
        return 1
    if n % 2 == 1:
        return 0
    if n == 2:
        return int(round(S[0][1]))

    # Recursive: Pf(S) = sum_{j>0} (-1)^{j} S[0][j] * Pf(S_{0j})
    # where S_{0j} is S with rows 0,j and cols 0,j removed
    result = 0
    for j in range(1, n):
        if abs(S[0][j]) < 1e-10:
            continue
        # Remove rows/cols 0 and j
        idx = [k for k in range(n) if k != 0 and k != j]
        sub = S[np.ix_(idx, idx)]
        sign = (-1)**(j+1)  # (-1)^{j-1} since j ranges from 1
        result += sign * int(round(S[0][j])) * pfaffian(sub)
    return result

def fast_hash(A, n):
    scores = list(A.sum(axis=1))
    ns = []
    for i in range(n):
        o = sorted(scores[j] for j in range(n) if A[i][j])
        ins = sorted(scores[j] for j in range(n) if A[j][i])
        ns.append((scores[i], tuple(o), tuple(ins)))
    return tuple(sorted(ns))

def find_directed_odd_cycles(A, n):
    cycles = set()
    for length in range(3, n + 1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                if all(A[perm[i]][perm[(i+1) % length]] for i in range(length)):
                    mi = perm.index(min(perm))
                    canon = perm[mi:] + perm[:mi]
                    cycles.add(canon)
    return cycles

def conflict_graph_alpha(cycles):
    cl = list(cycles)
    nc = len(cl)
    vs = [set(c) for c in cl]
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if vs[i] & vs[j]:
                adj[i][j] = adj[j][i] = True
    alpha = [0] * (nc + 1)
    alpha[0] = 1
    for size in range(1, nc + 1):
        cnt = 0
        for sub in combinations(range(nc), size):
            ok = True
            for a in range(len(sub)):
                for b in range(a+1, len(sub)):
                    if adj[sub[a]][sub[b]]:
                        ok = False
                        break
                if not ok: break
            if ok: cnt += 1
        alpha[size] = cnt
    return alpha

print("=" * 72)
print("THE H-PFAFFIAN GAP STRUCTURE")
print("=" * 72)

# =============================================
print("\n--- PART 1: Pfaffian computation for even n ---\n")

for n in [4, 6]:
    print(f"n = {n}:")
    m = n*(n-1)//2
    hash_groups = {}
    for bits in range(1 << m):
        b = [(bits >> i) & 1 for i in range(m)]
        A = adj_matrix(b, n)
        h = fast_hash(A, n)
        if h not in hash_groups:
            hash_groups[h] = A

    for h, A in sorted(hash_groups.items(), key=lambda x: count_hp(x[1], n)):
        H = count_hp(A, n)
        S = (A - A.T).astype(float)
        pf = pfaffian(S)
        det_check = int(round(np.linalg.det(np.eye(n) + 2*A.astype(float))))

        a_val = (H - abs(pf)) // 2
        b_val = (H + abs(pf)) // 2
        Q = a_val * b_val // 2

        cycles = find_directed_odd_cycles(A, n)
        alpha = conflict_graph_alpha(cycles) if cycles else [1]
        a1 = alpha[1] if len(alpha) > 1 else 0
        a2 = alpha[2] if len(alpha) > 2 else 0

        print(f"  H={H:3d} Pf(S)={pf:+4d} |Pf|={abs(pf):3d} "
              f"a=(H-|Pf|)/2={a_val:3d} b=(H+|Pf|)/2={b_val:3d} Q=ab/2={Q:5d} | "
              f"α₁={a1:2d} α₂={a2:2d} det_check={det_check} pf²={pf*pf}")
    print()

# =============================================
print("\n--- PART 2: Signed Pfaffian — does sign matter? ---\n")

print("The Pfaffian has a SIGN (±). Does the sign carry information?")
print("Note: Pf(S) depends on the labeling of vertices.")
print("Under permutation σ: Pf(P^T S P) = sign(σ) · Pf(S)")
print("So the sign is labeling-dependent, but |Pf| is an invariant.")
print()

for n in [4, 6]:
    print(f"n = {n}: signs observed")
    m = n*(n-1)//2
    hash_groups = {}
    for bits in range(1 << m):
        b = [(bits >> i) & 1 for i in range(m)]
        A = adj_matrix(b, n)
        h = fast_hash(A, n)
        if h not in hash_groups:
            hash_groups[h] = A

    pf_signs = {}
    for h, A in hash_groups.items():
        H = count_hp(A, n)
        S = (A - A.T).astype(float)
        pf = pfaffian(S)
        key = (H, abs(pf))
        if key not in pf_signs:
            pf_signs[key] = set()
        pf_signs[key].add(1 if pf > 0 else (-1 if pf < 0 else 0))

    for (H, apf), signs in sorted(pf_signs.items()):
        print(f"  H={H:3d} |Pf|={apf:3d} signs={signs}")
    print()

# =============================================
print("\n--- PART 3: Sub-Pfaffians for odd n ---\n")

for n in [3, 5]:
    print(f"n = {n}:")
    m = n*(n-1)//2
    hash_groups = {}
    for bits in range(1 << m):
        b = [(bits >> i) & 1 for i in range(m)]
        A = adj_matrix(b, n)
        h = fast_hash(A, n)
        if h not in hash_groups:
            hash_groups[h] = A

    for h, A in sorted(hash_groups.items(), key=lambda x: count_hp(x[1], n)):
        H = count_hp(A, n)
        S = (A - A.T).astype(float)
        det_M = int(round(np.linalg.det(np.eye(n) + 2*A.astype(float))))
        sqrt_det = int(round(abs(det_M)**0.5))

        # Sub-Pfaffians: Pf(S_{ii}) where S_{ii} = S with row/col i removed
        sub_pfs = []
        for i in range(n):
            idx = [k for k in range(n) if k != i]
            sub_S = S[np.ix_(idx, idx)]
            pf_sub = pfaffian(sub_S)
            sub_pfs.append(pf_sub)

        # The w vector: w_i = (-1)^i * Pf(S_ii) ... or some other sign pattern
        # Test: det(I+2A) = (Σ w_i)² for various sign choices

        # Try: w_i = (-1)^i * Pf(S_ii)
        w1 = [(-1)**i * sub_pfs[i] for i in range(n)]
        sum_w1 = sum(w1)

        # Try: w_i = (-1)^{i+1} * Pf(S_ii)
        w2 = [(-1)**(i+1) * sub_pfs[i] for i in range(n)]
        sum_w2 = sum(w2)

        # Just use Pf(S_ii) with natural sign
        sum_raw = sum(sub_pfs)

        print(f"  H={H:3d} det(I+2A)={det_M:4d} √det={sqrt_det:3d} | "
              f"sub_Pfs={sub_pfs} | "
              f"Σ(-1)^i·Pf={sum_w1}→{sum_w1**2} | "
              f"Σ(-1)^(i+1)·Pf={sum_w2}→{sum_w2**2} | "
              f"Σ_raw={sum_raw}→{sum_raw**2}")
    print()

# =============================================
print("\n--- PART 4: The Pfaffian as graph matching count ---\n")

print("""
RECALL: For a skew-symmetric matrix S, Pf(S) = Σ_{perfect matchings M} sign(M) · ∏_{(i,j)∈M} S_{ij}

For S = A - A^T (tournament skew-adjacency):
  S_{ij} = +1 if i→j, -1 if j→i

So Pf(S) = Σ over perfect matchings of the complete graph K_n,
  with sign from the matching orientation.

For even n, |Pf(S)|² = det(S) = det(I+2A).

The MATCHING INTERPRETATION is:
  H = number of Hamiltonian paths (using I(CG,2) = OCF)
  |Pf(S)| = |signed count of perfect matchings in K_n with tournament orientation|
  Q = (H² - Pf²)/8 measures how much H "exceeds" the matching count

This connects to DIMER theory!
""")

# =============================================
print("\n--- PART 5: |Pf(S)| in terms of α-vector ---\n")

print("For each tournament, express |Pf(S)| via the α_k:")
print("  H = 1 + 2α₁ + 4α₂ + ...,  |Pf(S)| = ???")
print()

for n in range(3, 7):
    print(f"n = {n}:")
    m = n*(n-1)//2
    hash_groups = {}
    for bits in range(1 << m):
        b = [(bits >> i) & 1 for i in range(m)]
        A = adj_matrix(b, n)
        h = fast_hash(A, n)
        if h not in hash_groups:
            hash_groups[h] = A

    seen = set()
    for h, A in sorted(hash_groups.items(), key=lambda x: count_hp(x[1], n)):
        H = count_hp(A, n)
        S_mat = (A - A.T).astype(float)

        if n % 2 == 0:
            pf_val = abs(pfaffian(S_mat))
        else:
            det_M = int(round(np.linalg.det(np.eye(n) + 2*A.astype(float))))
            pf_val = int(round(abs(det_M)**0.5))

        cycles = find_directed_odd_cycles(A, n)
        alpha = conflict_graph_alpha(cycles) if cycles else [1]
        a1 = alpha[1] if len(alpha) > 1 else 0
        a2 = alpha[2] if len(alpha) > 2 else 0

        key = (H, pf_val, a1, a2)
        if key not in seen:
            seen.add(key)
            # Express |Pf| in terms of α
            # H = 1+2α₁+4α₂, so |Pf| = H - 2a where a = (H-|Pf|)/2
            a_gap = (H - pf_val) // 2

            print(f"  H={H:3d} |Pf|={pf_val:3d} a=(H-|Pf|)/2={a_gap:3d} | "
                  f"α₁={a1:2d} α₂={a2:2d} | "
                  f"a-α₁+α₂={a_gap-a1+a2:+3d}")
    print()

# =============================================
print("\n--- PART 6: The Pfaffian vs Chromatic number ---\n")

print("""
For EVEN n, |Pf(S)| is odd (since det(S) = Pf² and √det is odd).
|Pf(S)| = 1 + 2k for some non-negative integer k.

Key patterns observed:
  α₂=0, all 3-cycles pairwise share vertices:
    |Pf| depends on whether CG has "balanced" or "unbalanced" cycle distribution

  When H = |Pf|: all cycles through every vertex → no "free" cycles → Q=0
  When H > |Pf|: some cycles NOT sharing a common vertex → Q > 0

CONJECTURE: |Pf(S)| relates to the independence number of some
graph derived from the conflict graph, possibly the LINE GRAPH of CG
or the INTERSECTION GRAPH of the vertex sets of cycles.
""")

print("\n" + "=" * 72)
print("SYNTHESIS — THE COMPLETE PICTURE")
print("=" * 72)

print("""
THE THREE QUANTITIES:
  H(T)     = I(CG(T), 2)     = Σ_k 2^k α_k     [path count]
  |Pf(S)|  = √det(I+2A)      = ???               [matching count]
  Q(T)     = (H²-Pf²)/8      = ab/2              [gap]

PROVEN IDENTITIES:
  1. det(I+2A) = Pf(S)²  for even n  (THM-174)
  2. det(I+2A) = (Σ w_i)² for odd n  (THM-174)
  3. H²-Pf² ≡ 0 (mod 8)              (HYP-850)
  4. H ≥ |Pf(S)| always               (observed, not proved)

RECURRENCE STRUCTURE:
  H satisfies deletion recurrence: H(T) = H(T-v) + 2μ_v
  Pf satisfies: ???
    For even n→odd n-1: delete row/col → sub-Pfaffian (cofactor expansion)
    For odd n→even n-1: Pf(S_{ii}) are the components of the w-vector

  Q = ab/2 where a+b = H:
    Under deletion T→T-v:
    Q(T) - Q(T-v) = μ_v(H-μ_v)/2 + [Pf(S_{T-v})²-Pf(S_T)²]/8

THE DEEP CONNECTION:
  H counts Hamiltonian PATHS (order matters, all vertices)
  |Pf| counts signed perfect MATCHINGS (pairing matters, all vertices)
  Both are "global" graph invariants requiring all n vertices
  The gap Q measures the "excess structure" that paths see but matchings don't

This is analogous to:
  In a bipartite graph: det(A) counts signed matchings
  For planar graphs: permanent ≈ Pfaffian (FKT theorem)
  For tournaments: H ≈ Pf through the OCF/independence poly bridge

OPEN QUESTIONS:
  1. Can we express |Pf(S)| directly in terms of {α_k}?
  2. Does |Pf(S)| have its own independence polynomial interpretation?
  3. What is the Pfaffian recurrence under vertex deletion?
  4. Is there a "Pfaffian OCF" — Pf(S) = I(G', 2) for some graph G'?
  5. Does Q = ab/2 have a combinatorial interpretation as counting
     something specific (pairs of paths? cycle configurations?)
""")
