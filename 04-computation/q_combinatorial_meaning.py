#!/usr/bin/env python3
"""
Q = (H² - det(I+2A))/8 — WHAT DOES IT COUNT?
opus-2026-03-13-S67k

From the deletion analysis, Q has these values at n=5:
  H= 1, det=  1: Q=  0
  H= 3, det=  9: Q=  0
  H= 5, det=  1: Q=  3
  H= 9, det=  9: Q=  9
  H= 9, det= 25: Q=  7
  H=11, det=  1: Q= 15
  H=13, det= 49: Q= 15
  H=15, det= 81: Q= 18
  H=15, det= 25: Q= 25

Key observation: Q = 0 iff H = √det (iff α₁ ≤ 1 from earlier).

Since H = 1 + 2α₁ + 4α₂ + ... and √det is also odd,
let S = √det. Then:
  Q = (H² - S²)/8 = (H-S)(H+S)/8

Both H-S and H+S are even (both odd), so (H-S)/2 and (H+S)/2 are integers.
Q = (H-S)(H+S)/8 = [(H-S)/2 · (H+S)/2] / 2

So Q is integer iff (H-S)/2 · (H+S)/2 is even, which we proved.

Let a = (H-S)/2, b = (H+S)/2. Then a+b = H, b-a = S.
Q = a·b/2. So Q = a·b/2 where a+b = H.

This means Q counts SOMETHING related to choosing from H-related sets!

Let's investigate: what are a and b for each tournament?
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict

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

def det_I2A(A, n):
    M = np.eye(n, dtype=float) + 2 * A.astype(float)
    return int(round(np.linalg.det(M)))

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

def conflict_graph_indpoly(cycles):
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
print("Q = (H² - det)/8 = a·b/2 WHERE a = (H-√det)/2, b = (H+√det)/2")
print("=" * 72)

print("\n--- PART 1: a, b decomposition ---\n")

data_all = []
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

    for h, A in sorted(hash_groups.items(), key=lambda x: count_hp(x[1], n)):
        H = count_hp(A, n)
        d = det_I2A(A, n)
        S = int(round(abs(d)**0.5))
        if S*S != abs(d):
            print(f"  WARNING: det={d} not perfect square!")
            continue

        a = (H - S) // 2
        b = (H + S) // 2
        Q = a * b // 2

        cycles = find_directed_odd_cycles(A, n)
        alpha = conflict_graph_indpoly(cycles) if cycles else [1]
        a1 = alpha[1] if len(alpha) > 1 else 0
        a2 = alpha[2] if len(alpha) > 2 else 0

        data_all.append((n, H, S, a, b, Q, a1, a2, tuple(alpha)))

        print(f"  H={H:3d} √det={S:3d} a=(H-S)/2={a:3d} b=(H+S)/2={b:3d} "
              f"Q=ab/2={Q:5d} α₁={a1:2d} α₂={a2:2d}")
    print()

print("\n--- PART 2: What is a = (H-√det)/2? ---\n")

print("Hypothesis: a = (H-√det)/2 counts something related to cycles.")
print("Since H = Σ 2^k α_k and √det = ???")
print()
print("Let's express √det in terms of α_k:")

for n, H, S, a, b, Q, a1, a2, alpha in data_all:
    # If α₂ = 0, then H = 1 + 2α₁, so S = H - 2a = 1+2α₁-2a
    # If S = 1 (the simplest case), then a = (H-1)/2 = α₁
    if a2 == 0:
        # When no α₂ channels: H = 1 + 2α₁
        # a = (H-S)/2
        s_from_a = H - 2*a
        print(f"  n={n} H={H:3d} α₁={a1:2d} a={a:3d} √det={S:3d} | "
              f"√det = H-2a = {H}-{2*a} = {s_from_a} | "
              f"a/α₁ = {a/a1:.3f}" if a1 > 0 else
              f"  n={n} H={H:3d} α₁={a1:2d} a={a:3d} √det={S:3d} | transitive")

print("\n--- PART 3: Is Q related to C(a₁, 2) or similar? ---\n")

print("Testing: Q vs C(a,2)=a(a-1)/2, Q vs a², Q vs C(α₁,2)=α₁(α₁-1)/2:")
for n, H, S, a, b, Q, a1, a2, alpha in data_all:
    ca2 = a*(a-1)//2
    a_sq = a*a
    ca1_2 = a1*(a1-1)//2
    ab_half = a * b // 2  # should equal Q

    print(f"  n={n} H={H:3d} a={a:3d} b={b:3d} Q={Q:5d} | "
          f"C(a,2)={ca2:5d} a²={a_sq:5d} C(α₁,2)={ca1_2:5d} | "
          f"Q-C(a,2)={Q-ca2:5d} Q-a²={Q-a_sq:5d}")

print("\n--- PART 4: a in terms of cycles ---\n")

print("What is a = (H-√det)/2 for tournaments with known √det?")
print()
print("For transitive tournament: H=1, det=1, S=1, a=0, b=1, Q=0")
print("For 3-cycle tournament: H=3, det=9, S=3, a=0, b=3, Q=0")
print("  So α₁=1 but a=0! This means a ≠ α₁ in general.")
print()
print("Wait — the unique tournament on 4 vertices with H=5:")
print("  H=5, det=1, S=1, a=2, b=3, Q=3")
print("  α₁=2, a=2, so a=α₁ here!")
print()

print("KEY PATTERN:")
for n, H, S, a, b, Q, a1, a2, alpha in data_all:
    diff_a_alpha1 = a - a1
    print(f"  n={n} H={H:3d} √det={S:3d} a={a:3d} α₁={a1:2d} a-α₁={diff_a_alpha1:+4d} | α₂={a2}")

print("\n--- PART 5: √det as independence polynomial of SOMETHING ---\n")

print("Recall H = I(CG, 2). Is √det also an independence polynomial?")
print()
print("Let's check: for each tournament, is √det of the form 1 + 2k for some k?")
for n, H, S, a, b, Q, a1, a2, alpha in data_all:
    if (S - 1) % 2 == 0:
        k = (S-1) // 2
        print(f"  n={n} H={H:3d} √det={S:3d} = 1+2·{k:2d}")
    else:
        print(f"  n={n} H={H:3d} √det={S:3d} NOT of form 1+2k!")

print("\n--- PART 6: √det and the COMPLEMENT conflict graph ---\n")

print("Idea: if H = I(CG, 2), maybe √det = I(CG_complement, 2)?")
print("Or √det = I(some_subgraph, 2)?")
print()
print("Let's check: √det as I(G, 2) for small cases")
print("  I(empty, 2) = 1 (zero vertices)")
print("  I(K_1, 2) = 3")
print("  I(K_2, 2) = 5")
print("  I(K_3, 2) = 7")
print("  I(P_2, 2) = 5")
print("  I(C_3, 2) = 7")
print()

# Possible values of I(G, 2) for small graphs
# G=empty: 1
# G=K_1: 3
# G=K_2: 5
# G=P_2: 5
# G=K_3: 7
# G=K_1+K_1: 9
# G=K_1+K_2: 15
# G=K_2+K_2: 25
# G=P_3: 11
# G=C_3: 7
# G=K_4: 9

print("Known I(G,2) values:")
print("  1 → empty graph")
print("  3 → K₁")
print("  5 → K₂ or P₂")
print("  7 → K₃ or C₃")
print("  9 → K₁+K₁ (= 3·3) or K₄ (1+8)")
print("  11 → P₃")
print("  15 → K₁+K₂ (= 3·5)")
print("  25 → K₂+K₂ (= 5·5) or P₄ (=21)? No, P₄ gives 21")
print("  25 → K₂+K₂ (= 5·5)")
print("  49 → K₃+K₃ (= 7·7) or C₃+C₃")
print("  81 → K₁+K₁+K₁+K₁ (= 3⁴)")
print()

# Check: which √det values match known I(G,2)?
print("√det values seen:")
sdet_values = sorted(set(S for _, _, S, _, _, _, _, _, _ in data_all))
for s in sdet_values:
    # Factor
    factors = []
    rem = s
    for p in [3, 5, 7, 11, 13, 17, 19, 23]:
        while rem % p == 0:
            factors.append(p)
            rem //= p
    if rem > 1:
        factors.append(rem)
    print(f"  √det = {s:4d} = {'·'.join(str(f) for f in factors) if factors else '1'} | "
          f"I(disjoint union of K_*)? {'YES' if all(f % 2 == 1 for f in factors) else 'NO'}")

print("\n--- PART 7: √det = product of I(K_{m_i}, 2) = ∏(2m_i + 1)? ---\n")

print("HYPOTHESIS: √det = ∏(2m_i + 1) for some multiset {m_i}")
print("This would mean √det = I(K_{m_1} ⊔ K_{m_2} ⊔ ..., 2)")
print()
print("Checking: √det always a product of odd numbers:")
for s in sdet_values:
    print(f"  {s} = ", end="")
    if s == 1:
        print("1 (empty)")
    elif s == 3:
        print("3 = I(K_1, 2)")
    elif s == 5:
        print("5 = I(K_2, 2)")
    elif s == 7:
        print("7 = I(K_3, 2)")
    elif s == 9:
        print("9 = 3·3 = I(K_1⊔K_1, 2)")
    elif s == 25:
        print("25 = 5·5 = I(K_2⊔K_2, 2)")
    elif s == 49:
        print("49 = 7·7 = I(K_3⊔K_3, 2)")
    elif s == 81:
        print("81 = 3⁴ = I(K_1⊔K_1⊔K_1⊔K_1, 2)")
    else:
        print(f"{s} = ???")

print("""
AMAZING! Every √det value we see is a PRODUCT of numbers of the form (2m+1),
which are exactly the values I(K_m, 2).

So √det = I(G, 2) for some graph G that is a disjoint union of complete graphs!

CONJECTURE: √det(I+2A) = I(CG_compressed, 2)
where CG_compressed is obtained from CG by contracting each clique to a single vertex
weighted by the clique size, and I(K_m, 2) = 2m+1.

Actually, more precisely:
  H = I(CG, 2) = Σ 2^k α_k
  √det = ∏ I(K_{c_i}, 2) = ∏ (2c_i + 1)
where {c_i} is related to the CLIQUE COVER of the conflict graph!
""")

print("\n--- PART 8: Testing clique cover hypothesis ---\n")

for n in range(3, 7):
    print(f"n = {n}:")
    m = n*(n-1)//2
    hash_groups = {}
    for bits in range(1 << m):
        bts = [(bits >> i) & 1 for i in range(m)]
        A = adj_matrix(bts, n)
        h = fast_hash(A, n)
        if h not in hash_groups:
            hash_groups[h] = A

    for h, A in sorted(hash_groups.items(), key=lambda x: count_hp(x[1], n)):
        H = count_hp(A, n)
        d = det_I2A(A, n)
        S = int(round(abs(d)**0.5))

        cycles = find_directed_odd_cycles(A, n)
        if not cycles:
            print(f"  H={H:3d} √det={S:3d} | no cycles → √det should be 1: {'✓' if S==1 else '✗'}")
            continue

        cl = list(cycles)
        nc = len(cl)
        vs = [set(c) for c in cl]

        # Build conflict graph adjacency
        adj = [[False]*nc for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if vs[i] & vs[j]:
                    adj[i][j] = adj[j][i] = True

        # Find cliques (maximal or minimum cover)
        # For small graphs, just find all maximal cliques
        # A clique in CG = set of pairwise-sharing cycles = cycles sharing a common vertex
        # Actually, through-v cycles for any vertex v form a clique!

        # Count cliques from vertex perspective
        vertex_cliques = {}
        for v in range(n):
            through_v = [i for i, c in enumerate(cl) if v in vs[i]]
            if through_v:
                vertex_cliques[v] = through_v

        # The clique cover number χ̄(CG) = minimum number of cliques to cover all vertices
        # This equals the chromatic number of the complement

        # For now, just check: ∏ (2|through_v| + 1) over v... no that's wrong
        # Let me think about what the clique structure looks like

        # Each vertex v of the tournament defines a clique in CG (cycles through v)
        # CG has at most n cliques of this form

        clique_sizes = [len(tv) for tv in vertex_cliques.values()]

        # Actually, let me just check: is √det related to the CLIQUE partition of CG?
        # A clique partition = partition of cycles into groups of pairwise-adjacent cycles

        # Since n is small, let's find the clique cover number
        # (minimum number of cliques covering all cycle-vertices)

        # For now just report
        print(f"  H={H:3d} √det={S:3d} | {nc} cycles, vertex clique sizes: {sorted(clique_sizes, reverse=True)}")
    print()

print("\n--- PART 9: ALTERNATIVE — √det and the skew-adjacency Pfaffian ---\n")

print("""
Actually, for skew-symmetric matrices, det(S) = Pf(S)² (Pfaffian squared).
And I+2A = J+S where J = all-ones, S = A-A^T (skew).

For ODD n: det(J+S) = det of (n×n) matrix with J rank-1 + skew S.
  J has eigenvalue n (once) and 0 (n-1 times).
  S has eigenvalues ±iλ_k and 0 (for odd n).

  det(J+S) = n · ∏ (1 + λ_k²) ... by matrix-tree type theorem?
  Actually det(J+S) where J=all-ones might factor as n · Pf'(something).

For EVEN n: det(J+S) with S skew n×n.
  S has eigenvalues ±iλ_k.
  J+S eigenvalues mix the rank-1 and skew parts.

In either case, det(I+2A) = det(J+S) is a perfect square.
The Pfaffian connection is: for a skew matrix S:
  det(xI + S) = Pf(xI + S)² (only for even n)
  For odd n, det(S) = 0.

But J+S is NOT of the form xI+S (since J ≠ xI).
""")

# Final: compute Pfaffian-related quantities
print("Pfaffian of S = A-A^T for even n:")
for n in [4, 6]:
    print(f"\n  n = {n}:")
    m_e = n*(n-1)//2
    hash_groups = {}
    for bits in range(1 << m_e):
        bts = [(bits >> i) & 1 for i in range(m_e)]
        A = adj_matrix(bts, n)
        h = fast_hash(A, n)
        if h not in hash_groups:
            hash_groups[h] = A

    for h, A in sorted(hash_groups.items(), key=lambda x: count_hp(x[1], n)):
        H = count_hp(A, n)
        S = A - A.T
        det_S = int(round(np.linalg.det(S.astype(float))))
        pf_S = int(round(abs(det_S)**0.5))
        d = det_I2A(A, n)
        sqrt_d = int(round(abs(d)**0.5))

        # Only print first few
        if H <= 15 or H >= 43:
            print(f"    H={H:3d} det(S)={det_S:6d} Pf(S)=±{pf_S:3d} det(J+S)={d:6d} √det={sqrt_d:3d}")

print("""

SUMMARY:
  Q = (H²-det)/8 = a·b/2 where a=(H-√det)/2, b=(H+√det)/2
  Q = 0 iff H = √det iff a = 0
  √det is always a product of numbers (2m+1) for non-negative integers m
  This suggests √det = I(G', 2) for some graph G' related to CG
  The exact relationship G' ↔ CG is the next puzzle.
""")
