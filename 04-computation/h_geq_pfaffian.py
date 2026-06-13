#!/usr/bin/env python3
"""
H ≥ |Pf(S)| — CAN WE PROVE THIS?
opus-2026-03-13-S67k

We know:
  H = I(CG, 2) = Σ 2^k α_k  (always ≥ 1, always odd)
  |Pf(S)| = √det(I+2A)       (always ≥ 1, always odd)
  Q = (H² - Pf²)/8 = ab/2    (observed ≥ 0)

For Q ≥ 0 we need H ≥ |Pf|, i.e., H² ≥ det(I+2A).

APPROACH 1: H = number of HP, |Pf| = |signed matching count|
  Since H counts unsigned objects and |Pf| counts signed objects:
  unsigned count ≥ |signed count| by triangle inequality... but these
  count DIFFERENT objects (paths vs matchings).

APPROACH 2: Algebraic. Both H and |Pf| are determined by A.
  H = I(CG, 2), what is Pf in terms of the conflict graph?

APPROACH 3: Inductive via deletion.
  H(T) = H(T-v) + 2μ_v
  Pf changes under deletion by cofactor expansion.

Let's explore all three.
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
    n = len(S)
    if n == 0: return 1
    if n % 2 == 1: return 0
    if n == 2: return int(round(S[0][1]))
    result = 0
    for j in range(1, n):
        if abs(S[0][j]) < 1e-10: continue
        idx = [k for k in range(n) if k != 0 and k != j]
        sub = S[np.ix_(idx, idx)]
        sign = (-1)**(j+1)
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

print("=" * 72)
print("H ≥ |Pf(S)| — TOWARD A PROOF")
print("=" * 72)

# =============================================
print("\n--- PART 1: Exhaustive verification up to n=7 ---\n")

import random
random.seed(42)

for n in range(3, 8):
    m = n*(n-1)//2
    max_bits = 1 << m

    if max_bits <= 2**21:
        hash_groups = {}
        for bits in range(max_bits):
            b = [(bits >> i) & 1 for i in range(m)]
            A = adj_matrix(b, n)
            h = fast_hash(A, n)
            if h not in hash_groups:
                hash_groups[h] = A
        sample_size = "exhaustive"
    else:
        hash_groups = {}
        for _ in range(10000):
            bits = random.randint(0, max_bits - 1)
            b = [(bits >> i) & 1 for i in range(m)]
            A = adj_matrix(b, n)
            h = fast_hash(A, n)
            if h not in hash_groups:
                hash_groups[h] = A
        sample_size = f"sampled ({len(hash_groups)} classes)"

    violations = 0
    min_ratio = float('inf')
    for h, A in hash_groups.items():
        H = count_hp(A, n)
        S = (A - A.T).astype(float)

        if n % 2 == 0:
            pf = abs(pfaffian(S))
        else:
            det_M = int(round(np.linalg.det(np.eye(n) + 2*A.astype(float))))
            pf = int(round(abs(det_M)**0.5))

        if H < pf:
            violations += 1
        ratio = H / pf if pf > 0 else float('inf')
        min_ratio = min(min_ratio, ratio)

    print(f"  n={n}: {sample_size}, {len(hash_groups)} classes, "
          f"violations={violations}, min(H/|Pf|)={min_ratio:.4f}")

# =============================================
print("\n--- PART 2: When is H/|Pf| closest to 1? ---\n")

print("Looking for tournaments where H ≈ |Pf| (ratio close to 1)")
print("These are the 'hardest' cases for the inequality.\n")

for n in range(3, 7):
    m = n*(n-1)//2
    hash_groups = {}
    for bits in range(1 << m):
        b = [(bits >> i) & 1 for i in range(m)]
        A = adj_matrix(b, n)
        h = fast_hash(A, n)
        if h not in hash_groups:
            hash_groups[h] = A

    # Find cases where H = |Pf|
    equal_cases = []
    close_cases = []
    for h, A in hash_groups.items():
        H = count_hp(A, n)
        S = (A - A.T).astype(float)
        if n % 2 == 0:
            pf = abs(pfaffian(S))
        else:
            det_M = int(round(np.linalg.det(np.eye(n) + 2*A.astype(float))))
            pf = int(round(abs(det_M)**0.5))

        if H == pf:
            equal_cases.append((H, pf))
        elif H <= 2 * pf:
            close_cases.append((H, pf, H/pf))

    print(f"n={n}: H = |Pf| in {len(equal_cases)} cases: "
          f"{sorted(set(h for h,_ in equal_cases))}")
    if close_cases:
        close_cases.sort(key=lambda x: x[2])
        for H, pf, r in close_cases[:3]:
            print(f"  Close: H={H}, |Pf|={pf}, ratio={r:.4f}")
    print()

# =============================================
print("\n--- PART 3: Algebraic approach ---\n")

print("""
For EVEN n:
  H = I(CG, 2) = Σ 2^k α_k
  |Pf(S)| = √det(S) where S = A - A^T

  We need: I(CG, 2)² ≥ det(S)

  ATTEMPT: Express det(S) in terms of the tournament structure.
  S is skew with entries ±1 off-diagonal.
  det(S) = Pf(S)² where Pf = signed matching count.

  For n=4: S is 4×4 skew with entries ±1.
    Pf(S) = S_{12}S_{34} - S_{13}S_{24} + S_{14}S_{23}
    Each S_{ij} = ±1, so Pf ∈ {-3,-1,+1,+3} (must be odd since each term is ±1).
    |Pf| ∈ {1, 3}.
    H ∈ {1, 3, 5}.
    H ≥ |Pf| holds: only need to check H=1→|Pf|≤1, H=3→|Pf|≤3.

  For n=6: S is 6×6 skew. Pf has 15 terms (number of perfect matchings of K_6).
    |Pf| ∈ {1,3,5,7,9,...}. H can be as small as 1.
    Need: when H=1 (transitive tournament), Pf(S)=±1.
    Check: transitive has S_{ij} = sign(j-i) for all i<j.
""")

# Check transitive tournament Pfaffian
for n in [4, 6, 8]:
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1  # transitive: i beats j for all i<j
    S = (A - A.T).astype(float)
    H = count_hp(A, n) if n <= 7 else "?"
    if n % 2 == 0:
        pf = pfaffian(S)
    else:
        pf = "N/A"
    print(f"  Transitive T_{n}: H={H}, Pf(S)={pf}")

# =============================================
print("\n--- PART 4: H=|Pf| characterization ---\n")

print("CONJECTURE: H = |Pf| iff the conflict graph has clique cover number ≤ 1")
print("(i.e., all cycles pairwise share a vertex, or the CG is a single clique)")
print()

for n in range(3, 7):
    m = n*(n-1)//2
    hash_groups = {}
    for bits in range(1 << m):
        b = [(bits >> i) & 1 for i in range(m)]
        A = adj_matrix(b, n)
        h = fast_hash(A, n)
        if h not in hash_groups:
            hash_groups[h] = A

    for h, A in hash_groups.items():
        H = count_hp(A, n)
        S_mat = (A - A.T).astype(float)
        if n % 2 == 0:
            pf = abs(pfaffian(S_mat))
        else:
            det_M = int(round(np.linalg.det(np.eye(n) + 2*A.astype(float))))
            pf = int(round(abs(det_M)**0.5))

        if H == pf:
            # Check: are all cycles pairwise sharing vertices?
            from itertools import combinations as combos
            cycles = set()
            for length in range(3, n + 1, 2):
                for verts in combos(range(n), length):
                    for perm in permutations(verts):
                        if all(A[perm[i]][perm[(i+1) % length]] for i in range(length)):
                            mi = perm.index(min(perm))
                            canon = perm[mi:] + perm[:mi]
                            cycles.add(canon)
            cl = list(cycles)
            vs = [set(c) for c in cl]

            # Check pairwise sharing
            all_share = True
            for i in range(len(cl)):
                for j in range(i+1, len(cl)):
                    if not (vs[i] & vs[j]):
                        all_share = False
                        break
                if not all_share:
                    break

            print(f"  n={n} H=|Pf|={H}: #cycles={len(cl)}, "
                  f"all_pairwise_share={all_share}")

# =============================================
print("\n--- PART 5: The H/|Pf| ratio as n grows ---\n")

print("For the regular tournament (all scores equal, n odd) and near-regular:")

for n in [3, 5, 7]:
    m = n*(n-1)//2
    max_bits = 1 << m

    max_ratio = 0
    max_H = 0
    max_pf = 0
    for _ in range(min(10000, max_bits)):
        if max_bits <= 2**15:
            bits = _
        else:
            bits = random.randint(0, max_bits - 1)
        b = [(bits >> i) & 1 for i in range(m)]
        A = adj_matrix(b, n)
        H = count_hp(A, n)
        det_M = int(round(np.linalg.det(np.eye(n) + 2*A.astype(float))))
        pf = int(round(abs(det_M)**0.5))
        if pf > 0:
            ratio = H / pf
            if ratio > max_ratio:
                max_ratio = ratio
                max_H = H
                max_pf = pf

    print(f"  n={n}: max(H/|Pf|) ≈ {max_ratio:.2f} (H={max_H}, |Pf|={max_pf})")

print("""
========================================================================
TOWARD A PROOF OF H ≥ |Pf(S)|
========================================================================

STRATEGY 1: Cauchy-Schwarz / AM-QM
  Both H and Pf are sums over combinatorial objects:
  H = Σ_{HP} 1  (count of Hamiltonian paths)
  Pf = Σ_{PM} ε(M) · ∏ S_{ij}  (signed matching count)
  |Pf| ≤ Σ |...| but this just gives |Pf| ≤ (n-1)!! (all matchings)
  Not useful since H can be 1 < (n-1)!!

STRATEGY 2: Use I+2A = J+S decomposition
  H² = det(I+2A) + 8Q where Q ≥ 0 (to prove)
  This is circular.

STRATEGY 3: Deletion induction
  H(T) = H(T-v) + 2μ_v ≥ H(T-v)
  |Pf(S_T)| vs |Pf(S_{T-v})|?

  For EVEN n → ODD n-1:
    √det(T) = |Pf(S)|
    √det(T-v) = |Σ (-1)^i Pf(S_{ii})|
    The sub-Pfaffians include Pf(S_{vv}) = Pf(S_{T-v})

  This is getting complicated. The inequality may require a
  more global argument.

STRATEGY 4: OCF bridge
  H = I(CG, 2) = Σ_k 2^k α_k
  |Pf| = √det(I+2A)

  For H ≥ |Pf|, we need I(CG, 2) ≥ √det(J+S).

  Maybe: I(CG, x) ≥ √det(J + x·S') for all x > 0?
  This would give the inequality as the special case x = 2.

  Or: there might be a "Pfaffian independence polynomial" P(x)
  such that Pf(S) = P(2) and I(CG, x) ≥ |P(x)| for x ≥ 0.

CONCLUSION: The inequality H ≥ |Pf(S)| is OBSERVED for all n ≤ 7
(exhaustive up to 6, sampled at 7). A proof likely requires either:
  (a) A direct algebraic identity relating I(CG,2) and Pf(S), or
  (b) An inductive argument using the deletion recurrence on both sides.

This is a STRONG open question worth pursuing.
""")
