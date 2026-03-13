#!/usr/bin/env python3
"""
matching_cycle_duality.py -- kind-pasteur-2026-03-13-S61

The MATCHING-CYCLE DUALITY in tournaments.

Key discovery: det(I + 2A) = (Pfaffian sum)^2 where the Pfaffian sum
involves MATCHINGS. But the hidden dimension beyond lambda is about
HAMILTONIAN CYCLES (c7_dir). This suggests a deep matching-cycle duality.

CONNECTIONS TO EXPLORE:
1. For even n: Pf(S) counts signed perfect matchings.
   det(I+2A) = Pf(S)^2.
   By contrast, H counts Hamiltonian paths.
   What's the relationship between Pf(S) and H?

2. For odd n: Pf(S\v) counts signed PMs of T\{v}.
   The "Pfaffian vector" w_v = (-1)^v Pf(S\v) lies in ker(S).
   What tournament structure does this null vector encode?

3. At n=7 for the ambiguous pair:
   w = [-7, 1, 3, -7, -1, -3, 1] vs [-9, 1, -3, 3, -9, -1, -1]
   These differ in specific vertices. Which vertices?
   The vertices with large |Pf| are "matching-central" vertices.

4. The VITALI INTERPRETATION:
   Lambda captures the cycle-overlap structure (W=2,1,0 on 3-cycles).
   The Pfaffian sum captures the matching structure.
   These are DUAL perspectives on the same tournament!
   The "non-measurable" content = what matchings see but cycles don't.

5. Connection to the so(n) Lie algebra (HYP-786):
   Tournament = choice of basis for so(n).
   Pfaffian = invariant of so(n) (via the invariant polynomial ring).
   The Pfaffian sum might be the "Euler class" of the tournament.

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations, permutations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def pfaffian(M):
    n = len(M)
    if n == 0:
        return 1
    if n == 1:
        return 0
    if n == 2:
        return M[0][1]
    result = 0
    for j in range(1, n):
        if M[0][j] == 0:
            continue
        indices = [i for i in range(n) if i != 0 and i != j]
        sub = [[M[indices[a]][indices[b]] for b in range(len(indices))] for a in range(len(indices))]
        result += ((-1) ** (j + 1)) * M[0][j] * pfaffian(sub)
    return result


def pfaffian_vector(A, n):
    """Compute the Pfaffian vector w_i = (-1)^i Pf(S\i)."""
    S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]
    w = []
    for i in range(n):
        remaining = [j for j in range(n) if j != i]
        S_del = [[S[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
        pf = pfaffian(S_del)
        w.append(((-1) ** i) * pf)
    return w


def enumerate_perfect_matchings(A, n, v_del=None):
    """Enumerate all perfect matchings of T\{v_del} with their signs.
    A matching is a set of n/2 disjoint edges covering all vertices except v_del.
    The sign comes from the Pfaffian convention."""
    if v_del is not None:
        vertices = [i for i in range(n) if i != v_del]
    else:
        vertices = list(range(n))
    k = len(vertices)
    if k % 2 == 1:
        return []  # No PMs for odd vertex count

    S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]

    # Enumerate all perfect matchings by recursive pairing
    matchings = []

    def find_pms(remaining, current_matching):
        if not remaining:
            matchings.append(current_matching[:])
            return
        first = remaining[0]
        for i in range(1, len(remaining)):
            partner = remaining[i]
            current_matching.append((first, partner))
            new_remaining = remaining[1:i] + remaining[i+1:]
            find_pms(new_remaining, current_matching)
            current_matching.pop()

    find_pms(vertices, [])

    # Compute signed weight for each matching
    results = []
    for pm in matchings:
        # Weight = product of S[i][j] for each edge (i,j) in the matching
        weight = 1
        for i, j in pm:
            weight *= S[i][j]

        # Pfaffian sign: the matching defines a permutation
        # sigma = (pm[0][0], pm[0][1], pm[1][0], pm[1][1], ...)
        # The sign of this permutation relative to the identity gives the Pf sign
        perm_list = []
        for i, j in pm:
            perm_list.extend([i, j])
        # Count inversions
        inv = sum(1 for a2 in range(len(perm_list)) for b2 in range(a2+1, len(perm_list))
                  if perm_list[a2] > perm_list[b2])
        pf_sign = (-1) ** inv

        results.append((pm, weight, pf_sign, pf_sign * weight))

    return results


# ========================================================================
# ANALYSIS 1: Perfect matchings of the ambiguous pair at n=7
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: PERFECT MATCHINGS OF T\\{v} FOR AMBIGUOUS PAIR")
print("=" * 70)

n = 7
targets = [(4728, "H=109, c7=8"), (4658, "H=111, c7=9")]

for bits, label in targets:
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    w = pfaffian_vector(A, n)
    S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]

    print(f"\n  {label}:")
    print(f"    Pfaffian vector: w = {w}")
    print(f"    Pf_sum = {sum(w)}")

    # For each deleted vertex, enumerate all PMs and show the signed count
    for v_del in range(n):
        pms = enumerate_perfect_matchings(A, n, v_del)
        total_signed = sum(r[3] for r in pms)
        total_pms = len(pms)

        # Count positive and negative contributions
        pos = sum(1 for r in pms if r[3] > 0)
        neg = sum(1 for r in pms if r[3] < 0)

        # The Pfaffian value
        pf_val = pfaffian([[S[i][j] for j in range(n) if j != v_del]
                           for i in range(n) if i != v_del])

        print(f"    T\\{v_del}: {total_pms} PMs, {pos} pos, {neg} neg, "
              f"signed_sum={total_signed}, Pf={pf_val}, w_{v_del}={w[v_del]}")


# ========================================================================
# ANALYSIS 2: What do the matchings LOOK LIKE?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: MATCHING ANATOMY FOR v_del=0")
print("=" * 70)

# For the vertex with the largest |Pf| (v=0 and v=3 for H=109, v=0 and v=4 for H=111)
for bits, label in targets:
    A = binary_to_tournament(bits, n)
    S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]

    print(f"\n  {label}, deleting vertex 0:")
    pms = enumerate_perfect_matchings(A, n, v_del=0)

    for pm, weight, pf_sign, contribution in sorted(pms, key=lambda x: -x[3]):
        # Show the matching edges with their tournament direction
        edges_str = []
        for i, j in pm:
            if A[i][j]:
                edges_str.append(f"{i}->{j}")
            else:
                edges_str.append(f"{j}->{i}")
        print(f"    {edges_str}: weight={weight:+d}, pf_sign={pf_sign:+d}, contrib={contribution:+d}")


# ========================================================================
# ANALYSIS 3: Matching-cycle connection at n=5
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: MATCHING-CYCLE DUALITY AT n=5 (EXHAUSTIVE)")
print("=" * 70)

n5 = 5
by_class = defaultdict(list)

for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    H = count_ham_paths(A, n5)
    w = pfaffian_vector(A, n5)
    scores = tuple(sorted(sum(A[v]) for v in range(n5)))

    # At n=5, deleting a vertex gives a 4-vertex tournament
    # PM of K_4 = 3 perfect matchings (always)
    # The Pfaffian counts signed PMs

    by_class[(scores, H)].append((bits, w))

print(f"  Score + H classes:")
for (sc, H_val), group in sorted(by_class.items()):
    w_set = set(tuple(w) for _, w in group)
    # Show the Pfaffian vectors
    sample_w = sorted(w_set)[:3]
    ps_vals = sorted(set(sum(w) for w in w_set))
    print(f"    scores={sc}, H={H_val}: {len(group)} tours, {len(w_set)} distinct w, Pf_sum={ps_vals}")

    # Key: does the Pfaffian vector contain MORE information than H?
    # At n=5, H determines the lambda class. So the Pfaffian gives
    # information BEYOND the lambda class.


# ========================================================================
# ANALYSIS 4: The Pfaffian vector as a tournament invariant
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: PFAFFIAN VECTOR SORTED (VERTEX-INDEPENDENT)")
print("=" * 70)

# The Pfaffian vector w_i depends on vertex labeling.
# A vertex-independent invariant: the SORTED absolute values |w_i|
# and the sorted w_i themselves (as a multiset).

for bits, label in targets:
    A = binary_to_tournament(bits, n)
    w = pfaffian_vector(A, n)

    abs_w = sorted(abs(wi) for wi in w)
    print(f"  {label}:")
    print(f"    w = {w}")
    print(f"    |w| sorted = {abs_w}")
    print(f"    w sorted = {sorted(w)}")

    # The "Pfaffian degree" of each vertex = |Pf(S\v)|
    # Vertices with large Pfaffian degree are "matching-central"
    # Vertices with small Pfaffian degree are "matching-peripheral"
    for v in range(n):
        print(f"    v={v}: score={sum(A[v])}, |Pf|={abs(w[v])}")


# ========================================================================
# ANALYSIS 5: The Pfaffian vector and the score sequence
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: PFAFFIAN VECTOR vs SCORE AT n=5")
print("=" * 70)

# Does |Pf(S\v)| depend on the score of v?
pf_vs_score = defaultdict(list)

for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    w = pfaffian_vector(A, n5)
    for v in range(n5):
        score_v = sum(A[v])
        pf_vs_score[score_v].append(abs(w[v]))

print(f"  Score -> |Pf(S\\v)| distribution at n=5:")
for sc in sorted(pf_vs_score.keys()):
    vals = sorted(set(pf_vs_score[sc]))
    avg = sum(pf_vs_score[sc]) / len(pf_vs_score[sc])
    print(f"    score={sc}: |Pf| in {vals}, avg={avg:.2f}")


# ========================================================================
# ANALYSIS 6: What is the algebraic identity behind det(I+2A) = (Pf_sum)^2?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 6: THE ALGEBRAIC MEANING")
print("=" * 70)

print("""
THEOREM (det(I+2A) = (Pf_sum)^2 for tournaments):

For a tournament T on n vertices:
  A + A^T = J - I  (tournament axiom)
  S = A - A^T  (skew-adjacency matrix)
  I + 2A = J + S  (key identity)

Case 1 (n even):
  det(J + S) = det(S) * det(I + S^{-1}J) = det(S) * (1 + 1^T S^{-1} 1)
  But S^{-1} is skew-symmetric, so 1^T S^{-1} 1 = 0.
  Therefore det(J + S) = det(S) = Pf(S)^2.

Case 2 (n odd):
  S is singular (rank n-1). adj(S) has rank 1 and is PSD.
  By the matrix determinant lemma for singular S:
    det(S + uv^T) = v^T adj(S) u  (when det(S) = 0)
  With u = v = 1 (all-ones vector):
    det(J + S) = 1^T adj(S) 1

  For skew-symmetric S of odd order, adj(S) = w * w^T where
  w_i = (-1)^i * Pf(S_ii) (S_ii = S with row i, col i deleted).

  Therefore: det(J + S) = (1^T w)^2 = (sum_i (-1)^i Pf(S_ii))^2.

THE MATCHING-CYCLE DUALITY:
  det(I + 2A) involves ALL minors of S (via expansion).
  The Pfaffian sum involves DELETION Pfaffians Pf(S_ii).
  Pf(S_ii) = signed count of perfect matchings of T\\{i}.

  The Hamiltonian path count H(T) = I(Omega(T), 2) involves all
  DIRECTED ODD CYCLES.

  The formula det(I+2A) = (Pf_sum)^2 is NOT a simple function of H.
  At n=5: both H=9 and H=15 have multiple Pf_sum values.
  At n=7: the ambiguous pair has Pf_sum = -13 vs -19, while H = 109 vs 111.

  So the Pfaffian sum encodes DIFFERENT information from H.
  It's a matching-based invariant that captures the oriented structure
  in a way that's complementary to the cycle-based H.

  CONJECTURE: The complete set of tournament invariants includes both:
  - H (cycle-based, via OCF)
  - sqrt(det(I+2A)) (matching-based, via Pfaffian)
  These two invariants together may capture ALL orientation information
  that the score sequence and lambda graph miss.
""")

# Test this conjecture: does (H, |Pf_sum|) resolve all ambiguities at n=5?
print(f"  Does (H, |Pf_sum|) uniquely determine score + c_k at n=5?")
by_H_ps = defaultdict(set)
for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    H = count_ham_paths(A, n5)
    w = pfaffian_vector(A, n5)
    ps = abs(sum(w))
    scores = tuple(sorted(sum(A[v]) for v in range(n5)))
    by_H_ps[(H, ps)].add(scores)

for (H_val, ps), score_set in sorted(by_H_ps.items()):
    if len(score_set) > 1:
        print(f"    (H={H_val}, |Pf_sum|={ps}): scores={sorted(score_set)} -- MULTIPLE SCORES")
    else:
        pass  # single score, fine


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
