#!/usr/bin/env python3
"""
n9_max_structure.py — Structural analysis of n=9 H=3357 maximizer

Questions:
1. Is this a Paley tournament on GF(9)?
2. How many isomorphism classes of maximizers exist?
3. What are the eigenvalues of the adjacency matrix?
4. How does c3=30 compare to the regular tournament prediction?

For GF(9): 9 = 3², so GF(9) = F_3[x]/(x²+1) or F_3[x]/(x²+x+2) etc.
The Paley tournament on GF(q) with q ≡ 3 mod 4 has connection set = QR(GF(q)).
q=9: 9 ≡ 1 mod 4?? No: 9 = 4*2+1, so 9 ≡ 1 mod 4. NOT 3 mod 4!
So Paley tournament on GF(9) does NOT exist!

Wait: q ≡ 3 mod 4. 9 mod 4 = 1. So indeed Paley doesn't work at q=9.
Then what IS the n=9 maximizer?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, random
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

def H_tournament(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full])

def random_regular_tournament(n):
    assert n % 2 == 1
    target = (n - 1) // 2
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    for _ in range(1000):
        scores = [sum(A[i]) for i in range(n)]
        if max(abs(s - target) for s in scores) == 0:
            break
        high_v = max(range(n), key=lambda v: scores[v])
        low_v = min(range(n), key=lambda v: scores[v])
        if high_v == low_v:
            break
        if A[high_v][low_v]:
            A[high_v][low_v] = 0
            A[low_v][high_v] = 1
        else:
            candidates = [w for w in range(n) if w != high_v and w != low_v
                         and A[high_v][w] and A[w][low_v]]
            if candidates:
                w = random.choice(candidates)
                A[high_v][w] = 0
                A[w][high_v] = 1
                A[w][low_v] = 0
                A[low_v][w] = 1
    return A

def is_self_converse(A, n):
    """Check if T is isomorphic to T^op (complement)."""
    Ac = [[1 - A[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
    # Simple check: same score sequence
    scores = sorted([sum(A[i]) for i in range(n)])
    scores_c = sorted([sum(Ac[i]) for i in range(n)])
    return scores == scores_c  # Necessary but not sufficient

def check_automorphism_size(A, n, max_perms=10000):
    """Estimate |Aut(T)| by sampling permutations."""
    from itertools import permutations
    count = 0
    total = 0
    if n <= 8:
        for perm in permutations(range(n)):
            total += 1
            if all(A[perm[i]][perm[j]] == A[i][j] for i in range(n) for j in range(n)):
                count += 1
    else:
        # Sample random permutations
        for _ in range(max_perms):
            perm = list(range(n))
            random.shuffle(perm)
            total += 1
            if all(A[perm[i]][perm[j]] == A[i][j] for i in range(n) for j in range(n)):
                count += 1
    return count, total

def adjacency_eigenvalues(A, n):
    """Eigenvalues of the adjacency matrix (not skew)."""
    Am = np.array(A, dtype=float)
    return sorted(np.linalg.eigvals(Am).real, reverse=True)

# ===== Find and analyze H=3357 maximizer =====
print("=" * 70)
print("n=9 MAXIMIZER STRUCTURE ANALYSIS")
print("=" * 70)

n = 9
t0 = time.time()

# Find multiple maximizers and check if they're all isomorphic
maximizers = []
eigenvalue_types = Counter()

print("\nSearching for H=3357 maximizers...")
for trial in range(50000):
    A = random_regular_tournament(n)
    scores = [sum(A[i]) for i in range(n)]
    if max(scores) - min(scores) > 0:
        continue
    H = H_tournament(A, n)
    if H == 3357:
        # Get skew eigenvalues as a signature
        S = np.array([[A[i][j] - A[j][i] for j in range(n)] for i in range(n)], dtype=float)
        evals = np.linalg.eigvals(S)
        pos_imag = sorted([e.imag for e in evals if e.imag > 0.01])
        sig = tuple(round(e, 3) for e in pos_imag)
        eigenvalue_types[sig] += 1
        if len(maximizers) < 20:
            maximizers.append([row[:] for row in A])

    if trial % 10000 == 0:
        print(f"  {trial}/50000, found {sum(eigenvalue_types.values())} maximizers")

print(f"\nSpectral types among H=3357 maximizers:")
for sig, cnt in eigenvalue_types.most_common():
    print(f"  eigenvalues {sig}: {cnt}")

# ===== Analyze first maximizer in detail =====
if maximizers:
    A = maximizers[0]
    print("\n" + "=" * 70)
    print("DETAILED ANALYSIS OF FIRST MAXIMIZER")
    print("=" * 70)

    # Adjacency matrix
    print("\nAdjacency matrix:")
    for row in A:
        print("  " + " ".join(map(str, row)))

    # Skew matrix S²
    S = np.array([[A[i][j] - A[j][i] for j in range(n)] for i in range(n)], dtype=float)
    S2 = S @ S
    print(f"\nS² diagonal: {[round(S2[i][i], 2) for i in range(n)]}")
    print(f"S² off-diagonal sample (0,1)={round(S2[0][1], 2)}, (0,2)={round(S2[0][2], 2)}")

    # Check if S² = -nI + J (conference matrix property)
    is_conference = True
    for i in range(n):
        for j in range(n):
            expected = -n if i == j else 1
            if abs(S2[i][j] - expected) > 0.01:
                is_conference = False
                break
        if not is_conference:
            break
    print(f"Conference matrix (S²=-nI+J): {is_conference}")

    # Adjacency eigenvalues
    adj_eigs = adjacency_eigenvalues(A, n)
    print(f"\nAdjacency eigenvalues: [{', '.join(f'{e:.4f}' for e in adj_eigs)}]")

    # c3, c5 counts
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                    c3 += 1
    print(f"c3 = {c3}")

    # Moon's formula for regular tournament: c3 = n(n²-1)/24 when uniform
    c3_moon = n * (n*n - 1) // 24  # For ALL regular tournaments
    print(f"Moon's c3 for regular n={n}: {c3_moon}")
    # Wait — Moon's formula: c3 = C(n,3)/4 - (n/8)sum(s_i - (n-1)/2)^2
    # For regular: sum = 0, so c3 = C(n,3)/4 = 84/4 = 21
    c3_moon_correct = n * (n-1) * (n-2) // 24  # C(n,3)/4
    print(f"Moon's c3 (corrected): C({n},3)/4 = {c3_moon_correct}")
    print(f"Actual c3 = {c3} vs Moon = {c3_moon_correct}: difference = {c3 - c3_moon_correct}")

    # Automorphism count
    print("\nChecking automorphism group size (sampling)...")
    aut_count, aut_total = check_automorphism_size(A, n, max_perms=50000)
    aut_size_est = aut_count * (362880) / aut_total  # 9! = 362880
    print(f"  Found {aut_count} automorphisms in {aut_total} random perms")
    print(f"  Estimated |Aut| = {aut_size_est:.1f}")

    # Check vertex-transitivity: does shifting by 1 preserve the tournament?
    # (Check if it's a circulant)
    is_circulant = True
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if A[i][j] != A[(i+1)%n][(j+1)%n]:
                is_circulant = False
                break
        if not is_circulant:
            break
    print(f"\nCirculant (shift-invariant): {is_circulant}")

    # Check if any circulant achieves H=3357
    print("\n" + "=" * 70)
    print("CIRCULANT TOURNAMENTS ON Z_9")
    print("=" * 70)

    # Connection set S ⊂ Z_9 with |S|=4, S∩(-S)=∅
    best_circ_H = 0
    best_circ_S = None
    from itertools import combinations
    elements = list(range(1, 9))  # 1..8
    # Pair (d, 9-d) for d=1..4
    pairs = [(1,8), (2,7), (3,6), (4,5)]
    # Choose one from each pair: 2^4 = 16 connection sets
    for mask in range(16):
        S = set()
        for bit, (a, b) in enumerate(pairs):
            if (mask >> bit) & 1:
                S.add(a)
            else:
                S.add(b)

        # Build circulant tournament
        A_circ = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j: continue
                if (j - i) % n in S:
                    A_circ[i][j] = 1

        # Verify it's a tournament
        is_tour = True
        for i in range(n):
            for j in range(i+1, n):
                if A_circ[i][j] + A_circ[j][i] != 1:
                    is_tour = False
                    break
            if not is_tour:
                break

        if not is_tour:
            continue

        H = H_tournament(A_circ, n)
        if H > best_circ_H:
            best_circ_H = H
            best_circ_S = S
        print(f"  S={sorted(S)}: H={H}")

    print(f"\nBest circulant H = {best_circ_H} (S={sorted(best_circ_S)})")
    print(f"vs global max H = 3357")

print(f"\nTotal time: {time.time()-t0:.1f}s")
