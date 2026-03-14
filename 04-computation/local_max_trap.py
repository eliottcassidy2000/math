"""
local_max_trap.py -- kind-pasteur-2026-03-14-S71
Investigate the H=37 local maximum at n=6 — why is it a trap?

From h_landscape_deep.py: at n=6, gradient ascent gets stuck at H=37
in ~10% of random starts. H=37 is NOT the global max (H=45).

Questions:
1. What tournaments have H=37 and are local maxima?
2. What is their structure (score, cycle counts, beta)?
3. How many of the 1200 local maxima are at H=37 vs H=45?
4. What is the "barrier" — min H on any path from H=37-max to H=45-max?
5. Can simulated annealing escape the H=37 trap?
6. Connection to the (2,3) structure: is 37 a "tournament number"?
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict
import sys, random

sys.stdout.reconfigure(encoding='utf-8')

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

def compute_H_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def main():
    print("=" * 70)
    print("LOCAL MAXIMUM TRAP ANALYSIS — H=37 at n=6")
    print("kind-pasteur-2026-03-14-S71")
    print("=" * 70)

    n = 6
    m = n * (n-1) // 2
    total = 2**m

    # ========================================
    # PART 1: Find and characterize all local maxima
    # ========================================
    print("\n" + "=" * 70)
    print("PART 1: ALL LOCAL MAXIMA AT n=6")
    print("=" * 70)

    local_maxima = []  # (bits, H)
    H_cache = {}

    for bits in range(total):
        A = bits_to_adj(bits, n)
        H = compute_H_dp(A, n)
        H_cache[bits] = H

    for bits in range(total):
        H = H_cache[bits]
        is_max = True
        for arc_bit in range(m):
            nbr = bits ^ (1 << arc_bit)
            if H_cache[nbr] > H:
                is_max = False
                break
        if is_max:
            local_maxima.append((bits, H))

    # Count by H value
    max_H_count = Counter(h for _, h in local_maxima)
    print(f"  Total local maxima: {len(local_maxima)}")
    print(f"  By H value: {dict(sorted(max_H_count.items()))}")

    # ========================================
    # PART 2: Structure of H=37 local maxima
    # ========================================
    print("\n" + "=" * 70)
    print("PART 2: STRUCTURE OF H=37 LOCAL MAXIMA")
    print("=" * 70)

    h37_maxima = [(b, h) for b, h in local_maxima if h == 37]
    h45_maxima = [(b, h) for b, h in local_maxima if h == 45]

    print(f"  H=37 local maxima: {len(h37_maxima)}")
    print(f"  H=45 local maxima: {len(h45_maxima)}")

    # Score sequences
    h37_scores = Counter()
    h45_scores = Counter()

    for bits, H in h37_maxima:
        A = bits_to_adj(bits, n)
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        h37_scores[scores] += 1

    for bits, H in h45_maxima:
        A = bits_to_adj(bits, n)
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        h45_scores[scores] += 1

    print(f"\n  H=37 score sequences:")
    for s, c in sorted(h37_scores.items()):
        print(f"    {s}: {c} tournaments")

    print(f"\n  H=45 score sequences:")
    for s, c in sorted(h45_scores.items()):
        print(f"    {s}: {c} tournaments")

    # Cycle counts
    print(f"\n  H=37 cycle structure:")
    for bits, H in h37_maxima[:5]:
        A = bits_to_adj(bits, n)
        c3 = int(np.trace(A @ A @ A)) // 3
        c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        score_var = np.var([sum(A[i]) for i in range(n)])
        print(f"    bits={bits}: c3={c3}, c5={c5}, scores={scores}, var={score_var:.2f}")

    print(f"\n  H=45 cycle structure:")
    for bits, H in h45_maxima[:5]:
        A = bits_to_adj(bits, n)
        c3 = int(np.trace(A @ A @ A)) // 3
        c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        score_var = np.var([sum(A[i]) for i in range(n)])
        print(f"    bits={bits}: c3={c3}, c5={c5}, scores={scores}, var={score_var:.2f}")

    # ========================================
    # PART 3: Barrier between H=37 and H=45
    # ========================================
    print("\n" + "=" * 70)
    print("PART 3: BARRIER BETWEEN H=37 AND H=45 LOCAL MAXIMA")
    print("  What is the minimum H on any path from a H=37 max to a H=45 max?")
    print("=" * 70)

    # BFS from a H=37 max, recording min H along path
    if h37_maxima:
        start = h37_maxima[0][0]
        # Find shortest path to any H=45 max
        from collections import deque
        visited = {start: (0, H_cache[start])}  # bits -> (dist, min_H_on_path)
        queue = deque([(start, H_cache[start])])  # (bits, min_H_so_far)
        found = False

        while queue and not found:
            curr, min_H = queue.popleft()
            dist = visited[curr][0]

            if dist > 6:  # limit search depth
                continue

            for arc_bit in range(m):
                nbr = bits ^ (1 << arc_bit)
                # Oops, should be curr not bits
                nbr = curr ^ (1 << arc_bit)
                new_min = min(min_H, H_cache[nbr])

                if nbr not in visited or visited[nbr][1] < new_min:
                    visited[nbr] = (dist + 1, new_min)
                    queue.append((nbr, new_min))

                    if H_cache[nbr] == 45:
                        # Check if nbr is a local max
                        is_max = True
                        for ab in range(m):
                            if H_cache[nbr ^ (1 << ab)] > H_cache[nbr]:
                                is_max = False
                                break
                        if is_max:
                            print(f"  Found path from H=37 max to H=45 max!")
                            print(f"    Distance: {dist + 1} flips")
                            print(f"    Min H on path: {new_min}")
                            found = True
                            break

        if not found:
            print(f"  No path found within 6 flips.")

    # ========================================
    # PART 4: What is special about H=37?
    # ========================================
    print("\n" + "=" * 70)
    print("PART 4: NUMBER-THEORETIC PROPERTIES OF H=37")
    print("=" * 70)

    print(f"  37 is prime")
    print(f"  37 = 36 + 1 = 6^2 + 1")
    print(f"  37 mod 4 = 1")
    print(f"  37 = C(6,2) + C(6,1) + 1 = 15 + 6 + 1 + ... no, 37 doesn't factor nicely into C(6,k)")
    print(f"  37 in OCF: H=37 means 1 + 2*alpha_1 + 4*alpha_2 = 37")
    print(f"  => alpha_1 + 2*alpha_2 = 18")
    print(f"  Possible: (18,0), (16,1), (14,2), ..., (0,9)")

    # Compute actual (alpha_1, alpha_2) for H=37 tournaments
    print(f"\n  Actual (alpha_1, alpha_2) for H=37:")
    alpha_count = Counter()
    for bits in range(total):
        if H_cache[bits] == 37:
            A = bits_to_adj(bits, n)
            c3 = int(np.trace(A @ A @ A)) // 3
            c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
            # alpha_1 = c3 + c5 (at n=6, no 7-cycles)
            alpha_1 = c3 + c5
            # alpha_2 = (H-1-2*alpha_1) / 4  (from OCF)
            alpha_2 = (37 - 1 - 2*alpha_1) // 4
            alpha_count[(alpha_1, alpha_2)] += 1

    for (a1, a2), cnt in sorted(alpha_count.items()):
        print(f"    (alpha_1={a1}, alpha_2={a2}): {cnt} tournaments")

    # Same for H=45
    print(f"\n  Actual (alpha_1, alpha_2) for H=45:")
    alpha_count45 = Counter()
    for bits in range(total):
        if H_cache[bits] == 45:
            A = bits_to_adj(bits, n)
            c3 = int(np.trace(A @ A @ A)) // 3
            c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
            alpha_1 = c3 + c5
            alpha_2 = (45 - 1 - 2*alpha_1) // 4
            alpha_count45[(alpha_1, alpha_2)] += 1

    for (a1, a2), cnt in sorted(alpha_count45.items()):
        print(f"    (alpha_1={a1}, alpha_2={a2}): {cnt} tournaments")

    # ========================================
    # PART 5: H=37 vs H=45 — what's the structural difference?
    # ========================================
    print("\n" + "=" * 70)
    print("PART 5: STRUCTURAL COMPARISON H=37 vs H=45 LOCAL MAXIMA")
    print("=" * 70)

    # Determinant
    det37 = set()
    det45 = set()
    for bits, H in h37_maxima:
        A = bits_to_adj(bits, n)
        det37.add(round(np.linalg.det(A.astype(float))))
    for bits, H in h45_maxima:
        A = bits_to_adj(bits, n)
        det45.add(round(np.linalg.det(A.astype(float))))

    print(f"  det(A) for H=37 local maxima: {sorted(det37)}")
    print(f"  det(A) for H=45 local maxima: {sorted(det45)}")

    # Self-complementary?
    sc37 = 0
    sc45 = 0
    for bits, H in h37_maxima:
        A = bits_to_adj(bits, n)
        A_op = 1 - A - np.eye(n, dtype=int)
        # Check if A is isomorphic to A_op (need permutation check)
        # Simplified: just check if same score sequence
        s1 = tuple(sorted([sum(A[i]) for i in range(n)]))
        s2 = tuple(sorted([sum(A_op[i]) for i in range(n)]))
        if s1 == s2:
            sc37 += 1
    for bits, H in h45_maxima:
        A = bits_to_adj(bits, n)
        A_op = 1 - A - np.eye(n, dtype=int)
        s1 = tuple(sorted([sum(A[i]) for i in range(n)]))
        s2 = tuple(sorted([sum(A_op[i]) for i in range(n)]))
        if s1 == s2:
            sc45 += 1

    print(f"  SC-score compatible: H=37: {sc37}/{len(h37_maxima)}, H=45: {sc45}/{len(h45_maxima)}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)

if __name__ == '__main__':
    main()
