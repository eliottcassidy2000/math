"""
ham_cycle_frontier.py -- kind-pasteur-2026-03-14-S83
Deep exploration of Ham cycle count C(T) and its relationship to H(T).

FROM S83:
  n=5: C(T) = 0 for H <= 5, C = 1 for H = 9,11, C = 2,3 for H = 13,15
  C(T) > 0 iff T is strongly connected (Moon's theorem)

QUESTIONS:
1. Is C(T) determined by H(T)? By (H, score)?
2. Is there a formula C(T) = f(H(T), alpha_1, alpha_2, ...)?
3. What is the maximum C(T) at each n? Does the maximizer of C = maximizer of H?
4. Does C(T) have an OCF-like formula in terms of Omega?
5. C(T)/H(T) as a function of tournament structure
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict
import sys, math

sys.stdout.reconfigure(encoding='utf-8')

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[j][i] = 1
            else: A[i][j] = 1
            idx += 1
    return A

def compute_H(A, n):
    dp = {}
    for v in range(n): dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = sum(dp.get((pm, u), 0) for u in range(n) if (pm & (1 << u)) and A[u][v])
                if t: dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def compute_C(A, n):
    """Count directed Hamiltonian cycles."""
    # Fix vertex 0, count Ham paths from 0 back to 0
    count = 0
    for perm in permutations(range(1, n)):
        path = (0,) + perm
        valid = all(A[path[i]][path[i+1]] for i in range(n-1))
        if valid and A[path[-1]][path[0]]:
            count += 1
    return count

def main():
    print("=" * 70)
    print("HAMILTONIAN CYCLE FRONTIER")
    print("kind-pasteur-2026-03-14-S83")
    print("=" * 70)

    # ========================================
    # PART 1: C(T) vs H(T) at n=3,4,5,6
    # ========================================
    for n in [3, 4, 5, 6]:
        m = n*(n-1)//2
        N = 2**m

        print(f"\n{'='*70}")
        print(f"n = {n}")
        print(f"{'='*70}")

        HC_data = []
        count = 0
        for bits in range(N):
            count += 1
            if n >= 6 and count > 10000: break
            A = bits_to_adj(bits, n)
            H = compute_H(A, n)
            C = compute_C(A.tolist(), n)
            scores = tuple(sorted(A.sum(axis=1).astype(int)))
            c3 = int(np.trace(A @ A @ A)) // 3
            HC_data.append({'H': H, 'C': C, 'scores': scores, 'c3': c3})

        # Joint distribution
        HC_dist = Counter((d['H'], d['C']) for d in HC_data)

        print(f"\n  (H, C) joint distribution:")
        print(f"  {'H':>4} {'C':>4} {'Count':>6} {'C/H':>8} {'C*n/H':>8}")
        for (h, c) in sorted(HC_dist.keys()):
            cnt = HC_dist[(h, c)]
            ratio = c / h if h > 0 else 0
            cn_ratio = c * n / h if h > 0 else 0
            print(f"  {h:4d} {c:4d} {cnt:6d} {ratio:8.4f} {cn_ratio:8.4f}")

        # Is C determined by H?
        C_by_H = defaultdict(set)
        for d in HC_data:
            C_by_H[d['H']].add(d['C'])

        C_det_by_H = all(len(v) == 1 for v in C_by_H.values())
        print(f"\n  C determined by H? {C_det_by_H}")
        if not C_det_by_H:
            for H, C_set in sorted(C_by_H.items()):
                if len(C_set) > 1:
                    print(f"    H={H}: C in {sorted(C_set)}")

        # Maximum C
        max_C = max(d['C'] for d in HC_data)
        max_C_H = [d['H'] for d in HC_data if d['C'] == max_C]
        max_H = max(d['H'] for d in HC_data)
        print(f"\n  max C = {max_C}, at H = {sorted(set(max_C_H))}")
        print(f"  max H = {max_H}")
        print(f"  C-maximizer = H-maximizer? {set(max_C_H) == {max_H}}")

        # C*n/H ratio analysis
        # For a REGULAR tournament on odd n: each vertex is equivalent.
        # C(T) = H(T) * (#cycles through vertex 0) / n
        # Wait: C = #{Ham cycles} = H_cycles. And each cycle has n edges.
        # Each Ham PATH of length n-1 can be closed to a cycle if the last
        # vertex connects back to the first.
        # C = sum over (start, end) pairs: M[start][end] * A[end][start]
        # where M is the transfer matrix.
        # C = sum_{v} M[v][v'] * A[v'][v] for all v' != v... no
        # Actually: C = (1/n) * sum over all Ham paths P: A[P[-1]][P[0]]
        # No wait, C is not 1/n of anything for general tournaments.

        # Relationship: for vertex-transitive (regular) tournaments:
        # C = (H / n) * (average number of "closing arcs")
        # Not obviously simple.

        # Check: is C*n / H close to an integer?
        if n <= 5:
            print(f"\n  C*n/H for each tournament:")
            for d in HC_data:
                if d['H'] > 0 and d['C'] > 0:
                    ratio = d['C'] * n / d['H']
                    if abs(ratio - round(ratio)) < 0.01:
                        pass  # print only non-trivial
                    else:
                        pass  # non-integer ratio exists

            # Aggregate
            ratios = set()
            for d in HC_data:
                if d['H'] > 0 and d['C'] > 0:
                    ratios.add(round(d['C'] * n / d['H'], 4))
            print(f"    Distinct C*n/H values: {sorted(ratios)[:10]}")

    # ========================================
    # PART 2: The C/H ratio and strong connectivity
    # ========================================
    print(f"\n{'='*70}")
    print("PART 2: STRONG CONNECTIVITY THRESHOLD")
    print("  C > 0 iff T is strongly connected (Moon)")
    print("  At what H does strong connectivity first appear?")
    print(f"{'='*70}")

    for n in [4, 5, 6]:
        m = n*(n-1)//2
        N = 2**m

        threshold_H = None
        count = 0
        for bits in range(N):
            count += 1
            if n >= 6 and count > 20000: break
            A = bits_to_adj(bits, n)
            H = compute_H(A, n)
            C = compute_C(A.tolist(), n)
            if C > 0:
                if threshold_H is None or H < threshold_H:
                    threshold_H = H

        print(f"  n={n}: smallest H with C>0: H={threshold_H}")

    # ========================================
    # PART 3: Does C have an OCF-like formula?
    # ========================================
    print(f"\n{'='*70}")
    print("PART 3: OCF-LIKE FORMULA FOR C?")
    print("  H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...")
    print("  Is there C = f(alpha_1, alpha_2, ...)?")
    print(f"{'='*70}")

    n = 5
    m = n*(n-1)//2
    N = 2**m

    # For n=5: alpha_2 = 0, so alpha_1 = (H-1)/2 determines H.
    # Does alpha_1 determine C?
    alpha_C = defaultdict(set)
    for bits in range(N):
        A = bits_to_adj(bits, n)
        H = compute_H(A, n)
        C = compute_C(A.tolist(), n)
        alpha_1 = (H - 1) // 2
        alpha_C[alpha_1].add(C)

    print(f"\n  n=5: alpha_1 -> C:")
    for a1 in sorted(alpha_C.keys()):
        print(f"    alpha_1={a1}: C in {sorted(alpha_C[a1])}")

    # Check: C = alpha_1 - (n-1)? or C = max(0, alpha_1 - threshold)?
    # From data: alpha_1 = 0->C=0, 1->C=0, 2->C=0, 3->C={0,1},
    #            4->C=1, 5->C=1, 6->C=2, 7->C={2,3}
    # Threshold at alpha_1 = 3 (first C > 0)
    # For alpha_1 >= 4: C = alpha_1 - 3? No: alpha_1=4->C=1, 5->C=1 (not 2)

    # Actually C is NOT determined by alpha_1 alone (alpha_1=3 gives C in {0,1})
    # and alpha_1=7 gives C in {2,3}

    # ========================================
    # PART 4: The H - n*C relationship
    # ========================================
    print(f"\n{'='*70}")
    print("PART 4: RELATIONSHIP H = n*C + R")
    print("  For each tournament, H = n*C + R where R = 'residual'")
    print("  R = number of Ham paths that DON'T close to a cycle")
    print(f"{'='*70}")

    for n in [5]:
        m = n*(n-1)//2
        N = 2**m

        R_by_H = defaultdict(list)
        for bits in range(N):
            A = bits_to_adj(bits, n)
            H = compute_H(A, n)
            C = compute_C(A.tolist(), n)
            R = H - n * C
            R_by_H[H].append(R)

        print(f"\n  n={n}: H = n*C + R decomposition:")
        for H in sorted(R_by_H.keys()):
            Rs = sorted(set(R_by_H[H]))
            Cs = sorted(set((H - r) // n for r in Rs))
            print(f"    H={H:3d}: C in {Cs}, R in {Rs}")

        # Is R always positive? (i.e., not all paths close to cycles)
        all_R_pos = all(r >= 0 for rs in R_by_H.values() for r in rs)
        print(f"\n  R >= 0 always? {all_R_pos}")

        # What is R for the maximizer?
        max_H = max(R_by_H.keys())
        R_max = sorted(set(R_by_H[max_H]))
        print(f"  H=max={max_H}: R in {R_max}")
        for r in R_max:
            c = (max_H - r) // n
            print(f"    C={c}, R={r}, check: {n}*{c}+{r} = {n*c+r}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
