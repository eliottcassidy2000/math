#!/usr/bin/env python3
"""
Verify EXACT w_0 formula at n=7 with corrected bc (alpha_2).

Previous run found: w_0 = 2*t3 - t5 + 2*t7 - bc_double - 17/4 (max err 0)
bc_double was double-counted. With alpha_2 = bc_double/2:
  w_0 = 2*t3 - t5 + 2*t7 - 2*alpha_2 - 17/4
  H = 1 + 2*(t3+t5+t7) + 4*alpha_2 (OCF verified)

Relationship: w_0 = H - 3*t5 - 6*alpha_2 - 21/4

At n=5: w_0 = H - 3*t3 = -t3 + 2*t5 + 1
At n=7: w_0 = H - 3*t5 - 6*alpha_2 - 21/4

Pattern analysis: what role does each term play?
kind-pasteur-2026-03-06-S25g
"""
from itertools import permutations, combinations
import random

def count_directed_k_cycles(A, n, k):
    count = 0
    for verts in combinations(range(n), k):
        seen = set()
        for p in permutations(verts):
            if all(A[p[i]][p[(i+1)%k]] == 1 for i in range(k)):
                start = min(range(k), key=lambda i: p[i])
                canon = tuple(p[start:] + p[:start])
                if canon not in seen:
                    seen.add(canon)
                    count += 1
    return count

def count_alpha2(A, n):
    """Count pairs of vertex-disjoint directed odd cycles (UNORDERED)."""
    # Get all directed odd cycles
    all_cycles = []
    for k in range(3, n+1, 2):
        for verts in combinations(range(n), k):
            seen = set()
            for p in permutations(verts):
                if all(A[p[i]][p[(i+1)%k]] == 1 for i in range(k)):
                    start = min(range(k), key=lambda i: p[i])
                    canon = tuple(p[start:] + p[:start])
                    if canon not in seen:
                        seen.add(canon)
                        all_cycles.append(frozenset(verts))

    # Count disjoint pairs (unordered)
    count = 0
    for i in range(len(all_cycles)):
        for j in range(i+1, len(all_cycles)):
            if all_cycles[i].isdisjoint(all_cycles[j]):
                count += 1
    return count

n = 7
random.seed(42)

print("=" * 70)
print(f"EXACT w_0 at n={n} with corrected alpha_2")
print("=" * 70)

for trial in range(15):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_directed_k_cycles(A, n, 3)
    t5 = count_directed_k_cycles(A, n, 5)
    t7 = count_directed_k_cycles(A, n, 7)
    a2 = count_alpha2(A, n)

    # H via DP
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for size in range(2, n+1):
        for S in range(1 << n):
            if bin(S).count('1') != size:
                continue
            for v in range(n):
                if not (S & (1 << v)):
                    continue
                S_prev = S ^ (1 << v)
                total = 0
                for u in range(n):
                    if (S_prev & (1 << u)) and A[u][v]:
                        total += dp.get((S_prev, u), 0)
                if total > 0:
                    dp[(S, v)] = total
    full = (1 << n) - 1
    H = sum(dp.get((full, v), 0) for v in range(n))

    # W(0)
    w0 = 0.0
    for p in permutations(range(n)):
        prod = 1.0
        for i in range(n-1):
            prod *= (A[p[i]][p[i+1]] - 0.5)
        w0 += prod

    # OCF check
    ocf_H = 1 + 2*(t3 + t5 + t7) + 4*a2

    # Predicted w0
    pred_w0 = 2*t3 - t5 + 2*t7 - 2*a2 - 17/4

    # Also: w0 = H - 3*t5 - 6*a2 - 21/4
    pred_w0_from_H = H - 3*t5 - 6*a2 - 21/4

    ocf_ok = "OK" if H == ocf_H else "FAIL"
    w0_ok = "OK" if abs(w0 - pred_w0) < 0.001 else "FAIL"

    print(f"  trial={trial:2d}: H={H:4d} ocf={ocf_ok} | t3={t3:2d} t5={t5:2d} t7={t7:2d} a2={a2:2d} | w0={w0:7.2f} pred={pred_w0:7.2f} [{w0_ok}]")

print(f"\n{'='*70}")
print("COEFFICIENT COMPARISON")
print(f"{'='*70}")
print(f"  n=5: w0 = -1*t3 + 2*t5 + 0*alpha_2 + 1")
print(f"  n=7: w0 =  2*t3 - 1*t5 + 2*t7 - 2*alpha_2 - 17/4")
print(f"")
print(f"  n=5: H - w0 = 3*t3")
print(f"  n=7: H - w0 = 3*t5 + 6*alpha_2 + 21/4")
print(f"")
print(f"  At n=5: H-w0 depends on t_3 only (simplest cycle)")
print(f"  At n=7: H-w0 depends on t_5 and alpha_2 (NOT t_3)")
print(f"  -> The 'penalty' shifts to HIGHER-ORDER cycles at larger n")
print(f"")
print(f"  W(r) expansion W(1/2)=H and W(0)=w0:")
print(f"  The difference H - w0 measures how W changes from r=0 to r=1/2")
print(f"  = contribution of the 'middle' W-polynomial terms")

print("\nDONE")
