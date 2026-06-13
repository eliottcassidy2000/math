#!/usr/bin/env python3
"""
Verify opus's c_0 formula at n=7 and connect to our w_0 computation.

Opus found: c_0 = 2*t_3 - t_5 + 2*t_7 - 2*bc + 253/4
But our w_0 = W(0) should equal c_0 (the constant term of W(r) expansion).

At n=5: w_0 = H - 3*t_3 = 1 - t_3 + 2*t_5  (since H = 1 + 2*(t_3+t_5))
And also: c_0 formula should match.

At n=7: The OCF gives H = 1 + 2*(t_3+t_5+t_7) + 4*bc
So: w_0 = H - something = 1 + 2*t_3 + 2*t_5 + 2*t_7 + 4*bc - something

If opus's c_0 = 2*t_3 - t_5 + 2*t_7 - 2*bc + 253/4 is correct,
then w_0 - H = c_0 - H = (2*t_3 - t_5 + 2*t_7 - 2*bc + 253/4) - (1 + 2*t_3 + 2*t_5 + 2*t_7 + 4*bc)
             = -3*t_5 - 6*bc + 249/4

Let me verify this numerically.
kind-pasteur-2026-03-06-S25g
"""
from itertools import permutations, combinations
from math import factorial
import random
import numpy as np

def count_directed_k_cycles(A, n, k):
    """Count directed k-cycles (each cycle counted once, rotation-canonical)."""
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

def count_bc(A, n):
    """Count disjoint 3-cycle pairs (bc = binomial convolution of t_3)."""
    count = 0
    # All pairs of disjoint 3-element subsets
    for v1 in combinations(range(n), 3):
        remaining = [x for x in range(n) if x not in v1]
        for v2 in combinations(remaining, 3):
            # Count directed 3-cycles on v1
            c1 = 0
            for p in permutations(v1):
                if all(A[p[i]][p[(i+1)%3]] == 1 for i in range(3)):
                    c1 += 1
            c1 //= 3  # rotation canonical
            # Count directed 3-cycles on v2
            c2 = 0
            for p in permutations(v2):
                if all(A[p[i]][p[(i+1)%3]] == 1 for i in range(3)):
                    c2 += 1
            c2 //= 3
            count += c1 * c2
    return count

n = 7
random.seed(42)

print("=" * 70)
print(f"w_0 exact formula at n={n}")
print("Opus claims: c_0 = 2*t_3 - t_5 + 2*t_7 - 2*bc + 253/4")
print("=" * 70)

data = []
for trial in range(20):
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
    bc = count_bc(A, n)

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

    # Opus formula
    opus_c0 = 2*t3 - t5 + 2*t7 - 2*bc + 253/4

    # OCF check
    ocf_H = 1 + 2*(t3 + t5 + t7) + 4*bc

    data.append({
        'H': H, 't3': t3, 't5': t5, 't7': t7, 'bc': bc,
        'w0': w0, 'opus_c0': opus_c0, 'ocf_H': ocf_H
    })

    err = abs(w0 - opus_c0)
    ocf_err = abs(H - ocf_H)
    if trial < 8:
        print(f"  trial={trial}: H={H}, ocf_H={ocf_H}, t3={t3}, t5={t5}, t7={t7}, bc={bc}")
        print(f"    w0={w0:.4f}, opus_c0={opus_c0:.4f}, err={err:.4f}, ocf_err={ocf_err}")

max_err = max(abs(d['w0'] - d['opus_c0']) for d in data)
max_ocf_err = max(abs(d['H'] - d['ocf_H']) for d in data)
print(f"\nMax |w0 - opus_c0| = {max_err:.6f}")
print(f"Max |H - ocf_H| = {max_ocf_err:.6f}")

if max_err < 0.01:
    print("\nOPUS FORMULA VERIFIED: w0 = 2*t3 - t5 + 2*t7 - 2*bc + 253/4")
else:
    print("\nOpus formula does NOT match. Fitting w0 = a*t3 + b*t5 + c*t7 + d*bc + e...")
    y = np.array([d['w0'] for d in data])
    X = np.array([[d['t3'], d['t5'], d['t7'], d['bc'], 1] for d in data])
    c = np.linalg.lstsq(X, y, rcond=None)[0]
    err_fit = max(abs(y - X @ c))
    print(f"  w0 = {c[0]:.4f}*t3 + {c[1]:.4f}*t5 + {c[2]:.4f}*t7 + {c[3]:.4f}*bc + {c[4]:.4f}")
    print(f"  Max error: {err_fit:.6f}")

    # Try with H instead
    X2 = np.array([[d['H'], d['t3'], d['t5'], d['t7'], d['bc'], 1] for d in data])
    c2 = np.linalg.lstsq(X2, y, rcond=None)[0]
    err_fit2 = max(abs(y - X2 @ c2))
    print(f"\n  w0 = {c2[0]:.4f}*H + {c2[1]:.4f}*t3 + {c2[2]:.4f}*t5 + {c2[3]:.4f}*t7 + {c2[4]:.4f}*bc + {c2[5]:.4f}")
    print(f"  Max error: {err_fit2:.6f}")

# Pattern comparison with n=5
print(f"\n{'='*70}")
print("PATTERN COMPARISON")
print(f"{'='*70}")
print(f"  n=5: w0 = -t3 + 2*t5 + 1")
print(f"        H  = 1 + 2*(t3 + t5)")
print(f"        w0 = H - 3*t3")
print(f"  n=7: w0 = a*t3 + b*t5 + c*t7 + d*bc + const")
print(f"        H  = 1 + 2*(t3 + t5 + t7) + 4*bc")
print(f"        w0 = H - ?*t3 - ?*t5 - ?*bc")

print("\nDONE")
