#!/usr/bin/env python3
"""
fibonacci_maximality.py -- Does the Fibonacci identity constrain H-maximization?

QUESTIONS:
1. Is prod(1+Q_k) = F_p ONLY for Interval, or for other orientations too?
2. What is prod(1+Q_k) for Paley?
3. Is there a monotonicity: larger prod(1+Q_k) => larger H?
4. What is P_m((p+1)/4) = the Interval polynomial at the Paley Q-value?

Author: kind-pasteur-2026-03-12-S59c
"""

import cmath
import math
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def fibonacci(n):
    if n <= 0:
        return 0
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def main():
    print("=" * 70)
    print("FIBONACCI IDENTITY AND H-MAXIMIZATION")
    print("=" * 70)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        omega = cmath.exp(2j * cmath.pi / p)

        pairs = [(s, p - s) for s in range(1, m + 1)]

        print(f"\n{'='*60}")
        print(f"p = {p}, m = {m}, F_p = {fibonacci(p)}")
        print(f"{'='*60}")

        # Compute all orientations
        H_dict = {}
        Q_dict = {}
        prod1Q_dict = {}

        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            A = build_adj(p, S)
            H = count_ham_paths(A, p)
            H_dict[sigma] = H

            Q_vals = []
            for k in range(1, m + 1):
                val = sum(omega ** (k * s) for s in S)
                Q_vals.append(abs(val)**2)
            Q_dict[sigma] = Q_vals

            prod_1Q = math.prod(1 + q for q in Q_vals)
            prod1Q_dict[sigma] = round(prod_1Q)

        # Part 1: Distribution of prod(1+Q_k)
        print(f"\n--- PART 1: prod(1+Q_k) DISTRIBUTION ---")
        prod1Q_groups = defaultdict(list)
        for sigma, prod1Q in prod1Q_dict.items():
            H = H_dict[sigma]
            prod1Q_groups[prod1Q].append((sigma, H))

        for prod1Q, items in sorted(prod1Q_groups.items(), key=lambda x: -x[0]):
            H_vals = sorted(set(item[1] for item in items), reverse=True)
            count = len(items)
            is_fib = (prod1Q == fibonacci(p))
            is_interval = any(s == tuple([1]*m) for s, _ in items)
            is_paley = any(s == tuple(1 if (i+1) in QR else -1 for i in range(m)) for s, _ in items)
            markers = []
            if is_interval:
                markers.append("INTERVAL")
            if is_paley:
                markers.append("PALEY")
            if is_fib:
                markers.append(f"= F_{p}")
            marker_str = " <-- " + ", ".join(markers) if markers else ""
            print(f"  prod(1+Q)={prod1Q:>10}, count={count:>3}, H={H_vals}{marker_str}")

        # Part 2: Correlation of prod(1+Q) with H
        print(f"\n--- PART 2: CORRELATION prod(1+Q) vs H ---")
        sigma_list = list(H_dict.keys())
        vals = [(prod1Q_dict[s], H_dict[s]) for s in sigma_list]
        mean_p = sum(v[0] for v in vals) / len(vals)
        mean_h = sum(v[1] for v in vals) / len(vals)
        cov = sum((v[0]-mean_p)*(v[1]-mean_h) for v in vals) / len(vals)
        std_p = (sum((v[0]-mean_p)**2 for v in vals) / len(vals)) ** 0.5
        std_h = (sum((v[1]-mean_h)**2 for v in vals) / len(vals)) ** 0.5
        if std_p > 0 and std_h > 0:
            corr = cov / (std_p * std_h)
        else:
            corr = 0
        print(f"  Rank correlation (Pearson): {corr:.6f}")
        print(f"  prod(1+Q) at max H: {prod1Q_dict[max(H_dict, key=H_dict.get)]}")
        print(f"  prod(1+Q) at min H: {prod1Q_dict[min(H_dict, key=H_dict.get)]}")

        # Part 3: Paley Q value analysis
        sigma_paley = tuple(1 if (i+1) in QR else -1 for i in range(m))
        sigma_int = tuple([1]*m)

        Q_paley = Q_dict.get(sigma_paley, [])
        Q_int = Q_dict.get(sigma_int, [])
        H_paley = H_dict.get(sigma_paley, 0)
        H_int = H_dict.get(sigma_int, 0)

        print(f"\n--- PART 3: PALEY vs INTERVAL ---")
        if Q_paley:
            print(f"  Paley: Q_k = [{', '.join(f'{q:.4f}' for q in Q_paley)}]")
            print(f"         H = {H_paley}")
            print(f"         prod(1+Q) = {prod1Q_dict.get(sigma_paley, 'N/A')}")
            print(f"         prod Q = {round(math.prod(Q_paley))}")
        print(f"  Interval: Q_k = [{', '.join(f'{q:.4f}' for q in Q_int)}]")
        print(f"            H = {H_int}")
        print(f"            prod(1+Q) = {prod1Q_dict.get(sigma_int, 'N/A')}")
        print(f"            prod Q = {round(math.prod(Q_int))}")

        # Part 4: P_m at (p+1)/4 (Paley value)
        q_paley = (p + 1) / 4
        P_at_paley = sum((-1)**j * math.comb(m+j, 2*j) * q_paley**(m-j)
                        for j in range(m + 1))
        print(f"\n--- PART 4: INTERVAL POLYNOMIAL AT PALEY Q ---")
        print(f"  P_m((p+1)/4) = P_{m}({q_paley}) = {P_at_paley:.6f}")
        print(f"  This should be zero only if Paley = Interval (impossible)")

        # Part 5: What IS the sorted Q-profile for each H class?
        print(f"\n--- PART 5: Q-PROFILE BY H VALUE ---")
        H_groups = defaultdict(list)
        for sigma in H_dict:
            H = H_dict[sigma]
            Q = tuple(sorted(Q_dict[sigma]))
            H_groups[H].append(Q)

        for H, Q_list in sorted(H_groups.items(), key=lambda x: -x[0]):
            unique_Q = set(tuple(round(q, 6) for q in Q) for Q in Q_list)
            count = len(Q_list)
            prod_1Q = round(math.prod(1 + q for q in Q_list[0]))
            prod_Q = round(math.prod(Q_list[0]))
            sum_Q = round(sum(Q_list[0]), 1)
            print(f"  H={H:>8}: count={count:>3}, #Q-profiles={len(unique_Q)}, "
                  f"sumQ={sum_Q}, prodQ={prod_Q}, prod(1+Q)={prod_1Q}")

    # Part 6: Does the product inequality work?
    print(f"\n{'='*60}")
    print("PART 6: PRODUCT INEQUALITY ANALYSIS")
    print(f"{'='*60}")
    print("Testing: Does max prod(1+Q_k) <=> max H?")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        pairs = [(s, p - s) for s in range(1, m + 1)]

        max_prod1Q = -1
        max_prod1Q_sigma = None
        max_H = -1
        max_H_sigma = None

        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            A = build_adj(p, S)
            H = count_ham_paths(A, p)

            Q_vals = []
            for k in range(1, m + 1):
                val = sum(omega ** (k * s) for s in S)
                Q_vals.append(abs(val)**2)

            prod_1Q = round(math.prod(1 + q for q in Q_vals))

            if H > max_H:
                max_H = H
                max_H_sigma = sigma

            if prod_1Q > max_prod1Q:
                max_prod1Q = prod_1Q
                max_prod1Q_sigma = sigma

        same = (max_H_sigma == max_prod1Q_sigma) or (max_H_sigma == tuple(-s for s in max_prod1Q_sigma))
        print(f"\n  p={p}: max_H_sigma = max_prod(1+Q)_sigma? {'YES' if same else 'NO'}")
        print(f"    max H = {max_H}, at prod(1+Q) = {max_prod1Q if same else '?'}")
        print(f"    max prod(1+Q) = {max_prod1Q}")

    # Part 7: prod(1+Q) for all primes (just Interval)
    print(f"\n{'='*60}")
    print("PART 7: FIBONACCI TABLE (INTERVAL)")
    print(f"{'='*60}")
    print(f"{'p':>4} {'m':>3} {'F_p':>12} {'prod(1+Q)':>12} {'match':>6}")
    print("-" * 40)
    for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        S = list(range(1, m + 1))
        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_vals.append(abs(val)**2)
        prod_1Q = math.prod(1 + q for q in Q_vals)
        fib = fibonacci(p)
        print(f"{p:>4} {m:>3} {fib:>12} {round(prod_1Q):>12} {abs(prod_1Q - fib) < 0.5:>6}")


if __name__ == '__main__':
    main()
