#!/usr/bin/env python3
"""
mod4_mechanism.py -- WHY the mod-4 dichotomy exists

At p = 3 mod 4: Paley maximizes H among all circulant tournaments.
At p = 1 mod 4: Interval maximizes H among all circulant tournaments.

This script investigates the mechanism by analyzing:
1. Spectral properties (eigenvalue magnitudes, IPR)
2. Additive energy E(S)
3. Higher-order alpha_k contributions
4. The "gap" between top and bottom H values as p grows

Key hypothesis: The maximizer is the tournament whose connection set has
the LOWEST additive energy (most spread-out Fourier spectrum).

Author: kind-pasteur-2026-03-12-S59c
"""

import cmath
from math import comb
from itertools import combinations


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths(A, n):
    dp = [[0] * n for _ in range(1 << n)]
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


def eigenvalues(p, S):
    """Compute eigenvalues lambda_r = sum_{s in S} omega^{rs}."""
    omega = cmath.exp(2j * cmath.pi / p)
    return [sum(omega ** (r * s) for s in S) for r in range(p)]


def additive_energy(p, S):
    """E(S) = #{(a,b,c,d) in S^4 : a+b = c+d mod p}."""
    S_set = set(S)
    count = 0
    for a in S_set:
        for b in S_set:
            for c in S_set:
                d = (a + b - c) % p
                if d in S_set:
                    count += 1
    return count


def spectral_ipr(lambdas, p):
    """Inverse Participation Ratio = sum |lambda_r|^4 / (sum |lambda_r|^2)^2."""
    l2 = [abs(l)**2 for l in lambdas]
    l4 = [abs(l)**4 for l in lambdas]
    return sum(l4) / sum(l2)**2


def main():
    print("=" * 70)
    print("MOD-4 DICHOTOMY MECHANISM")
    print("=" * 70)

    for p in [7, 11, 13, 17]:
        m = (p - 1) // 2

        # All circulant tournaments
        pairs = [(s, p - s) for s in range(1, m + 1)]
        all_S = []
        for bits in range(1 << m):
            S = []
            for i, (a, b) in enumerate(pairs):
                S.append(a if (bits & (1 << i)) else b)
            all_S.append(sorted(S))

        # Compute H and spectral data for each
        S_int = list(range(1, m + 1))
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        paley_exists = (p % 4 == 3) and len(S_qr) == m

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}, p mod 4 = {p % 4}")
        print(f"{'='*70}")

        data = []
        for S in all_S:
            A = build_adj(p, S)
            H = count_ham_paths(A, p)
            lam = eigenvalues(p, S)
            E = additive_energy(p, S)
            ipr = spectral_ipr(lam, p)

            # Spectral flatness: std dev of |lambda_r|^2 for r != 0
            l2 = [abs(lam[r])**2 for r in range(1, p)]
            mean_l2 = sum(l2) / len(l2)
            var_l2 = sum((x - mean_l2)**2 for x in l2) / len(l2)

            # L4/L2^2 ratio (Parseval normalized)
            L4 = sum(abs(lam[r])**4 for r in range(p))
            L2 = sum(abs(lam[r])**2 for r in range(p))

            is_int = (S == S_int)
            is_pal = paley_exists and (S == S_qr)

            data.append({
                'S': S, 'H': H, 'E': E, 'ipr': ipr,
                'var_l2': var_l2, 'L4': L4, 'L2': L2,
                'is_int': is_int, 'is_pal': is_pal
            })

        # Sort by H
        data.sort(key=lambda d: d['H'], reverse=True)

        max_H = data[0]['H']
        min_H = data[-1]['H']

        print(f"\nH range: {min_H} to {max_H} (gap = {max_H - min_H})")
        print(f"Relative gap: {(max_H-min_H)/max_H*100:.3f}%")

        print(f"\n{'S':>25} {'H':>12} {'E':>6} {'IPR':>8} {'Var(|l|^2)':>12} {'L4/L2^2':>8}")
        print("-" * 75)

        for d in data[:5]:
            label = ""
            if d['is_int']:
                label = " INT"
            if d['is_pal']:
                label = " PAL"
            print(f"{str(d['S']):>25} {d['H']:>12} {d['E']:>6} "
                  f"{d['ipr']:>8.4f} {d['var_l2']:>12.2f} "
                  f"{d['L4']/d['L2']**2:>8.4f}{label}")

        if len(data) > 10:
            print("  ...")
        for d in data[-3:]:
            label = ""
            if d['is_int']:
                label = " INT"
            if d['is_pal']:
                label = " PAL"
            print(f"{str(d['S']):>25} {d['H']:>12} {d['E']:>6} "
                  f"{d['ipr']:>8.4f} {d['var_l2']:>12.2f} "
                  f"{d['L4']/d['L2']**2:>8.4f}{label}")

        # ===== Key question: what characterizes the maximizer? =====
        print(f"\nCorrelation analysis:")

        # Rank correlation between H and various quantities
        H_vals = [d['H'] for d in data]
        E_vals = [d['E'] for d in data]
        ipr_vals = [d['ipr'] for d in data]
        var_vals = [d['var_l2'] for d in data]

        def rank_corr(xs, ys):
            """Spearman rank correlation."""
            n = len(xs)
            rx = [0]*n
            ry = [0]*n
            sx = sorted(range(n), key=lambda i: xs[i])
            sy = sorted(range(n), key=lambda i: ys[i])
            for rank, idx in enumerate(sx):
                rx[idx] = rank
            for rank, idx in enumerate(sy):
                ry[idx] = rank
            mean_r = (n-1)/2
            cov = sum((rx[i]-mean_r)*(ry[i]-mean_r) for i in range(n))
            var_x = sum((rx[i]-mean_r)**2 for i in range(n))
            var_y = sum((ry[i]-mean_r)**2 for i in range(n))
            if var_x * var_y == 0:
                return 0
            return cov / (var_x * var_y)**0.5

        print(f"  H vs E (additive energy): r = {rank_corr(H_vals, E_vals):.4f}")
        print(f"  H vs IPR:                 r = {rank_corr(H_vals, ipr_vals):.4f}")
        print(f"  H vs Var(|lambda|^2):     r = {rank_corr(H_vals, var_vals):.4f}")

        # Check: does max H correspond to max or min of each?
        max_E = data[0]['E']
        min_E = data[-1]['E']
        max_ipr = data[0]['ipr']
        min_ipr = data[-1]['ipr']
        print(f"\n  Maximizer properties:")
        print(f"    E = {max_E} (range {min_E}-{data[-1]['E']})")
        print(f"    IPR = {max_ipr:.4f} (range {data[-1]['ipr']:.4f}-{data[0]['ipr']:.4f})")

        # ===== Additive energy analysis =====
        # For Interval S = {1,...,m}: E = sum over all (a+b=c+d mod p)
        # This is related to the number of solutions to a+b=c+d in S.
        # For consecutive sets, this is a known quantity (related to the Fejér kernel).

        if paley_exists:
            pal_data = [d for d in data if d['is_pal']][0]
            int_data = [d for d in data if d['is_int']][0]
            print(f"\n  Paley vs Interval:")
            print(f"    Paley:    H={pal_data['H']}, E={pal_data['E']}, IPR={pal_data['ipr']:.4f}")
            print(f"    Interval: H={int_data['H']}, E={int_data['E']}, IPR={int_data['ipr']:.4f}")

    # ===== H VALUE DISTRIBUTION =====
    print(f"\n{'='*70}")
    print("H VALUE DISTRIBUTION BY PRIME")
    print("=" * 70)

    for p in [7, 11, 13, 17]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        H_vals = []
        for bits in range(1 << m):
            S = []
            for i, (a, b) in enumerate(pairs):
                S.append(a if (bits & (1 << i)) else b)
            A = build_adj(p, sorted(S))
            H_vals.append(count_ham_paths(A, p))

        distinct = sorted(set(H_vals), reverse=True)
        print(f"\n  p={p}: {len(distinct)} distinct H values out of {len(H_vals)} tournaments")
        for h in distinct:
            count = H_vals.count(h)
            print(f"    H={h:>15}: {count} tournaments ({count*100/len(H_vals):.1f}%)")


if __name__ == '__main__':
    main()
