#!/usr/bin/env python3
"""
omega_growth_rate.py — opus-2026-03-13-S70

MAJOR FINDING: For Paley tournaments, Omega_m = C(p,2) * (Q-1)^{m-1} for m<=4.
Now investigate:
1. Does this extend to general circulant tournaments?
2. What's the growth rate for interval tournaments?
3. Is there a general formula Omega_m ≈ C(n,2) * product of (Q_k - 1)?
"""

import numpy as np
from math import comb

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

def count_regular_m_paths(A, m):
    n = A.shape[0]
    count = 0
    def dfs(path_set, last, prev, depth):
        nonlocal count
        if depth == m:
            count += 1
            return
        for v in range(n):
            if v in path_set: continue
            if not A[last][v]: continue
            if depth >= 1 and not A[prev][v]: continue
            path_set.add(v)
            dfs(path_set, v, last, depth + 1)
            path_set.remove(v)
    for start in range(n):
        dfs({start}, start, -1, 0)
    return count

def fourier_magnitudes(n, S):
    """Compute Q_k = |Ŝ(k)|² for k=1,...,n-1."""
    Q = []
    for k in range(1, n):
        omega = np.exp(2j * np.pi * k / n)
        S_hat = sum(omega**s for s in S)
        Q.append(abs(S_hat)**2)
    return Q

# ============================================================
# Compare Omega growth rates across tournament types
# ============================================================
print("="*70)
print("OMEGA GROWTH RATE ANALYSIS")
print("="*70)

def analyze_tournament(name, n, S):
    A = circulant_tournament(n, S)
    Q = fourier_magnitudes(n, S)

    print(f"\n  {name} (n={n}, S={sorted(S)}):")
    print(f"    Q_k values: {[f'{q:.4f}' for q in Q]}")
    print(f"    Distinct Q: {sorted(set(round(q, 6) for q in Q))}")

    # Product formula: prod(Q_k - 1) for k where Q_k > 1
    prod_q_minus_1 = 1
    for q in Q:
        if q > 1.001:
            prod_q_minus_1 *= (q - 1)

    # Geometric mean of (Q_k - 1)
    q_minus_1_vals = [q - 1 for q in Q if q > 0.001]
    geo_mean = np.exp(np.mean(np.log([max(v, 0.001) for v in q_minus_1_vals])))

    print(f"    Geometric mean of (Q_k-1): {geo_mean:.4f}")

    profile = []
    cn2 = comb(n, 2)
    for mm in range(min(n, 10)):
        omega = count_regular_m_paths(A, mm)
        profile.append(omega)
        ratio = omega / cn2 if cn2 > 0 else 0
        # Expected from geometric model
        expected_ratio = geo_mean**(mm-1) if mm >= 1 else n/cn2
        print(f"    Omega_{mm} = {omega:10d}, ratio={ratio:10.4f}, "
              f"geo_pred={expected_ratio:10.4f}, "
              f"actual/pred={ratio/expected_ratio:.4f}" if mm >= 1 and expected_ratio > 0
              else f"    Omega_{mm} = {omega:10d}")

    # Compute successive ratios: Omega_{m+1}/Omega_m
    print(f"    Successive ratios:")
    for mm in range(1, len(profile)-1):
        if profile[mm] > 0:
            ratio = profile[mm+1] / profile[mm]
            print(f"      Omega_{mm+1}/Omega_{mm} = {ratio:.4f}")

# Paley tournaments (p ≡ 3 mod 4)
def legendre(a, p):
    if a % p == 0: return 0
    return pow(a, (p-1)//2, p)

for p in [7, 11]:
    QR = {a % p for a in range(1, p) if legendre(a, p) == 1}
    analyze_tournament(f"Paley p={p}", p, QR)

# Interval tournaments
for p in [7, 9, 11]:
    m = (p-1)//2
    S = set(range(1, m+1))
    analyze_tournament(f"Interval n={p}", p, S)

# ============================================================
# KEY TEST: Does Omega_2 = C(n,2) * mean(Q_k - 1)?
# ============================================================
print(f"\n{'='*70}")
print("Omega_2 vs sum/mean of (Q_k - 1)")
print("="*70)

for p in [5, 7, 9, 11, 13]:
    m = (p-1)//2
    S = set(range(1, m+1))
    A = circulant_tournament(p, S)
    Q = fourier_magnitudes(p, S)

    omega2 = count_regular_m_paths(A, 2)
    cn2 = comb(p, 2)
    cn3 = comb(p, 3)

    sum_Q = sum(Q)
    mean_Q = sum_Q / len(Q)
    sum_Q_minus_1 = sum(q - 1 for q in Q)
    mean_Q_minus_1 = sum_Q_minus_1 / len(Q)

    # Omega_2 = C(n,3) - t3
    t3 = cn3 - omega2

    print(f"\n  n={p}: Omega_2={omega2}, t3={t3}")
    print(f"    C(n,2)={cn2}, C(n,3)={cn3}")
    print(f"    mean(Q_k) = {mean_Q:.4f}")
    print(f"    mean(Q_k-1) = {mean_Q_minus_1:.4f}")
    print(f"    Omega_2/C(n,2) = {omega2/cn2:.4f}")
    print(f"    sum(Q_k-1) = {sum_Q_minus_1:.4f}")

    # Is Omega_2 = C(n,3) * something?
    if cn3 > 0:
        print(f"    Omega_2/C(n,3) = {omega2/cn3:.4f}")

    # Known: Omega_2 = C(n,3) - t3, and t3 = (n/24)(n²-12*sum_Q + 3(n-1))
    # Actually: t3 for circulant on Z_n with connection set S:
    # Every unordered triple {a,b,c} forms a 3-cycle iff...
    # Let me just check if t3 = n*sum_Q/8 or similar
    print(f"    t3*24/(n*(n-1)) = {t3*24/(p*(p-1)):.4f}")

# ============================================================
# Does Omega_3 = C(n,2) * mean((Q_k-1)^2)?
# ============================================================
print(f"\n{'='*70}")
print("Omega_m vs moments of (Q_k - 1)")
print("="*70)

for p in [7, 11]:
    m_half = (p-1)//2
    S_int = set(range(1, m_half+1))
    A_int = circulant_tournament(p, S_int)
    Q_int = fourier_magnitudes(p, S_int)

    if p % 4 == 3:
        QR = {a % p for a in range(1, p) if legendre(a, p) == 1}
        A_pal = circulant_tournament(p, QR)
        Q_pal = fourier_magnitudes(p, QR)
    else:
        Q_pal = None

    cn2 = comb(p, 2)

    for source, Q, A in [("Interval", Q_int, A_int)] + \
                          ([("Paley", Q_pal, A_pal)] if Q_pal else []):
        print(f"\n  {source} n={p}:")
        print(f"    Q values: {[f'{q:.2f}' for q in Q]}")

        for mm in range(1, min(p, 8)):
            omega = count_regular_m_paths(A, mm)
            # Moment: mean((Q_k - 1)^{m-1})
            moment = np.mean([(q-1)**(mm-1) for q in Q])
            predicted = cn2 * moment
            match = abs(omega - predicted) < 0.01
            print(f"    Omega_{mm}={omega:8d}, C(n,2)*E[(Q-1)^{mm-1}]={predicted:10.2f}, "
                  f"{'✓' if match else '✗'} ratio={omega/predicted:.6f}" if predicted != 0
                  else f"    Omega_{mm}={omega:8d}")

print("\nDONE.")
