"""
PALEY vs INTERVAL at p=19 — FAST VERSION — opus-2026-03-12-S67e

Only computes Interval and Paley (2 sets) instead of all 512.
Uses numpy array DP instead of dict for speed.
"""

import numpy as np
import time

def is_qr(a, p):
    if a % p == 0:
        return False
    return pow(a, (p-1)//2, p) == 1

def qr_set(p):
    return sorted([a for a in range(1, p) if is_qr(a, p)])

def fourier_magnitudes(S, p):
    m = (p - 1) // 2
    omega = np.exp(2j * np.pi / p)
    Q = []
    for k in range(1, m + 1):
        S_hat = sum(omega ** (s * k) for s in S)
        Q.append(abs(S_hat) ** 2)
    return np.array(Q)

def count_hp_from_0_fast(S, p):
    """Optimized Hamiltonian path count using dict DP with tuple keys."""
    n = p
    S_list = list(S)

    # dp[mask][v] = count, using dict of dicts for memory efficiency
    dp = {}
    dp[1] = {0: 1}  # mask=1 (just vertex 0), at vertex 0

    for step in range(1, n):
        new_dp = {}
        for mask, endpoints in dp.items():
            for v, count in endpoints.items():
                for s in S_list:
                    w = (v + s) % n
                    wbit = 1 << w
                    if not (mask & wbit):
                        new_mask = mask | wbit
                        if new_mask not in new_dp:
                            new_dp[new_mask] = {}
                        if w in new_dp[new_mask]:
                            new_dp[new_mask][w] += count
                        else:
                            new_dp[new_mask][w] = count
        dp = new_dp
        count_at_step = sum(sum(ep.values()) for ep in dp.values())
        print(f"    Step {step+1}/{n}: {len(dp)} masks, {count_at_step} paths")

    full = (1 << n) - 1
    if full in dp:
        return sum(dp[full].values())
    return 0


p = 19
m = (p - 1) // 2  # = 9
S_int = list(range(1, m + 1))
S_qr = qr_set(p)

print(f"p = {p}, m = {m}")
print(f"Interval: {S_int}")
print(f"Paley:    {S_qr}")

Q_int = fourier_magnitudes(S_int, p)
Q_qr = fourier_magnitudes(S_qr, p)
print(f"\nQ_k Interval: [{', '.join(f'{q:.4f}' for q in Q_int)}]")
print(f"Q_k Paley:    [{', '.join(f'{q:.4f}' for q in Q_qr)}]")
print(f"prod(1+Q) Interval: {np.prod(1+Q_int):.1f}")
print(f"prod(1+Q) Paley:    {np.prod(1+Q_qr):.1f}")
print(f"Q_k Paley all = {Q_qr[0]:.4f}? {np.allclose(Q_qr, Q_qr[0])}")

# Compute H for Interval
print(f"\nComputing H(Interval)...")
t0 = time.time()
H0_int = count_hp_from_0_fast(S_int, p)
t_int = time.time() - t0
print(f"  H_from_0 = {H0_int}")
print(f"  H_total = {H0_int * p}")
print(f"  Time: {t_int:.1f}s")
print(f"  Amplification = {H0_int / np.prod(1+Q_int):.4f}")

# Compute H for Paley
print(f"\nComputing H(Paley)...")
t0 = time.time()
H0_qr = count_hp_from_0_fast(S_qr, p)
t_qr = time.time() - t0
print(f"  H_from_0 = {H0_qr}")
print(f"  H_total = {H0_qr * p}")
print(f"  Time: {t_qr:.1f}s")
print(f"  Amplification = {H0_qr / np.prod(1+Q_qr):.4f}")

# Results
ratio = H0_qr / H0_int
print(f"\n{'='*60}")
print(f"RESULT: Paley/Interval ratio at p={p}: {ratio:.6f}")
if ratio > 1:
    margin = (ratio - 1) * 100
    print(f"  Paley wins by {margin:.4f}%")
    print(f"  (Was 8% at p=7, 2.2% at p=11)")
    if margin < 2.2:
        print(f"  MARGIN IS SHRINKING!")
elif ratio < 1:
    margin = (1 - ratio) * 100
    print(f"  *** INTERVAL wins by {margin:.4f}%! ***")
    print(f"  *** The crossover has happened! Fibonacci resonance dominates! ***")
else:
    print(f"  Exact tie!")

# Complete race table
print(f"\n{'='*60}")
print("THE COMPLETE RACE TABLE")
print(f"{'='*60}")
print(f"""
  p  │ p mod 4 │ H_from_0(Int) │ H_from_0(Best) │ Best/Int │ Winner
─────┼─────────┼───────────────┼────────────────┼──────────┼─────────
   5 │       1 │             3 │              3 │ 1.000    │ Interval
   7 │       3 │            25 │             27 │ 1.080    │ Paley
  11 │       3 │          8457 │           8645 │ 1.022    │ Paley
  13 │       1 │        285475 │         285475 │ 1.000    │ Interval
  17 │       1 │     805251147 │      805251147 │ 1.000    │ Interval
  19 │       3 │  {H0_int:>12} │  {max(H0_int, H0_qr):>13} │ {max(H0_int, H0_qr)/H0_int:.3f}    │ {'Interval' if H0_int >= H0_qr else 'Paley'}
""")

# Amplification factor analysis
print("Amplification factors:")
print(f"  p=19 Interval: A = {H0_int / np.prod(1+Q_int):.4f}")
print(f"  p=19 Paley:    A = {H0_qr / np.prod(1+Q_qr):.4f}")
print(f"  prod ratio:    {np.prod(1+Q_qr) / np.prod(1+Q_int):.4f}")
