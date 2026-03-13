"""
PALEY vs INTERVAL at p=23 — opus-2026-03-12-S67e

Confirming the trend: Interval should win by even more.
p=23 ≡ 3 mod 4, so Paley exists.
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
    n = p
    S_list = list(S)
    dp = {}
    dp[1] = {0: 1}

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
        if step % 3 == 0 or step >= n - 3:
            print(f"    Step {step+1}/{n}: {len(dp)} masks, {count_at_step} paths", flush=True)

    full = (1 << n) - 1
    if full in dp:
        return sum(dp[full].values())
    return 0


p = 23
m = (p - 1) // 2  # = 11
S_int = list(range(1, m + 1))
S_qr = qr_set(p)

print(f"p = {p}, m = {m}")
Q_int = fourier_magnitudes(S_int, p)
Q_qr = fourier_magnitudes(S_qr, p)
print(f"prod(1+Q) Interval: {np.prod(1+Q_int):.1f}")
print(f"prod(1+Q) Paley:    {np.prod(1+Q_qr):.1f}")
print(f"Spectral ratio:     {np.prod(1+Q_qr) / np.prod(1+Q_int):.1f}")

print(f"\nComputing H(Interval)...", flush=True)
t0 = time.time()
H0_int = count_hp_from_0_fast(S_int, p)
t_int = time.time() - t0
print(f"  H_from_0 = {H0_int}, time = {t_int:.1f}s")
A_int = H0_int / np.prod(1+Q_int)
print(f"  Amplification = {A_int:.4f}")

print(f"\nComputing H(Paley)...", flush=True)
t0 = time.time()
H0_qr = count_hp_from_0_fast(S_qr, p)
t_qr = time.time() - t0
print(f"  H_from_0 = {H0_qr}, time = {t_qr:.1f}s")
A_qr = H0_qr / np.prod(1+Q_qr)
print(f"  Amplification = {A_qr:.4f}")

ratio = H0_qr / H0_int
margin = (1 - ratio) * 100 if ratio < 1 else -(ratio - 1) * 100
print(f"\nPaley/Interval = {ratio:.6f}")
print(f"Interval {'wins' if ratio < 1 else 'loses'} by {abs(margin):.4f}%")

print(f"\n{'='*60}")
print("UPDATED RACE TABLE")
print(f"{'='*60}")
print(f"""
  p  │ p mod 4 │ Best/Int │ Winner       │ Int amp
─────┼─────────┼──────────┼──────────────┼────────────
   5 │       1 │ 1.000    │ Interval     │ 0.6
   7 │       3 │ 1.080    │ Paley        │ 1.9
  11 │       3 │ 1.022    │ Paley        │ 95
  13 │       1 │ 1.000    │ Interval     │ 1,225
  17 │       1 │ 1.000    │ Interval     │ 504,227
  19 │       3 │ 0.990    │ Interval     │ 14,907,197
  23 │       3 │ {ratio:.3f}    │ {'Interval' if ratio < 1 else 'Paley':13}│ {A_int:,.0f}
""")

if ratio < 1:
    print("CONFIRMED: Interval dominates for ALL tested p ≥ 13!")
    print(f"Amplification factor ratio: Int/Paley = {A_int/A_qr:.1f}×")
    print(f"Spectral base ratio: Paley/Int = {np.prod(1+Q_qr)/np.prod(1+Q_int):.1f}×")
    print(f"Net advantage: {H0_int/H0_qr:.6f}× for Interval")
