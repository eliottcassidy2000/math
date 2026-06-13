"""
PALEY vs INTERVAL RACE at p=19 — opus-2026-03-12-S67e

p=19 ≡ 3 mod 4, so the PALEY tournament exists.
This is the CRITICAL test: does Paley still beat Interval?
The margin was 8% at p=7, 2.2% at p=11 — is it shrinking to zero?
"""

import numpy as np
from itertools import combinations
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

def count_hp_from_0(S, p):
    """Count Hamiltonian paths starting from vertex 0 using bitmask DP."""
    n = p
    adj = {}
    for v in range(n):
        adj[v] = [(v + s) % n for s in S]

    full = (1 << n) - 1
    dp = {(1, 0): 1}

    for step in range(1, n):
        new_dp = {}
        for (mask, v), count in dp.items():
            for w in adj[v]:
                if not (mask & (1 << w)):
                    key = (mask | (1 << w), w)
                    if key in new_dp:
                        new_dp[key] += count
                    else:
                        new_dp[key] = count
        dp = new_dp

    return sum(dp.get((full, v), 0) for v in range(n))


p = 19
m = (p - 1) // 2  # = 9
print(f"p = {p}, m = {m}")

# Interval
S_int = list(range(1, m + 1))
print(f"\nInterval: S = {S_int}")
t0 = time.time()
H0_int = count_hp_from_0(S_int, p)
print(f"  H_from_0 = {H0_int}, time = {time.time()-t0:.1f}s")
Q_int = fourier_magnitudes(S_int, p)
prod_int = np.prod(1 + Q_int)
print(f"  prod(1+Q) = {prod_int:.1f}")
print(f"  Amplification = {H0_int / prod_int:.4f}")

# Paley (QR set)
S_qr = qr_set(p)
print(f"\nPaley: S = {S_qr}")
t0 = time.time()
H0_qr = count_hp_from_0(S_qr, p)
print(f"  H_from_0 = {H0_qr}, time = {time.time()-t0:.1f}s")
Q_qr = fourier_magnitudes(S_qr, p)
prod_qr = np.prod(1 + Q_qr)
print(f"  prod(1+Q) = {prod_qr:.1f}")
print(f"  Q_k flat? All = {Q_qr[0]:.4f}? {np.allclose(Q_qr, Q_qr[0])}")
print(f"  Amplification = {H0_qr / prod_qr:.4f}")

print(f"\nPaley/Interval ratio: {H0_qr / H0_int:.6f}")

# Enumerate ALL valid sets
print(f"\nEnumerating all {2**m} valid connection sets...")
valid_sets = []
elements = list(range(1, p))
for S in combinations(elements, m):
    S = list(S)
    ok = True
    for s in S:
        if (p - s) in S:
            ok = False
            break
    if not ok:
        continue
    valid_sets.append(S)

print(f"Total valid: {len(valid_sets)}")

results = []
t0 = time.time()
for i, S in enumerate(valid_sets):
    if i % 50 == 0:
        elapsed = time.time() - t0
        if i > 0:
            eta = elapsed / i * (len(valid_sets) - i)
            print(f"  {i}/{len(valid_sets)}, elapsed {elapsed:.1f}s, ETA {eta:.1f}s")
        else:
            print(f"  Starting...")

    H0 = count_hp_from_0(S, p)
    Q = fourier_magnitudes(S, p)
    prod_Q = np.prod(1 + Q)
    var_Q = np.var(Q)
    results.append((S, H0, prod_Q, var_Q))

elapsed = time.time() - t0
print(f"  Done in {elapsed:.1f}s")

# Sort by H
results.sort(key=lambda x: x[1], reverse=True)

print(f"\nTop 10 by H_from_0:")
print(f"  {'Rank':>4} {'H_from_0':>14} {'prod(1+Q)':>14} {'var(Q)':>10} {'amp':>12}")
for i, (S, H0, pq, vq) in enumerate(results[:10]):
    amp = H0 / pq
    is_int = " [INT]" if S == S_int else ""
    is_qr = " [QR]" if set(S) == set(S_qr) else ""
    print(f"  {i+1:>4} {H0:>14} {pq:>14.1f} {vq:>10.4f} {amp:>12.4f}{is_int}{is_qr}")

# Find Paley and Interval ranks
int_rank = [i for i, r in enumerate(results) if r[0] == S_int][0] + 1
qr_results = [i for i, r in enumerate(results) if set(r[0]) == set(S_qr)]
qr_rank = qr_results[0] + 1 if qr_results else "N/A"

print(f"\nInterval rank: {int_rank}/{len(results)}")
print(f"Paley rank: {qr_rank}/{len(results)}")

# H distribution
from collections import Counter
H_counts = Counter(r[1] for r in results)
print(f"\nH distribution ({len(H_counts)} distinct values):")
for h, count in sorted(H_counts.items(), reverse=True):
    rep = [r for r in results if r[1] == h][0]
    is_int = " <-- INTERVAL" if rep[0] == S_int else ""
    is_qr = " <-- PALEY" if set(rep[0]) == set(S_qr) else ""
    print(f"  H_from_0 = {h}: {count} sets, prod={rep[2]:.0f}, var={rep[3]:.2f}{is_int}{is_qr}")

# Correlations
H_vals = [r[1] for r in results]
prod_vals = [r[2] for r in results]
var_vals = [r[3] for r in results]
print(f"\nCorrelations:")
print(f"  Corr(H, prod(1+Q)) = {np.corrcoef(H_vals, prod_vals)[0,1]:.6f}")
print(f"  Corr(H, var(Q))    = {np.corrcoef(H_vals, var_vals)[0,1]:.6f}")

# The big race table
print(f"\n{'='*72}")
print("THE COMPLETE RACE TABLE")
print(f"{'='*72}")
print(f"""
  p  │ p mod 4 │ H_from_0(Int)│ H_from_0(Best)│ Best/Int │ Winner
─────┼─────────┼──────────────┼───────────────┼──────────┼─────────
   5 │       1 │            3 │             3 │ 1.000    │ Interval
   7 │       3 │           25 │            27 │ 1.080    │ Paley
  11 │       3 │         8457 │          8645 │ 1.022    │ Paley
  13 │       1 │       285475 │        285475 │ 1.000    │ Interval
  17 │       1 │    805251147 │     805251147 │ 1.000    │ Interval
  19 │       3 │   {H0_int:>10} │  {results[0][1]:>12} │ {results[0][1]/H0_int:.3f}    │ {'INTERVAL!' if results[0][1] == H0_int else 'Paley' if set(results[0][0]) == set(S_qr) else 'OTHER'}
""")

if results[0][1] == H0_int:
    print("*** INTERVAL WINS AT p=19! Paley advantage has VANISHED! ***")
    print("*** The Fibonacci resonance cascade DOMINATES for all tested p ≥ 13! ***")
elif H0_qr == results[0][1]:
    margin = (H0_qr - H0_int) / H0_int * 100
    print(f"Paley margin: {margin:.4f}% (was 8% at p=7, 2.2% at p=11)")
