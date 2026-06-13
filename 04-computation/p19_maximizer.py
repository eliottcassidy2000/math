#!/usr/bin/env python3
"""
Explore H-maximization at p=19: which tournaments beat Paley?

Key discovery: the cyclic interval S={1,...,9} and arithmetic progressions
all achieve H=1184212824763 > Paley's 1172695746915.

This script tests structured connection sets to understand the pattern.
"""

import numpy as np
import time

def adjacency_matrix(S, p):
    A = np.zeros((p,p), dtype=np.int8)
    for i in range(p):
        for j in range(p):
            if i!=j and (j-i)%p in S:
                A[i][j] = 1
    return A

def count_hp_fast(A):
    n = len(A)
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp[mask, v]
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v, u]:
                    dp[mask | (1 << u), u] += c
    return int(np.sum(dp[full]))

def is_qr(a, p):
    if a % p == 0: return False
    return pow(a, (p-1)//2, p) == 1

p = 19
m = 9
paley_S = frozenset(j for j in range(1,p) if is_qr(j,p))
print(f"p = {p}, Paley QR = {sorted(paley_S)}")

# Test specific structured connection sets
tests = {
    'Paley QR': paley_S,
    'Paley complement': frozenset(range(1,p)) - paley_S,
    'Cyclic interval {1,...,9}': frozenset(range(1, 10)),
    'Cyclic interval {10,...,18}': frozenset(range(10, 19)),
    'Even {2,4,...,18}': frozenset(range(2, 19, 2)),
    'Odd {1,3,...,17}': frozenset(range(1, 18, 2)),
    'Step3 {1,4,7,10,13,16,2,5,8}': frozenset({1,4,7,10,13,16,2,5,8}),
}

# Add multiplicative orbits of the interval
for k in range(1, 19):
    S = frozenset((k * j) % p for j in range(1, 10))
    name = f'k*{{1,...,9}} k={k}'
    # Check validity
    if all((p - j) % p not in S for j in S) and len(S) == 9:
        tests[name] = S

# Add some Galois orbits of Paley
for k in range(2, 19):
    S = frozenset((k * j) % p for j in paley_S)
    name = f'k*Paley k={k}'
    if all((p - j) % p not in S for j in S) and len(S) == 9:
        tests[name] = S

print(f"\nTesting {len(tests)} connection sets...")
print(f"{'Name':<35} {'H':>20} {'vs Paley':>12}")
print("-" * 70)

results = {}
paley_H = None

for name, S in sorted(tests.items(), key=lambda x: x[0]):
    A = adjacency_matrix(S, p)
    t0 = time.time()
    H = count_hp_fast(A)
    elapsed = time.time() - t0

    if name == 'Paley QR':
        paley_H = H

    results[name] = H

for name in sorted(results.keys()):
    H = results[name]
    diff = H - paley_H if paley_H else 0
    marker = "***" if diff > 0 else ""
    print(f"  {name:<35} {H:>20} {diff:>+12} {marker}")

# Group by H value
by_H = {}
for name, H in results.items():
    by_H.setdefault(H, []).append(name)

print(f"\n{'='*60}")
print(f"DISTINCT H VALUES ({len(by_H)} total)")
print(f"{'='*60}")
for H in sorted(by_H.keys(), reverse=True):
    print(f"\n  H = {H}:")
    for name in by_H[H]:
        S = tests[name]
        print(f"    {name}: S = {sorted(S)}")

# Check: which H values appear and what's the pattern?
# Are the multiplicative orbits of {1,...,9} giving different H's?
print(f"\n{'='*60}")
print(f"ORBIT STRUCTURE")
print(f"{'='*60}")

print("\nMultiplicative orbits of {1,...,9}:")
orbit_H = {}
for k in range(1, 19):
    S = frozenset((k * j) % p for j in range(1, 10))
    if all((p - j) % p not in S for j in S) and len(S) == 9:
        name = f'k={k}'
        H = results.get(f'k*{{1,...,9}} k={k}', None)
        if H:
            orbit_H.setdefault(H, []).append(k)
            print(f"  k={k}: S={sorted(S)}, H={H}")

print("\nMultiplicative orbits of Paley QR:")
for k in range(1, 19):
    S = frozenset((k * j) % p for j in paley_S)
    if all((p - j) % p not in S for j in S) and len(S) == 9:
        name = f'k*Paley k={k}'
        H = results.get(name, None)
        if H:
            print(f"  k={k}: S={sorted(S)}, H={H}")

# === UNDERSTANDING: Why does the interval beat Paley? ===
print(f"\n{'='*60}")
print(f"EIGENVALUE COMPARISON: Interval vs Paley")
print(f"{'='*60}")

def eigenvalues_circulant(S, p):
    omega = np.exp(2j*np.pi/p)
    return [sum(omega**(k*s) for s in S) for k in range(p)]

eigs_paley = eigenvalues_circulant(paley_S, p)
eigs_interval = eigenvalues_circulant(frozenset(range(1,10)), p)

print("\nPaley eigenvalues (|lambda_k|):")
paley_mags = [abs(eigs_paley[k]) for k in range(1, 10)]
print(f"  {['%.4f'%x for x in paley_mags]}")
print(f"  All equal? {max(paley_mags) - min(paley_mags) < 0.001}")
print(f"  |lambda| = sqrt((p+1)/4) = {np.sqrt((p+1)/4):.4f}")

print("\nInterval eigenvalues (|lambda_k|):")
interval_mags = [abs(eigs_interval[k]) for k in range(1, 10)]
print(f"  {['%.4f'%x for x in interval_mags]}")
print(f"  Range: [{min(interval_mags):.4f}, {max(interval_mags):.4f}]")
print(f"  Var: {np.var(interval_mags):.4f}")

# Closed form for interval eigenvalues
# S = {1,...,m}: lambda_k = sum_{s=1}^{m} omega^{ks}
# = omega^k * (1 - omega^{mk}) / (1 - omega^k) = sin(m*pi*k/p)/sin(pi*k/p) * e^{i*phase}
print("\nInterval eigenvalues (analytic):")
for k in range(1, 10):
    lam = eigs_interval[k]
    mag = abs(lam)
    phase = np.angle(lam) / np.pi
    print(f"  k={k}: |lambda| = {mag:.6f}, phase/pi = {phase:.6f}")

print("\nPaley eigenvalues (analytic):")
for k in range(1, 10):
    lam = eigs_paley[k]
    mag = abs(lam)
    phase = np.angle(lam) / np.pi
    print(f"  k={k}: |lambda| = {mag:.6f}, phase/pi = {phase:.6f}")

# Power sums comparison
print(f"\nPower sums p_j = sum |lambda_k|^j:")
for j in range(2, 8):
    ps_paley = sum(abs(eigs_paley[k])**j for k in range(1, p))
    ps_interval = sum(abs(eigs_interval[k])**j for k in range(1, p))
    print(f"  j={j}: Paley={ps_paley:.2f}, Interval={ps_interval:.2f}, "
          f"diff={ps_interval-ps_paley:+.2f}")
