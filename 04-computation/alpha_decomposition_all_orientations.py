#!/usr/bin/env python3
"""
alpha_decomposition_all_orientations.py -- kind-pasteur-2026-03-13-S60

Compare H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 across ALL 32 orientations
at p=11, to understand how each term contributes to H-maximization.

Key questions:
1. Does Paley maximize alpha_1 (total cycles)?
2. Does Paley minimize alpha_2 (disjoint pairs)?
3. How does alpha_3 vary?
4. What's the relative contribution of each term?
"""

from itertools import combinations
from collections import defaultdict

def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A

def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    start = 0
    dp = {(1 << start, start): 1}
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total

def compute_H_heldkarp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


p = 11
m = (p - 1) // 2

print(f"={'='*70}")
print(f"  ALPHA DECOMPOSITION: ALL {1 << m} ORIENTATIONS AT p={p}")
print(f"={'='*70}")

results = []

for bits in range(1 << m):
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j + 1)
        else:
            S.append(p - (j + 1))

    A = build_adj(p, S)
    H = compute_H_heldkarp(A, p)

    # Get all active vertex sets
    active_vsets = []
    for k in range(3, p+1, 2):
        for subset in combinations(range(p), k):
            fs = frozenset(subset)
            nc = count_ham_cycles(A, list(subset))
            if nc > 0:
                active_vsets.append((fs, k, nc))

    # Per-length cycle counts
    c_k = defaultdict(int)
    for _, k, nc in active_vsets:
        c_k[k] += nc

    N_total = sum(nc for _, _, nc in active_vsets)

    # alpha_2: disjoint pairs weighted by cycle counts
    alpha_2 = 0
    disjoint_pairs = []
    for i in range(len(active_vsets)):
        for j in range(i+1, len(active_vsets)):
            V1, k1, n1 = active_vsets[i]
            V2, k2, n2 = active_vsets[j]
            if not (V1 & V2):
                alpha_2 += n1 * n2
                disjoint_pairs.append((i, j))

    # alpha_3: disjoint triples weighted by cycle counts
    alpha_3 = 0
    for i, j in disjoint_pairs:
        V1, k1, n1 = active_vsets[i]
        V2, k2, n2 = active_vsets[j]
        used = V1 | V2
        for l in range(j+1, len(active_vsets)):
            V3, k3, n3 = active_vsets[l]
            if not (used & V3):
                alpha_3 += n1 * n2 * n3

    H_check = 1 + 2*N_total + 4*alpha_2 + 8*alpha_3

    results.append({
        'bits': bits, 'S': S, 'H': H,
        'N': N_total, 'alpha_2': alpha_2, 'alpha_3': alpha_3,
        'H_check': H_check, 'c_k': dict(sorted(c_k.items())),
        'n_active': len(active_vsets), 'n_disj_pairs': len(disjoint_pairs)
    })

    if bits % 4 == 0:
        print(f"  Progress: {bits+1}/{1 << m} orientations done...")

# Sort by H
results.sort(key=lambda r: r['H'], reverse=True)

print(f"\n{'='*70}")
print(f"  RESULTS (sorted by H descending)")
print(f"{'='*70}")
print(f"{'bits':>5} {'H':>8} {'N_total':>8} {'2N':>8} {'4*a2':>8} {'8*a3':>8} {'a1%':>6} {'a2%':>6} {'a3%':>6}")
print(f"{'-'*70}")

for r in results:
    term1 = 2 * r['N']
    term2 = 4 * r['alpha_2']
    term3 = 8 * r['alpha_3']
    total = term1 + term2 + term3  # without the +1
    pct1 = 100 * term1 / total if total > 0 else 0
    pct2 = 100 * term2 / total if total > 0 else 0
    pct3 = 100 * term3 / total if total > 0 else 0
    print(f"{r['bits']:>5} {r['H']:>8} {r['N']:>8} {term1:>8} {term2:>8} {term3:>8} {pct1:>5.1f}% {pct2:>5.1f}% {pct3:>5.1f}%")

# Statistics
Hs = [r['H'] for r in results]
Ns = [r['N'] for r in results]
a2s = [r['alpha_2'] for r in results]
a3s = [r['alpha_3'] for r in results]

print(f"\n{'='*70}")
print(f"  STATISTICS")
print(f"{'='*70}")
print(f"  H:  min={min(Hs)}, max={max(Hs)}, ratio={max(Hs)/min(Hs):.3f}")
print(f"  N:  min={min(Ns)}, max={max(Ns)}, ratio={max(Ns)/min(Ns):.3f}")
print(f"  a2: min={min(a2s)}, max={max(a2s)}, ratio={max(a2s)/min(a2s):.3f}")
print(f"  a3: min={min(a3s)}, max={max(a3s)}, ratio={'inf' if min(a3s)==0 else f'{max(a3s)/min(a3s):.3f}'}")

# Correlations
import math
def corr(xs, ys):
    n = len(xs)
    mx, my = sum(xs)/n, sum(ys)/n
    cov = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sx = math.sqrt(sum((x-mx)**2 for x in xs))
    sy = math.sqrt(sum((y-my)**2 for y in ys))
    if sx == 0 or sy == 0:
        return float('nan')
    return cov / (sx * sy)

print(f"\n  Correlations with H:")
print(f"    corr(H, N_total)  = {corr(Hs, Ns):.6f}")
print(f"    corr(H, alpha_2)  = {corr(Hs, a2s):.6f}")
print(f"    corr(H, alpha_3)  = {corr(Hs, a3s):.6f}")

# Cross correlations
print(f"\n  Cross correlations:")
print(f"    corr(N, alpha_2)  = {corr(Ns, a2s):.6f}")
print(f"    corr(N, alpha_3)  = {corr(Ns, a3s):.6f}")
print(f"    corr(alpha_2, alpha_3) = {corr(a2s, a3s):.6f}")

# Identify Paley
paley_bits = None
S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
for r in results:
    if sorted(r['S']) == S_qr:
        paley_bits = r['bits']
        break

print(f"\n  Paley (QR={S_qr}): bits={paley_bits}")
for r in results:
    if r['bits'] == paley_bits:
        print(f"    H={r['H']}, N={r['N']}, alpha_2={r['alpha_2']}, alpha_3={r['alpha_3']}")
        print(f"    Rank: H=#1, N={'#1' if r['N'] == max(Ns) else f'#{sorted(Ns, reverse=True).index(r[chr(78)]) + 1}'}")

# Who has max N but not max H?
max_N_bits = max(results, key=lambda r: r['N'])
max_H_bits = max(results, key=lambda r: r['H'])
print(f"\n  Max H: bits={max_H_bits['bits']}, H={max_H_bits['H']}, N={max_H_bits['N']}")
print(f"  Max N: bits={max_N_bits['bits']}, H={max_N_bits['H']}, N={max_N_bits['N']}")

# Regression: H = a + b*N + c*alpha_2 + d*alpha_3
# Since H = 1 + 2*N + 4*a2 + 8*a3, this should be exact
print(f"\n  Verification: all H = 1 + 2N + 4a2 + 8a3?")
all_match = all(r['H'] == r['H_check'] for r in results)
print(f"    {'YES — exact for all 32 orientations' if all_match else 'NO — some mismatch!'}")

# Term dominance analysis
print(f"\n{'='*70}")
print(f"  TERM DOMINANCE ANALYSIS")
print(f"{'='*70}")

# For each pair of orientations, which term drives H difference?
top = results[0]  # highest H (should be Paley)
# Find distinct H values
seen_H = set()
comparisons = []
for r in results:
    if r['H'] not in seen_H and r['H'] != top['H']:
        seen_H.add(r['H'])
        comparisons.append(r)

for r in comparisons:
    dH = top['H'] - r['H']
    dN = 2*(top['N'] - r['N'])
    da2 = 4*(top['alpha_2'] - r['alpha_2'])
    da3 = 8*(top['alpha_3'] - r['alpha_3'])
    print(f"\n  Paley vs bits={r['bits']} (H={r['H']}): delta_H = {dH}")
    print(f"    delta(2N)  = {dN:>+8} ({100*dN/dH:>+6.1f}%)")
    print(f"    delta(4a2) = {da2:>+8} ({100*da2/dH:>+6.1f}%)")
    print(f"    delta(8a3) = {da3:>+8} ({100*da3/dH:>+6.1f}%)")

# Already covered all distinct H values above

# Marginal analysis: if we keep alpha_1 fixed at Paley value, what alpha_2 range?
# And vice versa
print(f"\n{'='*70}")
print(f"  TRADE-OFF ANALYSIS")
print(f"{'='*70}")

# Group by c_k signature
sig_groups = defaultdict(list)
for r in results:
    sig = tuple(sorted(r['c_k'].items()))
    sig_groups[sig].append(r)

print(f"\n  Distinct c_k signatures: {len(sig_groups)}")
for sig, group in sorted(sig_groups.items(), key=lambda x: -x[1][0]['H']):
    Hs_g = [r['H'] for r in group]
    if len(group) > 1:
        print(f"    c_k={dict(sig)}: {len(group)} orientations, H in [{min(Hs_g)}, {max(Hs_g)}]")
    else:
        print(f"    c_k={dict(sig)}: 1 orientation, H={Hs_g[0]}")

print("\nDONE.")
