#!/usr/bin/env python3
"""
spectral_cycle_bridge.py -- kind-pasteur-2026-03-13-S60

Connect the spectral invariant sum|lambda|^{2k} to cycle counts.

For a tournament T with adjacency matrix A:
- tr(A^k) = number of CLOSED walks of length k
- For tournaments: these closed walks include k-cycles but also
  backtracking walks, so it's not just c_k.

Key: tr((AA^T)^k) = sum |lambda_j|^{2k}
- k=1: tr(AA^T) = sum deg_out(v) = p*(p-1)/2 = 55
- k=2: tr((AA^T)^2) = sum |lambda|^4

For a tournament: AA^T = matrix where (i,j) entry = |{v : i->v and j->v}|
= number of common out-neighbors. For regular tournament: diagonal = (p-1)/2.

Actually, AA^T for circulant: entry (i,j) depends only on i-j.
So tr((AA^T)^2) = sum_d |f(d)|^4 where f(d) = |S cap (S+d)|.
Hmm, this gets complicated.

Let's just compute the correlation empirically.
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
omega = cmath.exp(2j * cmath.pi / p)

print(f"SPECTRAL-CYCLE BRIDGE AT p={p}")
print(f"{'='*70}")

def bits_to_S(bits):
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j + 1)
        else:
            S.append(p - (j + 1))
    return S

# Collect both spectral and cycle data for all orientations
all_data = []
for bits in range(1 << m):
    S = bits_to_S(bits)
    A = build_adj(p, S)

    # Eigenvalues
    eigs = [sum(omega**(j*s) for s in S) for j in range(p)]
    sum_4 = sum(abs(e)**4 for e in eigs).real
    sum_6 = sum(abs(e)**6 for e in eigs).real
    prod_abs = 1
    for e in eigs:
        prod_abs *= abs(e)
    prod_abs = abs(prod_abs)

    # Cycle counts
    c_k = {}
    for k in range(3, p+1, 2):
        total = 0
        for subset in combinations(range(p), k):
            total += count_ham_cycles(A, list(subset))
        c_k[k] = total

    H = compute_H_heldkarp(A, p)
    N = sum(c_k.values())

    all_data.append({
        'bits': bits, 'H': H, 'N': N, 'c_k': c_k,
        'sum_4': sum_4, 'sum_6': sum_6, 'prod': prod_abs
    })

    if bits % 8 == 0:
        print(f"  Progress: {bits+1}/{1 << m}...")

# Correlations
def corr(xs, ys):
    n = len(xs)
    mx, my = sum(xs)/n, sum(ys)/n
    cov = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sx = math.sqrt(sum((x-mx)**2 for x in xs))
    sy = math.sqrt(sum((y-my)**2 for y in ys))
    if sx == 0 or sy == 0:
        return float('nan')
    return cov / (sx * sy)

Hs = [d['H'] for d in all_data]
Ns = [d['N'] for d in all_data]
s4s = [d['sum_4'] for d in all_data]
s6s = [d['sum_6'] for d in all_data]
prods = [d['prod'] for d in all_data]
c5s = [d['c_k'][5] for d in all_data]
c7s = [d['c_k'][7] for d in all_data]
c9s = [d['c_k'][9] for d in all_data]
c11s = [d['c_k'][11] for d in all_data]

print(f"\n{'='*70}")
print(f"  CORRELATIONS")
print(f"{'='*70}")

print(f"\n  With H:")
print(f"    corr(H, sum_4) = {corr(Hs, s4s):.6f}")
print(f"    corr(H, sum_6) = {corr(Hs, s6s):.6f}")
print(f"    corr(H, prod)  = {corr(Hs, prods):.6f}")
print(f"    corr(H, N)     = {corr(Hs, Ns):.6f}")
print(f"    corr(H, c5)    = {corr(Hs, c5s):.6f}")
print(f"    corr(H, c7)    = {corr(Hs, c7s):.6f}")
print(f"    corr(H, c9)    = {corr(Hs, c9s):.6f}")
print(f"    corr(H, c11)   = {corr(Hs, c11s):.6f}")

print(f"\n  sum_4 with cycles:")
print(f"    corr(sum_4, N)  = {corr(s4s, Ns):.6f}")
print(f"    corr(sum_4, c5) = {corr(s4s, c5s):.6f}")
print(f"    corr(sum_4, c7) = {corr(s4s, c7s):.6f}")
print(f"    corr(sum_4, c9) = {corr(s4s, c9s):.6f}")
print(f"    corr(sum_4, c11)= {corr(s4s, c11s):.6f}")

# Key relationship: is sum_4 = f(c_k)?
# sum_4 = tr((AA^T)^2) counts 4-step walks in the "common neighbor" graph
# For tournaments: (AA^T)_{ij} = |{v : i->v AND j->v}| = d+(i,j)
# So sum_4 = sum_{i,j} d+(i,j)^2

# For regular tournament: d+(i,i) = (p-1)/2.
# Off-diagonal: sum_j d+(i,j) = ... complicated.

# Let me compute (AA^T)_{ij} directly for one case
print(f"\n{'='*70}")
print(f"  COMMON NEIGHBOR MATRIX ANALYSIS")
print(f"{'='*70}")

for label, bits in [("Paley", 29), ("ClassC", 0)]:
    S = bits_to_S(bits)
    A = build_adj(p, S)

    # Compute AA^T
    AAT = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            AAT[i][j] = sum(A[i][v] * A[j][v] for v in range(p))

    # For circulant: AAT depends only on i-j
    aat_vals = {}
    for d in range(p):
        aat_vals[d] = AAT[0][d]

    print(f"\n  {label} (bits={bits}):")
    print(f"    (AA^T)_d: ", end="")
    for d in range(p):
        print(f"d={d}:{aat_vals[d]} ", end="")
    print()

    # sum_4 = p * sum_d aat_d^2 (for circulant)
    sum_4_from_aat = p * sum(v**2 for v in aat_vals.values())
    eigs = [sum(omega**(j*s) for s in S) for j in range(p)]
    sum_4_from_eigs = sum(abs(e)**4 for e in eigs).real

    print(f"    sum_4 from AAT: p * sum d^2 = {p} * {sum(v**2 for v in aat_vals.values())} = {sum_4_from_aat}")
    print(f"    sum_4 from eigenvalues: {sum_4_from_eigs:.1f}")
    print(f"    Match: {abs(sum_4_from_aat - sum_4_from_eigs) < 0.1}")

    # tr(A^3) = 3 * c_3(directed) (counts directed 3-cycles)
    trA3 = sum(sum(sum(A[i][j]*A[j][k]*A[k][i] for k in range(p)) for j in range(p)) for i in range(p))
    print(f"    tr(A^3) = {trA3} (expected: 3 * c3_directed = 3 * {2*55} = {3*2*55})")

# KEY STRUCTURAL RELATION: sum_4 relates to the 4-vertex subgraph counts
# (not directly to c_k for k > 4, but to lower-order invariants).

# At p=11: sum_4 = p * sum_d (AA^T)(0,d)^2
# = p * [(p-1)/2)^2 + sum_{d!=0} (# common out-neighbors of 0 and d)^2]

# For regular tournament: (AA^T)(0,0) = (p-1)/2 = 5 (out-degree)
# (AA^T)(0,d) = number of v with 0->v and d->v = number of common successors

# For Paley: by the 2-transitivity of the Paley tournament,
# all off-diagonal entries of AA^T are equal: each = (p-3)/4
# This gives the flattest possible common-neighbor distribution.

# Check:
print(f"\n{'='*70}")
print(f"  PALEY FLAT COMMON-NEIGHBOR PROPERTY")
print(f"{'='*70}")

S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
A = build_adj(p, S_qr)
AAT = [[0]*p for _ in range(p)]
for i in range(p):
    for j in range(p):
        AAT[i][j] = sum(A[i][v] * A[j][v] for v in range(p))

off_diag_vals = set()
for d in range(1, p):
    off_diag_vals.add(AAT[0][d])

print(f"  Paley off-diagonal AA^T values: {sorted(off_diag_vals)}")
print(f"  Expected: (p-3)/4 = {(p-3)/4}")
print(f"  Flat? {len(off_diag_vals) == 1}")

# For flat common-neighbor: sum_4 = p * [((p-1)/2)^2 + (p-1)*((p-3)/4)^2]
flat_sum4 = p * (((p-1)//2)**2 + (p-1) * ((p-3)/4)**2)
print(f"  Predicted sum_4 = {flat_sum4}")
print(f"  Actual sum_4 = {sum(abs(e)**4 for e in [sum(omega**(j*s) for s in S_qr) for j in range(p)]).real:.1f}")

# The deviation from flatness in AA^T determines sum_4 - sum_4(flat)
# More deviation => higher sum_4 => fewer cycles => lower H
# This is the SPECTRAL GAP interpretation of H-maximization!

print(f"\n  KEY INSIGHT: Paley minimizes sum|lam|^4 because it has the")
print(f"  flattest common-neighbor distribution (2-transitivity).")
print(f"  Flatter spectrum => more uniform cycle distribution => higher H.")

# Cross-orientation AAT analysis
print(f"\n{'='*70}")
print(f"  COMMON-NEIGHBOR VARIANCE ACROSS ORIENTATIONS")
print(f"{'='*70}")

for label, bits in [("Paley", 29), ("ClassB", 4), ("ClassC", 0), ("ClassD", 1)]:
    S = bits_to_S(bits)
    A = build_adj(p, S)
    AAT_vals = []
    for d in range(1, p):
        AAT_vals.append(sum(A[0][v] * A[d][v] for v in range(p)))

    mean_aat = sum(AAT_vals) / len(AAT_vals)
    var_aat = sum((v - mean_aat)**2 for v in AAT_vals) / len(AAT_vals)
    print(f"  {label}: AAT_off_diag = {sorted(AAT_vals)}, mean={mean_aat:.2f}, var={var_aat:.4f}")

print("\nDONE.")
