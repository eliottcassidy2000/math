#!/usr/bin/env python3
"""
paley_cocktail_general_87b.py — opus-2026-03-14-S87b

Does the QR_p → CP(α₁/2) cocktail party theorem generalize?

At p=7: QR₇ has 14 3-cycles, conflict graph = CP(7) = K₁₄ - 7K₂
        α₁ = 14 = 2 × 7 = 2p

At p=3: QR₃ has 1 3-cycle, conflict graph = K₁ (trivial CP(1/2)?)

At p=11: QR₁₁ has α₁ = 55 3-cycles. If the pattern holds:
         conflict graph = CP(55/2)? But 55 is odd!
         So maybe it's CP(27.5) — impossible.

Actually: for regular tournaments on n vertices,
α₁ = C(n,3) - n×C((n-1)/2, 2) = n(n-1)(n-2)/6 - n×(n-1)(n-3)/8
For p=7: α₁ = 35 - 21 = 14 = 2 × 7
For p=11: α₁ = 165 - 110 = 55 (ODD!)

So CP structure requires α₁ even. Let's check what happens at p=11.
"""

import numpy as np
from itertools import combinations
from collections import Counter

def build_circulant(n, conn_set):
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for d in conn_set:
            adj[i][(i+d)%n] = 1
    return adj

def find_3cycles(adj, n):
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cycles.append(((a,b,c), frozenset({a,b,c})))
        elif adj[a][c] and adj[c][b] and adj[b][a]:
            cycles.append(((a,c,b), frozenset({a,b,c})))
    return cycles

# ══════════════════════════════════════════════════════════════════
# General formula for α₁ in regular tournaments
# ══════════════════════════════════════════════════════════════════

print("=" * 70)
print("α₁ FOR REGULAR TOURNAMENTS ON n VERTICES")
print("=" * 70)

for p in [3, 5, 7, 9, 11, 13]:
    s = (p-1)//2
    alpha1 = p*(p-1)*(p-2)//6 - p*s*(s-1)//2
    print(f"  n={p:2d}: α₁ = {alpha1:4d}, α₁ mod 2 = {alpha1 % 2}, α₁/n = {alpha1/p:.1f}")

# α₁ = n(n²-3n+2)/6 - n(n-1)(n-3)/8
# = n[(n²-3n+2)/6 - (n-1)(n-3)/8]
# = n[(4(n²-3n+2) - 3(n-1)(n-3))/24]
# = n[(4n²-12n+8 - 3n²+12n-9)/24]
# = n[(n²-1)/24]
# = n(n²-1)/24

print("\nSimplified: α₁ = n(n²-1)/24")
for p in [3, 5, 7, 9, 11, 13]:
    alpha1 = p*(p*p-1)//24
    print(f"  n={p}: α₁ = {p}×{p*p-1}/24 = {p}×{(p*p-1)//24} = {alpha1}")

# When is n(n²-1)/24 even?
# n(n²-1) = n(n-1)(n+1) = product of 3 consecutive integers
# This is always divisible by 6. So α₁ = n(n-1)(n+1)/24.
# For this to be even: n(n-1)(n+1) ≡ 0 mod 48.
# n=3: 3×2×4=24, α₁=1 (odd)
# n=5: 5×4×6=120, α₁=5 (odd)
# n=7: 7×6×8=336, α₁=14 (even)
# n=9: 9×8×10=720, α₁=30 (even)
# n=11: 11×10×12=1320, α₁=55 (odd)
# n=13: 13×12×14=2184, α₁=91 (odd)

# Pattern: α₁ even iff 48 | n(n-1)(n+1)
# n=7: 48 | 336 ✓ (336/48=7)
# n=9: 48 | 720 ✓ (720/48=15)
# n=11: 1320/48 = 27.5 ✗

print("\nα₁ parity:")
for p in range(3, 20, 2):
    alpha1 = p*(p*p-1)//24
    is_even = alpha1 % 2 == 0
    print(f"  n={p:2d}: α₁={alpha1:4d}, even={is_even}")

# ══════════════════════════════════════════════════════════════════
# QR₁₁ conflict graph
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("QR₁₁: 3-CYCLE CONFLICT GRAPH")
print("=" * 70)

n = 11
qr11 = set()
for x in range(1, 11):
    qr11.add((x*x) % 11)
print(f"QR mod 11: {sorted(qr11)}")

adj_qr11 = build_circulant(11, qr11)
cycles_11 = find_3cycles(adj_qr11, n)
nc = len(cycles_11)
print(f"3-cycles: {nc}")

# Build conflict graph
print("Building conflict graph...")
cg = np.zeros((nc, nc), dtype=int)
for i in range(nc):
    for j in range(i+1, nc):
        if cycles_11[i][1] & cycles_11[j][1]:
            cg[i][j] = 1
            cg[j][i] = 1

degrees = [sum(cg[i]) for i in range(nc)]
deg_counter = Counter(degrees)
print(f"Degree distribution: {dict(sorted(deg_counter.items()))}")
is_regular = len(deg_counter) == 1
k = degrees[0] if is_regular else None
print(f"Regular graph: {is_regular}" + (f" with k={k}" if is_regular else ""))

# Compute complement degrees
comp = 1 - cg - np.eye(nc, dtype=int)
comp_degrees = [sum(comp[i]) for i in range(nc)]
comp_deg_counter = Counter(comp_degrees)
print(f"Complement degrees: {dict(sorted(comp_deg_counter.items()))}")

if is_regular:
    # Check SRG
    lam_vals = Counter()
    mu_vals = Counter()
    for i in range(nc):
        for j in range(i+1, nc):
            common = sum(1 for w in range(nc) if w != i and w != j
                         and cg[i][w] and cg[j][w])
            if cg[i][j]:
                lam_vals[common] += 1
            else:
                mu_vals[common] += 1

    print(f"λ values: {dict(lam_vals)}")
    print(f"μ values: {dict(mu_vals)}")

    if len(lam_vals) == 1 and len(mu_vals) == 1:
        lam = list(lam_vals.keys())[0]
        mu = list(mu_vals.keys())[0]
        print(f"\n✓ STRONGLY REGULAR: srg({nc}, {k}, {lam}, {mu})")

        # Is it a cocktail party graph?
        # CP(m) = srg(2m, 2m-2, 2m-4, 2m-2)
        m = nc // 2
        is_cp = (k == 2*m - 2 and lam == 2*m - 4 and mu == 2*m - 2) if nc % 2 == 0 else False
        print(f"  Cocktail party CP({m}): {is_cp}")
    else:
        print("NOT strongly regular")
else:
    # Spectrum anyway
    print("\nComputing spectrum...")
    eigenvalues = sorted(np.linalg.eigvalsh(cg.astype(float)), reverse=True)
    rounded = [round(e, 3) for e in eigenvalues]
    spec_counter = Counter(rounded)
    for val in sorted(spec_counter.keys(), reverse=True)[:15]:
        mult = spec_counter[val]
        print(f"  λ = {val:10.3f} (×{mult})")
    remaining = sum(v for k_, v in spec_counter.items()
                    if k_ not in list(sorted(spec_counter.keys(), reverse=True))[:15])
    if remaining > 0:
        print(f"  ... and {remaining} more")

# ══════════════════════════════════════════════════════════════════
# QR₅ — but p=5≡1 mod 4, so Paley doesn't exist
# ══════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("GENERAL PATTERN")
print("=" * 70)

print("""
RESULTS:
  n=3: α₁=1, trivial
  n=7: α₁=14, conflict = CP(7) = srg(14,12,10,12) ✓
  n=11: α₁=55, ...

For the cocktail party structure to hold:
  1. α₁ must be EVEN (so 55 at n=11 prevents CP)
  2. The conflict graph must be (n-2)-regular
  3. Every non-conflicting pair must be unique

The CP(7) theorem at n=7 may be specific to p=7:
  - p=7 is the smallest prime ≡ 3 mod 4 with p>3
  - α₁ = 14 = 2×7 (even, = 2p)
  - The Fano structure makes the matching work

For general Paley primes p ≡ 3 mod 4:
  α₁(QR_p) = p(p²-1)/24
  This is even iff 48 | p(p-1)(p+1)

  p=3: 24/48 → α₁=1 (odd)
  p=7: 336/48=7 → α₁=14 (even)
  p=11: 1320/48=27.5 → α₁=55 (odd)
  p=19: 6840/48=142.5 → α₁=285 (odd)
  p=23: 12144/48=253 → α₁=506 (even!)

  So CP structure can hold at p=7 and p=23 but not p=11 or p=19.
""")

# Let's check p=23 would have α₁ even
for p in [3, 7, 11, 19, 23, 31, 43, 47]:
    if pow(p-1, 1) % 4 == 2:  # p ≡ 3 mod 4
        alpha1 = p*(p*p-1)//24
        print(f"  p={p:2d} (≡3 mod 4): α₁={alpha1}, even={alpha1%2==0}")
