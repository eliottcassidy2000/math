#!/usr/bin/env python3
"""
Search for H=21 at n=6 — corrected Ω computation.
opus-2026-03-14-S71h

At n=5: max |Ω| = 7 (giving H=15). We need |Ω| large enough for I(Ω,2)=21.
At n=6: check all 2^15 = 32768 tournaments.
"""

from itertools import combinations, permutations
from collections import Counter

def find_all_directed_odd_cycles(n, adj, max_len=None):
    if max_len is None:
        max_len = n
    cycles = []
    seen = set()
    for length in range(3, max_len + 1, 2):
        for combo in combinations(range(n), length):
            for perm in permutations(combo):
                is_cycle = True
                for i in range(length):
                    nxt = (i + 1) % length
                    if perm[nxt] not in adj[perm[i]]:
                        is_cycle = False
                        break
                if is_cycle:
                    rotations = [perm[i:] + perm[:i] for i in range(length)]
                    canonical = min(rotations)
                    if canonical not in seen:
                        seen.add(canonical)
                        cycles.append((canonical, frozenset(combo)))
    return cycles

def tournament_from_bits(n, bits):
    adj = [set() for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i].add(j)
            else:
                adj[j].add(i)
            idx += 1
    return adj

def hp_count(n, adj):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in adj[v]:
                if not (mask & (1 << w)):
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[full])

def independence_poly_at_2(num_verts, edges):
    adj = set()
    for u, v in edges:
        adj.add((u, v))
        adj.add((v, u))
    total = 0
    for mask in range(1 << num_verts):
        verts = [i for i in range(num_verts) if mask & (1 << i)]
        k = len(verts)
        ok = True
        for a in range(k):
            for b in range(a+1, k):
                if (verts[a], verts[b]) in adj:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            total += 2**k
    return total

print("=" * 70)
print("H SPECTRUM AND Ω STRUCTURE AT n=6")
print("=" * 70)
print()

n = 6
num_edges = 15
total = 1 << num_edges  # 32768

h_counter = Counter()
h_omega_size = {}
h21_found = False

# This will be slow with full cycle enumeration including 5-cycles
# Let's first just compute H values to find the spectrum
print("Phase 1: H value spectrum at n=6")
for bits in range(total):
    adj = tournament_from_bits(n, bits)
    H = hp_count(n, adj)
    h_counter[H] += 1

    if H == 21:
        h21_found = True
        print(f"  H=21 FOUND at bits={bits}!")
        # Get cycle structure
        cycles = find_all_directed_odd_cycles(n, adj, max_len=5)
        nc = len(cycles)
        omega_edges = []
        for i in range(nc):
            for j in range(i+1, nc):
                if cycles[i][1] & cycles[j][1]:
                    omega_edges.append((i, j))
        I_omega = independence_poly_at_2(nc, omega_edges) if nc > 0 else 1
        print(f"    |Ω|={nc}, |E(Ω)|={len(omega_edges)}, I(Ω,2)={I_omega}")
        # Show cycles
        c3 = sum(1 for c, _ in cycles if len(c) == 3)
        c5 = sum(1 for c, _ in cycles if len(c) == 5)
        print(f"    3-cycles: {c3}, 5-cycles: {c5}")

print()
if not h21_found:
    print("H=21 NOT found at n=6!")
print()

print("H value distribution at n=6:")
for h in sorted(h_counter.keys()):
    pct = 100 * h_counter[h] / total
    marker = " ← FORBIDDEN!" if h in (7, 21) else ""
    if pct >= 0.5 or h <= 30 or h in (7, 21, 63):
        print(f"  H={h:3d}: {h_counter[h]:5d} ({pct:5.1f}%){marker}")

print()
print(f"Total distinct H values at n=6: {len(h_counter)}")
print(f"Range: {min(h_counter.keys())} to {max(h_counter.keys())}")

# Check: is 7 in the spectrum?
if 7 in h_counter:
    print(f"\nH=7 FOUND at n=6: {h_counter[7]} times!")
else:
    print(f"\nH=7 NOT in n=6 spectrum (as expected)")

if 21 in h_counter:
    print(f"H=21 FOUND at n=6: {h_counter[21]} times!")
else:
    print(f"H=21 NOT in n=6 spectrum")

print()
print("=" * 70)
print("KEY INSIGHT: What Ω structures arise at n=6?")
print("=" * 70)
print()

# For a few H values near 21, check Ω structure
for target_h in [17, 19, 21, 23, 25]:
    if target_h not in h_counter:
        print(f"H={target_h}: absent from n=6 spectrum")
        continue

    print(f"H={target_h}: checking Ω structures (first 3 examples)...")
    count = 0
    for bits in range(total):
        adj = tournament_from_bits(n, bits)
        H = hp_count(n, adj)
        if H != target_h:
            continue
        count += 1
        if count <= 3:
            cycles = find_all_directed_odd_cycles(n, adj, max_len=5)
            nc = len(cycles)
            omega_edges = []
            for i in range(nc):
                for j in range(i+1, nc):
                    if cycles[i][1] & cycles[j][1]:
                        omega_edges.append((i, j))
            c3 = sum(1 for c, _ in cycles if len(c) == 3)
            c5 = sum(1 for c, _ in cycles if len(c) == 5)
            print(f"  bits={bits}: |Ω|={nc} ({c3} tri + {c5} pent), |E|={len(omega_edges)}")
    print()
