#!/usr/bin/env python3
"""
pi_ip_fast_89c.py — IP structure of G(P_p), FAST version
opus-2026-03-14-S89c

Only enumerate 3-cycles and 5-cycles (tractable),
then compute first few IP coefficients.
"""

from itertools import combinations
from collections import Counter

def paley_tournament(p):
    qr = set()
    for a in range(1, p):
        qr.add((a*a) % p)
    adj = {i: set() for i in range(p)}
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                adj[i].add(j)
    return adj, qr

def find_directed_cycles_by_length(adj, n, length):
    """Find directed cycles of exactly given length, canonicalized by min vertex."""
    cycles = []
    for start in range(n):
        stack = [(start, [start], 1 << start)]
        while stack:
            v, path, mask = stack.pop()
            if len(path) == length:
                if start in adj[v]:
                    cycles.append(tuple(path))
                continue
            for u in adj[v]:
                if u < start:
                    continue
                if mask & (1 << u):
                    continue
                stack.append((u, path + [u], mask | (1 << u)))
    return cycles

print("=" * 70)
print("PART 1: Odd cycle counts by length")
print("=" * 70)

for p in [3, 7, 11]:
    adj, qr = paley_tournament(p)
    print(f"\n  P_{p}:")
    total_cycles = 0
    cycle_data = {}
    for length in range(3, p + 1, 2):
        cycles = find_directed_cycles_by_length(adj, p, length)
        cycle_data[length] = cycles
        total_cycles += len(cycles)
        print(f"    t_{length} = {len(cycles)}")

    print(f"    Total odd cycles: {total_cycles}")

    # For P_7: 14 three-cycles + 42 five-cycles + 24 seven-cycles = 80
    # For P_11: need to compute...

print()
print("=" * 70)
print("PART 2: G(P_7) independence polynomial (full)")
print("=" * 70)

p = 7
adj, qr = paley_tournament(p)
all_cycles_7 = []
for length in range(3, 8, 2):
    all_cycles_7.extend(find_directed_cycles_by_length(adj, p, length))

n_cyc = len(all_cycles_7)
print(f"\n  P_7: {n_cyc} odd cycles total")

# Build disjointness graph
g_adj = {i: set() for i in range(n_cyc)}
for i in range(n_cyc):
    for j in range(i+1, n_cyc):
        if not set(all_cycles_7[i]) & set(all_cycles_7[j]):
            g_adj[i].add(j)
            g_adj[j].add(i)

# Degree sequence
degrees = [len(g_adj[i]) for i in range(n_cyc)]
print(f"  G(P_7): {n_cyc} vertices, {sum(degrees)//2} edges")
print(f"  Degree range: {min(degrees)}-{max(degrees)}, mean {sum(degrees)/n_cyc:.1f}")

# Compute IP coefficients by enumerating independent sets
ip_coeffs = [0] * (n_cyc + 1)
for size in range(n_cyc + 1):
    count = 0
    for subset in combinations(range(n_cyc), size):
        indep = True
        for a in range(len(subset)):
            for b in range(a+1, len(subset)):
                if subset[b] in g_adj[subset[a]]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            count += 1
    ip_coeffs[size] = count
    if size > 0 and count == 0:
        ip_coeffs = ip_coeffs[:size]
        break

print(f"\n  IP(G(P_7), x) coefficients: {ip_coeffs}")
print(f"  IP polynomial: ", end="")
terms = []
for i, c in enumerate(ip_coeffs):
    if c > 0:
        if i == 0:
            terms.append(str(c))
        elif i == 1:
            terms.append(f"{c}x")
        else:
            terms.append(f"{c}x^{i}")
print(" + ".join(terms))

# Evaluate at various x
print(f"\n  Evaluations:")
for x in [-2, -1, 0, 1, 2, 3]:
    val = sum(c * x**i for i, c in enumerate(ip_coeffs))
    print(f"    IP(G, {x:2d}) = {val}")

# IP(G, 2) should equal H(P_7) = 189
ip_at_2 = sum(c * 2**i for i, c in enumerate(ip_coeffs))
print(f"\n  IP(G(P_7), 2) = {ip_at_2}, H(P_7) = 189, match: {'✓' if ip_at_2 == 189 else '✗'}")

# IP(G, -1): the alternating sum
ip_at_neg1 = sum(c * (-1)**i for i, c in enumerate(ip_coeffs))
print(f"  IP(G(P_7), -1) = {ip_at_neg1}")

# Zeros of IP
import numpy as np
if len(ip_coeffs) > 1:
    # Polynomial: ip_coeffs[0] + ip_coeffs[1]x + ...
    # numpy roots expects highest degree first
    coeffs_rev = list(reversed(ip_coeffs))
    roots = np.roots(coeffs_rev)
    print(f"\n  Zeros of IP(G(P_7), x):")
    for r in sorted(roots, key=lambda x: (x.real, x.imag)):
        if abs(r.imag) < 1e-10:
            print(f"    x = {r.real:.6f}")
        else:
            print(f"    x = {r.real:.6f} + {r.imag:.6f}i")

print()
print("=" * 70)
print("PART 3: P_7 orbit decomposition of G")
print("=" * 70)

# Group cycles into orbits under σ: i ↦ i+1 mod 7
def normalize_cycle(c):
    """Canonical form: tuple sorted by vertex."""
    return tuple(sorted(c))

def shift(c, p):
    return tuple(sorted([(v+1) % p for v in c]))

seen = set()
orbits = []
for c in all_cycles_7:
    nc = normalize_cycle(c)
    if nc in seen:
        continue
    orbit = []
    current = nc
    for _ in range(p):
        orbit.append(current)
        seen.add(current)
        current = shift(current, p)
    orbits.append(orbit)

print(f"\n  {len(orbits)} orbits of odd cycles under Z/7Z:")
for i, orbit in enumerate(orbits):
    rep = orbit[0]
    length = len(rep)
    print(f"    Orbit {i}: {len(orbit)} cycles of length {length}, rep = {rep}")

# In the quotient G/Z_p, orbits are vertices.
# Two orbits are "adjacent" if some cycle in one is disjoint from some cycle in the other.
# By symmetry: either ALL pairs between orbits are disjoint, or they follow a pattern.

print(f"\n  Quotient adjacency (orbit pairs with disjoint representatives):")
for i in range(len(orbits)):
    for j in range(i+1, len(orbits)):
        # Check if any representative from orbit i is disjoint from any in orbit j
        disjoint_count = 0
        total_pairs = len(orbits[i]) * len(orbits[j])
        for ci in orbits[i]:
            for cj in orbits[j]:
                if not set(ci) & set(cj):
                    disjoint_count += 1
        if disjoint_count > 0:
            print(f"    O_{i} -- O_{j}: {disjoint_count}/{total_pairs} disjoint pairs")

print()
print("=" * 70)
print("PART 4: First few IP coefficients for P_11")
print("=" * 70)

p = 11
adj, qr = paley_tournament(p)

# Only 3-cycles for P_11 (fast)
three_cycles = find_directed_cycles_by_length(adj, p, 3)
five_cycles = find_directed_cycles_by_length(adj, p, 5)
print(f"\n  P_11: {len(three_cycles)} 3-cycles, {len(five_cycles)} 5-cycles")

# Just using 3-cycles for a quick IP estimate:
# Disjoint 3-cycle pairs
disjoint_3 = 0
for i in range(len(three_cycles)):
    for j in range(i+1, len(three_cycles)):
        if not set(three_cycles[i]) & set(three_cycles[j]):
            disjoint_3 += 1
print(f"  Disjoint 3-cycle pairs: {disjoint_3}")

# Disjoint 3-cycle triples
triple_count = 0
for i in range(len(three_cycles)):
    si = set(three_cycles[i])
    for j in range(i+1, len(three_cycles)):
        if si & set(three_cycles[j]):
            continue
        sj = set(three_cycles[j])
        for k in range(j+1, len(three_cycles)):
            sk = set(three_cycles[k])
            if not (si & sk) and not (sj & sk):
                triple_count += 1
print(f"  Disjoint 3-cycle triples: {triple_count}")

# IP contribution from 3-cycles alone:
# 1 + 55×2 + disjoint_3×4 + triple_count×8
ip_3only = 1 + len(three_cycles)*2 + disjoint_3*4 + triple_count*8
print(f"\n  IP contribution from 3-cycles at x=2:")
print(f"    1 + {len(three_cycles)}×2 + {disjoint_3}×4 + {triple_count}×8")
print(f"    = {ip_3only}")
print(f"    H(P_11) = 95095, so 5-cycles and higher contribute {95095 - ip_3only}")

# 5-cycle disjoint pairs
print(f"\n  Computing 5-cycle disjointness...")
disjoint_5 = 0
for i in range(len(five_cycles)):
    for j in range(i+1, len(five_cycles)):
        if not set(five_cycles[i]) & set(five_cycles[j]):
            disjoint_5 += 1
print(f"  Disjoint 5-cycle pairs: {disjoint_5}")

# 3-cycle and 5-cycle disjoint pairs
cross_disjoint = 0
for c3 in three_cycles:
    s3 = set(c3)
    for c5 in five_cycles:
        if not s3 & set(c5):
            cross_disjoint += 1
print(f"  Disjoint (3-cycle, 5-cycle) pairs: {cross_disjoint}")

print(f"\n  Full IP_2 (at level 2, counting disjoint PAIRS of any type):")
total_pairs = disjoint_3 + disjoint_5 + cross_disjoint
print(f"    3-3: {disjoint_3}, 5-5: {disjoint_5}, 3-5: {cross_disjoint}")
print(f"    Total = {total_pairs}")

print()
print("=" * 70)
print("PART 5: Cycle counts mod p — orbits confirmation")
print("=" * 70)

for p in [7, 11]:
    adj, qr = paley_tournament(p)
    for length in [3, 5]:
        cycles = find_directed_cycles_by_length(adj, p, length)
        print(f"  P_{p}: t_{length} = {len(cycles)}, mod p = {len(cycles) % p}")

    # Hamiltonian cycles (length p)
    # Already computed: t_7(P_7) = 24, t_11(P_11) = 5505
    # t_7 mod 7 = 3 = (p-1)/2 (from THM-214)
    # t_11 mod 11 = 5 = (p-1)/2 (from THM-214)

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)

print("""
  THE INDEPENDENCE POLYNOMIAL VIEWPOINT

  H(P_p) = IP(G(P_p), 2) where G = odd-cycle disjointness graph

  G(P_7): 80 vertices (14+42+24 cycles of length 3,5,7)
          in 80/7 ≈ 11-12 orbits under Z/7Z

  IP(G(P_7), x) is a polynomial of degree = max packing number
  IP(G(P_7), 2) = H(P_7) = 189

  KEY: The IP viewpoint gives a COMBINATORIAL decomposition of H
  into contributions from k disjoint cycles (IP_k × 2^k).

  The π connection: IP zeros on the complex plane relate to
  Lee-Yang theory (statistical mechanics) where phase transitions
  occur at roots of the partition function Z(G, λ) = IP(G, λ).
  For H = IP(G, 2), we evaluate PAST the phase transition point.
""")

print("Done!")
