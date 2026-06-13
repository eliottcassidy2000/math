#!/usr/bin/env python3
"""
Quick check: does R-minimization hold at n=7?
n=7 has 2^21 = 2097152 tournaments. We compute H and sum H(T-v) for each.

kind-pasteur-2026-03-06-S18g
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count

MAX_H = {7: 189, 6: 45}

def delete_vertex(T, v):
    n = len(T)
    verts = [i for i in range(n) if i != v]
    return [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]

n = 7
m = n * (n - 1) // 2  # 21

min_r = float('inf')
min_r_h = 0
max_h = 0
count = 0
batch = 100000

print(f"Checking R-minimization at n={n} (2^{m} = {1 << m} tournaments)")

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    if h == 0:
        continue

    if h > max_h:
        max_h = h

    ds = sum(hamiltonian_path_count(delete_vertex(T, v)) for v in range(n))
    r = ds / h

    if r < min_r:
        min_r = r
        min_r_h = h

    count += 1
    if count % batch == 0:
        print(f"  {count} checked, min R so far = {min_r:.4f} (at H={min_r_h}), max H = {max_h}")

print(f"\nFinal: {count} tournaments with H>0")
print(f"Min R = {min_r:.6f}, achieved at H = {min_r_h}")
print(f"Max H = {max_h}")
print(f"R-minimizer = H-maximizer? {min_r_h == max_h}")
print(f"Expected: max H = {MAX_H[7]}, R = {7 * MAX_H[6] / MAX_H[7]:.6f} = 5/3")

print("\nDone.")
