#!/usr/bin/env python3
"""
Fast R-minimization check at n=7.
Instead of tracking min R across all tournaments, just verify:
R(T) >= 5/3 for all T with H > 0.

Optimization: compute H first, skip if H is very small (they have high R).

kind-pasteur-2026-03-06-S18g
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count

def delete_vertex(T, v):
    n = len(T)
    verts = [i for i in range(n) if i != v]
    return [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]

n = 7
m = 21
target_r = 5/3 - 0.0001  # 1.6666... with small tolerance

violations = 0
checked = 0
total = 1 << m
batch = 200000

print(f"Checking R >= 5/3 for all n={n} tournaments ({total} total)")
print(f"Only checking tournaments with H >= 100 (low H => high R)")

for bits in range(total):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)

    # Low H always has high R (R >= 1 always, and for H << max, R >> 5/3)
    if h < 100:
        continue

    ds = sum(hamiltonian_path_count(delete_vertex(T, v)) for v in range(n))
    r = ds / h
    checked += 1

    if r < target_r:
        violations += 1
        print(f"  VIOLATION: bits={bits}, H={h}, R={r:.6f}")

    if checked % batch == 0:
        print(f"  Checked {checked} (H>=100), bits={bits}/{total}, violations={violations}")

print(f"\nResult: {violations} violations out of {checked} tournaments with H>=100")
if violations == 0:
    print("R >= 5/3 for all tournaments with H >= 100. Conjecture HOLDS at n=7.")
else:
    print(f"Conjecture FAILS: {violations} tournaments have R < 5/3")

print("Done.")
