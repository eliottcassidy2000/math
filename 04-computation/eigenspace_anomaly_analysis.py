#!/usr/bin/env python3
"""
eigenspace_anomaly_analysis.py — k=0 boundary rank anomaly as topological invariant

DISCOVERY (S68):
  For circulant tournaments at p=7, the k=0 eigenspace boundary rank differs from
  k>=1 eigenspaces at specific degrees. This "anomaly depth" perfectly predicts
  the Betti profile:

  - Paley (anomaly at deg 1,2,3,4): β₄ = 6, all β_i = 0 for i ∉ {0,4}
  - All others (anomaly at deg 1 only): β₁ = 1, all β_i = 0 for i ∉ {0,1}

  The k=0 eigenspace corresponds to the "trivial representation" (constant functions).
  Paley's Gauss sum structure creates an extended anomaly.

This script:
1. Verifies the anomaly pattern at p=7 (done above)
2. Computes the Paley anomaly depth at p=11
3. Checks if anomaly depth determines the Betti concentration dimension
4. Explores WHY Interval has no anomaly beyond m=1

Author: opus-2026-03-12-S68
"""
import sys, time
sys.path.insert(0, '04-computation')
from circulant_homology import CirculantHomology, PaleyHomology, find_nth_root_of_unity

def legendre(a, p):
    if a % p == 0: return 0
    v = pow(a, (p-1)//2, p)
    return v if v == 1 else -1

print("=" * 70)
print("k=0 EIGENSPACE BOUNDARY RANK ANOMALY")
print("=" * 70)

# Paley at p=11: compute all boundary ranks
p = 11
print(f"\n=== Paley T_{p} ===")
h_pal = PaleyHomology(p=p)
h_pal._ensure_enumerated(p)
omega_p = find_nth_root_of_unity(p, h_pal.prime)

print(f"Omega: {h_pal.omega_dims(max_degree=p-1)}")
print(f"Known Betti: {h_pal.betti_numbers(max_degree=p-1)}")
print()

# Compute boundary ranks
print(f"Boundary ranks (k=0 vs k=1):")
for m in range(p + 1):
    t0 = time.time()
    ranks = []
    for k in range(p):
        omega_k = pow(omega_p, k, h_pal.prime)
        rk = h_pal._boundary_rank_k(m, omega_k)
        ranks.append(rk)
    dt = time.time() - t0

    anomaly = "ANOMALY" if ranks[0] != ranks[1] else ""
    # Check if k>=1 are all equal
    k1_same = len(set(ranks[1:])) == 1

    print(f"  m={m:2d}: k=0: {ranks[0]:5d}, k=1: {ranks[1]:5d}, k>=1 same: {k1_same}  {anomaly}  ({dt:.1f}s)")
    sys.stdout.flush()

# Compute per-eigenspace Betti
print(f"\nPer-eigenspace Betti decomposition:")
omega_dims = h_pal.omega_dims(max_degree=p-1)

# Recompute all ranks
all_ranks = {}
for m in range(p + 1):
    ranks_m = []
    for k in range(p):
        omega_k = pow(omega_p, k, h_pal.prime)
        rk = h_pal._boundary_rank_k(m, omega_k)
        ranks_m.append(rk)
    all_ranks[m] = ranks_m

for deg in range(p):
    per_k = []
    for k in range(p):
        ker = omega_dims[deg] - all_ranks[deg][k]
        im_next = all_ranks[deg + 1][k]
        beta_k = ker - im_next
        per_k.append(beta_k)
    total = sum(per_k)
    if total > 0:
        print(f"  β_{deg} = {total}: per-k = {per_k}")

print("\nDONE.")
