#!/usr/bin/env python3
"""
interval_betti_p11.py — Memory-efficient Betti computation for Interval T_11.

Computes one boundary rank at a time, freeing memory between degrees.

Author: opus-2026-03-12-S68
"""
import sys, gc, time
sys.path.insert(0, '04-computation')
from circulant_homology import CirculantHomology

p = 11
S_int = {1, 2, 3, 4, 5}
S_pal = set((a*a) % p for a in range(1, p))

print(f"Interval S = {sorted(S_int)}")
print(f"Paley    S = {sorted(S_pal)}")
print()

# Known Paley
print("Paley T_11 (cached): β = [1, 0, 0, 0, 0, 5, 15, 0, 0, 0, 0], χ = 11")
print()

# Compute Interval Betti degree-by-degree
h = CirculantHomology(n=p, S=S_int)

# Get Omega dims (cheap)
omega = h.omega_dims(max_degree=p-1, verbose=True)
print(f"\nInterval Omega: {omega}")
chi_omega = sum((-1)**m * d for m, d in enumerate(omega))
print(f"χ(Omega) = {chi_omega}")

# Now compute boundary ranks for each eigenspace and degree
from circulant_homology import find_nth_root_of_unity
omega_p = find_nth_root_of_unity(p, h.prime)

print(f"\nComputing boundary ranks (prime={h.prime})...")

# boundary_ranks[k][m] = rank(d_m restricted to Ω_m^(k))
all_ranks = {}
for m in range(p + 1):
    t0 = time.time()
    ranks_m = []
    for k in range(p):
        omega_k = pow(omega_p, k, h.prime)
        rk = h._boundary_rank_k(m, omega_k)
        ranks_m.append(rk)
    dt = time.time() - t0
    all_ranks[m] = ranks_m
    print(f"  m={m:2d}: ranks = {ranks_m}  ({dt:.1f}s)")
    sys.stdout.flush()
    gc.collect()

# Compute Betti numbers
print(f"\nBetti numbers:")
betti = []
for m in range(p):
    total = 0
    for k in range(p):
        ker = omega[m] - all_ranks[m][k]
        im_next = all_ranks[m+1][k]
        total += ker - im_next
    betti.append(total)

print(f"  β = {betti}")
chi = sum((-1)**i * b for i, b in enumerate(betti))
print(f"  χ = {chi}")

print(f"\nComparison:")
print(f"  Paley:    β = [1, 0, 0, 0, 0, 5, 15, 0, 0, 0, 0], χ = 11")
print(f"  Interval: β = {betti}, χ = {chi}")
print("DONE.")
