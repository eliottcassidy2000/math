#!/usr/bin/env python3
"""Quick test of alpha decomposition at p=7 to verify against known values."""
import sys
sys.path.insert(0, '.')
from alpha_decomp_p19 import build_adj, enumerate_cycles, compute_alpha_decomposition

p = 7
m = 3
S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
S_int = list(range(1, m + 1))

print(f"p={p}, QR={S_qr}, Int={S_int}")

for name, S in [("Paley", S_qr), ("Interval", S_int)]:
    A = build_adj(p, S)
    all_cycles = []
    for k in range(3, p + 1, 2):
        cycles_k = enumerate_cycles(A, p, k)
        print(f"  {name} c_{k} = {len(cycles_k)}")
        all_cycles.extend(cycles_k)

    print(f"  {name} alpha_1 = {len(all_cycles)}")
    alphas = compute_alpha_decomposition(all_cycles)

    H = sum(a * (2**j) for j, a in enumerate(alphas))
    print(f"  {name} alphas = {alphas}")
    print(f"  {name} H = {H}")
    print(f"  Expected: Paley H=189, alpha=[1,80,7]; Interval H=175, alpha=[1,59,14]")
    print()
