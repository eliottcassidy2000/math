#!/usr/bin/env python3
"""Alpha decomposition at p=11 — timing baseline for p=19."""
import sys, time
sys.path.insert(0, '.')
from alpha_decomp_p19 import build_adj, enumerate_cycles, compute_alpha_decomposition

p = 11
m = 5
S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
S_int = list(range(1, m + 1))

print(f"p={p}, QR={S_qr}, Int={S_int}")

for name, S in [("Paley", S_qr), ("Interval", S_int)]:
    A = build_adj(p, S)
    all_cycles = []
    total_t0 = time.time()

    for k in range(3, p + 1, 2):
        cycles_k = enumerate_cycles(A, p, k)
        print(f"  {name} c_{k} = {len(cycles_k)}")
        all_cycles.extend(cycles_k)

    print(f"  {name} alpha_1 = {len(all_cycles)}")
    enum_time = time.time() - total_t0
    print(f"  Cycle enumeration took {enum_time:.1f}s total")

    alphas = compute_alpha_decomposition(all_cycles)

    H = sum(a * (2**j) for j, a in enumerate(alphas))
    print(f"  {name} alphas = {alphas}")
    print(f"  {name} H (partial) = {H}")
    print(f"  {name} H (known) = {95095}")
    print(f"  Remaining = {95095 - H}")
    print()
