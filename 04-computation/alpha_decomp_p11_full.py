#!/usr/bin/env python3
"""Full alpha decomposition at p=11 including alpha_3, alpha_4, alpha_5."""
import sys, time
sys.path.insert(0, '.')
from alpha_decomp_p19 import build_adj, enumerate_cycles

p = 11
m = 5
S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
S_int = list(range(1, m + 1))

print(f"p={p}, QR={S_qr}, Int={S_int}")
print(f"Known: H(Paley)=95095, H(Interval)=93027")

for name, S in [("Paley", S_qr), ("Interval", S_int)]:
    A = build_adj(p, S)
    all_cycles = []

    for k in range(3, p + 1, 2):
        cycles_k = enumerate_cycles(A, p, k)
        all_cycles.extend(cycles_k)

    n = len(all_cycles)
    print(f"\n{'='*60}")
    print(f"{name}: alpha_1 = {n}")

    # Build compatibility (disjointness) structure
    print(f"  Building compatibility structure...", flush=True)
    t0 = time.time()
    compatible = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if not (all_cycles[i] & all_cycles[j]):
                compatible[i].append(j)
    elapsed = time.time() - t0
    print(f"  Done ({elapsed:.1f}s)", flush=True)

    # alpha_2
    alpha_2 = sum(len(compatible[i]) for i in range(n))
    print(f"  alpha_2 = {alpha_2}")

    # alpha_3: count triples (i < j < k) all mutually disjoint
    print(f"  Computing alpha_3...", flush=True)
    t0 = time.time()
    alpha_3 = 0
    for i in range(n):
        comp_i = compatible[i]  # already sorted, all > i
        for ji, j in enumerate(comp_i):
            # Find k in comp_i with k > j and k disjoint from j
            for k in comp_i[ji+1:]:
                if not (all_cycles[j] & all_cycles[k]):
                    alpha_3 += 1
    elapsed = time.time() - t0
    print(f"  alpha_3 = {alpha_3} ({elapsed:.1f}s)")

    # alpha_4
    print(f"  Computing alpha_4...", flush=True)
    t0 = time.time()
    alpha_4 = 0
    for i in range(n):
        comp_i = compatible[i]
        for ji, j in enumerate(comp_i):
            comp_ij = [k for k in comp_i[ji+1:] if not (all_cycles[j] & all_cycles[k])]
            for ki, k in enumerate(comp_ij):
                for l in comp_ij[ki+1:]:
                    if not (all_cycles[k] & all_cycles[l]):
                        alpha_4 += 1
    elapsed = time.time() - t0
    print(f"  alpha_4 = {alpha_4} ({elapsed:.1f}s)")

    # alpha_5: check if nonzero (need 5 mutually disjoint cycles, using 15+ vertices from 11)
    # This is impossible: 5 cycles of length 3 need 15 vertices but p=11.
    # Even 4 disjoint 3-cycles need 12 > 11 vertices. So alpha_4 of 3-cycles = 0.
    # alpha_4 can come from mixed lengths: e.g., 3+3+3+... nope, 4*3=12>11.
    # 3+3+5 = 11: possible! 3+3+3 = 9 <= 11: possible for alpha_3.
    # 3+3+3+3 = 12 > 11: impossible. So alpha_4 from all 3-cycles = 0.
    # But alpha_4 can include 5-cycles: 3+3+3+5 = 14 > 11. Nope.
    # Even 3+3+5 = 11: that's 3 cycles using 11 vertices. That's alpha_3!
    # alpha_4 would need 4 disjoint cycles using <= 11 vertices.
    # Minimum: 3+3+3+3 = 12 > 11. Impossible!
    # So alpha_4 = 0 and alpha_j = 0 for j >= 4.
    # Wait, but the computation above should verify this.

    # Verify H
    known_H = {"Paley": 95095, "Interval": 93027}[name]
    H = 1 + 2*n + 4*alpha_2 + 8*alpha_3 + 16*alpha_4
    print(f"\n  H = 1 + 2*{n} + 4*{alpha_2} + 8*{alpha_3} + 16*{alpha_4}")
    print(f"    = {H}")
    print(f"  H (known) = {known_H}")
    print(f"  Match: {H == known_H}")

    # Ising decomposition
    print(f"\n  Ising decomposition:")
    print(f"    alpha_0 = 1")
    print(f"    2*alpha_1 = {2*n}")
    print(f"    4*alpha_2 = {4*alpha_2}")
    print(f"    8*alpha_3 = {8*alpha_3}")
    print(f"    16*alpha_4 = {16*alpha_4}")
