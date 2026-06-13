#!/usr/bin/env python3
"""
Fast Hamiltonian path count for p=23 using layer-by-layer DP.
Confirms the Interval > Paley trend at the next p ≡ 3 mod 4.

opus-2026-03-12-S62d
"""
import numpy as np
import time

def make_tournament(p, S):
    A = np.zeros((p, p), dtype=np.int8)
    for i in range(p):
        for s in S:
            A[i][(i + s) % p] = 1
    return A

def count_H_layered(A, verbose=True):
    """Layer-by-layer DP for Hamiltonian paths."""
    n = len(A)
    full = (1 << n) - 1

    prev = {}
    for i in range(n):
        prev[1 << i] = np.zeros(n, dtype=np.int64)
        prev[1 << i][i] = 1

    for layer in range(2, n + 1):
        t0 = time.time()
        curr = {}
        for mask, counts in prev.items():
            for u in range(n):
                if counts[u] == 0:
                    continue
                for v in range(n):
                    if (mask & (1 << v)) == 0 and A[u][v]:
                        new_mask = mask | (1 << v)
                        if new_mask not in curr:
                            curr[new_mask] = np.zeros(n, dtype=np.int64)
                        curr[new_mask][v] += counts[u]

        t1 = time.time()
        if verbose:
            print(f"    Layer {layer}/{n}: {len(curr)} masks [{t1-t0:.1f}s]", flush=True)

        if layer == n:
            if full in curr:
                return int(curr[full].sum())
            return 0
        prev = curr
    return 0

def get_QR(p):
    return sorted(set(pow(a, 2, p) for a in range(1, p)) - {0})

p = 23
m = (p - 1) // 2
QR = get_QR(p)
S_int = list(range(1, m + 1))

print(f"p = {p}, m = {m}")
print(f"QR = {QR}")
print(f"Interval = {S_int}")
print(f"2^{p} = {2**p:,} states")
print()

results = {}
for name, S in [("Interval", S_int), ("Paley", QR)]:
    A = make_tournament(p, S)
    print(f"Computing H({name})...", flush=True)
    t0 = time.time()
    H = count_H_layered(A)
    t1 = time.time()
    results[name] = H
    print(f"  H({name}) = {H:,} ({t1-t0:.1f}s)")
    print(flush=True)

H_I = results["Interval"]
H_P = results["Paley"]
gap = H_I - H_P
print(f"SUMMARY p={p}:")
print(f"  H(Interval) = {H_I:,}")
print(f"  H(Paley)    = {H_P:,}")
print(f"  Gap: {gap:+,} ({gap/H_P*100:+.4f}%)")
print(f"  Interval {'WINS' if H_I > H_P else 'LOSES'}")
print()

# Quick local max check: just swap boundary elements
print("=" * 72)
print(f"LOCAL MAX CHECK (selected swaps from Interval)")
print("=" * 72)
for a_swap in [1, m//2, m]:
    b_swap = p - a_swap
    S_new = sorted([x for x in S_int if x != a_swap] + [b_swap])
    if len(set(S_new)) == m:
        A_new = make_tournament(p, S_new)
        t0 = time.time()
        H_new = count_H_layered(A_new, verbose=False)
        t1 = time.time()
        delta = H_new - H_I
        print(f"  Swap {a_swap}→{b_swap}: ΔH={delta:+,} ({delta/H_I*100:+.4f}%) [{t1-t0:.1f}s]")

print("\nDONE.")
