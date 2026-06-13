#!/usr/bin/env python3
"""
N_maximization_paley.py -- kind-pasteur-2026-03-13-S60

Does Paley ALWAYS maximize the total directed odd cycle count N
among all circulant tournaments on Z_p?

We know:
- c3 = p(p-1)(p+1)/24 is constant (depends only on regularity)
- c_p = H / (p-1)! * something... actually c_p = H/p (Hamiltonian cycles / p)
  No: c_p is the number of directed Hamiltonian cycles.
  For circulant T: H = c_p (number of directed Ham cycles through all p vertices)

So N = c3 + c5 + c7 + ... + c_p.
Paley maximizes EACH c_k? Or just the total?

Let's check per-k maximization.
"""

from itertools import combinations
from collections import defaultdict

def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A

def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    start = 0
    dp = {(1 << start, start): 1}
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


for p in [7, 11]:
    m = (p - 1) // 2
    print(f"\n{'='*70}")
    print(f"  PER-LENGTH CYCLE MAXIMIZATION AT p={p}")
    print(f"{'='*70}")

    all_ck = defaultdict(lambda: defaultdict(int))

    for bits in range(1 << m):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)

        for k in range(3, p+1, 2):
            total = 0
            for subset in combinations(range(p), k):
                total += count_ham_cycles(A, list(subset))
            all_ck[bits][k] = total

    # For each k, find who maximizes
    S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    paley_bits = None
    for bits in range(1 << m):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))
        if sorted(S) == S_qr:
            paley_bits = bits
            break

    print(f"  Paley bits = {paley_bits}")

    for k in range(3, p+1, 2):
        max_val = max(all_ck[bits][k] for bits in range(1 << m))
        paley_val = all_ck[paley_bits][k]
        maximizers = [bits for bits in range(1 << m) if all_ck[bits][k] == max_val]
        paley_is_max = paley_val == max_val

        min_val = min(all_ck[bits][k] for bits in range(1 << m))
        distinct_vals = len(set(all_ck[bits][k] for bits in range(1 << m)))

        print(f"\n  k={k}: Paley c_k={paley_val}, max c_k={max_val}, min={min_val}")
        print(f"    Distinct values: {distinct_vals}")
        print(f"    Paley maximizes? {paley_is_max}")
        if not paley_is_max:
            print(f"    Paley rank: {sorted(set(all_ck[bits][k] for bits in range(1 << m)), reverse=True).index(paley_val) + 1}")
        print(f"    # maximizers: {len(maximizers)}")

        # Show all distinct values
        vals = sorted(set(all_ck[bits][k] for bits in range(1 << m)), reverse=True)
        counts = [sum(1 for bits in range(1 << m) if all_ck[bits][k] == v) for v in vals]
        for v, c in zip(vals, counts):
            paley_mark = " <-- PALEY" if v == paley_val else ""
            print(f"      c_{k} = {v}: {c} orientations{paley_mark}")

    # Savchenko's theorem: For DRTs, c_k(Paley) > c_k(any other DRT)
    # But we're comparing across ALL circulant orientations, not just DRTs
    # (all circulant on Z_p are regular = doubly regular)

    # Wait: ALL circulant tournaments on Z_p are regular (each vertex has
    # out-degree (p-1)/2), hence doubly regular.
    # So Savchenko applies: Paley should maximize EVERY c_k.
    # Let's check this!

    print(f"\n  Savchenko check: Does Paley maximize EVERY c_k?")
    all_max = True
    for k in range(3, p+1, 2):
        max_val = max(all_ck[bits][k] for bits in range(1 << m))
        paley_val = all_ck[paley_bits][k]
        if paley_val != max_val:
            all_max = False
            print(f"    FAILS at k={k}: Paley={paley_val}, max={max_val}")
    if all_max:
        print(f"    YES -- Paley maximizes c_k for ALL k in {{3,5,...,{p}}}")

    # Also check: is the anti-Paley (complement) ALSO a maximizer?
    anti_paley = (1 << m) - 1 - paley_bits
    anti_bits_actual = None
    for bits in range(1 << m):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))
        S_set = set(S)
        comp = set(range(1, p)) - S_set
        if len(comp) == m:
            # This is the reversal of bits
            pass

    # Simpler: anti-Paley has S = NQR
    S_nqr = sorted(set(range(1, p)) - set(S_qr))
    for bits in range(1 << m):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))
        if sorted(S) == S_nqr:
            anti_paley_bits = bits
            break

    print(f"\n  Anti-Paley (S=NQR={S_nqr}): bits={anti_paley_bits}")
    for k in range(3, p+1, 2):
        ap_val = all_ck[anti_paley_bits][k]
        p_val = all_ck[paley_bits][k]
        print(f"    k={k}: anti-Paley c_k={ap_val}, Paley c_k={p_val}, same={ap_val==p_val}")

print("\nDONE.")
