#!/usr/bin/env python3
"""
Find the crossover prime where the interval tournament first beats Paley.
Check p = 3, 7, 11, 13, 17, 19 (all primes where we can do DP).

At p=13: QR_13 = {1,3,4,9,10,12} but p=13 is 1 mod 4.
For p not 3 mod 4, Paley is NOT a tournament (since -1 is a QR,
S and -S overlap). But we can still compare interval vs other circulants.

Focus on primes p = 3 mod 4: 3, 7, 11, 19, 23, ...
"""

import numpy as np
import time

def adjacency_matrix(S, p):
    A = np.zeros((p,p), dtype=np.int8)
    for i in range(p):
        for j in range(p):
            if i!=j and (j-i)%p in S:
                A[i][j] = 1
    return A

def count_hp_fast(A):
    n = len(A)
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp[mask, v]
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v, u]:
                    dp[mask | (1 << u), u] += c
    return int(np.sum(dp[full]))

def is_qr(a, p):
    if a % p == 0: return False
    return pow(a, (p-1)//2, p) == 1

# Also check non-Paley primes for the interval vs Paley-like comparison
print("CROSSOVER PRIME SEARCH")
print("=" * 60)

for p in [3, 5, 7, 11, 13, 17, 19]:
    m = (p-1)//2
    paley_mod = p % 4

    # Interval tournament
    S_interval = frozenset(range(1, m+1))
    A_int = adjacency_matrix(S_interval, p)
    H_interval = count_hp_fast(A_int)

    # Paley tournament (only valid for p = 3 mod 4)
    if paley_mod == 3:
        S_paley = frozenset(j for j in range(1,p) if is_qr(j,p))
        # Verify it's a valid tournament
        valid = all((p-j)%p not in S_paley for j in S_paley) and len(S_paley) == m
        if valid:
            A_pal = adjacency_matrix(S_paley, p)
            H_paley = count_hp_fast(A_pal)
        else:
            H_paley = None
    else:
        H_paley = None

    # Interval complement
    S_comp = frozenset(range(m+1, p))
    A_comp = adjacency_matrix(S_comp, p)
    H_comp = count_hp_fast(A_comp)

    print(f"\np = {p} (p mod 4 = {paley_mod}), m = {m}")
    print(f"  Interval S={{1,...,{m}}}: H = {H_interval}")
    print(f"  Interval complement S={{{m+1},...,{p-1}}}: H = {H_comp}")
    if H_paley is not None:
        diff = H_interval - H_paley
        print(f"  Paley QR: H = {H_paley}")
        print(f"  H(interval) - H(Paley) = {diff:+d}")
        if diff > 0:
            print(f"  *** INTERVAL BEATS PALEY by {diff/H_paley*100:.2f}% ***")
        elif diff == 0:
            print(f"  TIE")
        else:
            print(f"  Paley wins by {-diff/H_paley*100:.2f}%")

# Focus on p=3 mod 4 primes
print(f"\n{'='*60}")
print(f"SUMMARY: p = 3 mod 4 primes")
print(f"{'='*60}")
print(f"{'p':>4} {'H(Paley)':>20} {'H(Interval)':>20} {'Diff':>12} {'Winner':>10}")

for p in [3, 7, 11, 19]:
    m = (p-1)//2
    S_int = frozenset(range(1, m+1))
    S_pal = frozenset(j for j in range(1,p) if is_qr(j,p))

    A_int = adjacency_matrix(S_int, p)
    A_pal = adjacency_matrix(S_pal, p)

    t0 = time.time()
    H_int = count_hp_fast(A_int)
    H_pal = count_hp_fast(A_pal)
    elapsed = time.time() - t0

    diff = H_int - H_pal
    winner = "INTERVAL" if diff > 0 else ("TIE" if diff == 0 else "PALEY")
    print(f"{p:>4} {H_pal:>20} {H_int:>20} {diff:>+12} {winner:>10}")

# Check the non-standard primes too (p = 1 mod 4)
# For these, the QR set is NOT a tournament (since -1 is a QR).
# But we can still check how the interval competes.
print(f"\n{'='*60}")
print(f"NON-PALEY PRIMES (p = 1 mod 4)")
print(f"{'='*60}")
for p in [5, 13, 17]:
    m = (p-1)//2
    S_int = frozenset(range(1, m+1))
    S_comp = frozenset(range(m+1, p))

    A_int = adjacency_matrix(S_int, p)
    A_comp = adjacency_matrix(S_comp, p)

    H_int = count_hp_fast(A_int)
    H_comp = count_hp_fast(A_comp)

    print(f"\n  p={p}: H(interval) = {H_int}, H(complement) = {H_comp}")

    # The interval and its complement are related by T -> T^op
    # For p = 1 mod 4, the interval and complement are related by j -> p-j
    # which maps S to -S. Since -1 is a QR for p=1 mod 4, this gives an
    # AUTOMORPHISM of the tournament? No: -S = {p-1, p-2, ..., p-m} = comp.
    # So H(interval) = H(complement)? Check:
    if H_int == H_comp:
        print(f"  H(S) = H(comp) = {H_int} (as expected by symmetry)")
    else:
        print(f"  H(S) != H(comp)!")
