#!/usr/bin/env python3
"""
Verify the p=19 counterexample: does S={2,4,6,8,10,12,14,16,18} beat Paley?
Cross-validate the Held-Karp DP against known p=11 values.
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
    """Held-Karp DP for Hamiltonian path count."""
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

# === CROSS-VALIDATION at p=11 ===
print("=== CROSS-VALIDATION: p=11 ===")
p = 11
known = {
    frozenset({1,3,4,5,9}): 95095,   # Paley
    frozenset({1,2,4,5,8}): 93467,
    frozenset({1,2,3,4,5}): 93027,
    frozenset({2,3,4,5,10}): 92411,
}

for S, expected_H in known.items():
    A = adjacency_matrix(S, p)
    H = count_hp_fast(A)
    match = "OK" if H == expected_H else "MISMATCH!!!"
    print(f"  S={sorted(S)}: H={H}, expected={expected_H} {match}")

# Also check even residues at p=11
S_even_11 = frozenset({2,4,6,8,10})
A = adjacency_matrix(S_even_11, p)
H_even_11 = count_hp_fast(A)
print(f"\n  Even residues S={sorted(S_even_11)}: H={H_even_11}")
print(f"  (Paley H = 95095)")

# === MAIN TEST: p=19 ===
print("\n=== p=19 VERIFICATION ===")
p = 19

# Paley
S_paley = frozenset({1,4,5,6,7,9,11,16,17})
A_paley = adjacency_matrix(S_paley, p)
t0 = time.time()
H_paley = count_hp_fast(A_paley)
t1 = time.time()
print(f"Paley S={sorted(S_paley)}: H={H_paley} ({t1-t0:.1f}s)")

# Even residues
S_even = frozenset({2,4,6,8,10,12,14,16,18})
A_even = adjacency_matrix(S_even, p)
t0 = time.time()
H_even = count_hp_fast(A_even)
t1 = time.time()
print(f"Even  S={sorted(S_even)}: H={H_even} ({t1-t0:.1f}s)")

# Verify: is S_even a valid tournament (anti-complete)?
valid = all((p - j) % p not in S_even for j in S_even)
print(f"\nS_even valid (S ∩ -S = empty)? {valid}")
for j in sorted(S_even):
    comp = (p - j) % p
    print(f"  {j} ↔ {comp}: complement in S? {comp in S_even}")

print(f"\n*** RESULT: H_even - H_paley = {H_even - H_paley} ***")
if H_even > H_paley:
    print(f"*** PALEY DOES NOT MAXIMIZE H AT p=19! ***")
    print(f"*** Even residue tournament beats Paley by {H_even - H_paley} ***")
    print(f"*** Ratio: {H_even / H_paley:.6f} ***")
else:
    print(f"*** Paley maximizes H among these two ***")

# Test a few more connection sets to find the true maximum
print("\n=== TESTING MORE CONNECTION SETS AT p=19 ===")
# Cyclic interval: S = {1,2,...,9}
S_interval = frozenset(range(1, 10))
valid = all((p - j) % p not in S_interval for j in S_interval)
if valid:
    A = adjacency_matrix(S_interval, p)
    H = count_hp_fast(A)
    print(f"Interval S={sorted(S_interval)}: H={H}")
else:
    print(f"Interval S={sorted(S_interval)} INVALID")

# Complement of Paley
S_comp = frozenset(range(1,p)) - S_paley
valid = all((p - j) % p not in S_comp for j in S_comp)
if valid:
    A = adjacency_matrix(S_comp, p)
    H = count_hp_fast(A)
    print(f"Paley_comp S={sorted(S_comp)}: H={H}")

# Odd residues
S_odd = frozenset({1,3,5,7,9,11,13,15,17})
valid = all((p - j) % p not in S_odd for j in S_odd)
if valid:
    A = adjacency_matrix(S_odd, p)
    H = count_hp_fast(A)
    print(f"Odd S={sorted(S_odd)}: H={H}")
else:
    print(f"Odd residues INVALID")

# Some random connection sets
import random
random.seed(42)
best_H = max(H_paley, H_even)
best_S = S_even if H_even > H_paley else S_paley

for trial in range(10):
    elems = list(range(1, p))
    random.shuffle(elems)
    S = set()
    for j in elems:
        if j not in S and (p - j) % p not in S:
            S.add(j)
            if len(S) == 9:
                break
    if len(S) != 9:
        continue
    S = frozenset(S)
    A = adjacency_matrix(S, p)
    H = count_hp_fast(A)
    if H > best_H:
        best_H = H
        best_S = S
        print(f"  NEW MAX: S={sorted(S)}: H={H}")
    elif H == best_H:
        print(f"  TIE: S={sorted(S)}: H={H}")
    else:
        print(f"  S={sorted(S)}: H={H}")

print(f"\nBest H found: {best_H} for S={sorted(best_S)}")
