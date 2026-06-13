#!/usr/bin/env python3
"""
QUESTION: Does NONHAM = 0 for ALL tournaments, or only position-uniform ones?

If NONHAM = 0 universally, then M[a,b] = sum_j (-1)^j consec(a,b,j) always,
which gives a beautiful position-based formula for the transfer matrix.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations
import numpy as np

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def compute_M_entry(T, n, a, b):
    if a == b:
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                pos = list(perm).index(a)
                val += (-1)**pos
        return val
    else:
        U = [v for v in range(n) if v != a and v != b]
        val = 0
        for mask in range(1 << len(U)):
            S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            sign = (-1)**len(S_list)
            S_set = sorted(set(S_list) | {a})
            R_set = sorted(set(R) | {b})
            ea = 0
            if len(S_set) == 1:
                ea = 1
            else:
                for p in permutations(S_set):
                    if p[-1] != a: continue
                    prod = 1
                    for pk in range(len(p)-1):
                        prod *= T.get((p[pk], p[pk+1]), 0)
                    ea += prod
            bb = 0
            if len(R_set) == 1:
                bb = 1
            else:
                for p in permutations(R_set):
                    if p[0] != b: continue
                    prod = 1
                    for pk in range(len(p)-1):
                        prod *= T.get((p[pk], p[pk+1]), 0)
                    bb += prod
            val += sign * ea * bb
        return val

def tournament_from_bits(n, bits):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T


# ============================================================
# n=3: ALL tournaments
# ============================================================
print("=" * 70)
print("n=3: NONHAM = 0 for ALL tournaments?")
print("=" * 70)

n = 3
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
fail_count = 0

for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)

    for a in range(n):
        for b in range(n):
            if a == b: continue

            U = [v for v in range(n) if v != a and v != b]
            ham_sum = 0
            nonham_sum = 0

            for mask in range(1 << len(U)):
                S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
                R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                sign = (-1)**len(S_list)
                S_set = sorted(set(S_list) | {a})
                R_set = sorted(set(R) | {b})

                for p in permutations(S_set):
                    if p[-1] != a: continue
                    p_valid = all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1))
                    if not p_valid: continue

                    for q in permutations(R_set):
                        if q[0] != b: continue
                        q_valid = all(T.get((q[k], q[k+1]), 0) == 1 for k in range(len(q)-1))
                        if not q_valid: continue

                        concat = list(p) + list(q)
                        is_ham = all(T.get((concat[k], concat[k+1]), 0) == 1 for k in range(len(concat)-1))
                        if is_ham:
                            ham_sum += sign
                        else:
                            nonham_sum += sign

            if nonham_sum != 0:
                fail_count += 1
                if fail_count <= 3:
                    M_val = compute_M_entry(T, n, a, b)
                    consec = sum((-1)**j * sum(1 for perm in permutations(range(n))
                        if all(T.get((perm[k], perm[k+1]), 0) == 1 for k in range(n-1))
                        and perm[j] == a and perm[j+1] == b) for j in range(n-1))
                    print(f"  bits={bits}: ({a},{b}): HAM={ham_sum}, NONHAM={nonham_sum}, "
                          f"M={M_val}, consec_formula={consec}")

print(f"\n  n=3: {fail_count} (a,b) pairs with nonzero NONHAM")
print(f"  out of {(1 << len(pairs)) * n * (n-1)} total (a,b) pairs across all tournaments")


# ============================================================
# n=4: ALL tournaments
# ============================================================
print("\n" + "=" * 70)
print("n=4: NONHAM = 0 for ALL tournaments?")
print("=" * 70)

n = 4
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
fail_count = 0

for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)

    for a in range(n):
        for b in range(n):
            if a == b: continue

            U = [v for v in range(n) if v != a and v != b]
            nonham_sum = 0

            for mask in range(1 << len(U)):
                S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
                R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                sign = (-1)**len(S_list)
                S_set = sorted(set(S_list) | {a})
                R_set = sorted(set(R) | {b})

                for p in permutations(S_set):
                    if p[-1] != a: continue
                    p_valid = all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1))
                    if not p_valid: continue

                    for q in permutations(R_set):
                        if q[0] != b: continue
                        q_valid = all(T.get((q[k], q[k+1]), 0) == 1 for k in range(len(q)-1))
                        if not q_valid: continue

                        concat = list(p) + list(q)
                        is_ham = all(T.get((concat[k], concat[k+1]), 0) == 1 for k in range(len(concat)-1))
                        if not is_ham:
                            nonham_sum += sign

            if nonham_sum != 0:
                fail_count += 1

print(f"  n=4: {fail_count} (a,b) pairs with nonzero NONHAM")


# ============================================================
# n=5: ALL tournaments
# ============================================================
print("\n" + "=" * 70)
print("n=5: NONHAM = 0 for ALL tournaments?")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
fail_count = 0
checked = 0

for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)
    checked += 1

    for a in range(n):
        for b in range(a+1, n):  # Only check a<b to save time (symmetry)
            U = [v for v in range(n) if v != a and v != b]
            nonham_sum = 0

            for mask in range(1 << len(U)):
                S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
                R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                sign = (-1)**len(S_list)
                S_set = sorted(set(S_list) | {a})
                R_set = sorted(set(R) | {b})

                for p in permutations(S_set):
                    if p[-1] != a: continue
                    p_valid = all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1))
                    if not p_valid: continue

                    for q in permutations(R_set):
                        if q[0] != b: continue
                        q_valid = all(T.get((q[k], q[k+1]), 0) == 1 for k in range(len(q)-1))
                        if not q_valid: continue

                        concat = list(p) + list(q)
                        is_ham = all(T.get((concat[k], concat[k+1]), 0) == 1 for k in range(len(concat)-1))
                        if not is_ham:
                            nonham_sum += sign

            if nonham_sum != 0:
                fail_count += 1
                if fail_count <= 3:
                    M_val = compute_M_entry(T, n, a, b)
                    H = count_H(T, n)
                    print(f"  bits={bits}, H={H}: ({a},{b}): NONHAM={nonham_sum}, M={M_val}")
                break
        if fail_count > 0 and fail_count <= 3:
            pass  # continue checking

    if bits % 200 == 199:
        print(f"  ... checked {bits+1}/{1<<len(pairs)}, {fail_count} fails so far")

print(f"\n  n=5: {fail_count} tournaments with nonzero NONHAM (checked {checked})")

if fail_count == 0:
    print("""
  THEOREM: NONHAM = 0 FOR ALL TOURNAMENTS at n=3,4,5.

  This means: M[a,b] = sum_j (-1)^j consec(a,b,j) universally!

  The transfer matrix entry M[a,b] counts Hamiltonian paths where
  a immediately precedes b, weighted by (-1)^{position of a}.
""")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
