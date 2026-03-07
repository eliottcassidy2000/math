#!/usr/bin/env python3
"""
PROOF ANALYSIS: Why NONHAM=0 for position-uniform tournaments?

Key observation from nonham_vanish_general.py + nonham_vanish_uniform.py:

NONHAM(a,b) = 0 trivially when T[a,b]=1 (all split pairs are Ham paths).
NONHAM(a,b) = M[a,b] when T[a,b]=0 (all split pairs fail at junction).

So NONHAM=0 for all (a,b) <=> M[a,b]=0 whenever T[a,b]=0.

PROOF CHAIN (verified n=3,5,7):
1. For position-uniform T: M[a,b]=0 whenever T[a,b]=0  [NONHAM=0]
2. THM-030: M[a,b]=M[b,a] for all tournaments
3. For a!=b: exactly one of T[a,b], T[b,a] is 0
4. If T[a,b]=1: T[b,a]=0 => M[b,a]=0 by (1) => M[a,b]=0 by (2)
5. If T[a,b]=0: M[a,b]=0 directly by (1)
6. So M[a,b]=0 for ALL a!=b
7. M[a,a] = H/n from position uniformity
8. M = (H/n)*I  QED

The only unproved step is (1): WHY does M[a,b]=0 when T[a,b]=0
for position-uniform tournaments?

This script investigates the STRUCTURE of M[a,b] when T[a,b]=0
to understand the cancellation mechanism.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations
import numpy as np

def tournament_from_bits(n, bits):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def position_matrix(T, n):
    P = np.zeros((n, n), dtype=int)
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            for k in range(n):
                P[perm[k], k] += 1
    return P

def compute_M_full(T, n):
    """Compute full transfer matrix M via inclusion-exclusion."""
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        # Diagonal: sum over Ham paths of (-1)^{pos(a)}
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                pos = list(perm).index(a)
                val += (-1)**pos
        M[a,a] = val

        for b in range(n):
            if b == a: continue
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
            M[a,b] = val
    return M


# ============================================================
# Detailed analysis: M[a,b] subset-by-subset when T[a,b]=0
# ============================================================
print("=" * 70)
print("n=5: Subset-by-subset M[a,b] contributions when T[a,b]=0")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

# Find a position-uniform tournament
for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    P = position_matrix(T, n)
    is_uniform = H % n == 0 and all(P[v,k] == H // n for v in range(n) for k in range(n))
    if not is_uniform or H != 15:  # Want Paley-like
        continue

    # Find (a,b) with T[a,b]=0
    for a in range(n):
        for b in range(n):
            if a == b: continue
            if T.get((a,b),0) == 1: continue  # Want T[a,b]=0

            U = [v for v in range(n) if v != a and v != b]
            print(f"\n  bits={bits}, H={H}: ({a},{b}), T[{a},{b}]=0")
            print(f"  U = {U}")

            total = 0
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

                contrib = sign * ea * bb
                total += contrib
                if ea > 0 and bb > 0:
                    print(f"    S={sorted(S_list)}, |S|={len(S_list)}, sign={sign:+d}: "
                          f"E_a({S_set})={ea}, B_b({R_set})={bb}, contrib={contrib:+d}")

            print(f"    TOTAL M[{a},{b}] = {total}")

            # Only show first 2 examples
            if a + b > 2:
                break
        if a + b > 2:
            break
    break


# ============================================================
# Compare: same analysis for NON-uniform tournament with T[a,b]=0
# ============================================================
print("\n" + "=" * 70)
print("n=5: Same analysis for NON-uniform tournament")
print("=" * 70)

for bits in [0]:  # Transitive tournament
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    P = position_matrix(T, n)

    for a in range(n):
        for b in range(n):
            if a == b: continue
            if T.get((a,b),0) == 1: continue

            U = [v for v in range(n) if v != a and v != b]
            print(f"\n  bits={bits}, H={H}: ({a},{b}), T[{a},{b}]=0")
            print(f"  U = {U}")

            total = 0
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

                contrib = sign * ea * bb
                total += contrib
                if ea > 0 and bb > 0:
                    print(f"    S={sorted(S_list)}, |S|={len(S_list)}, sign={sign:+d}: "
                          f"E_a({S_set})={ea}, B_b({R_set})={bb}, contrib={contrib:+d}")

            print(f"    TOTAL M[{a},{b}] = {total}")
            if a == 1:
                break
        if a == 1:
            break


# ============================================================
# E_a and B_b by subset size for uniform vs non-uniform
# ============================================================
print("\n" + "=" * 70)
print("n=5: E_a and B_b totals by subset size |S|")
print("=" * 70)

for label, bits_val in [("Paley T_5 (uniform)", None), ("Transitive (non-uniform)", 0)]:
    if bits_val is None:
        # Find Paley
        for bits in range(1 << len(pairs)):
            T = tournament_from_bits(n, bits)
            H = count_H(T, n)
            P = position_matrix(T, n)
            is_uniform = H % n == 0 and all(P[v,k] == H // n for v in range(n) for k in range(n))
            if is_uniform and H == 15:
                bits_val = bits
                break
    T = tournament_from_bits(n, bits_val)
    H = count_H(T, n)

    print(f"\n  {label}: bits={bits_val}, H={H}")

    for a in range(min(n, 3)):
        for b in range(min(n, 3)):
            if a == b: continue
            if T.get((a,b),0) == 1: continue

            U = [v for v in range(n) if v != a and v != b]
            print(f"    ({a},{b}) T[{a},{b}]=0:")

            by_size = {}
            for mask in range(1 << len(U)):
                S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
                s = len(S_list)
                R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
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
                bb2 = 0
                if len(R_set) == 1:
                    bb2 = 1
                else:
                    for p in permutations(R_set):
                        if p[0] != b: continue
                        prod = 1
                        for pk in range(len(p)-1):
                            prod *= T.get((p[pk], p[pk+1]), 0)
                        bb2 += prod

                if s not in by_size:
                    by_size[s] = []
                by_size[s].append((sorted(S_list), ea, bb2, ea * bb2))

            for s in sorted(by_size):
                entries = by_size[s]
                sign = (-1)**s
                total_contrib = sign * sum(e * b for _, e, b, _ in entries)
                detail = ", ".join(f"E={e}*B={b}" for _, e, b, _ in entries)
                sum_eb = sum(e * b for _, e, b, _ in entries)
                print(f"      |S|={s}: sign={sign:+d}, sum(E*B)={sum_eb}, contrib={total_contrib:+d}  [{detail}]")

            total = sum((-1)**s * sum(e * b for _, e, b, _ in entries) for s, entries in by_size.items())
            print(f"      TOTAL = {total}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
