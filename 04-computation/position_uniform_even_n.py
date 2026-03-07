#!/usr/bin/env python3
"""
Do position-uniform tournaments exist at EVEN n?

At even n, position uniformity requires H/n integer AND P[v,k] = H/n for all v,k.
If no such tournaments exist, our proof chain only applies at odd n (as expected).

Also: analyze the e_s values (E*B sums by subset size) to understand
the NONHAM=0 cancellation pattern at n=5 and n=7.

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
            if T.get((perm[k], perm[k+1]), 0) == 0:
                prod = 0; break
        count += prod
    return count

def position_matrix(T, n):
    P = np.zeros((n, n), dtype=int)
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            if T.get((perm[k], perm[k+1]), 0) == 0:
                prod = 0; break
        if prod > 0:
            for k in range(n):
                P[perm[k], k] += 1
    return P

# ============================================================
# Check even n = 4 and n = 6
# ============================================================
for n in [3, 4, 5, 6]:
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    uniform_count = 0
    h_values = []

    for bits in range(1 << len(pairs)):
        T = tournament_from_bits(n, bits)
        H = count_H(T, n)
        if H % n != 0:
            continue
        P = position_matrix(T, n)
        is_uniform = all(P[v,k] == H // n for v in range(n) for k in range(n))
        if is_uniform:
            uniform_count += 1
            if H not in h_values:
                h_values.append(H)

    print(f"n={n}: {uniform_count} position-uniform tournaments out of {1 << len(pairs)}")
    print(f"  H values: {sorted(h_values)}")
    print(f"  H/n: {[H//n for H in sorted(h_values)]}")
    print()


# ============================================================
# e_s analysis for position-uniform n=5
# ============================================================
print("=" * 70)
print("e_s analysis: E*B sums by subset size")
print("=" * 70)

def E_paths(T, verts, a):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == a else 0
    count = 0
    for p in permutations(verts):
        if p[-1] != a: continue
        if all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1)):
            count += 1
    return count

def B_paths(T, verts, b):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == b else 0
    count = 0
    for p in permutations(verts):
        if p[0] != b: continue
        if all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1)):
            count += 1
    return count

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

# Find all position-uniform and a non-uniform for comparison
uniform_examples = []
nonuniform_example = None

for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    P = position_matrix(T, n)
    is_uniform = H % n == 0 and all(P[v,k] == H // n for v in range(n) for k in range(n))

    if is_uniform and len(uniform_examples) < 3:
        uniform_examples.append((bits, T, H))
    if not is_uniform and nonuniform_example is None and H > 1:
        nonuniform_example = (bits, T, H)

for label, (bits, T, H) in [("Uniform", uniform_examples[0])] + ([("Non-uniform", nonuniform_example)] if nonuniform_example else []):
    print(f"\n  {label}: bits={bits}, H={H}")

    for a in range(min(n, 3)):
        for b in range(min(n, 3)):
            if a == b: continue
            if T[(a,b)] == 1: continue  # Only non-edges

            U = [v for v in range(n) if v != a and v != b]
            e_s = [0] * (len(U) + 1)

            for mask in range(1 << len(U)):
                S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
                R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                s = len(S_list)
                S_set = sorted(set(S_list) | {a})
                R_set = sorted(set(R) | {b})

                ea = E_paths(T, S_set, a)
                bb = B_paths(T, R_set, b)
                e_s[s] += ea * bb

            alt_sum = sum((-1)**s * e_s[s] for s in range(len(e_s)))
            print(f"    ({a},{b}) T[{a},{b}]=0: e_s={e_s}, alt_sum={alt_sum}")

            # Check: e_s + e_{|U|-s} pattern for odd |U|
            U_len = len(U)
            if U_len % 2 == 1:
                pairings = []
                for s in range((U_len + 1) // 2):
                    s2 = U_len - s
                    if s == s2:
                        pairings.append(f"  e_{s}={e_s[s]} (unpaired)")
                    else:
                        pairings.append(f"  e_{s}={e_s[s]}, e_{s2}={e_s[s2]}, diff={e_s[s]-e_s[s2]}")
                print(f"      Pairing: {'; '.join(pairings)}")


# ============================================================
# n=7 circulant: e_s analysis
# ============================================================
print("\n" + "=" * 70)
print("n=7 circulant: e_s analysis")
print("=" * 70)

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

n = 7
# Paley T_7: QR mod 7 = {1,2,4}
T = circulant_tournament(n, {1, 2, 4})

for b in range(1, 4):
    if T[(0,b)] == 1: continue  # Only non-edges
    a = 0
    U = [v for v in range(n) if v != a and v != b]
    e_s = [0] * (len(U) + 1)

    for mask in range(1 << len(U)):
        S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
        R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
        s = len(S_list)
        S_set = sorted(set(S_list) | {a})
        R_set = sorted(set(R) | {b})

        ea = E_paths(T, S_set, a)
        bb = B_paths(T, R_set, b)
        e_s[s] += ea * bb

    alt_sum = sum((-1)**s * e_s[s] for s in range(len(e_s)))
    print(f"  (0,{b}) T[0,{b}]=0: e_s={e_s}, alt_sum={alt_sum}")

    # Check complement pairing: e_s vs e_{|U|-s}
    U_len = len(U)
    for s in range((U_len + 1) // 2):
        s2 = U_len - s
        if s != s2:
            print(f"    e_{s}={e_s[s]}, e_{s2}={e_s[s2]}, diff={e_s[s]-e_s[s2]}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
