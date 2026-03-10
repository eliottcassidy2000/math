# ⚠️ WARNING: This script uses QR mod p for p ≡ 1 (mod 4), which does NOT
# produce a tournament (S ∩ (-S) ≠ ∅ gives bidirectional edges).
# Results for those primes are INVALID. See MISTAKE-011b.
# Valid Paley tournaments require p ≡ 3 (mod 4).

#!/usr/bin/env python3
"""
AUGMENTED TOURNAMENT reformulation of NONHAM=0.

For T[a,b]=0, define T' = T + {a->b} (add the missing edge).
Then M[a,b] = sum over Ham paths pi in T' using edge (a,b),
              weighted by (-1)^{pos(a, pi)}.

For M[a,b]=0: the Ham paths in T' using edge a->b must have
positions of a balanced in alternating sum.

Key question: what does T' look like?
- T is a tournament on V
- T' has both a->b AND b->a (a "double edge")
- T' restricted to V\{a,b} is still a tournament
- T'[a,v] = T[a,v] for v != b
- T'[b,v] = T[b,v] for v != a
- T'[a,b] = T'[b,a] = 1 (bidirectional)

For position-uniform T: every vertex appears at each position
with frequency H/n. Does this constrain the positions of a
in Ham paths of T' using edge (a,b)?

Let consec'(a,b,j) = #{Ham paths in T' with a at position j, b at j+1}.
Then M[a,b] = sum_j (-1)^j consec'(a,b,j).

But T' has both a->b and b->a, so there are also Ham paths with b->a.
We only count those using a->b (where a precedes b consecutively).

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations
import numpy as np

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

# ============================================================
# n=5 Paley: Augmented tournament analysis
# ============================================================
print("=" * 70)
print("n=5 Paley: Augmented tournament T' = T + {a->b}")
print("=" * 70)

n = 5
T = circulant_tournament(n, {1, 4})

for a in range(3):
    for b in range(3):
        if a == b: continue
        if T[(a,b)] == 1: continue  # Only non-edges

        # Create T'
        Tp = dict(T)
        Tp[(a,b)] = 1  # Add edge a->b (b->a already exists)

        # Count Ham paths in T' using edge (a,b)
        # These are paths where a immediately precedes b
        consec = [0] * (n - 1)  # consec[j] = # paths with a at position j, b at j+1
        total_ham_Tp = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                if Tp.get((perm[k], perm[k+1]), 0) == 0:
                    prod = 0; break
            if prod > 0:
                total_ham_Tp += 1
                for j in range(n-1):
                    if perm[j] == a and perm[j+1] == b:
                        consec[j] += 1

        alt_sum = sum((-1)**j * consec[j] for j in range(n-1))

        # Also count Ham paths in T (original)
        total_ham_T = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                if T.get((perm[k], perm[k+1]), 0) == 0:
                    prod = 0; break
            if prod > 0:
                total_ham_T += 1

        print(f"\n  ({a},{b}): T[{a},{b}]=0, T'[{a},{b}]=1, T'[{b},{a}]=1")
        print(f"    H(T)={total_ham_T}, H(T')={total_ham_Tp}")
        print(f"    consec(a,b,j) = {consec}")
        print(f"    sum(-1)^j consec = {alt_sum} (= M[{a},{b}])")

        # The NEW Ham paths in T' (not in T) are those using edge a->b
        new_paths = 0
        new_consec = [0] * (n - 1)
        for perm in permutations(range(n)):
            prod_Tp = 1
            for k in range(n-1):
                if Tp.get((perm[k], perm[k+1]), 0) == 0:
                    prod_Tp = 0; break
            prod_T = 1
            for k in range(n-1):
                if T.get((perm[k], perm[k+1]), 0) == 0:
                    prod_T = 0; break

            if prod_Tp > 0 and prod_T == 0:
                # This is a NEW path (uses edge a->b)
                new_paths += 1
                for j in range(n-1):
                    if perm[j] == a and perm[j+1] == b:
                        new_consec[j] += 1

        new_alt = sum((-1)**j * new_consec[j] for j in range(n-1))
        print(f"    NEW paths (using a->b): {new_paths}")
        print(f"    new_consec = {new_consec}")
        print(f"    sum(-1)^j new_consec = {new_alt}")

        # Paths using ONLY original edges that happen to have a->b consecutively
        # (must use b->a somewhere else in the path)
        orig_consec = [consec[j] - new_consec[j] for j in range(n-1)]
        print(f"    orig_consec (paths in T with a,b consecutive via b->a): {orig_consec}")
        print(f"    sum(-1)^j orig_consec = {sum((-1)**j * orig_consec[j] for j in range(n-1))}")


# ============================================================
# Check: does M[a,b] = signed count of NEW Ham paths only?
# ============================================================
print("\n" + "=" * 70)
print("n=5: Does M[a,b] = signed count of NEW paths in T'?")
print("=" * 70)

for a in range(n):
    for b in range(n):
        if a == b: continue
        if T[(a,b)] == 1: continue

        U = [v for v in range(n) if v != a and v != b]
        M_ab = 0
        for mask in range(1 << len(U)):
            S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            sign = (-1)**len(S_list)
            S_set = sorted(set(S_list) | {a})
            R_set = sorted(set(R) | {b})
            ea = 0
            for p in permutations(S_set):
                if p[-1] != a: continue
                if all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1)):
                    ea += 1
            bb = 0
            for p in permutations(R_set):
                if p[0] != b: continue
                if all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1)):
                    bb += 1
            M_ab += sign * ea * bb

        # Compute new_alt_sum
        Tp = dict(T)
        Tp[(a,b)] = 1
        new_consec = [0] * (n-1)
        for perm in permutations(range(n)):
            prod_Tp = 1
            for k in range(n-1):
                if Tp.get((perm[k], perm[k+1]), 0) == 0:
                    prod_Tp = 0; break
            prod_T = 1
            for k in range(n-1):
                if T.get((perm[k], perm[k+1]), 0) == 0:
                    prod_T = 0; break
            if prod_Tp > 0 and prod_T == 0:
                for j in range(n-1):
                    if perm[j] == a and perm[j+1] == b:
                        new_consec[j] += 1
        new_alt = sum((-1)**j * new_consec[j] for j in range(n-1))

        match = "YES" if M_ab == new_alt else "NO"
        print(f"  ({a},{b}): M[a,b]={M_ab}, new_alt={new_alt}, match={match}")


# ============================================================
# n=7 Paley: same analysis
# ============================================================
print("\n" + "=" * 70)
print("n=7 Paley: Augmented tournament analysis")
print("=" * 70)

n = 7
T = circulant_tournament(n, {1, 2, 4})

# Only check one non-edge pair by circulant symmetry
for b in range(1, n):
    if T[(0,b)] == 1: continue
    a = 0

    Tp = dict(T)
    Tp[(a,b)] = 1

    consec = [0] * (n - 1)
    new_consec = [0] * (n - 1)
    total_T = 0
    total_Tp = 0

    for perm in permutations(range(n)):
        prod_Tp = 1
        for k in range(n-1):
            if Tp.get((perm[k], perm[k+1]), 0) == 0:
                prod_Tp = 0; break
        prod_T = 1
        for k in range(n-1):
            if T.get((perm[k], perm[k+1]), 0) == 0:
                prod_T = 0; break

        if prod_T > 0:
            total_T += 1
        if prod_Tp > 0:
            total_Tp += 1
            for j in range(n-1):
                if perm[j] == a and perm[j+1] == b:
                    consec[j] += 1
                    if prod_T == 0:
                        new_consec[j] += 1

    alt_sum = sum((-1)**j * consec[j] for j in range(n-1))
    new_alt = sum((-1)**j * new_consec[j] for j in range(n-1))

    print(f"  (0,{b}): H(T)={total_T}, H(T')={total_Tp}")
    print(f"    consec = {consec}")
    print(f"    new_consec = {new_consec}")
    print(f"    alt_sum(consec) = {alt_sum}")
    print(f"    alt_sum(new) = {new_alt}")
    break  # One example suffices

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
