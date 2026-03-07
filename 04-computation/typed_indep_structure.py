#!/usr/bin/env python3
"""
Explore the typed independence polynomial I_typed(Omega(T); y_3, y_5, y_7, ...)
and its relationship to the trivariate GF G_T(t,x).

Key facts:
- I(Omega, x) = I_typed(x, x, x, ...) [specialization]
- U_T = sum_{S indep} 2^|S| * prod p_{len(c)} * p_1^{n-sum_len} [GS-OCF bridge]
- G_T(0, x) = I(Omega, x) [THM-063]

Question: Can we extend G_T to separate cycle types?
  G_T(t, y_3, y_5, ...) = A_n(t) + sum_{I: indep sets} prod_{c in I} y_{len(c)} * g_I(t)

At t=0, this gives I_typed. At y_k = x for all k, gives G_T(t,x).

Let's compute I_typed at n=5 and n=7 and look for structure.

opus-2026-03-07-S36
"""
from itertools import permutations, combinations
import random
from collections import defaultdict

def tournament_from_seed(n, seed):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def find_directed_odd_cycles(T, n, vertices=None):
    if vertices is None:
        vertices = list(range(n))
    cycles = set()
    for length in range(3, len(vertices)+1, 2):
        for subset in combinations(vertices, length):
            first = subset[0]
            for perm in permutations(subset[1:]):
                cycle = (first,) + perm
                if all(T[cycle[i]][cycle[(i+1) % length]] for i in range(length)):
                    min_idx = cycle.index(min(cycle))
                    rotated = cycle[min_idx:] + cycle[:min_idx]
                    cycles.add(rotated)
    return list(cycles)

def main():
    # n=5: cycles are length 3 and 5
    n = 5
    m = n*(n-1)//2
    print(f"=== Typed Independence Polynomial at n={n} ===\n")

    # Collect all possible (alpha_3, alpha_5, alpha_33) patterns
    # alpha_3 = # indep sets of size 1 with a 3-cycle
    # alpha_5 = # indep sets of size 1 with a 5-cycle
    # alpha_33 = # indep sets of size 2 with two 3-cycles (must be disjoint)

    patterns = defaultdict(int)

    for bits_int in range(2**m):
        b = [(bits_int >> k) & 1 for k in range(m)]
        A = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if b[idx] == 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        cycles = find_directed_odd_cycles(A, n)
        cvsets = [frozenset(c) for c in cycles]

        c3 = [i for i in range(len(cycles)) if len(cycles[i]) == 3]
        c5 = [i for i in range(len(cycles)) if len(cycles[i]) == 5]

        # Count disjoint 3-cycle pairs
        disjoint_33 = 0
        for a in range(len(c3)):
            for b_idx in range(a+1, len(c3)):
                if not (cvsets[c3[a]] & cvsets[c3[b_idx]]):
                    disjoint_33 += 1

        pattern = (len(c3), len(c5), disjoint_33)
        patterns[pattern] += 1

    print("Pattern (c3, c5, disjoint_33): count")
    for p in sorted(patterns.keys()):
        c3, c5, d33 = p
        H = 1 + 2*(c3 + c5) + 4*d33
        print(f"  ({c3}, {c5}, {d33}): {patterns[p]:4d} tournaments, H={H}")

    print(f"\nI_typed(y3, y5) = 1 + c3*y3 + c5*y5 + d33*y3^2")
    print(f"  (no y3*y5 term: 3-cycle and 5-cycle always share vertex at n=5)")
    print(f"  (no y5^2 term: at most one 5-cycle can be independent at n=5)")

    # Now n=7: more interesting structure
    print(f"\n{'='*60}")
    n = 7
    print(f"=== Typed Independence at n={n} (random samples) ===\n")

    for seed in range(5):
        A = tournament_from_seed(n, seed)
        cycles = find_directed_odd_cycles(A, n)
        cvsets = [frozenset(c) for c in cycles]
        nc = len(cycles)

        c3 = [i for i in range(nc) if len(cycles[i]) == 3]
        c5 = [i for i in range(nc) if len(cycles[i]) == 5]
        c7 = [i for i in range(nc) if len(cycles[i]) == 7]

        # Build adjacency
        adj = {}
        for i in range(nc):
            adj[i] = set()
            for j in range(nc):
                if i != j and cvsets[i] & cvsets[j]:
                    adj[i].add(j)

        # Count typed independent sets
        # Type (k): single cycle of length k
        # Type (k1, k2): pair of non-adjacent cycles of lengths k1, k2
        # Type (k1, k2, k3): triple, etc.

        typed_counts = defaultdict(int)

        # Size 1
        for i in range(nc):
            typed_counts[(len(cycles[i]),)] += 1

        # Size 2
        for i in range(nc):
            for j in range(i+1, nc):
                if j not in adj[i]:
                    key = tuple(sorted([len(cycles[i]), len(cycles[j])]))
                    typed_counts[key] += 1

        # Size 3
        for i in range(nc):
            for j in range(i+1, nc):
                if j in adj[i]:
                    continue
                for k in range(j+1, nc):
                    if k in adj[i] or k in adj[j]:
                        continue
                    key = tuple(sorted([len(cycles[i]), len(cycles[j]), len(cycles[k])]))
                    typed_counts[key] += 1

        # Total H check
        total_alpha = defaultdict(int)
        for key, cnt in typed_counts.items():
            total_alpha[len(key)] += cnt

        H_check = 1 + sum(2**k * total_alpha[k] for k in total_alpha)
        H_direct = 0
        for mask in range(2**nc):
            selected = [i for i in range(nc) if (mask >> i) & 1]
            ok = True
            for idx1 in range(len(selected)):
                for idx2 in range(idx1+1, len(selected)):
                    if selected[idx2] in adj[selected[idx1]]:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                H_direct += 2**len(selected)

        print(f"seed={seed}: #cycles = {len(c3)}x3 + {len(c5)}x5 + {len(c7)}x7 = {nc} total")
        print(f"  H = {H_direct}")
        print(f"  Typed independent sets:")
        for key in sorted(typed_counts.keys()):
            print(f"    type {key}: {typed_counts[key]}")

        # Verify I_typed specialization
        H_from_typed = 1 + sum(2**len(key) * typed_counts[key] for key in typed_counts)
        print(f"  H from typed: {H_from_typed} (check: {H_from_typed == H_direct})")

        # S = sum(l_i - 1) for each type
        print(f"  S-values:")
        for key in sorted(typed_counts.keys()):
            S = sum(l - 1 for l in key)
            f = len(key)
            print(f"    type {key}: S={S}, f={f}, d=S+f={S+f}")
        print()

if __name__ == "__main__":
    main()
