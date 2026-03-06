#!/usr/bin/env python3
"""
Correct Schur expansion analysis using hardcoded character tables.
Instance: opus-2026-03-06-S10

Corrects the buggy chi_hook formula. Uses full character tables.
Focus: which Schur coefficients are negative? Are the HOOK coefficients
(for n>=4 where non-hook partitions exist) always positive?
"""

from itertools import permutations
from math import factorial
from fractions import Fraction
from collections import defaultdict, Counter

CHARACTER_TABLES = {
    3: {
        'partitions': [(3,), (2,1), (1,1,1)],
        'classes': [(1,1,1), (2,1), (3,)],
        'table': [
            [1, 1, 1],       # (3) = trivial
            [2, 0, -1],      # (2,1) = standard
            [1, -1, 1],      # (1,1,1) = sign
        ]
    },
    4: {
        'partitions': [(4,), (3,1), (2,2), (2,1,1), (1,1,1,1)],
        'classes': [(1,1,1,1), (2,1,1), (2,2), (3,1), (4,)],
        'table': [
            [1, 1, 1, 1, 1],
            [3, 1, -1, 0, -1],
            [2, 0, 2, -1, 0],
            [3, -1, -1, 0, 1],
            [1, -1, 1, 1, -1],
        ]
    },
    5: {
        'partitions': [(5,), (4,1), (3,2), (3,1,1), (2,2,1), (2,1,1,1), (1,1,1,1,1)],
        'classes': [(1,1,1,1,1), (2,1,1,1), (2,2,1), (3,1,1), (3,2), (4,1), (5,)],
        'table': [
            [1, 1, 1, 1, 1, 1, 1],        # (5)
            [4, 2, 0, 1, -1, 0, -1],       # (4,1)
            [5, 1, 1, -1, 1, -1, 0],       # (3,2)
            [6, 0, -2, 0, 0, 0, 1],        # (3,1,1)
            [5, -1, 1, -1, -1, 1, 0],      # (2,2,1)
            [4, -2, 0, 1, 1, 0, -1],       # (2,1,1,1)
            [1, -1, 1, 1, -1, -1, 1],      # (1,1,1,1,1)
        ]
    },
}

def z_lambda(partition):
    cnt = Counter(partition)
    result = 1
    for k, mk in cnt.items():
        result *= (k ** mk) * factorial(mk)
    return result

def is_hook(partition):
    """A hook partition (k, 1^{n-k}) has at most one part > 1."""
    return sum(1 for p in partition if p > 1) <= 1

def tournament_from_bits(bits, n):
    adj = {}
    for i in range(n):
        for j in range(n):
            if i != j: adj[(i,j)] = 0
    for i in range(n-1):
        adj[(i,i+1)] = 1; adj[(i+1,i)] = 0
    idx = 0
    for gap in range(2, n):
        for i in range(n - gap):
            j = i + gap
            if (bits >> idx) & 1: adj[(i,j)] = 1; adj[(j,i)] = 0
            else: adj[(j,i)] = 1; adj[(i,j)] = 0
            idx += 1
    return adj

def get_cycles(perm):
    n = len(perm); visited = [False]*n; cycles = []
    for start in range(n):
        if visited[start]: continue
        cycle = []; v = start
        while not visited[v]:
            visited[v] = True; cycle.append(v); v = perm[v]
        cycles.append(tuple(cycle))
    return cycles

def is_directed_cycle_in(cycle, adj):
    k = len(cycle)
    if k == 1: return True
    return all(adj.get((cycle[i], cycle[(i+1)%k]), 0) == 1 for i in range(k))

def compute_p_expansion(n, adj):
    adj_op = {(j,i): adj.get((i,j),0) for i in range(n) for j in range(n) if i!=j}
    coeffs = defaultdict(int)
    for perm_tuple in permutations(range(n)):
        perm = list(perm_tuple); cycles = get_cycles(perm)
        phi = 0; valid = True
        for cycle in cycles:
            if len(cycle) == 1: continue
            in_T = is_directed_cycle_in(cycle, adj)
            in_Top = is_directed_cycle_in(cycle, adj_op)
            if not in_T and not in_Top: valid = False; break
            if in_Top and not in_T: phi += len(cycle) - 1
        if not valid: continue
        partition = tuple(sorted([len(c) for c in cycles], reverse=True))
        coeffs[partition] += (-1)**phi
    return dict(coeffs)

def p_to_schur_full(n, p_coeffs):
    ct = CHARACTER_TABLES.get(n)
    if ct is None: return None
    parts = ct['partitions']; classes = ct['classes']; table = ct['table']
    class_idx = {c: i for i, c in enumerate(classes)}
    s_coeffs = {}
    for i, lam in enumerate(parts):
        coeff = Fraction(0)
        for mu, c_mu in p_coeffs.items():
            j = class_idx.get(mu)
            if j is None: continue
            coeff += Fraction(table[i][j]) * Fraction(c_mu) / Fraction(z_lambda(mu))
        s_coeffs[lam] = coeff
    return s_coeffs

def count_ham_paths(n, adj):
    dp = {}
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if adj.get((v, u), 0):
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def count_directed_3cycles(n, adj):
    count = 0
    for i in range(n):
        for j in range(n):
            if j == i: continue
            for k in range(n):
                if k == i or k == j: continue
                if adj.get((i,j),0) and adj.get((j,k),0) and adj.get((k,i),0):
                    if i < j and i < k: count += 1
    return count

def partition_str(p):
    return '(' + ','.join(map(str, p)) + ')'

def main():
    print("=== CORRECT SCHUR EXPANSION ANALYSIS ===\n")

    for n in [3, 4, 5]:
        ct = CHARACTER_TABLES[n]
        parts = ct['partitions']
        hooks = [p for p in parts if is_hook(p)]
        non_hooks = [p for p in parts if not is_hook(p)]

        print(f"\n{'='*80}")
        print(f"n = {n}")
        print(f"Hooks: {[partition_str(h) for h in hooks]}")
        print(f"Non-hooks: {[partition_str(h) for h in non_hooks]}")
        print(f"{'='*80}")

        m = n*(n-1)//2 - (n-1)
        seen = {}

        hook_pos_count = 0
        non_hook_neg_count = 0
        total = 0

        for bits in range(2**m):
            adj = tournament_from_bits(bits, n)
            H = count_ham_paths(n, adj)
            scores = tuple(sorted([sum(adj.get((i,j),0) for j in range(n) if j!=i) for i in range(n)]))
            c3 = count_directed_3cycles(n, adj)
            key = (scores, H, c3)
            if key in seen: continue
            seen[key] = bits; total += 1

            p_coeffs = compute_p_expansion(n, adj)
            s_coeffs = p_to_schur_full(n, p_coeffs)

            hooks_pos = all(s_coeffs[h] >= 0 for h in hooks)
            any_non_hook_neg = any(s_coeffs[h] < 0 for h in non_hooks) if non_hooks else False

            if hooks_pos: hook_pos_count += 1
            if any_non_hook_neg: non_hook_neg_count += 1

            # Print Schur expansion with hook/non-hook annotation
            print(f"\n  bits={bits}, H={H}, c3={c3}, scores={scores}")
            print(f"  Hooks:     ", end='')
            for h in hooks:
                c = s_coeffs[h]
                sign = '+' if c >= 0 else '-'
                print(f"  {sign}{abs(float(c)):.4f}*s{partition_str(h)}", end='')
            print()
            if non_hooks:
                print(f"  Non-hooks: ", end='')
                for h in non_hooks:
                    c = s_coeffs[h]
                    sign = '+' if c >= 0 else '-'
                    print(f"  {sign}{abs(float(c)):.4f}*s{partition_str(h)}", end='')
                print()
            print(f"  Hooks all >=0: {hooks_pos}, Non-hook negative: {any_non_hook_neg}")

            # Check conjugate symmetry: [s_lam] = [s_{lam'}]
            for lam in parts:
                lam_conj = conjugate_partition(lam, n)
                if lam_conj != lam:
                    if s_coeffs[lam] != s_coeffs[lam_conj]:
                        print(f"  !! Conjugate mismatch: [{partition_str(lam)}]={s_coeffs[lam]} != [{partition_str(lam_conj)}]={s_coeffs[lam_conj]}")

        print(f"\n  SUMMARY n={n}: hook-positive {hook_pos_count}/{total}, non-hook-negative {non_hook_neg_count}/{total}")

def conjugate_partition(p, n):
    """Compute conjugate partition."""
    if not p: return p
    max_val = p[0]
    conj = []
    for i in range(1, max_val + 1):
        conj.append(sum(1 for part in p if part >= i))
    return tuple(conj)

if __name__ == '__main__':
    main()
