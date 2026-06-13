#!/usr/bin/env python3
"""
Hook Schur positivity test at n=7.
Instance: opus-2026-03-06-S10

Uses the corrected MN hook character formula from hook_positivity_n6.py.
"""

import sys
sys.path.insert(0, '04-computation')
from hook_positivity_n6 import (
    chi_hook_mn, tournament_from_bits, get_cycles, is_directed_cycle_in,
    count_ham_paths, generate_partitions, z_lambda, hook_partitions, partition_str
)
from itertools import permutations
from fractions import Fraction
from collections import defaultdict, Counter

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

def main():
    n = 7
    print(f"=== HOOK POSITIVITY TEST AT n={n} ===\n")

    hooks = hook_partitions(n)
    hook_names = [f"s{partition_str(h)}" for h in hooks]
    print(f"Hooks: {', '.join(hook_names)}")

    all_parts = list(generate_partitions(n))
    m = n*(n-1)//2 - (n-1)
    print(f"Free bits: {m}, total: {2**m}")

    seen = {}
    count_tested = 0
    count_hook_pos = 0
    min_hook_val = Fraction(1)
    failures = []

    for bits in range(2**m):
        adj = tournament_from_bits(bits, n)
        H = count_ham_paths(n, adj)
        scores = tuple(sorted([sum(adj.get((i,j),0) for j in range(n) if j!=i) for i in range(n)]))

        key = (scores, H)
        if key in seen: continue
        seen[key] = bits
        count_tested += 1

        p_coeffs = compute_p_expansion(n, adj)

        hook_coeffs = {}
        for h in hooks:
            coeff = Fraction(0)
            for mu in all_parts:
                if mu not in p_coeffs: continue
                chi = chi_hook_mn(h, mu)
                coeff += Fraction(chi) * Fraction(p_coeffs[mu]) / Fraction(z_lambda(mu))
            hook_coeffs[h] = coeff

        is_pos = all(v >= 0 for v in hook_coeffs.values())
        if is_pos:
            count_hook_pos += 1
        else:
            failures.append((bits, H, scores, hook_coeffs))

        for v in hook_coeffs.values():
            if v < min_hook_val:
                min_hook_val = v

        if count_tested % 50 == 0:
            print(f"  tested {count_tested}, hook-positive so far: {count_hook_pos}")

    print(f"\n  RESULT: Hook-positive {count_hook_pos}/{count_tested}")
    print(f"  Min hook coefficient: {float(min_hook_val):.8f}")

    if failures:
        print(f"\n  FAILURES ({len(failures)} total):")
        for bits, H, scores, hc in failures[:5]:
            print(f"    bits={bits}, H={H}, scores={scores}")
            for h in hooks:
                if hc[h] < 0:
                    print(f"      [{partition_str(h)}] = {hc[h]} = {float(hc[h]):.8f}")
    else:
        print(f"\n  ALL TOURNAMENTS HOOK-POSITIVE AT n=7!")

if __name__ == '__main__':
    main()
