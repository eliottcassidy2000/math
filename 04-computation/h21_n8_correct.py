#!/usr/bin/env python3
"""
CORRECT analysis of (alpha_1, i_2) at n=8 with CORRECT cycle counting.

Previous analysis had a buggy 5-cycle counter that missed many 5-cycles.
This uses brute-force permutation enumeration for correctness.

Instance: kind-pasteur-2026-03-07-S33
"""

from itertools import combinations, permutations
from collections import Counter, defaultdict
import random
import time

def find_3cycles(adj, n):
    """Find all 3-cycle vertex sets (each has exactly 1 directed 3-cycle)."""
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.append(frozenset([i,j,k]))
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    cycles.append(frozenset([i,j,k]))
    return cycles

def count_directed_5cycles(adj, verts):
    """Count directed 5-cycles on given 5-vertex set."""
    v = list(verts)
    count = 0
    for perm in permutations(range(5)):
        ok = True
        for idx in range(5):
            if not adj[v[perm[idx]]][v[perm[(idx+1)%5]]]:
                ok = False
                break
        if ok:
            count += 1
    return count // 5  # cyclic rotations

def find_all_odd_cycles_correct(adj, n):
    """Find ALL directed odd cycles, returning (vertex_set, count) pairs.
    Each vertex_set maps to the number of DISTINCT directed cycles on it.
    In Omega(T), each directed cycle is a separate vertex.
    Multiple on same vertex set form a clique (all share all vertices).
    """
    cycles = []

    # 3-cycles: each vertex triple has 0 or 1 directed 3-cycle
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.append(frozenset([i,j,k]))
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    cycles.append(frozenset([i,j,k]))

    # 5-cycles
    if n >= 5:
        for verts in combinations(range(n), 5):
            vs = frozenset(verts)
            d = count_directed_5cycles(adj, verts)
            for _ in range(d):
                cycles.append(vs)

    # 7-cycles
    if n >= 7:
        for verts in combinations(range(n), 7):
            v = list(verts)
            dp = {}
            dp[(1, 0)] = 1
            for S in range(1, 128):
                for i in range(7):
                    if not (S & (1 << i)):
                        continue
                    if (S, i) not in dp:
                        continue
                    c = dp[(S, i)]
                    for j in range(7):
                        if S & (1 << j):
                            continue
                        if adj[v[i]][v[j]]:
                            key = (S | (1 << j), j)
                            dp[key] = dp.get(key, 0) + c
            count = 0
            for j in range(1, 7):
                if (127, j) in dp and adj[v[j]][v[0]]:
                    count += dp[(127, j)]
            d = count // 7
            vs = frozenset(verts)
            for _ in range(d):
                cycles.append(vs)

    return cycles


def compute_alpha_i2_from_cycles(cycles):
    """Compute alpha_1 (= len(cycles)) and i_2 (= # vertex-disjoint pairs)."""
    alpha1 = len(cycles)
    i2 = 0
    for a in range(alpha1):
        for b in range(a+1, alpha1):
            if not (cycles[a] & cycles[b]):
                i2 += 1
    return alpha1, i2


def independence_poly_at_2(cycles):
    """Compute I(Omega(T), 2) exactly using the cycle list."""
    nc = len(cycles)
    total = 0
    for mask in range(1 << nc):
        bits = [i for i in range(nc) if mask & (1 << i)]
        independent = True
        for a in range(len(bits)):
            for b in range(a+1, len(bits)):
                if cycles[bits[a]] & cycles[bits[b]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            total += 2**len(bits)
    return total


def held_karp(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for S in range(1, full + 1):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full])


def main():
    # First: verify OCF at n=6 with correct counting
    print("=== VERIFICATION at n=6 (exhaustive) ===")
    n = 6
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    ocf_ok = 0
    ocf_fail = 0
    h21_count = 0

    for bits in range(2**15):
        adj = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                adj[j][i] = 1
            else:
                adj[i][j] = 1

        cycles = find_all_odd_cycles_correct(adj, n)
        H = held_karp(adj, n)

        if len(cycles) <= 20:
            I_val = independence_poly_at_2(cycles)
            if H == I_val:
                ocf_ok += 1
            else:
                ocf_fail += 1
                if ocf_fail <= 3:
                    print(f"  OCF FAIL: H={H}, I={I_val}, #cycles={len(cycles)}")

        if H == 21:
            h21_count += 1

    print(f"  OCF verified: {ocf_ok} ok, {ocf_fail} fail")
    print(f"  H=21 count: {h21_count}")

    # Now n=7 sampling with correct counting
    print("\n=== n=7 SAMPLING (5000) with correct cycle counting ===")
    n = 7
    random.seed(42)
    alpha_i2 = defaultdict(lambda: defaultdict(int))
    h_vals = []

    for trial in range(5000):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        cycles = find_all_odd_cycles_correct(adj, n)
        alpha1 = len(cycles)

        i2 = 0
        for a in range(alpha1):
            for b in range(a+1, alpha1):
                if not (cycles[a] & cycles[b]):
                    i2 += 1

        alpha_i2[alpha1][i2] += 1

        H = held_karp(adj, n)
        if H == 21:
            print(f"  H=21 FOUND at trial {trial}!")
        h_vals.append(H)

    print(f"  H=21 count: {h_vals.count(21)}")
    print(f"  alpha_1 near 10:")
    for a1 in range(6, 15):
        if a1 in alpha_i2:
            i2_dict = alpha_i2[a1]
            i2_vals = sorted(i2_dict.keys())
            total = sum(i2_dict.values())
            needed = (10 - a1) // 2 if (10 - a1) >= 0 and (10 - a1) % 2 == 0 else None
            tag = ""
            if needed is not None:
                if needed in i2_vals:
                    tag = f" [i_2={needed} EXISTS]"
                else:
                    tag = f" [BLOCKED: need i_2={needed}]"
            print(f"    alpha_1={a1:2d}: i_2 in {{{','.join(str(v) for v in i2_vals)}}}, count={total}{tag}")

    # Now n=8 sampling
    print("\n=== n=8 SAMPLING (2000) with correct cycle counting ===")
    n = 8
    random.seed(123)
    alpha_i2 = defaultdict(lambda: defaultdict(int))
    h21_found = 0

    t0 = time.time()
    for trial in range(2000):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        cycles = find_all_odd_cycles_correct(adj, n)
        alpha1 = len(cycles)

        if alpha1 > 25:  # too many to compute i2 efficiently
            continue

        i2 = 0
        for a in range(alpha1):
            for b in range(a+1, alpha1):
                if not (cycles[a] & cycles[b]):
                    i2 += 1

        alpha_i2[alpha1][i2] += 1

        H = held_karp(adj, n)
        if H == 21:
            h21_found += 1
            print(f"  H=21 at trial {trial}! alpha_1={alpha1}, i_2={i2}")

    elapsed = time.time() - t0
    print(f"Done in {elapsed:.1f}s")
    print(f"H=21 found: {h21_found}")

    for a1 in range(6, 15):
        if a1 in alpha_i2:
            i2_dict = alpha_i2[a1]
            i2_vals = sorted(i2_dict.keys())
            total = sum(i2_dict.values())
            needed = (10 - a1) // 2 if (10 - a1) >= 0 and (10 - a1) % 2 == 0 else None
            tag = ""
            if needed is not None:
                if needed in i2_vals:
                    tag = f" [i_2={needed} EXISTS]"
                else:
                    tag = f" [BLOCKED: need i_2={needed}]"
            print(f"    alpha_1={a1:2d}: i_2 in {{{','.join(str(v) for v in i2_vals)}}}, count={total}{tag}")


if __name__ == "__main__":
    main()
