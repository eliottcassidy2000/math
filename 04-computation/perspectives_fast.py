#!/usr/bin/env python3
"""perspectives_fast.py -- Fast perspectives computation.

Session: kind-pasteur-2026-03-20-S4

Optimized: use score + H + c3 as invariant hash for grouping,
then check isomorphism only within groups. Skip full n! canon_form.
"""

from itertools import permutations
from collections import defaultdict, Counter
from math import factorial

def adj_from_bits(bits, n):
    adj = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    return adj

def held_karp_H(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))

def score_vector(adj, n):
    """Unsorted score vector."""
    return tuple(sum(adj[i][j] for j in range(n) if j != i) for i in range(n))

def count_3cycles(adj, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    c3 += 1
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    c3 += 1
    return c3

def are_isomorphic(adj1, adj2, n):
    """Check if two tournaments are isomorphic."""
    s1 = score_seq(adj1, n)
    s2 = score_seq(adj2, n)
    if s1 != s2:
        return False
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(i+1, n):
                if adj1[perm[i]][perm[j]] != adj2[i][j]:
                    match = False
                    break
            if not match:
                break
        if match:
            return True
    return False

def find_aut_group(adj, n):
    """Find automorphism group. Returns list of permutations."""
    auts = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(i+1, n):
                if adj[perm[i]][perm[j]] != adj[i][j]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            auts.append(perm)
    return auts

def vertex_orbit_count(adj, n):
    """Count vertex orbits under Aut(T)."""
    auts = find_aut_group(adj, n)
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry
    for aut in auts:
        for v in range(n):
            union(v, aut[v])
    return len(set(find(v) for v in range(n)))

def is_self_converse_fast(adj, n):
    """Check self-converse by looking for anti-automorphism."""
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(i+1, n):
                if adj[perm[i]][perm[j]] != (1 - adj[i][j]):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return True
    return False


def main():
    print("=" * 70)
    print("PERSPECTIVES AND CLASSES — FAST COMPUTATION")
    print("=" * 70)

    known_T = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}

    all_perspectives = {}
    all_class_counts = {}

    for n in range(1, 8):
        m = n * (n - 1) // 2
        N = 1 << m

        if n >= 7:
            # n=7 has 2M tournaments, too slow for full iso check
            all_class_counts[n] = known_T[n]
            print(f"\n  n={n}: using known T({n}) = {known_T[n]}")
            continue

        print(f"\n  n={n} ({N} labeled tournaments)")

        # Group by invariants for fast iso testing
        groups = defaultdict(list)
        for bits in range(N):
            adj = adj_from_bits(bits, n)
            sc = score_seq(adj, n)
            H = held_karp_H(adj, n)
            c3 = count_3cycles(adj, n)
            key = (sc, H, c3)
            groups[key].append(bits)

        # Within each group, find isomorphism classes
        classes = []
        for key, members in groups.items():
            reps = []  # representatives of iso classes in this group
            for bits in members:
                adj = adj_from_bits(bits, n)
                is_new = True
                for rep_bits in reps:
                    rep_adj = adj_from_bits(rep_bits, n)
                    if are_isomorphic(adj, rep_adj, n):
                        is_new = False
                        break
                if is_new:
                    reps.append(bits)

            for rep_bits in reps:
                adj = adj_from_bits(rep_bits, n)
                num_orbits = vertex_orbit_count(adj, n)
                aut_size = len(find_aut_group(adj, n))
                is_sc = is_self_converse_fast(adj, n)
                sc = score_seq(adj, n)
                H = held_karp_H(adj, n)

                classes.append({
                    'bits': rep_bits,
                    'score': sc,
                    'H': H,
                    'aut_size': aut_size,
                    'num_orbits': num_orbits,
                    'labeled': factorial(n) // aut_size,
                    'is_sc': is_sc,
                })

        num_classes = len(classes)
        total_perspectives = sum(c['num_orbits'] for c in classes)

        all_class_counts[n] = num_classes
        all_perspectives[n] = total_perspectives

        print(f"    Classes: {num_classes} (expected: {known_T.get(n, '?')})")
        print(f"    Total perspectives: {total_perspectives}")

        # Print table
        if num_classes <= 60:
            print(f"\n    {'Score':>24} {'|Aut|':>6} {'Orb':>4} {'SC':>3} {'H':>6} {'#Lab':>6}")
            for c in sorted(classes, key=lambda x: (x['score'], x['H'])):
                print(f"    {str(c['score']):>24} {c['aut_size']:>6} {c['num_orbits']:>4} "
                      f"{'Y' if c['is_sc'] else 'N':>3} {c['H']:>6} {c['labeled']:>6}")

        # SC breakdown
        sc_c = [c for c in classes if c['is_sc']]
        nsc_c = [c for c in classes if not c['is_sc']]
        print(f"\n    SC: {len(sc_c)} classes, {sum(c['num_orbits'] for c in sc_c)} perspectives")
        print(f"    Non-SC: {len(nsc_c)} classes ({len(nsc_c)//2} pairs), "
              f"{sum(c['num_orbits'] for c in nsc_c)} perspectives")

    # Summary comparison
    print("\n" + "=" * 70)
    print("SUMMARY: P(n) vs T(n+1)")
    print("=" * 70)

    print(f"\n  {'n':>4} {'T(n)':>8} {'P(n)':>8} {'T(n+1)':>8} {'P=T?':>6}")
    for n in sorted(all_perspectives.keys()):
        T_n = all_class_counts.get(n, '?')
        P_n = all_perspectives[n]
        T_next = known_T.get(n + 1, '?')
        match = 'YES' if isinstance(T_next, int) and P_n == T_next else ('NO' if isinstance(T_next, int) else '?')
        print(f"  {n:4d} {str(T_n):>8} {P_n:>8} {str(T_next):>8} {match:>6}")

    # The perspective sequence
    print(f"\n  Perspective sequence P(n):")
    for n in sorted(all_perspectives.keys()):
        print(f"    P({n}) = {all_perspectives[n]}")

    print(f"\n  Class sequence T(n):")
    for n in sorted(all_class_counts.keys()):
        print(f"    T({n}) = {all_class_counts[n]}")


if __name__ == "__main__":
    main()
