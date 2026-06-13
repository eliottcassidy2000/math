#!/usr/bin/env python3
"""alpha_constraint_analysis.py -- Investigate the alpha_1 + 2*alpha_2 constraint.

Session: kind-pasteur-2026-03-20-S1 (continued)

At n=6 regular, the SC Maximizer achieves alpha_1 + 2*alpha_2 = 22.
This generalizes: for max H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...,
the constraint is on the weighted sum.

Questions:
1. What is the maximum of alpha_1 + 2*alpha_2 at n=6 for each score class?
2. Does SC always achieve this maximum?
3. At n=7, what is the constraint? (alpha_3=0, so same form)
4. Is there a formula for the max of alpha_1 + 2*alpha_2 in terms of score?
"""

from itertools import combinations, permutations
from collections import defaultdict, Counter
from math import comb

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

def is_score_sc(score):
    n = len(score)
    s = sorted(score)
    return all(s[i] + s[n-1-i] == n-1 for i in range(n))

def count_directed_odd_cycles(adj, n):
    """Count ALL directed odd cycles. Return list of (frozenset, tuple)."""
    cycles = []
    for size in range(3, n+1, 2):
        for verts in combinations(range(n), size):
            sub_adj = [[adj[verts[i]][verts[j]] for j in range(size)] for i in range(size)]
            for perm in permutations(range(1, size)):
                path = [0] + list(perm)
                ok = True
                for idx in range(size):
                    if not sub_adj[path[idx]][path[(idx+1) % size]]:
                        ok = False
                        break
                if ok:
                    cycles.append((frozenset(verts), tuple(verts[path[i]] for i in range(size))))
    return cycles

def alpha_vector(adj, n):
    """Compute [alpha_0, alpha_1, alpha_2, ...]."""
    all_cycles = count_directed_odd_cycles(adj, n)
    num = len(all_cycles)
    if num == 0:
        return [1]

    vsets = [c[0] for c in all_cycles]
    conflict = [[False]*num for _ in range(num)]
    for i in range(num):
        for j in range(i+1, num):
            if vsets[i] & vsets[j]:
                conflict[i][j] = True
                conflict[j][i] = True

    alpha = [1]
    for k in range(1, num + 1):
        count = 0
        for subset in combinations(range(num), k):
            ok = True
            for a in range(len(subset)):
                for b in range(a+1, len(subset)):
                    if conflict[subset[a]][subset[b]]:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                count += 1
        if count == 0:
            break
        alpha.append(count)
    return alpha


def main():
    print("=" * 70)
    print("ALPHA_1 + 2*ALPHA_2 CONSTRAINT ANALYSIS")
    print("=" * 70)

    # Analysis at n=5
    print("\n--- n=5 ---")
    n = 5
    m = n * (n - 1) // 2
    by_score = defaultdict(list)

    for bits in range(1 << m):
        adj = adj_from_bits(bits, n)
        sc = score_seq(adj, n)
        H = held_karp_H(adj, n)
        by_score[sc].append({'bits': bits, 'H': H})

    for sc in sorted(by_score.keys()):
        entries = by_score[sc]
        max_H = max(e['H'] for e in entries)
        is_sc_score = is_score_sc(sc)
        # For max_H entry, compute alpha
        for e in entries:
            if e['H'] == max_H:
                adj = adj_from_bits(e['bits'], n)
                av = alpha_vector(adj, n)
                weighted = sum(av[k] * 2**k for k in range(len(av))) - 1  # H-1
                a1 = av[1] if len(av) > 1 else 0
                a2 = av[2] if len(av) > 2 else 0
                print(f"  score {sc} (SC={is_sc_score}): max H={max_H}, "
                      f"IP={av}, a1+2*a2={a1+2*a2}, weighted_sum={weighted}")
                break

    # Analysis at n=6 — ALL score classes
    print("\n--- n=6 ALL SCORE CLASSES ---")
    n = 6
    m = n * (n - 1) // 2
    by_score = defaultdict(list)

    for bits in range(1 << m):
        adj = adj_from_bits(bits, n)
        sc = score_seq(adj, n)
        H = held_karp_H(adj, n)
        by_score[sc].append({'bits': bits, 'H': H})

    print(f"  {len(by_score)} score classes")

    for sc in sorted(by_score.keys()):
        entries = by_score[sc]
        max_H = max(e['H'] for e in entries)
        is_sc_score = is_score_sc(sc)

        # Compute alpha for max_H entry
        for e in entries:
            if e['H'] == max_H:
                adj = adj_from_bits(e['bits'], n)
                av = alpha_vector(adj, n)
                a1 = av[1] if len(av) > 1 else 0
                a2 = av[2] if len(av) > 2 else 0
                print(f"  score {sc} (SC={is_sc_score}): max H={max_H}, "
                      f"IP={av}, a1+2*a2={a1+2*a2}")
                break

    # Special analysis: for SC score classes, check if max a1+2*a2 comes from SC tournament
    print("\n--- SC SCORE CLASSES: Is max alpha_1+2*alpha_2 always from SC? ---")
    for sc in sorted(by_score.keys()):
        if not is_score_sc(sc):
            continue
        entries = by_score[sc]

        # Compute alpha for ALL entries in this score class (expensive!)
        # Just do a sample for non-regular scores
        if sc == (2, 2, 2, 3, 3, 3):
            # Full analysis already done
            print(f"  score {sc}: see THM-255 for complete analysis")
            continue

        # For other SC scores, check max_H entries
        max_H = max(e['H'] for e in entries)
        max_entries = [e for e in entries if e['H'] == max_H]

        sample = max_entries[:5]
        for e in sample:
            adj = adj_from_bits(e['bits'], n)
            av = alpha_vector(adj, n)
            a1 = av[1] if len(av) > 1 else 0
            a2 = av[2] if len(av) > 2 else 0
            print(f"  score {sc}: max H={max_H}, IP={av}, a1+2*a2={a1+2*a2}")
            break


if __name__ == "__main__":
    main()
