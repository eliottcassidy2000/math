#!/usr/bin/env python3
"""
cycle_cover_deep_investigation.py — Unifying per(A), ham(A), I(Ω,2), and Q(T)

THREE PARTITION FUNCTIONS ON THE SAME TOURNAMENT:

1. per(A) = Σ over vertex-covering cycle partitions (any length ≥ 3) of 1
   Counts: cycle covers. Standard permanent of adjacency matrix.

2. I(Ω(T), 2) = H(T) = Σ_k α_k · 2^k
   Counts: vertex-disjoint ODD-cycle collections (partial or total), weighted 2^k.
   This is the Hamiltonian path count (by OCF = THM-002).

3. Q(T) = Σ over vertex-COVERING odd-cycle partitions of 2^k
   Counts: cycle covers using ONLY odd-length cycles, weighted by 2^(#cycles).
   NEW: the "odd-cycle permanent" — the covering part of I(Ω(T),2).

DECOMPOSITIONS:
   H(T) = Q(T) + R(T)    where R(T) = partial packing correction (non-covering terms)
   per(A) = per_odd(A) + per_even(A)   where per_even counts covers with ≥1 even cycle

KEY QUESTION: Does Q(T) or per(A) explain H(T) better?
              What does the ratio Q/H look like across classes?
              Is R(T) = H(T) - Q(T) structured?
"""

import sys, time
from itertools import permutations, combinations
from collections import defaultdict, Counter
from math import factorial, gcd
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

# ================================================================
# CORE: PERMUTATION / CYCLE UTILITIES
# ================================================================

def cycle_decomp(sigma):
    """Return list of cycles as lists of vertices."""
    n = len(sigma)
    visited = [False] * n
    cycles = []
    for i in range(n):
        if visited[i]: continue
        cyc = []
        j = i
        while not visited[j]:
            visited[j] = True
            cyc.append(j)
            j = sigma[j]
        cycles.append(cyc)
    return cycles

def is_derangement(sigma):
    return all(sigma[i] != i for i in range(len(sigma)))


# ================================================================
# TOURNAMENT ENUMERATION
# ================================================================

def all_tournaments_n(n):
    arcs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(arcs)
    perms_list = list(permutations(range(n)))

    def adj_from_bits(bits):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(arcs):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        return A

    def canon(A):
        best = None
        for p in perms_list:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best

    classes = {}
    for bits in range(1 << m):
        A = adj_from_bits(bits)
        c = canon(A)
        if c not in classes:
            classes[c] = A

    return sorted(classes.items()), perms_list


# ================================================================
# CYCLE COVER DECOMPOSITION
# ================================================================

def cycle_cover_analysis(adj, n):
    """Full cycle cover analysis of a tournament.

    Returns:
        per_A: permanent (total cycle covers)
        ham_A: Hamiltonian circuits (single-cycle covers)
        C_k: dict k → count of k-cycle covers
        per_odd: covers using ONLY odd-length cycles
        per_even: covers using ≥1 even-length cycle
        Q_T: Σ (2^k) over odd-only vertex-covering cycle partitions
        cover_types: detailed breakdown by cycle-length partition
    """
    C_k = defaultdict(int)
    cover_types = defaultdict(int)  # partition of n → count
    per_odd = 0
    per_even = 0
    Q_T = 0

    for sigma in permutations(range(n)):
        # Check if contributes to permanent
        prod = 1
        for i in range(n):
            prod *= adj[i][sigma[i]]
            if prod == 0: break
        if prod == 0: continue

        # Decompose into cycles
        cycles = cycle_decomp(list(sigma))
        k = len(cycles)
        lengths = tuple(sorted([len(c) for c in cycles], reverse=True))

        C_k[k] += 1
        cover_types[lengths] += 1

        all_odd = all(len(c) % 2 == 1 for c in cycles)
        if all_odd:
            per_odd += 1
            Q_T += (1 << k)  # 2^k
        else:
            per_even += 1

    per_A = sum(C_k.values())
    ham_A = C_k.get(1, 0)

    return {
        'per': per_A, 'ham': ham_A, 'C_k': dict(C_k),
        'per_odd': per_odd, 'per_even': per_even,
        'Q': Q_T, 'cover_types': dict(cover_types)
    }


# ================================================================
# ODD-CYCLE ENUMERATION & INDEPENDENCE POLYNOMIAL
# ================================================================

def find_all_odd_cycles(adj, n):
    """Find ALL directed odd cycles in the tournament."""
    cycles = []
    for length in range(3, n+1, 2):  # odd lengths only
        for vertices in combinations(range(n), length):
            # Try all cyclic orderings
            vlist = list(vertices)
            for perm in permutations(vlist[1:]):
                order = [vlist[0]] + list(perm)
                # Check if this is a directed cycle
                is_cycle = True
                for i in range(length):
                    if not adj[order[i]][order[(i+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    # Canonical form: start from minimum vertex, direction matters
                    min_pos = order.index(min(order))
                    canonical = tuple(order[min_pos:] + order[:min_pos])
                    cycles.append(canonical)

    # Deduplicate
    return list(set(cycles))


def independence_polynomial(cycles, n):
    """Compute independence polynomial I(Ω,x) of the odd-cycle conflict graph.

    Ω has vertices = odd directed cycles, edges = shared vertex.
    An independent set = vertex-disjoint collection of odd cycles.
    I(Ω, x) = Σ_k α_k x^k where α_k = # independent sets of size k.

    Returns: list alpha where alpha[k] = number of independent sets of size k.
    Also returns: Q_from_indep = Σ over COVERING independent sets of 2^k.
    """
    nc = len(cycles)

    # Build conflict graph
    # Two cycles conflict if they share a vertex
    cycle_vertices = [set(c) for c in cycles]
    conflicts = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycle_vertices[i] & cycle_vertices[j]:
                conflicts[i][j] = conflicts[j][i] = True

    # Enumerate all independent sets by brute force (feasible for small n)
    alpha = defaultdict(int)
    Q_from_indep = 0
    covering_sets = []

    all_vertices = set(range(n))

    for mask in range(1 << nc):
        # Check independence
        selected = [i for i in range(nc) if (mask >> i) & 1]
        k = len(selected)

        is_independent = True
        for i in range(len(selected)):
            for j in range(i+1, len(selected)):
                if conflicts[selected[i]][selected[j]]:
                    is_independent = False
                    break
            if not is_independent: break

        if not is_independent: continue

        alpha[k] += 1

        # Check if covering (vertices of selected cycles = all vertices)
        covered = set()
        for idx in selected:
            covered |= cycle_vertices[idx]

        if covered == all_vertices and k > 0:
            Q_from_indep += (1 << k)
            covering_sets.append((k, [cycles[i] for i in selected]))

    return dict(alpha), Q_from_indep, covering_sets


def count_H_dp(adj, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)) or dp[S][v] == 0: continue
            for w in range(n):
                if not (S & (1 << w)) and adj[v][w]:
                    dp[S | (1 << w)][w] += dp[S][v]
    return sum(dp[full])


def score_seq(adj, n):
    return tuple(sorted(sum(adj[i]) for i in range(n)))


def c3_count(adj, n):
    c = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if i != j and j != k and k != i:
                    if adj[i][j] and adj[j][k] and adj[k][i]:
                        c += 1
    return c // 3


# ================================================================
# MAIN INVESTIGATION
# ================================================================

print("DEEP CYCLE COVER INVESTIGATION")
print("Connecting per(A), ham(A), I(Ω,2)=H(T), Q(T), and cycle packing")
print("=" * 85)

for n in range(3, 7):
    t0 = time.time()
    print(f"\n{'='*85}")
    print(f"n = {n}")
    print(f"{'='*85}")

    class_list, perms_list = all_tournaments_n(n)

    results = []
    for ci, (canon, adj) in enumerate(class_list):
        H = count_H_dp(adj, n)
        scores = score_seq(adj, n)
        t3 = c3_count(adj, n)

        # Cycle cover analysis
        cc = cycle_cover_analysis(adj, n)

        # Independence polynomial of Ω(T)
        odd_cycles = find_all_odd_cycles(adj, n)
        alpha, Q_indep, covering = independence_polynomial(odd_cycles, n)

        # Verify H = I(Ω, 2)
        H_from_indep = sum(alpha.get(k, 0) * (1 << k) for k in range(max(alpha.keys())+1 if alpha else 1))

        # R(T) = H(T) - Q(T) = partial packing correction
        R = H - Q_indep

        results.append({
            'ci': ci, 'H': H, 'scores': scores, 't3': t3,
            'per': cc['per'], 'ham': cc['ham'], 'C_k': cc['C_k'],
            'per_odd': cc['per_odd'], 'per_even': cc['per_even'],
            'Q_cc': cc['Q'], 'Q_indep': Q_indep,
            'alpha': alpha, 'R': R,
            'n_odd_cycles': len(odd_cycles),
            'cover_types': cc['cover_types'],
            'H_check': H_from_indep,
            'covering': covering,
        })

    elapsed = time.time() - t0

    # Print detailed table
    print(f"\n{len(results)} classes, computed in {elapsed:.1f}s")
    print(f"\n{'#':>2} {'H':>4} {'t3':>3} {'per':>5} {'ham':>4} {'Q':>5} {'R':>5} "
          f"{'H=I?':>4} {'Q=Qc':>4} {'#odd':>4} {'alpha_k':>20} {'cover types'}")
    print("-" * 110)

    for d in sorted(results, key=lambda x: x['H']):
        alpha_str = ','.join(f'{k}:{v}' for k,v in sorted(d['alpha'].items()) if k > 0)
        ct_str = ' '.join(f'{k}:{v}' for k,v in sorted(d['cover_types'].items()))
        h_ok = '✓' if d['H'] == d['H_check'] else '✗'
        q_ok = '✓' if d['Q_cc'] == d['Q_indep'] else '✗'

        print(f"{d['ci']:2d} {d['H']:4d} {d['t3']:3d} {d['per']:5d} {d['ham']:4d} "
              f"{d['Q_indep']:5d} {d['R']:5d} {h_ok:>4} {q_ok:>4} {d['n_odd_cycles']:4d} "
              f"{alpha_str:>20} {ct_str}")

    # ================================================================
    # AGGREGATE ANALYSIS
    # ================================================================

    print(f"\n--- Aggregate Analysis for n={n} ---")

    # 1. Q/H ratio
    qh_ratios = [(d['Q_indep']/d['H'], d['H']) for d in results if d['H'] > 0]
    if qh_ratios:
        print(f"\nQ(T)/H(T) ratio (covering fraction of independence polynomial):")
        for ratio, H in sorted(qh_ratios):
            print(f"  H={H:4d}: Q/H = {ratio:.4f}  (R = {round((1-ratio)*H)})")

    # 2. per/H ratio
    ph = [(d['per']/d['H'] if d['H']>0 else 0, d['H'], d['per']) for d in results]
    print(f"\nper(A)/H(T) ratio (permanent vs Hamiltonian paths):")
    for ratio, H, per in sorted(ph):
        if H > 0:
            print(f"  H={H:4d}: per/H = {ratio:.4f}  (per={per})")

    # 3. Correlation between t3 and C_2
    t3s = [d['t3'] for d in results]
    c2s = [d['C_k'].get(2, 0) for d in results]
    if max(c2s) > 0:
        corr = np.corrcoef(t3s, c2s)[0, 1] if len(set(t3s)) > 1 else 0
        print(f"\nCorrelation(t3, C_2) = {corr:.4f}")
        print(f"  t3 values: {sorted(set(t3s))}")
        print(f"  C_2 values: {sorted(set(c2s))}")

        # Detailed t3 vs C_2
        for d in sorted(results, key=lambda x: x['t3']):
            if d['C_k'].get(2, 0) > 0 or d['t3'] > 0:
                print(f"    t3={d['t3']:2d}, C_2={d['C_k'].get(2,0):2d}, H={d['H']:4d}, per={d['per']:4d}")

    # 4. Even-cycle contributions
    even_count = sum(1 for d in results if d['per_even'] > 0)
    print(f"\nClasses with even-length cycle covers: {even_count}/{len(results)}")
    if even_count > 0:
        for d in results:
            if d['per_even'] > 0:
                print(f"  #{d['ci']}: per_even={d['per_even']}, per_odd={d['per_odd']}, "
                      f"cover_types={d['cover_types']}")

    # 5. The R(T) = H(T) - Q(T) decomposition
    print(f"\nH(T) = Q(T) + R(T) decomposition:")
    print(f"  Q(T) = contribution from vertex-covering disjoint odd-cycle collections")
    print(f"  R(T) = contribution from PARTIAL odd-cycle packings (+ empty set)")
    for d in sorted(results, key=lambda x: x['H']):
        if d['H'] > 0:
            print(f"  H={d['H']:4d} = Q={d['Q_indep']:4d} + R={d['R']:4d}  "
                  f"(Q/H={d['Q_indep']/d['H']:.3f})")

    # 6. Alpha distribution shape
    all_alpha = defaultdict(list)
    for d in results:
        for k, v in d['alpha'].items():
            all_alpha[k].append(v)
    print(f"\nAlpha_k (independent set counts) statistics:")
    for k in sorted(all_alpha.keys()):
        vals = all_alpha[k]
        print(f"  α_{k}: min={min(vals)}, max={max(vals)}, mean={sum(vals)/len(vals):.1f}")

print(f"""

{'='*85}
THEORETICAL SUMMARY
{'='*85}

Three partition functions on tournament T:

  per(A) = Σ (cycle covers of ALL vertices using cycles of length ≥ 3)
  H(T)   = I(Ω(T), 2) = Σ_k α_k · 2^k  (disjoint ODD-cycle collections, partial or total)
  Q(T)   = Σ (2^k over vertex-COVERING disjoint odd-cycle partitions)

Decompositions:
  H(T) = Q(T) + R(T)       where R counts partial odd-cycle packings
  per(A) = per_odd + per_even   where per_even uses ≥1 even-length cycle

Key observations:
  1. At n≤5: per(A) = ham(A) (no multi-cycle covers possible)
  2. At n=6: C_2 counts disjoint 3-cycle PAIRS (the first multi-cycle covers)
  3. Q(T) from cycle covers EQUALS Q(T) from independence polynomial (verified)
  4. R(T) = H(T) - Q(T) measures the "partial packing surplus"
  5. The possible cycle-count k in covers = partitions of n into parts ≥ 3
""")
