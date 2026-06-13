#!/usr/bin/env python3
"""
burnside_gf2_bridge_s263.py — Bridge between GF(2) decomposition and Burnside
opus-2026-03-23-S263

Key question: How does the GF(2) decomposition Edge = Cut ⊕ Cycle (odd n)
interact with the Burnside formula for tournaments and even graphs?

Specifically:
1. For each permutation type σ, how do arc orbits decompose into cut and cycle?
2. The fiber ratio V_tourn/V_even = average number of score assignments per even graph
3. Does the Burnside sum factor as Burnside(cut) × Burnside(cycle)?
4. Connection between edge-centric Burnside and GF(2) projection
5. Why does the double Burnside fail at n=5? Is it related to multi-wire?
"""

from math import comb, factorial, gcd, lcm
from collections import Counter, defaultdict
import numpy as np
from itertools import permutations

def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield []; return
    for f in range(min(n, mx), 0, -1):
        for r in partitions(n - f, f): yield [f] + r

def ccs(n, ct):
    c = Counter(ct); r = factorial(n)
    for l, k in c.items(): r //= (l ** k) * factorial(k)
    return r

def fix_tournament(ct):
    """Number of tournaments fixed by a permutation with cycle type ct."""
    for c in ct:
        if c % 2 == 0: return 0
    exp = sum((c - 1) // 2 for c in ct)
    for i in range(len(ct)):
        for j in range(i + 1, len(ct)):
            exp += gcd(ct[i], ct[j])
    return 2 ** exp

def edge_orbits(n, ct):
    """Number of edge orbits of K_n under permutation with cycle type ct."""
    total = 0
    for i in range(len(ct)):
        total += ct[i] // 2
        for j in range(i + 1, len(ct)):
            total += gcd(ct[i], ct[j])
    return total

def fix_even_graph(n, ct):
    """Number of even graphs (all even degrees) fixed by permutation with cycle type ct.
    = 2^(edge_orbits - rank(degree_constraint_matrix))"""
    # Build representative permutation
    perm = [0] * n; pos = 0
    for c in ct:
        for i in range(c):
            perm[pos + i] = pos + (i + 1) % c
        pos += c

    P = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(P)
    eidx = {e: i for i, e in enumerate(P)}

    # Find edge orbits
    visited = [False] * m
    orbits = []
    for idx in range(m):
        if visited[idx]: continue
        orbit = []
        i, j = P[idx]
        while True:
            e = (min(i,j), max(i,j))
            ei = eidx[e]
            if visited[ei]: break
            visited[ei] = True
            orbit.append(ei)
            i, j = perm[i], perm[j]
        orbits.append(orbit)

    num_orbits = len(orbits)

    # Degree constraint matrix A (n × num_orbits) over GF(2)
    A = np.zeros((n, num_orbits), dtype=np.int8)
    for oidx, orbit in enumerate(orbits):
        for ei in orbit:
            i, j = P[ei]
            A[i, oidx] = (A[i, oidx] + 1) % 2
            A[j, oidx] = (A[j, oidx] + 1) % 2

    # Rank of A over GF(2)
    M = A.copy()
    rows, cols = M.shape
    rank = 0
    for col in range(cols):
        found = -1
        for r in range(rank, rows):
            if M[r, col]: found = r; break
        if found == -1: continue
        M[[rank, found]] = M[[found, rank]]
        for r in range(rows):
            if r != rank and M[r, col]:
                M[r] = (M[r] + M[rank]) % 2
        rank += 1

    return 2 ** (num_orbits - rank), num_orbits, rank

def analyze_burnside_decomposition(n):
    """Analyze how Burnside terms decompose by cut/cycle space."""
    print(f"\n{'#'*70}")
    print(f"  n = {n}")
    print(f"{'#'*70}")

    m = comb(n, 2)
    nf = factorial(n)
    cut_dim = n - 1
    cycle_dim = m - cut_dim

    print(f"  m = {m}, cut_dim = {cut_dim}, cycle_dim = {cycle_dim}")

    # For each cycle type, compute:
    # - fix_tourn (= 2^{arc_orbits} for odd cycle types, 0 for even)
    # - fix_even (= 2^{edge_orbits - rank_A} for the even graph constraint)
    # - fix_graph (= 2^{edge_orbits})
    # - The "cut contribution" = fix_tourn / fix_even

    V_tourn = 0
    V_even = 0
    V_graph = 0

    print(f"\n  {'cycle_type':>15} {'|class|':>8} {'fix_T':>10} {'fix_E':>10} {'fix_G':>10} "
          f"{'e_orb':>6} {'rank':>5} {'null':>5} {'ratio_T/E':>10}")

    detail_rows = []

    for ct in partitions(n):
        ct_list = list(ct)
        cc = ccs(n, ct_list)
        ft = fix_tournament(ct_list)
        fe, eo, rk = fix_even_graph(n, ct_list)
        fg = 2 ** eo

        ratio = ft / fe if fe > 0 else float('inf')

        V_tourn += cc * ft
        V_even += cc * fe
        V_graph += cc * fg

        detail_rows.append({
            'ct': ct_list, 'cc': cc, 'ft': ft, 'fe': fe, 'fg': fg,
            'eo': eo, 'rk': rk, 'null': eo - rk, 'ratio': ratio
        })

        if n <= 8:
            print(f"  {str(ct_list):>15} {cc:>8} {ft:>10} {fe:>10} {fg:>10} "
                  f"{eo:>6} {rk:>5} {eo-rk:>5} {ratio:>10.3f}")

    V_tourn //= nf
    V_even //= nf
    V_graph //= nf

    print(f"\n  V_tourn = {V_tourn}, V_even = {V_even}, V_graph = {V_graph}")
    print(f"  V_tourn / V_even = {V_tourn / V_even:.4f}")
    print(f"  V_graph / V_tourn = {V_graph / V_tourn:.4f}")

    # Key question: does fix_tourn = fix_even × 2^{cut_free_bits}?
    # If so, the Burnside sum factors.
    print(f"\n  FACTORIZATION CHECK:")
    print(f"  Is fix_T = fix_E × 2^k for some k at each cycle type?")

    factors_consistently = True
    for d in detail_rows:
        if d['ft'] == 0 and d['fe'] > 0:
            # Even cycle type: tournament fixes nothing, even graph fixes something
            pass
        elif d['ft'] > 0 and d['fe'] > 0:
            ratio = d['ft'] / d['fe']
            import math
            if ratio > 0 and ratio == int(ratio) and (int(ratio) & (int(ratio)-1) == 0):
                k = int(math.log2(ratio))
                if n <= 8:
                    print(f"    {str(d['ct']):>15}: fix_T/fix_E = 2^{k} = {int(ratio)}")
            else:
                factors_consistently = False
                if n <= 8:
                    print(f"    {str(d['ct']):>15}: fix_T/fix_E = {ratio:.4f} (NOT a power of 2!)")

    if factors_consistently:
        print(f"  => YES, fix_T = fix_E × 2^k at every odd cycle type")
    else:
        print(f"  => NO, the ratio is not always a power of 2")

    # For odd cycle types only: analyze the exponent k
    print(f"\n  EXPONENT ANALYSIS (odd cycle types only):")
    for d in detail_rows:
        if d['ft'] == 0: continue
        import math
        k = int(round(math.log2(d['ratio']))) if d['ratio'] > 0 else 0
        arc_orbits_tourn = int(round(math.log2(d['ft'])))
        # arc_orbits = sum (c_i-1)/2 + sum gcd(c_i,c_j)
        intra = sum((c-1)//2 for c in d['ct'])
        inter = sum(gcd(d['ct'][i], d['ct'][j]) for i in range(len(d['ct'])) for j in range(i+1, len(d['ct'])))
        if n <= 8:
            print(f"    {str(d['ct']):>15}: arc_orbits={arc_orbits_tourn}, "
                  f"cycle_nullity={d['null']}, cut_free={k}, "
                  f"intra={intra}, inter={inter}")

    # Check: arc_orbits = cycle_nullity + cut_free_bits?
    print(f"\n  DECOMPOSITION CHECK: arc_orbits = cycle_null + cut_free?")
    all_match = True
    for d in detail_rows:
        if d['ft'] == 0: continue
        import math
        arc_orb = int(round(math.log2(d['ft'])))
        k = int(round(math.log2(d['ratio'])))
        check = (arc_orb == d['null'] + k)
        if n <= 8:
            print(f"    {str(d['ct']):>15}: {arc_orb} = {d['null']} + {k} → {'✓' if check else '✗'}")
        if not check: all_match = False

    if all_match:
        print(f"  => YES! Arc orbits = cycle nullity + cut free bits")
        print(f"  => The Burnside exponent FACTORS as cycle + cut contributions")
    else:
        print(f"  => NO, decomposition doesn't hold for all types")

    return V_tourn, V_even, V_graph

if __name__ == '__main__':
    print("=" * 70)
    print("  BURNSIDE ↔ GF(2) BRIDGE")
    print("  opus-2026-03-23-S263")
    print("=" * 70)

    results = []
    for n in range(3, 10):
        Vt, Ve, Vg = analyze_burnside_decomposition(n)
        results.append((n, Vt, Ve, Vg))

    print("\n" + "=" * 70)
    print("  SUMMARY TABLE")
    print("=" * 70)
    print(f"{'n':>3} {'V_tourn':>12} {'V_even':>10} {'V_graph':>12} {'T/E ratio':>10} {'G/T ratio':>10}")
    for n, Vt, Ve, Vg in results:
        print(f"{n:>3} {Vt:>12} {Ve:>10} {Vg:>12} {Vt/Ve:>10.3f} {Vg/Vt:>10.3f}")

    print("\nDONE.")
