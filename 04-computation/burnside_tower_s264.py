#!/usr/bin/env python3
"""
burnside_tower_s264.py — The Tower of Hierarchy with Burnside
opus-2026-03-23-S264

The Burnside factorization reveals a TOWER of mathematical objects,
each one a quotient of the previous, with explicit Burnside formulas
connecting them.

THE TOWER (from top to bottom):

Level 5: Labeled tournaments          2^{C(n,2)}
Level 4: Tournament iso classes       V_n = A000568
Level 3: Merged classes (÷ Z_2)       V_merged = (V_n + SC_n)/2
Level 2: Even graph iso classes       V_even = A002854
Level 1: Score sequences              #score_seqs
Level 0: A single number              1

Each transition is a quotient by a group action:
  5→4: quotient by S_n (vertex relabeling)
  4→3: quotient by Z_2 (complement)
  4→2: projection to cycle space (GF(2), odd n only)
  4→1: projection to cut space

The factorization theorem connects levels:
  Fix_tourn = Fix_even × 2^{k-1}
  |Level 4| × |Level 4→2 fiber| ≈ |Level 2|
"""

from math import comb, factorial, gcd, lcm
from collections import Counter

def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield []; return
    for f in range(min(n, mx), 0, -1):
        for r in partitions(n - f, f): yield [f] + r

def ccs_val(n, ct):
    c = Counter(ct); r = factorial(n)
    for l, k in c.items(): r //= (l ** k) * factorial(k)
    return r

def arc_orbits(ct):
    total = 0
    for i in range(len(ct)):
        total += (ct[i] - 1) // 2
        for j in range(i + 1, len(ct)):
            total += gcd(ct[i], ct[j])
    return total

def edge_orbits(ct):
    total = 0
    for i in range(len(ct)):
        total += ct[i] // 2
        for j in range(i + 1, len(ct)):
            total += gcd(ct[i], ct[j])
    return total

def fix_anti(ct, n):
    """Fixed points under complement × σ (for SC counting)."""
    # Build permutation
    perm = [0]*n; pos = 0
    for c in ct:
        for i in range(c):
            perm[pos + i] = pos + (i + 1) % c
        pos += c

    P = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(P)
    eidx = {e: i for i, e in enumerate(P)}

    visited = set(); free = 0
    for i in range(n):
        for j in range(n):
            if i == j or (i,j) in visited: continue
            orbit = []; a, b = i, j
            while (a, b) not in visited:
                visited.add((a, b)); orbit.append((a, b))
                a, b = perm[a], perm[b]
            L = len(orbit); orbit_set = set(orbit); rev = (j, i)
            if rev in orbit_set:
                k_idx = orbit.index(rev)
                if k_idx % 2 == 1: free += 1
                else: return 0
            else:
                if rev in visited: pass
                else:
                    if L % 2 == 0: free += 1
                    else: return 0
    return 2 ** free

def compute_tower(n):
    """Compute all levels of the tower at a given n."""
    nf = factorial(n)
    m = comb(n, 2)

    # ===== LEVEL 5: Labeled =====
    labeled_tournaments = 2 ** m
    labeled_graphs = 2 ** m
    labeled_even = 2 ** comb(n-1, 2)
    labeled_scores = 2 ** (n - 1)  # cut space dimension

    # ===== LEVEL 4: Iso classes (÷ S_n) =====
    V_tourn = 0
    SC_n = 0
    V_even_all = 0  # using closed form for all cycle types
    V_graph = 0

    # Burnside components
    T_transition = 0
    twin_SL_sum = 0

    # Per-level Burnside detail
    for ct in partitions(n):
        ct_list = list(ct)
        cc = ccs_val(n, ct_list)
        k = len(ct_list)
        f = ct_list.count(1)
        eo = edge_orbits(ct_list)
        order = 1
        for c in ct_list: order = lcm(order, c)

        # Graphs
        fg = 2 ** eo
        V_graph += cc * fg

        # Even graphs (closed form at odd n, general at all n)
        fe_closed = 2 ** (eo - k + 1)
        V_even_all += cc * fe_closed

        # Tournaments (odd cycle types only)
        if all(c % 2 == 1 for c in ct_list):
            ao = arc_orbits(ct_list)
            ft = 2 ** ao
            V_tourn += cc * ft

            # Transition orbits
            T_transition += cc * ft * comb(f, 2)

            # Twin SL
            fix_even = 2 ** (ao - k + 1)
            twin_per = 2 * comb(f, 2) * fix_even if f >= 2 else 0
            twin_SL_sum += cc * twin_per

        # SC (anti-automorphisms)
        fa = fix_anti(ct_list, n)
        SC_n += cc * fa

    V_tourn //= nf
    V_graph //= nf
    V_even_all //= nf
    SC_n //= nf
    T_transition //= nf
    twin_SL = twin_SL_sum // nf

    # ===== LEVEL 3: Merged (÷ Z_2) =====
    V_merged = (V_tourn + SC_n) // 2

    # ===== LEVEL 1: Score sequences =====
    # Number of distinct score sequences of tournaments on n vertices
    # This is a known sequence but let's estimate
    # For n vertices, scores are (s_0,...,s_{n-1}) with sum = C(n,2), each 0 ≤ s_i ≤ n-1
    # Number of SORTED score sequences = partitions of C(n,2) into n parts from 0 to n-1
    # Known values: 1, 1, 2, 4, 9, 22, 59, 167 for n=1..8 (A000571)

    return {
        'n': n, 'm': m,
        # Level 5
        'labeled_tourn': labeled_tournaments,
        'labeled_graph': labeled_graphs,
        'labeled_even': labeled_even,
        'labeled_score': labeled_scores,
        # Level 4
        'V_tourn': V_tourn,
        'V_graph': V_graph,
        'V_even': V_even_all,
        'SC': SC_n,
        # Level 3
        'V_merged': V_merged,
        # Transitions
        'T': T_transition,
        'twin_SL': twin_SL,
        # Ratios
        'fiber_even': V_tourn / V_even_all if V_even_all > 0 else 0,
        'fiber_merged': V_tourn / V_merged if V_merged > 0 else 0,
        'graph_to_tourn': V_graph / V_tourn if V_tourn > 0 else 0,
    }

if __name__ == '__main__':
    print("=" * 80)
    print("  THE TOWER OF HIERARCHY WITH BURNSIDE")
    print("  opus-2026-03-23-S264")
    print("=" * 80)

    print("""
  LEVEL 5: Labeled objects           2^{C(n,2)} tournaments, 2^{C(n,2)} graphs
  LEVEL 4: Iso classes (÷ S_n)      V_n tournaments, A000088 graphs, A002854 even
  LEVEL 3: Merged (÷ Z_2)           V_merged = (V_n + SC_n)/2
  LEVEL 2: Even graph classes        V_even (cycle space quotient)
  LEVEL 1: Score sequences           cut space quotient

  FACTORIZATION: Fix_tourn = Fix_even × 2^{k-1}
  Each level strips away one type of symmetry.
""")

    results = []
    for n in range(3, 16):
        r = compute_tower(n)
        results.append(r)

    # The Tower Table
    print(f"\n{'='*80}")
    print("  THE TOWER — ALL LEVELS")
    print(f"{'='*80}")
    print(f"{'n':>3} {'2^m':>14} {'V_graph':>12} {'V_tourn':>12} {'V_merged':>12} "
          f"{'V_even':>10} {'SC':>8} {'G/T':>6} {'T/E':>8}")

    for r in results:
        print(f"{r['n']:>3} {r['labeled_tourn']:>14} {r['V_graph']:>12} "
              f"{r['V_tourn']:>12} {r['V_merged']:>12} {r['V_even']:>10} "
              f"{r['SC']:>8} {r['graph_to_tourn']:>6.2f} {r['fiber_even']:>8.2f}")

    # The Fibonacci-like growth
    print(f"\n{'='*80}")
    print("  GROWTH RATIOS")
    print(f"{'='*80}")
    print(f"{'n':>3} {'V_tourn':>12} {'T/T_{n-1}':>10} {'V_even':>10} {'E/E_{n-1}':>10} "
          f"{'T/E':>8} {'(T/E)/(T/E)_{n-1}':>18}")

    for i in range(len(results)):
        r = results[i]
        if i > 0:
            prev = results[i-1]
            t_ratio = r['V_tourn'] / prev['V_tourn'] if prev['V_tourn'] > 0 else 0
            e_ratio = r['V_even'] / prev['V_even'] if prev['V_even'] > 0 else 0
            te_ratio = r['fiber_even']
            te_prev = prev['fiber_even'] if prev['fiber_even'] > 0 else 1
            te_growth = te_ratio / te_prev
            print(f"{r['n']:>3} {r['V_tourn']:>12} {t_ratio:>10.2f} {r['V_even']:>10} "
                  f"{e_ratio:>10.2f} {r['fiber_even']:>8.2f} {te_growth:>18.4f}")
        else:
            print(f"{r['n']:>3} {r['V_tourn']:>12} {'—':>10} {r['V_even']:>10} "
                  f"{'—':>10} {r['fiber_even']:>8.2f} {'—':>18}")

    # The Burnside Cascade
    print(f"\n{'='*80}")
    print("  THE BURNSIDE CASCADE — FROM LABELED TO SCORE")
    print(f"{'='*80}")
    print("""
  At each level, a Burnside quotient removes one layer of symmetry:

  2^{C(n,2)}  →  V_tourn  →  V_merged  →  V_even  →  V_score
     ÷ S_n      ÷ Z_2       project      project
                             cycle sp.    cut sp.

  The EXPONENT at each level:
    Labeled: C(n,2) bits total
    Tournaments: arc_orbits(σ) per σ
    Even graphs: arc_orbits(σ) - (k-1) per σ
    Scores: k - 1 per σ

  So: C(n,2) → arc_orbits → (cycle_null, cut_free) → score_seq
  Each step removes degrees of freedom.
""")

    # The 2^{k-1} structure across cycle types
    print(f"\n{'='*80}")
    print("  THE k-1 SPECTRUM: how many cut-free bits per cycle type")
    print(f"{'='*80}")

    for n in [5, 7, 9]:
        print(f"\n  n = {n}:")
        print(f"  {'cycle_type':>20} {'k':>3} {'k-1':>4} {'ao':>4} {'cyc_null':>9} "
              f"{'2^(k-1)':>8} {'Fix_even':>10} {'Fix_tourn':>12}")
        nf = factorial(n)
        for ct in partitions(n):
            ct_list = list(ct)
            if any(c % 2 == 0 for c in ct_list): continue
            k = len(ct_list)
            ao = arc_orbits(ct_list)
            cn = ao - k + 1
            print(f"  {str(ct_list):>20} {k:>3} {k-1:>4} {ao:>4} {cn:>9} "
                  f"{2**(k-1):>8} {2**cn:>10} {2**ao:>12}")

    # The CD tower intersection
    print(f"\n{'='*80}")
    print("  CAYLEY-DICKSON TOWER × BURNSIDE TOWER")
    print(f"{'='*80}")
    cd_levels = [
        (2, "R (real)"),
        (3, "C (complex)"),
        (5, "H (quaternion)"),
        (9, "O (octonion)"),
    ]
    print(f"  {'n':>3} {'CD level':>15} {'V_tourn':>12} {'V_even':>10} {'T/E':>8} "
          f"{'V_merged':>10} {'SC':>6} {'SC/V':>6}")
    for n, name in cd_levels:
        r = compute_tower(n)
        sc_frac = r['SC'] / r['V_tourn'] if r['V_tourn'] > 0 else 0
        print(f"  {n:>3} {name:>15} {r['V_tourn']:>12} {r['V_even']:>10} "
              f"{r['fiber_even']:>8.2f} {r['V_merged']:>10} {r['SC']:>6} {sc_frac:>6.3f}")

    # The identity dominance
    print(f"\n{'='*80}")
    print("  IDENTITY DOMINANCE IN THE BURNSIDE SUM")
    print(f"{'='*80}")
    print(f"  {'n':>3} {'identity_fix':>14} {'total_fix':>14} {'identity%':>10}")
    for n in range(3, 13):
        nf = factorial(n)
        identity_fix = 2 ** comb(n, 2)
        total_fix = 0
        for ct in partitions(n):
            ct_list = list(ct)
            if any(c % 2 == 0 for c in ct_list): continue
            total_fix += ccs_val(n, ct_list) * (2 ** arc_orbits(ct_list))
        pct = 100.0 * identity_fix / total_fix if total_fix > 0 else 0
        print(f"  {n:>3} {identity_fix:>14} {total_fix:>14} {pct:>9.4f}%")

    # The factorization as a product of two towers
    print(f"\n{'='*80}")
    print("  THE FACTORIZATION AS PRODUCT OF TWO TOWERS")
    print(f"{'='*80}")
    print("""
  TOWER A (Cycle Space / Even Graphs):
    Exponent per σ: arc_orbits - k + 1 = cycle_nullity
    Counts even graph classes at odd n

  TOWER B (Cut Space / Scores):
    Exponent per σ: k - 1 = cut_free
    Counts score-space orbits

  TOWER A × TOWER B ≠ TOURNAMENT TOWER (the product is wrong!)
  Because the S_n action entangles cut and cycle spaces.
  But the BURNSIDE EXPONENT factors: ao = cn + (k-1).
  The entanglement is in the conjugacy class sizes, not the exponents.
""")

    # Verify: V_tourn ≠ V_even × V_score (product doesn't work)
    print(f"  {'n':>3} {'V_tourn':>10} {'V_even':>8} {'V_score':>8} {'V_e×V_s':>10} {'match':>6}")
    # V_score = number of score sequence iso classes (A000571)
    a000571 = {3:2, 4:4, 5:9, 6:22, 7:59, 8:167, 9:490}
    for n in range(3, 10):
        r = compute_tower(n)
        vs = a000571.get(n, '?')
        if isinstance(vs, int):
            product = r['V_even'] * vs
            match = "✓" if product == r['V_tourn'] else "✗"
            print(f"  {n:>3} {r['V_tourn']:>10} {r['V_even']:>8} {vs:>8} {product:>10} {match:>6}")

    print("\nDONE.")
