#!/usr/bin/env python3
"""
opus-2026-04-05-S28: Design structure of directed 5-cycles in tournaments.

We know: 3-cycles of Paley P_p form a 2-(p, 3, (p+1)/4) BIBD.
Question: Do 5-cycles (and general odd k-cycles) also form designs?

For Paley: Aut(P_p) ≅ AGL(1,p) is 2-transitive → k-cycles automatically
form 2-designs. But:
  - What are the parameters λ_k?
  - Do they form t-designs for t ≥ 3?
  - What is the structure for NON-Paley tournaments?
  - Is there a formula for λ_k in terms of QR structure?

Key definitions:
  - A directed k-cycle: (v₁,v₂,...,vₖ) with v₁→v₂→...→vₖ→v₁, up to cyclic shift
  - The "block" of a directed k-cycle: the unordered k-subset {v₁,...,vₖ}
  - A 5-subset can support 0, 1, or more directed 5-cycles
  - For the design: each directed 5-cycle = one block (even if same 5-subset)
"""

from math import factorial, comb, gcd
from itertools import permutations, combinations
from collections import Counter, defaultdict
import time

# ═══════════════════════════════════════════════════════════════════
# Tournament construction
# ═══════════════════════════════════════════════════════════════════

def paley(p):
    """Paley tournament on p vertices (p ≡ 3 mod 4 prime)."""
    qr = {(a*a) % p for a in range(1, p)}
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                T[i][j] = 1
    return T

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i*i <= n:
        if n%i==0 or n%(i+2)==0: return False
        i += 6
    return True

# ═══════════════════════════════════════════════════════════════════
# Cycle counting
# ═══════════════════════════════════════════════════════════════════

def directed_k_cycles(T, k):
    """Find all directed k-cycles in tournament T.
    Returns list of tuples (v₁,...,vₖ) in canonical form (min vertex first)."""
    n = len(T)
    cycles = set()

    for subset in combinations(range(n), k):
        # Find Hamiltonian cycles in subtournament T[subset]
        verts = list(subset)
        # Fix first vertex to break cyclic symmetry
        v0 = verts[0]
        for perm in permutations(verts[1:]):
            path = [v0] + list(perm)
            ok = True
            for i in range(k):
                if T[path[i]][path[(i+1) % k]] != 1:
                    ok = False
                    break
            if ok:
                # Canonical form: rotate so minimum is first
                min_idx = path.index(min(path))
                canon = tuple(path[min_idx:] + path[:min_idx])
                cycles.add(canon)

    return list(cycles)

def pair_incidence(cycles, n, k):
    """For each unordered pair {u,v}, count how many k-cycles contain both."""
    counts = {}
    for i in range(n):
        for j in range(i+1, n):
            counts[(i,j)] = 0

    for cyc in cycles:
        verts = set(cyc)
        for i in range(len(cyc)):
            for j in range(i+1, len(cyc)):
                a, b = min(cyc[i], cyc[j]), max(cyc[i], cyc[j])
                counts[(a,b)] += 1

    return counts

def triple_incidence(cycles, n):
    """For each unordered triple {u,v,w}, count how many cycles contain all three."""
    counts = defaultdict(int)
    for cyc in cycles:
        verts = sorted(cyc)
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                for l in range(j+1, len(verts)):
                    counts[(verts[i], verts[j], verts[l])] += 1
    return dict(counts)

def quad_incidence(cycles, n):
    """For each unordered 4-subset, count how many cycles contain all four."""
    counts = defaultdict(int)
    for cyc in cycles:
        verts = sorted(cyc)
        for combo in combinations(verts, 4):
            counts[combo] += 1
    return dict(counts)

def check_t_design(incidence_counts, t_name):
    """Check if incidence counts are uniform (= t-design)."""
    vals = list(incidence_counts.values())
    if not vals:
        return False, 0, []
    unique = sorted(set(vals))
    is_design = len(unique) == 1
    return is_design, unique[0] if is_design else None, unique

# ═══════════════════════════════════════════════════════════════════
# Main analysis
# ═══════════════════════════════════════════════════════════════════

def analyze_tournament(T, name, max_k=None):
    """Analyze k-cycle designs for a tournament."""
    n = len(T)
    if max_k is None:
        max_k = min(n, 9)

    print(f"\n{'='*70}")
    print(f"TOURNAMENT: {name} (n={n})")
    print(f"{'='*70}")

    results = {}

    for k in range(3, max_k + 1, 2):  # odd k only
        if k > n:
            break

        t_start = time.time()
        cycles = directed_k_cycles(T, k)
        c_k = len(cycles)
        elapsed = time.time() - t_start

        print(f"\n  k={k}: {c_k} directed {k}-cycles (computed in {elapsed:.2f}s)")

        if c_k == 0:
            print(f"    No {k}-cycles found.")
            continue

        # 2-design check (pair incidence)
        pair_counts = pair_incidence(cycles, n, k)
        is_2, lam_2, unique_2 = check_t_design(pair_counts, "pair")

        print(f"    2-design: {'YES' if is_2 else 'NO'}", end="")
        if is_2:
            print(f", λ = {lam_2}")
            # Verify: c_k × C(k,2) / C(n,2) should = λ
            expected = c_k * comb(k, 2) / comb(n, 2)
            print(f"    Verify: c_{k}×C({k},2)/C({n},2) = {c_k}×{comb(k,2)}/{comb(n,2)} = {expected:.1f}")
        else:
            print(f", pair distribution: {Counter(pair_counts.values())}")

        # 3-design check (triple incidence)
        if k >= 5:  # only meaningful for k ≥ 5 (need triples in blocks of size k)
            triple_counts = triple_incidence(cycles, n)
            # Only check triples that could appear (i.e., within k-cycles)
            # For a t-design, ALL t-subsets must have the same count
            all_triples = {}
            for i in range(n):
                for j in range(i+1, n):
                    for l in range(j+1, n):
                        all_triples[(i,j,l)] = triple_counts.get((i,j,l), 0)

            is_3, lam_3, unique_3 = check_t_design(all_triples, "triple")
            print(f"    3-design: {'YES' if is_3 else 'NO'}", end="")
            if is_3:
                print(f", λ = {lam_3}")
            else:
                print(f", triple distribution: {Counter(all_triples.values())}")

        # 4-design check
        if k >= 5 and n <= 23:
            quad_counts = quad_incidence(cycles, n)
            all_quads = {}
            for combo in combinations(range(n), 4):
                all_quads[combo] = quad_counts.get(combo, 0)

            is_4, lam_4, unique_4 = check_t_design(all_quads, "quad")
            print(f"    4-design: {'YES' if is_4 else 'NO'}", end="")
            if is_4:
                print(f", λ = {lam_4}")
            else:
                dist = Counter(all_quads.values())
                print(f", quad distribution: {dist}")

        results[k] = {
            'count': c_k, 'is_2': is_2, 'lambda_2': lam_2,
            'pair_dist': Counter(pair_counts.values())
        }

    return results

def formula_exploration():
    """Explore the formula for λ_k across primes."""
    print(f"\n{'='*70}")
    print(f"DESIGN PARAMETER TABLE: λ_k for Paley tournaments")
    print(f"{'='*70}")
    print()

    all_results = {}

    for p in [7, 11, 19, 23]:
        if not (is_prime(p) and p % 4 == 3):
            continue
        T = paley(p)

        print(f"\n--- P_{p} ---")
        max_k = min(p, 7) if p <= 11 else 5
        if p >= 19:
            max_k = 5
        results = analyze_tournament(T, f"Paley P_{p}", max_k=max_k)
        all_results[p] = results

    # Summary table
    print(f"\n{'='*70}")
    print(f"SUMMARY: Design parameters λ_k for Paley P_p")
    print(f"{'='*70}")
    print()
    print(f"{'p':>4} {'k':>3} {'c_k':>10} {'2-design':>9} {'λ_2':>8} {'3-design':>9} {'formula?':>12}")
    print("-" * 60)

    for p in sorted(all_results.keys()):
        for k in sorted(all_results[p].keys()):
            r = all_results[p][k]
            is2 = "YES" if r['is_2'] else "no"
            lam2 = str(r['lambda_2']) if r['lambda_2'] is not None else "-"

            # Try formula: λ_k(p) = c_k × C(k,2) / C(p,2)
            formula = ""
            if r['lambda_2'] is not None:
                # Check if λ/(p-k)! relates to something nice
                pass

            print(f"{p:>4} {k:>3} {r['count']:>10} {is2:>9} {lam2:>8}")

def non_paley_comparison():
    """Compare 5-cycle design properties for Paley vs non-Paley tournaments."""
    print(f"\n{'='*70}")
    print(f"COMPARISON: 5-cycle designs in Paley vs non-Paley at n=7")
    print(f"{'='*70}")
    print()

    n = 7
    T_paley = paley(7)

    # Also construct the "rotational" tournament (circulant with S={1,2,3})
    # and the "interval" tournament (circulant with S={1,2,3} — same for n=7)
    T_interval = [[0]*n for _ in range(n)]
    for i in range(n):
        for d in [1, 2, 3]:  # out-neighbors at distance 1,2,3
            j = (i + d) % n
            T_interval[i][j] = 1

    # Check if interval = Paley at n=7
    # QR_7 = {1, 2, 4}, interval uses {1, 2, 3}
    print("P_7 uses QR_7 = {1, 2, 4}")
    print("Interval uses S = {1, 2, 3}")
    print()

    analyze_tournament(T_paley, "Paley P_7", max_k=7)
    analyze_tournament(T_interval, "Interval [1,2,3] mod 7", max_k=7)

def deep_structure():
    """Deeper analysis of the 5-cycle design at P_7."""
    print(f"\n{'='*70}")
    print(f"DEEP STRUCTURE: 5-cycles in P_7")
    print(f"{'='*70}")
    print()

    T = paley(7)
    cycles_5 = directed_k_cycles(T, 5)
    cycles_3 = directed_k_cycles(T, 3)
    cycles_7 = directed_k_cycles(T, 7)

    print(f"P_7 cycle census:")
    print(f"  3-cycles: {len(cycles_3)}")
    print(f"  5-cycles: {len(cycles_5)}")
    print(f"  7-cycles (Hamiltonian): {len(cycles_7)}")

    # For each 5-cycle, which pair of vertices is NOT in it?
    print(f"\n5-cycle complements (the 2-element set NOT in each 5-cycle):")
    complement_counts = Counter()
    for cyc in cycles_5:
        missing = tuple(sorted(set(range(7)) - set(cyc)))
        complement_counts[missing] += 1

    print(f"  Complement pair distribution: {dict(complement_counts)}")
    is_uniform = len(set(complement_counts.values())) == 1
    print(f"  Uniform: {is_uniform}")
    if is_uniform:
        print(f"  Each complement pair appears in {list(complement_counts.values())[0]} 5-cycles")
        print(f"  → The 5-cycles form a COMPLEMENTARY design:")
        print(f"    The COMPLEMENTS of 5-blocks are a 2-(7,2,λ') design")

    # Connection between 3-cycles and 5-cycles
    print(f"\nRelation between 3-cycles and 5-cycles:")
    print(f"  A 5-cycle on {'{'}a,b,c,d,e{'}'} avoids pair {'{'}f,g{'}'}.")
    print(f"  The pair {'{'}f,g{'}'} forms an arc f→g or g→f.")
    print(f"  Question: does the number of 5-cycles avoiding {'{'}f,g{'}'}")
    print(f"  relate to the number of 3-cycles THROUGH {'{'}f,g{'}'}?")

    pair_3 = pair_incidence(cycles_3, 7, 3)
    pair_5 = pair_incidence(cycles_5, 7, 5)

    print(f"\n  {'pair':>8} {'3-cyc thru':>11} {'5-cyc thru':>11} {'5-cyc avoid':>12} {'sum':>5}")
    print(f"  {'-'*50}")
    for i in range(7):
        for j in range(i+1, 7):
            c3 = pair_3[(i,j)]
            c5 = pair_5[(i,j)]
            avoid_5 = len(cycles_5) - c5
            # complement_counts gives 5-cycles NOT containing this pair
            print(f"  ({i},{j}){' ':>3} {c3:>11} {c5:>11} {avoid_5:>12} {c3+c5:>5}")

    # The OCF connection
    print(f"\n  OBSERVATION: c₃(pair) + c₅(pair) = ? Check if constant")
    sums = [pair_3[(i,j)] + pair_5[(i,j)] for i in range(7) for j in range(i+1, 7)]
    print(f"  Sums: {sorted(set(sums))} → {'CONSTANT' if len(set(sums))==1 else 'NOT constant'}")

    if len(set(sums)) == 1:
        print(f"  λ₃ + λ₅ = {sums[0]} = constant across all pairs!")
        print(f"  This means: odd-cycle incidence is CONSERVED.")
        print(f"  A pair in more 3-cycles has correspondingly fewer 5-cycles.")

# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    t0 = time.time()

    analyze_tournament(paley(7), "Paley P_7", max_k=7)
    deep_structure()
    non_paley_comparison()
    formula_exploration()

    print(f"\nTotal time: {time.time()-t0:.1f}s")
