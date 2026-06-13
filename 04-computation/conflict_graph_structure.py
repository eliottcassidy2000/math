#!/usr/bin/env python3
"""
conflict_graph_structure.py — opus-2026-03-13-S67k
Analyze the structure of odd-cycle conflict graphs in detail.

Key questions:
1. What is the chromatic number / clique cover number of CG?
2. How does CG decompose into cliques (cycles sharing vertex sets)?
3. What is the relationship between CG structure and Fibonacci?
4. Does the "two-channel" structure (THM-170) generalize?

THM-170 (kind-pasteur): delta_H = 2*dc7 + 4*di2 at n=8
Our HYP-818: H = I_{CG}(2)
These are the SAME observation from different angles.
"""

from itertools import permutations, combinations
from collections import defaultdict

def tournament_from_bits(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def canonical_form(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(i+1, n))
        if best is None or form < best:
            best = form
    return best

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def find_all_directed_odd_cycles(A, n):
    """Find all directed odd cycles, returning (normalized_cycle, vertex_set)."""
    cycles = []
    for length in range(3, n+1, 2):
        for combo in combinations(range(n), length):
            for perm in permutations(combo):
                is_cycle = True
                for k in range(length):
                    if not A[perm[k]][perm[(k+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    min_idx = perm.index(min(perm))
                    normalized = tuple(perm[min_idx:] + perm[:min_idx])
                    cycles.append((normalized, frozenset(combo)))
    # Deduplicate
    seen = set()
    unique = []
    for c, vs in cycles:
        if c not in seen:
            seen.add(c)
            unique.append((c, vs))
    return unique

print("=" * 70)
print("CONFLICT GRAPH DEEP STRUCTURE")
print("=" * 70)

for n in range(3, 7):
    m = n*(n-1)//2
    classes = defaultdict(list)
    for bits in range(1 << m):
        A = tournament_from_bits(n, bits)
        cf = canonical_form(A, n)
        classes[cf].append(bits)

    iso_classes = sorted(classes.keys())
    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")

    for cf in iso_classes:
        A = tournament_from_bits(n, classes[cf][0])
        H = count_ham_paths(A, n)
        cycles = find_all_directed_odd_cycles(A, n)
        nc = len(cycles)

        if nc == 0:
            continue
        if nc > 25:
            # Skip very large conflict graphs
            print(f"\n  H={H}: {nc} cycles (too large for detailed analysis)")
            continue

        # Group cycles by vertex set
        vs_groups = defaultdict(list)
        for c, vs in cycles:
            vs_groups[vs].append(c)

        # The clique structure: cycles on same vertex set form a clique
        print(f"\n  H={H}: {nc} directed odd cycles on {len(vs_groups)} vertex sets")

        # Count cycles by length
        by_length = defaultdict(int)
        for c, vs in cycles:
            by_length[len(c)] += 1
        print(f"    By length: {dict(sorted(by_length.items()))}")

        # Count cycles per vertex set
        vs_sizes = defaultdict(list)
        for vs, cycs in vs_groups.items():
            vs_sizes[len(vs)].append(len(cycs))
        for vlen in sorted(vs_sizes.keys()):
            mults = sorted(vs_sizes[vlen])
            print(f"    {vlen}-vertex sets: {len(mults)} sets, cycles/set = {mults}")

        # Build conflict graph edges by vertex set overlap
        cycle_list = list(cycles)
        vertex_sets = [vs for c, vs in cycle_list]

        # Key metric: which vertex sets are disjoint?
        disjoint_pairs = 0
        disjoint_vs_pairs = set()
        for i in range(len(vertex_sets)):
            for j in range(i+1, len(vertex_sets)):
                if not vertex_sets[i] & vertex_sets[j]:
                    disjoint_pairs += 1
                    vs_pair = frozenset([vertex_sets[i], vertex_sets[j]])
                    disjoint_vs_pairs.add(vs_pair)

        if disjoint_pairs > 0:
            print(f"    *** DISJOINT CYCLE PAIRS: {disjoint_pairs} (from {len(disjoint_vs_pairs)} disjoint VS pairs)")
            for dvp in sorted(disjoint_vs_pairs, key=lambda x: str(sorted(x))):
                vsets = list(dvp)
                sizes = sorted([len(vs) for vs in vsets])
                counts = [len(vs_groups[vs]) for vs in vsets]
                total_pairs = counts[0] * counts[1]
                print(f"      {set(vsets[0])} ∪ {set(vsets[1])} = sizes {sizes}, "
                      f"cycles = {counts[0]}×{counts[1]} = {total_pairs} disjoint pairs")

# Specific deep analysis for the n=6 H=45 classes
print(f"\n{'='*70}")
print(f"DEEP DIVE: n=6 CLASSES WITH α₂ > 0")
print(f"{'='*70}")

n = 6
m = n*(n-1)//2
classes6 = defaultdict(list)
for bits in range(1 << m):
    A = tournament_from_bits(n, bits)
    cf = canonical_form(A, n)
    classes6[cf].append(bits)

iso6 = sorted(classes6.keys())

alpha2_classes = []
for cf in iso6:
    A = tournament_from_bits(n, classes6[cf][0])
    H = count_ham_paths(A, n)
    cycles = find_all_directed_odd_cycles(A, n)
    vertex_sets = [vs for c, vs in cycles]

    # Count disjoint pairs
    dp = 0
    for i in range(len(vertex_sets)):
        for j in range(i+1, len(vertex_sets)):
            if not vertex_sets[i] & vertex_sets[j]:
                dp += 1
    if dp > 0:
        alpha2_classes.append((H, len(cycles), dp, cf))

print(f"\n{len(alpha2_classes)} classes with α₂ > 0:")
print(f"{'H':>4} {'#cyc':>5} {'α₂':>4}")
for H, nc, a2, cf in sorted(alpha2_classes):
    print(f"{H:4d} {nc:5d} {a2:4d}")

# Distribution of α₂
a2_vals = [a2 for H, nc, a2, cf in alpha2_classes]
print(f"\nα₂ distribution: {sorted(set(a2_vals))} with counts {[a2_vals.count(v) for v in sorted(set(a2_vals))]}")

# Key finding: what determines α₂?
print(f"\nα₂ = number of vertex-disjoint directed cycle pairs")
print(f"At n=6, disjoint pairs can only be (3-cycle, 3-cycle) since 3+3=6=n")
print(f"Two 3-cycles on disjoint vertex sets {'{a,b,c}'} and {'{d,e,f}'}")
print(f"Each vertex set has 0 or 1 directed 3-cycle (in a tournament)")
print()

# Check: for each class with α₂>0, what are the disjoint 3-cycle pairs?
for H, nc, a2, cf in sorted(alpha2_classes):
    A = tournament_from_bits(n, classes6[cf][0])
    cycles = find_all_directed_odd_cycles(A, n)

    # Find 3-cycles
    c3 = [(c, vs) for c, vs in cycles if len(c) == 3]

    # Find disjoint 3-cycle pairs
    disjoint_3c_pairs = []
    for i in range(len(c3)):
        for j in range(i+1, len(c3)):
            if not c3[i][1] & c3[j][1]:
                disjoint_3c_pairs.append((c3[i], c3[j]))

    # Also check 3+5 disjoint (impossible at n=6: 3+5=8>6)
    c5 = [(c, vs) for c, vs in cycles if len(c) == 5]
    disjoint_35_pairs = []
    for a in c3:
        for b in c5:
            if not a[1] & b[1]:
                disjoint_35_pairs.append((a, b))

    print(f"  H={H}: {len(c3)} 3-cycles, {len(c5)} 5-cycles, "
          f"disjoint (3,3) pairs: {len(disjoint_3c_pairs)}, "
          f"disjoint (3,5) pairs: {len(disjoint_35_pairs)}")

print(f"\n{'='*70}")
print(f"THE TWO-CHANNEL FORMULA GENERALIZED")
print(f"{'='*70}")
print("""
THM-170 (kind-pasteur) at n=8:
  delta_H = 2 * delta_c7_dir + 4 * delta_i2

Our framework (HYP-818):
  H = I_{CG}(2) = 1 + 2*α_1 + 4*α_2 + 8*α_3 + ...
  where α_k = independent sets of size k in CG

UNIFICATION:
  At n=8: α_1 = c3 + c5 + c7 (directed counts)
           α_2 = disjoint directed cycle pairs
           α_3 = 0 (need 9+ vertices)

  Under (1,1,2,2) reversal preserving lambda:
    delta_c3 = delta_c5 = 0 (THM-170 step 3)
    So delta_α_1 = delta_c7
    And delta_H = 2*delta_c7 + 4*delta_α_2 ✓

GENERAL FORMULA FOR ARBITRARY n:
  H = 1 + 2*(c3+c5+...+c_n_odd) + 4*α_2 + 8*α_3 + ...

  where c_{2k+1} = number of directed (2k+1)-cycles
  and α_j = number of independent sets of size j in CG
         = vertex-disjoint j-tuples of directed odd cycles

  The "channels" are:
    Channel 0: constant 1
    Channel 1: individual odd cycles (coefficient 2)
    Channel 2: disjoint pairs (coefficient 4)
    Channel k: disjoint k-tuples (coefficient 2^k)

  This is a MULTI-CHANNEL decomposition where each channel
  contributes independently to H.

FIBONACCI CONNECTION (refined):
  The independence polynomial I_G(x) satisfies:
    If G = P_m (path): I(x) = Fibonacci polynomial F_{m+2}(x)
    If G = K_m (complete): I(x) = 1 + mx (only channel 1)
    If G = mK_1 (independent): I(x) = (1+x)^m (all channels)

  Tournament conflict graphs are between K_m and mK_1:
    - Mostly complete (most cycles share vertices)
    - A few disjoint pairs (α_2 small relative to α_1)
    - Very few disjoint triples (α_3 tiny)

  The "Fibonacci" structure emerges because the conflict graph
  is CLOSE to a path graph when viewed at the vertex-set level
  (each vertex set can conflict with the next in a chain).
""")
