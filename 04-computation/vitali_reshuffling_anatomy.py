"""
vitali_reshuffling_anatomy.py -- kind-pasteur-2026-03-13-S61
Anatomy of the c3/c5 reshuffling under the Vitali atom at n=8.

KEY QUESTION: The overlap weight spectrum is preserved, but i2 changes.
How is this possible?

ANSWER: i2 counts pairs of DIRECTED cycles that are vertex-disjoint.
The overlap spectrum counts pairs of VERTEX SETS that overlap by a given amount.
These are DIFFERENT because multiple directed cycles can share the same vertex set.

When vertex set V_a is "lost" (has multiplicity 0 in T') and V_b is "gained"
(has multiplicity 0 in T), the overlap weight spectrum doesn't change if
V_a and V_b have the same size (which they do, since c3/c5/c7 counts preserved).
But the DIRECTED cycles that lived on V_a may have been disjoint from different
cycles than the ones that will live on V_b.

Actually, wait: at the vertex set level, the overlap spectrum IS preserved.
This means the number of disjoint VERTEX SET pairs is preserved.
But i2 counts disjoint DIRECTED CYCLE pairs, which includes multiplicity.

If V has m_A directed cycles and V' is disjoint from V with m'_A directed
cycles, then the number of disjoint pairs from V x V' is m_A * m'_A.
When multiplicities change, this product changes.

This is the HIDDEN DIMENSION: the independence polynomial depends not just
on the conflict graph structure (which vertex sets overlap) but on the
MULTIPLICITIES (how many directed cycles per vertex set).

TEST: Is delta_i2 explained by multiplicity changes on fixed vertex sets?
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def find_cycle_vertex_sets(A, n, length):
    sets = []
    for combo in combinations(range(n), length):
        for perm in permutations(combo[1:]):
            path = (combo[0],) + perm
            valid = True
            for k in range(length):
                if A[path[k]][path[(k+1) % length]] != 1:
                    valid = False
                    break
            if valid:
                sets.append(frozenset(combo))
                break
    return sets

def count_directed_on_set(A, vset):
    combo = tuple(sorted(vset))
    k = len(combo)
    count = 0
    for perm in permutations(combo[1:]):
        path = (combo[0],) + perm
        valid = True
        for i in range(k):
            if A[path[i]][path[(i+1) % k]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

n = 8
total_bits = n * (n - 1) // 2

print("=" * 70)
print("ANATOMY OF THE RESHUFFLING: WHY i2 CHANGES")
print("=" * 70)

np.random.seed(42)
examples = []

for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    lam = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(lam, lambda_graph(B, n)):
            continue

        H_A = count_ham_paths(A, n)
        H_B = count_ham_paths(B, n)
        if H_A == H_B:
            continue

        # Get all vertex sets and their multiplicities
        vset_mults_A = {}
        vset_mults_B = {}

        for length in [3, 5, 7]:
            vsets_A = find_cycle_vertex_sets(A, n, length)
            vsets_B = find_cycle_vertex_sets(B, n, length)
            all_vsets = set(map(frozenset, vsets_A)) | set(map(frozenset, vsets_B))
            for vs in all_vsets:
                m_A = count_directed_on_set(A, vs)
                m_B = count_directed_on_set(B, vs)
                if m_A > 0:
                    vset_mults_A[vs] = (length, m_A)
                if m_B > 0:
                    vset_mults_B[vs] = (length, m_B)

        # Compute i2 using DIRECTED cycles (vertex sets * multiplicity)
        def compute_i2_from_vsets(vset_mults):
            """Compute i2 = number of vertex-disjoint DIRECTED cycle pairs."""
            vsets = list(vset_mults.keys())
            total = 0
            for i in range(len(vsets)):
                for j in range(i+1, len(vsets)):
                    if len(vsets[i] & vsets[j]) == 0:
                        _, m_i = vset_mults[vsets[i]]
                        _, m_j = vset_mults[vsets[j]]
                        total += m_i * m_j
            # Also count pairs within same vertex set (but these are NEVER disjoint since same set)
            # So only cross-set pairs contribute
            return total

        i2_A = compute_i2_from_vsets(vset_mults_A)
        i2_B = compute_i2_from_vsets(vset_mults_B)
        delta_i2 = i2_B - i2_A

        # Also compute vertex-set-level disjoint pairs (ignoring multiplicity)
        def compute_vset_disjoint_pairs(vset_mults):
            vsets = list(vset_mults.keys())
            total = 0
            for i in range(len(vsets)):
                for j in range(i+1, len(vsets)):
                    if len(vsets[i] & vsets[j]) == 0:
                        total += 1
            return total

        vs_disjoint_A = compute_vset_disjoint_pairs(vset_mults_A)
        vs_disjoint_B = compute_vset_disjoint_pairs(vset_mults_B)
        delta_vs_disjoint = vs_disjoint_B - vs_disjoint_A

        # Lost and gained vertex sets
        lost_vs = set(vset_mults_A.keys()) - set(vset_mults_B.keys())
        gained_vs = set(vset_mults_B.keys()) - set(vset_mults_A.keys())
        common_vs = set(vset_mults_A.keys()) & set(vset_mults_B.keys())

        # Multiplicity changes on common vertex sets
        mult_delta = {}
        for vs in common_vs:
            l_A, m_A = vset_mults_A[vs]
            l_B, m_B = vset_mults_B[vs]
            if m_A != m_B:
                mult_delta[vs] = (l_A, m_A, m_B)

        examples.append({
            'delta_H': H_B - H_A,
            'delta_i2': delta_i2,
            'delta_vs_disjoint': delta_vs_disjoint,
            'i2_A': i2_A,
            'i2_B': i2_B,
            'lost': len(lost_vs),
            'gained': len(gained_vs),
            'common': len(common_vs),
            'mult_changes': len(mult_delta),
            'lost_vs': lost_vs,
            'gained_vs': gained_vs,
            'mult_delta': mult_delta,
            'vset_mults_A': vset_mults_A,
            'vset_mults_B': vset_mults_B,
        })
        break

    if trial % 400 == 0 and trial > 0:
        print(f"  ... {trial}/2000, found {len(examples)} examples")

print(f"\nTotal examples: {len(examples)}")

# Check: does delta_vs_disjoint = 0 always? (overlap spectrum preserved)
print(f"\n  delta_vs_disjoint distribution: {dict(sorted(Counter(e['delta_vs_disjoint'] for e in examples).items()))}")

# What DOES explain delta_i2?
# i2 = sum over disjoint VS pairs (i,j) of m_i * m_j
# delta_i2 = sum over disjoint VS pairs of (m_i' * m_j' - m_i * m_j)

# Decompose: for each disjoint VS pair, check if BOTH are common, one is common, or both changed
print(f"\n{'='*70}")
print(f"DECOMPOSITION OF delta_i2")
print(f"{'='*70}")

for idx, e in enumerate(examples[:5]):
    print(f"\n  Example {idx+1}: dH={e['delta_H']}, di2={e['delta_i2']}, dVS_disj={e['delta_vs_disjoint']}")
    print(f"    Lost VS: {e['lost']}, Gained VS: {e['gained']}, Common: {e['common']}, Mult changes: {e['mult_changes']}")

    # Compute delta_i2 contribution from different sources
    vA = e['vset_mults_A']
    vB = e['vset_mults_B']
    all_vs = set(vA.keys()) | set(vB.keys())

    # Group disjoint pairs by type
    delta_i2_both_common = 0
    delta_i2_one_changed = 0
    delta_i2_both_changed = 0

    all_vs_list = list(all_vs)
    for i in range(len(all_vs_list)):
        for j in range(i+1, len(all_vs_list)):
            vi, vj = all_vs_list[i], all_vs_list[j]
            if vi & vj:
                continue  # not disjoint

            mi_A = vA.get(vi, (0, 0))[1] if vi in vA else 0
            mi_B = vB.get(vi, (0, 0))[1] if vi in vB else 0
            mj_A = vA.get(vj, (0, 0))[1] if vj in vA else 0
            mj_B = vB.get(vj, (0, 0))[1] if vj in vB else 0

            delta = mi_B * mj_B - mi_A * mj_A

            if delta != 0:
                vi_changed = (mi_A != mi_B)
                vj_changed = (mj_A != mj_B)

                if vi_changed and vj_changed:
                    delta_i2_both_changed += delta
                elif vi_changed or vj_changed:
                    delta_i2_one_changed += delta
                else:
                    delta_i2_both_common += delta

    print(f"    delta_i2 breakdown:")
    print(f"      Both common, both unchanged: {delta_i2_both_common}")
    print(f"      One changed: {delta_i2_one_changed}")
    print(f"      Both changed: {delta_i2_both_changed}")
    print(f"      Total: {delta_i2_both_common + delta_i2_one_changed + delta_i2_both_changed} (should be {e['delta_i2']})")

    # Show the lost/gained vertex sets
    if e['lost'] <= 8:
        print(f"\n    Lost vertex sets:")
        for vs in sorted(e['lost_vs'], key=lambda x: (len(x), sorted(x))):
            l, m = vA[vs]
            print(f"      c{l} {sorted(vs)}: mult={m}")

        print(f"    Gained vertex sets:")
        for vs in sorted(e['gained_vs'], key=lambda x: (len(x), sorted(x))):
            l, m = vB[vs]
            print(f"      c{l} {sorted(vs)}: mult={m}")

    # Show multiplicity changes on common sets
    if e['mult_changes'] <= 15:
        print(f"\n    Multiplicity changes on common sets:")
        for vs, (l, mA, mB) in sorted(e['mult_delta'].items(), key=lambda x: (x[1][0], sorted(x[0]))):
            print(f"      c{l} {sorted(vs)}: {mA} -> {mB}")

# Summary statistics
print(f"\n{'='*70}")
print(f"SUMMARY STATISTICS")
print(f"{'='*70}")

# Is delta_vs_disjoint = 0 always?
all_zero = all(e['delta_vs_disjoint'] == 0 for e in examples)
print(f"\n  delta_vs_disjoint = 0 ALWAYS: {all_zero}")

# How often is delta_i2 driven by multiplicity changes vs vertex set swaps?
# Actually: since delta_vs_disjoint = 0, the disjoint pairs are the SAME vertex set pairs
# The change comes purely from multiplicity changes m_i' * m_j' != m_i * m_j

# Formula: delta_i2 = sum_{disjoint pairs (V_i,V_j)} (m_i'*m_j' - m_i*m_j)
# = sum_{disjoint pairs (V_i,V_j)} [(m_i'-m_i)*m_j' + m_i*(m_j'-m_j)]
# = sum_{disjoint pairs} delta_m_i * m_j' + m_i * delta_m_j

# This is a BILINEAR form in the multiplicity changes!

print(f"\n  For disjoint-pair-preserving changes (delta_vs_disjoint=0):")
print(f"  delta_i2 is a BILINEAR form in multiplicity changes")
print(f"  = sum over disjoint VS pairs of (delta_m_i * m_j' + m_i * delta_m_j)")

# What types of disjoint VS pairs exist at n=8?
# (c3, c5): need 3+5=8 = n vertices. Exactly n. EXISTS.
# (c3, c3): need 3+3=6 < n. EXISTS.
# (c3, c7): need 3+7=10 > n. IMPOSSIBLE.
# (c5, c5): need 5+5=10 > n. IMPOSSIBLE.
# (c5, c7): need 5+7=12 > n. IMPOSSIBLE.
# (c7, c7): need 7+7=14 > n. IMPOSSIBLE.

print(f"\n  Possible disjoint VS pair types at n=8:")
print(f"    (c3, c3): need 6 vertices, possible")
print(f"    (c3, c5): need 8 vertices = n, possible (tight!)")
print(f"    Others: impossible (need > 8 vertices)")

# Count disjoint pair types
for e in examples[:1]:
    vA = e['vset_mults_A']
    c35_pairs = 0
    c33_pairs = 0
    for vi in vA:
        for vj in vA:
            if vi >= vj:
                continue
            if vi & vj:
                continue
            li = vA[vi][0]
            lj = vA[vj][0]
            if (li, lj) == (3, 5) or (li, lj) == (5, 3):
                c35_pairs += 1
            elif li == 3 and lj == 3:
                c33_pairs += 1
    print(f"\n  Example 1 disjoint pair types:")
    print(f"    (c3, c3) pairs: {c33_pairs}")
    print(f"    (c3, c5) pairs: {c35_pairs}")

# Which pair type drives delta_i2?
print(f"\n{'='*70}")
print(f"WHICH DISJOINT PAIR TYPE DRIVES delta_i2?")
print(f"{'='*70}")

for idx, e in enumerate(examples[:10]):
    vA = e['vset_mults_A']
    vB = e['vset_mults_B']
    all_vs = set(vA.keys()) | set(vB.keys())
    all_vs_list = list(all_vs)

    delta_33 = 0
    delta_35 = 0

    for i in range(len(all_vs_list)):
        for j in range(i+1, len(all_vs_list)):
            vi, vj = all_vs_list[i], all_vs_list[j]
            if vi & vj:
                continue

            mi_A = vA.get(vi, (0, 0))[1] if vi in vA else 0
            mi_B = vB.get(vi, (0, 0))[1] if vi in vB else 0
            mj_A = vA.get(vj, (0, 0))[1] if vj in vA else 0
            mj_B = vB.get(vj, (0, 0))[1] if vj in vB else 0

            delta = mi_B * mj_B - mi_A * mj_A

            li = vA[vi][0] if vi in vA else vB[vi][0]
            lj = vA[vj][0] if vj in vA else vB[vj][0]

            pair_type = tuple(sorted([li, lj]))
            if pair_type == (3, 3):
                delta_33 += delta
            elif pair_type == (3, 5):
                delta_35 += delta

    print(f"  Ex {idx+1}: dH={e['delta_H']:>3}, di2={e['delta_i2']:>3}, d33={delta_33:>3}, d35={delta_35:>3}, sum={delta_33+delta_35:>3}")

print(f"\n{'='*70}")
print(f"INSIGHT: THE {'{'}2,1,0{'}'} OVERLAP WEIGHT STRUCTURE")
print(f"{'='*70}")
print("""
The three overlap weights W=2, W=1, W=0 play different roles:

W >= 1 (share vertex): CONFLICT in Omega. These pairs are adjacent.
  The Vitali atom preserves the NUMBER of such pairs (overlap spectrum preserved)
  but changes WHICH specific pairs are at each weight.

W = 0 (disjoint): INDEPENDENT in Omega. These contribute to i2.
  The Vitali atom preserves the number of disjoint VERTEX SET pairs
  but changes the MULTIPLICITIES on each vertex set.
  delta_i2 = sum over disjoint VS pairs of (delta_m_i * m_j' + m_i * delta_m_j)

At n=8, the only disjoint pair types are:
  (c3, c3): 3+3=6 <= 8. ALWAYS exists. Drives part of delta_i2.
  (c3, c5): 3+5=8 = n. Tight! Drives part of delta_i2.

The multiplicity delta_m on a vertex set V is determined by
how the Vitali reversal changes the INTERNAL TOURNAMENT on V.
The reversal changes arcs within the 4-subset S, which affects
the structure of any cycle that passes through S.

This is the HIDDEN HIGHER-DIMENSIONAL STRUCTURE:
- Level 0: Cycle counts (c3, c5, c7) — the "gross" statistics
- Level 1: Cycle identities (which vertex sets) — preserved in count but reshuffled
- Level 2: Directed multiplicities — how many directed cycles per vertex set
- Level 3: Disjoint pair multiplicities — products of multiplicities on disjoint sets

The Vitali atom acts at Level 2 (multiplicities) which propagates to Level 3 (i2).
Levels 0 and 1 are preserved (at n=8).
""")

print("Done.")
