"""
vitali_c5_multiplicity_mechanism.py -- kind-pasteur-2026-03-13-S61
WHY does the Vitali atom change c5 multiplicities?

The reversal flips all 6 arcs within the 4-subset S = {s1,s2,s3,s4}.
A 5-cycle on vertex set V = {v1,v2,v3,v4,v5} is affected iff V
contains at least 2 members of S (since the reversal only changes
arcs between S members).

At n=8 with |S|=4, a 5-vertex set V intersects S in:
- |V cap S| = 0: impossible (|V|+|S| = 9 > 8)
  Wait: 5+4=9 > 8, so |V cap S| >= 1. Actually: |V cap S| >= |V|+|S|-n = 1.
- |V cap S| = 1: V has 1 S-vertex + 4 ext-vertices (4 ext vertices exist)
  Only arcs within S that touch v are reversed. But a 5-cycle through V
  uses arcs between consecutive vertices. If only 1 is in S, at most 2 arcs
  in the cycle touch the S-vertex, and those arcs go between S-vertex and
  ext-vertices, which are NOT reversed.
  So |V cap S| = 1 => cycle is INVARIANT.
- |V cap S| = 2: V has 2 S-vertices. The single arc between them IS reversed.
  This can change whether a specific directed cycle exists.
- |V cap S| = 3: V has 3 S-vertices. Three arcs between them are reversed.
- |V cap S| = 4: V = S + 1 ext vertex. All 6 internal arcs reversed.

KEY INSIGHT: |V cap S| = 2 is the MARGINAL case. One reversed arc.
This can create or destroy exactly ONE directed cycle.

For disjoint (c3, c5) pairs at n=8:
- c3 uses 3 vertices, c5 uses 5 vertices, together = 8 = all vertices
- c3 cap S has size >= 3+4-8 = -1, so cap >= 0. Actually cap >= max(0, 3+4-8) = 0.
  And cap <= min(3,4) = 3.
- c5 cap S has size >= 5+4-8 = 1. And cap <= 4.
- Since c3 and c5 are disjoint: c3 cap S + c5 cap S = |S| = 4.
- So if c3 cap S = k, then c5 cap S = 4-k.
  k=0: c5 cap S = 4. V = S + 1 ext. All arcs reversed. BIG change.
  k=1: c5 cap S = 3. Three arcs reversed. MEDIUM change.
  k=2: c5 cap S = 2. One arc reversed. SMALL change.
  k=3: c5 cap S = 1. Zero reversed arcs in c5! INVARIANT.

For the c5 to change multiplicity, we need c5 cap S >= 2, i.e., k <= 2.
For the c3 to exist and be disjoint, the 3 vertices of c3 include k S-vertices
and 3-k ext-vertices.

This gives a beautiful DECOMPOSITION of the disjoint (c3,c5) pair changes:
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

n = 8
total_bits = n * (n - 1) // 2

print("=" * 70)
print("c5 MULTIPLICITY CHANGE MECHANISM")
print("=" * 70)

np.random.seed(314)
data = []

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

        S = frozenset(subset)
        ext = frozenset(range(n)) - S

        # Get c5 vertex sets and check multiplicity changes
        c5_A = find_cycle_vertex_sets(A, n, 5)
        c5_B = find_cycle_vertex_sets(B, n, 5)
        all_c5 = set(map(frozenset, c5_A)) | set(map(frozenset, c5_B))

        c5_changes = []
        for vs in all_c5:
            m_A = count_directed_on_set(A, vs)
            m_B = count_directed_on_set(B, vs)
            if m_A != m_B:
                overlap_S = len(vs & S)
                c5_changes.append({
                    'vs': vs,
                    'overlap_S': overlap_S,
                    'm_A': m_A,
                    'm_B': m_B,
                    'delta': m_B - m_A,
                })

        # Get c3 vertex sets
        c3_A = find_cycle_vertex_sets(A, n, 3)
        c3_B = find_cycle_vertex_sets(B, n, 3)
        lost_c3 = set(map(frozenset, c3_A)) - set(map(frozenset, c3_B))
        gained_c3 = set(map(frozenset, c3_B)) - set(map(frozenset, c3_A))

        # For each lost/gained c3, compute overlap with S
        c3_changes = []
        for vs in lost_c3:
            c3_changes.append({'vs': vs, 'overlap_S': len(vs & S), 'type': 'lost'})
        for vs in gained_c3:
            c3_changes.append({'vs': vs, 'overlap_S': len(vs & S), 'type': 'gained'})

        data.append({
            'delta_H': H_B - H_A,
            'c5_changes': c5_changes,
            'c3_changes': c3_changes,
            'subset': subset,
            'lost_c3': lost_c3,
            'gained_c3': gained_c3,
        })
        break

    if trial % 500 == 0 and trial > 0:
        print(f"  ... {trial}/2000, found {len(data)} examples")

print(f"\nTotal examples: {len(data)}")

# Analysis 1: c5 overlap with S distribution
print(f"\n{'='*70}")
print(f"c5 MULTIPLICITY CHANGES BY |V cap S|")
print(f"{'='*70}")

overlap_counts = Counter()
overlap_deltas = {}

for d in data:
    for ch in d['c5_changes']:
        ov = ch['overlap_S']
        overlap_counts[ov] += 1
        if ov not in overlap_deltas:
            overlap_deltas[ov] = Counter()
        overlap_deltas[ov][ch['delta']] += 1

for ov in sorted(overlap_counts.keys()):
    print(f"\n  |V cap S| = {ov}: {overlap_counts[ov]} changes")
    print(f"    delta distribution: {dict(sorted(overlap_deltas[ov].items()))}")

# Analysis 2: c3 changes by overlap with S
print(f"\n{'='*70}")
print(f"c3 VERTEX SET CHANGES BY |V cap S|")
print(f"{'='*70}")

c3_ov_counts = Counter()
for d in data:
    for ch in d['c3_changes']:
        c3_ov_counts[(ch['type'], ch['overlap_S'])] += 1

for key in sorted(c3_ov_counts.keys()):
    print(f"  {key[0]}, |V cap S| = {key[1]}: {c3_ov_counts[key]}")

# Analysis 3: The c3+c5=8 constraint for disjoint pairs
print(f"\n{'='*70}")
print(f"DISJOINT (c3,c5) PAIRS AND THE S-OVERLAP CONSTRAINT")
print(f"{'='*70}")
print(f"""
For disjoint c3 and c5: |c3 cap S| + |c5 cap S| = 4
  k=0, 4-k=4: c3 fully external, c5 = S + 1 ext. All 6 arcs reversed in c5.
  k=1, 4-k=3: c3 has 1 S-vert + 2 ext. c5 has 3 S-verts.
  k=2, 4-k=2: c3 has 2 S-verts + 1 ext. c5 has 2 S-verts.
  k=3, 4-k=1: c3 has 3 S-verts. c5 has 1 S-vert => c5 INVARIANT.

Only k=0,1,2 can cause c5 multiplicity change (need |c5 cap S| >= 2).
k=3 is the "protected" case: c3 absorbs most of S, protecting c5.
""")

# Verify: do we ever see c5 cap S = 1 changing?
c5_cap1_changes = sum(1 for d in data for ch in d['c5_changes'] if ch['overlap_S'] == 1)
print(f"  c5 with |V cap S| = 1 changing: {c5_cap1_changes} (should be 0)")

# Analysis 4: How the c3 swap affects disjoint (c3,c5) pair products
print(f"\n{'='*70}")
print(f"HOW THE c3 SWAP AFFECTS DISJOINT PAIR PRODUCTS")
print(f"{'='*70}")

for idx, d in enumerate(data[:5]):
    print(f"\n  Example {idx+1} (dH={d['delta_H']}):")
    print(f"    Subset S = {d['subset']}")

    # For each lost c3, find disjoint c5 sets and their multiplicity changes
    for c3_vs in sorted(d['lost_c3'], key=sorted):
        c5_disjoint = [ch for ch in d['c5_changes'] if len(ch['vs'] & c3_vs) == 0]
        if c5_disjoint:
            print(f"    Lost c3 {sorted(c3_vs)} (cap S = {len(c3_vs & frozenset(d['subset']))}):")
            for ch in c5_disjoint:
                print(f"      Disjoint c5 {sorted(ch['vs'])}: mult {ch['m_A']}->{ch['m_B']} (delta={ch['delta']})")

    for c3_vs in sorted(d['gained_c3'], key=sorted):
        c5_disjoint = [ch for ch in d['c5_changes'] if len(ch['vs'] & c3_vs) == 0]
        if c5_disjoint:
            print(f"    Gained c3 {sorted(c3_vs)} (cap S = {len(c3_vs & frozenset(d['subset']))}):")
            for ch in c5_disjoint:
                print(f"      Disjoint c5 {sorted(ch['vs'])}: mult {ch['m_A']}->{ch['m_B']} (delta={ch['delta']})")

# Analysis 5: The complete picture
print(f"\n{'='*70}")
print(f"THE COMPLETE PICTURE: VITALI SET ANALOGY")
print(f"{'='*70}")
print("""
The Vitali set V in measure theory is constructed by choosing one
representative from each equivalence class of R/Q. The key property:
V is non-measurable because Q-translations of V tile R without overlap,
but countably many disjoint translates can't cover a positive measure set.

The Vitali ATOM in tournaments plays an analogous role:
- The (1,1,2,2) reversal is like a "translation" of the tournament structure
- It preserves the "measure-like" quantities (lambda graph, cycle counts)
- But it changes the "non-measurable" quantity (Hamiltonian path count H)

The {2,1,0} overlap structure IS the measure-theoretic structure:
- W=2 (share edge): deterministic, preserved (like measurable sets)
- W=1 (share vertex): conflict boundary (like the Lebesgue measurable sets)
- W=0 (disjoint): independence (like the non-measurable sets)

The Vitali atom acts at the boundary between W=1 and W=0:
- It preserves the STATISTICAL distribution of overlaps
  (like a measure-preserving transformation)
- But it changes the SPECIFIC identities of cycles
  (like a non-measurable set that looks the same statistically)
- The i2 change is the "non-measurable" part:
  it depends on the EXACT placement of cycles, not just their statistics

At n=7, the Vitali atom can only change cycle COUNTS (c7), which is
"measurable" — you can detect it from gross statistics.

At n=8, the Vitali atom acts on the HIDDEN DIMENSION:
- It changes c5 multiplicities (how many directed cycles per vertex set)
- These multiplicities determine i2 via the bilinear form
  i2 = sum over disjoint (VS_i, VS_j) of m_i * m_j
- The change is "non-measurable" in the sense that:
  total directed counts are preserved (delta_c3_dir = delta_c5_dir = 0)
  overlap spectrum is preserved (delta_vs_disjoint = 0)
  but the PRODUCT structure (bilinear form) changes

This is a PERFECT analogy to Vitali non-measurability:
the set "looks the same" under every statistical test (measure zero/one),
but it has structure that is invisible to the sigma-algebra.
""")

print("Done.")
