"""
vitali_dc5_zero_proof.py -- kind-pasteur-2026-03-13-S61
WHY is dc5 = 0 for all n? (Verified n=8,9,10)

The Vitali atom reverses S = {s1,s2,s3,s4} with score (1,1,2,2).
A 5-cycle on vertex set V is affected iff |V ∩ S| >= 2.

Key: the reversal on S is the COMPLEMENT operation on the 4-vertex
tournament T[S]. If T[S] has score (1,1,2,2), then T'[S] also has
score (1,1,2,2) — the complement of a (1,1,2,2) tournament is
another (1,1,2,2) tournament.

For a 5-cycle vertex set V with |V ∩ S| = k (k=2,3,4):
  The sub-tournament on V has some of its arcs reversed (those within V ∩ S).
  The NUMBER of directed Hamiltonian cycles on V can change.

But does the TOTAL directed 5-cycle count across ALL vertex sets cancel?

Approach: For each vertex set V with |V ∩ S| = k, compute the
change in directed cycle count. Sum over all such V.

The complement operation on the k vertices of V ∩ S reverses
C(k,2) arcs within V. The other arcs of V are unchanged.

For k=2: 1 arc reversed. The 5-cycle uses this arc iff both endpoints
are consecutive in the cycle. There are exactly 2 positions in a
5-cycle where these 2 vertices could be consecutive.

For k=3: 3 arcs reversed.
For k=4: 6 arcs reversed.

Let's verify dc5=0 and try to understand the CANCELLATION mechanism.
Track which vertex sets have multiplicity changes and the SIGN pattern.
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
print(f"WHY dc5 = 0: DETAILED CANCELLATION ANALYSIS AT n={n}")
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

        S = frozenset(subset)

        # Check all 5-vertex sets for directed cycle changes
        c5_changes_by_overlap = {}
        total_dc5 = 0

        for combo in combinations(range(n), 5):
            vs = frozenset(combo)
            m_A = count_directed_on_set(A, vs)
            m_B = count_directed_on_set(B, vs)
            delta = m_B - m_A

            if delta != 0:
                k = len(vs & S)
                if k not in c5_changes_by_overlap:
                    c5_changes_by_overlap[k] = []
                c5_changes_by_overlap[k].append({
                    'vs': vs,
                    'm_A': m_A,
                    'm_B': m_B,
                    'delta': delta,
                    'S_vertices': sorted(vs & S),
                    'ext_vertices': sorted(vs - S),
                })
                total_dc5 += delta

        if c5_changes_by_overlap:
            data.append({
                'subset': subset,
                'changes': c5_changes_by_overlap,
                'total_dc5': total_dc5,
            })
        break

    if trial % 500 == 0 and trial > 0:
        print(f"  ... {trial}/2000, found {len(data)}")

print(f"\nTotal examples with c5 multiplicity changes: {len(data)}")

# Is dc5 always 0?
dc5_zero = all(d['total_dc5'] == 0 for d in data)
print(f"dc5 = 0 always? {dc5_zero}")

# For each overlap level, what's the NET change?
print(f"\nPer-overlap-level NET changes:")
for d in data[:5]:
    print(f"\n  Subset S = {d['subset']}:")
    for k in sorted(d['changes'].keys()):
        changes = d['changes'][k]
        net = sum(c['delta'] for c in changes)
        positives = sum(c['delta'] for c in changes if c['delta'] > 0)
        negatives = sum(c['delta'] for c in changes if c['delta'] < 0)
        print(f"    |V cap S|={k}: {len(changes)} sets changed, net={net} (+{positives}/{negatives})")
        for c in changes[:3]:
            print(f"      {sorted(c['vs'])}: S-verts={c['S_vertices']}, {c['m_A']}->{c['m_B']} (delta={c['delta']})")

# Check: is the cancellation WITHIN each overlap level or ACROSS levels?
print(f"\n{'='*70}")
print(f"CANCELLATION STRUCTURE")
print(f"{'='*70}")

for d in data:
    by_k = {}
    for k, changes in d['changes'].items():
        by_k[k] = sum(c['delta'] for c in changes)
    nonzero = {k: v for k, v in by_k.items() if v != 0}
    if nonzero:
        print(f"  Cross-level cancellation needed: {by_k}")

# Check: which S-PAIRS appear in the changed sets?
print(f"\n{'='*70}")
print(f"S-PAIR ANALYSIS (for |V∩S|=2)")
print(f"{'='*70}")

spair_counter = Counter()
spair_delta = {}
for d in data:
    if 2 not in d['changes']:
        continue
    S = frozenset(d['subset'])
    for c in d['changes'][2]:
        spair = frozenset(c['S_vertices'])
        spair_counter[tuple(sorted(spair))] += 1
        if spair not in spair_delta:
            spair_delta[spair] = []
        spair_delta[spair].append(c['delta'])

print(f"  S-pairs with changes: {len(spair_counter)}")
for spair in sorted(spair_counter.keys()):
    s = frozenset(spair)
    deltas = spair_delta[s]
    print(f"    {spair}: {spair_counter[spair]} changes, net={sum(deltas)}, vals={sorted(deltas)}")

# The (1,1,2,2) tournament has special structure in its S-pairs
# Score 1 vertices: low-degree. Score 2 vertices: high-degree.
# Which pairs are affected most?
print(f"\n{'='*70}")
print(f"SCORE STRUCTURE OF S-PAIRS")
print(f"{'='*70}")

for d in data[:3]:
    S = list(d['subset'])
    A_local = bits_to_adj(np.random.randint(0, 1 << total_bits), n)  # dummy, need original
    # We'd need the original A here, but let's just show the S-pair delta structure
    if 2 in d['changes']:
        spairs = set()
        for c in d['changes'][2]:
            spairs.add(frozenset(c['S_vertices']))
        print(f"  S={d['subset']}, affected S-pairs at |cap S|=2: {[sorted(s) for s in spairs]}")

print("\nDone.")
