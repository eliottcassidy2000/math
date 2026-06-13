"""
directed_multiplicity_vitali.py -- kind-pasteur-2026-03-13-S61
The Vitali atom changes DIRECTED CYCLE MULTIPLICITIES, not vertex set counts.

KEY INSIGHT from omega_ip_decomposition_n8.py:
At n=8, ALL H-changing Vitali reversals:
1. Preserve cycle vertex set counts (c3_vs, c5_vs, c7_vs all unchanged)
2. Preserve the overlap weight spectrum (W=0..5 pair counts identical)
3. But change directed cycle multiplicities on individual vertex sets

This means:
- A 5-vertex set {a,b,c,d,e} might carry 3 directed 5-cycles before reversal
  but only 2 after (or vice versa)
- The 7-vertex set {0..6} might carry 9 directed 7-cycles before but 8 after
- These multiplicity changes ARE the mechanism

The OCF counts DIRECTED cycles as Omega vertices.
Omega(T) vertex count = sum over vertex sets of multiplicity.
When multiplicity on a vertex set decreases, Omega loses vertices.
When multiplicity increases, Omega gains vertices.
ALL the gained/lost vertices on the SAME vertex set are pairwise adjacent
(they share ALL vertices), forming a clique in Omega.

This is a CLIQUE PERTURBATION: adding/removing vertices from existing cliques.

QUESTION: Does clique perturbation explain delta_H exactly?
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
    """Count directed Hamiltonian cycles on a vertex set."""
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
    """Find all vertex sets supporting at least one directed cycle."""
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
print("DIRECTED MULTIPLICITY ANALYSIS OF VITALI ATOM AT n=8")
print("=" * 70)

# Collect multiplicity change data
np.random.seed(2026)
data = []

for trial in range(3000):
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

        # Get vertex sets and their multiplicities
        mult_changes = {}
        total_dir_A = 0
        total_dir_B = 0

        for length in [3, 5, 7]:
            vsets_A = find_cycle_vertex_sets(A, n, length)
            vsets_B = find_cycle_vertex_sets(B, n, length)

            all_vsets = set(map(frozenset, vsets_A)) | set(map(frozenset, vsets_B))

            for vs in all_vsets:
                m_A = count_directed_on_set(A, vs)
                m_B = count_directed_on_set(B, vs)
                total_dir_A += m_A
                total_dir_B += m_B
                if m_A != m_B:
                    mult_changes[tuple(sorted(vs))] = (length, m_A, m_B)

        data.append({
            'delta_H': H_B - H_A,
            'total_dir_A': total_dir_A,
            'total_dir_B': total_dir_B,
            'delta_dir': total_dir_B - total_dir_A,
            'mult_changes': mult_changes,
            'bits': bits,
            'subset': subset,
        })
        break

    if trial % 500 == 0 and trial > 0:
        print(f"  ... {trial}/3000, found {len(data)} examples")

print(f"\nTotal H-changing examples: {len(data)}")

# Analysis 1: Is delta_H = 2 * delta_dir?
print(f"\n{'='*70}")
print(f"TEST: delta_H = 2 * delta_dir (directed cycle count)?")
print(f"{'='*70}")

match_2 = 0
for d in data:
    if d['delta_H'] == 2 * d['delta_dir']:
        match_2 += 1

print(f"  Match rate: {match_2}/{len(data)} ({100*match_2/len(data):.1f}%)")

if match_2 != len(data):
    print(f"\n  Non-matching examples:")
    for d in data:
        if d['delta_H'] != 2 * d['delta_dir']:
            print(f"    dH={d['delta_H']}, d_dir={d['delta_dir']}, 2*d_dir={2*d['delta_dir']}")
            if len(d['mult_changes']) <= 10:
                for vs, (length, m_A, m_B) in d['mult_changes'].items():
                    print(f"      c{length} on {list(vs)}: {m_A} -> {m_B}")

# Distribution
print(f"\n  delta_dir distribution:")
dir_dist = Counter(d['delta_dir'] for d in data)
for k in sorted(dir_dist.keys()):
    print(f"    delta_dir={k}: {dir_dist[k]}")

print(f"\n  delta_H distribution:")
dh_dist = Counter(d['delta_H'] for d in data)
for k in sorted(dh_dist.keys()):
    print(f"    delta_H={k}: {dh_dist[k]}")

# Analysis 2: Which cycle lengths change multiplicity?
print(f"\n{'='*70}")
print(f"WHICH CYCLE LENGTHS HAVE MULTIPLICITY CHANGES?")
print(f"{'='*70}")

len_change_count = Counter()
for d in data:
    lengths_changed = set()
    for vs, (length, m_A, m_B) in d['mult_changes'].items():
        lengths_changed.add(length)
    for l in lengths_changed:
        len_change_count[l] += 1

for l in sorted(len_change_count.keys()):
    print(f"  c{l} multiplicity changes in {len_change_count[l]}/{len(data)} examples ({100*len_change_count[l]/len(data):.1f}%)")

# Analysis 3: How many vertex sets have multiplicity changes per example?
print(f"\n{'='*70}")
print(f"NUMBER OF VERTEX SETS WITH MULTIPLICITY CHANGES")
print(f"{'='*70}")

nchanges_dist = Counter(len(d['mult_changes']) for d in data)
for k in sorted(nchanges_dist.keys()):
    print(f"  {k} vertex sets changed: {nchanges_dist[k]} examples")

# Analysis 4: Multiplicity change magnitudes
print(f"\n{'='*70}")
print(f"MULTIPLICITY CHANGE MAGNITUDES")
print(f"{'='*70}")

all_deltas = []
for d in data:
    for vs, (length, m_A, m_B) in d['mult_changes'].items():
        all_deltas.append((length, m_B - m_A))

delta_by_len = {}
for length, delta in all_deltas:
    if length not in delta_by_len:
        delta_by_len[length] = Counter()
    delta_by_len[length][delta] += 1

for l in sorted(delta_by_len.keys()):
    print(f"\n  c{l} multiplicity deltas:")
    for d in sorted(delta_by_len[l].keys()):
        print(f"    delta={d}: {delta_by_len[l][d]}")

# Analysis 5: Clique perturbation formula
print(f"\n{'='*70}")
print(f"CLIQUE PERTURBATION FORMULA")
print(f"{'='*70}")
print("""
When vertex set V changes multiplicity from m_A to m_B:
- The m_A directed cycles on V form a clique K_{m_A} in Omega_A
- The m_B directed cycles on V form a clique K_{m_B} in Omega_B
- This is a CLIQUE RESIZE operation

For the independence polynomial:
- Adding a vertex v to a clique doesn't change any independent set
  that already has a member of the clique (since v is adjacent to all of them)
- It CAN add to independent sets that have NO member of the clique
- So adding v to clique K increases alpha_j by the number of j-1 sized
  independent sets in Omega that avoid ALL vertices of K

Key: if V is a LARGE vertex set (e.g., all 8 vertices for c7),
then EVERY other cycle shares at least 1 vertex with V.
So every other Omega vertex is adjacent to all members of the K clique.
Adding/removing from K only changes alpha_0 (trivially) and alpha_1 (by 1 each).
This gives delta_H = 2 * (number added - number removed) = 2 * delta_mult.

But if V is a SMALL vertex set (e.g., 3 vertices for c3),
then cycles on disjoint vertex sets are NOT adjacent to the clique.
Adding a vertex to K can increase alpha_2 by the number of cycles
disjoint from V.
""")

# Test: For c7 changes (all 7 or 8 vertices), the clique is "universal"
# Every other cycle conflicts with c7, so delta_i_k = delta_mult for k=1 only
# For c3/c5 changes, disjoint cycles can create higher alpha changes

# Decompose: delta_H = 2*delta_c7_dir + ??? from c3/c5 changes
print(f"TESTING: delta_H = 2*delta_c7_dir + contribution from c3/c5 changes")
for d in data[:20]:
    delta_c7_dir = 0
    delta_c5_dir = 0
    delta_c3_dir = 0
    for vs, (length, m_A, m_B) in d['mult_changes'].items():
        if length == 7:
            delta_c7_dir += m_B - m_A
        elif length == 5:
            delta_c5_dir += m_B - m_A
        elif length == 3:
            delta_c3_dir += m_B - m_A

    # For c7: the 7-cycle uses 7 of 8 vertices, only 1 vertex free
    # A 3-cycle on 3 vertices: at most 1 vertex free from the 7-cycle
    # So a 3-cycle MUST share at least 2 vertices with the 7-set => always conflicts
    # Similarly for 5-cycles: must share at least 4 vertices
    # Therefore c7 clique is "universal" in Omega
    # delta from c7 change = 2 * delta_c7_dir

    c7_contribution = 2 * delta_c7_dir
    residual = d['delta_H'] - c7_contribution
    print(f"  dH={d['delta_H']:>3}, dc7_dir={delta_c7_dir:>2}, dc5_dir={delta_c5_dir:>2}, dc3_dir={delta_c3_dir:>2}, 2*dc7={c7_contribution:>3}, residual={residual:>3}")

# Further test: does 2*delta_dir_total work?
print(f"\n  Full test: delta_H vs 2*delta_dir_total")
for d in data[:20]:
    print(f"  dH={d['delta_H']:>3}, d_dir={d['delta_dir']:>2}, 2*d_dir={2*d['delta_dir']:>3}, match={'YES' if d['delta_H'] == 2*d['delta_dir'] else 'NO'}")

print("\nDone.")
