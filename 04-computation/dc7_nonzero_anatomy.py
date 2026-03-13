"""
dc7_nonzero_anatomy.py -- kind-pasteur-2026-03-13-S61

Specifically find dc7 != 0 cases at n=7 and analyze
the unpaired Hamiltonian cycles.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict

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

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

def count_directed_k_cycles(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

n = 7
total_bits = n * (n-1) // 2

print("=" * 60)
print("ANATOMY OF dc7 != 0 CASES AT n=7")
print("=" * 60)

# Find dc7 != 0 examples
np.random.seed(123)  # Different seed
nz_examples = []

for trial in range(50000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(L, lambda_graph(B, n)):
            continue

        c7_A = count_directed_k_cycles(A, n, 7)
        c7_B = count_directed_k_cycles(B, n, 7)
        dc7 = c7_B - c7_A

        if dc7 != 0:
            nz_examples.append({
                'bits': bits,
                'subset': subset,
                'dc7': dc7,
                'c7_A': c7_A,
                'c7_B': c7_B,
                'A': A.copy(),
                'B': B.copy(),
            })
        break

    if len(nz_examples) >= 20:
        break

print(f"Found {len(nz_examples)} dc7 != 0 examples")

for ex_idx in range(min(5, len(nz_examples))):
    ex = nz_examples[ex_idx]
    A = ex['A']
    B = ex['B']
    S = set(ex['subset'])
    S_sorted = sorted(S)
    outside = sorted(set(range(n)) - S)

    print(f"\n{'='*40}")
    print(f"Example {ex_idx}: dc7={ex['dc7']}, c7_A={ex['c7_A']}, c7_B={ex['c7_B']}")
    print(f"S={ex['subset']}, outside={outside}")

    # Enumerate ALL Hamiltonian cycles in A and B
    # Use canonical form: start from smallest vertex, 7 rotations per undirected cycle
    all_cycles_A = []
    all_cycles_B = []

    for perm in permutations(range(1, n)):
        cycle = (0,) + perm
        valid_A = all(A[cycle[i]][cycle[(i+1)%n]] for i in range(n))
        valid_B = all(B[cycle[i]][cycle[(i+1)%n]] for i in range(n))

        if valid_A:
            all_cycles_A.append(cycle)
        if valid_B:
            all_cycles_B.append(cycle)

    # Deduplicate by rotation (each directed cycle has 7 rotations)
    def canon(cycle, n):
        """Canonical form: rotate to start with 0."""
        idx = cycle.index(0)
        return tuple(cycle[(idx + j) % n] for j in range(n))

    unique_A = set()
    unique_B = set()
    for c in all_cycles_A:
        unique_A.add(canon(c, n))
    for c in all_cycles_B:
        unique_B.add(canon(c, n))

    print(f"  Unique cycles: A={len(unique_A)}, B={len(unique_B)}")
    print(f"  Shared: {len(unique_A & unique_B)}")

    # Show each cycle and its S-structure
    for label, cycle_set in [("A-only", unique_A - unique_B),
                              ("B-only", unique_B - unique_A),
                              ("Shared", unique_A & unique_B)]:
        if len(cycle_set) == 0:
            continue
        print(f"\n  {label} cycles ({len(cycle_set)}):")
        for cycle in sorted(cycle_set)[:10]:
            # S-structure
            in_S_markers = ''.join('S' if v in S else 'X' for v in cycle)
            # S-S arcs
            s_arcs = sum(1 for i in range(n) if cycle[i] in S and cycle[(i+1)%n] in S)
            # The specific S-arcs used
            s_arc_list = [(cycle[i], cycle[(i+1)%n]) for i in range(n)
                         if cycle[i] in S and cycle[(i+1)%n] in S]
            # The specific X-X arcs used
            x_arc_list = [(cycle[i], cycle[(i+1)%n]) for i in range(n)
                         if cycle[i] not in S and cycle[(i+1)%n] not in S]

            print(f"    {cycle} [{in_S_markers}] S-arcs={s_arc_list} X-arcs={x_arc_list}")

    # Try the S-reversal conjugation
    print(f"\n  Conjugation analysis:")
    conj_map = {}  # maps A-cycle -> conjugated version
    for cycle in sorted(unique_A):
        # Find S-segments and reverse each
        conj = list(cycle)
        in_S = [v in S for v in cycle]

        i = 0
        while i < n:
            if in_S[i]:
                j = i
                while j < n and in_S[j]:
                    j += 1
                # Reverse S-segment [i, j)
                seg = list(cycle[i:j])
                seg.reverse()
                for k, v in enumerate(seg):
                    conj[i + k] = v
                i = j
            else:
                i += 1

        conj = tuple(conj)
        # Canonicalize
        conj_canon = canon(conj, n)

        # Is the conjugated cycle valid in B?
        valid_in_B = all(B[conj[i]][conj[(i+1)%n]] for i in range(n))

        conj_map[cycle] = (conj_canon, valid_in_B)
        if not valid_in_B:
            # Which arcs fail?
            failed_arcs = [(conj[i], conj[(i+1)%n]) for i in range(n) if not B[conj[i]][conj[(i+1)%n]]]
            print(f"    A-cycle {cycle} -> conj {conj} : FAILS at arcs {failed_arcs}")
        else:
            print(f"    A-cycle {cycle} -> conj {conj_canon} : SUCCESS (in B)")

    # Also conjugate B-cycles back to A
    print(f"\n  Reverse conjugation (B -> A):")
    for cycle in sorted(unique_B):
        conj = list(cycle)
        in_S = [v in S for v in cycle]

        i = 0
        while i < n:
            if in_S[i]:
                j = i
                while j < n and in_S[j]:
                    j += 1
                seg = list(cycle[i:j])
                seg.reverse()
                for k, v in enumerate(seg):
                    conj[i + k] = v
                i = j
            else:
                i += 1

        conj = tuple(conj)
        conj_canon = canon(conj, n)
        valid_in_A = all(A[conj[i]][conj[(i+1)%n]] for i in range(n))

        if not valid_in_A:
            failed_arcs = [(conj[i], conj[(i+1)%n]) for i in range(n) if not A[conj[i]][conj[(i+1)%n]]]
            print(f"    B-cycle {cycle} -> conj {conj} : FAILS at arcs {failed_arcs}")
        else:
            print(f"    B-cycle {cycle} -> conj {conj_canon} : SUCCESS (in A)")

# Summary statistics
print(f"\n{'='*60}")
print("SUMMARY")
print(f"{'='*60}")

dc7_dist = Counter(ex['dc7'] for ex in nz_examples)
print(f"dc7 distribution among nonzero: {dict(sorted(dc7_dist.items()))}")

# c7 values
print(f"\nc7_A vs c7_B for dc7 != 0:")
for ex in nz_examples[:10]:
    print(f"  dc7={ex['dc7']:+d}: c7_A={ex['c7_A']}, c7_B={ex['c7_B']}")

# XX type for nonzero
for ex in nz_examples:
    outside = sorted(set(range(n)) - set(ex['subset']))
    x0, x1, x2 = outside
    A = ex['A']
    c3_xx = (A[x0][x1]*A[x1][x2]*A[x2][x0] + A[x0][x2]*A[x2][x1]*A[x1][x0])
    ex['xx_type'] = 'C3' if c3_xx > 0 else 'T3'

xx_dist = Counter(ex['xx_type'] for ex in nz_examples)
print(f"\nXX type among dc7 != 0: {dict(xx_dist)}")

print("\nDone.")
