"""
metagraph_recursive_gen.py — opus-2026-04-04-S12

RECURSIVE METAGRAPH GENERATION: Build G_n from G_{n-1} without
enumerating all 2^m tilings.

Algorithm:
1. Start with iso classes at n-1 (with Aut groups)
2. For each class C, enumerate extension orbits via Burnside
3. For each orbit, build a representative tournament at n
4. Group representatives by cheap invariants (H, score_seq)
5. Within groups, check isomorphism via canonical form
6. The MERGED set = V_n

This is O(V_{n-1} * 2^{n-1} / avg_Aut * n!) for step 5,
vs O(2^{C(n-1,2)} * n!) for brute force. MUCH faster when V_{n-1} is small.

For n=7: V_6 = 56 classes, 2^6 = 64 signatures → 56*64/avg_Aut ≈ 2000 pre-classes
vs 2^{15} = 32768 brute force. ~16x speedup.

For n=8: V_7 = 456 classes, 2^7 = 128 → ~35000 pre-classes
vs 2^{21} = 2M brute force. ~57x speedup.
"""
import numpy as np
from itertools import permutations
from collections import defaultdict, Counter
import time
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import count_hamiltonian_paths

def tournament_canonical(T, n):
    """Canonical form via lexicographic minimization over all permutations."""
    best = None
    for perm in permutations(range(n)):
        key = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(i+1, n))
        if best is None or key < best:
            best = key
    return best

def compute_aut_group(T, n):
    """Full automorphism group as list of permutation tuples."""
    auts = []
    for perm in permutations(range(n)):
        if all(T[perm[i]][perm[j]] == T[i][j] for i in range(n) for j in range(n)):
            auts.append(tuple(perm))
    return auts

def cycle_count_of_perm(perm, n):
    """Count cycles of a permutation."""
    visited = [False]*n
    cycles = 0
    for i in range(n):
        if not visited[i]:
            cycles += 1
            j = i
            while not visited[j]:
                visited[j] = True
                j = perm[j]
    return cycles

def score_sequence(T, n):
    """Sorted score sequence of tournament."""
    return tuple(sorted(int(sum(T[v])) for v in range(n)))

def cheap_invariants(T, n):
    """Fast invariants for grouping: (H, score_seq, c3)."""
    h = count_hamiltonian_paths(T, n)
    scores = score_sequence(T, n)
    c3 = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if T[a][b]+T[b][c]+T[c][a]==3 or T[a][c]+T[c][b]+T[b][a]==3:
                    c3 += 1
    return (int(h), scores, c3)

def recursive_generate(n_target):
    """Generate all iso classes at n_target from those at n_target-1."""
    n_sub = n_target - 1
    print(f"\n{'='*60}")
    print(f"RECURSIVE GENERATION: G_{n_target} from G_{n_sub}")
    print(f"{'='*60}")

    t_start = time.time()

    # Step 1: Get iso classes at n_sub
    print(f"\n  Step 1: Enumerate iso classes at n={n_sub}...")
    t0 = time.time()

    all_T_sub = []
    n_arcs = n_sub*(n_sub-1)//2
    for bits in range(1 << n_arcs):
        T = np.zeros((n_sub, n_sub), dtype=int)
        idx = 0
        for i in range(n_sub):
            for j in range(i+1, n_sub):
                if bits & (1 << idx):
                    T[i][j] = 1
                else:
                    T[j][i] = 1
                idx += 1
        all_T_sub.append(T)

    classes = defaultdict(list)
    for i, T in enumerate(all_T_sub):
        canon = tournament_canonical(T, n_sub)
        classes[canon].append(i)

    V_sub = len(classes)
    t1 = time.time()
    print(f"    V_{n_sub} = {V_sub} classes ({t1-t0:.1f}s)")

    # Step 2: Compute Aut groups
    print(f"\n  Step 2: Compute Aut groups...")
    t0 = time.time()

    class_data = {}
    for canon, indices in classes.items():
        T_rep = all_T_sub[indices[0]]
        auts = compute_aut_group(T_rep, n_sub)
        h_sub = count_hamiltonian_paths(T_rep, n_sub)
        class_data[canon] = {
            'T': T_rep, 'auts': auts, 'aut_size': len(auts), 'H': h_sub
        }

    t1 = time.time()
    print(f"    Done ({t1-t0:.1f}s)")

    # Step 3: Generate extension orbits
    print(f"\n  Step 3: Generate extension representatives...")
    t0 = time.time()

    pre_classes = []  # list of (T_extended, invariants, parent_canon)

    total_orbits = 0
    for canon, info in class_data.items():
        T_sub = info['T']
        auts = info['auts']

        # Burnside: #orbits = (1/|Aut|) sum 2^{cycles(alpha)}
        n_orbits = sum(2**cycle_count_of_perm(alpha, n_sub) for alpha in auts) // len(auts)
        total_orbits += n_orbits

        # Enumerate actual orbit representatives
        seen_sigmas = set()
        for sigma_int in range(1 << n_sub):
            sigma = tuple((sigma_int >> i) & 1 for i in range(n_sub))
            if sigma in seen_sigmas:
                continue

            # Compute orbit
            orbit = set()
            for alpha in auts:
                inv = [0]*n_sub
                for i in range(n_sub):
                    inv[alpha[i]] = i
                new_sigma = tuple(sigma[inv[j]] for j in range(n_sub))
                orbit.add(new_sigma)

            for s in orbit:
                seen_sigmas.add(s)

            # Build extended tournament
            T_ext = np.zeros((n_target, n_target), dtype=int)
            T_ext[:n_sub, :n_sub] = T_sub
            canon_sigma = min(orbit)
            for i in range(n_sub):
                if canon_sigma[i]:
                    T_ext[n_sub][i] = 1
                else:
                    T_ext[i][n_sub] = 1

            inv_data = cheap_invariants(T_ext, n_target)
            pre_classes.append((T_ext, inv_data, canon))

    t1 = time.time()
    print(f"    Total extension orbits: {total_orbits}")
    print(f"    Pre-classes generated: {len(pre_classes)} ({t1-t0:.1f}s)")

    # Step 4: Group by cheap invariants
    print(f"\n  Step 4: Group by (H, scores, c3)...")
    t0 = time.time()

    groups = defaultdict(list)
    for i, (T, inv, parent) in enumerate(pre_classes):
        groups[inv].append(i)

    n_groups = len(groups)
    max_group = max(len(v) for v in groups.values())
    t1 = time.time()
    print(f"    {n_groups} invariant groups, max group size {max_group} ({t1-t0:.1f}s)")

    # Step 5: Within groups, merge by canonical form
    print(f"\n  Step 5: Merge within groups (canonical form check)...")
    t0 = time.time()

    final_classes = {}  # canonical → representative tournament
    canon_checks = 0

    for inv, indices in groups.items():
        # Compute canonical form for each pre-class in this group
        local_canons = {}
        for i in indices:
            T = pre_classes[i][0]
            canon = tournament_canonical(T, n_target)
            canon_checks += 1
            if canon not in local_canons:
                local_canons[canon] = T
                if canon not in final_classes:
                    final_classes[canon] = {
                        'T': T, 'H': inv[0], 'scores': inv[1], 'c3': inv[2]
                    }

    V_n = len(final_classes)
    t1 = time.time()
    print(f"    Canonical form checks: {canon_checks}")
    print(f"    V_{n_target} = {V_n} ({t1-t0:.1f}s)")

    t_total = time.time() - t_start
    print(f"\n  TOTAL TIME: {t_total:.1f}s")
    print(f"  Merging ratio: {total_orbits}/{V_n} = {total_orbits/V_n:.2f}")
    print(f"  Speedup vs brute force ({1<<(n_target*(n_target-1)//2)} tilings): "
          f"~{(1<<(n_target*(n_target-1)//2))/len(pre_classes):.1f}x fewer canonical checks")

    # H distribution
    h_dist = Counter(d['H'] for d in final_classes.values())
    print(f"\n  H distribution at n={n_target}: {dict(sorted(h_dist.items()))}")

    return final_classes

if __name__ == "__main__":
    # Test at small n
    for n in [5, 6]:
        recursive_generate(n)

    # Try n=7
    print(f"\n\n{'#'*60}")
    print(f"ATTEMPTING n=7 RECURSIVE GENERATION")
    print(f"{'#'*60}")
    result = recursive_generate(7)
