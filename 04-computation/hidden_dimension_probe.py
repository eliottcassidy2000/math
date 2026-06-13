"""
hidden_dimension_probe.py -- kind-pasteur-2026-03-13-S61

The simplex profile (multiset of (sigma, lambda, delta) over all 21 pairs)
almost determines c7 — only 6 profiles are ambiguous, each by exactly ±1.

This script probes: WHAT is the hidden dimension that distinguishes them?

Key candidates:
1. The LABELED assignment (which pair gets which simplex point)
2. The orientation pattern of arcs within each simplex class
3. Presence/absence of specific Vitali atom sub-tournaments
4. The "phase" of the 3-cycle structure (handedness, chirality)
5. The parity of some combinatorial quantity (a Z/2 invariant)

Since the difference is always exactly 1 in c7, this hidden dimension
is likely a BINARY (Z/2) variable — a parity/chirality bit.
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

def count_directed_k_cycles(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

def get_simplex_data(A, n):
    L = lambda_graph(A, n)
    A2 = A @ A
    pairs = []
    for u in range(n):
        for v in range(u+1, n):
            sig = n - 2 - int(A2[u][v]) - int(A2[v][u])
            lam = int(L[u][v])
            delta = n - 2 - sig - lam
            pairs.append((sig, lam, delta))
    return pairs

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

def has_vitali_atom(A, n):
    """Check if tournament has a (1,1,2,2) sub-tournament that preserves lambda."""
    L = lambda_graph(A, n)
    atoms = []
    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        # Check lambda preservation
        B = A.copy()
        for i in subset:
            for j in subset:
                if i != j:
                    B[i][j] = A[j][i]
        L2 = lambda_graph(B, n)
        if np.array_equal(L, L2):
            atoms.append(subset)
    return atoms

n = 7
total_bits = n * (n-1) // 2

print("=" * 60)
print("HIDDEN DIMENSION PROBE AT n=7")
print("=" * 60)

np.random.seed(42)

# Collect data, focusing on ambiguous profiles
all_data = []
profile_groups = defaultdict(list)

for trial in range(8000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    pairs = get_simplex_data(A, n)
    c7 = count_directed_k_cycles(A, n, 7)
    c3 = count_directed_k_cycles(A, n, 3)
    c5 = count_directed_k_cycles(A, n, 5)
    scores = tuple(sorted([int(sum(A[i])) for i in range(n)]))
    profile = tuple(sorted(pairs))

    entry = {
        'bits': bits,
        'A': A.copy(),
        'profile': profile,
        'pairs': pairs,
        'c7': c7,
        'c3': c3,
        'c5': c5,
        'scores': scores,
    }
    all_data.append(entry)
    profile_groups[profile].append(entry)

# Find ambiguous profiles
ambig_profiles = {p: [e['c7'] for e in entries]
                  for p, entries in profile_groups.items()
                  if len(set(e['c7'] for e in entries)) > 1}

print(f"\nTotal profiles: {len(profile_groups)}")
print(f"Ambiguous profiles: {len(ambig_profiles)}")

# For each ambiguous profile, analyze what distinguishes the two c7 values
for prof_idx, (profile, _) in enumerate(sorted(ambig_profiles.items(),
                                                key=lambda x: -len(profile_groups[x[0]]))):
    entries = profile_groups[profile]
    c7_vals = sorted(set(e['c7'] for e in entries))

    if len(c7_vals) != 2:
        continue

    c7_lo, c7_hi = c7_vals
    lo_entries = [e for e in entries if e['c7'] == c7_lo]
    hi_entries = [e for e in entries if e['c7'] == c7_hi]

    print(f"\n{'='*50}")
    print(f"Profile {prof_idx+1}: c7 in {{{c7_lo}, {c7_hi}}}")
    print(f"  n_lo={len(lo_entries)}, n_hi={len(hi_entries)}")

    # Probe 1: Vitali atoms
    lo_has_atom = 0
    hi_has_atom = 0
    for e in lo_entries[:10]:
        atoms = has_vitali_atom(e['A'], n)
        if atoms:
            lo_has_atom += 1
    for e in hi_entries[:10]:
        atoms = has_vitali_atom(e['A'], n)
        if atoms:
            hi_has_atom += 1
    print(f"  Vitali atoms: lo={lo_has_atom}/{min(10,len(lo_entries))}, "
          f"hi={hi_has_atom}/{min(10,len(hi_entries))}")

    # Probe 2: The labeled pair assignment
    # For each tournament, compute the "labeled profile" = the actual (u,v) -> (s,l,d) map
    # Two tournaments with same multiset but different labeling may differ
    # But since we're looking at isomorphism-free properties, let's look at
    # the GRAPH structure induced by the simplex assignment

    # Probe 3: Parity of specific combinatorial quantities
    lo_parities = []
    hi_parities = []
    for e in lo_entries:
        A = e['A']
        # Compute number of transitive triples (delta_total)
        dt = sum(d for s, l, d in e['pairs'])
        # Compute parity of various things
        c3 = e['c3']
        # Parity of c7
        p_c7 = e['c7'] % 2
        # Compute det(A) mod 2
        # For adjacency matrix of tournament...
        # Number of 4-cycles
        A4 = np.linalg.matrix_power(A, 4)
        c4_trace = int(np.trace(A4))
        # Also: parity of number of (1,1,2,2) subsets
        n_1122 = 0
        for subset in combinations(range(n), 4):
            if sub_scores(A, n, list(subset)) == (1, 1, 2, 2):
                n_1122 += 1
        lo_parities.append((p_c7, c4_trace % 2, n_1122 % 2, dt % 2))

    for e in hi_entries:
        A = e['A']
        dt = sum(d for s, l, d in e['pairs'])
        c3 = e['c3']
        p_c7 = e['c7'] % 2
        A4 = np.linalg.matrix_power(A, 4)
        c4_trace = int(np.trace(A4))
        n_1122 = 0
        for subset in combinations(range(n), 4):
            if sub_scores(A, n, list(subset)) == (1, 1, 2, 2):
                n_1122 += 1
        hi_parities.append((p_c7, c4_trace % 2, n_1122 % 2, dt % 2))

    lo_par_dist = Counter(lo_parities)
    hi_par_dist = Counter(hi_parities)
    print(f"  Parities (c7%2, tr(A^4)%2, n_1122%2, delta_total%2):")
    print(f"    lo: {dict(lo_par_dist)}")
    print(f"    hi: {dict(hi_par_dist)}")

    # Probe 4: Check if the difference is EXACTLY a Vitali atom flip
    # Take one lo and one hi example, check if they're related by a 4-subset reversal
    if lo_entries and hi_entries:
        A_lo = lo_entries[0]['A']
        A_hi = hi_entries[0]['A']
        found_link = False
        for subset in combinations(range(n), 4):
            B = A_lo.copy()
            for i in subset:
                for j in subset:
                    if i != j:
                        B[i][j] = A_lo[j][i]
            if np.array_equal(B, A_hi):
                print(f"  DIRECT LINK: lo -> hi by reversing {subset}")
                found_link = True
                break
        if not found_link:
            # Check if they're related by relabeling + atom flip
            print(f"  No direct 4-subset reversal link (expected — different labelings)")

    # Probe 5: For the labeled profile, check the "pair graph" structure
    # Build a graph where vertices are the 21 pairs, and two pairs are connected
    # if they share a vertex AND have the same (sigma, lambda, delta)
    # This graph structure might differ between lo and hi

    # Actually, let's try something more concrete:
    # For each tournament, compute the number of 5-element subsets
    # that contain a (1,1,2,2) sub-tournament
    lo_contain = []
    hi_contain = []
    for e in lo_entries[:5]:
        A = e['A']
        cnt = 0
        for sub5 in combinations(range(n), 5):
            for sub4 in combinations(sub5, 4):
                if sub_scores(A, n, list(sub4)) == (1, 1, 2, 2):
                    cnt += 1
                    break
        lo_contain.append(cnt)
    for e in hi_entries[:5]:
        A = e['A']
        cnt = 0
        for sub5 in combinations(range(n), 5):
            for sub4 in combinations(sub5, 4):
                if sub_scores(A, n, list(sub4)) == (1, 1, 2, 2):
                    cnt += 1
                    break
        hi_contain.append(cnt)
    print(f"  5-subsets containing (1,1,2,2): lo={lo_contain}, hi={hi_contain}")

    if prof_idx >= 3:  # Only show first 4 ambiguous profiles
        break

# Probe 6: Global invariant search
# For each ambiguous profile, check a wide range of graph-theoretic invariants
print(f"\n\n{'='*60}")
print("GLOBAL INVARIANT SEARCH")
print("='*60")

for prof_idx, (profile, _) in enumerate(sorted(ambig_profiles.items(),
                                                key=lambda x: -len(profile_groups[x[0]]))):
    entries = profile_groups[profile]
    c7_vals = sorted(set(e['c7'] for e in entries))
    if len(c7_vals) != 2:
        continue
    c7_lo, c7_hi = c7_vals

    print(f"\nProfile {prof_idx+1}: c7={c7_lo} vs {c7_hi}")

    # For each tournament, compute a bunch of invariants
    for label, group in [("lo", [e for e in entries if e['c7'] == c7_lo][:3]),
                         ("hi", [e for e in entries if e['c7'] == c7_hi][:3])]:
        for e in group:
            A = e['A']
            scores = [int(sum(A[i])) for i in range(n)]

            # Number of (1,1,2,2) subsets
            n_1122 = sum(1 for s in combinations(range(n), 4)
                        if sub_scores(A, n, list(s)) == (1, 1, 2, 2))

            # Number of (0,1,2,3) subsets (transitive 4-tournaments)
            n_trans4 = sum(1 for s in combinations(range(n), 4)
                         if sub_scores(A, n, list(s)) == (0, 1, 2, 3))

            # Number of (1,1,1,3) subsets (near-regular 4-tournaments)
            n_1113 = sum(1 for s in combinations(range(n), 4)
                        if sub_scores(A, n, list(s)) == (1, 1, 1, 3))

            # 4-cycles (directed)
            c4 = count_directed_k_cycles(A, n, 4)

            # Number of 5-vertex sub-tournaments with score (2,2,2,2,2)
            n_reg5 = sum(1 for s in combinations(range(n), 5)
                        if sub_scores(A, n, list(s)) == (2, 2, 2, 2, 2))

            print(f"    {label}: c7={e['c7']}, c3={e['c3']}, c5={e['c5']}, "
                  f"n_1122={n_1122}, n_trans4={n_trans4}, n_1113={n_1113}, "
                  f"c4={c4}, n_reg5={n_reg5}")

    if prof_idx >= 3:
        break

# Probe 7: The Vitali atom connection
# For each ambiguous profile, find a SPECIFIC pair of tournaments
# with different c7 that are related by a Vitali atom
print(f"\n\n{'='*60}")
print("VITALI ATOM CONNECTION TO AMBIGUITY")
print("=" * 60)

for prof_idx, (profile, _) in enumerate(sorted(ambig_profiles.items(),
                                                key=lambda x: -len(profile_groups[x[0]]))):
    entries = profile_groups[profile]
    c7_vals = sorted(set(e['c7'] for e in entries))
    if len(c7_vals) != 2:
        continue
    c7_lo, c7_hi = c7_vals

    lo_entries = [e for e in entries if e['c7'] == c7_lo]
    hi_entries = [e for e in entries if e['c7'] == c7_hi]

    print(f"\nProfile {prof_idx+1}: c7={c7_lo} vs {c7_hi}")

    # For each lo tournament, apply all Vitali atoms and check if any produces a hi tournament
    # with the SAME simplex profile
    found = False
    for e_lo in lo_entries[:5]:
        A = e_lo['A']
        atoms = has_vitali_atom(A, n)
        for atom_subset in atoms:
            B = A.copy()
            for i in atom_subset:
                for j in atom_subset:
                    if i != j:
                        B[i][j] = A[j][i]
            c7_B = count_directed_k_cycles(B, n, 7)
            pairs_B = get_simplex_data(B, n)
            profile_B = tuple(sorted(pairs_B))

            if profile_B == profile and c7_B != e_lo['c7']:
                print(f"  FOUND: atom {atom_subset} maps c7={e_lo['c7']} -> c7={c7_B}")
                print(f"    Same simplex profile: YES")
                found = True
                break
            elif c7_B != e_lo['c7']:
                print(f"  Atom {atom_subset}: c7 {e_lo['c7']}->{c7_B}, "
                      f"profile same: {profile_B == profile}")
        if found:
            break

    if not found:
        print(f"  No Vitali atom connection found in first 5 lo examples")

    if prof_idx >= 5:
        break

print("\n\nDone.")
