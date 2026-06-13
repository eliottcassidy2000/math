"""
spectral_hidden_dim.py -- kind-pasteur-2026-03-13-S61

The simplex profile ambiguity (6 profiles with 2 c7 values each) is
NOT explained by lambda-preserving reversals. What IS the hidden dimension?

Candidates to test:
1. Skew-adjacency spectrum (eigenvalues of A - A^T)
2. Adjacency spectrum (eigenvalues of A + A^T)
3. Characteristic polynomial of the lambda graph
4. The "interaction graph" between 3-cycles (Omega graph structure)
5. Some homological invariant
"""

import numpy as np
from itertools import combinations
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

def get_3cycles(A, n):
    """Get all directed 3-cycles as frozensets of vertices."""
    cycles = set()
    for u in range(n):
        for v in range(u+1, n):
            for w in range(v+1, n):
                if A[u][v] and A[v][w] and A[w][u]:
                    cycles.add(frozenset([u,v,w]))
                elif A[u][w] and A[w][v] and A[v][u]:
                    cycles.add(frozenset([u,v,w]))
    return cycles

def omega_graph(cycles):
    """Build the conflict graph Omega: edges between overlapping 3-cycles."""
    cycle_list = list(cycles)
    adj = np.zeros((len(cycle_list), len(cycle_list)), dtype=int)
    for i in range(len(cycle_list)):
        for j in range(i+1, len(cycle_list)):
            if cycle_list[i] & cycle_list[j]:  # shared vertex
                adj[i][j] = 1
                adj[j][i] = 1
    return adj, cycle_list

n = 7
total_bits = n * (n-1) // 2

print("=" * 60)
print("SPECTRAL HIDDEN DIMENSION ANALYSIS AT n=7")
print("=" * 60)

np.random.seed(42)

# Collect data
all_data = []
profile_groups = defaultdict(list)

for trial in range(10000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    pairs = get_simplex_data(A, n)
    c7 = count_directed_k_cycles(A, n, 7)
    profile = tuple(sorted(pairs))

    entry = {
        'bits': bits,
        'A': A.copy(),
        'profile': profile,
        'c7': c7,
    }
    all_data.append(entry)
    profile_groups[profile].append(entry)

# Find ambiguous profiles
ambig_profiles = {p for p, entries in profile_groups.items()
                  if len(set(e['c7'] for e in entries)) > 1}

print(f"\nTotal profiles: {len(profile_groups)}")
print(f"Ambiguous: {len(ambig_profiles)}")

for prof_idx, profile in enumerate(sorted(ambig_profiles,
                                           key=lambda p: -len(profile_groups[p]))):
    entries = profile_groups[profile]
    c7_vals = sorted(set(e['c7'] for e in entries))
    if len(c7_vals) != 2:
        continue

    c7_lo, c7_hi = c7_vals
    lo_entries = [e for e in entries if e['c7'] == c7_lo]
    hi_entries = [e for e in entries if e['c7'] == c7_hi]

    print(f"\n{'='*50}")
    print(f"Profile {prof_idx+1}: c7={c7_lo} vs {c7_hi} (dc7={c7_hi-c7_lo})")
    print(f"  n_lo={len(lo_entries)}, n_hi={len(hi_entries)}")

    # 1. Skew-adjacency spectrum
    lo_spectra = []
    hi_spectra = []
    for e in lo_entries[:10]:
        S = e['A'] - e['A'].T  # Skew-symmetric
        eigs = sorted(np.linalg.eigvalsh(1j * S).real)  # Eigenvalues of iS are real
        lo_spectra.append(tuple(round(x, 6) for x in eigs))
    for e in hi_entries[:10]:
        S = e['A'] - e['A'].T
        eigs = sorted(np.linalg.eigvalsh(1j * S).real)
        hi_spectra.append(tuple(round(x, 6) for x in eigs))

    lo_spec_set = set(lo_spectra)
    hi_spec_set = set(hi_spectra)
    overlap = lo_spec_set & hi_spec_set
    print(f"  Skew spectra: lo has {len(lo_spec_set)} distinct, hi has {len(hi_spec_set)} distinct")
    print(f"    Overlap: {len(overlap)}")
    if len(overlap) == 0:
        print(f"    SKEW SPECTRUM DISTINGUISHES lo from hi!")
    else:
        print(f"    Skew spectrum does NOT distinguish (shared spectra)")

    # 2. Omega graph invariants
    lo_omega = []
    hi_omega = []
    for e in lo_entries[:5]:
        cycles = get_3cycles(e['A'], n)
        adj_omega, _ = omega_graph(cycles)
        # Omega spectrum
        eigs = sorted(np.linalg.eigvalsh(adj_omega.astype(float)))
        eigs_r = tuple(round(x, 4) for x in eigs)
        lo_omega.append(eigs_r)
    for e in hi_entries[:5]:
        cycles = get_3cycles(e['A'], n)
        adj_omega, _ = omega_graph(cycles)
        eigs = sorted(np.linalg.eigvalsh(adj_omega.astype(float)))
        eigs_r = tuple(round(x, 4) for x in eigs)
        hi_omega.append(eigs_r)

    lo_omega_set = set(lo_omega)
    hi_omega_set = set(hi_omega)
    omega_overlap = lo_omega_set & hi_omega_set
    print(f"  Omega spectra: lo has {len(lo_omega_set)} distinct, hi has {len(hi_omega_set)} distinct")
    print(f"    Overlap: {len(omega_overlap)}")
    if len(omega_overlap) == 0:
        print(f"    OMEGA SPECTRUM DISTINGUISHES lo from hi!")
    else:
        print(f"    Omega spectrum does NOT distinguish")

    # 3. Adjacency matrix eigenvalues (of A + A^T)
    lo_adj_spec = []
    hi_adj_spec = []
    for e in lo_entries[:10]:
        sym = e['A'] + e['A'].T
        eigs = sorted(np.linalg.eigvalsh(sym.astype(float)))
        lo_adj_spec.append(tuple(round(x, 6) for x in eigs))
    for e in hi_entries[:10]:
        sym = e['A'] + e['A'].T
        eigs = sorted(np.linalg.eigvalsh(sym.astype(float)))
        hi_adj_spec.append(tuple(round(x, 6) for x in eigs))

    lo_adj_set = set(lo_adj_spec)
    hi_adj_set = set(hi_adj_spec)
    adj_overlap = lo_adj_set & hi_adj_set
    print(f"  Adj spectrum (A+A^T): lo has {len(lo_adj_set)} distinct, hi has {len(hi_adj_set)}")
    print(f"    Overlap: {len(adj_overlap)}")

    # 4. The sigma graph spectrum (treating sigma as a weighted graph)
    lo_sig_spec = []
    hi_sig_spec = []
    for e in lo_entries[:10]:
        sig_mat = np.zeros((n, n))
        pidx = 0
        for u in range(n):
            for v in range(u+1, n):
                s, l, d = e['profile'][pidx]  # Wrong — this is sorted
                pidx += 1
        # Use actual pairs
        pairs = get_simplex_data(e['A'], n)
        sig_mat = np.zeros((n, n))
        pidx = 0
        for u in range(n):
            for v in range(u+1, n):
                sig_mat[u][v] = pairs[pidx][0]  # sigma
                sig_mat[v][u] = pairs[pidx][0]
                pidx += 1
        eigs = sorted(np.linalg.eigvalsh(sig_mat))
        lo_sig_spec.append(tuple(round(x, 4) for x in eigs))

    for e in hi_entries[:10]:
        pairs = get_simplex_data(e['A'], n)
        sig_mat = np.zeros((n, n))
        pidx = 0
        for u in range(n):
            for v in range(u+1, n):
                sig_mat[u][v] = pairs[pidx][0]
                sig_mat[v][u] = pairs[pidx][0]
                pidx += 1
        eigs = sorted(np.linalg.eigvalsh(sig_mat))
        hi_sig_spec.append(tuple(round(x, 4) for x in eigs))

    lo_sig_set = set(lo_sig_spec)
    hi_sig_set = set(hi_sig_spec)
    sig_overlap = lo_sig_set & hi_sig_set
    print(f"  Sigma graph spectrum: lo has {len(lo_sig_set)} distinct, hi has {len(hi_sig_set)}")
    print(f"    Overlap: {len(sig_overlap)}")
    if len(sig_overlap) == 0:
        print(f"    SIGMA SPECTRUM DISTINGUISHES!")

    # 5. det(A - A^T) mod small primes
    lo_dets = []
    hi_dets = []
    for e in lo_entries[:20]:
        S = e['A'] - e['A'].T
        d = int(round(np.linalg.det(S.astype(float))))
        lo_dets.append(d)
    for e in hi_entries[:20]:
        S = e['A'] - e['A'].T
        d = int(round(np.linalg.det(S.astype(float))))
        hi_dets.append(d)
    print(f"  det(A-A^T): lo={set(lo_dets)}, hi={set(hi_dets)}")

    # 6. Pfaffian / det mod small primes
    lo_det2 = Counter(d % 7 for d in lo_dets)
    hi_det2 = Counter(d % 7 for d in hi_dets)
    print(f"  det mod 7: lo={dict(lo_det2)}, hi={dict(hi_det2)}")

    if prof_idx >= 5:
        break

# Summary
print(f"\n\n{'='*60}")
print("DISTINGUISHING INVARIANT SUMMARY")
print("=" * 60)
print("""
For each ambiguous simplex profile (same multiset of (sigma,lambda,delta)):
  - c3, c5, c4 are IDENTICAL
  - n_1122, n_trans4, n_1113 are IDENTICAL
  - n_reg5 is IDENTICAL

Testing spectral invariants to find the hidden dimension...
""")

print("\nDone.")
