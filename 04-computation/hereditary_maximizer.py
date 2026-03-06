#!/usr/bin/env python3
"""
Hereditary Maximizer Chain: Does deleting a vertex from the n-maximizer
always give the (n-1)-maximizer?

Known: T_7 (H=189) -> delete v -> T_6 (H=45). [Paley deletion heredity]
Question: Is this true at every n? Does the chain go all the way down?

OEIS A038375: max H(T) = 1, 1, 3, 5, 15, 45, 189, 661, ...
                          n=  1  2  3  4   5   6    7    8

kind-pasteur-2026-03-06-S18f
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count, find_odd_cycles

# OEIS A038375 max H values
MAX_H = {1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189, 8: 661}

def all_vertex_deletions(T):
    """Return all (n-1)-subtournaments from deleting one vertex."""
    n = len(T)
    results = []
    for v in range(n):
        verts = [i for i in range(n) if i != v]
        sub = [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]
        results.append((v, sub))
    return results

def all_pair_deletions(T):
    """Return all (n-2)-subtournaments from deleting two vertices."""
    n = len(T)
    results = []
    for v1 in range(n):
        for v2 in range(v1+1, n):
            verts = [i for i in range(n) if i != v1 and i != v2]
            sub = [[T[verts[i]][verts[j]] for j in range(n-2)] for i in range(n-2)]
            results.append(((v1, v2), sub))
    return results

# ============================================================
# Phase 1: Exhaustive check for n = 3..7
# Find ALL maximizers at each n, then check vertex-deletion structure
# ============================================================
print("=" * 70)
print("HEREDITARY MAXIMIZER CHAIN")
print("=" * 70)

# Store all maximizers at each n
maximizers = {}

for n in range(3, 8):
    m = n * (n - 1) // 2
    total = 1 << m
    best_h = 0
    best_bits = []

    for bits in range(total):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h > best_h:
            best_h = h
            best_bits = [bits]
        elif h == best_h:
            best_bits.append(bits)

    maximizers[n] = (best_h, best_bits)
    assert best_h == MAX_H[n], f"n={n}: got {best_h}, expected {MAX_H[n]}"
    print(f"n={n}: max H = {best_h}, achieved by {len(best_bits)} tournaments")

# ============================================================
# Phase 2: Check hereditary property
# For each maximizer at n, check if ALL vertex-deletions give maximizer at n-1
# ============================================================
print(f"\n{'='*70}")
print("HEREDITARY CHECK: Does deleting ANY vertex from n-max give (n-1)-max?")
print("=" * 70)

for n in range(4, 8):
    best_h_n, best_bits_n = maximizers[n]
    best_h_n1, best_bits_n1 = maximizers[n-1]

    print(f"\nn={n}: max H={best_h_n} ({len(best_bits_n)} tournaments)")

    any_hereditary = False
    all_hereditary_count = 0

    for bits in best_bits_n:
        T = tournament_from_bits(n, bits)
        deletions = all_vertex_deletions(T)

        vertex_results = []
        all_give_max = True
        any_give_max = False

        for v, sub in deletions:
            h_sub = hamiltonian_path_count(sub)
            is_max = (h_sub == best_h_n1)
            vertex_results.append((v, h_sub, is_max))
            if is_max:
                any_give_max = True
            else:
                all_give_max = False

        if all_give_max:
            all_hereditary_count += 1

        if any_give_max:
            any_hereditary = True

        if len(best_bits_n) <= 20:  # print details for small counts
            status = "ALL MAX" if all_give_max else ("SOME MAX" if any_give_max else "NONE MAX")
            h_values = [h for _, h, _ in vertex_results]
            print(f"  bits={bits}: {status} | deletions H = {h_values}")

    print(f"  Summary: {all_hereditary_count}/{len(best_bits_n)} have ALL deletions giving max")
    print(f"  Any with at least one max deletion: {any_hereditary}")

# ============================================================
# Phase 3: Check n -> n-2 hereditary property (even-odd)
# ============================================================
print(f"\n{'='*70}")
print("n -> n-2 HEREDITARY CHECK (pair deletions)")
print("=" * 70)

for n in range(5, 8):
    best_h_n, best_bits_n = maximizers[n]
    best_h_n2, _ = maximizers[n-2]

    print(f"\nn={n} -> n-2={n-2}: max H({n})={best_h_n}, max H({n-2})={best_h_n2}")

    # Check first few maximizers
    for bits in best_bits_n[:5]:
        T = tournament_from_bits(n, bits)
        pair_dels = all_pair_deletions(T)

        max_count = 0
        total = 0
        for (v1, v2), sub in pair_dels:
            h_sub = hamiltonian_path_count(sub)
            if h_sub == best_h_n2:
                max_count += 1
            total += 1

        print(f"  bits={bits}: {max_count}/{total} pair-deletions give max H({n-2})")

# ============================================================
# Phase 4: The H-value spectrum of vertex deletions
# For EVERY n-maximizer, what is the multiset of deletion H-values?
# ============================================================
print(f"\n{'='*70}")
print("DELETION SPECTRUM: H-values of all vertex-deletions from maximizers")
print("=" * 70)

for n in range(4, 8):
    best_h_n, best_bits_n = maximizers[n]
    best_h_n1, _ = maximizers[n-1]
    diff = best_h_n - best_h_n1

    print(f"\nn={n}: max H={best_h_n}, max H({n-1})={best_h_n1}, diff={diff}")

    spectra = {}
    for bits in best_bits_n:
        T = tournament_from_bits(n, bits)
        deletions = all_vertex_deletions(T)
        h_vals = tuple(sorted([hamiltonian_path_count(sub) for _, sub in deletions], reverse=True))
        spectra[h_vals] = spectra.get(h_vals, 0) + 1

    print(f"  Distinct deletion spectra: {len(spectra)}")
    for spectrum, count in sorted(spectra.items(), key=lambda x: -x[1]):
        print(f"    {spectrum} x{count}")

# ============================================================
# Phase 5: Vertex-deletion H-sum identity
# Sum of H(T-v) over all v: is there a formula?
# ============================================================
print(f"\n{'='*70}")
print("VERTEX-DELETION H-SUM")
print("=" * 70)

for n in range(3, 8):
    m = n * (n - 1) // 2
    total = 1 << m if n <= 6 else min(1 << m, 1000)

    print(f"\nn={n}:")
    for bits in range(min(total, 10)):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        deletions = all_vertex_deletions(T)
        h_sum = sum(hamiltonian_path_count(sub) for _, sub in deletions)
        # Each edge of the original path appears in n-2 subtournaments
        # So sum of H(T-v) should be (n-2) * H(T)? No, that's wrong.
        # Actually, a Hamiltonian path in T-v extends to at most 2 paths in T
        # (by inserting v at one of n positions)
        # And each path in T, when we delete v_i, gives a path in T-v_i
        # that visits the other n-1 vertices. This path exists iff v_i is at
        # position 0 or n-1 (endpoints) OR the path restricted to T-v_i is
        # still Hamiltonian (it's always Hamiltonian: just remove v_i from
        # the sequence).
        # So: sum_{v} H(T-v) = sum over paths P sum over positions i of 1
        # = sum over paths P of n = n * H(T)? No...
        # Actually: deleting v from a Ham path P in T gives a sequence of
        # n-1 vertices. This IS a Ham path in T-v (all arcs exist).
        # But different paths P can give the same path in T-v.
        # If v is at position k in P, then P-v has a "gap" at position k
        # which is filled by direct arc from P[k-1] to P[k+1].
        # So each path P contributes exactly 1 to H(T-v) for each vertex v.
        # Wait: P has n vertices, so deleting any one gives n paths in
        # n different T-v's. So sum_v H(T-v) >= n * H(T).
        # But some of those T-v paths might be the same when derived from
        # different T-paths... No, each (P, v) pair gives a unique T-v path.
        # Actually: P -> P-v is a well-defined map from {Ham paths of T} to
        # {Ham paths of T-v}. Different P can map to the same T-v path.
        # P1 and P2 differ only in where v sits, and P1-v = P2-v is possible.
        # So sum_v H(T-v) >= H(T), but the exact relation depends on structure.
        ratio = h_sum / h if h > 0 else 0
        if bits < 10:
            print(f"  bits={bits}: H={h}, sum H(T-v)={h_sum}, ratio={ratio:.4f}")

    # Check ratio for maximizers
    if n >= 4:
        best_h, best_bits = maximizers[n]
        for bits in best_bits[:3]:
            T = tournament_from_bits(n, bits)
            h = hamiltonian_path_count(T)
            h_sum = sum(hamiltonian_path_count(sub) for _, sub in all_vertex_deletions(T))
            ratio = h_sum / h
            print(f"  MAX bits={bits}: H={h}, sum H(T-v)={h_sum}, ratio={ratio:.4f}")

print("\nDone.")
