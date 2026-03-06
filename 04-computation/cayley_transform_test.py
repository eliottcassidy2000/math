#!/usr/bin/env python3
"""
Irving-Omar Cayley Transform: test whether the matrix (I+A)(I-A)^{-1}
or related constructions give information about I(Omega(T), x).

From Irving-Omar (arXiv:2412.10572): The odd-cycle extraction from the
cycle generating function uses arctanh splitting:
    sum_{k odd} tr(A^k)/k = (1/2) tr(log((I+zA)/(I-zA)))

The matrix C(z) = (I+zA)(I-zA)^{-1} is the Cayley transform of zA.
For tournaments, the adjacency matrix A has entries +1 (arc i->j) and 0 (no arc).
The skew-symmetric matrix S = A - A^T has entries +1/-1.

Key question: Does any natural matrix derived from the Cayley transform
of S (or A) have eigenvalues related to I(Omega(T), x)?

Specifically, if we write the OCF as:
    H(T) = I(Omega(T), 2) = sum_{S in Ind(Omega)} 2^|S|

Then the independence polynomial coefficients alpha_k count independent
sets of size k. Can we connect these to trace formulas from the Cayley transform?

kind-pasteur-2026-03-06-S18f
"""
import sys
import os
import numpy as np
from itertools import combinations

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (tournament_from_bits, hamiltonian_path_count,
                             find_odd_cycles, conflict_graph)

def indep_poly_coeffs(adj):
    m = len(adj)
    if m == 0:
        return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            coeffs[bin(mask).count('1')] += 1
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def skew_adj(T):
    """Skew-symmetric adjacency matrix: S[i][j] = T[i][j] - T[j][i]."""
    n = len(T)
    S = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                S[i][j] = T[i][j] - T[j][i]
    return S

# ============================================================
# Test 1: Cayley transform eigenvalues
# C(z) = (I + z*S)(I - z*S)^{-1}
# For S skew-symmetric, eigenvalues of S are +/- i*sigma
# Eigenvalues of C(z) at z=1: (1 + i*sigma)/(1 - i*sigma) = e^{2i*arctan(sigma)}
# These lie on the unit circle.
# ============================================================
print("=" * 70)
print("CAYLEY TRANSFORM EIGENVALUES")
print("=" * 70)

for n in [4, 5, 6]:
    print(f"\nn={n}:")
    m = n * (n-1) // 2
    count = min(1 << m, 50)

    for bits in range(count):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        S = skew_adj(T)

        # Eigenvalues of S (purely imaginary)
        eig_S = np.linalg.eigvals(S)
        sigmas = sorted([abs(e.imag) for e in eig_S], reverse=True)

        # Cayley transform at z=1
        I_n = np.eye(n)
        try:
            C = (I_n + S) @ np.linalg.inv(I_n - S)
        except np.linalg.LinAlgError:
            continue

        eig_C = np.linalg.eigvals(C)

        # det(I - S) -- this is related to permanent/determinant
        det_I_minus_S = np.linalg.det(I_n - S)
        det_I_plus_S = np.linalg.det(I_n + S)

        # For skew-symmetric S:
        # det(I+S) = det(I-S) when n is even (since det(-S)=(-1)^n det(S)=det(S) for even n)
        # Pfaffian structure

        cycles = find_odd_cycles(T)
        num_cycles = len(cycles) if cycles else 0

        if bits < 10:
            print(f"  bits={bits}: H={h}, |cycles|={num_cycles}, "
                  f"sigmas={[f'{s:.2f}' for s in sigmas[:3]]}, "
                  f"det(I-S)={det_I_minus_S:.0f}, det(I+S)={det_I_plus_S:.0f}")

# ============================================================
# Test 2: Trace formulas and cycle counting
# tr(S^k) for odd k counts directed k-cycles with signs
# sum_{k odd} tr(S^k)/k = log(det(I+S)/det(I-S))/2 = arctanh trace
# ============================================================
print(f"\n{'='*70}")
print("TRACE FORMULA: sum tr(S^k)/k for odd k vs cycle counts")
print("=" * 70)

for n in [4, 5, 6]:
    print(f"\nn={n}:")
    m = n * (n-1) // 2
    count = min(1 << m, 30)

    for bits in range(count):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        S = skew_adj(T)

        # Compute tr(S^k) for k = 1, 3, 5, 7
        traces = {}
        Sk = S.copy()
        for k in range(1, min(n+1, 8)):
            tr = np.trace(Sk)
            if k % 2 == 1:
                traces[k] = tr
            Sk = Sk @ S

        # Arctanh trace sum
        arctanh_sum = sum(traces[k] / k for k in traces)

        # Actual cycle counts
        cycles = find_odd_cycles(T)
        if cycles:
            cg = conflict_graph(cycles)
            ip = indep_poly_coeffs(cg)
            c3 = sum(1 for c in cycles if len(c) == 3)
            c5 = sum(1 for c in cycles if len(c) == 5)
            c7 = sum(1 for c in cycles if len(c) == 7)
        else:
            ip = [1]
            c3 = c5 = c7 = 0

        if bits < 10:
            print(f"  bits={bits}: H={h}, tr(S^k)/k={{3:{traces.get(3,0)/3:.0f}, "
                  f"5:{traces.get(5,0)/5:.0f}}}, "
                  f"c3={c3}, c5={c5}, arctanh_sum={arctanh_sum:.1f}")

        # KEY RELATION: tr(S^3)/3 should count directed 3-cycles
        # In a tournament, tr(A^3)/3 = number of directed 3-cycles
        # (but here S = A - A^T, so S^3 has different trace)
        # tr(S^3) = tr((A-A^T)^3)
        # = tr(A^3) - 3*tr(A^2*A^T) + 3*tr(A*A^T^2) - tr(A^T^3)
        # For tournaments, A^T = J - I - A, so this gets complicated

# ============================================================
# Test 3: The KEY relation — tr(A^k) for ADJACENCY matrix A
# For tournament adjacency matrix A[i][j] = 1 if i->j, 0 otherwise:
# tr(A^3) = 3 * (number of directed 3-cycles) since each 3-cycle
# contributes 1 to each of its 3 cyclic rotations.
# Similarly tr(A^5) = 5 * (number of directed 5-cycles)... but overcounting!
# tr(A^5)/5 counts 5-cycles BUT includes "backtracking" walks.
# ============================================================
print(f"\n{'='*70}")
print("ADJACENCY MATRIX: tr(A^k) vs cycle counts")
print("=" * 70)

for n in [4, 5, 6, 7]:
    print(f"\nn={n}:")
    m = n * (n-1) // 2
    count = min(1 << m, 30) if n <= 6 else 10

    for bits in range(count):
        T = tournament_from_bits(n, bits)
        A = np.array(T, dtype=float)
        h = hamiltonian_path_count(T)

        # tr(A^k) for k = 3, 5, 7
        traces = {}
        Ak = A.copy()
        for k in range(1, min(n+1, 8)):
            Ak = Ak @ A if k > 1 else A
            if k >= 3 and k % 2 == 1:
                traces[k] = int(round(np.trace(Ak)))

        # Actual cycle counts
        cycles = find_odd_cycles(T)
        c3 = sum(1 for c in cycles if len(c) == 3) if cycles else 0
        c5 = sum(1 for c in cycles if len(c) == 5) if cycles else 0

        if bits < 8:
            tr3_div = traces.get(3, 0) // 3 if 3 in traces else 0
            tr5_div = traces.get(5, 0) // 5 if 5 in traces else 0
            # For directed cycles: each k-cycle is counted k times in tr(A^k)
            # But tr(A^5)/5 includes BOTH 5-cycles and products of shorter cycles
            # Actually no: tr(A^k) counts CLOSED WALKS of length k, which
            # for a tournament (no 2-cycles) decomposes into:
            # - genuine k-cycles (each counted k times)
            # - composite walks (e.g., 3+2 for k=5, but 2-cycles don't exist)
            # Wait, in a tournament there ARE no 2-cycles (no i->j and j->i).
            # So tr(A^3)/3 = number of directed 3-cycles (exact!)
            # tr(A^5)/5 = number of directed 5-cycles + (terms involving 3-cycles)?
            # Actually: a closed walk of length 5 on a tournament is either:
            # - A 5-cycle: i->j->k->l->m->i (5 rotations)
            # - A walk that revisits vertices (but in a tournament, from any vertex
            #   there's exactly one arc direction to each other vertex, so
            #   a walk i->j->k->i->j->k is impossible since it would need
            #   both i->j and j->k->i->j which is fine... actually in tournaments
            #   closed walks of length 5 can revisit vertices)
            # So tr(A^5)/5 is NOT just 5-cycles.
            print(f"  bits={bits}: H={h}, tr(A^3)/3={tr3_div}=c3?{tr3_div==c3}, "
                  f"tr(A^5)/5={tr5_div}, c5={c5}")

# ============================================================
# Test 4: The Pfaffian connection
# For even n, det(S) = Pf(S)^2 where S is skew-symmetric.
# The Pfaffian counts perfect matchings of the "tournament graph."
# Is Pf(S) related to H(T)?
# ============================================================
print(f"\n{'='*70}")
print("PFAFFIAN CONNECTION")
print("=" * 70)

for n in [4, 6]:
    print(f"\nn={n}:")
    m = n * (n-1) // 2
    count = min(1 << m, 50)

    for bits in range(count):
        T = tournament_from_bits(n, bits)
        S = skew_adj(T)
        h = hamiltonian_path_count(T)

        det_S = np.linalg.det(S)
        pf_sq = det_S  # det(S) = Pf(S)^2 for skew-symmetric
        pf = np.sqrt(abs(pf_sq)) * (1 if pf_sq >= 0 else 1j)

        if bits < 15:
            print(f"  bits={bits}: H={h}, det(S)={det_S:.0f}, |Pf|={abs(pf):.1f}")

# ============================================================
# Test 5: The REAL question — A*A^T eigenvalues and H
# For tournaments, A*A^T has interesting structure.
# Off-diagonal: (A*A^T)[i][j] = |{k : i->k, k->j}| = common successors
# Diagonal: (A*A^T)[i][i] = out-degree of i
# ============================================================
print(f"\n{'='*70}")
print("A*A^T EIGENVALUES AND H")
print("=" * 70)

for n in [5, 6]:
    print(f"\nn={n}:")
    m = n * (n-1) // 2

    data = []
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        A = np.array(T, dtype=float)
        AAT = A @ A.T
        h = hamiltonian_path_count(T)
        eigs = sorted(np.linalg.eigvalsh(AAT), reverse=True)
        data.append((h, eigs, bits))

    # Sort by H
    data.sort(key=lambda x: -x[0])

    print(f"  Top 5 by H:")
    for h, eigs, bits in data[:5]:
        print(f"    H={h:3d}, eigs(AA^T)={[f'{e:.2f}' for e in eigs[:4]]}")

    print(f"  Bottom 5 by H:")
    for h, eigs, bits in data[-5:]:
        print(f"    H={h:3d}, eigs(AA^T)={[f'{e:.2f}' for e in eigs[:4]]}")

    # Check: does lambda_1(AA^T) correlate with H?
    h_vals = [d[0] for d in data]
    lam1_vals = [d[1][0] for d in data]
    corr = np.corrcoef(h_vals, lam1_vals)[0, 1]
    print(f"  Correlation(H, lambda_1(AA^T)) = {corr:.4f}")

    # Check spectral gap
    gaps = [d[1][0] - d[1][1] for d in data]
    corr_gap = np.corrcoef(h_vals, gaps)[0, 1]
    print(f"  Correlation(H, spectral_gap(AA^T)) = {corr_gap:.4f}")

print("\nDone.")
