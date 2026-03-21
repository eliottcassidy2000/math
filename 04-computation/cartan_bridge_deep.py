#!/usr/bin/env python3
"""
cartan_bridge_deep.py -- kind-pasteur-2026-03-21-S13

Deep investigation of the Cartan bridge between tournament theory and
the decomposition gl(n,R) = so(n) + p + R.

KEY THREADS:
1. S_T^2 structure: For Paley, S_T^2 = -p*I + J (uniform off-diagonal).
   What does S_T^2 look like for other tournaments?

2. Transfer matrix M(T) in the Cartan decomposition:
   M is symmetric (THM-030), tr(M) = H.
   M = (H/n)*I + M_traceless.  How does ||M_traceless|| relate to H?

3. Powers of S_T: S_T^{2k} is symmetric, S_T^{2k+1} is antisymmetric.
   The alternation between Cartan sectors encodes tournament dynamics.

4. The "Casimir deviation" of S_T^2: D(T) = S_T^2 - (-n*I + J) measures
   how far T is from being doubly regular (Paley-like).

5. Can we express H(T) in terms of eigenvalues of S_T?

Author: kind-pasteur-2026-03-21-S13
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict


def binary_to_tournament(bits, n):
    """Convert bit encoding to adjacency matrix."""
    A = np.zeros((n, n), dtype=int)
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def signed_adjacency(A):
    """
    Signed adjacency matrix S_T:
    S_T[i,j] = +1 if i->j, -1 if j->i, 0 if i=j.
    """
    n = A.shape[0]
    S = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            S[i][j] = 2 * A[i][j] - 1
    return S


def count_ham_paths(A):
    n = A.shape[0]
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[(mask, v)] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))


def transfer_matrix(A):
    """
    Compute M[a,b] = number of Hamiltonian paths from a to b.
    """
    n = A.shape[0]
    M = np.zeros((n, n), dtype=int)
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[(mask, v)] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << n) - 1
    for a in range(n):
        for b in range(n):
            # Count paths from a to b
            M[a][b] = dp.get((full, b), 0)
    # Wait - the DP above mixes all starting points. Need separate DP per start.
    # Redo properly.
    M = np.zeros((n, n), dtype=int)
    for start in range(n):
        dp2 = defaultdict(int)
        dp2[(1 << start, start)] = 1
        for mask in range(1, 1 << n):
            if not (mask & (1 << start)):
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                if dp2[(mask, v)] == 0:
                    continue
                for w in range(n):
                    if mask & (1 << w):
                        continue
                    if A[v][w]:
                        dp2[(mask | (1 << w), w)] += dp2[(mask, v)]
        full = (1 << n) - 1
        for end in range(n):
            M[start][end] = dp2.get((full, end), 0)
    return M


def paley_tournament(p):
    """Paley tournament on p vertices."""
    qr = set()
    for k in range(1, p):
        qr.add((k*k) % p)
    A = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                A[i][j] = 1
    return A


# ============================================================
# INVESTIGATION 1: S_T^2 Structure Theorem
# ============================================================

def investigate_ST2(n=7, max_tournaments=None):
    """
    For each tournament on n vertices, compute S_T^2 and analyze.

    For Paley: S_T^2 = -n*I + J (uniform off-diagonal = 1).
    Question: How much does S_T^2 vary for non-Paley tournaments?

    The off-diagonal variation of S_T^2 measures "doubly regular deviation."
    """
    print("=" * 70)
    print(f"INVESTIGATION 1: S_T^2 STRUCTURE (n={n})")
    print("=" * 70)

    m = n * (n - 1) // 2
    total = 1 << m
    if max_tournaments is None:
        max_tournaments = total if total <= 100000 else 10000

    # Collect statistics
    off_diag_stats = defaultdict(int)  # (min, max) of off-diagonal elements
    S2_types = defaultdict(list)  # signature -> list of H values

    np.random.seed(42)
    for trial in range(min(max_tournaments, total)):
        if total <= 100000:
            bits = trial
        else:
            bits = np.random.randint(0, total)

        A = binary_to_tournament(bits, n)
        S = signed_adjacency(A)
        S2 = S @ S

        # Off-diagonal elements
        off_diag = []
        for i in range(n):
            for j in range(n):
                if i != j:
                    off_diag.append(S2[i][j])

        od_min = min(off_diag)
        od_max = max(off_diag)
        od_set = tuple(sorted(set(off_diag)))
        H = count_ham_paths(A)

        off_diag_stats[(od_min, od_max)] += 1
        # Signature: sorted list of unique off-diagonal values
        S2_types[od_set].append(H)

    print(f"\nOff-diagonal range distribution of S_T^2:")
    for (mn, mx), count in sorted(off_diag_stats.items()):
        print(f"  range [{mn}, {mx}]: {count} tournaments")

    print(f"\nOff-diagonal value signature types:")
    for sig, H_list in sorted(S2_types.items(), key=lambda x: -len(x[1])):
        print(f"  values {sig}: {len(H_list)} tournaments, "
              f"H range [{min(H_list)}, {max(H_list)}], "
              f"mean H = {np.mean(H_list):.1f}")

    # For Paley specifically
    if n in [3, 5, 7, 11, 13]:
        if all(n % i != 0 for i in range(2, int(n**0.5)+1)):
            A_paley = paley_tournament(n)
            S_paley = signed_adjacency(A_paley)
            S2_paley = S_paley @ S_paley
            print(f"\nPaley T_{n} specifically:")
            print(f"  S_T^2 diagonal: all {S2_paley[0,0]}")
            print(f"  S_T^2 off-diagonal: all {S2_paley[0,1]}")
            expected = -n * np.eye(n) + np.ones((n, n))
            if np.allclose(S2_paley, expected):
                print(f"  CONFIRMED: S_T^2 = -{n}*I + J")
            else:
                print(f"  S_T^2 differs from -n*I + J!")


# ============================================================
# INVESTIGATION 2: Transfer Matrix in Cartan
# ============================================================

def investigate_transfer_cartan(n=5):
    """
    The transfer matrix M(T) is symmetric (THM-030) with tr(M) = H.
    Decompose M = (H/n)*I + M_traceless.

    Question: How does ||M_traceless|| / ||M|| relate to tournament structure?
    For Paley: M_traceless = 0 (M is scalar).
    For non-Paley: M_traceless != 0.

    The "scalarness" of M measures how close the tournament is to
    distributing Hamiltonian paths uniformly across all start-end pairs.
    """
    print("\n" + "=" * 70)
    print(f"INVESTIGATION 2: TRANSFER MATRIX CARTAN DECOMPOSITION (n={n})")
    print("=" * 70)

    m = n * (n - 1) // 2
    total = 1 << m

    results = []
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        H = count_ham_paths(A)
        M = transfer_matrix(A)

        # Verify symmetry (THM-030)
        sym_err = np.max(np.abs(M - M.T))
        if sym_err > 0:
            print(f"  WARNING: M not symmetric at bits={bits}, err={sym_err}")

        # Cartan decomposition of M
        scalar = np.trace(M) / n  # = H/n
        M_traceless = M - scalar * np.eye(n)
        norm_traceless = np.linalg.norm(M_traceless, 'fro')
        norm_total = np.linalg.norm(M, 'fro')
        scalarness = 1.0 - (norm_traceless**2 / norm_total**2) if norm_total > 0 else 1.0

        # Eigenvalues of M
        eigvals = sorted(np.linalg.eigvalsh(M), reverse=True)

        # Score sequence
        scores = tuple(sorted(A.sum(axis=1)))

        results.append({
            'bits': bits, 'H': H, 'scalarness': scalarness,
            'eigvals': eigvals, 'scores': scores,
            'norm_traceless': norm_traceless
        })

    # Sort by H
    results.sort(key=lambda r: r['H'])

    print(f"\nTransfer Matrix scalarness by H value:")
    print(f"  {'H':>5s} {'scalarness':>11s} {'eigvals':>30s} {'scores':>20s}")
    print(f"  {'-'*5} {'-'*11} {'-'*30} {'-'*20}")

    # Group by H
    by_H = defaultdict(list)
    for r in results:
        by_H[r['H']].append(r)

    for H in sorted(by_H.keys()):
        group = by_H[H]
        avg_scalar = np.mean([r['scalarness'] for r in group])
        # Show one representative
        rep = group[0]
        eigvals_str = ','.join(f"{e:.0f}" for e in rep['eigvals'])
        print(f"  {H:5d} {avg_scalar:11.4f} [{eigvals_str:>28s}] {str(rep['scores']):>20s}")

    print(f"\n  scalarness = 1 means M = (H/n)*I (perfectly uniform path distribution)")
    print(f"  scalarness < 1 means some start-end pairs have more paths than others")

    # Correlation between scalarness and H
    H_vals = [r['H'] for r in results]
    scalar_vals = [r['scalarness'] for r in results]
    corr = np.corrcoef(H_vals, scalar_vals)[0, 1]
    print(f"\n  Correlation(H, scalarness) = {corr:.4f}")

    return results


# ============================================================
# INVESTIGATION 3: Powers of S_T (Cartan sector alternation)
# ============================================================

def investigate_power_evolution(n=7):
    """
    S_T^k alternates between Cartan sectors:
    - k odd: S_T^k is antisymmetric (in so(n))
    - k even: S_T^k is symmetric (in p)

    But the SPECTRAL structure evolves. Track how eigenvalues change.
    For Paley: eigenvalues of S_T are {0, +/-i*sqrt(p)} (multiplicity (p-1)/2 each).
    So S_T^k has eigenvalues {0, (-p)^{k/2}} for k even, {0, +/-i*(-p)^{(k-1)/2}*sqrt(p)} for k odd.
    """
    print("\n" + "=" * 70)
    print(f"INVESTIGATION 3: POWER EVOLUTION OF S_T (n={n})")
    print("=" * 70)

    # Paley tournament
    if n in [3, 5, 7, 11, 13] and all(n % i != 0 for i in range(2, int(n**0.5)+1)):
        A = paley_tournament(n)
    else:
        # Use a random regular tournament if n isn't prime
        np.random.seed(42)
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i+1, n):
                if np.random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

    S = signed_adjacency(A).astype(float)
    H = count_ham_paths(A) if n <= 12 else None

    print(f"\nTournament: {'Paley T_' + str(n) if n in [3,5,7,11,13] else 'n=' + str(n)}")
    if H is not None:
        print(f"H = {H}")

    # Eigenvalues of S_T
    eigvals_S = np.linalg.eigvals(S)
    print(f"\nEigenvalues of S_T (sorted by |lambda|):")
    for ev in sorted(eigvals_S, key=lambda x: abs(x)):
        if abs(ev.imag) < 1e-10:
            print(f"  {ev.real:+.4f}")
        else:
            print(f"  {ev.real:+.4f} {ev.imag:+.4f}i  |lambda| = {abs(ev):.4f}")

    # Verify S_T^2 = -n*I + J for Paley
    S2 = S @ S
    S2_expected = -n * np.eye(n) + np.ones((n, n))
    if np.allclose(S2, S2_expected, atol=0.01):
        print(f"\n  CONFIRMED: S_T^2 = -{n}*I + J for this tournament")

    # Powers evolution
    print(f"\nPower evolution (Frobenius norm and symmetry type):")
    print(f"  {'k':>3s} {'||S^k||_F':>12s} {'type':>12s} {'tr(S^k)/n':>12s} {'max|offdiag|':>14s}")
    S_k = np.eye(n)
    for k in range(1, min(n+2, 15)):
        S_k = S_k @ S
        norm_fro = np.linalg.norm(S_k, 'fro')
        # Check symmetry type
        asym_err = np.linalg.norm(S_k + S_k.T, 'fro')
        sym_err = np.linalg.norm(S_k - S_k.T, 'fro')
        if asym_err < 0.01:
            stype = "antisymmetric"
        elif sym_err < 0.01:
            stype = "symmetric"
        else:
            stype = "mixed"
        tr_norm = np.trace(S_k) / n
        max_off = max(abs(S_k[i][j]) for i in range(n) for j in range(n) if i != j)
        print(f"  {k:3d} {norm_fro:12.2f} {stype:>12s} {tr_norm:12.2f} {max_off:14.2f}")

    # Key observation: for Paley with eigenvalues {0, +/-i*sqrt(p)}:
    # S^2 has eigenvalues {0, -p} -> all real
    # S^3 has eigenvalues {0, +/-i*p*sqrt(p)} -> imaginary
    # S^4 has eigenvalues {0, p^2} -> real positive
    # S^(p-1) has eigenvalues {0, (-p)^{(p-1)/2}} = {0, p^{(p-1)/2} * (-1)^{(p-1)/2}}
    # For p=7: S^6 has eigenvalue 7^3 * (-1)^3 = -343
    # For p=7: S^7 should be... S_T^7 = (-7)^3 * S_T for Paley? Let me check.
    #
    # S_T has min poly: x(x^2 + p) = 0 (since eigenvalues 0, +/-i*sqrt(p))
    # So S_T^3 = -p * S_T. This is the key recurrence!
    if n in [3, 5, 7, 11, 13]:
        S3 = S @ S @ S
        check = S3 + n * S
        if np.max(np.abs(check)) < 0.01:
            print(f"\n  KEY IDENTITY VERIFIED: S_T^3 = -{n} * S_T")
            print(f"  Minimal polynomial of S_T: x(x^2 + {n}) = 0")
            print(f"  Eigenvalues: 0, +/- i*sqrt({n})")
        else:
            print(f"\n  S_T^3 != -{n} * S_T (tournament is not doubly regular)")
            # What IS the minimal polynomial?
            # Cayley-Hamilton: S satisfies its characteristic polynomial
            # For S with eigenvalues lambda_k: char poly = product (x - lambda_k)
            pass


# ============================================================
# INVESTIGATION 4: Casimir Deviation and H
# ============================================================

def investigate_casimir_deviation(n=5):
    """
    For a tournament T, define the Casimir deviation:
    D(T) = S_T^2 - (-n*I + J)

    For doubly regular (Paley): D = 0.
    For general: D measures off-diagonal non-uniformity.

    Question: Does ||D(T)|| predict H(T)?
    """
    print("\n" + "=" * 70)
    print(f"INVESTIGATION 4: CASIMIR DEVIATION AND H (n={n})")
    print("=" * 70)

    m = n * (n - 1) // 2
    total = 1 << m

    results = []
    target = -n * np.eye(n) + np.ones((n, n))

    for bits in range(total):
        A = binary_to_tournament(bits, n)
        S = signed_adjacency(A).astype(float)
        S2 = S @ S
        D = S2 - target
        H = count_ham_paths(A)
        scores = tuple(sorted(A.sum(axis=1)))

        # D should have zero diagonal (since S2[i,i] = -(n-1) and target[i,i] = -(n-1))
        norm_D = np.linalg.norm(D, 'fro')

        # Off-diagonal of D
        off_D = []
        for i in range(n):
            for j in range(n):
                if i != j:
                    off_D.append(D[i][j])

        results.append({
            'bits': bits, 'H': H, 'norm_D': norm_D,
            'scores': scores, 'D_range': (min(off_D), max(off_D)),
            'D_variance': np.var(off_D)
        })

    # Correlation
    H_vals = [r['H'] for r in results]
    D_vals = [r['norm_D'] for r in results]
    corr = np.corrcoef(H_vals, D_vals)[0, 1]

    print(f"\n  Correlation(H, ||D||) = {corr:.4f}")
    print(f"  D = 0 means doubly regular (Paley-like)")
    print(f"  Negative correlation would mean: more regular -> more H paths")

    # Group by scores
    by_scores = defaultdict(list)
    for r in results:
        by_scores[r['scores']].append(r)

    print(f"\n  Results by score sequence:")
    print(f"  {'scores':>20s} {'count':>6s} {'mean H':>8s} {'mean ||D||':>10s}")
    for scores in sorted(by_scores.keys()):
        group = by_scores[scores]
        print(f"  {str(scores):>20s} {len(group):6d} "
              f"{np.mean([r['H'] for r in group]):8.1f} "
              f"{np.mean([r['norm_D'] for r in group]):10.4f}")

    # Within regular tournaments (if any)
    if n % 2 == 1:
        reg_score = tuple([n//2] * n) if n % 2 == 1 else None
        if reg_score and reg_score in by_scores:
            reg_group = by_scores[reg_score]
            print(f"\n  REGULAR tournaments ({len(reg_group)} total):")
            for r in sorted(reg_group, key=lambda x: -x['H']):
                print(f"    H={r['H']:5d}, ||D||={r['norm_D']:.4f}, "
                      f"D range={r['D_range']}")

    return results


# ============================================================
# INVESTIGATION 5: H from eigenvalues of S_T
# ============================================================

def investigate_H_from_spectrum(n=5):
    """
    Can we express H(T) as a function of the eigenvalues of S_T?

    For Paley T_p: eigenvalues are {0, +/-i*sqrt(p)}, H(T_p) is known.
    For general T: eigenvalues are complex conjugate pairs + possibly 0.

    The characteristic polynomial of S_T encodes the cycle structure
    of T (by Newton identities). H is related to Hamiltonian paths
    which are NOT directly cycle counts, but OCF connects them via
    independent sets of cycles.

    If H = f(eigenvalues of S_T), this would be a SPECTRAL formula for H.
    """
    print("\n" + "=" * 70)
    print(f"INVESTIGATION 5: H FROM SPECTRUM OF S_T (n={n})")
    print("=" * 70)

    m = n * (n - 1) // 2
    total = 1 << m

    results = []
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        S = signed_adjacency(A).astype(float)
        H = count_ham_paths(A)
        eigvals = sorted(np.linalg.eigvals(S), key=lambda x: (abs(x), x.real))
        # Store magnitudes squared (since eigenvalues come in conjugate pairs)
        mags_sq = sorted([abs(ev)**2 for ev in eigvals], reverse=True)

        results.append({
            'H': H, 'eigvals': eigvals, 'mags_sq': tuple(round(m, 4) for m in mags_sq),
            'det_S': round(np.real(np.linalg.det(S)), 4),
            'tr_S2': round(np.real(np.trace(S @ S)), 4),
            'tr_S4': round(np.real(np.trace(np.linalg.matrix_power(S, 4))), 4),
        })

    # Group by eigenvalue magnitudes
    by_spec = defaultdict(list)
    for r in results:
        by_spec[r['mags_sq']].append(r['H'])

    print(f"\n  Distinct spectra: {len(by_spec)}")
    print(f"\n  Do isospectral tournaments have the same H?")
    agree = 0
    disagree = 0
    for spec, H_list in by_spec.items():
        if len(set(H_list)) == 1:
            agree += 1
        else:
            disagree += 1

    print(f"    Spectral classes with unique H: {agree}")
    print(f"    Spectral classes with multiple H: {disagree}")

    if disagree > 0:
        print(f"\n  CONCLUSION: H is NOT a spectral invariant (counterexamples exist)")
        print(f"  Isospectral tournaments with different H:")
        for spec, H_list in by_spec.items():
            if len(set(H_list)) > 1:
                print(f"    |lambda|^2 = {spec}: H values = {sorted(set(H_list))}")
                if disagree > 10:
                    break
    else:
        print(f"\n  CONCLUSION: H IS a spectral invariant at n={n}!")
        print(f"  (Every isospectral class has a unique H value)")

    # Newton identity analysis: power sums vs H
    print(f"\n  Power sum analysis:")
    print(f"  {'tr(S^2)':>10s} {'tr(S^4)':>10s} {'det(S)':>10s} {'H range':>15s}")
    by_powersums = defaultdict(list)
    for r in results:
        key = (r['tr_S2'], r['tr_S4'], r['det_S'])
        by_powersums[key].append(r['H'])

    for key in sorted(by_powersums.keys()):
        H_list = by_powersums[key]
        tr2, tr4, det = key
        H_range = f"[{min(H_list)}, {max(H_list)}]" if min(H_list) != max(H_list) else f"{min(H_list)}"
        print(f"  {tr2:10.1f} {tr4:10.1f} {det:10.1f} {H_range:>15s}")

    return results


# ============================================================
# INVESTIGATION 6: The Grand Cartan Theorem
# ============================================================

def grand_cartan_theorem(n=7):
    """
    CONJECTURE (The Grand Cartan Theorem):

    For a tournament T on n vertices, the transfer matrix M(T) satisfies:

    M(T) = f(S_T^2)

    where f is a POLYNOMIAL and S_T is the signed adjacency matrix.

    Why this might be true:
    - S_T^2 is symmetric (like M)
    - S_T^2 encodes 2-path structure (like M which encodes n-1 path structure)
    - For Paley: S_T^2 = -p*I + J and M = (H/p)*I, so f(x) = (H/p)*something
    - The transfer matrix M counts Hamiltonian paths, which are related to
      the permanent of S_T, which is related to S_T^k products

    Actually for Paley: S_T^3 = -p * S_T, so the entire ring generated by S_T
    is 3-dimensional: {I, S_T, S_T^2}. And M = 27*I for T_7.
    So M = f(S_T^2) = 27*I, and f(-7*I + J) = 27*I.

    For non-Paley: the ring generated by S_T is larger (full matrix algebra).
    But M might still be expressible as a polynomial in S_T^2.

    Let's test this.
    """
    print("\n" + "=" * 70)
    print(f"INVESTIGATION 6: GRAND CARTAN THEOREM TEST (n={n})")
    print("=" * 70)

    if n > 8:
        print("  n too large for exhaustive test, using sampling")
        return

    m = n * (n - 1) // 2
    total = 1 << m

    # For each tournament, check if M can be expressed as polynomial in S^2
    agrees = 0
    disagrees = 0

    num_test = min(total, 5000)
    np.random.seed(42)

    for trial in range(num_test):
        if total <= 5000:
            bits = trial
        else:
            bits = np.random.randint(0, total)

        A = binary_to_tournament(bits, n)
        S = signed_adjacency(A).astype(float)
        S2 = S @ S
        M = transfer_matrix(A).astype(float)

        # Check if M is in the algebra generated by S2
        # The algebra generated by an n×n matrix X has dimension at most n
        # (by Cayley-Hamilton). So M is in Alg(S2) iff M is a linear
        # combination of {I, S2, S2^2, ..., S2^{n-1}}.
        #
        # Build the basis vectors (flattened)
        basis = []
        X_k = np.eye(n)
        for k in range(n):
            basis.append(X_k.flatten())
            X_k = X_k @ S2

        basis = np.array(basis).T  # (n^2) x n
        M_flat = M.flatten()

        # Solve basis @ coeffs = M_flat (least squares)
        try:
            coeffs, residual, rank, sv = np.linalg.lstsq(basis, M_flat, rcond=None)
            reconstruction = basis @ coeffs
            error = np.linalg.norm(M_flat - reconstruction)
            if error < 0.01:
                agrees += 1
            else:
                disagrees += 1
                if disagrees <= 3:
                    print(f"  Counter: bits={bits}, H={int(np.trace(M))}, "
                          f"error={error:.4f}")
        except Exception as e:
            print(f"  Error at bits={bits}: {e}")
            disagrees += 1

    print(f"\n  Results: {agrees}/{num_test} tournaments have M in Alg(S_T^2)")
    print(f"  Failures: {disagrees}")

    if disagrees == 0:
        print(f"  CONJECTURE SUPPORTED at n={n}!")
        print(f"  M(T) IS a polynomial in S_T^2 for all tested tournaments")
    else:
        print(f"  CONJECTURE FAILS at n={n}")
        print(f"  M(T) is NOT always in Alg(S_T^2)")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    # Investigation 1: S_T^2 structure
    investigate_ST2(n=5)
    investigate_ST2(n=7, max_tournaments=10000)

    # Investigation 2: Transfer matrix Cartan
    investigate_transfer_cartan(n=5)

    # Investigation 3: Power evolution (Paley T_7)
    investigate_power_evolution(n=7)

    # Investigation 4: Casimir deviation
    investigate_casimir_deviation(n=5)
    investigate_casimir_deviation(n=7)

    # Investigation 5: H from spectrum
    investigate_H_from_spectrum(n=5)

    # Investigation 6: Grand Cartan Theorem
    grand_cartan_theorem(n=5)
    grand_cartan_theorem(n=7)
