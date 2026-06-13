"""
spectral_tournament_knot.py -- kind-pasteur-2026-03-14-S69
Spectral theory of tournament adjacency matrices.

KEY CONNECTIONS:
1. Tournament adjacency A has eigenvalues related to H(T) how?
2. The skew-symmetric part S = A - A^T has purely imaginary eigenvalues
3. For even n, det(S) = Pfaffian(S)^2, and Pfaffian is always odd
4. Trace(A^k) = #{directed k-walks} = k * c_k(T) for k odd
5. The characteristic polynomial det(xI - A) -- how does it relate to H?

CREATIVE IDEAS:
- Define "tournament Jones" as a specialization of det(xI - A)?
- Connect eigenvalue distribution to H via OCF
- Use the "spectral gap" to predict H-maximization
- Random matrix theory for tournaments
"""

import numpy as np
from collections import Counter, defaultdict
import sys

sys.stdout.reconfigure(encoding='utf-8')

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

def compute_H_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def char_poly_coeffs(A):
    """Compute characteristic polynomial coefficients of A."""
    n = len(A)
    # Use numpy eigenvalues and reconstruct
    # det(xI - A) = sum_k c_k x^k
    eigs = np.linalg.eigvals(A.astype(float))
    # Coefficients from roots
    from numpy.polynomial import polynomial as P
    coeffs = np.array([1.0])
    for e in eigs:
        coeffs = np.convolve(coeffs, [1.0, -e])
    # Return as real rounded (exact integers for integer matrices)
    return [round(c.real) for c in coeffs]

def main():
    print("=" * 70)
    print("SPECTRAL TOURNAMENT THEORY")
    print("kind-pasteur-2026-03-14-S69")
    print("=" * 70)

    # ========================================
    # PART 1: Characteristic polynomial vs H
    # ========================================
    print("\n" + "=" * 70)
    print("PART 1: CHARACTERISTIC POLYNOMIAL det(xI - A) vs H(T)")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        charpoly_by_H = defaultdict(list)
        det_by_H = defaultdict(list)

        for bits in range(2**(n*(n-1)//2)):
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)

            # Eigenvalues
            eigs = np.linalg.eigvals(A.astype(float))
            eigs_sorted = sorted(eigs, key=lambda x: -x.real)

            det_A = round(np.linalg.det(A.astype(float)))
            tr_A = int(np.trace(A))  # = sum of out-degrees / 2... no, = sum of diagonal = 0

            charpoly_by_H[H].append(tuple(round(e.real, 3) + round(e.imag, 3)*1j for e in eigs_sorted))
            det_by_H[H].append(det_A)

        for H in sorted(charpoly_by_H.keys()):
            dets = sorted(set(det_by_H[H]))
            eig_patterns = len(set(charpoly_by_H[H]))
            # Max real eigenvalue
            max_re = max(max(e.real for e in eigs) for eigs in charpoly_by_H[H])
            min_re = min(min(e.real for e in eigs) for eigs in charpoly_by_H[H])
            print(f"  H={H:3d}: det(A) in {dets}, {eig_patterns} eig patterns, max_re={max_re:.3f}")

    # ========================================
    # PART 2: Spectral gap and H-maximization
    # ========================================
    print("\n" + "=" * 70)
    print("PART 2: SPECTRAL GAP vs H")
    print("  'Spectral gap' = largest real eigenvalue of A")
    print("  Expected: H-maximizers have largest spectral gap")
    print("=" * 70)

    for n in [5, 6]:
        print(f"\n--- n = {n} ---")
        data = []

        count = 0
        for bits in range(2**(n*(n-1)//2)):
            count += 1
            if n >= 6 and count > 5000:
                break

            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)

            eigs = np.linalg.eigvals(A.astype(float))
            spectral_radius = max(abs(e) for e in eigs)
            max_real = max(e.real for e in eigs)

            # Also compute: permanent approximation via eigenvalues
            # van der Waerden: perm(A) >= n!/n^n for doubly stochastic
            # For tournaments: A is 0/1 with row sums = out-degrees

            data.append((H, spectral_radius, max_real))

        H_vals = np.array([d[0] for d in data])
        sr_vals = np.array([d[1] for d in data])
        mr_vals = np.array([d[2] for d in data])

        corr_sr = np.corrcoef(H_vals, sr_vals)[0, 1]
        corr_mr = np.corrcoef(H_vals, mr_vals)[0, 1]

        print(f"  Correlation(H, spectral_radius) = {corr_sr:.6f}")
        print(f"  Correlation(H, max_real_eigenvalue) = {corr_mr:.6f}")

        # H-maximizer eigenvalues
        max_H = max(H_vals)
        max_data = [(d[1], d[2]) for d in data if d[0] == max_H]
        print(f"  H-maximizer (H={max_H}): spectral_radius in {sorted(set(d[0] for d in max_data))[:3]}")

        # H-minimizer eigenvalues
        min_H = min(H_vals)
        min_data = [(d[1], d[2]) for d in data if d[0] == min_H]
        print(f"  H-minimizer (H={min_H}): spectral_radius in {sorted(set(d[0] for d in min_data))[:3]}")

    # ========================================
    # PART 3: Skew-symmetric spectrum
    # ========================================
    print("\n" + "=" * 70)
    print("PART 3: SKEW-SYMMETRIC SPECTRUM of S = A - A^T")
    print("  S has purely imaginary eigenvalues (skew-sym)")
    print("  For even n: Pfaffian(S) is always odd (HYP-1171)")
    print("=" * 70)

    for n in [4, 5, 6]:
        print(f"\n--- n = {n} ---")
        pfaff_by_H = defaultdict(list)
        imag_eigs_by_H = defaultdict(list)

        count = 0
        for bits in range(2**(n*(n-1)//2)):
            count += 1
            if n >= 6 and count > 3000:
                break

            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            S = A - A.T  # skew-symmetric, entries +-1

            eigs_S = np.linalg.eigvals(S.astype(float))
            # Should be purely imaginary
            imag_parts = sorted([abs(e.imag) for e in eigs_S if abs(e.imag) > 0.001], reverse=True)

            if n % 2 == 0:
                pfaff_sq = round(np.linalg.det(S.astype(float)))
                pfaff = round(pfaff_sq**0.5) if pfaff_sq > 0 else 0
                pfaff_by_H[H].append(pfaff)

            imag_eigs_by_H[H].append(tuple(round(x, 2) for x in imag_parts[:3]))

        if n % 2 == 0:
            print(f"  H -> Pfaffian(A-A^T):")
            for H in sorted(pfaff_by_H.keys()):
                vals = sorted(set(pfaff_by_H[H]))
                print(f"    H={H:3d}: Pfaff in {vals}")

        print(f"\n  H -> largest imaginary eigenvalues of A-A^T:")
        for H in sorted(imag_eigs_by_H.keys()):
            vals = sorted(set(imag_eigs_by_H[H]))
            print(f"    H={H:3d}: top imag eigs: {vals[:3]}")

    # ========================================
    # PART 4: NEW INVARIANT — "Tournament Determinant"
    # ========================================
    print("\n" + "=" * 70)
    print("PART 4: 'TOURNAMENT DETERMINANT' det(A) vs H")
    print("  For a knot, the determinant = |Alexander(-1)| = |det(V-V^T)|")
    print("  For tournament: det(A) is our 'tournament determinant'")
    print("=" * 70)

    for n in [3, 4, 5, 6]:
        print(f"\n--- n = {n} ---")
        det_by_H = defaultdict(list)

        count = 0
        for bits in range(2**(n*(n-1)//2)):
            count += 1
            if n >= 6 and count > 5000:
                break

            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            det_A = round(np.linalg.det(A.astype(float)))
            det_by_H[H].append(det_A)

        print(f"  H -> det(A):")
        for H in sorted(det_by_H.keys()):
            vals = sorted(set(det_by_H[H]))
            count_vals = len(det_by_H[H])
            print(f"    H={H:3d}: det in {vals[:10]}{'...' if len(vals)>10 else ''}, count={count_vals}")

        # Does det(A) determine H?
        det_to_H = defaultdict(set)
        for H, dets in det_by_H.items():
            for d in dets:
                det_to_H[d].add(H)
        ambiguous = sum(1 for hs in det_to_H.values() if len(hs) > 1)
        print(f"\n  det(A) ambiguous for H: {ambiguous}/{len(det_to_H)} det values")

    # ========================================
    # PART 5: Eigenvalue statistics of random tournaments
    # ========================================
    print("\n" + "=" * 70)
    print("PART 5: EIGENVALUE STATISTICS (Random Matrix Theory)")
    print("  Tournament A is a random {0,1} matrix with A+A^T = J-I")
    print("  What is the eigenvalue distribution?")
    print("=" * 70)

    np.random.seed(42)
    n = 8
    num_samples = 1000

    all_real_eigs = []
    all_imag_eigs = []
    H_vals = []

    for _ in range(num_samples):
        bits = np.random.randint(0, 2**(n*(n-1)//2))
        A = bits_to_adj(bits, n)
        H = compute_H_dp(A, n)
        H_vals.append(H)

        eigs = np.linalg.eigvals(A.astype(float))
        for e in eigs:
            all_real_eigs.append(e.real)
            all_imag_eigs.append(e.imag)

    print(f"  n={n}, {num_samples} samples")
    print(f"  Eigenvalue real parts: mean={np.mean(all_real_eigs):.4f}, "
          f"std={np.std(all_real_eigs):.4f}, "
          f"range=[{min(all_real_eigs):.4f}, {max(all_real_eigs):.4f}]")
    print(f"  Eigenvalue imag parts: mean={np.mean(all_imag_eigs):.4f}, "
          f"std={np.std(all_imag_eigs):.4f}, "
          f"range=[{min(all_imag_eigs):.4f}, {max(all_imag_eigs):.4f}]")

    print(f"\n  Expected: mean real = (n-1)/2 = {(n-1)/2:.1f} (row sum average)")
    print(f"  Note: A has row sums = out-degrees, trace = 0")
    print(f"  So sum(eigenvalues) = trace(A) = 0, but eigenvalues not centered at 0")

    # Actually: trace(A) = 0 (diagonal is 0), so sum of eigenvalues = 0
    # But the mean real part can be nonzero because complex eigenvalues come in conjugate pairs
    # and contribute 0 to imaginary part but nonzero to real part

    # Correlation between spectral properties and H
    spec_rad = []
    max_real = []
    for _ in range(num_samples):
        bits = np.random.randint(0, 2**(n*(n-1)//2))
        A = bits_to_adj(bits, n)
        eigs = np.linalg.eigvals(A.astype(float))
        spec_rad.append(max(abs(e) for e in eigs))
        max_real.append(max(e.real for e in eigs))

    print(f"\n  Correlation(H, spectral_radius) at n={n}: {np.corrcoef(H_vals[:num_samples], spec_rad)[0,1]:.4f}")
    print(f"  Correlation(H, max_real_eig) at n={n}: {np.corrcoef(H_vals[:num_samples], max_real)[0,1]:.4f}")
    print(f"  Mean H = {np.mean(H_vals):.2f}, n!/2^(n-1) = {np.math.factorial(n)/2**(n-1):.2f}")

    # ========================================
    # PART 6: det(A) for Paley tournaments
    # ========================================
    print("\n" + "=" * 70)
    print("PART 6: det(A) FOR PALEY TOURNAMENTS")
    print("=" * 70)

    for p in [3, 7, 11]:
        qr = set()
        for x in range(1, p):
            qr.add((x*x) % p)
        A = np.zeros((p, p), dtype=int)
        for i in range(p):
            for j in range(p):
                if i != j and ((j-i) % p) in qr:
                    A[i][j] = 1

        det_A = np.linalg.det(A.astype(float))
        eigs = np.linalg.eigvals(A.astype(float))
        eigs_sorted = sorted(eigs, key=lambda x: -abs(x))

        print(f"\n  Paley T_{p}:")
        print(f"    det(A) = {round(det_A)}")
        print(f"    Eigenvalues: {[round(e.real, 4) + round(e.imag, 4)*1j for e in eigs_sorted]}")

        # Spectral radius
        sr = max(abs(e) for e in eigs)
        print(f"    Spectral radius = {sr:.4f}")
        print(f"    (p-1)/2 = {(p-1)/2}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)

if __name__ == '__main__':
    import math
    np.math = math  # Fix numpy math attribute
    main()
