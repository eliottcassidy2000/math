"""
dihedral_ladder_analysis.py — The complete dihedral ladder for tournament maximization.

THE LADDER (odd n only — even n has no circulant tournaments):
  n=3:  D_6,  1 spectral param,   2 circulants
  n=5:  D_10, 2 spectral params,  4 circulants
  n=7:  D_14, 3 spectral params,  8 circulants
  n=9:  D_18, 4 spectral params, 16 circulants  [COMPOSITE: 9=3^2]
  n=11: D_22, 5 spectral params, 32 circulants
  n=13: D_26, 6 spectral params, 64 circulants

KEY QUESTION: Does QR_n (when valid) always maximize H?
  - n prime, p=3 mod 4: QR_p is Paley, confirmed maximizer
  - n=9 (composite, -1 not in QR_9): QR_9 = {1,4,7} — does it maximize?
  - n prime, p=1 mod 4: QR_p is NOT valid (self-paired), cyclic interval wins

EVEN n IMPOSSIBILITY: For even n, the element n/2 satisfies n/2 = n - n/2,
so i→(i+n/2) AND (i+n/2)→i simultaneously. No circulant tournament exists.
This is WHY the dihedral groups at even n "interlace" differently.

Author: kind-pasteur-2026-03-12-S56
"""
import sys
import cmath
import math
from collections import defaultdict

sys.path.insert(0, '04-computation')
from satake_ndrt_h import all_circulant_H


def circulant_eigenvalues(n, S):
    omega = cmath.exp(2j * cmath.pi / n)
    return [sum(omega ** (k * s) for s in S) for k in range(n)]


def qr_set_general(n):
    """Quadratic residues mod n (coprime to n only)."""
    from math import gcd
    residues = set()
    for a in range(1, n):
        if gcd(a, n) == 1:
            residues.add(pow(a, 2, n))
    return sorted(residues - {0})


def is_valid_tournament_set(n, S):
    """Check if S is a valid tournament connection set for Z_n."""
    S_set = set(S)
    for s in S:
        neg_s = (n - s) % n
        if neg_s == 0:
            continue
        if neg_s in S_set and neg_s != s:
            # Both s and -s in S: not a tournament (bidirectional)
            return False
        if neg_s == s:
            # Self-paired: s = n-s, so 2s = n. For odd n, impossible.
            # For even n: s = n/2, creates bidirectional edge.
            return False
    # Check completeness: every pair {s, n-s} has exactly one in S
    for s in range(1, n):
        neg_s = (n - s) % n
        if neg_s == s:
            return False  # even n problem
        if (s in S_set) == (neg_s in S_set):
            return False  # both or neither
    return True


def main():
    print("=" * 70)
    print("THE DIHEDRAL LADDER — FULL SPECTRAL ANALYSIS")
    print("=" * 70)

    # ==================================================================
    # SECTION 1: Even n impossibility
    # ==================================================================
    print("\nSECTION 1: EVEN n CIRCULANT TOURNAMENT IMPOSSIBILITY")
    print("-" * 50)
    print("For even n, element n/2 is self-paired (n/2 = n - n/2).")
    print("This forces bidirectional or missing edge => no tournament.")
    for n in [2, 4, 6, 8]:
        print(f"  n={n}: element {n // 2} is self-paired. "
              f"D_{{{2 * n}}} acts on {n}-gon but NO circulant tournament exists.")

    # ==================================================================
    # SECTION 2: Full ladder for odd n = 3, 5, 7, 9, 11, 13
    # ==================================================================
    print(f"\n\n{'=' * 70}")
    print("SECTION 2: THE ODD LADDER")
    print("=" * 70)

    ladder_data = []

    for n in [3, 5, 7, 9, 11, 13]:
        print(f"\n{'=' * 60}")
        half = (n - 1) // 2
        n_tours = 2 ** half
        is_prime = all(n % i != 0 for i in range(2, n))

        # Check QR structure
        qr = qr_set_general(n)
        neg1 = n - 1
        neg1_is_qr = neg1 in qr
        qr_valid = not neg1_is_qr and is_valid_tournament_set(n, qr) if qr else False

        print(f"n={n}: D_{{{2 * n}}} ({2 * n} elements), {half} spectral params, "
              f"{n_tours} circulants")
        print(f"  {'PRIME' if is_prime else f'COMPOSITE ({n}=?)'}, "
              f"n mod 4 = {n % 4}, n mod 8 = {n % 8}")
        print(f"  QR_{n} = {qr}, -1 mod {n} = {neg1} "
              f"{'IN QR' if neg1_is_qr else 'NOT IN QR'}")
        if qr_valid:
            print(f"  => QR_{n} is a VALID tournament set (Paley-like)")
        else:
            print(f"  => QR_{n} is NOT a valid tournament set")

        # Compute all circulant H values
        results = all_circulant_H(n)
        H_max = results[0][0]
        H_min = results[-1][0]
        n_distinct_H = len(set(H for H, _ in results))

        # Find maximizer spectral data
        S_max = results[0][1]
        eigs_max = circulant_eigenvalues(n, S_max)
        y2_max = [eigs_max[k].imag ** 2 for k in range(1, half + 1)]
        spread_max = max(abs(eigs_max[k]) for k in range(1, n)) - min(abs(eigs_max[k]) for k in range(1, n))

        # Find flattest tournament
        best_spread = float('inf')
        S_flat = None
        H_flat = None
        for H, S in results:
            eigs = circulant_eigenvalues(n, S)
            mags = [abs(eigs[k]) for k in range(1, n)]
            sp = max(mags) - min(mags)
            if sp < best_spread:
                best_spread = sp
                S_flat = S
                H_flat = H

        # Find QR tournament H (if valid)
        H_qr = None
        if qr_valid:
            for H, S in results:
                if sorted(S) == sorted(qr):
                    H_qr = H
                    break

        # Cyclic interval
        ci = list(range(n - half, n))
        H_ci = None
        for H, S in results:
            if sorted(S) == sorted(ci):
                H_ci = H
                break

        print(f"\n  H range: [{H_min}, {H_max}], {n_distinct_H} distinct values")
        print(f"  Maximizer: S={S_max}, H={H_max}")
        print(f"    y^2 values: {['%.4f' % y for y in y2_max]}")
        print(f"    Spread: {spread_max:.4f}")
        print(f"  Flattest: S={S_flat}, H={H_flat}, spread={best_spread:.4f}")
        if H_qr is not None:
            qr_rank = sorted(set(H for H, _ in results), reverse=True).index(H_qr) + 1
            print(f"  QR tournament: S={sorted(qr)}, H={H_qr}, "
                  f"rank={qr_rank}/{n_distinct_H} "
                  f"{'*** MAXIMIZER ***' if H_qr == H_max else ''}")
        print(f"  Cyclic interval: S={ci}, H={H_ci}, "
              f"{'*** MAXIMIZER ***' if H_ci == H_max else ''}")

        # Check: flat = max?
        flat_is_max = H_flat == H_max
        qr_is_max = (H_qr == H_max) if H_qr is not None else False
        print(f"\n  FLATTEST = MAX H: {'YES' if flat_is_max else 'NO'}")
        if qr_valid:
            print(f"  QR = MAX H: {'YES' if qr_is_max else 'NO'}")

        ladder_data.append({
            'n': n, 'is_prime': is_prime, 'n_mod4': n % 4,
            'H_max': H_max, 'H_flat': H_flat, 'H_qr': H_qr, 'H_ci': H_ci,
            'S_max': S_max, 'flat_is_max': flat_is_max, 'qr_valid': qr_valid,
            'qr_is_max': qr_is_max, 'spread_max': spread_max,
            'best_spread': best_spread, 'n_distinct_H': n_distinct_H,
        })

    # ==================================================================
    # SECTION 3: Summary table
    # ==================================================================
    print(f"\n\n{'=' * 70}")
    print("SECTION 3: LADDER SUMMARY TABLE")
    print("=" * 70)
    print(f"\n{'n':>3} {'|D|':>4} {'prime':>5} {'mod4':>4} {'mod8':>4} "
          f"{'n_H':>4} {'H_max':>10} {'flat=max':>8} {'QR_valid':>8} "
          f"{'QR=max':>6} {'CI=max':>6} {'winner':<15}")
    print("-" * 95)

    for d in ladder_data:
        winner = "?"
        if d['qr_is_max']:
            winner = "PALEY/QR"
        elif d['H_ci'] == d['H_max']:
            winner = "CYCLIC INT"
        elif d['flat_is_max']:
            winner = "FLATTEST"
        else:
            winner = "OTHER"

        print(f"{d['n']:>3} {2 * d['n']:>4} {str(d['is_prime']):>5} "
              f"{d['n_mod4']:>4} {d['n'] % 8:>4} {d['n_distinct_H']:>4} "
              f"{d['H_max']:>10} {str(d['flat_is_max']):>8} "
              f"{str(d['qr_valid']):>8} "
              f"{str(d['qr_is_max']):>6} "
              f"{str(d['H_ci'] == d['H_max']):>6} "
              f"{winner:<15}")

    # ==================================================================
    # SECTION 4: n=9 deep dive (first composite case)
    # ==================================================================
    print(f"\n\n{'=' * 70}")
    print("SECTION 4: n=9 DEEP DIVE (COMPOSITE: 9 = 3^2)")
    print("=" * 70)

    n = 9
    results = all_circulant_H(n)
    half = (n - 1) // 2  # 4

    print(f"\nAll {len(results)} circulant tournaments on Z_9:")
    print(f"  {'S':<25} {'H':>8} {'spread':>8} {'y^2 values':<50}")
    print("  " + "-" * 93)

    for H, S in results:
        eigs = circulant_eigenvalues(n, S)
        mags = [abs(eigs[k]) for k in range(1, n)]
        spread = max(mags) - min(mags)
        y2 = [eigs[k].imag ** 2 for k in range(1, half + 1)]
        S_str = str(S)
        y2_str = str(['%.3f' % y for y in y2])

        # Check if QR
        is_qr = sorted(S) == sorted(qr_set_general(9))
        label = " QR" if is_qr else ""
        if sorted(S) == sorted(list(range(n - half, n))):
            label += " CI"

        print(f"  {S_str:<25} {H:>8} {spread:>8.4f} {y2_str:<50}{label}")

    # Orbit decomposition at n=9
    print(f"\n  Z_9^* orbit decomposition:")
    from math import gcd
    Z9_star = [a for a in range(1, 9) if gcd(a, 9) == 1]  # {1,2,4,5,7,8}
    print(f"  Z_9^* = {Z9_star} (order {len(Z9_star)})")

    H_to_sets = defaultdict(list)
    for H, S in results:
        H_to_sets[H].append(S)

    for H in sorted(H_to_sets.keys(), reverse=True):
        sets = H_to_sets[H]
        rep = sets[0]
        # Compute orbit under Z_9^*
        orbit = set()
        for a in Z9_star:
            aS = tuple(sorted((a * s) % 9 for s in rep))
            orbit.add(aS)
        print(f"  H={H}: {len(sets)} tournaments, orbit size {len(orbit)}")

    # ==================================================================
    # SECTION 5: QR maximization pattern
    # ==================================================================
    print(f"\n\n{'=' * 70}")
    print("SECTION 5: THE QR MAXIMIZATION THEOREM")
    print("=" * 70)
    print("""
    PATTERN FROM THE LADDER:
      n=3  (prime, mod4=3):  QR valid, QR = max       (Paley)
      n=5  (prime, mod4=1):  QR NOT valid              (all tied)
      n=7  (prime, mod4=3):  QR valid, QR = max        (Paley)
      n=9  (composite):      QR valid (since -1 not QR_9), QR = ???
      n=11 (prime, mod4=3):  QR valid, QR = max        (Paley)
      n=13 (prime, mod4=1):  QR NOT valid              (cyclic int wins)

    The KEY QUESTION: at n=9, does QR_9 = {1,4,7} maximize H?
    If YES: QR maximization extends to composites where -1 not in QR.
    If NO: Paley maximization is purely a prime phenomenon.
""")

    # Check n=9 QR
    qr9 = sorted(qr_set_general(9))
    H_qr9 = None
    for H, S in results:
        if sorted(S) == qr9:
            H_qr9 = H
            break
    H_max9 = results[0][0]
    print(f"  n=9: QR_9 = {qr9}, H(QR) = {H_qr9}, H_max = {H_max9}")
    if H_qr9 is None:
        print(f"  QR_9 NOT FOUND among circulant tournaments!")
        print(f"  REASON: |QR_9| = {len(qr9)} but need |S| = {(9-1)//2}")
        print(f"  QR tournaments ONLY work at PRIMES (|QR_p| = (p-1)/2 = required |S|)")
        print(f"  At n=p^k: |QR_n| = phi(n)/2 != (n-1)/2 in general.")
    elif H_qr9 == H_max9:
        print(f"  QR = MAX: YES!")
    else:
        rank = sorted(set(H for H, _ in results), reverse=True).index(H_qr9) + 1
        print(f"  QR = MAX: NO, rank {rank}/{len(set(H for H, _ in results))}")
        print(f"  Gap: {H_max9 - H_qr9}")

    # Check spectral flatness of QR at n=9
    eigs_qr = circulant_eigenvalues(9, qr9)
    y2_qr = [eigs_qr[k].imag ** 2 for k in range(1, half + 1)]
    mags_qr = [abs(eigs_qr[k]) for k in range(1, 9)]
    spread_qr = max(mags_qr) - min(mags_qr)
    print(f"  QR spectrum: y^2 = {['%.4f' % y for y in y2_qr]}")
    print(f"  QR spread = {spread_qr:.4f}")
    is_flat = len(set(round(y, 6) for y in y2_qr)) == 1
    print(f"  QR is spectrally FLAT: {is_flat}")

    # ==================================================================
    # SECTION 6: The unified principle
    # ==================================================================
    print(f"\n\n{'=' * 70}")
    print("SECTION 6: THE UNIFIED SPECTRAL PRINCIPLE")
    print("=" * 70)
    print("""
    Combining all evidence:

    1. FLAT SPECTRUM => MAX H (when achievable):
       - n=3,7,11: Paley achieves flat, Paley = max. [CONFIRMED]
       - n=9: QR_9 achieves flat (if it does) => QR = max? [TEST ABOVE]

    2. WHEN FLAT NOT ACHIEVABLE:
       - n=5: All spectra identical (degenerate case), all tied. [TRIVIAL]
       - n=13: No flat spectrum. Concentrated (Dirichlet kernel) = max. [CONFIRMED]

    3. THE PRINCIPLE: H is maximized at the "spectral center" (flat point).
       When center is achievable: flat wins (Paley/QR).
       When center is not achievable: landscape reverses, concentration wins.

    4. ACHIEVABILITY: flat spectrum <=> -1 not in QR_n AND S=QR_n works.
       For prime p: -1 not in QR_p iff p = 3 mod 4.
       For composite n=p^2: -1 not in QR_n iff p = 3 mod 4.

    5. THE DIHEDRAL INTERLACING:
       EVEN n: No circulant tournaments exist (self-paired n/2 element).
       ODD n: Circulant tournaments exist, D_{2n} acts.
       The even-order dihedral groups D_{2n} provide the spectral framework.
       The odd-order rotation groups Z_n provide the eigenvalue decomposition.
       Between each tournament size (odd n), the "interlaced" group mediates.
""")


if __name__ == '__main__':
    main()
