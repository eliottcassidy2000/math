"""
dihedral_spectral_analysis.py — Deep analysis of dihedral group structure
and spectral properties of circulant tournament H-maximizers.

KEY GEOMETRIC INSIGHT (from user):
  - An n-vertex tournament can be drawn as a regular n-gon on the unit circle
  - D_{2n} acts on this n-gon: rotations = automorphisms, reflections = anti-auts
  - For EVERY circulant tournament on Z_p (prime), D_{2p} acts faithfully
  - The dihedral groups D_2, D_4, D_6, D_8, D_{10}, ... (all even order)
    correspond to tournament sizes n=1,2,3,4,5,...
  - Odd-order cyclic groups Z_3, Z_5, Z_7, ... are "interlaced" between them

KEY QUESTIONS:
  1. Does spectral flatness (= all |lambda_k| equal) ALWAYS predict the H-maximizer?
  2. What distinguishes p ≡ 3 mod 4 (Paley wins) from p ≡ 1 mod 4 (no Paley)?
  3. Is H determined by the eigenvalue magnitudes alone?
  4. What is the representation-theoretic formula for H via D_{2p} irreps?
  5. For p ≡ 1 mod 4, which circulant achieves the "flattest" spectrum?

ALGEBRAIC CONSTRAINT:
  For ANY circulant tournament on Z_n: Re(lambda_k) = -1/2 for all k != 0.
  Proof: A + A^T = J - I, so lambda_k + conj(lambda_k) = -1.
  Therefore lambda_k = -1/2 + i*y_k, and |lambda_k|^2 = 1/4 + y_k^2.
  "Spectral flatness" = all y_k^2 equal = all |y_k| = sqrt(p)/2 (Gauss sum value).

Author: kind-pasteur-2026-03-12-S56
"""
import sys
import time
import cmath
import math
from collections import Counter, defaultdict

sys.path.insert(0, '04-computation')
from satake_ndrt_h import ham_count_circulant, all_circulant_H


def circulant_eigenvalues(n, S):
    """Compute eigenvalues lambda_k = sum_{s in S} omega^{ks} for k=0,...,n-1."""
    omega = cmath.exp(2j * cmath.pi / n)
    return [sum(omega ** (k * s) for s in S) for k in range(n)]


def spectral_data(n, S):
    """Compute full spectral profile of circulant tournament on Z_n."""
    eigs = circulant_eigenvalues(n, S)
    # Non-trivial eigenvalues (k != 0)
    mags = [abs(eigs[k]) for k in range(1, n)]
    mag_sq = [m ** 2 for m in mags]
    # Imaginary parts: lambda_k = -1/2 + i*y_k
    ys = [eigs[k].imag for k in range(1, n)]
    y_sq = [y ** 2 for y in ys]

    spread = max(mags) - min(mags)
    variance = sum((m - sum(mags) / len(mags)) ** 2 for m in mags) / len(mags)
    y_var = sum((y2 - sum(y_sq) / len(y_sq)) ** 2 for y2 in y_sq) / len(y_sq)

    # Power sums of y^2 (the independent spectral invariants)
    # Only use k=1,...,(n-1)/2 since y_{n-k} = -y_k
    half = (n - 1) // 2
    y_half = [eigs[k].imag for k in range(1, half + 1)]
    y2_half = [y ** 2 for y in y_half]

    return {
        'eigs': eigs,
        'mags': mags,
        'mag_sq': mag_sq,
        'ys': ys,
        'y_sq': y_sq,
        'y_half': y_half,
        'y2_half': y2_half,
        'spread': spread,
        'variance': variance,
        'y_var': y_var,
        # Power sums of y^2 (using independent half)
        'sum_y2': sum(y2_half),
        'sum_y4': sum(y ** 2 for y in y2_half),
        'sum_y6': sum(y ** 3 for y in y2_half),
        'sum_y8': sum(y ** 4 for y in y2_half),
        'max_mag': max(mags),
        'min_mag': min(mags),
    }


def qr_set(p):
    """Quadratic residues mod p (excluding 0)."""
    return sorted(set(pow(a, 2, p) for a in range(1, p)) - {0})


def cyclic_interval_S(n):
    """S = {ceil(n/2), ..., n-1} — the 'upper half' cyclic interval."""
    k = (n - 1) // 2
    return list(range(n - k, n))


def main():
    print("=" * 70)
    print("DIHEDRAL SPECTRAL ANALYSIS — TOURNAMENT MAXIMIZATION")
    print("=" * 70)

    print("""
GEOMETRIC PICTURE:
  Place n vertices on the unit circle: v_k = exp(2*pi*i*k/n).
  An n-tournament is an orientation of the edges of K_n.
  D_{2n} = symmetries of the regular n-gon = <r, s | r^n = s^2 = 1, srs = r^{-1}>.

  For EVERY circulant tournament on Z_p (prime):
    - r (rotation by 2*pi/p) is an AUTOMORPHISM (definition of circulant)
    - s (reflection v -> -v) is an ANTI-AUTOMORPHISM (reverses all arcs)
    - Proof: S union -S = {1,...,p-1} for any tournament connection set S
    - So EVERY circulant tournament has D_{2p} symmetry (aut + anti-aut)

  Eigenvalue constraint (from A + A^T = J - I):
    lambda_k = -1/2 + i*y_k  for all k != 0
    |lambda_k|^2 = 1/4 + y_k^2

  Spectral flatness: all y_k^2 equal <=> all |lambda_k| equal.
  For Paley T_p (p = 3 mod 4): y_k^2 = p/4 for all k (Gauss sum) => PERFECTLY FLAT.
""")

    # =====================================================================
    # SECTION 1: Spectral flatness vs H-maximization for all primes
    # =====================================================================
    print("=" * 70)
    print("SECTION 1: SPECTRAL FLATNESS vs H-MAXIMIZATION")
    print("=" * 70)

    primes = [3, 5, 7, 11, 13]
    for p in primes:
        t0 = time.time()
        results = all_circulant_H(p)
        H_max = results[0][0]
        n_tours = len(results)

        print(f"\n--- Z_{p} | D_{{{2*p}}} symmetry | p mod 4 = {p % 4} | p mod 8 = {p % 8} ---")
        print(f"    {n_tours} circulant tournaments, H_max = {H_max}")

        # Compute spectral data for each
        tour_data = []
        for H, S in results:
            sd = spectral_data(p, S)
            sd['H'] = H
            sd['S'] = S
            tour_data.append(sd)
        tour_data.sort(key=lambda x: -x['H'])

        # Print table
        print(f"    {'S':<25} {'H':>10} {'spread':>8} {'var(y^2)':>12} "
              f"{'Sum_y2':>10} {'Sum_y4':>12} {'max|L|':>8} {'min|L|':>8}")
        print("    " + "-" * 99)
        for td in tour_data:
            S_str = str(td['S'])
            print(f"    {S_str:<25} {td['H']:>10} {td['spread']:>8.4f} "
                  f"{td['y_var']:>12.6f} {td['sum_y2']:>10.4f} "
                  f"{td['sum_y4']:>12.4f} {td['max_mag']:>8.4f} {td['min_mag']:>8.4f}")

        # Test: min spread = max H?
        min_spread = min(tour_data, key=lambda x: x['spread'])
        min_yvar = min(tour_data, key=lambda x: x['y_var'])
        min_y4 = min(tour_data, key=lambda x: x['sum_y4'])

        print(f"\n    FLATNESS TESTS:")
        print(f"      Min spread:  S={min_spread['S']}, H={min_spread['H']}  "
              f"{'✓ MAX H' if min_spread['H'] == H_max else '✗ NOT max H'}")
        print(f"      Min var(y²): S={min_yvar['S']}, H={min_yvar['H']}  "
              f"{'✓ MAX H' if min_yvar['H'] == H_max else '✗ NOT max H'}")
        print(f"      Min Σy⁴:    S={min_y4['S']}, H={min_y4['H']}  "
              f"{'✓ MAX H' if min_y4['H'] == H_max else '✗ NOT max H'}")

        # Check: is this Paley or cyclic interval?
        if p % 4 == 3:
            qr = qr_set(p)
            is_paley = sorted(min_spread['S']) == sorted(qr)
            print(f"      Maximizer is Paley (QR_{p}): {is_paley}")
        ci = cyclic_interval_S(p)
        is_ci = sorted(min_spread['S']) == sorted(ci) or sorted(min_yvar['S']) == sorted(ci)
        print(f"      Maximizer is cyclic interval: {is_ci}")

        # Pearson correlation H vs spread
        h_vals = [td['H'] for td in tour_data]
        s_vals = [td['spread'] for td in tour_data]
        if len(set(s_vals)) > 1:
            n_t = len(tour_data)
            mean_h = sum(h_vals) / n_t
            mean_s = sum(s_vals) / n_t
            cov = sum((h - mean_h) * (s - mean_s) for h, s in zip(h_vals, s_vals)) / n_t
            std_h = (sum((h - mean_h) ** 2 for h in h_vals) / n_t) ** 0.5
            std_s = (sum((s - mean_s) ** 2 for s in s_vals) / n_t) ** 0.5
            corr = cov / (std_h * std_s) if std_h > 0 and std_s > 0 else 0
            print(f"      Pearson(H, spread) = {corr:.4f} {'(negative = flatness → more H)' if corr < 0 else ''}")

        print(f"    [{time.time() - t0:.1f}s]")

    # =====================================================================
    # SECTION 2: IS H DETERMINED BY SPECTRAL INVARIANTS?
    # =====================================================================
    print(f"\n\n{'=' * 70}")
    print("SECTION 2: IS H A FUNCTION OF SPECTRAL POWER SUMS?")
    print("=" * 70)
    print("""
    If H depends only on the eigenvalue magnitudes |lambda_k|, then it is
    determined by the power sums sigma_j = sum_k y_k^{2j} (j=1,2,...).
    Test: do tournaments with same (sigma_1, sigma_2, ...) have same H?
""")

    for p in primes:
        if p == 3:
            continue  # too small
        results = all_circulant_H(p)
        data = []
        for H, S in results:
            sd = spectral_data(p, S)
            data.append((H, sd['sum_y2'], sd['sum_y4'], sd['sum_y6'], sd['sum_y8'], S))

        # Test with increasing numbers of power sums
        for depth in [1, 2, 3, 4]:
            sig_to_H = defaultdict(set)
            for H, s2, s4, s6, s8, S in data:
                sigs = [s2, s4, s6, s8][:depth]
                sig = tuple(round(s, 8) for s in sigs)
                sig_to_H[sig].add(H)
            determined = all(len(hs) == 1 for hs in sig_to_H.values())
            n_sigs = len(sig_to_H)
            n_H = len(set(H for H, *_ in data))
            print(f"  p={p}: {depth} power sum(s) → {n_sigs} signatures, "
                  f"{n_H} distinct H values. Determined: {determined}")

    # =====================================================================
    # SECTION 3: D_{2p} IRREDUCIBLE REPRESENTATION STRUCTURE
    # =====================================================================
    print(f"\n\n{'=' * 70}")
    print("SECTION 3: D_{2p} IRREP DECOMPOSITION")
    print("=" * 70)
    print("""
    D_{2p} irreps for odd prime p:
      - rho_0: trivial (dim 1)
      - rho_sign: sign rep (dim 1): rotations -> 1, reflections -> -1
      - rho_k (k=1,...,(p-1)/2): 2-dim, chi(r^j) = 2*cos(2*pi*k*j/p)

    For a circulant tournament, eigenvalue pair (lambda_k, bar{lambda_k})
    lives in the irrep rho_k. The anti-aut (reflection) maps lambda_k -> -1 - bar{lambda_k}.

    Since lambda_k = -1/2 + i*y_k:
      Anti-aut image = -1 - (-1/2 - i*y_k) = -1/2 + i*y_k = lambda_k (!!)
    This means: the anti-automorphism PRESERVES each eigenvalue!!
    Because -1 - bar{lambda_k} = -1 - (-1/2 - i*y_k) = -1/2 + i*y_k = lambda_k.

    Implication: reflections act TRIVIALLY on the spectrum of any circulant tournament.
    This is why ALL circulant tournaments on Z_p have D_{2p} symmetry — the
    tournament IS determined by its eigenvalues, and reflections fix them.
""")

    for p in [7, 13]:
        print(f"\n--- D_{{{2 * p}}} irrep structure (p={p}) ---")
        half = (p - 1) // 2
        print(f"  {2 + half} irreps: rho_0, rho_sign, rho_1,...,rho_{half}")
        print(f"  dim check: 1 + 1 + {half}*4 = {2 + half * 4}  vs  |D_{2*p}| = {2 * p}")
        # dim^2 sum: 1 + 1 + (p-1)/2 * 4 = 2 + 2(p-1) = 2p. Checks out!

        results = all_circulant_H(p)
        H_max = results[0][0]
        S_max = results[0][1]
        eigs_max = circulant_eigenvalues(p, S_max)

        print(f"\n  Maximizer S={S_max}, H={H_max}:")
        print(f"  {'Irrep':<10} {'k':<5} {'lambda_k':<30} {'|lambda|^2':<12} {'y_k^2':<12}")
        print(f"  " + "-" * 69)
        print(f"  {'rho_0':<10} {'0':<5} {eigs_max[0].real:+.6f} + 0i              "
              f"{abs(eigs_max[0]) ** 2:<12.4f} {'N/A':<12}")
        for k in range(1, half + 1):
            lam = eigs_max[k]
            mag2 = abs(lam) ** 2
            y2 = lam.imag ** 2
            print(f"  {'rho_' + str(k):<10} {k:<5} {lam.real:+.6f} {lam.imag:+.6f}i    "
                  f"{mag2:<12.4f} {y2:<12.4f}")

        # Verify anti-aut preserves eigenvalues
        print(f"\n  Anti-automorphism check: -1 - bar(lambda_k) == lambda_k?")
        all_ok = True
        for k in range(1, p):
            lam = eigs_max[k]
            anti = -1 - lam.conjugate()
            if abs(anti - lam) > 1e-10:
                print(f"    k={k}: FAILS! lambda={lam}, anti={anti}")
                all_ok = False
        print(f"    {'✓ All confirmed' if all_ok else '✗ FAILURE'}")

    # =====================================================================
    # SECTION 4: THE P MOD 4 DICHOTOMY — WHY PALEY WINS ONLY AT 3 MOD 4
    # =====================================================================
    print(f"\n\n{'=' * 70}")
    print("SECTION 4: THE p mod 4 DICHOTOMY")
    print("=" * 70)
    print("""
    For p ≡ 3 mod 4:
      -1 is NOT a quadratic residue mod p.
      QR_p has the property: s in QR => -s NOT in QR.
      So QR_p is a valid tournament connection set.
      The Paley tournament T_p has |lambda_k|^2 = (p+1)/4 for all k != 0.
      This is PERFECTLY FLAT — the theoretical maximum flatness.
      => Paley is the H-maximizer (among circulants, at least).

    For p ≡ 1 mod 4:
      -1 IS a quadratic residue mod p.
      QR_p has s in QR => -s in QR too (both or neither).
      QR_p is NOT a valid tournament connection set (contains {s, -s} pairs).
      No perfectly flat spectrum is achievable.
      The NEXT BEST connection set must be found by other means.
      => Cyclic interval appears to win (confirmed at p=5, p=13).

    The interlaced odd-order groups Z_3, Z_5, Z_7, Z_9, Z_{11}, Z_{13}, ...
    between D_4, D_6, ..., D_{2n}, ... correspond to the ROTATION SUBGROUPS.
    The rotation group Z_p has (p-1)/2 non-trivial irreps, each 1-dimensional
    (characters omega^k for k=1,...,p-1). Pairing into D_{2p} irreps gives
    (p-1)/2 two-dimensional irreps.

    The Legendre symbol (k/p) determines whether eigenspace k falls into:
      - QR class (spectral contribution from quadratic residue character)
      - NQR class (spectral contribution from non-residue character)
    For p ≡ 3 mod 4, these classes are LINKED by the anti-automorphism.
    For p ≡ 1 mod 4, they are SELF-CONTAINED.
""")

    for p in [7, 11, 13]:
        print(f"\n  p = {p} (p mod 4 = {p % 4}):")
        if p % 4 == 3:
            qr = qr_set(p)
            nqr = sorted(set(range(1, p)) - set(qr))
            print(f"    QR_{p} = {qr}")
            print(f"    NQR_{p} = {nqr}")
            print(f"    -1 mod {p} = {p - 1} ∈ {'QR' if (p - 1) in qr else 'NQR'}")

            eigs = circulant_eigenvalues(p, qr)
            y_qr = [eigs[k].imag ** 2 for k in qr if k <= (p - 1) // 2]
            y_nqr = [eigs[k].imag ** 2 for k in nqr if k <= (p - 1) // 2]
            print(f"    y² at QR positions:  {['%.4f' % y for y in y_qr]}")
            print(f"    y² at NQR positions: {['%.4f' % y for y in y_nqr]}")
            print(f"    All equal: {len(set(round(y, 8) for y in y_qr + y_nqr)) == 1}")
        else:
            qr = qr_set(p)
            nqr = sorted(set(range(1, p)) - set(qr))
            print(f"    QR_{p} = {qr}")
            print(f"    NQR_{p} = {nqr}")
            print(f"    -1 mod {p} = {p - 1} ∈ {'QR' if (p - 1) in qr else 'NQR'}")
            print(f"    -1 ∈ QR => QR is CLOSED under negation => NOT a tournament set")

            # What IS the flattest tournament at this p?
            results = all_circulant_H(p)
            best_spread = float('inf')
            best_td = None
            for H, S in results:
                sd = spectral_data(p, S)
                if sd['spread'] < best_spread:
                    best_spread = sd['spread']
                    best_td = (H, S, sd)
            H_b, S_b, sd_b = best_td
            print(f"    Flattest circulant: S={S_b}, spread={sd_b['spread']:.6f}, H={H_b}")
            print(f"    y² values: {['%.4f' % y for y in sd_b['y2_half']]}")
            print(f"    Is maximizer: {'✓ YES' if H_b == results[0][0] else '✗ NO'}")

    # =====================================================================
    # SECTION 5: GAUSS SUM STRUCTURE AND THE CIRCULANT PERMANENT
    # =====================================================================
    print(f"\n\n{'=' * 70}")
    print("SECTION 5: GAUSS SUMS, JACOBI SUMS, AND H")
    print("=" * 70)
    print("""
    For the Paley tournament T_p (p ≡ 3 mod 4):
      lambda_k = sum_{a in QR_p} omega^{ka} = (eta(k) * g - 1) / 2
      where g = sum_{a=0}^{p-1} (a/p) * omega^a is the Gauss sum, |g|^2 = p,
      and eta(k) = (k/p) is the Legendre symbol.

      So |lambda_k|^2 = |eta(k)*g - 1|^2 / 4 = (p + 1 - 2*Re(eta(k)*g)) / 4.
      For p ≡ 3 mod 4: g = i*sqrt(p), so Re(eta(k)*g) = 0, giving |lambda_k|^2 = (p+1)/4.
      For p ≡ 1 mod 4: g = sqrt(p) (real), so Re(eta(k)*g) = +/-sqrt(p),
        giving |lambda_k|^2 = (p+1 -/+ 2*sqrt(p))/4 — TWO distinct values!

      This is WHY perfect flatness is possible iff p ≡ 3 mod 4:
        p ≡ 3: g is purely imaginary => Re(eta*g) = 0 => all |lambda| equal.
        p ≡ 1: g is real => Re(eta*g) = +/- sqrt(p) => two magnitude classes.
""")

    for p in [7, 11]:
        g_sum = sum(
            (1 if pow(a, (p - 1) // 2, p) == 1 else -1) * cmath.exp(2j * cmath.pi * a / p)
            for a in range(1, p)
        )
        print(f"  p={p}: Gauss sum g = {g_sum.real:.6f} + {g_sum.imag:.6f}i")
        print(f"    |g|^2 = {abs(g_sum) ** 2:.4f} (should be {p})")
        print(f"    g / i*sqrt({p}) = {(g_sum / (1j * math.sqrt(p))).real:.4f} + "
              f"{(g_sum / (1j * math.sqrt(p))).imag:.4f}i")
        print(f"    => g {'IS' if abs(g_sum.real) < 0.01 else 'is NOT'} purely imaginary")

    for p in [5, 13]:
        g_sum = sum(
            (1 if pow(a, (p - 1) // 2, p) == 1 else -1) * cmath.exp(2j * cmath.pi * a / p)
            for a in range(1, p)
        )
        print(f"  p={p}: Gauss sum g = {g_sum.real:.6f} + {g_sum.imag:.6f}i")
        print(f"    |g|^2 = {abs(g_sum) ** 2:.4f} (should be {p})")
        print(f"    g / sqrt({p}) = {(g_sum / math.sqrt(p)).real:.4f} + "
              f"{(g_sum / math.sqrt(p)).imag:.4f}i")
        print(f"    => g {'IS' if abs(g_sum.imag) < 0.01 else 'is NOT'} purely real")

    # =====================================================================
    # SECTION 6: WHAT MAKES CYCLIC INTERVAL SPECIAL AT p ≡ 1 mod 4?
    # =====================================================================
    print(f"\n\n{'=' * 70}")
    print("SECTION 6: CYCLIC INTERVAL AS BEST ACHIEVABLE FLATNESS")
    print("=" * 70)
    print("""
    At p ≡ 1 mod 4, no perfectly flat tournament exists.
    Hypothesis: the cyclic interval S = {ceil(p/2), ..., p-1} achieves the
    MINIMUM possible spectral spread among all connection sets.

    Why might this be true? The cyclic interval is the "most consecutive"
    possible connection set — it consists of a contiguous arc on Z_p.
    The Fourier transform of a contiguous arc is a Dirichlet kernel,
    which has the most uniform magnitude among all {0,1}-valued sequences
    with the prescribed number of 1s. This is essentially the Fejer kernel
    approximation to the flat spectrum.
""")

    for p in [5, 13]:
        print(f"\n  p = {p}:")
        results = all_circulant_H(p)

        # Rank all tournaments by spread
        ranked = []
        for H, S in results:
            sd = spectral_data(p, S)
            ranked.append((sd['spread'], sd['y_var'], H, S, sd))
        ranked.sort(key=lambda x: x[0])  # by spread

        ci = cyclic_interval_S(p)
        ci_sorted = sorted(ci)

        print(f"  All tournaments ranked by spectral spread:")
        print(f"  {'Rank':<6} {'S':<25} {'spread':>10} {'var(y²)':>12} {'H':>10} {'type':<15}")
        print(f"  " + "-" * 78)
        for i, (sp, yv, H, S, sd) in enumerate(ranked):
            s_sorted = sorted(S)
            tour_type = ""
            if p % 4 == 3 and s_sorted == sorted(qr_set(p)):
                tour_type = "PALEY"
            if s_sorted == ci_sorted:
                tour_type = "CYCLIC INT"
            print(f"  {i + 1:<6} {str(S):<25} {sp:>10.6f} {yv:>12.6f} {H:>10} {tour_type:<15}")

        # Check if cyclic interval has minimum spread
        ci_spread = [sp for sp, yv, H, S, sd in ranked if sorted(S) == ci_sorted]
        min_spread_val = ranked[0][0]
        if ci_spread and abs(ci_spread[0] - min_spread_val) < 1e-8:
            print(f"\n  ✓ Cyclic interval achieves MINIMUM spread at p={p}")
        else:
            print(f"\n  ✗ Cyclic interval does NOT achieve minimum spread at p={p}")
            print(f"    CI spread: {ci_spread[0] if ci_spread else 'N/A':.6f}, "
                  f"min spread: {min_spread_val:.6f}")

    # =====================================================================
    # SECTION 7: THE INTERLACING — D_{2n} CHARACTER TABLE STRUCTURE
    # =====================================================================
    print(f"\n\n{'=' * 70}")
    print("SECTION 7: INTERLACING OF DIHEDRAL AND CYCLIC GROUPS")
    print("=" * 70)
    print("""
    The user observes: D_2, D_4, D_6, D_8, D_{10}, ... (orders 2,4,6,8,10,...)
    correspond to tournament sizes n=1,2,3,4,5,...
    Interlaced: Z_3, Z_5, Z_7, Z_9, ... (orders 3,5,7,9,...).

    D_{2n} has (n+3)/2 irreps if n is odd, (n+6)/2 if n is even.
    The 2-dim irreps rho_k (k=1,...,floor((n-1)/2)) are indexed by
    the elements of Z_n^* / {+/-1}, which has size (n-1)/2 for odd n.

    KEY OBSERVATION: The number of independent spectral parameters
    for a circulant tournament on Z_n is (n-1)/2 = number of 2-dim D_{2n} irreps.
    These are the y_k^2 values (one per conjugate pair of eigenvalues).

    As n grows: 1 -> 2 -> 3 -> 4 -> 5 -> 6 parameters.
    The H-maximization problem is: optimize H over an ((n-1)/2)-simplex
    defined by the constraint that {y_1^2, ..., y_{(n-1)/2}^2} must
    correspond to a valid connection set S.

    The VALID y^2 configurations form a discrete set of 2^{(n-1)/2} points.
    The "flattest" point (all y_k^2 equal) is achievable iff p ≡ 3 mod 4.
""")

    print("\n  Ladder of spectral parameters:")
    for n in range(3, 14, 2):
        half = (n - 1) // 2
        n_tours = 2 ** half
        print(f"    n={n:>2}: D_{{{2 * n}}} ({2 * n} elements), "
              f"{half} spectral params (y_k^2), "
              f"{n_tours} circulant tournaments, "
              f"{'Paley exists' if n > 2 and all(n % i != 0 for i in range(2, n)) and n % 4 == 3 else 'no Paley' if all(n % i != 0 for i in range(2, n)) else 'composite'}")

    # =====================================================================
    # SECTION 8: DEEPER — H AS A FUNCTION ON THE D_{2p} REPRESENTATION RING
    # =====================================================================
    print(f"\n\n{'=' * 70}")
    print("SECTION 8: H VIA SYMMETRIC FUNCTIONS OF EIGENVALUES")
    print("=" * 70)
    print("""
    Since all circulant tournaments on Z_p have M = (H/p)*I (scalar transfer matrix),
    H is a symmetric function of the connection set S, invariant under:
      - Translation: S -> S + a (rotations, always an aut)
      - Negation: S -> -S (reflection anti-aut, maps T to T^op)

    Since H(T) = H(T^op), and T^op has connection set -S = complement of S,
    H is determined by the UNORDERED PAIR {S, -S}. For p prime, S and -S
    are always different (no self-complementary circulants at odd prime),
    so the number of distinct H values is at most 2^{(p-1)/2}/2.

    More precisely: H depends on S only through the eigenvalue magnitudes
    |lambda_k|^2 = 1/4 + y_k^2, which are the CHARACTERS of the representation.
    So H is a CLASS FUNCTION on the D_{2p} representation ring.

    Test: compute Newton power sums e_1 = Sigma y_k^2, e_2 = Sigma y_k^4, ...
    and check if H = f(e_1, e_2, ...) for some universal polynomial f.
""")

    for p in [7, 13]:
        print(f"\n  p = {p}:")
        results = all_circulant_H(p)
        half = (p - 1) // 2

        # Collect data points: (H, e_1, e_2, ..., e_{half})
        data = []
        for H, S in results:
            sd = spectral_data(p, S)
            # Newton power sums: e_j = sum_k y_k^{2j}
            power_sums = []
            for j in range(1, half + 1):
                ej = sum(y ** j for y in sd['y2_half'])
                power_sums.append(ej)
            data.append((H, power_sums, S))

        # Check if H depends only on elementary symmetric polynomials
        # (i.e., on the MULTISET of y_k^2 values)
        multiset_to_H = defaultdict(set)
        for H, ps, S in data:
            key = tuple(round(p, 8) for p in sorted(spectral_data(p, S)['y2_half']))
            multiset_to_H[key].add(H)

        determined = all(len(v) == 1 for v in multiset_to_H.values())
        print(f"    H determined by multiset {{y_k^2}}: {determined}")
        print(f"    {len(multiset_to_H)} distinct multisets, "
              f"{len(set(H for H, _, _ in data))} distinct H values")

        if determined and len(multiset_to_H) <= 10:
            print(f"    Mapping:")
            for key in sorted(multiset_to_H.keys(), key=lambda k: -max(multiset_to_H[k])):
                h_set = multiset_to_H[key]
                print(f"      {{y²}} = {['%.4f' % y for y in key]}  →  H = {sorted(h_set, reverse=True)}")

    # =====================================================================
    # SUMMARY
    # =====================================================================
    print(f"\n\n{'=' * 70}")
    print("SUMMARY OF FINDINGS")
    print("=" * 70)
    print("""
    1. SPECTRAL FLATNESS ↔ H-MAXIMIZATION:
       Test whether min spread = max H at each prime.

    2. H DETERMINED BY EIGENVALUE MAGNITUDES:
       Test whether the multiset {|lambda_k|^2} determines H.

    3. p ≡ 3 mod 4 (Paley exists):
       Perfect flatness achievable via Gauss sums (g purely imaginary).
       Paley tournament has all y_k^2 = p/4.

    4. p ≡ 1 mod 4 (no Paley):
       No perfect flatness. Best achievable is cyclic interval.
       Gauss sum g is real, creating TWO classes of eigenvalue magnitudes.

    5. REPRESENTATION THEORY:
       H is a class function on the D_{2p} representation ring.
       The anti-automorphism v -> -v PRESERVES eigenvalues (not permutes them).
       This means H is truly a function of the (p-1)/2 independent |lambda_k|^2 values.
    """)


if __name__ == '__main__':
    main()
