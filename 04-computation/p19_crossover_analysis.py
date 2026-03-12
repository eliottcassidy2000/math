"""
p19_crossover_analysis.py — WHY does the interval beat Paley at p=19?

FACTS (THM-135):
  - Paley wins at p=3 (tie), 7 (+7.4%), 11 (+2.2%)
  - Interval wins at p=19 (+1.0%)
  - Paley margin SHRINKS monotonically: 7.4% -> 2.2% -> -1.0%

SPECTRAL COMPARISON:
  - Paley: |lambda_k| = sqrt(p)/2 for all k != 0 (FLAT)
  - Interval: |lambda_1| ~ p/pi (CONCENTRATED), rest O(1)

THIS SCRIPT INVESTIGATES:
  1. Cycle count comparison (c_3, c_5 via trace; structure of advantage)
  2. How the Paley margin varies with p (analytical formula?)
  3. The dihedral/geometric interpretation
  4. Eigenvalue phase structure at p=7,11,19

Author: kind-pasteur-2026-03-12-S56c
"""

import math
import cmath
import sys

sys.path.insert(0, '04-computation')


def is_qr(a, p):
    if a % p == 0:
        return False
    return pow(a, (p - 1) // 2, p) == 1


def eigenvalues(S, p):
    """Compute eigenvalues of circulant tournament adjacency matrix."""
    omega = cmath.exp(2j * cmath.pi / p)
    return [sum(omega ** (k * s) for s in S) for k in range(p)]


def trace_Ak(eigs, k):
    """Compute tr(A^k) = sum lambda_j^k."""
    return sum(e ** k for e in eigs).real


def main():
    print("=" * 70)
    print("P=19 CROSSOVER ANALYSIS: WHY INTERVAL BEATS PALEY")
    print("=" * 70)

    # ================================================================
    # SECTION 1: Spectral comparison across Paley primes
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 1: SPECTRAL STRUCTURE COMPARISON")
    print(f"{'=' * 60}")

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
        S_interval = frozenset(range(1, m + 1))

        eigs_p = eigenvalues(S_paley, p)
        eigs_i = eigenvalues(S_interval, p)

        mags_p = sorted([abs(e) for e in eigs_p[1:]], reverse=True)
        mags_i = sorted([abs(e) for e in eigs_i[1:]], reverse=True)

        print(f"\n  p={p}, m={m}:")
        print(f"    Paley eigenvalue magnitudes: "
              f"[{', '.join(f'{m:.3f}' for m in mags_p[:5])}...]")
        print(f"    Interval eigenvalue magnitudes: "
              f"[{', '.join(f'{m:.3f}' for m in mags_i[:5])}...]")
        print(f"    Paley: all |lambda| = {mags_p[0]:.4f} = sqrt({p})/2 = {math.sqrt(p)/2:.4f}")
        print(f"    Interval: |lambda_1| = {mags_i[0]:.4f}, "
              f"|lambda_2| = {mags_i[1]:.4f}")

        # Dirichlet kernel formula: |lambda_1| = sin(m*pi/p) / sin(pi/p)
        dk = abs(math.sin(m * math.pi / p) / math.sin(math.pi / p))
        print(f"    Dirichlet kernel |lambda_1| = {dk:.4f}")

        # Ratio of dominant eigenvalue
        ratio = mags_i[0] / mags_p[0]
        print(f"    Ratio |lambda_1(interval)| / |lambda(Paley)| = {ratio:.4f}")

    # ================================================================
    # SECTION 2: Trace / cycle count comparison
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 2: CYCLE COUNTS (TRACE-SIMPLE RANGE)")
    print(f"{'=' * 60}")

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
        S_interval = frozenset(range(1, m + 1))

        eigs_p = eigenvalues(S_paley, p)
        eigs_i = eigenvalues(S_interval, p)

        print(f"\n  p={p}:")
        print(f"  {'k':>4} {'tr(A^k)/k Paley':>20} {'tr(A^k)/k Interval':>20} {'Diff':>12}")

        for k in range(3, min(p + 1, 12), 2):
            tr_p = trace_Ak(eigs_p, k)
            tr_i = trace_Ak(eigs_i, k)
            ck_p = tr_p / k
            ck_i = tr_i / k
            diff = ck_p - ck_i
            # Note: for k >= 7, these are NOT exact cycle counts (non-simple walks)
            simple = "exact" if k <= 5 else "approx"
            print(f"  {k:>4} {ck_p:>20.1f} {ck_i:>20.1f} {diff:>+12.1f}  ({simple})")

        # c_3 should be constant
        c3_p = trace_Ak(eigs_p, 3) / 3
        c3_i = trace_Ak(eigs_i, 3) / 3
        print(f"    c_3: Paley={c3_p:.0f}, Interval={c3_i:.0f} "
              f"({'EQUAL' if abs(c3_p - c3_i) < 0.5 else 'DIFFERENT!'})")

        # c_5 comparison
        c5_p = trace_Ak(eigs_p, 5) / 5
        c5_i = trace_Ak(eigs_i, 5) / 5
        if c5_p > c5_i:
            print(f"    c_5: Paley WINS ({c5_p:.0f} > {c5_i:.0f})")
        else:
            print(f"    c_5: Interval WINS ({c5_i:.0f} > {c5_p:.0f})")

    # ================================================================
    # SECTION 3: Analytical margin formula
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 3: ANALYTICAL MARGIN ANALYSIS")
    print(f"{'=' * 60}")

    print("""
    For Paley T_p (p = 3 mod 4):
      All |lambda_k| = sqrt(p)/2 for k = 1,...,p-1

    For Interval C_p (S = {1,...,m}):
      lambda_k = sum_{s=1}^{m} omega^{ks} = (omega^k - omega^{k(m+1)}) / (1 - omega^k)
      |lambda_k| = |sin(pi*k*m/p) / sin(pi*k/p)|

    The TRACE COMPARISON at each power k:
      tr(A^k) = m^k + sum_{j=1}^{p-1} lambda_j^k

    For Paley: sum_{j=1}^{p-1} lambda_j^k = ?
    By Gauss sum theory: lambda_j = (chi(j)*g - 1)/2 where g = Gauss sum, |g|=sqrt(p)

    For k=3:
      tr(A^3)_Paley = m^3 + sum (chi(j)*g - 1)^3 / 8
      = m^3 + (1/8) * sum [chi(j)^3 * g^3 - 3*chi(j)^2 * g^2 + 3*chi(j)*g - 1]
      = m^3 + (1/8) * [g^3 * sum chi(j)^3 - 3*g^2 * sum chi(j)^2
                        + 3*g * sum chi(j) - (p-1)]

    Since chi is the Legendre symbol: chi(j)^2 = 1 for j != 0.
    sum_{j=1}^{p-1} chi(j) = 0
    sum_{j=1}^{p-1} chi(j)^2 = p-1
    sum_{j=1}^{p-1} chi(j)^3 = sum chi(j) = 0

    So: tr(A^3)_Paley = m^3 + (1/8) * [0 - 3*g^2*(p-1) + 0 - (p-1)]
                       = m^3 + (1/8) * (p-1) * [-3*g^2 - 1]

    With g^2 = chi(-1)*p = -p (for p = 3 mod 4):
      = m^3 + (1/8) * (p-1) * [3p - 1]
      = m^3 + (p-1)(3p-1)/8

    For Interval: tr(A^3)_Interval = same (c_3 is constant for regular circulants)

    So the difference comes from k >= 5. Let's compute tr(A^5):
    """)

    # Compute tr(A^5) analytically for both
    for p in [7, 11, 19, 23, 29, 31]:
        if p % 4 != 3:
            continue
        m = (p - 1) // 2
        S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
        S_interval = frozenset(range(1, m + 1))

        eigs_p = eigenvalues(S_paley, p)
        eigs_i = eigenvalues(S_interval, p)

        c5_p = trace_Ak(eigs_p, 5) / 5
        c5_i = trace_Ak(eigs_i, 5) / 5

        # Also compute higher traces
        c7_tr_p = trace_Ak(eigs_p, 7) / 7
        c7_tr_i = trace_Ak(eigs_i, 7) / 7

        c9_tr_p = trace_Ak(eigs_p, 9) / 9
        c9_tr_i = trace_Ak(eigs_i, 9) / 9

        # Who wins each power?
        wins_p = 0
        wins_i = 0
        for k in range(3, p + 1, 2):
            ck_p = trace_Ak(eigs_p, k) / k
            ck_i = trace_Ak(eigs_i, k) / k
            if ck_p > ck_i + 0.5:
                wins_p += 1
            elif ck_i > ck_p + 0.5:
                wins_i += 1

        print(f"  p={p}: c_5 Paley={c5_p:.0f}, Interval={c5_i:.0f}, "
              f"diff={c5_p-c5_i:+.0f}")
        print(f"         c_7/7 Paley={c7_tr_p:.0f}, Interval={c7_tr_i:.0f}, "
              f"diff={c7_tr_p-c7_tr_i:+.0f}")
        print(f"         Trace wins: Paley {wins_p}, Interval {wins_i} "
              f"(of {(p-1)//2} odd k values)")

    # ================================================================
    # SECTION 4: Eigenvalue phase structure
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 4: EIGENVALUE PHASE STRUCTURE")
    print(f"{'=' * 60}")

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
        S_interval = frozenset(range(1, m + 1))

        eigs_p = eigenvalues(S_paley, p)
        eigs_i = eigenvalues(S_interval, p)

        print(f"\n  p={p}: Paley eigenvalue phases (k=1,...,{m}):")
        for k in range(1, m + 1):
            e = eigs_p[k]
            mag = abs(e)
            phase = cmath.phase(e)
            print(f"    k={k}: |lambda|={mag:.4f}, arg={phase:.4f} ({phase/cmath.pi:.4f}*pi)")

        print(f"\n  p={p}: Interval eigenvalue phases (k=1,...,{m}):")
        for k in range(1, m + 1):
            e = eigs_i[k]
            mag = abs(e)
            phase = cmath.phase(e)
            print(f"    k={k}: |lambda|={mag:.4f}, arg={phase:.4f} ({phase/cmath.pi:.4f}*pi)")

    # ================================================================
    # SECTION 5: The crossover as competition of growth rates
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 5: CROSSOVER GROWTH RATE ANALYSIS")
    print(f"{'=' * 60}")

    print("""
    The H-advantage of Paley vs Interval depends on two competing effects:

    1. SHORT CYCLES (k=5): Paley's flat spectrum gives sum Re(lam^5) proportional
       to (p-1) * (sqrt(p)/2)^5 * (phase coherence factor).
       Interval's concentrated spectrum: ~ |lam_1|^5 * cos(5*arg(lam_1)).

    2. LONG CYCLES (k ~ p): For k near p, we're computing H itself (directed
       Hamiltonian cycles/paths). The interval's strong directional bias
       creates a "highway" effect.

    The crossover happens when the long-cycle advantage of the interval
    overcomes the short-cycle advantage of Paley.

    Growth rate comparison:
    - Paley tr(A^k): ~ (p-1) * (sqrt(p)/2)^k * O(1)  [phases may help or hurt]
    - Interval tr(A^k): ~ |lam_1|^k * cos(k*arg) ~ (p/pi)^k * O(1)

    Ratio: (p/pi)^k / (sqrt(p)/2)^k = (2/pi * sqrt(p))^k

    For p = 7: 2*sqrt(7)/pi ~ 1.68 > 1, so interval's dominant eigenvalue
    grows FASTER, but the (p-1)=6 other eigenvalues of Paley compensate.

    For large p: the SINGLE dominant eigenvalue (p/pi)^k overwhelms
    the (p-1) Paley terms (sqrt(p)/2)^k each.

    Crossover when: (p/pi)^k >> (p-1) * (sqrt(p)/2)^k
    i.e., (2*sqrt(p)/pi)^k >> p-1
    i.e., k > log(p) / log(2*sqrt(p)/pi)

    For k = p (the H computation), this becomes:
    (2*sqrt(p)/pi)^p >> p
    """)

    # Numerical crossover estimate
    for p in [7, 11, 19, 23, 29, 31, 37, 41, 43]:
        if p % 4 != 3:
            continue
        m = (p - 1) // 2
        # Ratio of dominant eigenvalue
        lam_interval = math.sin(m * math.pi / p) / math.sin(math.pi / p)
        lam_paley = math.sqrt(p) / 2
        ratio = lam_interval / lam_paley
        # At power k=p, single dominant term vs (p-1) equal terms
        # Dominant: lam_interval^p vs (p-1) * lam_paley^p
        log_ratio_p = p * math.log(ratio) - math.log(p - 1)
        winner = "INTERVAL" if log_ratio_p > 0 else "PALEY"
        print(f"  p={p}: |lam_1/lam_P| = {ratio:.4f}, "
              f"p*log(ratio) - log(p-1) = {log_ratio_p:.2f} -> {winner}")

    # ================================================================
    # SECTION 6: Geometric / tournament structure comparison
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 6: TOURNAMENT STRUCTURE COMPARISON")
    print(f"{'=' * 60}")

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        S_paley = set(j for j in range(1, p) if is_qr(j, p))
        S_interval = set(range(1, m + 1))

        # How many CONSECUTIVE pairs (s, s+1) are both in S?
        consec_p = sum(1 for s in range(1, p - 1) if s in S_paley and s + 1 in S_paley)
        consec_i = sum(1 for s in range(1, p - 1) if s in S_interval and s + 1 in S_interval)

        # How "clustered" is S?
        # Interval: S = {1,...,m} is maximally clustered
        # Paley: S is "random-looking"
        sorted_S = sorted(S_paley)
        gaps_p = [sorted_S[i + 1] - sorted_S[i] for i in range(len(sorted_S) - 1)]
        max_gap_p = max(gaps_p) if gaps_p else 0
        min_gap_p = min(gaps_p) if gaps_p else 0

        print(f"\n  p={p}:")
        print(f"    Paley QR = {sorted(S_paley)}")
        print(f"    Interval = {sorted(S_interval)}")
        print(f"    Consecutive pairs: Paley={consec_p}, Interval={m-1}")
        print(f"    Paley gaps: min={min_gap_p}, max={max_gap_p}, pattern={gaps_p}")
        print(f"    Interval: perfectly clustered (1 block of {m} consecutive elements)")


if __name__ == '__main__':
    main()
