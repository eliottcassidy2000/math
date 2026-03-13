#!/usr/bin/env python3
"""
resonance_cascade.py -- Resonance cascade in degree-2 Walsh of D^{2n}

DISCOVERY from true_deg2_walsh.py:
- D^{2n} degree-2 Walsh is nonzero only at "resonant" pairs
- D^4: 3-resonance (3a = +/-b mod p)
- D^6: 3-resonance + 5-resonance (5a = +/-b mod p)
- D^8: 3+5+7-resonance
- In general: D^{2n} picks up odd-resonances up to (2n-1)

This script:
1. Verifies the resonance cascade pattern
2. Analyzes the sign at each resonance level
3. Connects signs to Legendre characters chi(2k-1)
4. Tests whether the product law sign(h_hat)=chi(ab) holds at p=19

Author: kind-pasteur-2026-03-12-S60
"""

import math
from collections import defaultdict


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def compute_D2n_for_all_sigma(p, max_n):
    """Compute sum_{k=1}^{p-1} D_k^{2n} for all sigma in {+/-1}^m."""
    m = (p - 1) // 2
    results = {}
    for bits in range(1 << m):
        sigma = [(1 if bits & (1 << i) else -1) for i in range(m)]
        D_vals = {}
        for k in range(1, p):
            D_vals[k] = sum(
                sigma[l] * math.sin(2 * math.pi * k * (l+1) / p)
                for l in range(m)
            )
        power_sums = {}
        for n in range(1, max_n + 1):
            power_sums[n] = sum(D_vals[k] ** (2*n) for k in range(1, p))
        results[bits] = power_sums
    return results


def walsh_deg2(values, m, i, j):
    """True degree-2 Walsh coefficient at {i,j}."""
    total = 0
    for bits in range(1 << m):
        si = 1 if bits & (1 << i) else -1
        sj = 1 if bits & (1 << j) else -1
        total += values[bits] * si * sj
    return total / (1 << m)


def classify_resonance(a, b, p):
    """Classify the resonance type of pair (a,b).
    Returns list of (2k-1, type) where type is 'pos' or 'neg'
    meaning (2k-1)*a = b or (2k-1)*a = -b mod p.
    """
    resonances = []
    for k in range(1, p):
        q = 2*k - 1
        if q >= p:
            break
        if (q*a - b) % p == 0:
            resonances.append((q, f"{q}a=b"))
        if (q*a + b) % p == 0:
            resonances.append((q, f"{q}a=-b"))
        if (a - q*b) % p == 0 and q != 1:  # avoid double-counting q=1
            resonances.append((q, f"a={q}b"))
        if (a + q*b) % p == 0 and q != 1:
            resonances.append((q, f"a=-{q}b"))
    return resonances


def main():
    print("=" * 70)
    print("RESONANCE CASCADE IN DEGREE-2 WALSH")
    print("=" * 70)

    # PART 1: Classify all pairs by their resonance level
    print("\n--- PART 1: RESONANCE CLASSIFICATION ---")

    for p in [7, 11, 19, 23]:
        m = (p - 1) // 2
        print(f"\n  p={p}, m={m}:")

        for a in range(1, m+1):
            for b in range(a+1, m+1):
                res = classify_resonance(a, b, p)
                chi_ab = legendre(a*b, p)
                if res:
                    print(f"    ({a},{b}): chi={chi_ab:+d}, resonances: "
                          f"{', '.join(f'{q}:{t}' for q,t in res)}")
                else:
                    print(f"    ({a},{b}): chi={chi_ab:+d}, NO resonance (a*inv(b) not odd)")

    # PART 2: Verify cascade pattern
    print("\n--- PART 2: CASCADE PATTERN VERIFICATION ---")
    print("D^{2n} deg-2 nonzero iff pair has q-resonance with q <= 2n-1")

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        max_n = min(5, m)
        print(f"\n  p={p}:")

        all_D2n = compute_D2n_for_all_sigma(p, max_n)

        for a_idx in range(m):
            for b_idx in range(a_idx+1, m):
                a, b = a_idx+1, b_idx+1
                res = classify_resonance(a, b, p)
                min_q = min(q for q, t in res) if res else float('inf')

                onset = None
                for n in range(1, max_n+1):
                    fn = {bits: all_D2n[bits][n] for bits in range(1 << m)}
                    w = walsh_deg2(fn, m, a_idx, b_idx)
                    if abs(w) > 0.01 and onset is None:
                        onset = n

                expected_onset = None
                if res:
                    # First nonzero at D^{2n} where 2n-1 >= min_q, i.e., n >= (min_q+1)/2
                    expected_onset = (min_q + 1) // 2 + (1 if (min_q + 1) % 2 else 0)
                    # Actually: D^{2n} picks up q-resonance if q appears as a partial sum
                    # of the multinomial. For D^{2n}, the relevant q values are 1,3,5,...,2n-1
                    # So onset is at n = (min_q+1)/2
                    expected_onset = (min_q + 1) // 2

                match_str = ""
                if onset is not None and expected_onset is not None:
                    match_str = f" {'MATCH' if onset == expected_onset else 'MISMATCH'}"
                elif onset is None and not res:
                    match_str = " (never nonzero, no resonance)"

                print(f"    ({a},{b}): min_q={min_q if res else 'inf':>3}, "
                      f"onset_n={onset}, expected={expected_onset}{match_str}")

    # PART 3: Sign analysis at each resonance level
    print("\n--- PART 3: SIGN AT EACH RESONANCE LEVEL ---")
    print("For q-resonant pair (a,b) with qa=b: sign(D^{2n}) vs chi(ab)")
    print("chi(ab) = chi(qa^2) = chi(q)*chi(a)^2 = chi(q)")
    print("For qa=-b: chi(ab) = chi(-qa^2) = chi(-q)")

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        max_n = min(5, m)
        print(f"\n  p={p} ({p%4} mod 4):")
        print(f"    chi(3)={legendre(3,p):+d}, chi(5)={legendre(5,p):+d}, "
              f"chi(7)={legendre(7,p):+d}")
        print(f"    chi(-3)={legendre(-3,p):+d}, chi(-5)={legendre(-5,p):+d}, "
              f"chi(-7)={legendre(-7,p):+d}")

        all_D2n = compute_D2n_for_all_sigma(p, max_n)

        # Group pairs by their FIRST resonance
        by_q = defaultdict(list)
        for a_idx in range(m):
            for b_idx in range(a_idx+1, m):
                a, b = a_idx+1, b_idx+1
                res = classify_resonance(a, b, p)
                if res:
                    min_q = min(q for q, t in res)
                    first_res = [(q, t) for q, t in res if q == min_q]
                    by_q[min_q].append((a_idx, b_idx, a, b, first_res))

        for q in sorted(by_q):
            print(f"\n    q={q} resonance:")
            onset_n = (q + 1) // 2

            for a_idx, b_idx, a, b, first_res in by_q[q]:
                chi_ab = legendre(a*b, p)
                line = f"      ({a},{b}) [{first_res[0][1]}] chi={chi_ab:+d}: "

                for n in range(onset_n, max_n+1):
                    fn = {bits: all_D2n[bits][n] for bits in range(1 << m)}
                    w = walsh_deg2(fn, m, a_idx, b_idx)
                    sign_w = 1 if w > 0.01 else (-1 if w < -0.01 else 0)
                    ratio = w / (chi_ab * p) if chi_ab != 0 else 0
                    line += f" D^{2*n}:{sign_w:+d}(R={ratio:>8.2f})"

                print(line)

    # PART 4: The KEY question — is the sign at each q determined by chi(q)?
    print("\n--- PART 4: SIGN = f(chi(q)) HYPOTHESIS ---")
    print("Hypothesis: at q-resonant pair with qa=+/-b,")
    print("sign(D^{2n} deg-2) = chi(q) * chi(ab) * (-1)^{something}")

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        max_n = min(5, m)
        all_D2n = compute_D2n_for_all_sigma(p, max_n)

        print(f"\n  p={p}:")

        for a_idx in range(m):
            for b_idx in range(a_idx+1, m):
                a, b = a_idx+1, b_idx+1
                chi_ab = legendre(a*b, p)
                res = classify_resonance(a, b, p)
                if not res:
                    continue

                min_q = min(q for q, t in res)
                onset_n = (min_q + 1) // 2

                fn = {bits: all_D2n[bits][onset_n] for bits in range(1 << m)}
                w = walsh_deg2(fn, m, a_idx, b_idx)
                sign_w = 1 if w > 0.01 else (-1 if w < -0.01 else 0)

                # What determines the sign?
                # For qa=b type: chi(ab) = chi(qa^2) = chi(q)
                # For qa=-b type: chi(ab) = chi(-qa^2) = chi(-q) = chi(-1)*chi(q)
                # For a=qb type: chi(ab) = chi(qb^2) = chi(q)
                # For a=-qb type: chi(ab) = chi(-qb^2) = chi(-q) = chi(-1)*chi(q)

                first_type = res[0][1]
                if "=-" in first_type:
                    # Negative type: qa=-b or a=-qb
                    predicted_chi_from_q = legendre(-min_q, p)
                else:
                    predicted_chi_from_q = legendre(min_q, p)

                sign_times_chi = sign_w * chi_ab if sign_w != 0 else 0

                print(f"    ({a},{b}): q={min_q}, type={first_type}, "
                      f"chi(ab)={chi_ab:+d}, sign={sign_w:+d}, "
                      f"sign*chi={sign_times_chi:+d}, chi(q)={legendre(min_q,p):+d}, "
                      f"chi(-q)={legendre(-min_q,p):+d}")

    # PART 5: Value pattern — is |D^{2n} deg-2| = f(q) * p for q-onset pairs?
    print("\n--- PART 5: MAGNITUDE AT ONSET ---")
    print("At the onset power n=(q+1)/2, what is |D^{2n} deg-2|?")

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        max_n = min(5, m)
        all_D2n = compute_D2n_for_all_sigma(p, max_n)

        print(f"\n  p={p}:")

        for a_idx in range(m):
            for b_idx in range(a_idx+1, m):
                a, b = a_idx+1, b_idx+1
                res = classify_resonance(a, b, p)
                if not res:
                    continue

                min_q = min(q for q, t in res)
                onset_n = (min_q + 1) // 2

                if onset_n > max_n:
                    continue

                fn = {bits: all_D2n[bits][onset_n] for bits in range(1 << m)}
                w = walsh_deg2(fn, m, a_idx, b_idx)

                ratio_to_p = abs(w) / p
                print(f"    ({a},{b}): q={min_q}, onset D^{2*onset_n}, "
                      f"|W|={abs(w):.4f}, |W|/p={ratio_to_p:.4f}")

    # PART 6: Full h_hat at p=19 — does product law hold?
    print("\n--- PART 6: PRODUCT LAW AT p=19 ---")
    print("Computing H for all 512 orientations at p=19...")
    print("(This takes a few minutes)")

    p = 19
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]

    def count_ham_paths(A, n):
        dp = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                if dp[mask][v] == 0:
                    continue
                for w in range(n):
                    if mask & (1 << w):
                        continue
                    if A[v][w]:
                        dp[mask | (1 << w)][w] += dp[mask][v]
        full = (1 << n) - 1
        return sum(dp[full][v] for v in range(n))

    H_vals = {}
    import time
    t0 = time.time()
    for bits in range(1 << m):
        sigma = [(1 if bits & (1 << i) else -1) for i in range(m)]
        S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
        A = [[0]*p for _ in range(p)]
        for v in range(p):
            for s in S:
                A[v][(v + s) % p] = 1
        H_vals[bits] = count_ham_paths(A, p)

        if bits % 64 == 0 and bits > 0:
            elapsed = time.time() - t0
            pct = bits / (1 << m) * 100
            eta = elapsed / bits * ((1 << m) - bits)
            print(f"  Progress: {pct:.0f}% ({bits}/{1<<m}), ETA: {eta:.0f}s")

    t1 = time.time()
    print(f"  Computed {1 << m} orientations in {t1-t0:.1f}s")

    # Degree-2 Walsh of H
    print(f"\n  Degree-2 Walsh of H at p=19:")
    match_count = 0
    total = 0

    for i in range(m):
        for j in range(i+1, m):
            a, b = i+1, j+1
            h_hat = walsh_deg2(H_vals, m, i, j)
            chi_ab = legendre(a*b, p)
            sign_h = 1 if h_hat > 0 else (-1 if h_hat < 0 else 0)
            match = (sign_h == chi_ab)

            total += 1
            if match:
                match_count += 1

            res = classify_resonance(a, b, p)
            min_q = min(q for q, t in res) if res else 'inf'

            print(f"    ({a},{b}): h_hat={h_hat:>14.2f}, sign={sign_h:+d}, "
                  f"chi={chi_ab:+d}, match={match}, q={min_q}")

    print(f"\n  PRODUCT LAW: sign=chi in {match_count}/{total} pairs")
    print(f"  {'HOLDS' if match_count == total else 'FAILS'} at p=19")


if __name__ == '__main__':
    main()
