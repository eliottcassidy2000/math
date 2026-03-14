#!/usr/bin/env python3
"""
grand_theorem_implications_89.py -- opus-2026-03-14-S89

Implications of the Grand Fourier Level Theorem:
  E_{2k}/E_0 = 2(n-2k)^k / P(n,2k)

1. Product representation of 1 - Var/Mean^2?
2. Generating function structure
3. Connection to forbidden values
4. Asymptotic analysis
5. Connection to tribonacci
"""

from fractions import Fraction
from math import comb, factorial, log, exp
from functools import reduce


def var_ratio(n):
    """Compute exact Var/Mean^2 using the Grand Theorem."""
    K = (n - 1) // 2
    total = Fraction(0)
    for k in range(1, K + 1):
        pn2k = 1
        for i in range(2 * k):
            pn2k *= (n - i)
        term = Fraction(2 * (n - 2*k)**k, pn2k)
        total += term
    return total


def falling_factorial(n, k):
    """P(n,k) = n(n-1)...(n-k+1)."""
    result = 1
    for i in range(k):
        result *= (n - i)
    return result


def main():
    print("="*70)
    print("GRAND THEOREM IMPLICATIONS")
    print("opus-2026-03-14-S89")
    print("="*70)

    # Compute Var/Mean^2 for large range
    print(f"\n{'='*70}")
    print("PART 1: VAR/MEAN^2 SEQUENCE")
    print("="*70)

    print(f"\n  {'n':>3} {'Var/Mean^2':>20} {'decimal':>14} {'1/3 - ratio':>14}")
    for n in range(3, 25):
        vr = var_ratio(n)
        dev = Fraction(1, 3) - vr
        print(f"  {n:3d} {str(vr):>20} {float(vr):14.10f} {float(dev):14.10f}")

    # Is 1/3 - Var/Mean^2 = some nice function?
    print(f"\n{'='*70}")
    print("PART 2: DEVIATION FROM 1/3")
    print("="*70)

    print(f"\n  1/3 - Var/Mean^2 = sum of level correction terms")
    print(f"\n  Actually: 1/3 - E2/E0 = 1/3 - 2(n-2)/(n(n-1))")
    print(f"    = [n(n-1) - 6(n-2)] / (3n(n-1))")
    print(f"    = [n^2-n-6n+12] / (3n(n-1))")
    print(f"    = [n^2-7n+12] / (3n(n-1))")
    print(f"    = [(n-3)(n-4)] / (3n(n-1))")

    print(f"\n  So: 1/3 - E2/E0 = (n-3)(n-4) / (3n(n-1))")
    for n in range(3, 10):
        m = n * (n - 1) // 2
        actual = Fraction(1, 3) - Fraction(n-2, m)
        predicted = Fraction((n-3)*(n-4), 3*n*(n-1))
        print(f"    n={n}: (n-3)(n-4)/(3n(n-1)) = {predicted}, actual = {actual}, match = {predicted == actual}")

    print(f"\n  For n=3,4: this is 0 (Var/Mean^2 = 1/3 exactly)")
    print(f"  For n>=5: level-4+ corrections bring it closer to 1/3 (but don't reach it)")

    # What fraction of deviation comes from E4?
    print(f"\n  Deviation budget:")
    for n in range(5, 10):
        e2 = Fraction(2*(n-2), n*(n-1))
        K = (n - 1) // 2
        e4_plus = Fraction(0)
        for k in range(2, K + 1):
            e4_plus += Fraction(2 * (n-2*k)**k, falling_factorial(n, 2*k))

        print(f"    n={n}: E2={e2}, E4+={e4_plus}, total={e2+e4_plus}, 1/3-total={Fraction(1,3)-e2-e4_plus}")

    # Product representation
    print(f"\n{'='*70}")
    print("PART 3: IS THERE A PRODUCT FORMULA?")
    print("="*70)

    # Var/Mean^2 = sum 2(n-2k)^k / P(n,2k)
    # Can this be written as a product?

    # Let's compute 1 - 3*Var/Mean^2 (which is 0 for n=3,4)
    print(f"\n  1 - 3*Var/Mean^2:")
    for n in range(3, 15):
        vr = var_ratio(n)
        x = 1 - 3 * vr
        print(f"    n={n}: {x} = {float(x):.10f}")

    # Factor these as (n-3)(n-4)/something
    print(f"\n  (n-3)(n-4) / [3n(n-1) * (1 - 3*Var)] = ?")
    for n in range(5, 12):
        vr = var_ratio(n)
        x = 1 - 3 * vr
        if x != 0:
            y = Fraction((n-3)*(n-4), 3*n*(n-1)) / x
            print(f"    n={n}: ratio = {y} = {float(y):.6f}")

    print(f"\n  1 - 3*Var/Mean^2 sequence:")
    devs = {}
    for n in range(3, 15):
        vr = var_ratio(n)
        x = 1 - 3 * vr
        devs[n] = x
        print(f"    n={n}: {x}")

    # Try: is (1-3V) = (n-3)(n-4)*f(n)?
    print(f"\n  (1-3V) / [(n-3)(n-4)]:")
    for n in range(5, 15):
        x = devs[n]
        y = x / Fraction((n-3)*(n-4), 1)
        print(f"    n={n}: {y} = {float(y):.10f}")

    # The MEAN is n!/2^{n-1}. The denominators of Var/Mean^2:
    print(f"\n  Denominators of Var/Mean^2:")
    for n in range(3, 15):
        vr = var_ratio(n)
        print(f"    n={n}: {vr.denominator}")

    # Connection to tribonacci
    print(f"\n{'='*70}")
    print("PART 4: TRIBONACCI CONNECTION")
    print("="*70)

    trib = [0, 0, 1]
    for _ in range(30):
        trib.append(trib[-1] + trib[-2] + trib[-3])

    print(f"\n  Tribonacci numbers: {[t for t in trib if 0 < t < 10000]}")

    for n in range(3, 15):
        vr = var_ratio(n)
        d = vr.denominator
        in_trib = d in trib
        print(f"    n={n}: denom = {d}, tribonacci? {in_trib}")

    # 504 = T(14) is tribonacci! (n=7 denominator)

    # Connection to n=7 value 131/504
    print(f"\n  The number 131:")
    print(f"    131 = 2^7 + 3 = 128 + 3")
    print(f"    131 is the 32nd prime")
    print(f"    131 in the Var/Mean^2 numerator at n=7")

    # Asymptotic behavior
    print(f"\n{'='*70}")
    print("PART 5: ASYMPTOTIC ANALYSIS")
    print("="*70)

    print(f"\n  As n -> inf:")
    print(f"    E2/E0 ~ 2/n -> 0")
    print(f"    E4/E0 ~ 2n^2/(n^4) = 2/n^2 -> 0 (faster)")
    print(f"    E_{{2k}}/E0 ~ 2n^k/n^{{2k}} = 2/n^k -> 0")
    print(f"  ")
    print(f"    Var/Mean^2 ~ sum_k 2/n^k = 2n/(n-1) * 1/n = 2/(n-1)")
    print(f"    More precisely: Var/Mean^2 = 2/n + 2/n^2 + 2/n^3 + ... = 2/(n-1)")
    print(f"    But this is a HEURISTIC. The actual terms have (n-2k)^k != n^k")

    print(f"\n  Exact vs 2/(n-1) approximation:")
    for n in range(3, 25):
        vr = var_ratio(n)
        approx = Fraction(2, n-1) if n > 1 else 0
        diff = float(vr - approx)
        print(f"    n={n:2d}: exact={float(vr):.8f}, 2/(n-1)={float(approx):.8f}, diff={diff:.8f}")

    # Better approximation: 2/(n-1) + correction
    print(f"\n  Better approximation: Var/Mean^2 vs 2(n-2)/(n(n-1)):")
    for n in range(3, 12):
        vr = var_ratio(n)
        e2 = Fraction(2*(n-2), n*(n-1))
        ratio = vr / e2 if e2 > 0 else 0
        print(f"    n={n}: ratio Var/E2 = {ratio} = {float(ratio):.8f}")

    # For the forbidden value connection
    print(f"\n{'='*70}")
    print("PART 6: FORBIDDEN VALUES AND FOURIER STRUCTURE")
    print("="*70)

    # The Grand Theorem tells us about the VARIANCE of H, not individual values.
    # But the forbidden values must satisfy Fourier constraints.
    # H(T) = sum_S c_hat_S chi_S(T) where chi_S(T) = (-1)^{S.T}
    # H is always odd (proved). H ranges from 1 to max_H.

    # For H = 63 = 111111 in binary:
    # 63 = Mean(7) - 15.75 = 78.75 - 15.75
    # Deviation from mean: 78.75 - 63 = 15.75 = 63/4
    # In units of sigma: |63 - 78.75|/sqrt(Var) ~ 15.75/40.15 ~ 0.39
    # So 63 is only 0.39 sigma below the mean -- NOT extreme!

    mean_7 = Fraction(315, 4)
    var_7 = Fraction(206325, 128)
    sd_7 = float(var_7) ** 0.5
    z_63 = float(mean_7 - 63) / sd_7
    print(f"\n  n=7: Mean = {mean_7} = {float(mean_7)}")
    print(f"  Var = {var_7} = {float(var_7):.2f}")
    print(f"  SD = {sd_7:.4f}")
    print(f"  Z-score of H=63: {z_63:.4f}")
    print(f"  Z-score of H=7: {float(mean_7 - 7)/sd_7:.4f}")
    print(f"  Z-score of H=21: {float(mean_7 - 21)/sd_7:.4f}")
    print(f"  Z-score of H=189 (max): {float(189 - mean_7)/sd_7:.4f}")

    print(f"\n  The forbidden values are NOT extreme! They're near the bulk.")
    print(f"  The obstruction must be STRUCTURAL, not statistical.")

    # Level-2 analysis for H=63
    # Mean = 315/4 = 78.75
    # Level 2: each c_hat is +/-15/4 = +/-3.75, with 105 nonzero
    # Level 4: c_hat is +/-3/16 (1260 of these) or +/-3/8 (630 of these)
    # Level 6: c_hat is +/-1/32, with 2520

    # H(T) = 78.75 + sum(level 2) + sum(level 4) + sum(level 6)
    # For H = 63: need level2 + level4 + level6 = -15.75
    # Each level-2 term is +/-3.75. If p are + and q=105-p are -:
    #   sum = (2p-105)*3.75
    #   Need (2p-105)*3.75 ~ -15.75 -> 2p-105 ~ -4.2 -> p ~ 50.4
    # So need p = 50 or 51 level-2 terms positive. Perfectly feasible!

    print(f"\n  Level-2 analysis for H=63:")
    print(f"  Need level-2 sum = -(78.75 - 63) - (level 4+6 corrections) ~ -15.75")
    print(f"  Level-2 terms are +/-{float(Fraction(15,4))}")
    print(f"  Need ~50.4 out of 105 terms positive -> perfectly feasible")
    print(f"  The obstruction is NOT a magnitude constraint but a PARITY/LATTICE constraint")

    # H must be an odd integer. What lattice does the Fourier sum live on?
    # Total H = 315/4 + 15/4 * S2 + 3/16 * S4a + 3/8 * S4b + 1/32 * S6
    # where S2 in {-105, -103, ..., 103, 105}
    # S4a in {-1260, ..., 1260}, S4b in {-630, ..., 630}
    # S6 in {-2520, ..., 2520}

    # Multiply everything by 32 to clear denominators:
    # 32H = 2520 + 120*S2 + 6*S4a + 12*S4b + S6

    print(f"\n  Lattice analysis (multiply by 32):")
    print(f"  32H = 2520 + 120*S2 + 6*S4a + 12*S4b + S6")
    print(f"  where S2, S4a, S4b, S6 are sums of +/-1 variables")
    print(f"  32*63 = {32*63} = 2016")
    print(f"  32*7 = {32*7} = 224")
    print(f"  32*21 = {32*21} = 672")

    # For H=63: 32*63 = 2016 = 2520 + 120*S2 + 6*S4a + 12*S4b + S6
    # -> 120*S2 + 6*S4a + 12*S4b + S6 = -504

    # Modular analysis:
    # mod 8: 120*S2 = 0 (120=8*15), 6*S4a, 4*S4b (12 mod 8=4), S6
    # 6*S4a + 4*S4b + S6 = -504 mod 8 = 0 mod 8
    print(f"\n  Mod 8 constraint: 6*S4a + 4*S4b + S6 = 0 mod 8")
    print(f"  (same for all H values)")

    # mod 3: everything divisible by 3 except S6
    # S6 = -504 mod 3 = 0 mod 3
    print(f"  Mod 3: S6 = 0 mod 3 for all H")

    # mod 7: 120=1, 6=6, 12=5, 1=1 (mod 7)
    # S2 + 6*S4a + 5*S4b + S6 = -504 mod 7 = 0 mod 7
    print(f"  Mod 7: S2 + 6*S4a + 5*S4b + S6 = 0 mod 7")
    print(f"  This constrains which H values are reachable mod 7!")

    # PART 7: GCD structure of H-spectrum
    print(f"\n{'='*70}")
    print("PART 7: GCD AND LATTICE STRUCTURE OF H-SPECTRUM")
    print("="*70)

    # At n=7, GCD of all H values is 16 (from the deep spectrum computation).
    # So H = 16*h for some odd integer h.
    # Equivalently, 32H = 32*16*h = 512*h.
    # In our lattice equation: 2520 + 120*S2 + 6*S4a + 12*S4b + S6 = 512*h
    # Is 512 | (2520 + 120*S2 + 6*S4a + 12*S4b + S6)?

    # GCD constraint mod 512:
    # 2520 mod 512 = 2520 - 4*512 = 2520 - 2048 = 472
    # So we need 120*S2 + 6*S4a + 12*S4b + S6 = -472 mod 512

    print(f"\n  At n=7: GCD of H-spectrum = 16")
    print(f"  So 32H = 512h for odd h, meaning all reachable lattice points are 0 mod 512")
    print(f"  2520 mod 512 = {2520 % 512}")
    print(f"  Need: 120*S2 + 6*S4a + 12*S4b + S6 = {(-2520) % 512} mod 512")

    # Actually let's verify: for H divisible by 16,
    # 32*16 = 512, so 32H mod 512 = 0 for all achievable H.
    # But wait -- H is NOT always divisible by 16. H can be 1 at n=7?
    # No, from the spectrum data, minimum H at n=7 is 1... wait.
    # Let me re-check. The spectrum has values 1, 3, 5, ... (odd).
    # GCD = 16 would mean all values are multiples of 16. But 1 is in the spectrum?
    # No -- at n=7, the spectrum was: {1, 3, 5, ...} but divided by... let me recheck.

    # From the deep spectrum computation, H values at n=7 range from 1 to 189.
    # The GCD was reported as 16 for some specific case, or maybe that was wrong.
    # Let me just verify: minimum H at n=7.

    # Actually: for n=7, minimum H is 1 (transitive tournament has H=1? No!)
    # Transitive tournament T_7 has MANY paths. Let me compute.

    # The transitive tournament on {0,...,6} has i->j iff i<j.
    # H(transitive T_n) = 1 (only one Hamiltonian path: 0->1->2->...->n-1).
    # Wait, that IS 1 for n=7.
    # But from our spectrum analysis, GCD was... let me reconsider.
    # From n7_deep_spectrum: H_min=1, H_max=189.
    # GCD of {1, 3, 5, ...} = 1 since 1 is in the set.
    # I previously reported GCD=16 which must have been an error, or for a different quantity.

    # Let's analyze the actual H values more carefully.
    # All H are odd (from complement symmetry: H(T) + H(T-bar) = 2*H(T) -> H always integer;
    # and the oddness follows from Redei's theorem).

    print(f"\n  Re-checking: H(transitive T_7) = 1 (only path: 0->1->...->6)")
    print(f"  So GCD of H-spectrum at n=7 includes 1, hence GCD = 1")
    print(f"  The previous GCD=16 report may have been for a DIFFERENT quantity")
    print(f"  (perhaps the level-2 energy sublattice or coefficient denominators)")

    # What IS structured: the H values mod small primes.
    # Since all H are odd, H mod 2 = 1 always.
    # What about H mod 4? H mod 8?

    # PART 8: Deeper asymptotic -- generating function for Var/Mean^2
    print(f"\n{'='*70}")
    print("PART 8: GENERATING FUNCTION FOR VAR/MEAN^2")
    print("="*70)

    # V(n) = sum_{k=1}^{floor((n-1)/2)} 2(n-2k)^k / P(n,2k)
    # Let's look at V(n) * n! / 2:
    print(f"\n  V(n) * n! / 2:")
    for n in range(3, 15):
        vr = var_ratio(n)
        scaled = vr * factorial(n) / 2
        print(f"    n={n}: V*n!/2 = {scaled}")

    # Or V(n) * P(n,2) = V(n) * n*(n-1):
    print(f"\n  V(n) * n*(n-1):")
    for n in range(3, 15):
        vr = var_ratio(n)
        scaled = vr * n * (n-1)
        print(f"    n={n}: V*n*(n-1) = {scaled} = {float(scaled):.6f}")

    # V * n*(n-1) = 2(n-2) + sum_{k>=2} 2(n-2k)^k * n*(n-1) / P(n,2k)
    # The first term is 2(n-2), and the rest are correction terms.

    # More interesting: partial sums up to level 2k
    print(f"\n  Partial sums V_K(n) = sum_{{k=1}}^K E_{{2k}}/E0:")
    for n in [5, 7, 9, 11]:
        K = (n - 1) // 2
        print(f"  n={n} (K_max={K}):")
        running = Fraction(0)
        for k in range(1, K + 1):
            term = Fraction(2 * (n-2*k)**k, falling_factorial(n, 2*k))
            running += term
            print(f"    k={k}: E_{{2{k}}}/E0 = {term}, cumulative = {running} = {float(running):.8f}")

    # PART 9: The "almost 1/3" phenomenon
    print(f"\n{'='*70}")
    print("PART 9: THE 'ALMOST 1/3' PHENOMENON")
    print("="*70)

    # For n=3,4: Var/Mean^2 = 1/3 EXACTLY
    # For n>=5: Var/Mean^2 < 1/3
    # WHY is E2/E0 = (n-2)/m = 2(n-2)/(n(n-1)) so close to 1/3?
    # Because (n-2)/m -> 2/n -> 0 as n->inf.
    # Actually 1/3 - (n-2)/m = (n-3)(n-4)/(3n(n-1)), which -> 1/3 as n->inf!
    # So E2/E0 is NOT close to 1/3 for large n. It's close only for small n.

    # The interesting thing is:
    # V(3) = V(4) = 1/3 (exact)
    # V(5) = 19/60 ~ 0.317 (close to 1/3)
    # V(6) = 13/45 ~ 0.289
    # V(7) = 131/504 ~ 0.260
    # It's DECREASING towards 0 (like 2/n).

    # The fact that V(3)=V(4)=1/3 is remarkable. WHY?
    # V(3) = E2/E0 = 2*1/(3*2) = 1/3. Only one level.
    # V(4) = E2/E0 = 2*2/(4*3) = 1/3. Also only one level (E4 = 0 since (n-4)^2 = 0).
    # So V(4) = 1/3 because (n-2)/m = 2/6 = 1/3, AND there are no higher levels.

    # For n=4: floor((4-1)/2) = 1, so only k=1 contributes. Good.
    # For n=5: k=1,2 contribute. E4/E0 = 2*1/(5*4*3*2) = 2/120 = 1/60.
    # V(5) = (n-2)/m + E4/E0 = 3/10 + 1/60 = 18/60 + 1/60 = 19/60.

    print(f"\n  Why V(3) = V(4) = 1/3:")
    print(f"    n=3: E2/E0 = 2*1/(3*2) = 1/3, only level")
    print(f"    n=4: E2/E0 = 2*2/(4*3) = 1/3, E4/E0 = 2*0^2/P(4,4) = 0")
    print(f"    n=5: E2/E0 = 3/10, E4/E0 = 1/60")
    print(f"    n=6: E2/E0 = 4/15, E4/E0 = 2*4/(6*5*4*3) = 8/360 = 1/45")

    # Is there a deeper reason? The n=4 case: E4 vanishes because (n-4)^2 = 0.
    # This is the SAME reason that E_{n-1} vanishes for even n: (n-2k)^k = 0 when 2k=n.

    # PART 10: Telescope / cancellation structure
    print(f"\n{'='*70}")
    print("PART 10: TELESCOPING AND PARTIAL FRACTION DECOMPOSITION")
    print("="*70)

    # V(n) = sum 2(n-2k)^k / P(n,2k)
    # Can we decompose each term as a partial fraction in n?

    # For k=1: 2(n-2)/P(n,2) = 2(n-2)/(n(n-1)) = 2/(n-1) - 2/(n(n-1))
    #   = 2/(n-1) - 2[1/(n-1) - 1/n] = 2/(n-1) - 2/(n-1) + 2/n = 2/n
    # Wait: 2(n-2)/(n(n-1)) = 2[n/(n(n-1)) - 2/(n(n-1))] = 2/(n-1) - 4/(n(n-1))
    # Hmm, let me just do it: (n-2)/(n(n-1)) by partial fractions.
    # A/n + B/(n-1) = (A(n-1) + Bn)/(n(n-1)) = ((A+B)n - A)/(n(n-1))
    # Need A+B = 1, -A = -2 -> A = 2, B = -1
    # So (n-2)/(n(n-1)) = 2/n - 1/(n-1)
    # Therefore E2/E0 = 2*(2/n - 1/(n-1)) = 4/n - 2/(n-1)
    # Check: 4/3 - 2/2 = 4/3 - 1 = 1/3. For n=3. YES!
    # 4/4 - 2/3 = 1 - 2/3 = 1/3. For n=4. YES!
    # 4/5 - 2/4 = 0.8 - 0.5 = 0.3 = 3/10. For n=5. YES!

    print(f"\n  Partial fraction for E2/E0:")
    print(f"  E2/E0 = 2(n-2)/(n(n-1)) = 4/n - 2/(n-1)")
    for n in range(3, 10):
        e2 = Fraction(2*(n-2), n*(n-1))
        pf = Fraction(4, n) - Fraction(2, n-1)
        print(f"    n={n}: E2={e2}, 4/n - 2/(n-1) = {pf}, match={e2 == pf}")

    # For k=2: 2(n-4)^2 / P(n,4) = 2(n-4)^2 / (n(n-1)(n-2)(n-3))
    # This is harder. But the pattern of the first term is suggestive.
    # 4/n - 2/(n-1) is a "telescoping" decomposition.

    # Let's see if the TOTAL V(n) has a nice partial fraction form.
    print(f"\n  V(n) as rational function of n (exact values):")
    print(f"  Looking for V(n) = A/n + B/(n-1) + C/(n-2) + ...")
    # V(n) is not a rational function of n because K = floor((n-1)/2) depends on n.
    # But for fixed parity of n, it might be.

    # For ODD n (K = (n-1)/2): these include the highest level.
    # For EVEN n (K = (n-2)/2 = n/2 - 1): one fewer level.

    # Let's focus on even n:
    print(f"\n  Even n sequence:")
    for n in range(4, 20, 2):
        vr = var_ratio(n)
        print(f"    n={n}: V = {vr} = {float(vr):.10f}")

    print(f"\n  Odd n sequence:")
    for n in range(3, 20, 2):
        vr = var_ratio(n)
        print(f"    n={n}: V = {vr} = {float(vr):.10f}")

    # PART 11: Ratio of consecutive V(n)
    print(f"\n{'='*70}")
    print("PART 11: RATIO ANALYSIS")
    print("="*70)

    print(f"\n  V(n)/V(n-1) ratio:")
    prev = None
    for n in range(3, 20):
        vr = var_ratio(n)
        if prev is not None:
            ratio = vr / prev
            print(f"    V({n})/V({n-1}) = {ratio} = {float(ratio):.8f}")
        prev = vr

    # V(n+1)/V(n) should approach 1 - 2/n or similar for large n.
    print(f"\n  1 - V(n+1)/V(n) (fractional decrease per step):")
    prev = None
    for n in range(3, 20):
        vr = var_ratio(n)
        if prev is not None:
            diff = 1 - vr / prev
            pred = Fraction(1, n-1)  # Guess: decrease ~ 1/(n-1)?
            print(f"    n={n}: 1-V(n)/V(n-1) = {float(diff):.8f}, 1/(n-1) = {float(pred):.8f}")
        prev = vr

    print(f"\n{'='*70}")
    print("DONE -- GRAND THEOREM IMPLICATIONS")
    print("="*70)


if __name__ == "__main__":
    main()
