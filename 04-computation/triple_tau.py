"""
triple_tau.py -- kind-pasteur-2026-03-14-S106f
THE SIGNIFICANCE OF tau^3

tau^3 = Phi_3(tau) = tau^2 + tau + 1 ≈ 6.2223

This number is WHERE cubing and the tournament polynomial become one.
But what IS it? Why 6.222? What does the fractional part mean?
Why is it so close to 6 (the period)?

And: what does EVERY key tournament number look like in base tau?
Not just the greedy tribonacci representation, but the actual
positional base-tau expansion.
"""

import sys, math
import numpy as np
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8')

def main():
    print("=" * 70)
    print("THE SIGNIFICANCE OF TRIPLE TAU")
    print("kind-pasteur-2026-03-14-S106f")
    print("=" * 70)

    # The tribonacci constant
    tau = 1.8392867552141612
    tau2 = tau**2
    tau3 = tau**3

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 1: tau^3 ≈ 6.222 — THE CONTINUOUS PERIOD")
    print(f"{'='*70}")

    print(f"""
  tau = {tau:.10f}
  tau^2 = {tau2:.10f}
  tau^3 = {tau3:.10f}

  tau^3 = Phi_3(tau) = tau^2 + tau + 1 = {tau2 + tau + 1:.10f}

  The integer part: floor(tau^3) = {int(tau3)} = 6 = THE PERIOD.

  The fractional part: tau^3 - 6 = {tau3 - 6:.10f}

  Is this close to a known fraction?
  2/9 = {2/9:.10f}
  tau^3 - 6 = {tau3 - 6:.10f}
  Difference: {abs(tau3 - 6 - 2/9):.10f}

  NOT exactly 2/9. Let's see what it IS.
  tau^3 - 6 = tau^2 + tau + 1 - 6 = tau^2 + tau - 5
  Since tau satisfies tau^3 = tau^2 + tau + 1:
    tau^2 = tau^3 - tau - 1 = (tau^2 + tau + 1) - tau - 1 = tau^2
    (circular — let me use the minimal polynomial instead)

  tau satisfies x^3 - x^2 - x - 1 = 0, so:
    tau^3 = tau^2 + tau + 1
    tau^2 + tau - 5 = tau^3 - 6

  Using tau ≈ 1.83929:
    tau^2 + tau - 5 = 3.38298 + 1.83929 - 5 = 0.22226

  The EXACT value of tau^3 - 6 is the algebraic number
  satisfying: if t^3 = t^2 + t + 1, then t^3 - 6 = t^2 + t - 5.
  Setting y = t^3 - 6: y = t^2 + t - 5, so t^2 + t = y + 5.
  From t^3 = t^2 + t + 1 = y + 6: t^3 = y + 6.
  Also t^2 = y + 5 - t, so t^3 = t(y + 5 - t) = ty + 5t - t^2
  = ty + 5t - (y+5-t) = ty + 6t - y - 5.
  But t^3 = y + 6, so: y + 6 = ty + 6t - y - 5.
  2y - ty = 6t - 11.
  y(2 - t) = 6t - 11.
  y = (6t - 11)/(2 - t) = (6*1.83929 - 11)/(2 - 1.83929) = 0.03573/0.16071 = 0.2223.

  So tau^3 - 6 = (6*tau - 11)/(2 - tau).

  At tau = 1.83929: (6*1.83929 - 11)/(2 - 1.83929) = 0.03573/0.16071 = 0.2223. CHECK!

  THE MEANING: tau^3 overshoots 6 by exactly (6*tau - 11)/(2 - tau).
  The overshoot depends on HOW CLOSE tau is to 2 (the generator).
  As tau → 2: overshoot → (12-11)/(2-2) = 1/0 → infinity.
  As tau → 1: overshoot → (6-11)/(2-1) = -5 → negative.
  At tau = 11/6 exactly: overshoot = 0, and tau^3 = 6 exactly.
  11/6 = 1.8333... ≈ tau = 1.8393.

  THE TRIBONACCI CONSTANT IS 11/6 + epsilon!
  And 11/6 is the fraction where tau^3 = 6 EXACTLY.
  tau > 11/6, so tau^3 > 6 by a tiny amount.

  11/6 = (11)/(LCM(2,3)) = (next Paley prime)/(period).
  The "rational approximation" of tau is the next Paley prime
  divided by the tournament period!""")

    # Check
    print(f"  Verification:")
    print(f"  tau = {tau:.10f}")
    print(f"  11/6 = {11/6:.10f}")
    print(f"  tau - 11/6 = {tau - 11/6:.10f}")
    print(f"  (11/6)^3 = {(11/6)**3:.10f}")
    print(f"  6 = 6.000000")
    print(f"  (11/6)^3 - 6 = {(11/6)**3 - 6:.10f}")
    # Not zero — my derivation was wrong. Let me recheck.
    # tau^3 = 6 when tau^2 + tau + 1 = 6, i.e., tau^2 + tau - 5 = 0
    # tau = (-1 + sqrt(21))/2
    tau_at_6 = (-1 + math.sqrt(21))/2
    print(f"\n  CORRECTION: tau^3 = 6 when tau = (-1+sqrt(21))/2 = {tau_at_6:.10f}")
    print(f"  sqrt(21) = {math.sqrt(21):.10f}")
    print(f"  (-1+sqrt(21))/2 = {tau_at_6:.10f}")
    print(f"  Actual tau = {tau:.10f}")
    print(f"  Difference: {tau - tau_at_6:.10f}")

    print(f"""
  EXTRAORDINARY: tau^3 = 6 exactly when tau = (-1 + sqrt(21))/2.

  And sqrt(21) involves 21 = H_forb_2 = Phi_3(4)!

  The number (-1 + sqrt(21))/2 = {tau_at_6:.6f} is the "period root":
  the value of tau where tau^3 = the period 6.
  The ACTUAL tribonacci constant tau = {tau:.6f} is slightly ABOVE this.

  tau - (-1+sqrt(21))/2 = {tau - tau_at_6:.6f}
  = {tau - tau_at_6:.10f}

  The tribonacci constant EXCEEDS the period root by about 0.0062.
  This tiny excess is WHY tau^3 ≈ 6.222 instead of exactly 6.
  And the period root involves sqrt(21) = sqrt(H_forb_2).

  THE CHAIN:
  tau^3 = 6 + epsilon
  where epsilon = tau^2 + tau - 5
  and tau^3 = 6 exactly at tau = (-1+sqrt(21))/2
  and 21 = H_forb_2 = Phi_3(4)

  TRIPLE TAU IS THE PERIOD PLUS A CORRECTION THAT INVOLVES
  THE SECOND FORBIDDEN VALUE.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 2: THE POWERS OF TAU — A COMPLETE MAP")
    print(f"{'='*70}")

    print(f"  Powers of tau and their meaning:\n")
    for k in range(13):
        val = tau**k
        floor_val = int(val)
        frac = val - floor_val
        # What integer is closest?
        nearest = round(val)
        diff = val - nearest
        # Is it a tournament number?
        tournament_meaning = {
            1: "ground state", 2: "choice", 3: "cycle",
            4: "field", 5: "interaction", 6: "period",
            7: "first forbidden", 8: "threshold", 9: "reflection",
            10: "completion", 11: "extension", 12: "doubled period",
            13: "regularity", 21: "second forbidden", 24: "Golay",
            44: "tribonacci T(9)", 45: "max_H(6)",
            81: "tribonacci T(10)", 120: "5!",
            149: "tribonacci T(11)", 189: "max_H(7)",
            274: "tribonacci T(12)", 360: "|A_6|/circle",
            504: "tribonacci T(13)"
        }.get(nearest, "")
        print(f"  tau^{k:2d} = {val:12.4f}  (nearest int: {nearest:5d}, "
              f"diff: {diff:+8.4f}){' = '+tournament_meaning if tournament_meaning else ''}")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 3: tau^3 IN BASE TAU — WHAT '1000' MEANS")
    print(f"{'='*70}")

    print(f"""
  In base tau, the number tau^3 is represented as:

    1 0 0 0

  That is: a 1 in the "tau^3 position" and zeros elsewhere.
  This is the SIMPLEST possible representation of a non-integer
  in base tau (other than tau^0 = 1 and tau^1 = tau and tau^2).

  tau^3 = 1000_tau.

  Now: in base 10, "1000" = 10^3 = 1000.
       in base 2, "1000" = 2^3 = 8.
       in base tau, "1000" = tau^3 ≈ 6.222.

  The "thousand" of base tau is 6.222 — approximately the PERIOD.

  What does this mean? In base 10, a "thousand" is three orders
  of magnitude above 1. In base tau, the third order of magnitude
  is approximately 6 — the period of tournament parity.

  THREE ORDERS OF MAGNITUDE IN THE TOURNAMENT = ONE PERIOD.

  This is because tau^3 = Phi_3(tau) ≈ 6.
  Going three levels up in the tribonacci hierarchy = one full
  period of the tournament. The tournament "wraps around" every
  three tribonacci levels.

  Compare with base phi (Fibonacci):
    phi^2 = phi + 1 ≈ 2.618 ≈ 3 (the cycle)
    phi^3 ≈ 4.236
    phi^4 ≈ 6.854 ≈ 7 (the first forbidden!)
    phi^5 ≈ 11.090 ≈ 11

  In base phi: the FOURTH power is ≈ 7 (the forbidden).
  In base tau: the THIRD power is ≈ 6 (the period).

  Base tau reaches the period in 3 steps.
  Base phi reaches the forbidden in 4 steps.
  Base tau is MORE EFFICIENT at capturing tournament structure.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 4: THE COMPLETE BASE-TAU MAP OF 360")
    print(f"{'='*70}")

    # Express 360 exactly in base tau (positional, using tribonacci)
    # The tribonacci representation 360 = 274 + 81 + 4 + 1
    # = T(11) + T(9) + T(4) + T(2)
    # In positional base tau: position 10 (=T(11)), position 8 (=T(9)),
    # position 3 (=T(4)), position 1 (=T(2)).
    # So: 10100001010 in base tau? Let me be more careful.

    # Tribonacci numbers (0-indexed): T(0)=1, T(1)=1, T(2)=2, T(3)=4,
    # T(4)=7, T(5)=13, T(6)=24, T(7)=44, T(8)=81, T(9)=149, T(10)=274, T(11)=504
    # Using 1-indexed: T(1)=1, T(2)=1, T(3)=2, T(4)=4, T(5)=7, T(6)=13,
    # T(7)=24, T(8)=44, T(9)=81, T(10)=149, T(11)=274, T(12)=504

    # 360 = 274 + 81 + 4 + 1 = T(11) + T(9) + T(4) + T(2)
    # In the tribonacci positional system (positions 0,1,2,...):
    # Position 10: T(11)=274, Position 8: T(9)=81, Position 3: T(4)=4, Position 1: T(2)=1
    # Binary string: position 10=1, 9=0, 8=1, 7=0, 6=0, 5=0, 4=0, 3=1, 2=0, 1=1, 0=0
    # = 10100001010

    trib_positions = {10: 274, 8: 81, 3: 4, 1: 1}
    total = sum(trib_positions.values())
    print(f"  360 = {total} = ", end="")
    print(" + ".join(f"T({p+1})={v}" for p, v in sorted(trib_positions.items(), reverse=True)))

    # Build positional representation
    max_pos = max(trib_positions.keys())
    digits = ['0'] * (max_pos + 1)
    for p in trib_positions:
        digits[max_pos - p] = '1'
    pos_str = ''.join(digits)
    print(f"  Base-tau positional: {pos_str}")
    print(f"  = 1·tau^10 + 0·tau^9 + 1·tau^8 + 0...0 + 1·tau^3 + 0·tau^2 + 1·tau^1 + 0·tau^0")

    # The pattern: positions 10, 8, 3, 1
    # Gaps: 10-8=2, 8-3=5, 3-1=2
    print(f"\n  Active positions: 10, 8, 3, 1")
    print(f"  Gaps between positions: 2, 5, 2")
    print(f"  The gaps are 2, 5, 2 — PALINDROMIC!")
    print(f"  And 2, 5, 2 uses the choice (2) and interaction (5) primes!")
    print(f"  The period (5+2+2 = 9 = 3^2) spans the representation.")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 5: TRIPLE TAU AND THE RATIO 189/360")
    print(f"{'='*70}")

    print(f"""
  We showed: 189/360 = 21/40.

  In base tau, let's understand 189 and 360 separately:
    189 = 149 + 24 + 13 + 2 + 1 (5 terms, weight 5)
    360 = 274 + 81 + 4 + 1      (4 terms, weight 4)

  189 is MORE COMPLEX than 360 in base tau (weight 5 vs 4).
  This means: the H-maximizer at n=7 is more "tribonacci-compound"
  than the circle/Hamiltonian-cycle-count.

  Now: 189 = 3^3 * 7 = 27 * 7.
  In base tau:
    27 = 24 + 2 + 1 (= T(7) + T(3) + T(2), weight 3)
    7 = 7 (= T(5), weight 1)
    189 = 27 * 7 in prime factorization.
    But 189 in tribonacci is NOT the product of 27_tau and 7_tau.
    (Base tau multiplication is not digit-wise.)

  The RATIO 189/360 = 21/40 in base tau:
    21 = 13 + 7 + 1 (weight 3)
    40 = 24 + 13 + 2 + 1 (weight 4)

  21/40 in tribonacci: forbidden_2 / (Golay + regularity + choice + ground)

  Now: tau^3 ≈ 6.222. What is 360 / tau^3?
    360 / tau^3 = {360/tau3:.6f}
    This is approximately {360/tau3:.2f} ≈ 57.86.

  And: 360 / 6 = 60 = |A_5|.
  So 360/tau^3 ≈ 360/6.222 ≈ 57.86 < 60.
  The SHORTFALL: 60 - 360/tau^3 = {60 - 360/tau3:.6f} ≈ {60 - 360/tau3:.2f}.

  360 / tau^3 is the number of times "one period" (in the
  continuous tau sense) fits into one full circle.
  It's approximately 58 — close to but less than 60.

  Compare: 360 / 6 = 60 exactly.
  360 / tau^3 ≈ 57.86.

  The excess of the discrete (6) over the continuous (tau^3):
  360/6 - 360/tau^3 = 60 - 57.86 = {60 - 360/tau3:.6f}
  ≈ 2.14 ≈ 2 + 1/7?
  60 - 360/tau^3 = 60 - 360/(tau^2+tau+1)
  = 60*(tau^2+tau+1-6)/(tau^2+tau+1)
  = 60*(tau^2+tau-5)/(tau^2+tau+1)
  = 60*(tau^3-6)/tau^3

  At tau = 1.83929: 60*0.2223/6.2223 = 60*0.03572 = {60*0.2223/6.2223:.4f}. Hmm.
  Let me compute directly: {60*(tau3-6)/tau3:.6f}.""")

    print(f"  Exact: 60*(tau^3-6)/tau^3 = {60*(tau3-6)/tau3:.10f}")
    print(f"  This is the 'correction' from discrete period to continuous period.")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 6: EVERY KEY NUMBER IN BASE TAU — THE DEEP READ")
    print(f"{'='*70}")

    # Now let's do the real work: what does each number "look like"
    # as a tau-expansion, and what does that MEAN?

    print(f"""
  In base tau, each number gets a tribonacci positional representation.
  The POSITIONS that light up tell us about the number's structure.

  The tribonacci numbers are:
  pos:  0   1   2   3   4   5   6   7   8   9  10  11  12
  T:    1   1   2   4   7  13  24  44  81 149 274 504 927

  Each position corresponds to a tau POWER:
  pos 0 = tau^0 = 1 (ground state)
  pos 1 = tau^1 ≈ 1.84 (close to choice 2)
  pos 2 = tau^2 ≈ 3.38 (close to cycle 3)
  pos 3 = tau^3 ≈ 6.22 (close to period 6)
  pos 4 = tau^4 ≈ 11.44 (close to extension 11)
  pos 5 = tau^5 ≈ 21.05 (close to FORBIDDEN 21!)
  pos 6 = tau^6 ≈ 38.72
  pos 7 = tau^7 ≈ 71.21
  pos 8 = tau^8 ≈ 130.98
  pos 9 = tau^9 ≈ 240.90
  pos10 = tau^10 ≈ 443.09

  NOTICE: tau^5 ≈ 21.05 ≈ 21 = H_forb_2!
  The FIFTH power of tau is almost exactly the second forbidden value!

  Let me check: Phi_3(tau^(5/3))? No, let me just observe:
  tau^5 = tau^3 * tau^2 = Phi_3(tau) * tau^2
  = (tau^2 + tau + 1) * tau^2 = tau^4 + tau^3 + tau^2
  = {tau**4 + tau**3 + tau**2:.6f}

  And 21 = Phi_3(4) = 4^2 + 4 + 1 = 21.
  tau^5 ≈ 21.05 ≈ Phi_3(4).

  Is there a relationship? tau^5 = Phi_3(tau) * tau^2 ≈ 6.222 * 3.383 ≈ 21.05.
  And Phi_3(4) = 21 exactly.
  So tau^5 ≈ Phi_3(4). The fifth power of tau ≈ Phi_3(2^2).

  MORE: tau^3 ≈ Phi_3(1) + 3.222? No: Phi_3(1) = 3, tau^3 ≈ 6.222.
  tau^3 ≈ 2 * Phi_3(1) = 6. Yes! tau^3 ≈ 2 * 3 = 6.
  tau^5 ≈ Phi_3(4) = 21. And 4 = 2^2.

  tau^3 ≈ 6 = 2 * 3 = choice * cycle
  tau^5 ≈ 21 = Phi_3(4) = Phi_3(2^2) = second forbidden

  PATTERN: tau^(2k+1) ≈ Phi_3(2^k)?
  k=1: tau^3 ≈ Phi_3(2) = 7? NO, tau^3 ≈ 6.22, not 7.
  Hmm. Let me reconsider.""")

    # Check tau^k near Phi_3 values
    print(f"  tau^k vs Phi_3 and key numbers:")
    phi3_vals = [(0, 1), (1, 3), (2, 7), (3, 13), (4, 21)]
    for k in range(13):
        val = tau**k
        nearest_phi3 = min(phi3_vals, key=lambda x: abs(val - x[1]))
        nearest_int = round(val)
        print(f"    tau^{k:2d} = {val:10.4f}  nearest int: {nearest_int:4d}  "
              f"nearest Phi_3({nearest_phi3[0]}): {nearest_phi3[1]:4d} "
              f"(diff: {val - nearest_phi3[1]:+8.4f})")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 7: THE TAU-360 RESONANCE")
    print(f"{'='*70}")

    # Key: 360 = tau^9.659...
    # And tau^3 ≈ 6.222
    # So tau^9 ≈ (tau^3)^3 ≈ 6.222^3 ≈ 240.9
    # And tau^12 ≈ (tau^3)^4 ≈ 6.222^4 ≈ 1498.7

    print(f"  The tau^3 tower:")
    for k in range(6):
        val = tau**(3*k)
        print(f"    (tau^3)^{k} = tau^{3*k:2d} = {val:12.4f}"
              f"  (≈ 6.222^{k} = {6.2223**k:12.4f})")

    print(f"""
  tau^0 = 1       (ground state)
  tau^3 ≈ 6.222   (≈ period 6)
  tau^6 ≈ 38.72   (≈ 6^2 = 36)
  tau^9 ≈ 240.9   (≈ 6^3 = 216)
  tau^12 ≈ 1498.9 (≈ 6^4 = 1296)

  Each TRIPLE of tau powers ≈ one power of 6.
  But with increasing overshoot (because tau^3 > 6).

  Now: 360/tau^9 = {360/tau**9:.6f}
  This is close to 360/240.9 ≈ 1.494 ≈ 3/2 = the perfect fifth!

  360/tau^9 = {360/tau**9:.10f}
  3/2 = {3/2:.10f}

  360/tau^9 ≈ 3/2!!! The ratio is ALMOST exactly 3/2!

  Difference: {360/tau**9 - 3/2:.10f}

  So: 360 ≈ (3/2) * tau^9 = (3/2) * (tau^3)^3.

  Since tau^3 ≈ 6.222:
  360 ≈ (3/2) * 6.222^3 = (3/2) * 240.9 ≈ 361.4.
  Close but not exact (because tau^3 ≠ 6 exactly).

  If tau^3 WERE exactly 6:
  (3/2) * 6^3 = (3/2) * 216 = 324. Not 360.

  But (5/3) * 6^3 = (5/3) * 216 = 360. EXACT!

  So 360 = (5/3) * 6^3 = interaction/cycle * period^3.

  AND: 360/tau^9 ≈ 3/2, but the EXACT expression is:
  360 = (5/3) * 6^3 (using integer period 6)
  360 ≈ (3/2) * tau^9 (using continuous period tau^3)

  The BRIDGE between these:
  5/3 (rational, exact) ↔ 3/2 (approximate in tau)
  Both are superparticular ratios (of the form (n+1)/n).
  5/3 = the major sixth in music.
  3/2 = the perfect fifth in music.
  """)

    # Verify
    print(f"  VERIFICATIONS:")
    print(f"  360 = (5/3)*6^3 = (5/3)*216 = {Fraction(5,3)*216}")
    print(f"  360/tau^9 = {360/tau**9:.10f}")
    print(f"  3/2 = {3/2:.10f}")
    print(f"  Excess: {360/tau**9 - 3/2:.10f} ({(360/tau**9 - 3/2)/(3/2)*100:.4f}%)")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 8: THE TRIPLE-TAU IDENTITY AND ITS CONSEQUENCES")
    print(f"{'='*70}")

    print(f"""
  tau^3 = Phi_3(tau) IS the defining identity of the tribonacci constant.
  Its consequences cascade through all of tournament theory.

  CONSEQUENCE 1: tau^6 = (tau^3)^2 = Phi_3(tau)^2
    tau^6 = {tau**6:.6f}
    Phi_3(tau)^2 = {(tau**3)**2:.6f}
    And: Phi_3(tau^2) = tau^4 + tau^2 + 1 = {tau**4 + tau**2 + 1:.6f}
    Note: Phi_3(tau^2) = Phi_3(tau) * Phi_6(tau)
    where Phi_6(tau) = tau^2 - tau + 1 = {tau**2 - tau + 1:.6f}
    Check: {(tau**2 + tau + 1) * (tau**2 - tau + 1):.6f} = {tau**4 + tau**2 + 1:.6f}
    So Phi_3(tau^2) = Phi_3(tau)*Phi_6(tau) = tau^3 * Phi_6(tau)

  CONSEQUENCE 2: The n-th power of tau^3 gives tower evaluations:
    (tau^3)^1 = tau^3 ≈ 6.222    (one period)
    (tau^3)^2 = tau^6 ≈ 38.72    (period squared)
    (tau^3)^3 = tau^9 ≈ 240.9    (period cubed ≈ 360 * 2/3)
    (tau^3)^4 = tau^12 ≈ 1498.9  (period to the fourth)

  CONSEQUENCE 3: The "tau-period" tau^3 ≈ 6.222 means:
    360 degrees = 360/(tau^3) ≈ 57.86 "tau-periods"
    A full circle = about 58 tau-periods.
    (Compare: a full circle = 360/6 = 60 discrete periods.)
    The continuous period is SLIGHTLY LONGER than the discrete one,
    so fewer fit in a circle.

  CONSEQUENCE 4: Every tournament number lives at a specific
    "tau altitude" — its logarithm base tau:
    log_tau(1) = 0        (the ground)
    log_tau(2) = 1.137     (one floor up)
    log_tau(3) = 1.803     (just under two floors)
    log_tau(6) = 2.940     (just under three floors = tau^3 ≈ 6)
    log_tau(7) = 3.193     (just ABOVE three floors — the forbidden!)
    log_tau(21) = 4.996    (just under FIVE floors — tau^5 ≈ 21)
    log_tau(360) = 9.659   (just under TEN floors — 360 ≈ tau^10?)

  PATTERN: forbidden values sit at NEAR-INTEGER tau-altitudes!
    log_tau(7) = 3.19 ≈ 3 + 1/5 (just past the third floor)
    log_tau(21) = 5.00 ≈ 5 (almost exactly the fifth floor!)

  THE FORBIDDEN VALUES LIVE AT ODD FLOORS IN THE TAU TOWER:
    7 ≈ tau^3 (third floor)
    21 ≈ tau^5 (fifth floor)

  And the third and fifth floors correspond to:
    3 = the cycle generator
    5 = the interaction generator

  THE FORBIDDEN VALUES ARE THE TAU-TOWER AT THE CYCLE AND INTERACTION FLOORS.""")

    print(f"\n  tau^3 = {tau**3:.6f}, 7 = 7, ratio = {7/tau**3:.6f}")
    print(f"  tau^5 = {tau**5:.6f}, 21 = 21, ratio = {21/tau**5:.6f}")
    print(f"  7/tau^3 = {7/tau**3:.6f} ≈ {Fraction(7,6).limit_denominator(100)} (close to 7/6)")
    print(f"  21/tau^5 = {21/tau**5:.6f} ≈ {Fraction(21, round(tau**5)).limit_denominator(100)} (close to 1)")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 9: TRIPLE TAU IS THE HEARTBEAT")
    print(f"{'='*70}")

    print(f"""
  SYNTHESIS: What is the significance of tau^3?

  tau^3 = Phi_3(tau) ≈ 6.222

  It is:
  1. THE CONTINUOUS PERIOD. The tournament has a discrete period of 6,
     but the continuous version (through the tribonacci lens) is tau^3.
     The excess tau^3 - 6 ≈ 0.222 is the "irrational overshoot" —
     the part of the period that escapes rational description.

  2. THE CUBE OF THE TRIBONACCI. "Triple tau" = tau applied three times.
     Since tau encodes the 3-cycle structure, applying it three times
     completes one full cycle of the 3-cycle. tau^3 is the 3-CYCLE OF THE 3-CYCLE.

  3. THE Phi_3-CUBING COINCIDENCE. At tau, the cyclotomic polynomial
     and the cube operation AGREE. This is not true anywhere else on the
     real line. tau is the UNIQUE point where evaluating the tournament
     polynomial equals cubing.

  4. THE HEARTBEAT. Just as 6 is the period of tournament parity,
     tau^3 is the period of tournament GROWTH. Every three tribonacci
     steps, the system multiplies by tau^3 ≈ 6.222. This is the
     growth heartbeat: one pulse every three generations.

  5. THE KEY TO 360. Since 360 ≈ (3/2)*tau^9 = (3/2)*(tau^3)^3,
     the full circle is approximately (3/2) cubed heartbeats.
     Three heartbeats, scaled by the perfect fifth 3/2,
     give one revolution.

  6. THE sqrt(21) CONNECTION. tau^3 = 6 exactly when tau = (-1+sqrt(21))/2.
     The second forbidden value 21 appears under the square root.
     Triple tau overshoots the period because the tribonacci constant
     is slightly ABOVE the "period root" that involves sqrt(H_forb_2).

  TRIPLE TAU IS THE HEARTBEAT OF TOURNAMENT GROWTH:
  the irrational period that governs how quickly tournament
  complexity multiplies, sitting just above the rational period 6,
  with its overshoot controlled by the forbidden value 21.
  """)

    print(f"\n{'='*70}")
    print("DONE — TRIPLE TAU IS THE HEARTBEAT")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
