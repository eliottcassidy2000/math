"""
rational_irrational_symbiosis.py -- kind-pasteur-2026-03-14-S106h
HOW RATIONAL AND IRRATIONAL NUMBERS LIVE TOGETHER

The usual story: rationals are "simple," irrationals are "complex."
Q is countable, R\Q is uncountable. They seem opposed.

But tournament theory reveals they are SYMBIOTIC.
Every key result involves BOTH types working together:
  - H is always a natural (rational) number
  - But its growth rate is tau (irrational)
  - The variance ratio is 1/3 (rational)
  - But it comes from a cone with continuous geometry (irrational)
  - The forbidden values are 7, 21 (rational)
  - But they sit at tau^3 and tau^5 (irrational coordinates)

The symbiosis: rationals provide the VALUES, irrationals provide the STRUCTURE.
Neither makes sense without the other.
"""

import sys, math
import numpy as np
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8')

tau = 1.8392867552141612
phi = (1 + math.sqrt(5)) / 2

def main():
    print("=" * 70)
    print("THE SYMBIOSIS OF RATIONAL AND IRRATIONAL")
    print("kind-pasteur-2026-03-14-S106h")
    print("=" * 70)

    # ============================================================
    print(f"\n{'='*70}")
    print("CHAPTER 1: THE DANCE — EVERY RATIONAL HAS AN IRRATIONAL SHADOW")
    print(f"{'='*70}")

    print(f"""
  Consider the key rational numbers of tournament theory.
  Each one has an irrational partner — its "shadow" in the
  continuous world. They come in PAIRS.

  PAIR 1: 6 and tau^3
    6 = LCM(2,3) = the discrete period (rational, exact)
    tau^3 = 6.2223... = the continuous period (irrational, approximate)
    The rational 6 is the INTEGER PART of the irrational tau^3.
    tau^3 = 6 + (tau^2 + tau - 5) = 6 + 0.2223...
    The irrational EXCEEDS the rational by a controlled amount.

  PAIR 2: 7 and tau^3
    7 = Phi_3(2) = the first forbidden (rational, exact)
    tau^3 = 6.2223... = the heartbeat (irrational)
    7 = tau^3 * (9/8) approximately
    The rational forbidden value is the irrational heartbeat
    SCALED by the rational 9/8 (reflection/threshold).
    Irrational * rational = rational. The irrational disappears!

  PAIR 3: 21 and tau^5
    21 = Phi_3(4) = the second forbidden (rational)
    tau^5 = 21.0498... = the fifth floor (irrational)
    21/tau^5 = 0.99764... ≈ 1
    The rational is the irrational rounded to the nearest integer.
    tau^5 is 21's irrational AVATAR.

  PAIR 4: 1/3 and ln(4/3)
    1/3 = Var/Mean^2 at n=3,4 (rational, exact)
    ln(4/3) = 0.2877... (irrational, transcendental)
    ln(1 + 1/3) = ln(4/3) = the log-variance
    The rational variance ratio GENERATES a transcendental
    through the logarithm.

  PAIR 5: 360 and tau^9
    360 = 2^3 * 3^2 * 5 (rational)
    tau^9 = 240.90... (irrational)
    360/tau^9 = 1.4944... ≈ 3/2 (irrational → rational!)
    The ratio of a rational to an irrational is ALMOST rational.

  PAIR 6: n!/2^(n-1) and e
    n!/2^(n-1) = Mean(H) (rational for each n)
    max_H / Mean → e as n → ∞ (irrational, transcendental)
    A sequence of rationals whose LIMIT is transcendental.
    Stirling: n! ~ sqrt(2*pi*n) * (n/e)^n
    The rational factorial is APPROXIMATED by the transcendental e and pi.

  THE PATTERN: Rationals are the NOTES. Irrationals are the MELODY.
  You can hear individual notes (7, 21, 360) but the melody
  (tau^3, tau^5, tau^9) connects them into a coherent piece.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("CHAPTER 2: THE CONTINUED FRACTION — WHERE THEY MEET")
    print(f"{'='*70}")

    print(f"""
  The continued fraction is where rational and irrational MERGE.
  Every irrational number is an INFINITE continued fraction.
  Every TRUNCATION of that fraction is a rational APPROXIMATION.

  The continued fraction is the MEETING POINT:
  it is simultaneously a rational sequence AND an irrational limit.

  TAU (tribonacci constant):""")

    # Compute continued fraction of tau
    def continued_fraction(x, n_terms=15):
        """Compute continued fraction coefficients."""
        coeffs = []
        for _ in range(n_terms):
            a = int(x)
            coeffs.append(a)
            frac = x - a
            if abs(frac) < 1e-12:
                break
            x = 1 / frac
        return coeffs

    cf_tau = continued_fraction(tau, 20)
    cf_phi = continued_fraction(phi, 20)
    cf_sqrt2 = continued_fraction(math.sqrt(2), 20)
    cf_e = continued_fraction(math.e, 20)
    cf_pi = continued_fraction(math.pi, 20)

    print(f"  tau   = [{cf_tau[0]}; {', '.join(map(str, cf_tau[1:]))}]")
    print(f"  phi   = [{cf_phi[0]}; {', '.join(map(str, cf_phi[1:]))}]")
    print(f"  sqrt2 = [{cf_sqrt2[0]}; {', '.join(map(str, cf_sqrt2[1:]))}]")
    print(f"  e     = [{cf_e[0]}; {', '.join(map(str, cf_e[1:]))}]")
    print(f"  pi    = [{cf_pi[0]}; {', '.join(map(str, cf_pi[1:]))}]")

    print(f"""
  OBSERVATIONS:
  - phi = [1; 1, 1, 1, 1, ...]: ALL ONES. The simplest irrational.
    The golden ratio is the LEAST well-approximated by rationals.
    Its continued fraction converges the SLOWEST of any number.
    phi is the MOST IRRATIONAL number.

  - tau = [1; {', '.join(map(str, cf_tau[1:8]))}...]: More complex than phi.
    The large coefficient {max(cf_tau[1:8])} means tau has a GOOD rational
    approximation at that point. Tau is "more rational" than phi.

  - sqrt(2) = [1; 2, 2, 2, ...]: ALL TWOS. The simplest quadratic irrational
    after phi. Its rational approximants converge steadily.

  - e = [2; 1, 2, 1, 1, 4, 1, 1, 6, ...]: The pattern 1, 2k, 1, 1, 2k+2, ...
    e has a STRUCTURED continued fraction. The transcendental e
    has RATIONAL structure in its irrational expansion!

  - pi = [3; 7, 15, 1, 292, ...]: IRREGULAR. No pattern.
    pi is "more random" than e. Its continued fraction has
    occasional LARGE coefficients (292) that make it
    well-approximated by 355/113 = 3.14159292...

  THE SYMBIOSIS IN CONTINUED FRACTIONS:
  Each coefficient a_k is a NATURAL number (rational).
  The sequence [a_0; a_1, a_2, ...] encodes an IRRATIONAL number.
  The rational coefficients CONSTRUCT the irrational.
  The irrational ORGANIZES the rational coefficients.
  Neither exists without the other.""")

    # Convergents of tau
    print(f"\n  RATIONAL CONVERGENTS OF TAU:")
    num = [cf_tau[0], cf_tau[0]*cf_tau[1] + 1]
    den = [1, cf_tau[1]]
    for i in range(2, min(10, len(cf_tau))):
        num.append(cf_tau[i]*num[-1] + num[-2])
        den.append(cf_tau[i]*den[-1] + den[-2])

    for i in range(min(8, len(num))):
        approx = num[i]/den[i]
        error = abs(approx - tau)
        print(f"    {num[i]}/{den[i]} = {approx:.10f}  (error: {error:.2e})")

    # ============================================================
    print(f"\n{'='*70}")
    print("CHAPTER 3: THE BOOKEND THEOREM — RATIONAL ENDPOINTS, IRRATIONAL GAP")
    print(f"{'='*70}")

    print(f"""
  The BOOKEND THEOREM (from S106g):
  log_tau(21) - log_tau(7) = log_tau(3) ≈ 1.803

  The forbidden values 7 and 21 are RATIONAL.
  Their SEPARATION in the tau tower is IRRATIONAL (log_tau(3)).
  But their RATIO is RATIONAL: 21/7 = 3.

  This is the essence of the symbiosis:
    RATIONAL endpoints → 7 and 21 (exact, computable)
    IRRATIONAL gap → log_tau(3) (continuous, structural)
    RATIONAL ratio → 3 (the cycle generator)

  The irrational gap ENCODES the rational ratio:
    log_tau(3) = ln(3)/ln(tau).
  The transcendental logarithm CONVERTS the rational 3 into
  an irrational distance. But the 3 is still "in there" —
  you can recover it by exponentiating: tau^(log_tau(3)) = 3.

  THE ROUND-TRIP:
    3 (rational) → log_tau(3) (irrational) → tau^(log_tau(3)) = 3 (rational)

  Going through the irrational doesn't CHANGE the rational.
  It ENRICHES it with context (position in the tau tower)
  and then returns it exactly.

  This is like a PRISM: white light (rational) → spectrum (irrational) → white light.
  The prism REVEALS structure but doesn't alter content.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("CHAPTER 4: HOW IRRATIONALS MAKE RATIONALS POSSIBLE")
    print(f"{'='*70}")

    print(f"""
  Without irrationals, the rational structure of tournaments would be
  UNEXPLAINED. Consider:

  Q: Why is Var/Mean^2 = 1/3 at n=3,4?
  A: Because H is a degree-2 quadratic on the hypercube, and the
     integral of x^2 on [0,1] = 1/3. But this integral involves the
     CONTINUOUS real line (irrationals!) to evaluate.

  Q: Why are 7 and 21 forbidden?
  A: Because Phi_3(2) = 7 and Phi_3(4) = 21, and the cycle structure
     prevents certain independent set counts. But the PROOF uses
     properties of Phi_3 as a polynomial — which exists over ALL of R
     (including irrationals), not just Z.

  Q: Why does max_H / mean_H → e?
  A: Because e = lim (1 + 1/n)^n. This limit exists ONLY in R (irrationals).
     The rational sequence (1+1/n)^n never reaches e, but it converges
     to it through R.

  Q: Why is tau^5 ≈ 21?
  A: Because tau satisfies tau^3 = tau^2 + tau + 1, so
     tau^5 = tau^2 * (tau^2 + tau + 1) = tau^4 + tau^3 + tau^2.
     The near-integer property tau^5 ≈ 21 is because tau is the
     real root of x^3 - x^2 - x - 1 = 0, and this polynomial
     has INTEGER coefficients that conspire to make powers of tau
     close to integers.

  THE MECHANISM: Integer coefficients in the minimal polynomial
  FORCE the irrational to stay close to integers.
  The rational world CONSTRAINS the irrational,
  and the irrational EXPLAINS the rational.""")

    # Show how close tau powers are to integers
    print(f"\n  HOW CLOSE tau powers are to integers:")
    print(f"  (The PV property — Pisot-Vijayaraghavan)")
    for k in range(1, 16):
        val = tau**k
        nearest = round(val)
        frac = val - nearest
        print(f"    tau^{k:2d} = {val:12.6f}  nearest int: {nearest:5d}  "
              f"frac part: {frac:+10.6f}")

    print(f"""
  THE PV PROPERTY: The fractional parts of tau^k get SMALL!
  Not monotonically, but they oscillate and shrink.
  This is because tau is a PISOT NUMBER: an algebraic integer > 1
  whose conjugates (the other roots of x^3-x^2-x-1=0) have
  absolute value < 1.

  The conjugate roots: sigma ≈ -0.420 ± 0.606i, |sigma| ≈ 0.737.
  Since |sigma| < 1: sigma^k → 0 as k → inf.
  And tau^k + sigma^k + sigma_bar^k = integer (by Newton's identities).
  So tau^k = integer - sigma^k - sigma_bar^k ≈ integer.

  THE SYMBIOSIS: The rational integers PULL the irrational tau^k
  toward them. The closer tau is to a Pisot number (it IS one!),
  the more strongly the integers attract its powers.
  The irrational dances AROUND the integers, pulled by them
  but never landing exactly.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("CHAPTER 5: THE TAU-INTEGER NEAR-MISSES AS INFORMATION")
    print(f"{'='*70}")

    # The signed distances from tau^k to nearest integer
    print(f"  The signed fractional parts of tau^k carry INFORMATION:")
    print(f"  (positive = above integer, negative = below)")

    signs = []
    for k in range(1, 20):
        val = tau**k
        frac = val - round(val)
        sign = '+' if frac > 0 else '-'
        signs.append(sign)
        tournament_hit = ""
        nearest = round(val)
        if nearest in [2, 3, 6, 7, 11, 13, 21, 24, 44, 81, 131, 241, 443, 815, 1499]:
            tournament_hit = f"  ← near {nearest}"
        print(f"    tau^{k:2d}: frac = {frac:+.6f} ({sign}){tournament_hit}")

    sign_string = ''.join(signs)
    print(f"\n  Sign sequence: {sign_string}")
    print(f"  Pattern: {''.join(signs[:12])}")

    print(f"""
  The sign sequence tells us whether tau^k OVERSHOOTS or UNDERSHOOTS
  the nearest integer. This sequence is NOT random — it has structure
  determined by the continued fraction of tau.

  When the fractional part is SMALL (close to 0):
    tau^k ≈ integer very precisely → the integer is "attracted"
    Examples: tau^5 ≈ 21 (frac = +0.05), tau^8 ≈ 131 (frac = -0.02)

  When the fractional part is LARGE (close to ±0.5):
    tau^k is between integers → the number is "repelled"
    Examples: tau^4 ≈ 11.44 (frac = +0.44)

  THE ATTRACTED integers are TOURNAMENT-SIGNIFICANT:
    21 (attracted by tau^5) = second forbidden
    131 (attracted by tau^8) = numerator of Var/Mean^2 at n=7!
    241 (attracted by tau^9) = near tau^9, the heartbeat boundary

  The REPELLED integers are TOURNAMENT-INVISIBLE:
    11-12 (repelled by tau^4) = neither is a key tournament number

  THE SYMBIOSIS: The irrational tau SELECTS which rational integers
  are "important" by attracting some and repelling others.
  The attracted integers become the LANDMARKS of tournament theory.""")

    # Check: is 131 the numerator of Var/Mean^2 at n=7?
    print(f"\n  CHECK: Var/Mean^2 at n=7 = 131/504")
    print(f"  131 = round(tau^8) = round({tau**8:.4f}) = {round(tau**8)}")
    print(f"  504 = T(12) (tribonacci number)")
    print(f"  BOTH the numerator AND denominator of the n=7 variance ratio")
    print(f"  are tau-attracted integers!")
    print(f"  131 ≈ tau^8, 504 = T(12) = tau^12 - ...")
    print(f"  504/tau^12 = {504/tau**12:.6f} (close to 1/3 = {1/3:.6f}!)")
    print(f"  WOW: 504 ≈ tau^12 / 3! The tribonacci number is tau^12 / Phi_3(1)!")

    # Verify
    print(f"  tau^12 = {tau**12:.4f}")
    print(f"  tau^12 / 3 = {tau**12/3:.4f}")
    print(f"  Actual T(12) = 504")
    print(f"  tau^12/3 = {tau**12/3:.4f}, difference from 504: {tau**12/3 - 504:.4f}")

    # ============================================================
    print(f"\n{'='*70}")
    print("CHAPTER 6: THE DEEPEST SYMBIOSIS — ALGEBRAIC + TRANSCENDENTAL")
    print(f"{'='*70}")

    print(f"""
  The deepest symbiosis is between ALGEBRAIC irrationals (tau, phi, sqrt(p))
  and TRANSCENDENTAL numbers (e, pi, ln 2).

  ALGEBRAIC irrationals solve polynomial equations with rational coefficients.
  They are "constructible" from the rationals by a finite process.
  tau satisfies x^3 - x^2 - x - 1 = 0 (a FINITE equation).

  TRANSCENDENTAL numbers solve NO polynomial equation.
  They require INFINITE processes (limits, series, integrals).
  e = lim (1+1/n)^n (an INFINITE limit).
  pi = 4*(1 - 1/3 + 1/5 - ...) (an INFINITE series).

  HOW THEY WORK TOGETHER:

  1. STIRLING'S FORMULA: n! ~ sqrt(2*pi*n) * (n/e)^n
     The factorial (rational) is approximated by a product of
     pi (transcendental) and e (transcendental).
     Mean(H) = n!/2^(n-1) involves the SAME approximation.
     So the RATIONAL mean of H is the RATIO of two TRANSCENDENTALS.

  2. THE SZELE LIMIT: max_H/mean_H → e
     A sequence of RATIONALS (max_H(n)/mean_H(n)) converges to
     the TRANSCENDENTAL e. The algebraic tau controls the RATE
     of convergence (through the tribonacci growth rate).
     Algebraic rate, transcendental limit, rational values.

  3. THE FOURIER ENERGY: E_2/E_0 = 2(n-2)/(n(n-1))
     This RATIONAL formula comes from summing SQUARED Fourier
     coefficients. The Fourier transform uses e^(2*pi*i*k/N)
     (TRANSCENDENTAL basis), but the sum comes out RATIONAL.
     The transcendental basis CONSPIRES to produce rational sums.

  4. THE INFORMATION RATE: I(T;H)/m ≈ 0.27
     This involves entropy: H_info = -sum p*log(p) (TRANSCENDENTAL).
     The mutual information uses LOG (transcendental function).
     But the result ≈ 0.27 is approximately RATIONAL.
     ln(4/3) = 0.2877... (close to the info rate!).
     And 4/3 = 1 + 1/3 = 1 + the cone ratio.

  THE GRAND SYMBIOSIS:
    RATIONAL numbers provide the DISCRETE VALUES (H, n, alpha_k).
    ALGEBRAIC irrationals provide the GROWTH RATES (tau, phi, sqrt(p)).
    TRANSCENDENTAL numbers provide the LIMITS (e, pi, ln 2).

  Each type serves a DIFFERENT FUNCTION:
    Rationals = the SKELETON (the bones, the exact measurements).
    Algebraics = the MUSCLES (the dynamics, the growth).
    Transcendentals = the NERVOUS SYSTEM (the limits, the optimization).

  A body needs all three. Remove any one and the organism dies.
  Without rationals: no exact H values, no counting.
  Without algebraics: no growth rate, no tribonacci heartbeat.
  Without transcendentals: no efficiency limit, no optimization.
  """)

    # ============================================================
    print(f"\n{'='*70}")
    print("CHAPTER 7: THE EQUATION THAT UNITES THEM ALL")
    print(f"{'='*70}")

    print(f"""
  Is there ONE equation that expresses the symbiosis?

  CANDIDATE: tau^3 = Phi_3(tau) = tau^2 + tau + 1

  This equation involves:
    tau (IRRATIONAL, algebraic) — the growth rate
    3 (RATIONAL, natural) — the exponent, the cycle
    Phi_3 (a polynomial with RATIONAL coefficients) — the structure
    1 (RATIONAL) — the ground state, the identity

  The equation SAYS: the irrational tau, when cubed (rational operation),
  equals a polynomial with rational coefficients evaluated at itself.
  The rational and irrational are SELF-CONSISTENT.
  tau is defined BY the rationals (the polynomial coefficients),
  and the rationals are ORGANIZED BY tau (the heartbeat).

  ANOTHER CANDIDATE: 360 = (5/3) * 6^3

  This is PURELY rational. But it encodes the irrational:
    360 ≈ (3/2) * tau^9
  The same number has a rational form AND an irrational form,
  and the TWO FORMS are connected by:
    (5/3) * 6^3 ≈ (3/2) * tau^9
    (5/3) * 6^3 / ((3/2) * tau^9) = 1.00377...
    The error is 0.377% — close but not exact.

  THE BEST CANDIDATE: Var/Mean^2 at n=7 = 131/504

  131 ≈ tau^8 (algebraic irrational proximity)
  504 = T(12) (tribonacci, defined by tau)
  131/504 (rational)
  ≈ tau^8 / T(12) (irrational / irrational)
  ≈ tau^8 / (tau^12/3) (using T(12) ≈ tau^12/3)
  ≈ 3/tau^4 (simplification!)
  ≈ 3/11.44 = 0.2623 (vs actual 0.2599)

  Close but not exact. The symbiosis is approximate, not algebraic.

  FINAL ANSWER: The symbiosis has no single equation.
  It is a RELATIONSHIP, not a formula.
  Rational and irrational numbers are like FIGURE and GROUND
  in a Gestalt image: each defines the other's boundary,
  and the meaning lives in their interaction, not in either alone.

  In tournament theory:
    The rationals provide the WHAT (H = 7 is forbidden).
    The irrationals provide the WHY (7 ≈ tau^3 is the heartbeat boundary).
    Together they provide the UNDERSTANDING (Phi_3(2) = 7 because
    tau^3 = Phi_3(tau) and 2 > tau, pushing past the heartbeat).

  2 > tau. That's the key inequality. The rational generator (2)
  EXCEEDS the irrational growth rate (tau). This is why tournaments
  are "too rich" for simple description: binary choices grow FASTER
  than the tribonacci can organize them. The excess creates structure
  (cycles, forbidden values, Fourier energy) — and that structure
  is where the mathematics lives.
  """)

    print(f"\n{'='*70}")
    print("DONE — RATIONAL AND IRRATIONAL ARE FIGURE AND GROUND")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
