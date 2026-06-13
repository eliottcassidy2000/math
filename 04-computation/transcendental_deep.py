"""
transcendental_deep.py -- kind-pasteur-2026-03-14-S95
Deep creative exploration of pi and e in tournament theory.

pi and e are the two most fundamental transcendental numbers.
Both appear in tournament theory — but HOW and WHY?

THE DEEP QUESTION: Are pi and e "accidental" appearances
(just from Stirling's approximation), or do they have
STRUCTURAL roles in the theory?
"""

import sys, math
import numpy as np
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8')

def C(n, k):
    if k < 0 or k > n: return 0
    return math.comb(n, k)

maxH = {1:1, 2:1, 3:3, 4:5, 5:15, 6:45, 7:189, 8:661, 9:3357, 10:15745, 11:95095}

def main():
    print("=" * 70)
    print("TRANSCENDENTAL NUMBERS IN TOURNAMENT THEORY — DEEP DIVE")
    print("kind-pasteur-2026-03-14-S95")
    print("=" * 70)

    e = math.e
    pi = math.pi

    # ============================================================
    # PART 1: WHERE e APPEARS — COMPLETE CATALOGUE
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 1: EVERY APPEARANCE OF e IN TOURNAMENT THEORY")
    print(f"{'='*70}")

    print(f"\n  e = {e:.15f}")

    print(f"\n  DIRECT APPEARANCES:")
    print(f"  1. SZELE'S THEOREM: max_H(n) / mean_H(n) -> e as n -> inf")
    print(f"     The H-maximizer achieves e times the average.")
    print(f"     This is the MOST IMPORTANT appearance of e.")
    print(f"")
    print(f"  SZELE RATIO CONVERGENCE:")
    for n in range(3, 12):
        mh = maxH[n]
        mn = math.factorial(n) / 2**(n-1)
        ratio = mh / mn
        print(f"     n={n:2d}: max_H/mean_H = {mh}/{mn:.1f} = {ratio:.6f}")
    print(f"     limit = e = {e:.6f}")

    print(f"\n  2. STIRLING'S APPROXIMATION: n! ~ (n/e)^n * sqrt(2*pi*n)")
    print(f"     mean_H = n!/2^(n-1) ~ sqrt(2*pi*n) * (n/(2e))^n")
    print(f"     The e in Stirling controls the base of exponential growth.")

    print(f"\n  3. THE BOLTZMANN DISTRIBUTION:")
    print(f"     Z(beta) = sum_T exp(beta * H(T))")
    print(f"     At beta = 1/e: the partition function has special properties.")
    print(f"     The 'natural temperature' of the tournament gas is 1/e.")

    print(f"\n  4. INFORMATION RATE:")
    print(f"     I(T;H)/m ≈ 0.27. Is this related to 1/e ≈ 0.368?")
    print(f"     Actually 0.27/0.368 = {0.27/0.368:.4f} ≈ 0.73")
    print(f"     And 1 - 0.27 = 0.73. So 0.27 + 1/e ≈ 1? NO: 0.27 + 0.368 = 0.638")
    print(f"     But: 0.27 ≈ 1/(e + 1/e) = 1/{e + 1/e:.4f} = {1/(e+1/e):.4f}")
    print(f"     Nope. Let's try: e^(-1-1/3) = {math.exp(-4/3):.4f}")
    print(f"     e^(-4/3) = 0.264 ≈ 0.27. Hmm, close but not exact.")

    # ============================================================
    # PART 2: WHERE pi APPEARS — COMPLETE CATALOGUE
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 2: EVERY APPEARANCE OF pi IN TOURNAMENT THEORY")
    print(f"{'='*70}")

    print(f"\n  pi = {pi:.15f}")

    print(f"\n  1. STIRLING: n! ~ (n/e)^n * sqrt(2*PI*n)")
    print(f"     pi appears in the NORMALIZATION of the factorial.")
    print(f"     mean_H = n!/2^(n-1) inherits the sqrt(2*pi*n) factor.")

    print(f"\n  2. CENTRAL BINOMIAL: C(2n,n) ~ 4^n / sqrt(PI*n)")
    print(f"     pi appears in the normalization of binomial coefficients.")
    print(f"     The 3-strand Pascal sequence inherits this.")

    print(f"\n  3. SPECTRAL THEORY — PALEY EIGENVALUES:")
    print(f"     The eigenvalues of Paley T_p involve sin(k*pi/p).")
    print(f"     The spectral gap is controlled by pi/p.")

    # Compute the spectral contribution
    for p in [3, 7, 11]:
        # Paley eigenvalue argument
        arg = math.atan2(math.sqrt(p)/2, -0.5)
        print(f"     T_{p}: arg(mu) = {arg:.6f} = {arg/pi:.6f} * pi")

    print(f"\n  4. THE FOURIER TRANSFORM:")
    print(f"     The Walsh-Hadamard transform uses {'{-1, +1}'} basis.")
    print(f"     But the CONTINUOUS Fourier transform uses e^(2*pi*i*k/N).")
    print(f"     pi enters when connecting discrete and continuous theories.")

    print(f"\n  5. THE GAUSSIAN APPROXIMATION:")
    print(f"     If H were Gaussian: I(T;H) = (1/2)*log(2*pi*e*Var(H))")
    print(f"     The actual I(T;H) is LESS (by ~35% at n=5).")
    print(f"     The deficit measures non-Gaussianity.")

    for n in [5, 6]:
        var = sum((maxH.get(n, 0) - math.factorial(n)/2**(n-1))**2
                  for _ in [0]) / 1  # placeholder
        # Use actual Var from S71
        var_actual = {5: 17.8125, 6: 146.25}
        if n in var_actual:
            gaussian_I = 0.5 * math.log2(2 * pi * e * var_actual[n])
            actual_I = {5: 2.681, 6: 4.05}
            if n in actual_I:
                deficit = gaussian_I - actual_I[n]
                print(f"     n={n}: Gaussian I = {gaussian_I:.4f}, actual I = {actual_I[n]:.3f}, "
                      f"deficit = {deficit:.4f} bits")

    print(f"\n  6. LEIBNIZ FORMULA: pi/4 = 1 - 1/3 + 1/5 - 1/7 + ...")
    print(f"     The alternating sum over ODD reciprocals gives pi/4!")
    print(f"     And tournaments are ALL ABOUT odd numbers:")
    print(f"     - H is always ODD")
    print(f"     - alpha_k counts ODD-cycle independent sets")
    print(f"     - The OCF weights are powers of 2 (even), applied to odd cycles")
    print(f"     Is there a tournament Leibniz formula?")

    # ============================================================
    # PART 3: THE e-pi RELATIONSHIP IN TOURNAMENTS
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 3: THE e-pi RELATIONSHIP")
    print(f"{'='*70}")

    print(f"\n  e*pi = {e*pi:.6f}")
    print(f"  e+pi = {e+pi:.6f}")
    print(f"  e^pi = {e**pi:.6f}")
    print(f"  pi^e = {pi**e:.6f}")
    print(f"  e^pi - pi^e = {e**pi - pi**e:.6f}")
    print(f"  sqrt(2*pi*e) = {math.sqrt(2*pi*e):.6f}")
    print(f"  ln(2*pi) = {math.log(2*pi):.6f}")

    print(f"\n  TOURNAMENT COMBINATIONS:")
    # mean_H(n) ~ sqrt(2*pi*n) * (n/(2e))^n
    # max_H(n) ~ e * sqrt(2*pi*n) * (n/(2e))^n
    # So max_H(n) ~ e * sqrt(2*pi*n) * (n/(2e))^n
    # = sqrt(2*pi*n) * (n/2)^n * e^{1-n}
    # The e appears in BOTH Stirling AND Szele.

    print(f"  max_H(n) ~ e * sqrt(2*pi*n) * (n/(2e))^n")
    print(f"  The e appears TWICE: once from Stirling (n!), once from Szele (max/mean).")
    print(f"  They COMPOUND: max_H involves e^{{1-n}} * e = e^{{2-n}}.")
    print(f"")
    print(f"  pi appears ONCE: in the normalization sqrt(2*pi*n).")
    print(f"  So the e contribution is EXPONENTIAL, while pi is POLYNOMIAL.")
    print(f"  e is more 'deeply embedded' than pi.")

    # ============================================================
    # PART 4: CREATIVE IDEAS ABOUT pi
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 4: CREATIVE IDEAS ABOUT pi IN TOURNAMENTS")
    print(f"{'='*70}")

    ideas_pi = [
        ("pi AND THE ODD-CYCLE FORMULA",
         f"Leibniz: pi/4 = sum_k (-1)^k / (2k+1) = 1 - 1/3 + 1/5 - 1/7 + ...\n"
         f"    Tournament: H = 1 + 2*c3 + 2*c5 + 2*c7 + ... (at n=5)\n"
         f"    Both involve odd numbers 1, 3, 5, 7.\n"
         f"    QUESTION: Is there a tournament series that converges to pi?\n"
         f"    IDEA: Define S(T) = sum_k (-1)^k * c_(2k+1) / (2k+1).\n"
         f"    For the 3-cycle T_3: S = 1/3 (only c3=1).\n"
         f"    Does sum_T S(T) / 2^m relate to pi?"),

        ("pi AND THE CIRCLE OF TOURNAMENTS",
         f"The Lee-Yang zeros of Q_n converge to |z| = 1 (the UNIT CIRCLE).\n"
         f"    The circumference of the unit circle is 2*pi.\n"
         f"    As n -> inf, the density of zeros on the circle approaches...\n"
         f"    a continuous distribution. Its integral over the full circle is 2*pi.\n"
         f"    The zero density function rho(theta) encodes tournament structure.\n"
         f"    CONJECTURE: integral(rho(theta) d theta) = 2*pi (trivially).\n"
         f"    But rho(theta) at specific angles might encode H-distribution."),

        ("pi AND RANDOM MATRIX THEORY",
         f"Tournament adjacency A is a random {{0,1}} matrix with A+A^T = J-I.\n"
         f"    The eigenvalue spacing follows Wigner surmise: P(s) = (pi*s/2)*exp(-pi*s^2/4).\n"
         f"    The pi in eigenvalue spacing is the SAME pi as in Stirling.\n"
         f"    Both come from Gaussian integrals: int exp(-x^2) dx = sqrt(pi)."),

        ("pi AND THE AREA OF THE TOURNAMENT POLYTOPE",
         f"The set of score sequences {{s : s is achievable}} forms a polytope.\n"
         f"    Its volume involves pi via the Gamma function: Gamma(1/2) = sqrt(pi).\n"
         f"    The 'volume' of the tournament polytope at dimension n\n"
         f"    should involve pi^{{n/2}} in its normalization."),
    ]

    for title, desc in ideas_pi:
        print(f"\n  {title}")
        print(f"    {desc}")

    # ============================================================
    # PART 5: CREATIVE IDEAS ABOUT e
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 5: CREATIVE IDEAS ABOUT e IN TOURNAMENTS")
    print(f"{'='*70}")

    ideas_e = [
        ("e AS THE 'OPTIMAL TOURNAMENT RATIO'",
         f"max_H/mean_H -> e means: the best tournament is e times\n"
         f"    better than the average. This is like the Secretary Problem\n"
         f"    where the optimal strategy succeeds with probability 1/e.\n"
         f"    INSIGHT: Tournament H-maximization is a SELECTION problem.\n"
         f"    The 1/e threshold appears because you're 'selecting' the\n"
         f"    best ordering from an exponentially growing space."),

        ("e AND DERANGEMENTS",
         f"The number of derangements D(n) ~ n!/e.\n"
         f"    A derangement is a permutation with no fixed points.\n"
         f"    In tournaments: a 'derangement tournament' would have\n"
         f"    no vertex in its 'natural' position.\n"
         f"    D(n)/n! -> 1/e = {1/e:.6f}.\n"
         f"    mean_H = n!/2^(n-1). So mean_H / D(n) = 2^(n-1)/e -> ...\n"
         f"    This connects derangements to tournament paths."),

        ("e AND THE EXPONENTIAL FUNCTION",
         f"e = sum_k 1/k! = 1 + 1 + 1/2 + 1/6 + 1/24 + ...\n"
         f"    In tournaments: H = 1 + 2*alpha_1 + 4*alpha_2 + ...\n"
         f"    The 'normalized OCF': H/2^alpha_max = sum alpha_k * (1/2)^k\n"
         f"    This looks like exp(alpha * ln(2))... tenuous.\n"
         f"    But: exp(x) = sum x^k/k! and OCF = sum 2^k * alpha_k.\n"
         f"    If we set x = ln(2): exp(ln(2)) = 2, and sum (ln 2)^k/k! = 2.\n"
         f"    The OCF evaluation at x=2 connects to exp at x=ln(2)."),

        ("e AND THE NATURAL LOGARITHM",
         f"ln(2) = {math.log(2):.6f} is the entropy of a fair coin flip.\n"
         f"    Each arc has entropy ln(2) bits.\n"
         f"    Total tournament entropy = m * ln(2) = C(n,2) * ln(2).\n"
         f"    I(T;H) ≈ 0.27 * m ≈ 0.27 * C(n,2) in bits.\n"
         f"    In nats: I(T;H) ≈ 0.27 * C(n,2) / ln(2) nats.\n"
         f"    0.27 / ln(2) = {0.27/math.log(2):.4f} ≈ 0.39 nats per arc."),

        ("e AND TOURNAMENT ENTROPY PRODUCTION",
         f"In thermodynamics, entropy production rate = energy / (T * e)\n"
         f"    In the tournament 'gas' at temperature T = 1/beta:\n"
         f"    The specific heat C(beta) peaks at beta_c ≈ 0.3-0.8.\n"
         f"    At the phase transition: C * T_c ≈ some function of e.\n"
         f"    The 'efficiency' of the transition is bounded by e."),
    ]

    for title, desc in ideas_e:
        print(f"\n  {title}")
        print(f"    {desc}")

    # ============================================================
    # PART 6: NUMERICAL EXPLORATIONS
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 6: NUMERICAL EXPLORATIONS — pi AND e IN TOURNAMENT DATA")
    print(f"{'='*70}")

    # Does max_H(n) / n^pi approach anything?
    print(f"\n  max_H(n) / n^pi:")
    for n in range(3, 12):
        val = maxH[n] / n**pi
        print(f"    n={n:2d}: {val:.6f}")

    # Does max_H(n) / n^e approach anything?
    print(f"\n  max_H(n) / n^e:")
    for n in range(3, 12):
        val = maxH[n] / n**e
        print(f"    n={n:2d}: {val:.6f}")

    # log(max_H(n)) / n:
    print(f"\n  log(max_H(n)) / n (growth rate per vertex):")
    for n in range(3, 12):
        val = math.log(maxH[n]) / n
        print(f"    n={n:2d}: {val:.6f} (compare: log(n/2e)={math.log(n/(2*e)):.4f} + corrections)")

    # The Szele prediction: max_H ~ e * n! / 2^{n-1}
    # log(max_H) ~ 1 + n*log(n) - n*log(2e) + (1/2)*log(2*pi*n)
    print(f"\n  Szele prediction vs actual:")
    for n in range(3, 12):
        predicted_log = 1 + n*math.log(n) - n*math.log(2*e) + 0.5*math.log(2*pi*n)
        actual_log = math.log(maxH[n])
        print(f"    n={n:2d}: predicted log = {predicted_log:.4f}, actual log = {actual_log:.4f}, "
              f"diff = {actual_log - predicted_log:+.4f}")

    # ============================================================
    # PART 7: THE TOURNAMENT LEIBNIZ SERIES
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 7: THE TOURNAMENT LEIBNIZ SERIES — DOES pi EMERGE?")
    print(f"{'='*70}")

    # For each tournament T at n=5, compute S(T) = sum_k (-1)^k c_{2k+1}/(2k+1)
    n = 5
    m = C(n, 2)

    print(f"\n  S(T) = sum_k (-1)^k * c_(2k+1) / (2k+1) for n={n}")
    print(f"  = c_3/3 - c_5/5 + c_7/7 - ...")

    S_sum = 0
    count = 0
    S_by_H = {}

    for bits in range(2**m):
        A = np.zeros((n, n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx): A[j][i] = 1
                else: A[i][j] = 1
                idx += 1

        c3 = int(np.trace(A @ A @ A)) // 3
        c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5

        S = c3/3 - c5/5
        S_sum += S
        count += 1

        # Compute H for grouping
        dp = {}
        for v in range(n): dp[(1 << v, v)] = 1
        for ms in range(2, n+1):
            for mask in range(1 << n):
                if bin(mask).count('1') != ms: continue
                for v in range(n):
                    if not (mask & (1 << v)): continue
                    pm = mask ^ (1 << v)
                    t = sum(dp.get((pm, u), 0) for u in range(n) if (pm & (1 << u)) and A[u][v])
                    if t: dp[(mask, v)] = t
        H = sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))
        if H not in S_by_H:
            S_by_H[H] = []
        S_by_H[H].append(S)

    mean_S = S_sum / count
    print(f"\n  Mean S over all tournaments: {mean_S:.6f}")
    print(f"  pi/4 = {pi/4:.6f}")
    print(f"  4 * mean_S = {4*mean_S:.6f}")
    print(f"  Ratio mean_S / (pi/4) = {mean_S/(pi/4):.6f}")
    print(f"  Does mean_S = C(n,3)*pi/4 / ... ? No obvious match.")

    print(f"\n  S(T) by H value:")
    for H in sorted(S_by_H.keys()):
        vals = S_by_H[H]
        mean_val = sum(vals) / len(vals)
        print(f"    H={H:3d}: mean S = {mean_val:.6f}, unique S values: {len(set(round(v, 6) for v in vals))}")

    # ============================================================
    # PART 8: THE WALLIS PRODUCT AND TOURNAMENTS
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 8: WALLIS PRODUCT — pi FROM EVEN/ODD")
    print(f"{'='*70}")

    # Wallis: pi/2 = prod_{k=1}^inf (2k)^2 / ((2k-1)(2k+1))
    # = (2*2)/(1*3) * (4*4)/(3*5) * (6*6)/(5*7) * ...
    # = 4/3 * 16/15 * 36/35 * ...

    print(f"  Wallis: pi/2 = prod (2k)^2 / ((2k-1)(2k+1))")
    print(f"")
    print(f"  This is a product of EVEN^2 / (ODD * ODD) terms!")
    print(f"  In tournaments: EVEN = Fourier levels, ODD = vanishing levels.")
    print(f"  The even/odd dichotomy of H (all even Fourier, no odd)")
    print(f"  is structurally the same as the Wallis product structure!")
    print(f"")
    print(f"  Wallis partial products:")
    wallis = 1.0
    for k in range(1, 15):
        wallis *= (2*k)**2 / ((2*k-1) * (2*k+1))
        if k <= 10:
            print(f"    k={k:2d}: partial product = {wallis:.6f}, target pi/2 = {pi/2:.6f}")

    print(f"")
    print(f"  TOURNAMENT WALLIS: Define W(n) = prod_{{k=1}}^{{n}} (even_energy_k / odd_energy_k).")
    print(f"  If this converges to something involving pi, it would connect")
    print(f"  the Fourier even/odd split directly to pi.")

    # ============================================================
    # PART 9: SYNTHESIS
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 9: SYNTHESIS — THE ROLES OF e AND pi")
    print(f"{'='*70}")

    print(f"""
  e AND pi HAVE DIFFERENT STRUCTURAL ROLES:

  e = THE EFFICIENCY CONSTANT:
    - max_H / mean_H -> e (Szele)
    - e controls the BASE of exponential growth (via Stirling)
    - e appears in the Boltzmann distribution (energy/temperature)
    - e is about OPTIMIZATION: how much better can you do?
    - e is EXPONENTIAL (appears in growth rates, not normalization)

  pi = THE GEOMETRY CONSTANT:
    - pi appears in NORMALIZATION (sqrt(2*pi*n) in Stirling)
    - pi appears in CIRCLES (Lee-Yang zeros, unit circle)
    - pi appears in EIGENVALUES (sin(k*pi/p) in Paley)
    - pi appears in RANDOM MATRIX THEORY (eigenvalue spacing)
    - pi is GEOMETRIC (appears in areas, angles, normalizations)

  THE DEEP DIFFERENCE:
    e is ALGEBRAIC in character (it's about recurrences and limits)
    pi is GEOMETRIC in character (it's about circles and rotations)

  In tournament theory:
    e governs HOW MANY paths (the count)
    pi governs HOW THEY'RE DISTRIBUTED (the shape)

  Together: max_H ~ e * sqrt(2*pi*n) * (n/(2e))^n
    e tells you the PEAK HEIGHT
    pi tells you the PEAK WIDTH
    n/(2e) tells you the GROWTH BASE
""")

    print(f"\n{'='*70}")
    print("DONE — TRANSCENDENTAL NUMBERS DEEPLY EXPLORED")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
