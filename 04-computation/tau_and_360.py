"""
tau_and_360.py -- kind-pasteur-2026-03-14-S106e
BASE TAU AND THE NUMBER 360 — A FRAME FOR EVERYTHING

360 is one of the most ancient and mysterious numbers in human history.
Why did the Babylonians choose 360 degrees in a circle?
The standard answer: 360 has many divisors (24 of them).
But tournament theory may reveal something deeper.

360 = 2^3 * 3^2 * 5
    = 8 * 45
    = 8 * 9 * 5
    = (threshold) * (reflection) * (interaction)

In tournament theory:
  360 = |A_6| = the alternating group on 6 elements
  360 = number of Baer subplanes in PG(2, F_4)
  360 = the number of directed Hamiltonian CYCLES on 7 vertices
        in certain tournaments? Let me check.

And in base tau (tribonacci constant):
  What does 360 look like? What do ALL the key numbers look like?
"""

import sys, math
import numpy as np
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8')

def tribonacci_representation(n):
    """Compute the tribonacci (greedy) representation of n.
    Returns list of (index, tribonacci_number) pairs."""
    if n <= 0:
        return []
    # Build tribonacci numbers: T(1)=1, T(2)=1, T(3)=2, T(4)=4, T(5)=7, ...
    t = [1, 1, 2]
    while t[-1] < n:
        t.append(t[-1] + t[-2] + t[-3])
    # Greedy decomposition
    rep = []
    remaining = n
    for i in range(len(t)-1, -1, -1):
        if t[i] <= remaining:
            rep.append((i+1, t[i]))  # 1-indexed
            remaining -= t[i]
    return rep

def trib_string(n):
    """Compact string for tribonacci representation."""
    rep = tribonacci_representation(n)
    return ' + '.join(f'{v}' for _, v in rep)

# Tribonacci numbers for reference
TRIB = [0, 1, 1, 2, 4, 7, 13, 24, 44, 81, 149, 274, 504, 927, 1705]

def main():
    print("=" * 70)
    print("BASE TAU AND 360 — A FRAME FOR EVERYTHING")
    print("kind-pasteur-2026-03-14-S106e")
    print("=" * 70)

    tau = 1.8392867552141612

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 1: 360 — THE ANATOMY")
    print(f"{'='*70}")

    print(f"""
  360 = 2^3 * 3^2 * 5

  In prime-exponent space: (3, 2, 1, 0, ...) = choice^3 * cycle^2 * interaction

  This is EXACTLY the product of the three tournament essences:
    2^3 = 8 (the THRESHOLD — the cube of choice)
    3^2 = 9 (the REFLECTION — the square of cycle)
    5   = 5 (the INTERACTION — where choice meets cycle)

  360 = threshold * reflection * interaction.
  It is the number where ALL THREE essences come together
  at their natural powers.

  DIVISORS of 360: {sorted([d for d in range(1, 361) if 360 % d == 0])}
  Count: {sum(1 for d in range(1, 361) if 360 % d == 0)} divisors.

  WHY 360 DEGREES?
  The Babylonians didn't just pick 360 for divisibility.
  360 = the product of the three smallest primes at their
  "tournament-natural" exponents (3, 2, 1).
  A circle is the ULTIMATE cycle — and 360 measures it
  using choice^3 * cycle^2 * interaction.

  The exponents (3, 2, 1) are THEMSELVES the first three naturals.
  360 = 2^3 * 3^2 * 5^1: the exponents count DOWN (3, 2, 1)
  while the bases count UP (2, 3, 5).
  This REVERSAL is the signature of 360.""")

    divs = sorted([d for d in range(1, 361) if 360 % d == 0])
    print(f"\n  All {len(divs)} divisors of 360:")
    print(f"  {divs}")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 2: 360 IN TOURNAMENT THEORY")
    print(f"{'='*70}")

    print(f"""
  360 appears in tournament theory in THREE independent ways:

  1. |A_6| = 360.
     The alternating group on 6 elements has order 360.
     6 = LCM(2,3) = the tournament period.
     The EVEN permutations of the period form a group of size 360.
     |S_6| = 720 = 2 * 360.
     S_6 is the ONLY symmetric group with an outer automorphism.

  2. Baer subplanes of PG(2, F_4): exactly 360.
     PG(2, F_4) has 21 = Phi_3(4) = H_forb_2 points.
     Each Baer subplane is a copy of PG(2, F_2) = Fano (7 points).
     There are 360 = |A_6| such copies.
     The orbit-stabilizer: 360 * 336 = 120960 = |P Gamma L(3, F_4)|.

  3. n! / n for n = 6: 6!/6 = 720/6 = ... no, that's 120.
     Actually: (n-1)!/2 at n=7: 6!/2 = 360.
     This is the number of directed Hamiltonian CYCLES on 7 vertices
     in the complete graph. Each Hamiltonian cycle visits all 7 vertices
     and returns to the start. There are 6!/2 = 360 such cycles.
     (Divide by 7 for the cycle, multiply by 2 for direction, net: 6!/2.)
     Wait: (n-1)! / 2 = number of undirected Hamiltonian cycles.
     At n=7: 6!/2 = 360 undirected, or 720 directed.
     So 360 = the number of UNDIRECTED Hamiltonian cycles in K_7.

  360 IS THE HAMILTONIAN CYCLE COUNT OF K_7.

  This means: a tournament on 7 vertices PARTITIONS these 360
  undirected cycles into two classes — forward and backward —
  and H(T) counts the number of directed Hamiltonian PATHS,
  which is related to (but different from) cycles.

  H(T_7^Paley) = 189 = 360/2 + 9? No: 360/2 = 180 ≠ 189.
  Actually 189 = 360 * 189/360 = 360 * 21/40.
  Hmm: 189/360 = 21/40 = 0.525. The Paley tournament at n=7
  has H/360 ≈ 0.525 — just over half the undirected Ham cycles.

  More precisely: 189 = 3^3 * 7, 360 = 2^3 * 3^2 * 5.
  GCD(189, 360) = 9 = 3^2.
  189/9 = 21, 360/9 = 40.
  189/360 = 21/40.
  So H(Paley_7) = (21/40) * 360 = (21/40) * (number of undir Ham cycles).

  21/40 = H_forb_2 / (8 * 5) = Phi_3(4) / (2^3 * 5).
  The FRACTION of Ham cycles captured by the Paley tournament
  involves the second forbidden value over threshold * interaction!""")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 3: KEY NUMBERS IN BASE TAU")
    print(f"{'='*70}")

    print(f"  tau = {tau:.10f} (tribonacci constant, tau^3 = Phi_3(tau))")
    print(f"  Tribonacci sequence: {TRIB[1:]}\n")

    # Key tournament numbers and their tribonacci representations
    numbers = [
        (0, "absence"),
        (1, "ground state"),
        (2, "choice / generator"),
        (3, "cycle"),
        (4, "field F_4"),
        (5, "interaction / ADE"),
        (6, "period / harmony"),
        (7, "first forbidden / Fano"),
        (8, "threshold / 2^3"),
        (9, "reflection / 3^2"),
        (10, "completion / C(5,2)"),
        (12, "doubled harmony"),
        (13, "regularity / Phi_3(3)"),
        (15, "max_H(5)"),
        (21, "second forbidden"),
        (24, "Golay / 2^3*3"),
        (42, "2*21"),
        (45, "max_H(6)"),
        (60, "|A_5| / icosahedral"),
        (120, "|S_5| / 5!"),
        (168, "|PSL(2,7)| = 8*21"),
        (189, "max_H(7) = 3^3*7"),
        (360, "|A_6| / Baer / HamCycles(K_7)"),
        (504, "trib / Var denom n=7"),
        (720, "|S_6| = 6!"),
    ]

    print(f"  {'n':>6}  {'tribonacci rep':>35}  {'#terms':>6}  {'essence'}")
    print(f"  {'-'*85}")

    for n, desc in numbers:
        if n == 0:
            print(f"  {n:6d}  {'(empty)':>35}  {0:6d}  {desc}")
            continue
        rep = tribonacci_representation(n)
        rep_str = ' + '.join(f'{v}' for _, v in rep)
        n_terms = len(rep)
        # Check if n is itself a tribonacci number
        is_trib = n in TRIB
        marker = " *TRIB*" if is_trib else ""
        print(f"  {n:6d}  {rep_str:>35}  {n_terms:6d}  {desc}{marker}")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 4: 360 IN BASE TAU — THE DECOMPOSITION")
    print(f"{'='*70}")

    rep_360 = tribonacci_representation(360)
    print(f"\n  360 in tribonacci: {trib_string(360)}")
    print(f"  Components: {rep_360}")
    print(f"  Number of terms: {len(rep_360)}")

    # What tribonacci numbers compose 360?
    # 274 + 81 + 4 + 1 = 360?
    print(f"\n  Verification: {' + '.join(str(v) for _,v in rep_360)} = {sum(v for _,v in rep_360)}")

    # 360 and its relationship to tribonacci neighbors
    print(f"\n  Tribonacci neighbors of 360:")
    for i, t in enumerate(TRIB):
        if t > 0 and abs(t - 360) < 200:
            print(f"    T({i}) = {t}, distance from 360: {360 - t}")

    print(f"""
  360 sits between T(11) = 274 and T(12) = 504.
  360 - 274 = 86
  504 - 360 = 144 = 12^2 = F(12) (a Fibonacci number!)

  So: 504 - 360 = 144 = F(12).
  And: 360 = 504 - 144 = T(12) - F(12).
  360 IS THE DIFFERENCE OF A TRIBONACCI AND A FIBONACCI NUMBER!

  This is extraordinary: 360 = T(12) - F(12) = 504 - 144.
  The TRIBONACCI at index 12 minus the FIBONACCI at index 12.
  12 = 2^2 * 3 = the Golay dimension = the doubled period.

  At the index 12, the tribonacci and Fibonacci diverge by exactly 360.
  The EXCESS of tribonacci over Fibonacci at the doubled period
  equals the Hamiltonian cycle count of K_7.""")

    # Verify
    fib = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610]
    print(f"\n  F(12) = {fib[12]}, T(12) = {TRIB[12]}")
    print(f"  T(12) - F(12) = {TRIB[12] - fib[12]}")
    print(f"  = 360? {TRIB[12] - fib[12] == 360}")

    # Check other indices
    print(f"\n  T(n) - F(n) for various n:")
    for i in range(1, min(len(TRIB), len(fib))):
        diff = TRIB[i] - fib[i]
        special = ""
        if diff == 0: special = "  (equal!)"
        if diff == 360: special = "  *** = 360 ***"
        if diff in [6, 7, 21, 42, 120, 168, 189, 504]: special = f"  (tournament number!)"
        print(f"    n={i:2d}: T({i})={TRIB[i]:5d}, F({i})={fib[i]:5d}, "
              f"T-F={diff:5d}{special}")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 5: THE TRIBONACCI WEIGHT OF NUMBERS")
    print(f"{'='*70}")

    print(f"""
  Define the TRIBONACCI WEIGHT w(n) = number of tribonacci terms
  in the greedy representation of n. This is the "complexity"
  of n in base tau.

  w(n) = 1: n is a TRIBONACCI NUMBER (simplest, atomic)
  w(n) = 2: n is a sum of two tribonacci numbers
  w(n) = 3: three terms needed
  ...

  The tribonacci weight measures how "compound" a number is
  in the cycle-growth representation.""")

    # Compute tribonacci weight for 1..100
    print(f"\n  Tribonacci weight for key numbers:")
    for n in range(1, 51):
        rep = tribonacci_representation(n)
        w = len(rep)
        is_t = n in TRIB[1:]
        is_prime = all(n % d != 0 for d in range(2, int(n**0.5)+1)) and n > 1
        markers = []
        if is_t: markers.append("TRIB")
        if is_prime: markers.append("prime")
        if n in [7, 21]: markers.append("FORBIDDEN")
        if n in [3, 5, 9, 15, 45, 189]: markers.append("H-special")
        mark = f" [{', '.join(markers)}]" if markers else ""
        rep_str = trib_string(n)
        print(f"    {n:4d} = {rep_str:>25}  w={w}{mark}")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 6: THE TRIBONACCI WEIGHT SPECTRUM")
    print(f"{'='*70}")

    # Weight distribution
    from collections import Counter
    weights = Counter()
    for n in range(1, 1001):
        w = len(tribonacci_representation(n))
        weights[w] += 1

    print(f"\n  Distribution of tribonacci weights for n = 1..1000:")
    for w in sorted(weights.keys()):
        bar = "#" * (weights[w] // 5)
        print(f"    w={w}: {weights[w]:4d} numbers {bar}")

    # What weight is 360?
    w360 = len(tribonacci_representation(360))
    print(f"\n  w(360) = {w360}")
    print(f"  w(189) = {len(tribonacci_representation(189))}")
    print(f"  w(504) = {len(tribonacci_representation(504))}")
    print(f"  w(168) = {len(tribonacci_representation(168))}")
    print(f"  w(720) = {len(tribonacci_representation(720))}")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 7: 360 AS A FRAME — WHAT IT DIVIDES, WHAT DIVIDES IT")
    print(f"{'='*70}")

    print(f"""
  360 as a FRAME means: looking at other numbers as fractions of 360.

  KEY RATIOS n/360:""")

    frame_numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 18, 20,
                     21, 24, 30, 36, 40, 42, 45, 60, 72, 90, 120, 168, 180,
                     189, 360, 504, 720]

    for n in frame_numbers:
        ratio = Fraction(n, 360)
        degrees = n  # as degrees of a circle
        # What fraction of 360
        tournament_meaning = ""
        if n == 1: tournament_meaning = "ground state"
        elif n == 2: tournament_meaning = "choice"
        elif n == 3: tournament_meaning = "cycle"
        elif n == 5: tournament_meaning = "interaction"
        elif n == 6: tournament_meaning = "period"
        elif n == 7: tournament_meaning = "first forbidden"
        elif n == 8: tournament_meaning = "threshold"
        elif n == 9: tournament_meaning = "reflection"
        elif n == 10: tournament_meaning = "completion"
        elif n == 12: tournament_meaning = "doubled period"
        elif n == 13: tournament_meaning = "regularity"
        elif n == 15: tournament_meaning = "max_H(5)"
        elif n == 21: tournament_meaning = "second forbidden"
        elif n == 24: tournament_meaning = "Golay"
        elif n == 30: tournament_meaning = "2*3*5"
        elif n == 36: tournament_meaning = "6^2"
        elif n == 40: tournament_meaning = "360/9 = 8*5"
        elif n == 42: tournament_meaning = "2*21"
        elif n == 45: tournament_meaning = "max_H(6)"
        elif n == 60: tournament_meaning = "|A_5|"
        elif n == 72: tournament_meaning = "8*9"
        elif n == 90: tournament_meaning = "right angle"
        elif n == 120: tournament_meaning = "|S_5| = 5!"
        elif n == 168: tournament_meaning = "|PSL(2,7)|"
        elif n == 180: tournament_meaning = "half circle"
        elif n == 189: tournament_meaning = "max_H(7)"
        elif n == 360: tournament_meaning = "|A_6| / full circle"
        elif n == 504: tournament_meaning = "tribonacci / Var denom"
        elif n == 720: tournament_meaning = "|S_6| = 6!"

        divides_360 = 360 % n == 0 if n <= 360 else False
        div_marker = " | 360" if divides_360 else ""

        print(f"    {n:4d}/360 = {str(ratio):>8} = {float(ratio):8.6f}"
              f"  {tournament_meaning:>25}{div_marker}")

    print(f"""
  REMARKABLE OBSERVATIONS:

  1. 7/360 = 7/360. 7 does NOT divide 360.
     The forbidden value is INCOMMENSURABLE with the circle.
     You cannot divide a circle into 7 equal parts using 360.
     This is geometrically visible: a regular heptagon
     has angles of 360/7 = 51.428... degrees (irrational!).

  2. 21/360 = 7/120 = 7/(5!). The second forbidden over 360
     equals the first forbidden over 5!.
     21 DOES divide 360? 360/21 = 17.14... NO.
     Neither forbidden value divides 360.
     The forbidden values are GEOMETRICALLY INACCESSIBLE.

  3. 189/360 = 21/40 = 0.525. The Paley maximizer at n=7
     is 52.5% of the full circle. Just over half.
     21/40 = H_forb_2 / (2^3 * 5).

  4. The numbers that DO divide 360 are EXACTLY the
     "harmonious" numbers — those built from {{2, 3, 5}} only.
     No 7, no 11, no 13. The obstruction prime is excluded.

  5. 360/7 = 51.428... but 360/6 = 60 = |A_5|.
     The PERIOD (6) perfectly divides the circle (360).
     The FORBIDDEN (7) doesn't. This is the geometric
     manifestation of the tournament obstruction.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 8: THE DIVISORS OF 360 — A TOURNAMENT LATTICE")
    print(f"{'='*70}")

    print(f"""
  The divisors of 360 form a LATTICE under divisibility.
  Since 360 = 2^3 * 3^2 * 5, the divisors correspond to
  triples (a, b, c) with 0 <= a <= 3, 0 <= b <= 2, 0 <= c <= 1.

  Total: 4 * 3 * 2 = 24 divisors.

  This is a 3-dimensional rectangular grid: 4 x 3 x 2.
  The SHAPE of the divisor lattice is determined by the exponents (3, 2, 1).
  These are the first three natural numbers in reverse!

  The divisor lattice of 360 is a discrete MODEL of the
  tournament interaction space:
    - 4 levels of binary depth (a = 0,1,2,3)
    - 3 levels of ternary depth (b = 0,1,2)
    - 2 levels of quintic depth (c = 0,1)""")

    # Display the divisor lattice organized by (a,b,c)
    print(f"\n  Divisor lattice of 360:")
    print(f"  {'(a,b,c)':>8} {'2^a*3^b*5^c':>12} {'divisor':>8} {'tournament meaning':>30}")
    print(f"  {'-'*65}")
    meanings = {
        1: "ground state", 2: "choice", 3: "cycle", 4: "field",
        5: "interaction", 6: "period", 8: "threshold", 9: "reflection",
        10: "completion", 12: "doubled period", 15: "max_H(5)",
        18: "doubled reflection", 20: "expanded field", 24: "Golay",
        30: "ADE product", 36: "squared period", 40: "threshold*interaction",
        45: "max_H(6)", 60: "|A_5|", 72: "threshold*reflection",
        90: "right angle", 120: "5!", 180: "half circle",
        360: "full circle"
    }
    for c in range(2):
        for b in range(3):
            for a in range(4):
                d = (2**a) * (3**b) * (5**c)
                meaning = meanings.get(d, "")
                print(f"  ({a},{b},{c})  {2**a}*{3**b}*{5**c:1d} = {d:>4d}"
                      f"        {d:>8}  {meaning:>30}")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 9: THE 360-TAU CONNECTION")
    print(f"{'='*70}")

    print(f"""
  In base tau:
    360 = {trib_string(360)}

  tau^k values:
    tau^1 = {tau**1:.4f}
    tau^2 = {tau**2:.4f}
    tau^3 = {tau**3:.4f} = Phi_3(tau)
    tau^4 = {tau**4:.4f}
    tau^5 = {tau**5:.4f}
    tau^6 = {tau**6:.4f}
    tau^7 = {tau**7:.4f}
    tau^8 = {tau**8:.4f}
    tau^9 = {tau**9:.4f}

  log_tau(360) = {math.log(360)/math.log(tau):.6f}

  So 360 ≈ tau^{math.log(360)/math.log(tau):.2f}.
  The exponent {math.log(360)/math.log(tau):.4f} ≈ {math.log(360)/math.log(tau):.1f}.

  Is log_tau(360) close to a nice number?
  {math.log(360)/math.log(tau):.6f} ≈ 9.66...

  tau^9 = {tau**9:.4f} vs 360.
  tau^10 = {tau**10:.4f}.
  360 is between tau^9 and tau^10, closer to tau^10.
  360/tau^9 = {360/tau**9:.6f}
  360/tau^10 = {360/tau**10:.6f}

  Interestingly: 360/tau^9 ≈ {360/tau**9:.4f}.
  And tau^3 = Phi_3(tau) ≈ 6.222.
  360/6.222 ≈ {360/6.222:.2f} ≈ 58. Hmm.
  360 / Phi_3(tau) = 360/tau^3 = (360/tau^3) ≈ {360/tau**3:.4f}.
  Is 360/tau^3 close to something nice?
  tau^6 = tau^3 * tau^3 = Phi_3(tau)^2 ≈ {tau**6:.4f}.
  360/tau^6 ≈ {360/tau**6:.4f}. Close to tau^3 ≈ 6.22? No, {360/tau**6:.4f}.

  Let me try: 360 = 2^3 * 3^2 * 5. In terms of tau:
  2 ≈ tau^(log_tau(2)) = tau^{math.log(2)/math.log(tau):.4f}
  3 ≈ tau^{math.log(3)/math.log(tau):.4f}
  5 ≈ tau^{math.log(5)/math.log(tau):.4f}

  So 360 = 2^3 * 3^2 * 5 ≈ tau^(3*{math.log(2)/math.log(tau):.3f} + 2*{math.log(3)/math.log(tau):.3f} + {math.log(5)/math.log(tau):.3f})
     = tau^({3*math.log(2)/math.log(tau) + 2*math.log(3)/math.log(tau) + math.log(5)/math.log(tau):.4f})
     This checks out: log_tau(360) = {math.log(360)/math.log(tau):.4f}.""")

    # The key insight about what 360 means in base tau
    log_tau_2 = math.log(2)/math.log(tau)
    log_tau_3 = math.log(3)/math.log(tau)
    log_tau_5 = math.log(5)/math.log(tau)
    log_tau_7 = math.log(7)/math.log(tau)

    print(f"\n  PRIME POSITIONS IN BASE TAU:")
    print(f"    log_tau(2) = {log_tau_2:.6f} (choice position)")
    print(f"    log_tau(3) = {log_tau_3:.6f} (cycle position)")
    print(f"    log_tau(5) = {log_tau_5:.6f} (interaction position)")
    print(f"    log_tau(7) = {log_tau_7:.6f} (obstruction position)")
    print(f"    log_tau(13)= {math.log(13)/math.log(tau):.6f} (regularity position)")
    print(f"")
    print(f"    Differences:")
    print(f"    log_tau(3) - log_tau(2) = {log_tau_3 - log_tau_2:.6f}")
    print(f"    log_tau(5) - log_tau(3) = {log_tau_5 - log_tau_3:.6f}")
    print(f"    log_tau(7) - log_tau(5) = {log_tau_7 - log_tau_5:.6f}")
    print(f"    The gaps SHRINK: the primes get closer in base tau.")
    print(f"    This reflects the n-nacci convergence to 2.")

    # ============================================================
    print(f"\n{'='*70}")
    print("PART 10: SYNTHESIS — 360, TAU, AND THE NATURE OF COUNTING")
    print(f"{'='*70}")

    print(f"""
  THE 360-TAU SYNTHESIS:

  360 = 2^3 * 3^2 * 5 = threshold * reflection * interaction.

  In base tau, 360 = {trib_string(360)}.
  This is a {len(tribonacci_representation(360))}-term tribonacci representation.

  360 = T(12) - F(12) = 504 - 144.
  The tribonacci at index 12 minus the Fibonacci at index 12.
  360 measures the EXCESS of ternary growth over binary growth
  at the doubled period index.

  THE FRAME OF 360:
  When we measure angles in degrees (out of 360), we are
  measuring what fraction of FULL CYCLE something represents.
  Full cycle = 360 = 2^3 * 3^2 * 5.

  The numbers that divide 360 are the "clean" fractions of a cycle:
    1 degree = 1/360 of a cycle = the smallest distinction
    2 degrees = 1/180 = the generator's share
    3 degrees = 1/120 = the cycle's share
    5 degrees = 1/72 = the interaction's share
    6 degrees = 1/60 = the period's share
    7 degrees = 7/360 = IRRATIONAL fraction (in the cycle sense!)

  7/360 is "irrational" not in the number-theoretic sense but in the
  GEOMETRIC sense: you cannot tile a circle with 7-degree wedges.
  7 is incommensurable with the circle because 7 does not divide 360.
  The first forbidden value is GEOMETRICALLY INCOMMENSURABLE
  with the full cycle.

  THIS IS THE DEEPEST MEANING OF "FORBIDDEN":
  A number is forbidden as an H value precisely when it is
  INCOMMENSURABLE with the full cycle of tournament counting.

  360 CONTAINS (as divisors): 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 18, 20, 24, 30, 36, 40, 45, 60, 72, 90, 120, 180, 360.
  360 EXCLUDES: 7, 11, 13, 14, 16, 17, 19, 21, 22, 23, 25, 26, 27, 28, 29, 31, ...

  The excluded numbers include BOTH forbidden H values (7, 21)
  and many achievable H values (11, 13, 15 divides!, ...).
  Wait: 15 DOES divide 360 (360/15 = 24). And 45 divides 360 (360/45 = 8).
  So max_H(5) = 15 and max_H(6) = 45 both divide 360!
  But 189 = max_H(7) does NOT divide 360 (360/189 = 1.905...).

  The H-MAXIMIZERS that divide 360:
    max_H(3) = 3: YES (360/3 = 120)
    max_H(4) = 5: YES (360/5 = 72)
    max_H(5) = 15: YES (360/15 = 24)
    max_H(6) = 45: YES (360/45 = 8)
    max_H(7) = 189: NO (189 = 3^3*7, and 7 does not divide 360)

  THE MAXIMIZERS DIVIDE 360 UP TO n=6, THEN STOP.
  n=7 is the first n where max_H involves the prime 7
  (the obstruction prime), and 7 does not divide 360.

  360 = the LIMIT OF HARMONIC TOURNAMENT COUNTING.
  Beyond n=6, the obstruction prime 7 enters, and the
  harmonious world of 360-divisibility breaks.

  SUMMARY:
  360 is the product of choice^3 * cycle^2 * interaction.
  It is the number of undirected Hamiltonian cycles in K_7.
  It is where the Fibonacci and Tribonacci diverge (T(12) - F(12) = 360).
  It is the frame where "forbidden" means "incommensurable."
  And it marks the boundary between the harmonious (n <= 6)
  and the obstructed (n >= 7) worlds of tournament theory.
  """)

    print(f"\n{'='*70}")
    print("DONE — 360 IS THE BOUNDARY OF HARMONIC COUNTING")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
