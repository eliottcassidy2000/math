"""
pascal_fibonacci_alignment.py -- kind-pasteur-2026-03-14-S88
The 3-strand Pascal sequence aligned with Fibonacci.

THE USER'S INSIGHT:
Both sequences advance at the same rate (n incrementing by 1).
Despite 3 strands, the Pascal sequence has period-2 symmetry
(odd/even rows), giving a period of 2×3 = 6.

This 6-periodicity is FUNDAMENTAL. Let me explore why.

The Fibonacci sequence at position n: F_n
The 3-strand Pascal at position n: P_n = C(floor(2n/3)+1, ...)

Line them up and look for relationships.
"""

import sys, math
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding='utf-8')

def C(n, k):
    if k < 0 or k > n: return 0
    return math.comb(n, k)

def fib(n):
    """Fibonacci: F_0=0, F_1=1, F_2=1, F_3=2, ..."""
    if n <= 0: return 0
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a

def pascal_3strand(k):
    """The 3-strand Pascal sequence at position k."""
    n = k // 3
    r = k % 3
    if r == 0: return C(2*n+1, n)
    elif r == 1: return C(2*n+2, n)
    else: return C(2*n+2, n+1)

def main():
    print("=" * 70)
    print("PASCAL-FIBONACCI ALIGNMENT AND PERIOD-6 STRUCTURE")
    print("kind-pasteur-2026-03-14-S88")
    print("=" * 70)

    # ============================================================
    # PART 1: Side-by-side alignment
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 1: FIBONACCI AND 3-STRAND PASCAL SIDE BY SIDE")
    print(f"{'='*70}")

    print(f"\n  {'k':>3} {'F_k':>8} {'P_k':>8} {'P/F':>10} {'k%3':>4} {'k%6':>4} {'k%2':>4}")
    print(f"  {'-'*50}")

    for k in range(30):
        f = fib(k)
        p = pascal_3strand(k)
        ratio = p / f if f > 0 else float('inf')
        print(f"  {k:3d} {f:8d} {p:8d} {ratio:10.4f} {k%3:4d} {k%6:4d} {k%2:4d}")

    # ============================================================
    # PART 2: The period-6 structure
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 2: PERIOD-6 STRUCTURE — P/F RATIOS BY k MOD 6")
    print(f"{'='*70}")

    # Group ratios by k mod 6
    mod6_ratios = defaultdict(list)
    for k in range(6, 30):
        f = fib(k)
        p = pascal_3strand(k)
        if f > 0:
            mod6_ratios[k % 6].append(p / f)

    for r in range(6):
        ratios = mod6_ratios[r]
        if ratios:
            print(f"\n  k ≡ {r} (mod 6): P/F ratios = {[f'{x:.4f}' for x in ratios]}")
            if len(ratios) >= 2:
                # Are the ratios converging?
                diffs = [ratios[i+1] - ratios[i] for i in range(len(ratios)-1)]
                print(f"    Differences: {[f'{d:+.4f}' for d in diffs]}")
                print(f"    Converging? {all(abs(d) < abs(diffs[0]) for d in diffs[1:]) if len(diffs) > 1 else 'N/A'}")

    # ============================================================
    # PART 3: Growth rate comparison
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 3: GROWTH RATES — BOTH GROW EXPONENTIALLY")
    print(f"{'='*70}")

    phi = (1 + math.sqrt(5)) / 2
    print(f"\n  Fibonacci growth rate: phi = {phi:.6f}")
    print(f"  Pascal 3-strand growth rate per step: 4^(1/3) = {4**(1/3):.6f}")
    print(f"  Pascal 3-strand growth rate per 3 steps: 4 (doubling of central binomial)")
    print(f"")
    print(f"  Fibonacci per 3 steps: phi^3 = {phi**3:.6f}")
    print(f"  Pascal per 3 steps: 4")
    print(f"  Ratio: 4/phi^3 = {4/phi**3:.6f}")
    print(f"")
    print(f"  Fibonacci per 6 steps: phi^6 = {phi**6:.6f}")
    print(f"  Pascal per 6 steps: 4^2 = 16")
    print(f"  Ratio: 16/phi^6 = {16/phi**6:.6f}")

    # The ratio P_k / F_k grows because Pascal grows faster (4^{1/3} > phi ≈ 1.618)
    # Actually 4^{1/3} ≈ 1.587 < phi ≈ 1.618!
    # So Fibonacci grows FASTER per step!

    print(f"\n  4^(1/3) = {4**(1/3):.6f} vs phi = {phi:.6f}")
    print(f"  Fibonacci per step is FASTER! phi > 4^(1/3)")
    print(f"  So P/F ratio should DECREASE over time.")

    # Verify
    for k in [6, 12, 18, 24]:
        f = fib(k)
        p = pascal_3strand(k)
        print(f"  k={k}: P/F = {p/f:.6f}")

    # ============================================================
    # PART 4: The "6-cycle" in ratios and tournament connection
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 4: THE 6-CYCLE AND TOURNAMENT THEORY")
    print(f"{'='*70}")
    print(f"""
  WHY 6?
  The 3-strand Pascal has period 3 (from the 3 strands).
  The Fibonacci has period 2 (even/odd symmetry from n -> -n involution).
  Their LCM is 6.

  In tournament theory:
  - 3 = the fundamental odd cycle (3-cycle)
  - 2 = the arc binary choice (and the OCF fugacity)
  - 6 = LCM(3,2) = the "full period" of the combined structure

  This 6 appears as:
  - h(G2) = 6 (Coxeter number of the exceptional Lie algebra G2)
  - 6 = C(4,2) = number of arcs at n=4
  - 6 = 3! = number of labeled transitive tournaments on 3 vertices
  - The GS DOF formula has period related to floor((n-1)/2)

  THE DEEP REASON:
  3 and 2 are the GENERATORS of the theory.
  3-cycles generate tournament complexity.
  2-fold arc choice generates tournament variety.
  Their interaction has period LCM(3,2) = 6.
""")

    # ============================================================
    # PART 5: Log-ratio analysis — the precise alignment
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 5: LOG-RATIO ANALYSIS")
    print("  log(P_k) / log(F_k) should approach log(4^{1/3}) / log(phi)")
    print(f"{'='*70}")

    target_ratio = math.log(4**(1/3)) / math.log(phi)
    print(f"\n  Target: log(4^{{1/3}}) / log(phi) = {target_ratio:.6f}")

    for k in range(3, 30):
        f = fib(k)
        p = pascal_3strand(k)
        if f > 1 and p > 1:
            lr = math.log(p) / math.log(f)
            print(f"  k={k:2d}: log(P)/log(F) = {lr:.6f}, target = {target_ratio:.6f}, "
                  f"diff = {lr - target_ratio:+.6f}")

    # ============================================================
    # PART 6: The THREE subsequences of Fibonacci × THREE of Pascal
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 6: 3×2 = 6 CROSS-ALIGNMENT")
    print("  Fibonacci has 2 strands (even/odd index)")
    print("  Pascal has 3 strands (mod 3)")
    print("  Cross-product: 6 sub-sequences at positions k mod 6")
    print(f"{'='*70}")

    for r6 in range(6):
        fib_strand = "even" if r6 % 2 == 0 else "odd"
        pascal_strand = r6 % 3

        print(f"\n  k ≡ {r6} (mod 6): Fib strand = {fib_strand}, Pascal strand = {pascal_strand}")

        vals_f = []
        vals_p = []
        for k in range(r6, 30, 6):
            vals_f.append(fib(k))
            vals_p.append(pascal_3strand(k))

        print(f"    Fibonacci: {vals_f}")
        print(f"    Pascal:    {vals_p}")

        # Ratio
        ratios = [p/f if f > 0 else 0 for f, p in zip(vals_f, vals_p)]
        print(f"    P/F ratios: {[f'{r:.4f}' for r in ratios if r > 0]}")

    # ============================================================
    # PART 7: The tournament at n=6 — where 3 and 2 meet
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 7: n=6 — THE MEETING POINT OF 3 AND 2")
    print(f"{'='*70}")
    print(f"""
  n=6 is the FIRST n where ALL of these things happen:
  - H-landscape becomes MULTIMODAL (phase transition)
  - alpha_2 first becomes nonzero (disjoint cycle pairs appear)
  - The lex product T_3 lex T_2 achieves max_H = 45
  - The blue line skeleton becomes NON-bipartite
  - Blueself tilings first appear (even n required)

  n = 6 = LCM(3, 2) = 2 * 3

  The tournament at n=6 is where the 3-cycle structure (period 3)
  and the binary arc structure (period 2) SYNCHRONIZE for the first time.
  This synchronization creates the richest new phenomena.

  ANALOGY TO THE SEQUENCE:
  The 3-strand Pascal sequence at position 6 completes its first
  full 6-cycle (one complete period of both the 3-fold and 2-fold
  structures). Similarly, n=6 is the first tournament size where
  the 3-cycle and 2-fold structures both fully manifest.

  The values at k=0..5 (one full period):
    k=0: P=1, F=0   (genesis)
    k=1: P=1, F=1   (unit)
    k=2: P=2, F=1   (the 2 = binary choice)
    k=3: P=3, F=2   (the 3 = 3-cycle)
    k=4: P=4, F=3   (4 = 2^2 = arc count at n=3)
    k=5: P=6, F=5   (6 = LCM(2,3) = the meeting point)
""")

    # ============================================================
    # PART 8: Divisibility patterns mod 6
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 8: DIVISIBILITY PATTERNS")
    print(f"{'='*70}")

    print(f"\n  P_k mod 2 (period?): ", end="")
    mods2 = [pascal_3strand(k) % 2 for k in range(24)]
    print(mods2)

    print(f"  P_k mod 3 (period?): ", end="")
    mods3 = [pascal_3strand(k) % 3 for k in range(24)]
    print(mods3)

    print(f"  P_k mod 6 (period?): ", end="")
    mods6 = [pascal_3strand(k) % 6 for k in range(24)]
    print(mods6)

    print(f"\n  F_k mod 2 (period 3): ", end="")
    fmods2 = [fib(k) % 2 for k in range(24)]
    print(fmods2)

    print(f"  F_k mod 3 (period 8): ", end="")
    fmods3 = [fib(k) % 3 for k in range(24)]
    print(fmods3)

    # Pisano periods: period of F_k mod m
    # mod 2: period 3
    # mod 3: period 8
    # mod 6: period 24

    print(f"\n  Pisano period of Fibonacci mod 2: 3")
    print(f"  Pisano period of Fibonacci mod 3: 8")
    print(f"  Pisano period of Fibonacci mod 6: 24")
    print(f"  LCM(3, 8) = 24 ✓")

    # What about Pascal 3-strand mod 2?
    # P_k mod 2 = C(something, something) mod 2
    # By Lucas' theorem, C(n,k) mod 2 = product of C(n_i, k_i) mod 2
    # where n_i, k_i are binary digits.

    # Period of P_k mod 2:
    print(f"\n  Looking for period of P_k mod 2:")
    for period in range(1, 25):
        is_period = all(mods2[i] == mods2[i + period] for i in range(min(24 - period, 12)))
        if is_period:
            print(f"    Period {period}: YES")
            break

    print(f"\n  Looking for period of P_k mod 3:")
    for period in range(1, 25):
        is_period = all(mods3[i] == mods3[i + period] for i in range(min(24 - period, 12)))
        if is_period:
            print(f"    Period {period}: YES")
            break

    # ============================================================
    # PART 9: The "Fibonacci-Pascal resonance"
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 9: THE FIBONACCI-PASCAL RESONANCE")
    print(f"{'='*70}")

    # When do F_k and P_k share a common factor > 1?
    print(f"\n  GCD(F_k, P_k) for k=1..29:")
    for k in range(1, 30):
        f = fib(k)
        p = pascal_3strand(k)
        g = math.gcd(f, p)
        if g > 1:
            print(f"    k={k:2d}: F={f}, P={p}, GCD={g}")

    # When does P_k / F_k happen to be an integer?
    print(f"\n  P_k / F_k integer?")
    for k in range(1, 30):
        f = fib(k)
        p = pascal_3strand(k)
        if f > 0 and p % f == 0:
            print(f"    k={k:2d}: P/F = {p//f}")

    # ============================================================
    # PART 10: The independence polynomial at the 6-cycle
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 10: INDEPENDENCE POLYNOMIAL OF THE 6-CYCLE")
    print("  The 6-cycle C_6 is the simplest graph with period 6.")
    print("  I(C_6, x) = ? and how does it relate to tournaments?")
    print(f"{'='*70}")

    # I(C_n, x) for cycle graphs:
    # I(C_n, x) = L_n(x) + 2*(-x)^n / (1+x)... actually:
    # I(C_n, x) = F_n(x) + F_{n-2}(x) where F_k is Fibonacci polynomial
    # Actually: I(C_n, x) = sum_{k=0}^{floor(n/2)} (n/(n-k)) * C(n-k, k) * x^k

    def indpoly_cycle(n, x):
        """Independence polynomial of cycle C_n at x."""
        total = 0
        for k in range(n // 2 + 1):
            if n - k > 0:
                coeff = n * C(n - k, k) // (n - k)
            else:
                coeff = 1
            total += coeff * x**k
        return total

    for n in [3, 4, 5, 6, 7, 8]:
        I1 = indpoly_cycle(n, 1)
        I2 = indpoly_cycle(n, 2)
        print(f"  I(C_{n}, 1) = {I1} (Lucas number L_{n})")
        print(f"  I(C_{n}, 2) = {I2}")

    # I(C_6, 2) = ?
    I_C6_2 = indpoly_cycle(6, 2)
    print(f"\n  I(C_6, 2) = {I_C6_2}")
    print(f"  Is this a tournament H value? ", end="")

    # C_6 as Omega: would need a tournament whose conflict graph IS C_6
    # C_6 has 6 vertices (= 6 odd cycles in the tournament)
    # Each pair of adjacent cycles shares a vertex
    # Each pair of non-adjacent cycles is vertex-disjoint

    # If Omega = C_6: H = I(C_6, 2) = {I_C6_2}
    # Is this achievable? Need 6 odd cycles arranged in a hexagonal
    # conflict pattern. This would require n >= 6 at minimum.

    print(f"H = {I_C6_2} would require Omega = C_6")
    print(f"  This needs 6 odd cycles in a hexagonal conflict pattern.")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
