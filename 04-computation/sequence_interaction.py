"""
sequence_interaction.py -- kind-pasteur-2026-03-14-S90 (continued)
How do a(n) = C(n,1)+C(n,3), Fibonacci, triangular, and 3-strand Pascal
INTERACT when aligned at the same index?

Focus on: GCD patterns, divisibility, ratios, and where they predict
each other's values.
"""

import sys, math
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding='utf-8')

def C(n, k):
    if k < 0 or k > n: return 0
    return math.comb(n, k)

def fib(n):
    if n <= 0: return 0
    a, b = 0, 1
    for _ in range(n): a, b = b, a + b
    return a

def tri(n): return n*(n+1)//2
def a_seq(n): return n + n*(n-1)*(n-2)//6 if n >= 0 else n + n*(n-1)*(n-2)//6

def pascal_3strand(k):
    n = k // 3
    r = k % 3
    if r == 0: return C(2*n+1, n)
    elif r == 1: return C(2*n+2, n)
    else: return C(2*n+2, n+1)

maxH = [0, 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095]

def main():
    print("=" * 70)
    print("SEQUENCE INTERACTION — HOW THEY MOVE TOGETHER")
    print("kind-pasteur-2026-03-14-S90 (continued)")
    print("=" * 70)

    # ============================================================
    # INTERACTION 1: GCD patterns
    # ============================================================
    print(f"\n{'='*70}")
    print("INTERACTION 1: GCD PATTERNS")
    print(f"{'='*70}")

    print(f"\n  {'n':>3} {'a(n)':>6} {'F(n)':>6} {'T(n)':>6} {'GCD(a,F)':>8} {'GCD(a,T)':>8} {'GCD(F,T)':>8}")
    for n in range(1, 16):
        an = a_seq(n)
        fn = fib(n)
        tn = tri(n)
        print(f"  {n:3d} {an:6d} {fn:6d} {tn:6d} "
              f"{math.gcd(an, fn):8d} {math.gcd(an, tn):8d} {math.gcd(fn, tn):8d}")

    # ============================================================
    # INTERACTION 2: Divisibility — when does one divide another?
    # ============================================================
    print(f"\n{'='*70}")
    print("INTERACTION 2: DIVISIBILITY PATTERNS")
    print(f"{'='*70}")

    print(f"\n  F(n) | a(n)? (Fibonacci divides a(n)):")
    for n in range(1, 20):
        an = a_seq(n)
        fn = fib(n)
        if fn > 0 and an % fn == 0:
            print(f"    n={n:2d}: a(n)/F(n) = {an}/{fn} = {an//fn}")

    print(f"\n  T(n) | a(n)? (Triangular divides a(n)):")
    for n in range(1, 20):
        an = a_seq(n)
        tn = tri(n)
        if tn > 0 and an % tn == 0:
            print(f"    n={n:2d}: a(n)/T(n) = {an}/{tn} = {an//tn}")

    print(f"\n  a(n) | T(n)? (a(n) divides triangular):")
    for n in range(1, 20):
        an = a_seq(n)
        tn = tri(n)
        if an > 0 and tn % an == 0:
            print(f"    n={n:2d}: T(n)/a(n) = {tn}/{an} = {tn//an}")

    # ============================================================
    # INTERACTION 3: The "prediction" game
    # ============================================================
    print(f"\n{'='*70}")
    print("INTERACTION 3: CAN ONE SEQUENCE PREDICT ANOTHER?")
    print(f"{'='*70}")

    # Does a(n) predict max_H(n)?
    print(f"\n  max_H(n) / a(n):")
    for n in range(3, min(12, len(maxH))):
        an = a_seq(n)
        mh = maxH[n]
        ratio = mh / an if an > 0 else 0
        print(f"    n={n:2d}: maxH={mh:6d}, a(n)={an:6d}, ratio={ratio:.4f}")

    # a(n) grows as n^3/6, max_H grows super-exponentially.
    # The ratio grows rapidly.
    # But at n=5: ratio = 1.0000 (exact match!)

    # Does a(n) * something = max_H(n)?
    print(f"\n  max_H(n) / a(n) pattern:")
    print(f"    n=3: 3/4 = 0.75")
    print(f"    n=4: 5/8 = 0.625")
    print(f"    n=5: 15/15 = 1.000  <-- EXACT!")
    print(f"    n=6: 45/26 = 1.731")
    print(f"    n=7: 189/42 = 4.500")
    print(f"    n=8: 661/64 = 10.328")
    print(f"    n=9: 3357/93 = 36.097")

    # The ratio maxH/a grows roughly as max_H grows exponentially.
    # But note: 189/42 = 4.5 = 9/2. And 45/26 is not clean.

    # ============================================================
    # INTERACTION 4: The "even-odd interleaving" of a(n)
    # ============================================================
    print(f"\n{'='*70}")
    print("INTERACTION 4: EVEN-ODD INTERLEAVING OF a(n)")
    print("  a(n) has a beautiful symmetry: a(n) + a(-n) = -n^2")
    print("  The 'symmetric part' s(n) = (a(n) - a(-n))/2 and")
    print("  'antisymmetric part' d(n) = (a(n) + a(-n))/2 = -n^2/2")
    print(f"{'='*70}")

    for n in range(0, 10):
        an = a_seq(n)
        amn = a_seq(-n)
        sym = (an - amn) // 2 if (an - amn) % 2 == 0 else (an - amn) / 2
        anti = (an + amn) // 2 if (an + amn) % 2 == 0 else (an + amn) / 2
        print(f"  n={n:2d}: a(n)={an:6d}, a(-n)={amn:6d}, "
              f"sym=(a-a')/2={sym}, anti=(a+a')/2={anti}")

    # The symmetric part: (a(n) - a(-n))/2
    # = (n + C(n,3) - (-n + C(-n,3))) / 2
    # = (2n + C(n,3) - C(-n,3)) / 2
    # C(-n,3) = (-n)(-n-1)(-n-2)/6 = -n(n+1)(n+2)/6
    # C(n,3) - C(-n,3) = n(n-1)(n-2)/6 + n(n+1)(n+2)/6 = n[(n-1)(n-2)+(n+1)(n+2)]/6
    # = n[2n^2+4]/6 = n(n^2+2)/3
    # sym = (2n + n(n^2+2)/3) / 2 = n(1 + (n^2+2)/6) = n(n^2+8)/6

    print(f"\n  Symmetric part: s(n) = n(n^2+8)/6")
    for n in range(0, 10):
        predicted = n * (n**2 + 8) // 6
        an = a_seq(n)
        amn = a_seq(-n)
        actual = (an - amn) // 2
        print(f"    n={n}: predicted={predicted}, actual={actual}, match={predicted==actual}")

    # The antisymmetric part: d(n) = (a(n)+a(-n))/2 = -n^2/2
    # This is ALWAYS -n^2/2 (proved in Part 3).

    print(f"\n  DECOMPOSITION: a(n) = s(n) + d(n)")
    print(f"  where s(n) = n(n^2+8)/6 (symmetric, s(n) = s(-n)... wait)")

    # Actually: s(n) = (a(n) - a(-n))/2, so s(-n) = (a(-n) - a(n))/2 = -s(n)
    # s is ANTISYMMETRIC! And d is symmetric (d(n) = d(-n) = -n^2/2).

    print(f"  CORRECTION: s(n) is ANTISYMMETRIC (s(-n) = -s(n))")
    print(f"  d(n) = -n^2/2 is SYMMETRIC (d(-n) = d(n))")
    print(f"")
    print(f"  a(n) = s(n) + d(n)")
    print(f"  a(-n) = -s(n) + d(n)")
    print(f"  Adding: a(n) + a(-n) = 2*d(n) = -n^2 ✓")
    print(f"  Subtracting: a(n) - a(-n) = 2*s(n) = n(n^2+8)/3")

    # ============================================================
    # INTERACTION 5: The Fibonacci-a(n) interaction at tournament sizes
    # ============================================================
    print(f"\n{'='*70}")
    print("INTERACTION 5: WHERE FIBONACCI MEETS a(n) = C(n,1)+C(n,3)")
    print(f"{'='*70}")

    # At what n does a(n) first exceed Fibonacci?
    for n in range(1, 20):
        an = a_seq(n)
        fn = fib(n)
        if an > fn and a_seq(n-1) <= fib(n-1):
            print(f"  a(n) first exceeds Fib at n={n}: a={an}, F={fn}")

    # At what n does a(n) first exceed C(n,2)?
    for n in range(1, 20):
        an = a_seq(n)
        cn2 = C(n, 2)
        if an > cn2 and a_seq(n-1) <= C(n-1, 2):
            print(f"  a(n) first exceeds C(n,2) at n={n}: a={an}, C(n,2)={cn2}")

    # At n=5: a(5) = 15, C(5,2) = 10, Fib(5) = 5
    # a crosses C(n,2) between n=4 and n=5? a(4)=8, C(4,2)=6. Already > !
    # a(3) = 4, C(3,2) = 3. Also >
    # a(2) = 2, C(2,2) = 1. Also >
    # a is ALWAYS > C(n,2) for n >= 2:
    # a(n) = n + C(n,3) vs C(n,2) = n(n-1)/2
    # a(n) - C(n,2) = n + C(n,3) - C(n,2) = n + C(n,3) - C(n,2)
    # = n + n(n-1)(n-2)/6 - n(n-1)/2 = n(1 + (n-1)(n-2)/6 - (n-1)/2)
    # = n(6 + (n-1)(n-2) - 3(n-1))/6 = n(6 + n^2 - 5n + 4)/6
    # = n(n^2 - 5n + 10)/6

    print(f"\n  a(n) - C(n,2) = n(n^2 - 5n + 10)/6")
    print(f"  Discriminant of n^2-5n+10: 25-40 = -15 < 0")
    print(f"  So n^2-5n+10 > 0 for ALL n! Hence a(n) > C(n,2) for all n >= 1.")

    # ============================================================
    # INTERACTION 6: The "layer" interpretation
    # ============================================================
    print(f"\n{'='*70}")
    print("INTERACTION 6: LAYERS OF PASCAL AND TOURNAMENT STRUCTURE")
    print(f"{'='*70}")

    print(f"""
  The sequence a(n) = C(n,1) + C(n,3) is the sum of Pascal layers 1 and 3.

  Each layer has a tournament interpretation:
    Layer 0: C(n,0) = 1 (the empty set — the "constant" term)
    Layer 1: C(n,1) = n (vertices — 0-dimensional structure)
    Layer 2: C(n,2) = n(n-1)/2 (edges/arcs — 1-dimensional structure)
    Layer 3: C(n,3) = n(n-1)(n-2)/6 (triangles — 2-dimensional structure)
    Layer 4: C(n,4) (tetrahedra — 3-dimensional)
    ...
    Layer k: C(n,k) (k-element subsets — (k-1)-dimensional)

  The OCF uses powers of 2 to weight cycle counts:
    H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
    = C(x,0) + 2*alpha_1 + 4*alpha_2 + ... where x counts at fugacity 2

  The connection:
    H counts weighted INDEPENDENT SETS (of the cycle graph)
    a(n) counts UNWEIGHTED subsets of sizes 1 and 3 (of the vertex set)
    Both are partial sums of a "layer decomposition"
    But the layers are different:
      H: layers are sizes of independent cycle sets (0,1,2,3,...)
      a(n): layers are sizes of vertex subsets (1,3)

  THE BRIDGE:
  At n=5: alpha_1 = (H-1)/2 ranges from 0 to 7.
  The max alpha_1 = 7 = C(5,1) + C(5,3)/... no, 7 ≠ 5+10.
  But a(5) = 15 = max_H(5) = the TOTAL structure count.
  And 15 = C(5,1) + C(5,3) = 5 + 10.
  The maximum Ham path count EQUALS the vertex + triangle count!
""")

    # ============================================================
    # INTERACTION 7: What is the sequence at n=6?
    # ============================================================
    print(f"\n{'='*70}")
    print("INTERACTION 7: THE SEQUENCE AT n=6 (THE CRITICAL TRANSITION)")
    print(f"{'='*70}")

    n = 6
    print(f"\n  n=6 values:")
    print(f"    a(6) = C(6,1) + C(6,3) = 6 + 20 = 26")
    print(f"    Tri(6) = 21 = FORBIDDEN H value")
    print(f"    Fib(6) = 8")
    print(f"    C(6,2) = 15 = max_H(5)")
    print(f"    max_H(6) = 45 = T(9)")
    print(f"    Jacobsthal(6) = 21 = FORBIDDEN!")
    print(f"")
    print(f"  At n=6: the 'forbidden' number 21 appears as BOTH")
    print(f"  Tri(6) AND Jacobsthal(6) AND C(7,2).")
    print(f"  While a(6) = 26 passes through, giving NO special alignment.")
    print(f"  But C(6,2) = 15 = max_H(5) = a(5) = Tri(5)!")
    print(f"  The 'echo' of n=5 appears in n=6's arc count.")

    # ============================================================
    # FINAL: The meta-pattern
    # ============================================================
    print(f"\n{'='*70}")
    print("THE META-PATTERN")
    print(f"{'='*70}")
    print(f"""
  The sequences align at SPECIFIC "resonance points":
    n=1: everything = 1 (trivial)
    n=3: a(3)=4=2^2, Fib(3)=2, Tri(3)=6, maxH(3)=3
    n=5: a(5)=Tri(5)=maxH(5)=15 (TRIPLE COINCIDENCE)
    n=7: a(7)=42=Catalan(5), Tri(7)=28=C(8,2), maxH(7)=189=3^3*7
    n=10: Fib(10)=Tri(10)=55, a(10)=130, C(10,2)=45=maxH(6)
    n=11: maxH(11)/Tri(10)=1729

  The RESONANCE POINTS are at n = 1, 5, 10 (and weakly 3, 7, 11).
  These are n ≡ 0, 1 (mod 5)!
  Or: the Fibonacci-indexed resonances.

  At each resonance, multiple sequences agree on a value that has
  TOURNAMENT SIGNIFICANCE (a max_H value, a forbidden value, or
  a structural quantity like arc count).

  The 2nd differences of a(n) are the natural numbers.
  The 2nd differences of Fibonacci are... Fibonacci itself!
  The 2nd differences of triangular are... 1, 1, 1, ... (constant)!

  HIERARCHY OF SIMPLICITY:
    Triangular: constant 2nd diff (simplest polynomial)
    a(n) = n+C(n,3): natural number 2nd diff (next simplest)
    Fibonacci: self-similar 2nd diff (simplest exponential)
    max_H: no simple diff pattern (most complex)
""")

    print(f"\n{'='*70}")
    print("DONE — DEEP INTERACTION ANALYSIS")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
