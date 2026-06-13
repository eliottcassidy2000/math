"""
mystery_sequence_deep.py -- kind-pasteur-2026-03-14-S90
LONG SESSION. Identify and deeply explore the sequence:
...-62,-40,-24,-13,-6,-2,0,1,2,4,8,15,26,42,64,93...

Also: continue the meta-alignment work from S89, connecting
ALL fundamental sequences at the same index rate.
"""

import sys, math
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding='utf-8')

def C(n, k):
    if k < 0 or k > n: return 0
    return math.comb(n, k)

def main():
    print("=" * 70)
    print("MYSTERY SEQUENCE + DEEP META-ALIGNMENT")
    print("kind-pasteur-2026-03-14-S90")
    print("=" * 70)

    # The given sequence (centered around 0)
    seq_neg = [-62, -40, -24, -13, -6, -2, 0]
    seq_pos = [0, 1, 2, 4, 8, 15, 26, 42, 64, 93]
    # Full sequence indexed from some negative number
    # Let's say a(0) = 0, a(1) = 1, a(2) = 2, a(3) = 4, a(4) = 8, a(5) = 15, ...
    # And a(-1) = -2, a(-2) = -6, a(-3) = -13, a(-4) = -24, a(-5) = -40, a(-6) = -62

    # ============================================================
    # PART 1: IDENTIFY THE SEQUENCE
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 1: IDENTIFY THE SEQUENCE")
    print(f"{'='*70}")

    seq = seq_pos  # Focus on non-negative part first
    print(f"\n  Positive part: {seq}")

    # First differences
    d1 = [seq[i+1] - seq[i] for i in range(len(seq)-1)]
    print(f"  1st differences: {d1}")

    # Second differences
    d2 = [d1[i+1] - d1[i] for i in range(len(d1)-1)]
    print(f"  2nd differences: {d2}")

    # Third differences
    d3 = [d2[i+1] - d2[i] for i in range(len(d2)-1)]
    print(f"  3rd differences: {d3}")

    # The first differences are: 1, 1, 2, 4, 7, 11, 16, 22, 29
    # These are 1, 1, 2, 4, 7, 11, 16, 22, 29
    # Second differences: 0, 1, 2, 3, 4, 5, 6, 7
    # THESE ARE JUST 0, 1, 2, 3, 4, 5, 6, 7 !!!

    print(f"\n  *** 2nd differences are 0, 1, 2, 3, 4, 5, 6, 7 ***")
    print(f"  *** = the natural numbers! ***")

    # So: d2(n) = n for n >= 0
    # d1(n) = sum_{k=0}^{n} k = n(n+1)/2 = triangular(n)
    # Wait: d1 = 1, 1, 2, 4, 7, 11, 16, 22, 29
    # Cumsum of d2 = [0, 0+1, 0+1+2, 0+1+2+3, ...] = [0, 1, 3, 6, 10, 15, 21, 28]
    # But d1 starts at 1, so d1(n) = 1 + T(n-1) where T is triangular?
    # d1(0) = 1, d1(1) = 1, d1(2) = 2, d1(3) = 4, d1(4) = 7
    # T(0)=0, T(1)=1, T(2)=3, T(3)=6, T(4)=10
    # d1(n) = 1 + T(n-1)? d1(0)=1+T(-1)=1+0=1 ✓, d1(1)=1+T(0)=1+0=1 ✓,
    # d1(2)=1+T(1)=1+1=2 ✓, d1(3)=1+T(2)=1+3=4 ✓, d1(4)=1+T(3)=1+6=7 ✓!

    print(f"\n  1st differences: d1(n) = 1 + T(n-1) = 1 + (n-1)n/2")
    print(f"  Verification:")
    for n in range(len(d1)):
        predicted = 1 + (n * (n-1)) // 2 if n >= 1 else 1
        if n == 0: predicted = 1
        else: predicted = 1 + n*(n-1)//2
        print(f"    n={n}: predicted = {predicted}, actual = {d1[n]}, match = {predicted == d1[n]}")

    # Hmm, not matching for n=0. Let me re-index.
    # d1 = [1, 1, 2, 4, 7, 11, 16, 22, 29]
    # d2 = [0, 1, 2, 3, 4, 5, 6, 7]
    # So d2(k) = k for k = 0, 1, 2, ...
    # d1(k) = d1(0) + sum_{j=0}^{k-1} d2(j) = 1 + sum_{j=0}^{k-1} j = 1 + k(k-1)/2
    print(f"\n  Corrected: d1(k) = 1 + k(k-1)/2")
    for k in range(len(d1)):
        predicted = 1 + k*(k-1)//2
        print(f"    k={k}: predicted = {predicted}, actual = {d1[k]}, match = {predicted == d1[k]}")

    # And the sequence itself: a(n) = a(0) + sum_{k=0}^{n-1} d1(k)
    # = 0 + sum_{k=0}^{n-1} (1 + k(k-1)/2)
    # = n + sum_{k=0}^{n-1} k(k-1)/2
    # = n + (1/2) * sum_{k=0}^{n-1} (k^2 - k)
    # = n + (1/2) * ((n-1)n(2n-1)/6 - (n-1)n/2)
    # = n + (n-1)n/2 * ((2n-1)/6 - 1/2)
    # = n + (n-1)n/2 * (2n-1-3)/6
    # = n + (n-1)n(2n-4)/12
    # = n + n(n-1)(n-2)/6
    # = n * (1 + (n-1)(n-2)/6)
    # = n * (6 + n^2 - 3n + 2) / 6
    # = n * (n^2 - 3n + 8) / 6

    # Wait, let me just compute: a(n) = sum_{k=0}^{n-1} (1 + k(k-1)/2)
    # = n + (1/2) sum_{k=2}^{n-1} k(k-1)
    # sum_{k=2}^{n-1} k(k-1) = sum_{k=0}^{n-1} k(k-1) = 2*C(n,3)... hmm
    # Actually sum_{k=0}^{m} k(k-1) = 2*C(m+1, 3) by hockey stick

    # sum_{k=0}^{n-1} k(k-1)/2 = C(n, 3)

    print(f"\n  a(n) = n + C(n, 3)")
    print(f"  Verification:")
    for n in range(len(seq)):
        predicted = n + C(n, 3)
        print(f"    n={n}: predicted = {predicted}, actual = {seq[n]}, match = {predicted == seq[n]}")

    # CHECK: a(n) = n + C(n,3) = n + n(n-1)(n-2)/6
    # = n(1 + (n-1)(n-2)/6) = n(6 + n^2 - 3n + 2)/6 = n(n^2 - 3n + 8)/6

    # Verify with more values:
    print(f"\n  Extended: a(n) = n + C(n,3) for n = -6..15:")
    for n in range(-6, 16):
        # C(n, 3) for negative n: C(n,3) = n(n-1)(n-2)/6
        cn3 = n * (n-1) * (n-2) // 6
        val = n + cn3
        print(f"    a({n:3d}) = {n:4d} + {cn3:6d} = {val:6d}")

    # Check against the given negative values:
    # a(-1) should be -2: -1 + C(-1,3) = -1 + (-1)(-2)(-3)/6 = -1 + (-6)/6 = -1 + (-1) = -2 ✓
    # a(-2) = -2 + C(-2,3) = -2 + (-2)(-3)(-4)/6 = -2 + (-24)/6 = -2 + (-4) = -6 ✓
    # a(-3) = -3 + C(-3,3) = -3 + (-3)(-4)(-5)/6 = -3 + (-60)/6 = -3 + (-10) = -13 ✓
    print(f"\n  Matches with given negative values:")
    given_neg = {-1: -2, -2: -6, -3: -13, -4: -24, -5: -40, -6: -62}
    for n, expected in sorted(given_neg.items()):
        cn3 = n * (n-1) * (n-2) // 6
        val = n + cn3
        print(f"    a({n}) = {val}, expected = {expected}, match = {val == expected}")

    print(f"\n  *** THE SEQUENCE IS a(n) = n + C(n, 3) ***")
    print(f"  *** = n + n(n-1)(n-2)/6 ***")
    print(f"  *** = n(n^2 - 3n + 8)/6 ***")

    # ============================================================
    # PART 2: WHAT IS n + C(n,3) ?
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 2: WHAT IS n + C(n,3)?")
    print(f"{'='*70}")

    # n + C(n,3) = C(n,1) + C(n,3)
    # This is the sum of ODD binomial coefficients up to some point!
    # Actually: C(n,1) + C(n,3) is a partial sum of the odd-indexed
    # binomial coefficients.

    # The FULL sum of odd-indexed: C(n,1)+C(n,3)+C(n,5)+... = 2^{n-1}
    # So n + C(n,3) = C(n,1) + C(n,3) = partial sum.

    print(f"  n + C(n,3) = C(n,1) + C(n,3)")
    print(f"  This is the sum of the FIRST TWO odd binomial coefficients!")
    print(f"")
    print(f"  Full odd sum: C(n,1) + C(n,3) + C(n,5) + ... = 2^(n-1)")
    print(f"  Our sequence = the PARTIAL sum C(n,1) + C(n,3)")
    print(f"")
    print(f"  Compare with:")
    print(f"  C(n,0) + C(n,2) + C(n,4) + ... = 2^(n-1) (even sum)")
    print(f"")
    print(f"  So a(n)/2^(n-1) -> 0 as n -> infinity (only 2 of ~n/2 terms)")

    # But there's a BEAUTIFUL interpretation:
    # In tournaments: 2^{n-1} = 2^{C(n,2)}/... no.
    # Actually: the OCF formula H = 1 + 2*alpha_1 + 4*alpha_2 + ...
    # = C(0,0)*2^0 + C(alpha_1,1)*2^1 + ... if we think of it
    # as a "binomial expansion in 2".

    # The sequence C(n,1) + C(n,3) = n + n(n-1)(n-2)/6
    # appears in combinatorics as the number of ways to choose
    # 1 or 3 items from n.

    # This connects to TOURNAMENT THEORY:
    # alpha_1 = #{odd cycles of size 1 or 3}... no, alpha_1 = total independent cycles.
    # But the OCF is H = sum 2^k * alpha_k.
    # If we set alpha_1 = n, alpha_3 = C(n,3), then
    # H = 1 + 2n + 4*0 + 8*C(n,3) = 1 + 2n + 8*C(n,3)... no, that's different.

    # ============================================================
    # PART 3: THE SEQUENCE AND STRANDS
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 3: STRANDS OF a(n) = C(n,1) + C(n,3)")
    print(f"{'='*70}")

    # Splitting into even/odd and mod-3 strands:
    print(f"\n  Even-index strand: a(0), a(2), a(4), a(6), ...")
    for n in range(0, 16, 2):
        val = n + C(n, 3)
        print(f"    a({n:2d}) = {val:6d}")

    print(f"\n  Odd-index strand: a(1), a(3), a(5), a(7), ...")
    for n in range(1, 16, 2):
        val = n + C(n, 3)
        print(f"    a({n:2d}) = {val:6d}")

    # The even/odd symmetry:
    # a(-n) = -n + C(-n, 3) = -n + (-n)(-n-1)(-n-2)/6
    # = -n - n(n+1)(n+2)/6
    # a(n) = n + n(n-1)(n-2)/6
    # a(n) + a(-n) = n(n-1)(n-2)/6 - n(n+1)(n+2)/6
    # = n/6 * [(n-1)(n-2) - (n+1)(n+2)]
    # = n/6 * [n^2 - 3n + 2 - n^2 - 3n - 2]
    # = n/6 * (-6n) = -n^2

    print(f"\n  SYMMETRY: a(n) + a(-n) = -n^2")
    print(f"  Verification:")
    for n in range(1, 7):
        an = n + C(n, 3)
        amn = -n + (-n)*(-n-1)*(-n-2)//6
        print(f"    a({n}) + a({-n}) = {an} + {amn} = {an + amn}, -n^2 = {-n**2}, "
              f"match = {an + amn == -n**2}")

    print(f"\n  So: a(n) = (-n^2 - a(-n))/1... or better:")
    print(f"  a(n) - a(-n) = 2n + n(n-1)(n-2)/6 + n(n+1)(n+2)/6")
    print(f"  = 2n + n[(n-1)(n-2) + (n+1)(n+2)]/6")
    print(f"  = 2n + n[2n^2 + 4]/6")
    print(f"  = 2n + n(n^2+2)/3")

    # ============================================================
    # PART 4: CONNECT TO TOURNAMENT THEORY
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 4: TOURNAMENT CONNECTIONS OF C(n,1) + C(n,3)")
    print(f"{'='*70}")

    # The number C(n,1) + C(n,3) at specific tournament values:
    for n in range(1, 12):
        val = n + C(n, 3)
        m = C(n, 2)  # arcs
        ratio_str = f"{val/m:.4f}" if m > 0 else "inf"
        print(f"  n={n:2d}: a(n)={val:6d}, C(n,2)={m:3d}, ratio={ratio_str}")

    # At n=3: a(3) = 3 + 1 = 4. And C(3,2) = 3.
    # At n=4: a(4) = 4 + 4 = 8 = 2^3 = number of tournaments on 3 vertices!
    # At n=5: a(5) = 5 + 10 = 15 = max_H(5)!!!
    # At n=6: a(6) = 6 + 20 = 26.
    # At n=7: a(7) = 7 + 35 = 42 = Catalan(5) = C(10,5)/6

    print(f"\n  *** a(5) = 15 = max_H(5) !!! ***")
    print(f"  *** a(7) = 42 = Catalan(5) ***")
    print(f"  *** a(4) = 8 = 2^3 = total tournaments on n=3 ***")
    print(f"  *** a(3) = 4 = 2^C(2,1) = total tournaments on n=2 ***")

    # The sequence C(n,1) + C(n,3) encodes "1st + 3rd layer of Pascal."
    # In the tournament context: n arcs (1st order) plus C(n,3) triangles (3rd order).
    # The number of ARC PLUS TRIANGLE structures.

    # ============================================================
    # PART 5: THE "PASCAL LAYER" SEQUENCES
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 5: PASCAL LAYER SEQUENCES")
    print("  Layer k: L_k(n) = C(n, 2k+1) (odd binomial coefficients)")
    print("  Our sequence = L_0 + L_1 = C(n,1) + C(n,3)")
    print(f"{'='*70}")

    for k in range(5):
        layer = [C(n, 2*k+1) for n in range(12)]
        print(f"  Layer {k} = C(n, {2*k+1}): {layer}")

    # Sum of layers 0..K:
    for K in range(4):
        partial = [sum(C(n, 2*k+1) for k in range(K+1)) for n in range(12)]
        print(f"  Sum layers 0..{K}: {partial}")

    # Full odd sum = 2^{n-1}:
    full_odd = [2**(n-1) if n > 0 else 0 for n in range(12)]
    print(f"  Full odd sum (2^(n-1)): {full_odd}")

    # The "ODD INDEPENDENCE NUMBER" interpretation:
    # C(n, 2k+1) counts ways to choose an odd number of items from n.
    # In tournament theory: choosing odd subsets of vertices...
    # Each odd cycle has an odd number of vertices.
    # C(n, 3) = number of 3-element subsets = potential 3-cycles.
    # C(n, 5) = number of 5-element subsets = potential 5-cycles.
    # The sequence sum C(n, 2k+1) for k=0,...,K is a "truncated odd partition function."

    print(f"\n  TOURNAMENT INTERPRETATION:")
    print(f"  C(n, 1) = n = number of vertices (trivial 1-cycles)")
    print(f"  C(n, 3) = number of potential 3-cycles (vertex triples)")
    print(f"  C(n, 5) = number of potential 5-cycles (vertex quintuples)")
    print(f"  C(n, 1) + C(n, 3) = vertices + potential triangles")
    print(f"  This counts the 'first two levels of cycle complexity'")

    # ============================================================
    # PART 6: RATE COMPARISON WITH FIBONACCI AND TRIANGULAR
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 6: GROWTH RATE OF a(n) = n + C(n,3)")
    print("  a(n) ~ n^3/6 (cubic polynomial)")
    print("  Compare: Tri(n) ~ n^2/2 (quadratic)")
    print("  Compare: Fib(n) ~ phi^n (exponential)")
    print(f"{'='*70}")

    def fib(n):
        if n <= 0: return 0
        a, b = 0, 1
        for _ in range(n): a, b = b, a + b
        return a

    def tri(n): return n*(n+1)//2

    print(f"\n  {'n':>3} {'a(n)':>8} {'Tri':>8} {'Fib':>8} {'a/Tri':>8} {'a/Fib':>8} {'Tri/Fib':>8}")
    for n in range(2, 16):
        an = n + C(n, 3)
        tn = tri(n)
        fn = fib(n)
        print(f"  {n:3d} {an:8d} {tn:8d} {fn:8d} "
              f"{an/tn if tn else 0:8.4f} "
              f"{an/fn if fn else 0:8.4f} "
              f"{tn/fn if fn else 0:8.4f}")

    # a(n)/Tri(n) -> n/3 (grows linearly)
    # a(n)/Fib(n) -> infinity (polynomial vs exponential)
    # Tri(n)/Fib(n) -> 0 (polynomial vs exponential, but slowly)

    print(f"\n  Key crossings:")
    print(f"    a(n) = Tri(n) at n where n+C(n,3) = n(n+1)/2")
    print(f"    C(n,3) = n(n+1)/2 - n = n(n-1)/2 = C(n,2)")
    print(f"    C(n,3) = C(n,2) when n(n-1)(n-2)/6 = n(n-1)/2")
    print(f"    (n-2)/3 = 1, so n = 5")
    print(f"    Verify: a(5) = 15, Tri(5) = 15 ✓ !!!")

    print(f"\n  *** a(n) = Tri(n) at n=5, and a(5) = 15 = max_H(5) ***")
    print(f"  *** Three sequences meet at the tournament maximizer! ***")

    # Where does a(n) cross Fibonacci?
    for n in range(2, 20):
        an = n + C(n, 3)
        fn = fib(n)
        if an <= fn and n + C(n+1, 3) > fib(n+1):
            print(f"\n  a(n) crosses Fib between n={n} and n={n+1}")
            print(f"    a({n})={an}, F({n})={fn}")
            print(f"    a({n+1})={n+1+C(n+1,3)}, F({n+1})={fib(n+1)}")

    # ============================================================
    # PART 7: THE GRAND MULTI-SEQUENCE ALIGNMENT
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 7: GRAND ALIGNMENT — ALL SEQUENCES AT KEY POINTS")
    print(f"{'='*70}")

    maxH = [0, 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095]

    for n in [3, 5, 6, 7, 10, 11]:
        an = n + C(n, 3)
        tn = tri(n)
        fn = fib(n)
        cn = C(2*n, n)
        mh = maxH[n] if n < len(maxH) else '?'
        arcs = C(n, 2)

        print(f"\n  n = {n}:")
        print(f"    a(n) = C(n,1)+C(n,3) = {an}")
        print(f"    Tri(n) = {tn}")
        print(f"    Fib(n) = {fn}")
        print(f"    C(2n,n) = {cn}")
        print(f"    max_H(n) = {mh}")
        print(f"    C(n,2) = {arcs}")

        # Check equalities
        checks = []
        if an == tn: checks.append("a(n) = Tri(n)")
        if an == fn: checks.append("a(n) = Fib(n)")
        if an == mh: checks.append("a(n) = max_H(n)")
        if an == arcs: checks.append("a(n) = C(n,2)")
        if tn == fn: checks.append("Tri(n) = Fib(n)")
        if tn == mh: checks.append("Tri(n) = max_H(n)")
        if fn == arcs: checks.append("Fib(n) = C(n,2)")
        if isinstance(mh, int) and mh == tn: checks.append("max_H = Tri(n)")
        if isinstance(mh, int) and mh == an: checks.append("max_H = a(n)")

        if checks:
            for c in checks:
                print(f"    *** {c} ***")

    # ============================================================
    # PART 8: THE n=5 TRIPLE COINCIDENCE
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 8: THE n=5 TRIPLE COINCIDENCE")
    print("  At n=5: a(5) = Tri(5) = max_H(5) = 15")
    print("  THREE different sequences meet at the tournament maximizer!")
    print(f"{'='*70}")

    print(f"""
  a(5) = 5 + C(5,3) = 5 + 10 = 15
  Tri(5) = 5*6/2 = 15
  max_H(5) = 15

  WHY?
  a(n) = Tri(n) iff C(n,3) = C(n,2) iff n = 5 (unique solution for n >= 3).
  max_H(5) = 15 = n!/2^(n-2) = 120/8 = 15 (regular tournament formula).

  So n=5 is the UNIQUE point where:
  - The "odd Pascal partial sum" C(n,1)+C(n,3) equals the triangular number
  - Both equal the max Hamiltonian path count
  - Both equal 15 = the 5th triangular number = C(6,2)

  This makes n=5 a UNIVERSAL FIXED POINT of the sequence alignment.
  It's also the LAST n where:
  - Omega(T) is always claw-free
  - The H-landscape is unimodal
  - alpha_2 = 0 (no disjoint cycle pairs)
  - The 3-strand Pascal aligns with Jacobsthal for max_H
""")

    # ============================================================
    # PART 9: WHAT IS C(n,1) + C(n,3) IN THE INDEPENDENCE WORLD?
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 9: C(n,1)+C(n,3) AS PARTIAL INDEPENDENCE POLYNOMIAL")
    print(f"{'='*70}")

    # I(K_n, x) = (1+x)^n (independence polynomial of the EMPTY graph on n vertices)
    # The odd part: sum_{k odd} C(n,k) x^k = ((1+x)^n - (1-x)^n) / 2
    # At x=1: 2^{n-1}

    # C(n,1) + C(n,3) is I_odd(K_n, x) truncated to degree 3, at x=1.
    # = n + n(n-1)(n-2)/6

    # In the TOURNAMENT framework:
    # I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2 + ...
    # The odd-indexed terms: alpha_1*x + alpha_3*x^3 + ...
    # At x=1: alpha_1 + alpha_3 + ... = total odd independent sets of odd size

    # But C(n,1) + C(n,3) is from the EMPTY graph (= K_n complement = no edges).
    # For the empty graph: every subset is independent.
    # C(n,1) + C(n,3) = #singletons + #triples from n elements.

    # DEEP CONNECTION:
    # For a tournament with alpha_1 = n (impossible for small n, but conceptually):
    # H = 1 + 2n. And alpha_1 = n means n independent odd cycles.
    # If additionally alpha_3 = C(n,3): all triples of these cycles are independent.
    # Then I(Omega, x) = (1+x)^n (the empty graph!).
    # H = I(Omega, 2) = 3^n.

    print(f"  If Omega = empty graph on n vertices: I(Omega, 2) = 3^n")
    print(f"  3^n sequence: {[3**n for n in range(10)]}")
    print(f"  max_H values: {maxH[1:]}")
    print(f"")
    print(f"  max_H(n) vs 3^{{alpha_1_max}}:")
    # At n=5: alpha_1 max = 7, 3^7 = 2187 >> 15 = max_H
    # The conflict graph is FAR from empty!
    # C(n,1)+C(n,3) is the "what if Omega were empty up to triples" counting.

    # ============================================================
    # PART 10: GOING FURTHER — a(n) as a "complexity measure"
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 10: a(n) = n + C(n,3) AS TOURNAMENT COMPLEXITY")
    print(f"{'='*70}")

    print(f"""
  INTERPRETATION:
  a(n) = n + C(n,3) counts the total number of "simple structures"
  in a tournament on n vertices:
    - n vertices (0-dimensional structure)
    - C(n,3) potential 3-cycles (1-dimensional structure)

  This is the EULER CHARACTERISTIC of the "potential cycle complex":
    chi = V - E + F = n - C(n,2) + C(n,3)
    = n - n(n-1)/2 + n(n-1)(n-2)/6

  Let me compute this:
""")

    for n in range(1, 12):
        V = n
        E = C(n, 2)
        F = C(n, 3)
        chi = V - E + F
        an = n + C(n, 3)
        print(f"  n={n:2d}: V={V:3d}, E={E:3d}, F={F:4d}, "
              f"chi=V-E+F={chi:5d}, a(n)=V+F={an:5d}")

    # chi = n - C(n,2) + C(n,3) is actually C(n-1, 0) - C(n-1,1) + C(n-1,2) = ...
    # By the alternating sum: sum (-1)^k C(n,k) = 0 for n >= 1.
    # So chi = C(n,0) - C(n,1) + C(n,2) - C(n,3) + ... = 0 (for n >= 1)
    # But here we truncate: V-E+F = C(n,1)-C(n,2)+C(n,3)
    # = n - n(n-1)/2 + n(n-1)(n-2)/6 = n[1 - (n-1)/2 + (n-1)(n-2)/6]
    # = n[(6 - 3(n-1) + (n-1)(n-2))/6]
    # = n[(6 - 3n + 3 + n^2 - 3n + 2)/6]
    # = n[(n^2 - 6n + 11)/6]

    print(f"\n  chi = n(n^2 - 6n + 11)/6")
    print(f"  a(n) = n(n^2 - 3n + 8)/6")
    print(f"  Difference: a(n) - chi = n(3n - 3)/6 = n(n-1)/2 = C(n,2)")
    print(f"  So: a(n) = chi + C(n,2)")
    print(f"  Or: a(n) = (V - E + F) + E = V + F = n + C(n,3)")

    print(f"\n  This is just the UNSIGNED version:")
    print(f"  a(n) = #{'{'}vertices{'}'} + #{'{'}triangles{'}'}")
    print(f"  chi(n) = #{'{'}vertices{'}'} - #{'{'}edges{'}'} + #{'{'}triangles{'}'}")

    print(f"\n{'='*70}")
    print("DONE — DEEP MULTI-SEQUENCE EXPLORATION")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
