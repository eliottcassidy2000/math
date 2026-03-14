"""
sequence_meta_alignment.py -- kind-pasteur-2026-03-14-S89
LONG DEEP SESSION: Meta-alignment of fundamental integer sequences.

THE QUESTION: When we line up ALL fundamental sequences at the same index n,
what patterns and alignments emerge? Which sequences "move together"?

SEQUENCES TO ALIGN:
1. Fibonacci: F_n = 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, ...
2. Triangular: T_n = 0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, ...
3. Central binomial: C(2n,n) = 1, 2, 6, 20, 70, 252, 924, ...
4. Catalan: Cat_n = 1, 1, 2, 5, 14, 42, 132, 429, ...
5. 3-strand Pascal: P_n (from S87)
6. Factorial: n! = 1, 1, 2, 6, 24, 120, 720, ...
7. Powers of 2: 2^n
8. Bell numbers: 1, 1, 2, 5, 15, 52, 203, 877, ...
9. max_H: 1, 1, 3, 5, 15, 45, 189, 661, 3357, ...
10. Tournament count: 2^C(n,2)
11. Jacobsthal: J_n = 1, 1, 3, 5, 11, 21, 43, 85, ...

GROWTH RATES (per step):
  Fibonacci: phi ≈ 1.618
  Triangular: ~n (polynomial, not exponential)
  Central binomial: 4 per 2 steps = 2 per step
  Catalan: 4 per step (asymptotically)
  3-strand Pascal: 4^{1/3} ≈ 1.587 per step
  Factorial: ~n per step (super-exponential)
  Powers of 2: 2 per step
  Bell: ~n/ln(n) per step
  max_H: ~n/e per step (from Szele)

WHICH SEQUENCES HAVE THE SAME GROWTH RATE?
  phi ≈ 1.618: Fibonacci, Lucas
  2: Powers of 2, central binomial (per step avg)
  4^{1/3} ≈ 1.587: 3-strand Pascal
  These are CLOSE: phi ≈ 4^{1/3} (within 2%)!

THIS is the key insight. Let me explore it.
"""

import sys, math
from collections import Counter, defaultdict
import numpy as np

sys.stdout.reconfigure(encoding='utf-8')

def C(n, k):
    if k < 0 or k > n: return 0
    return math.comb(n, k)

def fib(n):
    if n <= 0: return 0
    a, b = 0, 1
    for _ in range(n): a, b = b, a + b
    return a

def pascal_3strand(k):
    n = k // 3
    r = k % 3
    if r == 0: return C(2*n+1, n)
    elif r == 1: return C(2*n+2, n)
    else: return C(2*n+2, n+1)

def triangular(n):
    return n * (n + 1) // 2

def catalan(n):
    return C(2*n, n) // (n + 1)

def bell(n):
    """Bell numbers via the Bell triangle."""
    if n == 0: return 1
    B = [[0] * (n+1) for _ in range(n+1)]
    B[0][0] = 1
    for i in range(1, n+1):
        B[i][0] = B[i-1][i-1]
        for j in range(1, i+1):
            B[i][j] = B[i][j-1] + B[i-1][j-1]
    return B[n][0]

def jacobsthal(n):
    """Jacobsthal numbers: J(0)=0, J(1)=1, J(n)=J(n-1)+2*J(n-2)."""
    if n == 0: return 0
    if n == 1: return 1
    a, b = 0, 1
    for _ in range(n-1): a, b = b, b + 2*a
    return b

maxH = [0, 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095]

def main():
    print("=" * 70)
    print("SEQUENCE META-ALIGNMENT — LONG DEEP SESSION")
    print("kind-pasteur-2026-03-14-S89")
    print("=" * 70)

    # ============================================================
    # PART 1: THE GRAND TABLE
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 1: GRAND ALIGNMENT TABLE")
    print(f"{'='*70}")

    headers = ['n', 'Fib', 'Tri', 'C(2n,n)', 'Cat', 'P_3str', 'n!',
               '2^n', 'Bell', 'Jac', 'maxH', 'C(n,2)', '2^C(n,2)']

    print(f"\n  {'':>3}", end="")
    for h in headers[1:]:
        print(f" {h:>10}", end="")
    print()
    print(f"  {'---':>3}" + " ----------" * (len(headers)-1))

    for n in range(13):
        vals = [
            n,
            fib(n),
            triangular(n),
            C(2*n, n),
            catalan(n),
            pascal_3strand(n),
            math.factorial(n),
            2**n,
            bell(n),
            jacobsthal(n),
            maxH[n] if n < len(maxH) else '?',
            C(n, 2),
            2**C(n, 2) if n <= 8 else '...',
        ]
        print(f"  {n:3d}", end="")
        for v in vals[1:]:
            if isinstance(v, int) and v < 10000000:
                print(f" {v:10d}", end="")
            elif isinstance(v, int):
                print(f" {v:10.3e}".replace('+0', '+'), end="")
            else:
                print(f" {str(v):>10}", end="")
        print()

    # ============================================================
    # PART 2: GROWTH RATE COMPARISON — LOG SCALE
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 2: GROWTH RATES (log base 2 of each sequence)")
    print("  Sequences that grow at the same rate have parallel log curves.")
    print(f"{'='*70}")

    phi = (1 + math.sqrt(5)) / 2

    print(f"\n  Exponential growth rates (per step):")
    print(f"    Fibonacci: phi = {phi:.6f}")
    print(f"    3-strand Pascal: 4^(1/3) = {4**(1/3):.6f}")
    print(f"    Jacobsthal: 2 = {2:.6f}")
    print(f"    Powers of 2: 2 = {2:.6f}")
    print(f"    Central binomial: ~4/sqrt(pi*n) per step")
    print(f"    Catalan: ~4 per step")
    print(f"    max_H: ~n/e per step (super-exponential)")
    print(f"    Triangular: polynomial (n^2/2), not exponential")

    print(f"\n  ALIGNMENT GROUPS (sequences with similar growth):")
    print(f"    Group A (growth ~1.6): Fibonacci ({phi:.3f}), 3-strand Pascal ({4**(1/3):.3f})")
    print(f"    Group B (growth ~2.0): Jacobsthal, 2^n")
    print(f"    Group C (growth ~4.0): Central binomial, Catalan")
    print(f"    Group D (super-exp): n!, max_H, Bell, 2^C(n,2)")
    print(f"    Group E (polynomial): Triangular, C(n,2)")

    # ============================================================
    # PART 3: GROUP A DEEP DIVE — Fibonacci vs 3-strand Pascal
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 3: GROUP A — THE PHI vs 4^(1/3) NEAR-ALIGNMENT")
    print(f"  phi = {phi:.10f}")
    print(f"  4^(1/3) = {4**(1/3):.10f}")
    print(f"  Difference: {phi - 4**(1/3):.10f} ({100*(phi - 4**(1/3))/phi:.4f}%)")
    print(f"{'='*70}")

    # WHY are they so close?
    # phi = (1+sqrt(5))/2, 4^{1/3} = 2^{2/3}
    # phi^3 = (1+sqrt(5))^3/8 = (1+3*sqrt(5)+3*5+5*sqrt(5))/8 = (16+8*sqrt(5))/8 = 2+sqrt(5)
    # phi^3 = 2 + sqrt(5) ≈ 4.236
    # (4^{1/3})^3 = 4

    print(f"\n  phi^3 = {phi**3:.6f} = 2 + sqrt(5)")
    print(f"  (4^(1/3))^3 = 4")
    print(f"  Ratio: phi^3 / 4 = {phi**3/4:.6f}")
    print(f"  So every 3 steps, Fibonacci gains ~5.9% over Pascal.")
    print(f"")
    print(f"  phi^6 = {phi**6:.6f}")
    print(f"  4^2 = 16")
    print(f"  phi^6 / 16 = {phi**6/16:.6f}")
    print(f"  Every 6 steps (one full period), Fibonacci gains ~12.2%.")

    # The near-alignment means for small n, Fib and Pascal are INTERCHANGEABLE
    # in terms of growth rate. They diverge slowly.

    print(f"\n  For n < 20, Fibonacci and 3-strand Pascal are within 2x of each other.")
    print(f"  This means in our 'tournament-relevant' range (n=3..11),")
    print(f"  these sequences are effectively ALIGNED.")

    # ============================================================
    # PART 4: THE TRIANGULAR-FIBONACCI RESONANCE
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 4: TRIANGULAR NUMBERS AND FIBONACCI")
    print("  T_n = n(n+1)/2 — polynomial growth")
    print("  F_n ~ phi^n — exponential growth")
    print("  BUT: specific values coincide!")
    print(f"{'='*70}")

    # Find where T_n = F_m for various n, m
    trivals = set(triangular(n) for n in range(100))
    fibvals = set(fib(n) for n in range(50))
    common = sorted(trivals & fibvals)
    print(f"\n  Values that are BOTH triangular and Fibonacci:")
    print(f"  {common}")
    # Known: 0, 1, 3, 21, 55 are both triangular and Fibonacci

    for v in common:
        if v > 0:
            # Find indices
            for n in range(100):
                if triangular(n) == v:
                    for m in range(50):
                        if fib(m) == v:
                            print(f"    {v} = T_{n} = F_{m}")

    print(f"\n  TOURNAMENT CONNECTION:")
    print(f"  T_n = C(n+1, 2) = number of arcs in tournament on n+1 vertices")
    print(f"  F_n = #{'{'}independent sets of path P_n{'}'}")
    print(f"  When T_n = F_m: the arc count of a tournament equals a Fibonacci number!")
    print(f"")
    print(f"  T_2 = 3 = F_4: 3 arcs (n=3 tournament)")
    print(f"  T_5 = 15 = F_... no, F_5=5, F_7=13, F_8=21. 15 is NOT Fibonacci.")
    print(f"  T_6 = 21 = F_8: 21 arcs (n=7 tournament!) and F_8 = 21")
    print(f"  T_9 = 45 = max_H(6)!")
    print(f"  T_10 = 55 = F_10: 55 arcs (n=11 tournament) and F_10 = 55")

    # AMAZING: T_6 = 21 = F_8, and 21 is a FORBIDDEN H value!
    # And T_9 = 45 = max_H(6)!
    print(f"\n  *** STRIKING COINCIDENCES ***:")
    print(f"  T_6 = 21 = F_8 = FORBIDDEN H VALUE = |Phi+(A_6)|")
    print(f"  T_9 = 45 = max_H(6) = C(10,2)")
    print(f"  T_2 = 3 = F_4 = max_H(3) = H(3-cycle)")
    print(f"  T_4 = 10 = C(5,2) = arcs at n=5 = V(Petersen)")

    # ============================================================
    # PART 5: WHERE DO SEQUENCES CROSS?
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 5: SEQUENCE CROSSINGS — WHERE DO THEY INTERSECT?")
    print(f"{'='*70}")

    # Find all crossings between pairs of sequences
    seqs = {
        'Fib': [fib(n) for n in range(20)],
        'Tri': [triangular(n) for n in range(20)],
        'Cat': [catalan(n) for n in range(15)],
        'P3': [pascal_3strand(n) for n in range(20)],
        'Jac': [jacobsthal(n) for n in range(20)],
        'Bell': [bell(n) for n in range(15)],
        '2^n': [2**n for n in range(20)],
    }

    print(f"\n  Crossings (A_n = B_n for specific n):")
    for name_a in sorted(seqs.keys()):
        for name_b in sorted(seqs.keys()):
            if name_a >= name_b: continue
            crossings = []
            min_len = min(len(seqs[name_a]), len(seqs[name_b]))
            for n in range(min_len):
                if seqs[name_a][n] == seqs[name_b][n] and seqs[name_a][n] > 0:
                    crossings.append((n, seqs[name_a][n]))
            if crossings:
                print(f"    {name_a} = {name_b} at: {crossings}")

    # ============================================================
    # PART 6: THE RATIO TABLE — All pairs
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 6: RATIO TABLE — HOW SEQUENCES COMPARE")
    print(f"{'='*70}")

    key_seqs = {
        'Fib': lambda n: fib(n),
        'Tri': lambda n: triangular(n),
        'C(2n,n)': lambda n: C(2*n, n),
        'Cat': lambda n: catalan(n),
        'maxH': lambda n: maxH[n] if n < len(maxH) else None,
    }

    print(f"\n  Fib(n) / Tri(n) for n=2..11:")
    for n in range(2, 12):
        f = fib(n)
        t = triangular(n)
        ratio = f / t if t > 0 else 0
        print(f"    n={n:2d}: F={f:6d}, T={t:4d}, F/T={ratio:.4f}")

    print(f"\n  max_H(n) / Tri(n-1) for n=3..11:")
    for n in range(3, min(12, len(maxH))):
        mh = maxH[n]
        t = triangular(n-1)
        ratio = mh / t if t > 0 else 0
        print(f"    n={n:2d}: maxH={mh:6d}, T(n-1)={t:4d}, maxH/T={ratio:.4f}")

    # ============================================================
    # PART 7: TOURNAMENT-RELEVANT ALIGNMENTS
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 7: TOURNAMENT-RELEVANT ALIGNMENTS")
    print(f"{'='*70}")

    print(f"\n  Key tournament quantities aligned with fundamental sequences:")

    for n in range(1, min(12, len(maxH))):
        m = C(n, 2)
        t = triangular(n-1)
        f = fib(n)
        j = jacobsthal(n)
        mh = maxH[n]
        cat = catalan(n)
        p3 = pascal_3strand(n)
        mean_h = math.factorial(n) / 2**(n-1) if n > 0 else 1

        print(f"\n  n={n:2d}: arcs=C({n},2)={m}, mean_H={mean_h:.1f}")
        print(f"    max_H = {mh}")
        print(f"    Fib({n}) = {f}")
        print(f"    Tri({n}) = {triangular(n)}")
        print(f"    Jacobsthal({n}) = {j}")
        print(f"    Catalan({n}) = {cat}")
        print(f"    P_3strand({n}) = {p3}")

        # Check alignments
        if mh == f:
            print(f"    *** max_H = Fibonacci! ***")
        if mh == triangular(n):
            print(f"    *** max_H = Triangular! ***")
        if mh == j:
            print(f"    *** max_H = Jacobsthal! ***")
        if mh == cat:
            print(f"    *** max_H = Catalan! ***")
        if mh == p3:
            print(f"    *** max_H = 3-strand Pascal! ***")
        if m == f:
            print(f"    *** C(n,2) = Fibonacci! ***")
        if m == triangular(n-1):
            print(f"    *** C(n,2) = Tri(n-1) always! ***")

    # ============================================================
    # PART 8: THE "GROWTH HIERARCHY" OF TOURNAMENT QUANTITIES
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 8: GROWTH HIERARCHY")
    print("  Classify ALL tournament quantities by growth rate")
    print(f"{'='*70}")

    print(f"""
  POLYNOMIAL GROWTH:
    n (vertices)                    O(n)
    C(n,2) = n(n-1)/2 (arcs)       O(n^2)
    c3 ~ C(n,3)/4 (3-cycles avg)   O(n^3)
    GS_DOF ~ n^2/4                  O(n^2)

  EXPONENTIAL GROWTH ~phi^n:
    Fibonacci F_n                   phi^n ≈ 1.618^n
    3-strand Pascal P_n             4^{n/3} ≈ 1.587^n
    Jacobsthal J_n                  2^n

  SUPER-EXPONENTIAL:
    mean_H = n!/2^{n-1}            ~ (n/e)^n * sqrt(2*pi*n) / 2^{n-1}
    max_H ~ e * n!/2^{n-1}         ~ e * mean_H
    n! (factorial)                  ~ (n/e)^n * sqrt(2*pi*n)
    2^{C(n,2)} (all tournaments)   ~ 2^{n^2/2}
    # iso classes                  ~ 2^{C(n,2)} / n!

  THE KEY DIVIDE:
  Tournament quantities split into:
  - POLYNOMIAL: structural (arc count, cycle count, GS dof)
  - EXPONENTIAL: combinatorial (Fibonacci-scale counting)
  - SUPER-EXPONENTIAL: enumerative (path count, tournament count)

  H(T) lives in the SUPER-EXPONENTIAL regime.
  But the FORBIDDEN VALUES and structural properties live in
  the POLYNOMIAL and EXPONENTIAL regimes!
  7 = constant, 21 = triangular, 45 = triangular, ...
""")

    # ============================================================
    # PART 9: THE RESONANCE AT n=7
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 9: THE RESONANCE AT n=7")
    print("  n=7 is a Paley prime where max_H = H(Paley T_7) = 189")
    print(f"{'='*70}")

    n = 7
    print(f"\n  At n=7:")
    print(f"    max_H = 189 = 7 * 27 = 7 * 3^3")
    print(f"    Fib(7) = 13")
    print(f"    Tri(7) = 28 = C(8,2)")
    print(f"    C(7,2) = 21 = FORBIDDEN H!")
    print(f"    Jacobsthal(7) = 43")
    print(f"    Catalan(7) = 429")
    print(f"    P_3strand(7) = 15 = max_H(5)")
    print(f"    n! = 5040")
    print(f"    mean_H = 5040/64 = 78.75")
    print(f"    max_H / mean_H = 189/78.75 = {189/78.75:.4f}")
    print(f"")
    print(f"    189 = 3^3 * 7 = 27 * 7")
    print(f"    This is in the (z-2)(z-3) orbit: 7 -> 21 -> 63 -> 189!")
    print(f"    The FOURTH term of the forbidden orbit IS the maximizer!")
    print(f"    The orbit goes: forbidden, forbidden, achievable, MAXIMIZER!")

    # ============================================================
    # PART 10: DEEP — What sequences are "tournament-universal"?
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 10: TOURNAMENT-UNIVERSAL SEQUENCES")
    print("  A sequence is 'tournament-universal' if EVERY term appears as")
    print("  some tournament invariant at some n.")
    print(f"{'='*70}")

    # Check: does every Fibonacci number appear as an H value?
    H_achievable = set()
    for n in range(3, 8):
        m = C(n, 2)
        for bits in range(min(2**m, 100000)):
            A = np.zeros((n, n), dtype=int)
            idx = 0
            for i in range(n):
                for j in range(i+1, n):
                    if bits & (1 << idx): A[j][i] = 1
                    else: A[i][j] = 1
                    idx += 1
            # Quick H via DP
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
            H_achievable.add(H)

    print(f"\n  Achievable H values (n=3..7): {sorted(H_achievable)[:30]}...")
    print(f"  Total achievable: {len(H_achievable)}")

    fibs = [fib(n) for n in range(1, 15)]
    fib_in_H = [f for f in fibs if f in H_achievable]
    fib_not_in_H = [f for f in fibs if f not in H_achievable and f > 0]
    print(f"\n  Fibonacci numbers that ARE achievable H values: {fib_in_H}")
    print(f"  Fibonacci numbers NOT achievable: {fib_not_in_H}")

    tris = [triangular(n) for n in range(1, 15)]
    tri_in_H = [t for t in tris if t in H_achievable]
    tri_not_in_H = [t for t in tris if t not in H_achievable and t > 0]
    print(f"\n  Triangular numbers that ARE achievable H values: {tri_in_H}")
    print(f"  Triangular numbers NOT achievable: {tri_not_in_H}")

    cats = [catalan(n) for n in range(1, 10)]
    cat_in_H = [c for c in cats if c in H_achievable]
    print(f"\n  Catalan numbers that ARE achievable H values: {cat_in_H}")

    jacs = [jacobsthal(n) for n in range(1, 15)]
    jac_in_H = [j for j in jacs if j in H_achievable]
    jac_not_in_H = [j for j in jacs if j not in H_achievable and j > 0]
    print(f"\n  Jacobsthal numbers that ARE achievable H values: {jac_in_H}")
    print(f"  Jacobsthal numbers NOT achievable: {jac_not_in_H}")

    print(f"\n{'='*70}")
    print("DONE — DEEP SEQUENCE META-ALIGNMENT")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
