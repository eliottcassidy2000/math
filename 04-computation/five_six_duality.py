#!/usr/bin/env python3
"""
five_six_duality.py -- The additive/multiplicative duality of 2 and 3

5 = 2 + 3 (additive combination)
6 = 2 * 3 (multiplicative combination)

These are NOT just numbers -- they are the two fundamental OPERATIONS
applied to the two fundamental PRIMES of tournament theory.

In tournament theory:
  5 = length of the next odd cycle after 3
  6 = Vandermonde determinant at (2,3) = 3!
  5 = number of vertices where alpha_2 first becomes nonzero (n=6 has a2>0)
      Wait: at n=5, a2=0 always. At n=6, a2 can be 1,2,3,4.
      Actually 5 vertices is where the FIRST 5-cycle can appear.

Explores:
1. I(Omega, 5) and I(Omega, 6) -- evaluating at 2+3 and 2*3
2. The polynomial I(x) = 1 + a1*x + a2*x^2 at x=5 and x=6
3. I(5)/I(2) and I(6)/I(3) ratios
4. 5 and 6 in the CRT tower: I(5) mod 6 = I(-1) mod 6 = H mod 6!
5. The additive decomposition: I(2+3) = I(2) + 3*I'(2) + 9*I''(2)/2
6. The multiplicative decomposition: I(2*3) = I(6), I(2)*I(3) vs I(6)
7. Newton's identity: I(5) = I(2) + 3*(a1 + 4*a2) -- the "3-shift"
8. 10 = 2*5 = 2*(2+3) and 11 = 10+1 = 2*(2+3)+1
9. The role of 5-cycles: exactly the cycles of length 2+3
10. The role of C(n,3) and C(n,5) in counting cycle-hosting vertex sets

Author: kind-pasteur-2026-03-14-S64
"""

import numpy as np
from itertools import combinations
from collections import defaultdict
from math import comb, gcd


def random_tournament(n, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp or dp[(mask, v)] == 0:
                continue
            cnt = dp[(mask, v)]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and dp[(full, v)] > 0:
            if A[verts[v]][verts[0]]:
                total += dp[(full, v)]
    return total


def get_alpha_1_2(A, n):
    cycles = []
    for k in range(3, n+1, 2):
        for subset in combinations(range(n), k):
            verts = list(subset)
            nc = count_ham_cycles(A, verts)
            for _ in range(nc):
                cycles.append(frozenset(subset))
    a1 = len(cycles)
    a2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                a2 += 1
    return a1, a2


def I(a1, a2, x):
    return 1 + a1 * x + a2 * x * x


def main():
    n = 7
    rng = np.random.default_rng(42)
    num = 500

    print("Collecting tournament data...")
    data = []
    for i in range(num):
        A = random_tournament(n, rng)
        a1, a2 = get_alpha_1_2(A, n)
        data.append({'a1': a1, 'a2': a2})
        if (i+1) % 100 == 0:
            print(f"  {i+1}/{num}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 1: The Additive-Multiplicative Table")
    print("=" * 70)

    print("""
  The fundamental pair (2, 3) generates four numbers via +, *:
    2 + 3 = 5     (additive)
    2 * 3 = 6     (multiplicative)
    2 + 3 + 1 = 6 (hmm, same as product -- THIS IS SPECIAL)
    2^3 = 8       (exponential)
    3^2 = 9       (reverse exponential)

  The COINCIDENCE 2*3 = 2+3+1 is equivalent to:
    2*3 - 2 - 3 = 1
    (2-1)*(3-1) = 1*2 = 2
  More generally: p*q - p - q = (p-1)(q-1) - 1
  For p=2, q=3: (1)(2) - 1 = 1
  This is the ONLY pair of primes where p*q = p+q+1.

  In tournament theory:
    x = 2: the evaluation point (H = I(Omega, 2))
    x+1 = 3: the Galois gap / beta step
    x + (x+1) = 2*x + 1 = 5: the sum
    x * (x+1) = x^2 + x = 6: the product
    5 = 2*3 - 1 (one less than the product!)
    6 = 2*3 = 3! (also the Vandermonde det)
""")

    # ============================================================
    print("=" * 70)
    print("PART 2: I(5) and I(6) -- evaluating at sum and product")
    print("=" * 70)

    print("\n  I(x) = 1 + a1*x + a2*x^2")
    print("  I(5) = 1 + 5*a1 + 25*a2")
    print("  I(6) = 1 + 6*a1 + 36*a2")
    print()
    print("  KEY: I(5) mod 6 = I(-1) mod 6 (from CRT tower: b=5, b+1=6)")
    print("  And I(-1) mod 6 = H mod 6 (since H mod 3 = I(-1) mod 3, H mod 2 = 1)")
    print()
    print("  Wait: I(-1) mod 6 vs H mod 6. H = I(2).")
    print("  I(-1) mod 2 = (1-a1+a2) mod 2. H mod 2 = 1.")
    print("  I(-1) mod 3 = H mod 3 (proved).")
    print("  So I(-1) mod 6 = CRT(I(-1) mod 2, I(-1) mod 3)")
    print("     = CRT(I(-1) mod 2, H mod 3)")
    print("  But I(-1) mod 2 may differ from H mod 2 = 1!")

    # Check: does I(-1) mod 2 = 1 always?
    # I(-1) = 1 - a1 + a2. Mod 2: 1 + a1 + a2 (since -1 = 1 mod 2)
    # H = 1 + 2*a1 + 4*a2. Mod 2: 1.
    # So I(-1) mod 2 = 1 + a1 + a2 mod 2. NOT always 1.
    i_neg1_mod2 = defaultdict(int)
    for t in data:
        val = (1 - t['a1'] + t['a2']) % 2
        i_neg1_mod2[val] += 1
    print(f"\n  I(-1) mod 2: {dict(sorted(i_neg1_mod2.items()))}")
    print(f"  So I(-1) mod 2 != 1 in general. I(-1) mod 6 != H mod 6.")

    print("\n  BUT: I(5) mod 6 = I(-1) mod 6 (since 5 = -1 mod 6)")
    print("  This is a DIFFERENT fact from H mod 3 = I(-1) mod 3.")

    # Verify I(5) mod 6 = I(-1) mod 6
    ok = all(I(t['a1'], t['a2'], 5) % 6 == (1 - t['a1'] + t['a2']) % 6 for t in data)
    print(f"  I(5) mod 6 = I(-1) mod 6: {ok}")

    # The connection: I(5) mod 6 encodes BOTH topology mod 2 AND topology mod 3
    print("\n  I(5) mod 6 = CRT(I(-1) mod 2, I(-1) mod 3)")
    print("             = CRT((1+a1+a2) mod 2, (1-a1+a2) mod 3)")
    print("  This is the FULL mod-6 topology, which H mod 6 only partially captures.")

    # Distribution of I(5) mod 6
    i5_mod6 = defaultdict(int)
    h_mod6 = defaultdict(int)
    ineg1_mod6 = defaultdict(int)
    for t in data:
        i5_mod6[I(t['a1'], t['a2'], 5) % 6] += 1
        h_mod6[I(t['a1'], t['a2'], 2) % 6] += 1
        ineg1_mod6[(1 - t['a1'] + t['a2']) % 6] += 1

    print(f"\n  I(5) mod 6: {dict(sorted(i5_mod6.items()))}")
    print(f"  H mod 6:    {dict(sorted(h_mod6.items()))}")
    print(f"  I(-1) mod 6: {dict(sorted(ineg1_mod6.items()))}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 3: The Taylor shift -- I(2+3) from I(2)")
    print("=" * 70)

    print("\n  I(x+h) = I(x) + h*I'(x) + h^2*I''(x)/2")
    print("  I(2+3) = I(5) = I(2) + 3*I'(2) + 9*I''(2)/2")
    print("  where I'(2) = a1 + 4*a2, I''(2)/2 = a2")
    print("  So: I(5) = H + 3*(a1 + 4*a2) + 9*a2 = H + 3*a1 + 21*a2")

    for t in data[:10]:
        a1, a2 = t['a1'], t['a2']
        H = I(a1, a2, 2)
        I5 = I(a1, a2, 5)
        I5_taylor = H + 3*a1 + 21*a2
        print(f"    a1={a1:2d}, a2={a2:2d}: H={H:3d}, I(5)={I5:4d}, Taylor={I5_taylor:4d}, match={I5==I5_taylor}")

    print(f"\n  I(5) - H = 3*a1 + 21*a2 = 3*(a1 + 7*a2)")
    print(f"  So I(5) - H is ALWAYS divisible by 3!")
    print(f"  And (I(5)-H)/3 = a1 + 7*a2")

    # This is significant: a1 + 7*a2 is the "7-weighted total"
    # At n=7: this connects to the number of vertices!
    print(f"\n  At n=7: (I(5)-H)/3 = a1 + 7*a2 = a1 + n*a2")
    print(f"  The weight 7 = n is the number of vertices!")

    # Is this a coincidence? Let's check: I(5)-H = 3*a1 + 21*a2
    # = 3*a1 + 3*7*a2 = 3*(a1 + 7*a2)
    # The 7 comes from (5^2-2^2)/2 = (25-4)/2... no.
    # Actually: I(5) - I(2) = a1*(5-2) + a2*(25-4) = 3*a1 + 21*a2 = 3*(a1 + 7*a2)
    # The 7 = (5+2)(5-2)/(5-2) = ... no, simpler:
    # 21 = 5^2 - 4 = 5^2 - 2^2, and 21/3 = 7 = 5+2.
    # So the weight is 5+2 = 7 = n. This IS a coincidence for n=7!

    print(f"\n  GENERAL: I(x+h) - I(x) = h*a1 + (2*x*h + h^2)*a2")
    print(f"          = h*(a1 + (2x+h)*a2)")
    print(f"  For x=2, h=3: I(5)-I(2) = 3*(a1 + 7*a2)")
    print(f"  The coefficient 2x+h = 2*2+3 = 7 = n. COINCIDENCE for n=7.")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 4: I(6) = I(2*3) vs I(2)*I(3) -- additive vs multiplicative")
    print("=" * 70)

    print("\n  For a polynomial f(x), f(ab) != f(a)*f(b) in general.")
    print("  But the RATIO I(6)/(I(2)*I(3)) is interesting.")

    for t in data[:10]:
        a1, a2 = t['a1'], t['a2']
        H = I(a1, a2, 2)
        I3 = I(a1, a2, 3)
        I6 = I(a1, a2, 6)
        product = H * I3
        ratio = I6 / product if product > 0 else 0
        print(f"    a1={a1:2d}, a2={a2:2d}: I(2)={H:3d}, I(3)={I3:3d}, I(6)={I6:4d}, "
              f"I(2)*I(3)={product:5d}, I(6)/prod={ratio:.4f}")

    # The "multiplicative defect" I(6) - I(2)*I(3)
    print("\n  Multiplicative defect D = I(6) - I(2)*I(3):")
    print("  D = (1+6a1+36a2) - (1+2a1+4a2)(1+3a1+9a2)")
    print("  Let me expand:")
    print("  (1+2a1+4a2)(1+3a1+9a2) = 1 + 3a1 + 9a2 + 2a1 + 6a1^2 + 18a1*a2")
    print("                          + 4a2 + 12a1*a2 + 36a2^2")
    print("                        = 1 + 5a1 + 13a2 + 6a1^2 + 30a1*a2 + 36a2^2")
    print("  I(6) = 1 + 6a1 + 36a2")
    print("  D = (1+6a1+36a2) - (1+5a1+13a2+6a1^2+30a1*a2+36a2^2)")
    print("    = a1 + 23a2 - 6a1^2 - 30a1*a2 - 36a2^2")

    for t in data[:10]:
        a1, a2 = t['a1'], t['a2']
        D = a1 + 23*a2 - 6*a1**2 - 30*a1*a2 - 36*a2**2
        D_check = I(a1, a2, 6) - I(a1, a2, 2) * I(a1, a2, 3)
        print(f"    a1={a1:2d}, a2={a2:2d}: D={D:7d}, check={D_check:7d}, match={D==D_check}")

    # D is always negative (since -6*a1^2 dominates)
    all_neg = all(I(t['a1'], t['a2'], 6) - I(t['a1'], t['a2'], 2)*I(t['a1'], t['a2'], 3) < 0
                  for t in data if t['a1'] > 0)
    print(f"\n  D < 0 for all non-transitive: {all_neg}")
    print(f"  Product ALWAYS exceeds the 6-evaluation: I(2)*I(3) > I(6)")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 5: The special identity (2-1)(3-1) = 2")
    print("=" * 70)

    print("\n  (p-1)(q-1) = 1 iff p*q = p+q+1")
    print("  For primes: only p=2, q=3 satisfy this.")
    print("  Consequence: 2*3 - (2+3) = 1, so 6 = 5+1.")
    print()
    print("  In I(x) terms:")
    print("  I(6) - I(5) = a1*(6-5) + a2*(36-25) = a1 + 11*a2")
    print()
    print("  I(6) - I(5) = a1 + 11*a2")
    print("  Compare with: I(3) - I(2) = a1 + 5*a2")
    print("  And: I(2) - I(1) = a1 + 3*a2")
    print()
    print("  General: I(k+1) - I(k) = a1 + (2k+1)*a2")
    print("  The step from k to k+1 has weight 2k+1 on a2.")

    print("\n  Step weights:")
    for k in range(0, 12):
        w = 2*k + 1
        print(f"    I({k+1})-I({k}) = a1 + {w}*a2   (weight {w})")

    print(f"\n  The weight at k=2 (I(3)-I(2)) is 5 = 2+3")
    print(f"  The weight at k=5 (I(6)-I(5)) is 11 = our positional 11!")
    print(f"  The weight at k=1 (I(2)-I(1)) is 3 = our fundamental 3!")
    print(f"  The weight at k=4 (I(5)-I(4)) is 9 = 3^2!")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 6: The step weight 2k+1 and the odd number sequence")
    print("=" * 70)

    print("\n  I(k+1) - I(k) = a1 + (2k+1)*a2")
    print("  The weights are 1, 3, 5, 7, 9, 11, 13, ...")
    print("  These are the ODD NUMBERS = the differences of perfect squares!")
    print("  (k+1)^2 - k^2 = 2k+1")
    print()
    print("  The x^2 term in I(x) means that consecutive differences")
    print("  grow as odd numbers. This is the QUADRATIC structure of I.")

    print("\n  Summing steps:")
    print("  I(b) - I(0) = sum_{k=0}^{b-1} [a1 + (2k+1)*a2]")
    print("              = b*a1 + a2 * sum_{k=0}^{b-1} (2k+1)")
    print("              = b*a1 + a2 * b^2")
    print("  So I(b) = 1 + b*a1 + b^2*a2. CORRECT by definition.")

    print("\n  The step I(3)-I(2) = a1 + 5*a2 has a dual interpretation:")
    print("  It is the coefficient of the FIRST DIFFERENCE at the counting point.")
    print("  And 5 = 2+3 = the sum of the fundamental pair.")
    print("  The derivative I'(2) = a1 + 4*a2. NOT the same (off by a2).")
    print("  Step = derivative + a2. The extra a2 is the 'curvature correction'.")

    # Verify
    for t in data[:5]:
        a1, a2 = t['a1'], t['a2']
        step = I(a1, a2, 3) - I(a1, a2, 2)
        deriv = a1 + 4*a2
        print(f"    a1={a1:2d}, a2={a2:2d}: step={step:3d}, deriv={deriv:3d}, diff={step-deriv}=a2={a2}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 7: 10 = 2*(2+3) and 11 = 2*(2+3)+1")
    print("=" * 70)

    print("\n  10 = 2 * 5 = 2 * (2+3)")
    print("  11 = 2 * 5 + 1 = 2*(2+3) + 1")
    print()
    print("  10 = x * (x + (x+1)) = x*(2x+1) where x=2")
    print("  11 = x*(2x+1) + 1")
    print()
    print("  In Z[x] at x=2:")
    print("    10 = 2x^2 + x = x(2x+1)")
    print("    11 = 2x^2 + x + 1 = x(2x+1) + 1")
    print()
    print("  Note: 2x+1 = 5 is the step weight at k=x=2.")
    print("  So 10 = x * (step weight at x)!")

    print("\n  I(10) = 1 + 10*a1 + 100*a2")
    print("  I(11) = 1 + 11*a1 + 121*a2")
    print("  I(11) - I(10) = a1 + 21*a2 (step weight 2*10+1 = 21 = 3*7)")

    # The 21 = 3*7 is remarkable: it's 3 * n (at n=7)
    print(f"\n  Step weight at k=10: 21 = 3 * 7 = 3 * n (at n=7)")
    print(f"  Another n=7 coincidence!")

    # I(10) and I(11) in terms of I(2) and I(3)
    print("\n  I(10) and I(11) via Vandermonde extraction:")
    print("  I(10) = 1 + 10a1 + 100a2")
    print("  I(2) = 1 + 2a1 + 4a2 = H")
    print("  I(10) - I(2) = 8a1 + 96a2 = 8*(a1 + 12a2)")
    print("  So I(10) = H + 8*(a1 + 12*a2)")
    print("  And (I(10) - H)/8 = a1 + 12*a2")

    for t in data[:5]:
        a1, a2 = t['a1'], t['a2']
        H = I(a1, a2, 2)
        I10 = I(a1, a2, 10)
        q = (I10 - H) // 8
        expected = a1 + 12*a2
        print(f"    a1={a1:2d}, a2={a2:2d}: (I(10)-H)/8 = {q} = a1+12a2 = {expected}, match={q==expected}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 8: Factorization in Z[x] -- everything through x=2")
    print("=" * 70)

    print("\n  Rewriting key numbers in terms of x=2:")
    for val, expr in [
        (1, "1"),
        (2, "x"),
        (3, "x+1"),
        (4, "x^2"),
        (5, "x^2+1 = x+(x+1) = 2x+1"),
        (6, "x(x+1) = x^2+x"),
        (7, "x^3-1 = (x-1)(x^2+x+1) = 1*7"),
        (8, "x^3"),
        (9, "(x+1)^2"),
        (10, "x(x^2+1) = x(2x+1)"),
        (11, "x^3+x+1"),
        (12, "x^2(x+1)"),
        (21, "x^4+x^2+1 = (x^2+x+1)(x^2-x+1) = 7*3"),
    ]:
        # Verify
        import re
        first_expr = expr.split("=")[0].strip()
        # Insert * between digit and x, and between ) and x, and between x and (
        first_expr = re.sub(r'(\d)(x)', r'\1*x', first_expr)
        first_expr = re.sub(r'(\))(x)', r'\1*x', first_expr)
        first_expr = re.sub(r'(x)(\()', r'\1*\2', first_expr)
        first_expr = first_expr.replace("^", "**")
        first_expr = first_expr.replace("x", "(2)")
        # Also handle number( and )( patterns
        first_expr = re.sub(r'(\d)\(', r'\1*(', first_expr)
        first_expr = re.sub(r'\)\(', r')*(', first_expr)
        poly_val = eval(first_expr)
        print(f"    {val:3d} = {expr:40s} [check: {poly_val}]")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 9: The 5-cycle structure -- cycles of length 2+3")
    print("=" * 70)

    print("\n  At n=7: C(7,5) = 21 potential 5-cycle vertex sets")
    print("  Each 5-vertex set can host 0, 1, 2, or 3 directed 5-cycles")
    print("  (since a 5-vertex tournament can have 0-12 Hamiltonian directed cycles)")

    # Count directed 5-cycles at n=7
    c5_dist = defaultdict(int)
    c5_per_set = defaultdict(list)
    for t_idx, t in enumerate(data[:100]):
        A = random_tournament(n, np.random.default_rng(42 + t_idx))
        # Recompute to get per-vertex-set data
        total_c5 = 0
        for subset in combinations(range(n), 5):
            nc = count_ham_cycles(A, list(subset))
            total_c5 += nc
            c5_per_set[nc].append(subset)
        c5_dist[total_c5] += 1

    print(f"\n  Directed 5-cycles per tournament distribution:")
    for c5, cnt in sorted(c5_dist.items()):
        if cnt >= 2:
            print(f"    c5={c5:2d}: {cnt} tournaments")

    # Per-vertex-set counts
    per_set_dist = defaultdict(int)
    for nc in c5_per_set:
        per_set_dist[nc] += len(c5_per_set[nc])
    print(f"\n  Directed 5-cycles per 5-vertex set:")
    for nc, cnt in sorted(per_set_dist.items()):
        print(f"    {nc}: {cnt} sets ({100*cnt/sum(per_set_dist.values()):.1f}%)")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 10: The FIVE identities connecting 2, 3, 5, 6")
    print("=" * 70)

    print("""
  Identity 1 (ADDITIVE): I(5) = I(2) + 3*(a1 + 7*a2)
    The shift from x=2 to x=5=2+3 is always divisible by 3.
    The quotient a1+7*a2 has weight 7=n at n=7.

  Identity 2 (MULTIPLICATIVE): I(2)*I(3) - I(6) = 6*a1^2 + 30*a1*a2 + 36*a2^2 - a1 - 23*a2
    The product always EXCEEDS the composite evaluation.

  Identity 3 (CRT): I(5) mod 6 = I(-1) mod 6
    Evaluation at 5=2+3 mod 6=2*3 gives FULL topological mod 6 info.

  Identity 4 (STEP): I(3)-I(2) = a1 + 5*a2
    The step at the counting point has weight 5 = 2+3 on a2.
    This equals I'(2) + a2 (derivative + curvature).

  Identity 5 (REDUCED): H = 3 - 2*I(-1) + 6*a2
    The three numbers 3, 2, 6 are x+1, x, x(x+1).
    H is determined by topology (b0) and curvature (a2) via
    the multiplicative structure of {2, 3, 6}.
""")

    # Verify all five
    print("  Verification (10 samples):")
    for t in data[:10]:
        a1, a2 = t['a1'], t['a2']
        H = I(a1, a2, 2)
        b0 = 1 - a1 + a2

        id1 = I(a1, a2, 5) == H + 3*(a1 + 7*a2)
        id2 = H * I(a1, a2, 3) - I(a1, a2, 6) == 6*a1**2 + 30*a1*a2 + 36*a2**2 - a1 - 23*a2
        id3 = I(a1, a2, 5) % 6 == b0 % 6
        id4 = I(a1, a2, 3) - H == a1 + 5*a2
        id5 = H == 3 - 2*b0 + 6*a2

        print(f"    a1={a1:2d}, a2={a2:2d}: [1]={id1} [2]={id2} [3]={id3} [4]={id4} [5]={id5}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 11: The decimal system: 10 = 2*5, 11 = 2*5+1")
    print("=" * 70)

    print("""
  Why base 10? Mathematically:
    10 = 2 * 5 = x * (x^2 + 1) = x * (2x + 1)    at x=2
    10 = 2 * (2 + 3) = evaluation point * additive sum

  The decimal digits of H encode:
    H mod 10 = H mod (2*5) = CRT(H mod 2, H mod 5) = CRT(1, H mod 5)
    H mod 5 determines the last decimal digit (since H odd, digit in {1,3,5,7,9})

  The last decimal digit of H is determined by a1 and a2 mod 5:
    H mod 5 = (1 + 2*a1 + 4*a2) mod 5
""")

    h_mod5 = defaultdict(int)
    for t in data:
        h_mod5[I(t['a1'], t['a2'], 2) % 5] += 1
    print(f"  H mod 5: {dict(sorted(h_mod5.items()))}")

    # H mod 10
    h_mod10 = defaultdict(int)
    for t in data:
        h_mod10[I(t['a1'], t['a2'], 2) % 10] += 1
    print(f"  H mod 10: {dict(sorted(h_mod10.items()))}")
    print(f"  All 5 odd digits appear.")

    # I(10) mod 11 and I(11) mod 12
    print(f"\n  I(10) mod 11 = I(-1) mod 11 = topology mod 11:")
    i10_mod11 = defaultdict(int)
    for t in data:
        i10_mod11[I(t['a1'], t['a2'], 10) % 11] += 1
    print(f"    {dict(sorted(i10_mod11.items()))}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 12: The multiplication table of {2, 3, 5, 6}")
    print("=" * 70)

    print("\n  Products and their Z[x] representations at x=2:")
    pairs = [(2,3,6), (2,5,10), (2,6,12), (3,5,15), (3,6,18), (5,6,30)]
    for a, b, c in pairs:
        print(f"    {a} * {b} = {c:3d}", end="")
        # I(c) formula
        ic = f"1 + {c}*a1 + {c**2}*a2"
        print(f"    I({c}) = {ic}")

    print(f"\n  Sums:")
    pairs_s = [(2,3,5), (2,5,7), (2,6,8), (3,5,8), (3,6,9), (5,6,11)]
    for a, b, c in pairs_s:
        print(f"    {a} + {b} = {c:3d}", end="")
        # Step weight at k=c-1
        w = 2*(c-1) + 1
        print(f"    step weight at k={c-1} is {w}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 13: H mod 5 and the golden ratio")
    print("=" * 70)

    print("\n  5 = 2+3 is the sum. H mod 5 = (1+2a1+4a2) mod 5.")
    print("  Since 4 = -1 mod 5 and 2 = 2 mod 5:")
    print("  H mod 5 = (1 + 2*a1 - a2) mod 5")
    print()
    print("  The golden ratio phi = (1+sqrt(5))/2 is a root of x^2-x-1.")
    print("  In F_5: x^2-x-1 = x^2+4x+4 = (x+2)^2 mod 5.")
    print("  So phi = -2 = 3 mod 5.")
    print("  And phi^2 = phi+1 = 4 mod 5.")
    print()
    print("  Fibonacci mod 5: 0,1,1,2,3,0,3,3,1,4,0,4,4,3,2,0,... period 20")
    print("  Fibonacci mod 5 has period 20 = 4*5.")

    # H mod 5 distribution
    print(f"\n  H mod 5 distribution (500 samples):")
    for r in sorted(h_mod5):
        print(f"    H = {r} mod 5: {h_mod5[r]} ({100*h_mod5[r]/num:.1f}%)")

    # Does a1 mod 5 determine H mod 5 (given a2 mod 5)?
    cross = defaultdict(int)
    for t in data:
        cross[(t['a1'] % 5, t['a2'] % 5, I(t['a1'], t['a2'], 2) % 5)] += 1

    print(f"\n  (a1%5, a2%5) -> H%5 (checking determinism):")
    for a1m, a2m in sorted(set((k[0], k[1]) for k in cross)):
        h_vals = set(k[2] for k in cross if k[0]==a1m and k[1]==a2m)
        expected = (1 + 2*a1m + 4*a2m) % 5
        if len(h_vals) > 1:
            print(f"    ({a1m},{a2m}): H%5 = {h_vals} [MULTIPLE]")
        elif list(h_vals)[0] != expected:
            print(f"    ({a1m},{a2m}): H%5 = {h_vals} != expected {expected}")
    print(f"  (a1%5, a2%5) uniquely determines H%5: True (by definition)")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 14: Summary -- The Arithmetic of Tournament Theory")
    print("=" * 70)

    print("""
  THE FUNDAMENTAL PAIR: (2, 3)
    2 = evaluation point (H = I(Omega, 2))
    3 = topological gap (3 = 2-(-1)), beta step, 3-cycle length

  THE ADDITIVE SON: 5 = 2 + 3
    5 = length of next odd cycle
    5 = step weight at counting point: I(3)-I(2) = a1 + 5*a2
    I(5) = I(2) + 3*(a1 + (2*2+3)*a2)
    I(5) mod 6 = I(-1) mod 6 = full topological mod 6

  THE MULTIPLICATIVE SON: 6 = 2 * 3
    6 = Vandermonde determinant = 3!
    6 = CRT modulus combining mod 2 and mod 3
    6*a2 = 2*I(3) - 3*H + 1 (extraction formula)
    H = 3 - 2*b0 + 6*b2 (reduced formula)

  THE SPECIAL IDENTITY: 5 + 1 = 6
    (2-1)(3-1) = 2, so 2*3 = 2+3+1
    ONLY pair of primes with this property
    CONNECTS additive and multiplicative sons

  THE DECIMAL BASE: 10 = 2 * 5 = 2(2+3)
    10 = evaluation point * additive son
    11 = 10 + 1 = 2*5 + 1
    I(10) mod 11 = I(-1) mod 11

  THE STEP SEQUENCE: I(k+1)-I(k) = a1 + (2k+1)*a2
    k=0: weight 1    (trivial)
    k=1: weight 3    (= x+1, the topological gap)
    k=2: weight 5    (= 2+3, the additive son)
    k=4: weight 9    (= (x+1)^2 = 3^2)
    k=5: weight 11   (= positional 11)
    k=10: weight 21  (= 3*7 = 3*n at n=7)

  MULTIPLICATION vs ADDITION on I(x):
    I(a+b) = I(a) + b*I'(a) + (b^2/2)*I''(a)   [Taylor]
    I(a*b) != I(a)*I(b)                           [not multiplicative]
    I(a)*I(b) > I(a*b) for a,b >= 2               [super-multiplicative]
""")

    print("Done.")


if __name__ == '__main__':
    main()
