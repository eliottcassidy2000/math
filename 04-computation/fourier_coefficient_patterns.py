#!/usr/bin/env python3
"""
Fourier Coefficient Patterns in W(r) — Analysis of constant terms, t3 coefficients,
alternating structure, and connections to classical number sequences.

Tasks:
  1. Verify w_k formulas at n=9 by sampling random tournaments
  2. Investigate constant terms — Bernoulli/Euler/Eulerian connection
  3. Trace t3 coefficients across Fourier degrees
  4. Verify w_{n-3} = 2*(n-2)! * (t3 - C(n,3)/4)
  5. Alternating sign structure and Euler zigzag numbers

kind-pasteur-2026-03-07
"""

from itertools import combinations, permutations
from math import comb, factorial, gcd
from fractions import Fraction
import random
import sys


# ─────────────────────────────────────────────────────────────────────
# Core tournament utilities (exact arithmetic with Fraction where needed)
# ─────────────────────────────────────────────────────────────────────

def random_tournament(n, seed):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def compute_W_poly(A, n):
    """
    Compute W(r) as exact polynomial coefficients using Fraction arithmetic.
    W(r) = sum over Hamiltonian paths of product_{consecutive (u,v)} (r + A[u][v] - 1/2).
    Returns dict: degree -> Fraction coefficient.
    """
    half = Fraction(1, 2)
    # dp[(mask, v)] = polynomial as dict {degree: Fraction}
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = {0: Fraction(1)}

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            poly_v = dp.get((mask, v))
            if poly_v is None:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                # Multiply by (r + A[v][u] - 1/2)
                c = Fraction(A[v][u]) - half  # constant part
                new_mask = mask | (1 << u)
                key = (new_mask, u)
                if key not in dp:
                    dp[key] = {}
                target = dp[key]
                for deg, coeff in poly_v.items():
                    # coeff * r -> degree deg+1
                    target[deg + 1] = target.get(deg + 1, Fraction(0)) + coeff
                    # coeff * c -> degree deg
                    target[deg] = target.get(deg, Fraction(0)) + coeff * c

    full = (1 << n) - 1
    result = {}
    for v in range(n):
        poly_v = dp.get((full, v), {})
        for deg, coeff in poly_v.items():
            result[deg] = result.get(deg, Fraction(0)) + coeff
    return result


def count_H(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[(mask | (1 << u), u)] = dp.get((mask | (1 << u), u), 0) + c
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_cycles(A, n, L):
    """Count directed cycles of length L (counting each cycle once per vertex set,
    but summing over all cycle orientations)."""
    total = 0
    for verts in combinations(range(n), L):
        sub = [[A[verts[i]][verts[j]] for j in range(L)] for i in range(L)]
        dp = [[0] * L for _ in range(1 << L)]
        dp[1][0] = 1
        for m in range(1, 1 << L):
            for v in range(L):
                if not (m & (1 << v)) or dp[m][v] == 0:
                    continue
                for u in range(L):
                    if m & (1 << u):
                        continue
                    if sub[v][u]:
                        dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << L) - 1
        total += sum(dp[full][v] for v in range(1, L) if sub[v][0])
    return total


def count_alpha_disjoint(A, n, k):
    """Count k-tuples of pairwise vertex-disjoint directed 3-cycles."""
    tri = []
    for i, j, l in combinations(range(n), 3):
        if A[i][j] and A[j][l] and A[l][i]:
            tri.append(frozenset([i, j, l]))
        if A[i][l] and A[l][j] and A[j][i]:
            tri.append(frozenset([i, j, l]))
    if k == 2:
        count = 0
        for a in range(len(tri)):
            for b in range(a + 1, len(tri)):
                if len(tri[a] & tri[b]) == 0:
                    count += 1
        return count
    elif k == 3:
        count = 0
        for a in range(len(tri)):
            for b in range(a + 1, len(tri)):
                if tri[a] & tri[b]:
                    continue
                for c in range(b + 1, len(tri)):
                    if not (tri[a] & tri[c]) and not (tri[b] & tri[c]):
                        count += 1
        return count
    return 0


def count_alpha_35(A, n):
    """Count disjoint (3-cycle, 5-cycle) pairs."""
    tri = []
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            tri.append(frozenset([i, j, k]))
        if A[i][k] and A[k][j] and A[j][i]:
            tri.append(frozenset([i, j, k]))
    pent = []
    for verts in combinations(range(n), 5):
        sub = [[A[verts[i]][verts[j]] for j in range(5)] for i in range(5)]
        dp = [[0] * 5 for _ in range(1 << 5)]
        dp[1][0] = 1
        for m in range(1, 1 << 5):
            for v in range(5):
                if not (m & (1 << v)) or dp[m][v] == 0:
                    continue
                for u in range(5):
                    if m & (1 << u):
                        continue
                    if sub[v][u]:
                        dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << 5) - 1
        cc = sum(dp[full][v] for v in range(1, 5) if sub[v][0])
        if cc > 0:
            pent.append((frozenset(verts), cc))
    count = 0
    for t in tri:
        for pv, pc in pent:
            if len(t & pv) == 0:
                count += pc
    return count


# ─────────────────────────────────────────────────────────────────────
# TASK 1: Verify n=9 Fourier coefficients on random tournaments
# ─────────────────────────────────────────────────────────────────────

def formula_w_coeffs_n9(t3, t5, t7, t9, a33, a35, a333):
    """Return dict {degree: Fraction} for each w_k at n=9 using the exact formulas."""
    w = {}
    w[8] = Fraction(362880)
    w[6] = Fraction(-211680) + Fraction(10080) * Fraction(t3)
    w[4] = Fraction(40320) - Fraction(4200) * Fraction(t3) + Fraction(240) * Fraction(t5) + Fraction(480) * Fraction(a33)
    w[2] = Fraction(-2640) + Fraction(462) * Fraction(t3) - Fraction(60) * Fraction(t5) + Fraction(12) * Fraction(t7) \
         - Fraction(120) * Fraction(a33) + Fraction(24) * Fraction(a35) + Fraction(48) * Fraction(a333)
    w[0] = Fraction(31) - Fraction(17, 2) * Fraction(t3) + Fraction(2) * Fraction(t5) - Fraction(1) * Fraction(t7) \
         + Fraction(4) * Fraction(a33) - Fraction(2) * Fraction(a35) - Fraction(4) * Fraction(a333) + Fraction(2) * Fraction(t9)
    return w


def task1_verify_n9(N=20):
    print("=" * 70)
    print("TASK 1: Verify n=9 Fourier coefficients on 20 random tournaments")
    print("=" * 70)
    n = 9
    all_ok = True
    for seed in range(N):
        A = random_tournament(n, seed)
        # Brute-force polynomial
        poly = compute_W_poly(A, n)

        # Compute invariants
        t3 = count_directed_cycles(A, n, 3)
        t5 = count_directed_cycles(A, n, 5)
        t7 = count_directed_cycles(A, n, 7)
        t9 = count_directed_cycles(A, n, 9)
        a33 = count_alpha_disjoint(A, n, 2)
        a35 = count_alpha_35(A, n)
        a333 = count_alpha_disjoint(A, n, 3)

        # Formula
        wf = formula_w_coeffs_n9(t3, t5, t7, t9, a33, a35, a333)

        ok = True
        for deg in [0, 2, 4, 6, 8]:
            bf = poly.get(deg, Fraction(0))
            fm = wf[deg]
            if bf != fm:
                ok = False
                print(f"  MISMATCH seed={seed} deg={deg}: brute_force={bf} formula={fm}")

        if ok:
            H = count_H(A, n)
            print(f"  seed={seed:2d}: H={H:6d}, t3={t3:3d}, t5={t5:4d}, t7={t7:5d}, t9={t9:6d}, "
                  f"a33={a33:3d}, a35={a35:4d}, a333={a333:3d} -- OK")
        else:
            all_ok = False

    print(f"\nAll {N} samples match: {all_ok}")
    return all_ok


# ─────────────────────────────────────────────────────────────────────
# TASK 2: Investigate constant terms
# ─────────────────────────────────────────────────────────────────────

def task2_constant_terms():
    print("\n" + "=" * 70)
    print("TASK 2: Constant terms in Fourier coefficients")
    print("=" * 70)

    # Known constant terms from the exact formulas:
    # n=3: w_2 = 6 = 3!, w_0 = -1/2
    # n=5: w_4 = 120 = 5!, w_2 = -30, w_0 = 1
    # n=7: w_6 = 5040 = 7!, w_4 = -2100, w_2 = 231, w_0 = -17/4
    # n=9: w_8 = 362880 = 9!, w_6 = -211680, w_4 = 40320, w_2 = -2640, w_0 = 31

    # Organize: const_term[n][k] = constant in w_k
    const = {
        3: {2: Fraction(6), 0: Fraction(-1, 2)},
        5: {4: Fraction(120), 2: Fraction(-30), 0: Fraction(1)},
        7: {6: Fraction(5040), 4: Fraction(-2100), 2: Fraction(231), 0: Fraction(-17, 4)},
        9: {8: Fraction(362880), 6: Fraction(-211680), 4: Fraction(40320), 2: Fraction(-2640), 0: Fraction(31)},
    }

    print("\nConstant terms table:")
    print(f"{'n':>3} | {'w_{n-1}':>12} | {'w_{n-3}':>12} | {'w_{n-5}':>12} | {'w_{n-7}':>12} | {'w_{n-9}':>12}")
    print("-" * 75)
    for n in [3, 5, 7, 9]:
        row = f"{n:3d} |"
        for j in range(5):
            k = n - 1 - 2 * j
            if k >= 0 and k in const[n]:
                row += f" {str(const[n][k]):>12s} |"
            else:
                row += f" {'':>12s} |"
        print(row)

    # Check: w_{n-1} = n!
    print("\nw_{n-1} = n! check:")
    for n in [3, 5, 7, 9]:
        print(f"  n={n}: {const[n][n-1]} = {factorial(n)}! = {factorial(n)}  {'OK' if const[n][n-1] == factorial(n) else 'FAIL'}")

    # Check: constant in w_{n-3} = -2*(n-2)! * C(n,3)/4 = -(n-2)! * C(n,3) / 2
    print("\nw_{n-3} constant = -2*(n-2)! * C(n,3)/4 = -(n-2)!*C(n,3)/2:")
    for n in [3, 5, 7, 9]:
        expected = -Fraction(factorial(n-2) * comb(n, 3), 2)
        actual = const[n][n - 3]
        print(f"  n={n}: expected={expected}, actual={actual}  {'OK' if expected == actual else 'FAIL'}")

    # Now look for patterns in all constant terms
    # Euler numbers E_n: 1, 0, -1, 0, 5, 0, -61, 0, 1385, ...
    # Tangent numbers T_n: 1, 2, 16, 272, 7936, ...
    # Bernoulli numbers B_n: 1, -1/2, 1/6, 0, -1/30, ...
    # Eulerian numbers A(n,k)

    print("\n--- Exploring ratios and patterns ---")

    # Look at the LEADING constant w_{n-1} = n!
    # Second constant w_{n-3}:
    print("\nw_{n-3} / (n-2)!:")
    for n in [3, 5, 7, 9]:
        ratio = const[n][n-3] / factorial(n-2)
        print(f"  n={n}: {ratio} = {float(ratio):.6f}")

    # The ratio is -C(n,3)/2: n=3: -1/2, n=5: -5, n=7: -35/2, n=9: -42
    # Check: -C(n,3)/2
    print("\n  Expected -C(n,3)/2:")
    for n in [3, 5, 7, 9]:
        print(f"    n={n}: -C({n},3)/2 = {-Fraction(comb(n,3), 2)}")

    # Third constant w_{n-5}:
    print("\nw_{n-5} constant (j=2):")
    for n in [5, 7, 9]:
        val = const[n][n-5]
        print(f"  n={n}: {val}")

    # n=5: w_0 = 1; n=7: w_2 = 231; n=9: w_4 = 40320
    # Ratios to (n-4)!:
    print("\nw_{n-5} / (n-4)!:")
    for n in [5, 7, 9]:
        ratio = const[n][n-5] / factorial(n-4)
        print(f"  n={n}: {ratio} = {float(ratio):.6f}")
    # n=5: 1/1=1, n=7: 231/6=77/2, n=9: 40320/120=336

    # Check connection to Bernoulli numbers
    # B_0=1, B_1=-1/2, B_2=1/6, B_4=-1/30, B_6=1/42
    print("\n--- Bernoulli number connection ---")
    bernoulli = {
        0: Fraction(1), 1: Fraction(-1, 2), 2: Fraction(1, 6),
        4: Fraction(-1, 30), 6: Fraction(1, 42), 8: Fraction(-1, 30),
        10: Fraction(5, 66), 12: Fraction(-691, 2730)
    }

    # w_0 values: -1/2, 1, -17/4, 31
    print("\nw_0 values across n:")
    w0_vals = {3: Fraction(-1, 2), 5: Fraction(1), 7: Fraction(-17, 4), 9: Fraction(31)}
    for n, v in w0_vals.items():
        print(f"  n={n}: w_0 = {v} = {float(v):.6f}")

    # Check: are these related to Euler numbers?
    # E_0=1, E_1=0, E_2=-1, E_3=0, E_4=5, E_5=0, E_6=-61, ...
    # Tangent numbers (odd index Euler): T_1=1, T_3=2, T_5=16, T_7=272, T_9=7936
    # Secant numbers (even index): S_0=1, S_2=1, S_4=5, S_6=61, S_8=1385

    # Zigzag numbers (up/down alternating perms): 1, 1, 1, 2, 5, 16, 61, 272, 1385, 7936
    # A_n for n=0,1,2,...
    zigzag = [1, 1, 1, 2, 5, 16, 61, 272, 1385, 7936, 50521]

    print("\n--- Euler zigzag numbers (A_n) ---")
    print("  A_0=1, A_1=1, A_2=1, A_3=2, A_4=5, A_5=16, A_6=61, A_7=272, A_8=1385, A_9=7936")

    # w_0 * something = zigzag?
    # n=3: w_0=-1/2. A_2=1. w_0 = -A_2/2?
    # n=5: w_0=1. A_4=5. Nope.
    # n=7: w_0=-17/4. A_6=61. Nope.
    # n=9: w_0=31. A_8=1385. Nope.

    # Try w_0 * 2^{(n-1)/2}:
    print("\nw_0 * 2^{(n-1)/2}:")
    for n in [3, 5, 7, 9]:
        val = w0_vals[n] * 2**((n-1)//2)
        print(f"  n={n}: {val} = {float(val):.2f}")
    # n=3: -1/2 * 2 = -1. n=5: 1*4=4. n=7: -17/4*8=-34. n=9: 31*16=496.

    # Try Eulerian numbers A(n,k)
    # A(n,k) = number of permutations of {1,...,n} with exactly k ascents
    def eulerian(nn, kk):
        """Eulerian number <n, k>."""
        return sum((-1)**j * comb(nn + 1, j) * (kk + 1 - j)**nn for j in range(kk + 2))

    print("\n--- Eulerian numbers <n,k> ---")
    for nn in range(1, 10):
        vals = [eulerian(nn, k) for k in range(nn)]
        print(f"  n={nn}: {vals}")

    # Check: alternating sum of Eulerian numbers
    # sum_k (-1)^k A(n,k) x^k evaluated somewhere?

    # Look at the constant terms more carefully
    # The Fourier decomposition has W(r) = sum_k w_k * r^k
    # w_k(T) = const_k + (coeff of t3 in w_k)*t3 + ...
    # So const_k = w_k when t3=t5=...=0, i.e., the transitive tournament (no odd cycles)

    print("\n--- KEY INSIGHT: constant = W(r) for transitive tournament ---")
    print("The transitive tournament T_trans has t_k = 0 for all odd k >= 3.")
    print("So const terms = Fourier coefficients of W_trans(r).")
    print("Let's verify by computing W(r) for T_trans directly.")

    for n in [3, 5, 7, 9]:
        # Transitive tournament: A[i][j] = 1 iff i < j (vertices ranked 0 < 1 < ... < n-1)
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                A[i][j] = 1
        poly = compute_W_poly(A, n)
        print(f"\n  n={n}: W_trans(r) coefficients:")
        for deg in sorted(poly.keys(), reverse=True):
            if poly[deg] != 0:
                print(f"    r^{deg}: {poly[deg]} = {float(poly[deg]):.6f}")

        # The transitive tournament has H = 1 (exactly one Hamiltonian path)
        H = count_H(A, n)
        print(f"    H(T_trans) = {H}")
        W_at_half = sum(poly.get(d, Fraction(0)) / Fraction(2)**d for d in range(n))
        print(f"    W(1/2) should = H? No, W(r) sums with r edge weights.")
        # Actually H = sum w_k / 2^k
        H_check = sum(poly.get(d, Fraction(0)) / Fraction(2)**d for d in range(n))
        print(f"    sum w_k/2^k = {H_check} = {float(H_check):.4f}")

    # Now let's look at whether the transitive tournament polynomial
    # has a known closed form
    print("\n--- W_trans(r) = number of permutations weighted by (r +/- 1/2) ---")
    print("For T_trans (i beats j iff i<j), W_trans(r) counts perms pi")
    print("where each consecutive pair (pi_k, pi_{k+1}) contributes")
    print("  r + 1/2 if pi_k < pi_{k+1} (ascent)")
    print("  r - 1/2 if pi_k > pi_{k+1} (descent)")
    print()
    print("So W_trans(r) = sum_{pi in S_n} (r+1/2)^{asc(pi)} * (r-1/2)^{des(pi)}")
    print("where asc + des = n-1.")
    print()
    print("Using Eulerian polynomials: A_n(t) = sum_k A(n,k) t^k")
    print("W_trans(r) = sum_k A(n,k) (r+1/2)^k (r-1/2)^{n-1-k}")

    # Verify this formula
    for n in [3, 5, 7]:
        A_trans = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                A_trans[i][j] = 1
        poly_bf = compute_W_poly(A_trans, n)

        # Compute via Eulerian formula
        half = Fraction(1, 2)
        poly_eul = {}
        for k in range(n):
            ek = eulerian(n, k)
            # (r+1/2)^k * (r-1/2)^{n-1-k}
            # Expand using binomial theorem
            for i in range(k + 1):
                for j in range(n - k):
                    deg = i + j
                    coeff = ek * comb(k, i) * comb(n - 1 - k, j) * half**(k - i) * ((-half)**(n - 1 - k - j))
                    poly_eul[deg] = poly_eul.get(deg, Fraction(0)) + Fraction(coeff).limit_denominator(10**15)

        match = True
        for deg in range(n):
            bf = poly_bf.get(deg, Fraction(0))
            eu = poly_eul.get(deg, Fraction(0))
            if bf != eu:
                match = False
                print(f"  n={n} deg={deg}: MISMATCH bf={bf} eul={eu}")
        print(f"  n={n}: Eulerian formula {'MATCHES' if match else 'FAILS'}")

    # Now analyze constant terms via the Eulerian formula
    print("\n--- Constant term via Eulerian numbers ---")
    print("w_0 = W_trans(0) = sum_k A(n,k) * (1/2)^k * (-1/2)^{n-1-k}")
    print("     = (1/2)^{n-1} * sum_k A(n,k) * (-1)^{n-1-k}")
    print("     = (1/2)^{n-1} * (-1)^{n-1} * sum_k A(n,k) * (-1)^k  [since (-1)^{-k}=(-1)^k]")
    print("     Wait, let me be more careful...")

    for n in [3, 5, 7, 9]:
        half = Fraction(1, 2)
        # w_0 = sum_k A(n,k) * (1/2)^k * (-1/2)^{n-1-k}
        w0 = Fraction(0)
        for k in range(n):
            ek = eulerian(n, k)
            term = Fraction(ek) * half**k * (-half)**(n - 1 - k)
            w0 += term
        print(f"  n={n}: w_0 via Eulerian = {w0} = {float(w0):.6f}")

    # For odd n, n-1 is even, so (-1/2)^{n-1-k} = (-1)^{n-1-k} * (1/2)^{n-1-k}
    # w_0 = (1/2)^{n-1} * sum_k A(n,k) * (-1)^{n-1-k}
    print("\n  Simplified: w_0 = (1/2)^{n-1} * sum_k A(n,k)*(-1)^{n-1-k}")
    print("  = (1/2)^{n-1} * A_n(-1) * (-1)^{n-1}  [where A_n(t) = sum A(n,k) t^k]")
    print()
    print("  But A_n(-1) for the Eulerian polynomial:")
    for n in [3, 5, 7, 9]:
        An_neg1 = sum(eulerian(n, k) * (-1)**k for k in range(n))
        print(f"    A_{n}(-1) = {An_neg1}")

    # These should be related to 2^n * B_n (Bernoulli) or similar
    # Actually the alternating Eulerian sum is:
    # sum_k A(n,k) (-1)^k = (-1)^{n-1} * n! * B_n  ... no
    # Let me just look at the numbers

    # A_3(-1) = A(3,0) - A(3,1) + A(3,2) = 1 - 4 + 1 = -2
    # w_0 for n=3 = (1/4) * (-1)^2 * (-2) = -1/2. Check!

    # A_5(-1) = 1 - 26 + 66 - 26 + 1 = 16
    # w_0 for n=5 = (1/16) * (-1)^4 * 16 = 1. Check!

    # Look up: the alternating Eulerian sum relates to TANGENT/SECANT numbers
    # Actually sum_k (-1)^k A(n,k) = E_n (Euler number) * something?

    # Let's check connection to tangent numbers
    # T_n = |E_n| for odd n: T_1=1, T_3=2, T_5=16, T_7=272, T_9=7936

    print("\n--- Connection to tangent numbers ---")
    for n in [3, 5, 7, 9]:
        An_neg1 = sum(eulerian(n, k) * (-1)**k for k in range(n))
        half_pow = Fraction(1, 2)**(n - 1)
        w0 = half_pow * (-1)**(n - 1) * An_neg1
        print(f"  n={n}: A_n(-1) = {An_neg1}, w_0 = {w0}")
        # Tangent number T_n
        T_n = zigzag[n]
        print(f"    Tangent/zigzag number A_{n} = {T_n}")
        print(f"    A_n(-1) / n! = {Fraction(An_neg1, factorial(n))}")
        # Check: is |A_n(-1)| = 2^n * |B_n| * n! / n  or something?

    # Actually, the correct identity is:
    # sum_k A(n,k) t^k = sum_{j>=0} (j+1)^n t^j * (1-t)^{n+1}  ... no
    # Let's try a different angle.
    # We know w_0 = W_trans(0).
    # For T_trans, every permutation pi contributes product of (0 + A[pi_k][pi_{k+1}] - 1/2)
    # = product of (+1/2 or -1/2) depending on ascent/descent
    # = (1/2)^{n-1} * (-1)^{des(pi)}
    # So w_0 = (1/2)^{n-1} * sum_pi (-1)^{des(pi)}

    # And sum_pi (-1)^{des(pi)} = sum_k A(n,k) (-1)^k = A_n(-1)

    # For the general Fourier coefficient at r^{2j}:
    # w_{2j} in W_trans = coefficient of r^{2j} in
    # sum_k A(n,k) (r+1/2)^k (r-1/2)^{n-1-k}

    # Let's compute w_{n-3} constant = coeff of r^{n-3} in W_trans
    print("\n--- General constant terms ---")
    for n in [3, 5, 7, 9]:
        A_trans = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                A_trans[i][j] = 1
        poly = compute_W_poly(A_trans, n)
        print(f"\n  n={n}: all even-degree coefficients of W_trans:")
        for deg in range(n - 1, -1, -2):
            val = poly.get(deg, Fraction(0))
            # Also compute (n-1)! * ratio
            if val != 0:
                ratio = val / factorial(n)
                print(f"    w_{deg} = {val} = {float(val):.6f}  (ratio to n!: {ratio})")

    # Let's check if the constant terms satisfy a recurrence
    # w_0 values: n=3: -1/2, n=5: 1, n=7: -17/4, n=9: 31
    # Multiply by 2^{(n-1)/2}: -2, 4, -34, 496
    # Divide by (-1)^{(n-1)/2}: 2, 4, 34, 496  [signs: (-1)^1, (-1)^2, (-1)^3, (-1)^4]
    # Hmm: 2, 4, 34, 496. Not obviously standard.

    # Let's try: w_0 * 4^{(n-1)/2} = w_0 * 2^{n-1}
    print("\n--- w_0 * 2^{n-1} ---")
    for n in [3, 5, 7, 9]:
        val = w0_vals[n] * 2**(n - 1)
        print(f"  n={n}: {val}")
    # -2, 16, -272, 7936
    # These are (-1)^{(n-1)/2} * tangent numbers!
    # T_1=1... no. Let's check:
    # n=3: -2. T_3 = 2. (-1)^1 * 2 = -2. YES!
    # n=5: 16. T_5 = 16. (-1)^2 * 16 = 16. YES!
    # n=7: -272. T_7 = 272. (-1)^3 * 272 = -272. YES!
    # n=9: 7936. T_9 = 7936. (-1)^4 * 7936 = 7936. YES!

    print("\n  *** MATCH: w_0 * 2^{n-1} = (-1)^{(n-1)/2} * T_n ***")
    print("  where T_n is the n-th tangent number (zigzag number A_n for odd n)")
    print("  Equivalently: w_0 = (-1)^{(n-1)/2} * T_n / 2^{n-1}")
    print()

    for n in [3, 5, 7, 9]:
        T_n = zigzag[n]
        expected_w0 = Fraction((-1)**((n - 1) // 2) * T_n, 2**(n - 1))
        actual_w0 = w0_vals[n]
        print(f"  n={n}: T_{n}={T_n}, expected w_0={expected_w0}, actual={actual_w0}  "
              f"{'OK' if expected_w0 == actual_w0 else 'FAIL'}")

    # Now let's check ALL constant terms, not just w_0
    print("\n--- All constant terms vs tangent/secant numbers ---")
    # For the transitive tournament, W_trans(r) = sum_k A(n,k) (r+1/2)^k (r-1/2)^{n-1-k}
    # This is a polynomial in r. The constant terms at each degree
    # should relate to derivatives of this at r=0.

    # Actually, let's express W_trans(r) differently.
    # Let s = r + 1/2, d = r - 1/2, so s - d = 1, s + d = 2r.
    # W_trans(r) = sum_k A(n,k) s^k d^{n-1-k} = A_n(s/d) * d^{n-1}  ... messy.

    # Better: substitute r = i*x/2 (imaginary) and use tangent/secant generating function?
    # Actually let's just compute the pattern numerically.

    # We know w_0 = (-1)^{(n-1)/2} * T_n / 2^{n-1}
    # Let's check higher constant terms.

    # For n=9, the constant terms are:
    # w_8 = 362880 = 9!
    # w_6 = -211680
    # w_4 = 40320
    # w_2 = -2640
    # w_0 = 31

    # Ratios to (degree+1)!:
    print("\nConstant term of w_k divided by (k+1)!:")
    for n in [3, 5, 7, 9]:
        A_trans = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                A_trans[i][j] = 1
        poly = compute_W_poly(A_trans, n)
        for deg in range(n - 1, -1, -2):
            val = poly.get(deg, Fraction(0))
            if deg >= 0:
                ratio = val / factorial(deg + 1) if deg >= 0 else val
                print(f"  n={n}, w_{deg} / (deg+1)! = {val}/{factorial(deg+1)} = {ratio} = {float(ratio):.8f}")

    # Check the w_{n-3} formula constant = -(n-2)! * C(n,3) / 2
    print("\n--- w_{n-3} constant = -C(n,3) * (n-2)! / 2: verified above ---")


# ─────────────────────────────────────────────────────────────────────
# TASK 3: t3 coefficient patterns
# ─────────────────────────────────────────────────────────────────────

def task3_t3_coefficients():
    print("\n" + "=" * 70)
    print("TASK 3: t3 coefficient patterns across Fourier degrees")
    print("=" * 70)

    # Known t3 coefficients in each w_k:
    # n=3: w_0 has t3 coeff = 2
    # n=5: w_2 has t3 coeff = 12, w_0 has t3 coeff = -1
    # n=7: w_4 has t3 coeff = 240, w_2 has t3 coeff = -60, w_0 has t3 coeff = 2
    # n=9: w_6 has t3 coeff = 10080, w_4 has t3 coeff = -4200, w_2 has t3 coeff = 462, w_0 has t3 coeff = -17/2

    t3_coeffs = {
        3: {0: Fraction(2)},
        5: {2: Fraction(12), 0: Fraction(-1)},
        7: {4: Fraction(240), 2: Fraction(-60), 0: Fraction(2)},
        9: {6: Fraction(10080), 4: Fraction(-4200), 2: Fraction(462), 0: Fraction(-17, 2)},
    }

    # TOP t3 coefficient (in w_{n-3}):
    print("\nTOP t3 coefficient (in w_{n-3}):")
    print("  Formula: 2*(n-2)!")
    for n in [3, 5, 7, 9]:
        top_deg = n - 3
        actual = t3_coeffs[n][top_deg]
        expected = 2 * factorial(n - 2)
        print(f"  n={n}: coeff={actual}, 2*(n-2)!={expected}  {'OK' if actual == expected else 'FAIL'}")

    # SECOND t3 coefficient (in w_{n-5}):
    print("\nSECOND t3 coefficient (in w_{n-5}):")
    for n in [5, 7, 9]:
        deg = n - 5
        val = t3_coeffs[n][deg]
        print(f"  n={n}: coeff of t3 in w_{deg} = {val}")
        # Ratio to (n-4)!
        ratio = val / factorial(n - 4)
        print(f"    ratio to (n-4)! = {ratio} = {float(ratio):.6f}")
    # n=5: -1/1 = -1. n=7: -60/6 = -10. n=9: -4200/120 = -35.
    # Pattern: -1, -10, -35. These are -C(n,3) * something?
    # -C(3,2)=-3? No. -1 = -C(2,2). -10 = -C(5,2). -35 = -C(7,2).
    # Wait: for n=5 the ratio is -1, n=7: -10, n=9: -35
    # -C(2,2)=-1, -C(5,2)=-10, -C(7,2)=-21... no -35 != -21
    # -C(n-2, 2)? n=5: -C(3,2)=-3. No.
    # Try: -1 = -1*1. -10 = -2*5. -35 = -5*7. Hmm.
    # -C(n,3)/4 * something? -C(5,3)/4 = -10/4 = -5/2. No.
    # Actually: -1, -10, -35.  Differences: -9, -25. Nah.
    # -C(2*1, 1)/2 = -1. -C(4,2) = -6. No.
    # Let me just check: are these -C(n-2, 2) * ... ?
    # n=5: val=-1, n=7: val=-10, n=9: val=-35
    # Ratios: -35/-10 = 3.5, -10/-1=10.
    # Or: val / (-C(n,3)/4) = n=5: -1/(-10/4) = 4/10 = 2/5
    #                          n=7: -10/(-35/4) = 40/35 = 8/7
    #                          n=9: -35/(-84/4) = -35/(-21) = 5/3

    # Let me try something else. Look at the full coefficient matrix.
    print("\n--- Full t3 coefficient table ---")
    print(f"{'n':>3} | {'w_{n-3}':>10} | {'w_{n-5}':>10} | {'w_{n-7}':>10} | {'w_{n-9}':>10}")
    print("-" * 55)
    for n in [3, 5, 7, 9]:
        row = f"{n:3d} |"
        for j in range(4):
            deg = n - 3 - 2 * j
            if deg >= 0 and deg in t3_coeffs[n]:
                row += f" {str(t3_coeffs[n][deg]):>10s} |"
            elif deg >= 0:
                row += f" {'0':>10s} |"
            else:
                row += f" {'':>10s} |"
        print(row)

    # Now look at ratios: coeff_t3 in w_{n-3-2j} / coeff_t3 in w_{n-3}
    print("\n--- Ratios: coeff_t3(w_{n-3-2j}) / coeff_t3(w_{n-3}) ---")
    for n in [3, 5, 7, 9]:
        top = t3_coeffs[n][n - 3]
        for j in range(4):
            deg = n - 3 - 2 * j
            if deg >= 0 and deg in t3_coeffs[n]:
                ratio = t3_coeffs[n][deg] / top
                print(f"  n={n}, j={j}: ratio = {ratio} = {float(ratio):.8f}")

    # j=0: always 1 (trivially)
    # j=1: n=5: -1/12. n=7: -60/240=-1/4. n=9: -4200/10080=-5/12
    # j=2: n=7: 2/240=1/120. n=9: 462/10080=11/240
    # j=3: n=9: (-17/2)/10080 = -17/20160

    # Ratios at j=1: -1/12, -1/4, -5/12. Pattern: -(j*(2j-1))/12? j=1: -1/12. Hmm.
    # Actually: -1/12 = -1/12, -1/4 = -3/12, -5/12 = -5/12
    # Numerators: -1, -3, -5. Arithmetic! So ratio at j=1 = -(2*(n-3)/2-1+1)/(12)...
    # n=5: -(2*1-1)/12 = -1/12. n=7: -(2*2-1)/12 = -3/12 = -1/4. n=9: -(2*3-1)/12 = -5/12.
    # So ratio = -(n-4)/something. n=5: n-4=1. n=7: n-4=3. n=9: n-4=5. Yes!
    # ratio at j=1 = -(n-4)/12

    print("\n  At j=1: ratio = -(n-4)/12 ?")
    for n in [5, 7, 9]:
        expected = -Fraction(n - 4, 12)
        actual = t3_coeffs[n][n - 5] / t3_coeffs[n][n - 3]
        print(f"    n={n}: expected={expected}, actual={actual}  {'OK' if expected == actual else 'FAIL'}")

    # ratio at j=1 = -(n-4)/12 means:
    # coeff_t3(w_{n-5}) = 2*(n-2)! * (-(n-4)/12) = -(n-4)*(n-2)!/6
    # n=5: -(1)(6)/6 = -1. CHECK.
    # n=7: -(3)(120)/6 = -60. CHECK.
    # n=9: -(5)(5040)/6 = -4200. CHECK.
    # So coeff_t3(w_{n-5}) = -(n-4)*(n-2)!/6 = -(n-4)!/3 ... no
    # = -(n-2)!*(n-4)/6

    print("\n  FORMULA: coeff_t3 in w_{n-5} = -(n-4)*(n-2)!/6")

    # Can also write: -(n-4)*(n-2)!/6 = -(n-2)!*(n-4)/6 = -C(n-2,2)*(n-4)!/1 ... hmm
    # More naturally: = -C(n-4,1) * (n-2)! / 6

    # What about j=2? coeff_t3(w_{n-7}):
    # n=7: 2, n=9: 462
    # ratio to 2*(n-2)!: n=7: 2/240=1/120. n=9: 462/10080=11/240
    # Guess: involves C(n-4,2)/something?
    # n=7: C(3,2)=3. 1/120 = 1/120. n=9: C(5,2)=10. 11/240.
    # Not obvious. Let me compute directly.

    print("\n  coeff_t3 in w_{n-7}:")
    for n in [7, 9]:
        deg = n - 7
        val = t3_coeffs[n][deg]
        print(f"    n={n}: {val}")
        # Factor out (n-4)!
        if n == 7:
            # (n-4)! = 6. val = 2. ratio = 1/3
            print(f"      / (n-6)! = {val / factorial(n-6)}")
        if n == 9:
            # (n-6)! = 6. val = 462. ratio = 77
            print(f"      / (n-6)! = {val / factorial(n-6)}")

    # j=3: n=9: coeff = -17/2
    # This equals w_0's t3 coefficient.

    # Let me also verify by computing W(r) for several tournaments with known t3
    # and extracting the t3 coefficient numerically
    print("\n--- Numerical verification of t3 coefficients at n=7 ---")
    n = 7
    for seed in range(5):
        A = random_tournament(n, seed)
        poly = compute_W_poly(A, n)
        t3 = count_directed_cycles(A, n, 3)

        # From formula: w_4 = -2100 + 240*t3
        w4_formula = -2100 + 240 * t3
        w4_actual = poly.get(4, Fraction(0))
        print(f"  seed={seed}: t3={t3}, w_4_formula={w4_formula}, w_4_actual={w4_actual}  "
              f"{'OK' if Fraction(w4_formula) == w4_actual else 'FAIL'}")


# ─────────────────────────────────────────────────────────────────────
# TASK 4: Verify w_{n-3} = 2*(n-2)! * (t3 - C(n,3)/4)
# ─────────────────────────────────────────────────────────────────────

def task4_wn3_formula():
    print("\n" + "=" * 70)
    print("TASK 4: Verify w_{n-3} = 2*(n-2)! * (t3 - C(n,3)/4)")
    print("=" * 70)

    # constant = 2*(n-2)! * (-C(n,3)/4) = -(n-2)!*C(n,3)/2
    print("\nConstant term = -(n-2)! * C(n,3) / 2:")
    for n in [3, 5, 7, 9]:
        expected = -Fraction(factorial(n - 2) * comb(n, 3), 2)
        # From known data
        known = {3: Fraction(-1, 2), 5: Fraction(-30), 7: Fraction(-2100), 9: Fraction(-211680)}
        actual = known[n]
        print(f"  n={n}: -(n-2)!*C(n,3)/2 = -{factorial(n-2)}*{comb(n,3)}/2 = {expected}, actual={actual}  "
              f"{'OK' if expected == actual else 'FAIL'}")

    # Full formula verification on random tournaments
    print("\nFull w_{n-3} verification on random tournaments:")
    for n in [3, 5, 7, 9]:
        print(f"\n  n={n}:")
        N_samples = 10 if n <= 7 else 5
        all_ok = True
        for seed in range(N_samples):
            A = random_tournament(n, seed)
            poly = compute_W_poly(A, n)
            t3 = count_directed_cycles(A, n, 3)
            expected = 2 * factorial(n - 2) * (Fraction(t3) - Fraction(comb(n, 3), 4))
            actual = poly.get(n - 3, Fraction(0))
            ok = expected == actual
            if not ok:
                all_ok = False
            if seed < 3 or not ok:
                print(f"    seed={seed}: t3={t3}, formula={expected}, actual={actual}  {'OK' if ok else 'FAIL'}")
        if all_ok:
            print(f"    All {N_samples} samples match.")


# ─────────────────────────────────────────────────────────────────────
# TASK 5: Alternating structure and Euler zigzag connection
# ─────────────────────────────────────────────────────────────────────

def task5_alternating_structure():
    print("\n" + "=" * 70)
    print("TASK 5: Alternating sign structure and Euler zigzag numbers")
    print("=" * 70)

    zigzag = [1, 1, 1, 2, 5, 16, 61, 272, 1385, 7936, 50521]
    # Tangent numbers T_n (odd n): T_1=1, T_3=2, T_5=16, T_7=272, T_9=7936
    # Secant numbers S_n (even n): S_0=1, S_2=1, S_4=5, S_6=61, S_8=1385

    w0_vals = {3: Fraction(-1, 2), 5: Fraction(1), 7: Fraction(-17, 4), 9: Fraction(31)}

    print("\n--- w_0 and tangent numbers ---")
    print("THEOREM: w_0 = (-1)^{(n-1)/2} * T_n / 2^{n-1}")
    print("where T_n is the n-th tangent number (for odd n).")
    print()
    for n in [3, 5, 7, 9]:
        T_n = zigzag[n]
        sign = (-1)**((n - 1) // 2)
        expected = Fraction(sign * T_n, 2**(n - 1))
        actual = w0_vals[n]
        print(f"  n={n}: T_{n}={T_n}, sign=({'-' if sign < 0 else '+'}1), "
              f"formula={expected}, actual={actual}  {'OK' if expected == actual else 'FAIL'}")

    # Now check if ALL constant terms relate to zigzag/tangent/secant numbers
    # The constant terms are the coefficients of W_trans(r)
    # W_trans(r) = sum_k A(n,k) (r+1/2)^k (r-1/2)^{n-1-k}
    # = sum_k A(n,k) * sum_{i,j} C(k,i) C(n-1-k,j) r^{i+j} (1/2)^{k-i} (-1/2)^{n-1-k-j}

    # The coefficient of r^m in W_trans is:
    # sum_k A(n,k) sum_{i+j=m} C(k,i) C(n-1-k,j) (1/2)^{k-i} (-1/2)^{n-1-k-j}
    # = (1/2)^{n-1-m} sum_k A(n,k) sum_{i+j=m} C(k,i) C(n-1-k,j) (-1)^{n-1-k-j}

    # For the transitive tournament we also know:
    # W_trans(r) = sum_{pi in S_n} product_{edges} (r + s_edge - 1/2)
    # where s_edge = 1 for ascent, 0 for descent
    # = sum_{pi} (r+1/2)^{asc(pi)} (r-1/2)^{des(pi)}

    # Let's verify the alternating sign pattern
    print("\n--- Sign pattern of constant terms ---")
    for n in [3, 5, 7, 9]:
        A_trans = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                A_trans[i][j] = 1
        poly = compute_W_poly(A_trans, n)
        print(f"\n  n={n}:")
        for deg in range(n - 1, -1, -2):
            val = poly.get(deg, Fraction(0))
            sign_char = '+' if val > 0 else '-'
            print(f"    w_{deg} = {sign_char}{abs(val)}")

    # Signs at n=9: w_8=+, w_6=-, w_4=+, w_2=-, w_0=+
    # Alternating starting from +! Pattern: (-1)^j for w_{n-1-2j}

    print("\n--- Sign pattern: w_{n-1-2j} has sign (-1)^j ---")
    for n in [3, 5, 7, 9]:
        A_trans = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                A_trans[i][j] = 1
        poly = compute_W_poly(A_trans, n)
        ok = True
        for j in range((n + 1) // 2):
            deg = n - 1 - 2 * j
            val = poly.get(deg, Fraction(0))
            expected_sign = (-1)**j
            actual_sign = 1 if val > 0 else -1
            if expected_sign != actual_sign:
                ok = False
        print(f"  n={n}: alternating {'OK' if ok else 'FAIL'}")

    # Connection between ALL constant terms and tangent numbers via generating function
    print("\n--- Generating function connection ---")
    print("W_trans(r) = sum_pi (r+1/2)^{asc} (r-1/2)^{des}")
    print("Let u = r+1/2, v = r-1/2, so u-v=1.")
    print("W_trans = sum_{k=0}^{n-1} A(n,k) u^k v^{n-1-k}")
    print("This is the Eulerian polynomial in u/v times v^{n-1}.")
    print()
    print("The exponential generating function of A_n(t) is:")
    print("  sum_{n>=1} A_n(t) x^n/n! = (t-1)/(t - exp((t-1)x))")
    print()
    print("Setting t = u/v and multiplying by v^{n-1}/n!:")
    print("  W_trans(r)/n! = coeff of x^n in v^{-1}(u/v - 1)/(u/v - exp((u/v-1)x))")
    print("                = coeff of x^n in (u-v)/(u - v*exp((u-v)x/v))")
    print("Since u-v=1: = coeff of x^n in 1/(u - v*exp(x/v))... ")
    print()
    print("For r=0: u=1/2, v=-1/2.")
    print("  W_trans(0)/n! = coeff of x^n in 1/(1/2 + 1/2*exp(-2x))")
    print("                = coeff of x^n in 2/(1 + exp(-2x))")
    print("                = coeff of x^n in 1/cosh(x) * something?")
    print()

    # 2/(1+e^{-2x}) = 2e^x/(e^x + e^{-x}) = e^x/cosh(x) = (sinh(x)+cosh(x))/cosh(x)
    #                = 1 + tanh(x)

    print("  2/(1+exp(-2x)) = 1 + tanh(x)")
    print()
    print("  So W_trans(0) = n! * [x^n] (1 + tanh(x))")
    print("  = n! * [x^n] tanh(x)  for n >= 1")
    print()
    print("  The Taylor series of tanh(x) = sum_{n odd} T_n * x^n / n!")
    print("  where T_n are tangent numbers (with appropriate signs).")
    print("  tanh(x) = x - x^3/3 + 2x^5/15 - 17x^7/315 + 62x^9/2835 - ...")
    print()

    # Let's verify: tanh(x) coefficients times n!
    # tanh(x) = sum B_{2n} * 4^n * (4^n - 1) / (2n)! * x^{2n-1}  ... complicated
    # Actually tanh(x) = sum_{n>=1, n odd} T_n / n! * x^n with
    # T_1 = 1, T_3 = -1/3, ...no

    # Let's just directly check n!*[x^n]tanh(x) vs w_0:
    # tanh(x) = x - x^3/3 + 2x^5/15 - 17x^7/315 + 62x^9/2835 - ...
    # Coefficients: a_1=1, a_3=-1/3, a_5=2/15, a_7=-17/315, a_9=62/2835
    # Multiply by n!:
    # n=1: 1*1=1. n=3: 6*(-1/3)=-2. n=5: 120*(2/15)=16. n=7: 5040*(-17/315)=-272.
    # n=9: 362880*(62/2835)=7936.

    # But our w_0 values are: n=3: -1/2, n=5: 1, n=7: -17/4, n=9: 31
    # And w_0 * 2^{n-1}: n=3: -2, n=5: 16, n=7: -272, n=9: 7936

    # So w_0 * 2^{n-1} = n! * [x^n] tanh(x)
    # => w_0 = n! * [x^n] tanh(x) / 2^{n-1}

    print("  VERIFICATION: n! * [x^n] tanh(x) vs w_0 * 2^{n-1}:")
    tanh_coeffs = {1: Fraction(1), 3: Fraction(-1, 3), 5: Fraction(2, 15),
                   7: Fraction(-17, 315), 9: Fraction(62, 2835)}
    for n in [3, 5, 7, 9]:
        nfact_times_coeff = factorial(n) * tanh_coeffs[n]
        w0_times_2 = w0_vals[n] * 2**(n - 1)
        print(f"    n={n}: n!*[x^n]tanh = {nfact_times_coeff}, w_0*2^{n-1} = {w0_times_2}  "
              f"{'OK' if nfact_times_coeff == w0_times_2 else 'FAIL'}")

    print()
    print("CONFIRMED: w_0 = n! * [x^n] tanh(x) / 2^{n-1}")
    print()
    print("More generally, the FULL transitive polynomial satisfies:")
    print("  W_trans(r) / n! = [x^n] (1 + tanh(x)) with x -> (r+1/2)x, etc.")
    print("  The constant terms are controlled by the Taylor series of tanh(x).")

    # Now let's investigate whether ALL constant terms (not just w_0) have
    # a clean formula via the tanh generating function
    print("\n--- ALL constant terms via tanh generating function ---")
    print("Since W_trans(r)/n! = [x^n] sum_k A(n,k) u^k v^{n-1-k} / n!")
    print("and this equals [x^n](1+tanh((u-v)x/2)) when u-v=1...")
    print()
    print("Actually, let's derive the general formula differently.")
    print("For the transitive tournament, W(r) = sum_pi prod_edges (r + epsilon)")
    print("where epsilon = +1/2 for ascent, -1/2 for descent.")
    print()
    print("This means coeff of r^m in W_trans = sum_pi C(n-1, m) * ...")
    print("Actually, let's just tabulate all constant terms / n! and look for patterns.")
    print()

    for n in [3, 5, 7, 9]:
        A_trans = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                A_trans[i][j] = 1
        poly = compute_W_poly(A_trans, n)
        print(f"  n={n}:")
        for deg in range(n - 1, -1, -2):
            val = poly.get(deg, Fraction(0))
            ratio_nfact = val / factorial(n)
            # Also check val / ((deg+1)! or (n-1 choose deg)... )
            ratio_binom = val / comb(n - 1, deg) if comb(n-1, deg) > 0 else None
            print(f"    w_{deg}/n! = {ratio_nfact} = {float(ratio_nfact):.10f}")

    # The ratios w_{n-1-2j}/n! at fixed j across n:
    print("\n--- Ratio w_{n-1-2j}/n! at fixed j across n ---")
    all_polys = {}
    for n in [3, 5, 7, 9]:
        A_trans = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                A_trans[i][j] = 1
        all_polys[n] = compute_W_poly(A_trans, n)

    for j in range(5):
        print(f"\n  j={j} (w_{{n-1-{2*j}}}):")
        for n in [3, 5, 7, 9]:
            deg = n - 1 - 2 * j
            if deg >= 0:
                val = all_polys[n].get(deg, Fraction(0))
                ratio = val / factorial(n)
                print(f"    n={n}: w_{deg}/n! = {ratio} = {float(ratio):.10f}")

    # j=0: always 1 (w_{n-1} = n!)
    # j=1: n=3: -1/12, n=5: -1/4, n=7: -5/12, n=9: -7/12
    #   = -(n-2)/12? n=3: -1/12. n=5: -3/12=-1/4. n=7: -5/12. n=9: -7/12. YES!
    print("\n  j=1: w_{n-3}/n! = -(n-2)/12 ?")
    for n in [3, 5, 7, 9]:
        expected = -Fraction(n - 2, 12)
        actual = all_polys[n].get(n - 3, Fraction(0)) / factorial(n)
        print(f"    n={n}: expected={expected}, actual={actual}  {'OK' if expected == actual else 'FAIL'}")
    # But wait: w_{n-3}/n! = -(n-2)!*C(n,3)/(2*n!) = -C(n,3)/(2*n*(n-1))
    # = -(n(n-1)(n-2)/6) / (2n(n-1)) = -(n-2)/12. Confirmed!

    print("\n  FORMULA: w_{n-3} constant / n! = -(n-2)/12 = -C(n,3)/(2*n*(n-1))")

    # j=2: n=5: 1/120, n=7: 231/(5040) = 11/240, n=9: 40320/362880 = 1/9
    print("\n  j=2: w_{n-5}/n!:")
    for n in [5, 7, 9]:
        val = all_polys[n].get(n - 5, Fraction(0)) / factorial(n)
        print(f"    n={n}: {val} = {float(val):.10f}")
    # 1/120, 11/240, 1/9
    # = 1/120, 11/240, 1/9
    # Common denominator form: 1/120 = 2/240 = 0.00833..
    # Try (n-2)(n-4) / (12*something)?
    # n=5: 3*1/x = 1/120 => x=360. n=7: 5*3/x = 11/240 => x = 15*240/11... no.
    # Try: (n-2)^2*(n-4) / (1440)?
    # n=5: 9*1/1440 = 1/160. No.
    # Just leave this open.

    # Summary of w_0 connection
    print("\n" + "=" * 70)
    print("SUMMARY OF FINDINGS")
    print("=" * 70)
    print()
    print("1. CONSTANT TERMS = W_trans(r) coefficients (transitive tournament)")
    print()
    print("2. w_0 = (-1)^{(n-1)/2} * T_n / 2^{n-1}")
    print("   where T_n is the n-th tangent number: T_3=2, T_5=16, T_7=272, T_9=7936")
    print("   Equivalently: w_0 * 2^{n-1} = n! * [x^n] tanh(x)")
    print()
    print("3. w_{n-1} = n!  (trivial)")
    print()
    print("4. w_{n-3} constant = -(n-2)! * C(n,3) / 2 = -n! * (n-2) / 12")
    print("   Confirmed for n = 3, 5, 7, 9.")
    print()
    print("5. Signs alternate: w_{n-1-2j} constant has sign (-1)^j")
    print("   This follows from tanh(x) having alternating signs.")
    print()
    print("6. t3 TOP coefficient: 2*(n-2)! in w_{n-3}  (confirmed)")
    print("   t3 SECOND coefficient: -(n-4)*(n-2)!/6 in w_{n-5}  (confirmed n=5,7,9)")
    print()
    print("7. The full constant-term generating function is")
    print("   sum_{n odd} W_trans(r) x^n / n! = f(r, x)")
    print("   where f(0, x) = tanh(x) (for n >= 1 terms).")
    print("   The transitive tournament polynomial is a 'tangent polynomial'.")


# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    ok1 = task1_verify_n9(N=20)
    task2_constant_terms()
    task3_t3_coefficients()
    task4_wn3_formula()
    task5_alternating_structure()
