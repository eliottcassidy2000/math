#!/usr/bin/env python3
"""
DEFORMED EULERIAN PROPERTIES -- analysis of inflated Eulerian coefficients c_k^{(f,d)}.

Tasks:
1. Compute c_k^{(f,d)} for d=4,6,8,10 and all valid f.
2. Check real roots of c_k^{(f,d)} as polynomial in k (or transform variable).
3. Check log-concavity of |c_k^{(f,d)}|.
4. Positivity polytope for n=7 forward-edge distribution.
5. Compare polytope to actual invariant ranges.
6. Minimum of a_k(T) over all n=7 tournaments.

Uses exact rational arithmetic where possible.

opus-2026-03-07
"""
from fractions import Fraction
from math import comb, factorial
from itertools import combinations, permutations
from collections import defaultdict
import random

# ====================================================================
# CORE FUNCTIONS
# ====================================================================

def eulerian_number(n, k):
    """A(n,k) = number of permutations of [n] with exactly k ascents (0-indexed)."""
    if k < 0 or k >= n:
        return 0
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))

def inflated_eulerian(f, d, k):
    """
    c_k^{(f,d)} -- inflated Eulerian coefficient.
    Exact integer arithmetic.
    """
    if k < 0 or k > d:
        return 0
    total = 0
    for j in range(max(0, k - (d - f)), min(f, k) + 1):
        sign = (-1) ** (d - f - k + j)
        total += eulerian_number(f + 1, j) * comb(d - f, k - j) * sign
    return total


# ====================================================================
# TASK 1: Compute c_k^{(f,d)} for d=4,6,8,10
# ====================================================================

def task1():
    print("=" * 78)
    print("TASK 1: Inflated Eulerian coefficients c_k^{(f,d)}")
    print("=" * 78)

    for d in [4, 6, 8, 10]:
        print(f"\n  d = {d}:")
        print(f"  {'f':>3} | " + " ".join(f"{'k='+str(k):>10}" for k in range(d+1)))
        print(f"  " + "-" * (6 + 11*(d+1)))
        for f in range(0, d+1, 2):  # f same parity as d for tournaments (odd cycles)
            coeffs = [inflated_eulerian(f, d, k) for k in range(d+1)]
            print(f"  {f:>3} | " + " ".join(f"{c:>10}" for c in coeffs))
            # Verify palindromy
            pal = all(coeffs[k] == coeffs[d-k] for k in range(d+1))
            # Verify zero-sum for f < d
            zs = sum(coeffs)
            if f < d:
                assert zs == 0, f"Zero-sum failed for f={f}, d={d}: sum={zs}"
            if not pal:
                print(f"  *** PALINDROMY FAILS for f={f}, d={d}")

        # Also show all valid f (not just even)
        print(f"\n  All f values for d={d}:")
        for f in range(d+1):
            coeffs = [inflated_eulerian(f, d, k) for k in range(d+1)]
            s = sum(coeffs)
            pal = all(coeffs[k] == coeffs[d-k] for k in range(d+1))
            print(f"    f={f}: {coeffs}  sum={s}  palindromic={pal}")


# ====================================================================
# TASK 2: Real roots analysis
# ====================================================================

def poly_discriminant_sign(coeffs):
    """
    For degree <= 4, check real roots exactly using discriminant.
    Returns True if all roots are real, False if not, None if degree too high.
    """
    # Strip leading/trailing zeros
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs = coeffs[:-1]
    while len(coeffs) > 1 and coeffs[0] == 0:
        coeffs = coeffs[1:]
    deg = len(coeffs) - 1
    if deg <= 1:
        return True
    if deg == 2:
        a, b, c = Fraction(coeffs[0]), Fraction(coeffs[1]), Fraction(coeffs[2])
        disc = b*b - 4*a*c  # note: reversed order for ax^2+bx+c
        # Actually coeffs[0] + coeffs[1]*x + coeffs[2]*x^2
        # disc = coeffs[1]^2 - 4*coeffs[2]*coeffs[0]
        disc = Fraction(coeffs[1])**2 - 4*Fraction(coeffs[2])*Fraction(coeffs[0])
        return disc >= 0
    return None  # degree too high for simple discriminant

def find_polynomial_roots_numpy(coeffs):
    """Find roots of polynomial with given coefficients [a_0, a_1, ..., a_d]
    interpreted as a_0 + a_1*x + a_2*x^2 + ... + a_d*x^d."""
    try:
        import numpy as np
        # numpy.roots expects [a_d, ..., a_1, a_0] (highest degree first)
        np_coeffs = list(reversed(coeffs))
        # Remove leading zeros
        while len(np_coeffs) > 1 and np_coeffs[0] == 0:
            np_coeffs.pop(0)
        if len(np_coeffs) <= 1:
            return [], True
        roots = np.roots(np_coeffs)
        # Use adaptive tolerance: for palindromic polynomials with roots near 1,
        # numerical errors can be large. Use 1e-6 as threshold.
        all_real = all(abs(r.imag) < 1e-6 for r in roots)
        return sorted(roots, key=lambda r: r.real), all_real
    except ImportError:
        return None, None

def task2():
    print("\n" + "=" * 78)
    print("TASK 2: Real roots of c_k^{(f,d)} as polynomial in x")
    print("        (coefficients of x^0, x^1, ..., x^d)")
    print("=" * 78)

    try:
        import numpy as np
        has_numpy = True
    except ImportError:
        has_numpy = False
        print("  numpy not available -- skipping numerical root finding")
        print("  Will check via discriminant / algebraic methods where possible")

    # First check: Eulerian numbers A(n,k) as polynomial in k
    # Actually, c_k^{(f,d)} is a SEQUENCE indexed by k, not naturally a polynomial.
    # We interpret it as coefficients of a generating polynomial:
    #   P_{f,d}(x) = sum_k c_k^{(f,d)} * x^k
    # and check if this polynomial has all real roots.

    print("\n  Generating polynomial P_{f,d}(x) = sum_k c_k^{(f,d)} x^k")
    print("  (Eulerian polynomials A_n(x) = sum_k A(n,k) x^k have all real roots)")
    print()

    results = {}
    for d in [4, 6, 8, 10]:
        print(f"  d = {d}:")
        for f in range(d+1):
            coeffs = [inflated_eulerian(f, d, k) for k in range(d+1)]
            if has_numpy:
                roots, all_real = find_polynomial_roots_numpy(coeffs)
                if roots is not None:
                    real_roots = [r.real for r in roots if abs(r.imag) < 1e-8]
                    complex_roots = [(r.real, r.imag) for r in roots if abs(r.imag) >= 1e-8]
                    results[(f, d)] = all_real
                    status = "ALL REAL" if all_real else f"COMPLEX ({len(complex_roots)} complex)"
                    print(f"    f={f}: {status}")
                    if not all_real and len(complex_roots) <= 6:
                        for cr in complex_roots:
                            print(f"           complex root: {cr[0]:.6f} + {cr[1]:.6f}i")
            else:
                # Without numpy, just report coefficients
                print(f"    f={f}: coeffs = {coeffs}")
        print()

    # Exact analysis for special cases
    print("  EXACT ANALYSIS of special cases:")
    print()
    print("  f=0: c_k^{(0,d)} = (-1)^{d-k} C(d,k)")
    print("    => P_{0,d}(x) = sum_k (-1)^{d-k} C(d,k) x^k = (-1)^d (1-x)^d")
    print("    => root x=1 with multiplicity d. ALL REAL (trivially).")
    print("    (numpy reports spurious imaginary parts due to repeated root instability)")
    print()
    print("  f=d: c_k^{(d,d)} = A(d+1, k) [Eulerian numbers]")
    print("    => P_{d,d}(x) = A_{d+1}(x), the Eulerian polynomial.")
    print("    => ALL REAL NEGATIVE roots (classical result).")
    print()
    print("  Odd f: c_k^{(f,d)} is ANTI-palindromic (c_k = -c_{d-k})")
    print("    => P_{f,d}(x) has factor (1+x) if d even, or (1-x) if d odd.")
    print("    => Always has x=-1 (or x=1) as a root.")
    print()

    # Refined analysis: for f=0, the polynomial is (1-x)^d, roots all real.
    # Mark these correctly.
    print("  CORRECTED SUMMARY (accounting for numerical artifacts):")
    for d in [4, 6, 8, 10]:
        print(f"    d={d}:")
        for f in range(d+1):
            coeffs = [inflated_eulerian(f, d, k) for k in range(d+1)]
            if f == 0:
                # (1-x)^d, all real
                status = "ALL REAL (exact: (1-x)^d)"
            elif f == d:
                # Eulerian polynomial, all real negative
                status = "ALL REAL (Eulerian polynomial A_{d+1}(x))"
            elif has_numpy:
                roots, all_real = find_polynomial_roots_numpy(coeffs)
                if roots is not None:
                    # Check more carefully: are "complex" roots just numerical noise near x=1?
                    max_imag = max(abs(r.imag) for r in roots) if roots else 0
                    if all_real:
                        status = "ALL REAL (numerical)"
                    elif max_imag < 0.05:
                        status = f"LIKELY ALL REAL (max |Im| = {max_imag:.2e}, near repeated root)"
                    else:
                        status = f"GENUINELY COMPLEX (max |Im| = {max_imag:.6f})"
                else:
                    status = "could not compute"
            else:
                status = "numpy unavailable"
            print(f"      f={f}: {status}")


# ====================================================================
# TASK 3: Log-concavity of |c_k^{(f,d)}|
# ====================================================================

def task3():
    print("\n" + "=" * 78)
    print("TASK 3: Log-concavity of |c_k^{(f,d)}|")
    print("=" * 78)
    print("  Checking: |c_k|^2 >= |c_{k-1}| * |c_{k+1}| for all k")
    print()

    for d in [4, 6, 8, 10]:
        print(f"  d = {d}:")
        for f in range(d+1):
            coeffs = [inflated_eulerian(f, d, k) for k in range(d+1)]
            abs_coeffs = [abs(c) for c in coeffs]

            log_concave = True
            violations = []
            for k in range(1, d):
                lhs = abs_coeffs[k] ** 2
                rhs = abs_coeffs[k-1] * abs_coeffs[k+1]
                if lhs < rhs:
                    log_concave = False
                    violations.append((k, lhs, rhs))

            status = "LOG-CONCAVE" if log_concave else f"FAILS at k={[v[0] for v in violations]}"
            print(f"    f={f}: |c_k| = {abs_coeffs}  {status}")
            if not log_concave:
                for k, lhs, rhs in violations:
                    print(f"           k={k}: |c_{k}|^2 = {lhs} < |c_{k-1}|*|c_{k+1}| = {rhs}")
        print()


# ====================================================================
# TASK 4: Positivity polytope for n=7
# ====================================================================

def task4():
    print("\n" + "=" * 78)
    print("TASK 4: Positivity polytope for n=7 (d=6)")
    print("=" * 78)
    print()
    print("  a_k(T) = A(7,k) + 2*c_k^{(4,6)}*t3 + 2*c_k^{(2,6)}*t5")
    print("                   + 2*c_k^{(0,6)}*t7 + 4*c_k^{(2,6)}*bc")
    print()
    print("  Constraint: a_k >= 0 for all k=0,...,6")
    print("  By palindromy, k and 6-k give same constraint.")
    print("  So independent constraints from k=0,1,2,3.")
    print()

    n = 7
    d = 6

    # Coefficients: A(7,k) + 2*c4_k*t3 + 2*c2_k*t5 + 2*c0_k*t7 + 4*c2_k*bc
    # where c4_k = c_k^{(4,6)}, c2_k = c_k^{(2,6)}, c0_k = c_k^{(0,6)}

    # Using Fraction for exact arithmetic
    constraints = []
    for k in range(d+1):
        Ak = Fraction(eulerian_number(n, k))
        c4 = Fraction(inflated_eulerian(4, d, k))
        c2 = Fraction(inflated_eulerian(2, d, k))
        c0 = Fraction(inflated_eulerian(0, d, k))
        # a_k = Ak + 2*c4*t3 + 2*c2*t5 + 2*c0*t7 + 4*c2*bc >= 0
        # i.e., 2*c4*t3 + 2*c2*t5 + 2*c0*t7 + 4*c2*bc >= -Ak
        # or:   2*c4*t3 + (2*c2)*t5 + 2*c0*t7 + (4*c2)*bc >= -Ak
        constraints.append({
            'k': k,
            'A': Ak,
            'coeff_t3': 2 * c4,
            'coeff_t5': 2 * c2,
            'coeff_t7': 2 * c0,
            'coeff_bc': 4 * c2,
        })

    print("  Linear constraints a_k >= 0:")
    print(f"  {'k':>3}  {'A(7,k)':>8}  {'2*c4':>6}*t3 + {'2*c2':>6}*t5 + {'2*c0':>6}*t7 + {'4*c2':>6}*bc >= 0")
    for c in constraints:
        print(f"  {c['k']:>3}  {int(c['A']):>8}  {int(c['coeff_t3']):>6}*t3 + {int(c['coeff_t5']):>6}*t5 + {int(c['coeff_t7']):>6}*t7 + {int(c['coeff_bc']):>6}*bc >= {int(-c['A'])}")

    print()
    print("  By palindromy, constraints for k and 6-k are identical.")
    print("  Independent constraints: k = 0, 1, 2, 3.")
    print()

    # The distinct constraints:
    print("  DISTINCT CONSTRAINTS:")
    for k in range(4):
        c = constraints[k]
        print(f"    k={k}: {int(c['A'])} + {int(c['coeff_t3'])}*t3 + {int(c['coeff_t5'])}*t5 + {int(c['coeff_t7'])}*t7 + {int(c['coeff_bc'])}*bc >= 0")

    print()
    print("  Note: k=0 gives 1 + 2*t3 + 2*t5 + 2*t7 + 4*bc >= 0")
    print("  Since all invariants (t3, t5, t7, bc) are non-negative, k=0 is always satisfied.")
    print("  Similarly k=6 (palindrome of k=0) is always satisfied.")
    print()

    # Check which constraints are non-trivial
    print("  NON-TRIVIAL CONSTRAINTS (those with negative coefficients):")
    for k in range(4):
        c = constraints[k]
        neg_terms = []
        for name, coeff in [('t3', c['coeff_t3']), ('t5', c['coeff_t5']),
                            ('t7', c['coeff_t7']), ('bc', c['coeff_bc'])]:
            if coeff < 0:
                neg_terms.append(f"{int(coeff)}*{name}")
        if neg_terms:
            print(f"    k={k}: negative terms: {', '.join(neg_terms)}")
        else:
            print(f"    k={k}: all coefficients non-negative => always satisfied")

    return constraints


# ====================================================================
# TASK 5 & 6: Actual invariant ranges and min a_k for n=7
# ====================================================================

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

def count_t3(A, n):
    """Count directed 3-cycles."""
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            t3 += 1
        if A[i][k] and A[k][j] and A[j][i]:
            t3 += 1
    return t3

def count_directed_cycles(A, n, cl):
    """Count directed cycles of length cl using DP."""
    if n < cl:
        return 0
    total = 0
    for verts in combinations(range(n), cl):
        sub = [[A[verts[i]][verts[j]] for j in range(cl)] for i in range(cl)]
        dp = [[0]*cl for _ in range(1 << cl)]
        dp[1][0] = 1
        for m in range(1, 1 << cl):
            for v in range(cl):
                if not (m & (1 << v)) or dp[m][v] == 0:
                    continue
                for u in range(cl):
                    if m & (1 << u):
                        continue
                    if sub[v][u]:
                        dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << cl) - 1
        total += sum(dp[full][v] for v in range(1, cl) if sub[v][0])
    return total

def count_bc(A, n):
    """Count pairs of vertex-disjoint directed 3-cycles (boundary cycles)."""
    cyc3_sets = []
    for t in combinations(range(n), 3):
        if (A[t[0]][t[1]] and A[t[1]][t[2]] and A[t[2]][t[0]]) or \
           (A[t[0]][t[2]] and A[t[2]][t[1]] and A[t[1]][t[0]]):
            cyc3_sets.append(set(t))
    return sum(1 for i in range(len(cyc3_sets))
               for j in range(i+1, len(cyc3_sets))
               if cyc3_sets[i].isdisjoint(cyc3_sets[j]))

def compute_ak(n, t3, t5, t7, bc):
    """Compute a_k for n=7 using the closed form."""
    d = n - 1
    result = []
    for k in range(n):
        ak = eulerian_number(n, k)
        ak += 2 * inflated_eulerian(4, d, k) * t3
        ak += 2 * inflated_eulerian(2, d, k) * t5
        ak += 2 * inflated_eulerian(0, d, k) * t7
        ak += 4 * inflated_eulerian(2, d, k) * bc
        result.append(ak)
    return result

def forward_edge_dist_dp(A, n):
    """Compute forward-edge distribution via DP."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            for fwd in range(n):
                c = dp.get((mask, v, fwd), 0)
                if c == 0:
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    new_fwd = fwd + A[v][u]
                    key = (mask | (1 << u), u, new_fwd)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    dist = defaultdict(int)
    for v in range(n):
        for fwd in range(n):
            dist[fwd] += dp.get((full, v, fwd), 0)
    return dict(dist)

def task5_6():
    print("\n" + "=" * 78)
    print("TASKS 5-6: Actual invariant ranges and min a_k for n=7")
    print("=" * 78)

    n = 7
    d = 6

    # Exhaustive enumeration of n=7 tournaments is 2^21 = 2M -- feasible but slow.
    # Use large random sample instead.
    num_samples = 5000
    print(f"\n  Sampling {num_samples} random n=7 tournaments...")

    all_invs = []
    min_ak = [float('inf')] * n
    max_ak = [0] * n
    min_ak_examples = [None] * n

    # Track which tournaments achieve a_k = 0
    zero_ak = [0] * n

    for trial in range(num_samples):
        A = random_tournament(n, trial * 7 + 13)
        t3 = count_t3(A, n)
        t5 = count_directed_cycles(A, n, 5)
        t7 = count_directed_cycles(A, n, 7)
        bc = count_bc(A, n)

        all_invs.append((t3, t5, t7, bc))

        ak = compute_ak(n, t3, t5, t7, bc)

        for k in range(n):
            if ak[k] < min_ak[k]:
                min_ak[k] = ak[k]
                min_ak_examples[k] = (t3, t5, t7, bc)
            if ak[k] > max_ak[k]:
                max_ak[k] = ak[k]
            if ak[k] == 0:
                zero_ak[k] += 1

    # Also check special tournaments
    # Transitive tournament: all invariants = 0
    ak_trans = compute_ak(n, 0, 0, 0, 0)
    print(f"\n  Transitive tournament (t3=t5=t7=bc=0):")
    print(f"    a_k = {ak_trans}")
    print(f"    = A(7,k) = {[eulerian_number(7,k) for k in range(7)]}")

    # Report invariant ranges
    t3_vals = sorted(set(inv[0] for inv in all_invs))
    t5_vals = sorted(set(inv[1] for inv in all_invs))
    t7_vals = sorted(set(inv[2] for inv in all_invs))
    bc_vals = sorted(set(inv[3] for inv in all_invs))

    print(f"\n  Invariant ranges (from {num_samples} samples):")
    print(f"    t3: min={min(t3_vals)}, max={max(t3_vals)}, distinct={len(t3_vals)}")
    print(f"    t5: min={min(t5_vals)}, max={max(t5_vals)}, distinct={len(t5_vals)}")
    print(f"    t7: min={min(t7_vals)}, max={max(t7_vals)}, distinct={len(t7_vals)}")
    print(f"    bc: min={min(bc_vals)}, max={max(bc_vals)}, distinct={len(bc_vals)}")

    print(f"\n  Forward-edge distribution ranges:")
    print(f"  {'k':>3}  {'min a_k':>10}  {'max a_k':>10}  {'#(a_k=0)':>10}  min example (t3,t5,t7,bc)")
    for k in range(n):
        ex = min_ak_examples[k]
        print(f"  {k:>3}  {min_ak[k]:>10}  {max_ak[k]:>10}  {zero_ak[k]:>10}  {ex}")

    print(f"\n  PALINDROMY CHECK: min_ak[k] == min_ak[{d}-k]?")
    for k in range(n):
        match = min_ak[k] == min_ak[d-k]
        print(f"    k={k}: min={min_ak[k]}, min[{d-k}]={min_ak[d-k]}, match={match}")

    # TASK 6: Is minimum exactly 0 or strictly positive?
    print(f"\n  TASK 6: Minimum of a_k(T)")
    print(f"  {'k':>3}  {'min a_k':>10}  {'status':>20}")
    for k in range(n):
        if min_ak[k] == 0:
            status = "ZERO achieved"
        elif min_ak[k] > 0:
            status = "STRICTLY POSITIVE"
        else:
            status = "NEGATIVE (BUG!)"
        print(f"  {k:>3}  {min_ak[k]:>10}  {status}")

    # Positivity constraint tightness
    print(f"\n  TIGHTNESS of positivity constraints:")
    print(f"  For each k with negative coefficients in the formula,")
    print(f"  the min a_k measures how tight the constraint is.")
    print(f"  If min a_k = 0, the constraint is TIGHT (active at some tournament).")
    print(f"  If min a_k > 0, there is SLACK in the constraint.")

    # Now compute the polytope constraints more carefully
    print(f"\n  POLYTOPE ANALYSIS:")
    print(f"  The non-trivial constraints are from k=1,2,3:")
    for k in [1, 2, 3]:
        Ak = eulerian_number(n, k)
        c4 = 2 * inflated_eulerian(4, d, k)
        c2_t5 = 2 * inflated_eulerian(2, d, k)
        c0 = 2 * inflated_eulerian(0, d, k)
        c2_bc = 4 * inflated_eulerian(2, d, k)
        print(f"    k={k}: {Ak} + {c4}*t3 + {c2_t5}*t5 + {c0}*t7 + {c2_bc}*bc >= 0")

    # The binding constraint: k=3 has the largest negative coefficients
    # k=3: 2416 - 160*t3 + 32*t5 - 40*t7 + 64*bc >= 0
    # i.e., 160*t3 + 40*t7 <= 2416 + 32*t5 + 64*bc
    print(f"\n  KEY CONSTRAINT (k=3):")
    print(f"    2416 - 160*t3 + 32*t5 - 40*t7 + 64*bc >= 0")
    print(f"    equivalently: 160*t3 + 40*t7 <= 2416 + 32*t5 + 64*bc")
    print(f"    or: 4*t3 + t7 <= 60.4 + 0.8*t5 + 1.6*bc")

    # Find tournaments closest to violating k=3
    print(f"\n  Tournaments closest to violating k=3 constraint:")
    scored = []
    for inv in all_invs:
        t3, t5, t7, bc = inv
        ak3 = 2416 - 160*t3 + 32*t5 - 40*t7 + 64*bc
        scored.append((ak3, inv))
    scored.sort()
    print(f"    Smallest a_3 values:")
    seen = set()
    for val, inv in scored[:10]:
        if inv not in seen:
            seen.add(inv)
            print(f"      a_3 = {val}, (t3, t5, t7, bc) = {inv}")

    # Extremal invariant combinations
    print(f"\n  EXTREMAL INVARIANT COMBINATIONS:")
    # Max t3
    max_t3_inv = max(all_invs, key=lambda x: x[0])
    print(f"    Max t3: {max_t3_inv}")
    ak = compute_ak(n, *max_t3_inv)
    print(f"      a_k = {ak}")
    # Max t7
    max_t7_inv = max(all_invs, key=lambda x: x[2])
    print(f"    Max t7: {max_t7_inv}")
    ak = compute_ak(n, *max_t7_inv)
    print(f"      a_k = {ak}")

    return all_invs


# ====================================================================
# ADDITIONAL: Verify a few against brute force
# ====================================================================

def verify_formula():
    """Quick brute-force verification at n=7."""
    print("\n" + "=" * 78)
    print("VERIFICATION: Brute-force check of a_k formula at n=7")
    print("=" * 78)

    n = 7
    d = 6
    all_match = True
    num_checks = 10

    for trial in range(num_checks):
        A = random_tournament(n, trial * 31 + 7)
        t3 = count_t3(A, n)
        t5 = count_directed_cycles(A, n, 5)
        t7 = count_directed_cycles(A, n, 7)
        bc = count_bc(A, n)

        predicted = compute_ak(n, t3, t5, t7, bc)
        actual = forward_edge_dist_dp(A, n)

        for k in range(n):
            if predicted[k] != actual.get(k, 0):
                print(f"  FAIL: trial={trial}, k={k}, pred={predicted[k]}, actual={actual.get(k,0)}")
                all_match = False

    print(f"  {num_checks} tournaments x {n} coefficients = {num_checks * n} checks: {'ALL PASS' if all_match else 'FAIL'}")
    return all_match


# ====================================================================
# MAIN
# ====================================================================

if __name__ == '__main__':
    task1()
    task2()
    task3()
    constraints = task4()

    # Verify formula first
    verify_formula()

    # Then compute ranges and min
    all_invs = task5_6()

    print("\n" + "=" * 78)
    print("DONE")
    print("=" * 78)
