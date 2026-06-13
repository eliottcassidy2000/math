import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
c6_from_trace.py
kind-pasteur-2026-03-07-S39b

CONJECTURE: tr(A^6) = 6*c_6 + 3*c_3^2 for any tournament.

If true, then c_6 = (tr(A^6) - 3*c_3^2) / 6, computable in O(n^3)
since both tr(A^6) and c_3 are O(n^3) via matrix multiplication and
Moon's formula respectively.

This extends the trace-cycle identity to k=6 with a correction term.

PROOF IDEA: Non-simple closed walks of length 6 decompose into (3,3)-compound
walks: two 3-cycles sharing a junction vertex. The total number of such walks
equals 3*c_3^2. Here's why:

Each pair of directed 3-cycles (C1, C2) (not necessarily vertex-disjoint)
generates 6 non-simple closed walk sequences (6 rotations of the compound walk).
The number of ORDERED pairs of directed 3-cycles is (c_3)^2 (since each
directed 3-cycle has a unique direction in a tournament, and we count
ordered pairs). Wait, need to be more careful...

Actually: there are c_3 undirected 3-cycles, each with exactly 1 direction
in a tournament. So 1 directed 3-cycle per vertex triple. An ordered pair
of directed 3-cycles is (C1, C2). But a compound walk requires them to
share a specific junction vertex v and be "composed" at v.

Let me just verify the formula first.
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations
from math import comb
import random


def count_directed_k_cycles_dp(T, k):
    n = len(T)
    if k > n or k < 3:
        return 0
    count = 0
    for verts in combinations(range(n), k):
        v = list(verts)
        dp = [[0] * k for _ in range(1 << k)]
        dp[1][0] = 1
        for mask in range(1, 1 << k):
            for last in range(k):
                if dp[mask][last] == 0 or not (mask & (1 << last)):
                    continue
                for nxt in range(1, k):
                    if mask & (1 << nxt):
                        continue
                    if T[v[last]][v[nxt]]:
                        dp[mask | (1 << nxt)][nxt] += dp[mask][last]
        full = (1 << k) - 1
        for last in range(1, k):
            if T[v[last]][v[0]]:
                count += dp[full][last]
    return count


def matrix_power_trace(T, k):
    n = len(T)
    Ak = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(k):
        new = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for l in range(n):
                    new[i][j] += Ak[i][l] * T[l][j]
        Ak = new
    return sum(Ak[i][i] for i in range(n))


def c3_from_score(T):
    n = len(T)
    scores = [sum(T[i]) for i in range(n)]
    return comb(n, 3) - sum(comb(s, 2) for s in scores)


def c6_via_trace(T):
    """Compute c_6(T) = (tr(A^6) - 3*c_3^2) / 6 in O(n^3)."""
    n = len(T)
    if n < 6:
        return 0
    # c_3 via Moon's formula: O(n^2)
    scores = [sum(T[i]) for i in range(n)]
    c3 = comb(n, 3) - sum(comb(s, 2) for s in scores)

    # tr(A^6) via matrix multiplication: O(n^3)
    A2 = [[sum(T[i][k]*T[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    A3 = [[sum(A2[i][k]*T[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    A6 = [[sum(A3[i][k]*A3[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    tr6 = sum(A6[i][i] for i in range(n))

    correction = 3 * c3 * c3
    assert (tr6 - correction) % 6 == 0, f"tr6={tr6}, correction={correction}, remainder={(tr6 - correction) % 6}"
    return (tr6 - correction) // 6


# ============================================================
# Exhaustive verification
# ============================================================
print("=" * 70)
print("CONJECTURE: tr(A^6) = 6*c_6 + 3*c_3^2")
print("=" * 70)

for n in range(5, 8):
    m = n * (n - 1) // 2
    total = 1 << m
    if n >= 7:
        random.seed(42)
        sample = [random.randint(0, total - 1) for _ in range(2000)]
        sample_type = "sampled (2000)"
    else:
        sample = range(total)
        sample_type = "exhaustive"

    matches = 0
    tested = 0

    for bits in sample:
        T = tournament_from_bits(n, bits)
        tr6 = matrix_power_trace(T, 6)
        c3 = c3_from_score(T)
        c6 = count_directed_k_cycles_dp(T, 6) if n >= 6 else 0

        predicted = 6 * c6 + 3 * c3 * c3
        if tr6 == predicted:
            matches += 1
        elif tested < 5:
            print(f"  MISMATCH bits={bits}: tr6={tr6}, predicted={predicted}, "
                  f"c3={c3}, c6={c6}")
        tested += 1

    print(f"n={n} ({sample_type}, {tested}): "
          f"tr(A^6) = 6*c_6 + 3*c_3^2: {matches}/{tested} "
          f"({'PASS' if matches == tested else 'FAIL'})")


# ============================================================
# c_6 via trace: verification
# ============================================================
print("\n" + "=" * 70)
print("c_6 via trace formula: c_6 = (tr(A^6) - 3*c_3^2) / 6")
print("=" * 70)

for n in [6]:
    m = n * (n - 1) // 2
    all_match = True
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        c6_dp = count_directed_k_cycles_dp(T, 6)
        c6_tr = c6_via_trace(T)
        if c6_dp != c6_tr:
            print(f"  n={n}, bits={bits}: MISMATCH dp={c6_dp}, trace={c6_tr}")
            all_match = False
    if all_match:
        print(f"  n={n}: c6_via_trace matches for ALL {1 << m} tournaments")

# n=7 sampled
random.seed(42)
n = 7
m = n * (n - 1) // 2
sample = [random.randint(0, (1 << m) - 1) for _ in range(500)]
all_match = True
for bits in sample:
    T = tournament_from_bits(n, bits)
    c6_dp = count_directed_k_cycles_dp(T, 7 if n >= 7 else 6)
    # Hmm, need c6 at n=7
    c6_dp = count_directed_k_cycles_dp(T, 6)
    c6_tr = c6_via_trace(T)
    if c6_dp != c6_tr:
        print(f"  n=7, bits={bits}: MISMATCH dp={c6_dp}, trace={c6_tr}")
        all_match = False
if all_match:
    print(f"  n=7: c6_via_trace matches for ALL 500 sampled tournaments")


# ============================================================
# Now try to extend to k=7 and k=8
# ============================================================
print("\n" + "=" * 70)
print("EXTENDING: tr(A^7) correction formula")
print("=" * 70)

# For k=7: compound walks are (3,4)-type.
# Try: tr(A^7) = 7*c_7 + alpha * c_3 * c_4
# for some constant alpha.

for n in [7]:
    random.seed(42)
    m = n * (n - 1) // 2
    sample = [random.randint(0, (1 << m) - 1) for _ in range(300)]

    # Collect (excess, c3, c4) triples
    data = []
    for bits in sample:
        T = tournament_from_bits(n, bits)
        tr7 = matrix_power_trace(T, 7)
        c3 = c3_from_score(T)
        c4 = count_directed_k_cycles_dp(T, 4)
        c5 = count_directed_k_cycles_dp(T, 5)
        c7 = count_directed_k_cycles_dp(T, 7)
        excess = tr7 - 7 * c7
        data.append({'excess': excess, 'c3': c3, 'c4': c4, 'c5': c5, 'c7': c7})

    # Test: excess = alpha * c3 * c4
    # Find the ratio excess / (c3 * c4) for each tournament
    ratios = set()
    for d in data:
        if d['c3'] * d['c4'] > 0:
            ratios.add(d['excess'] / (d['c3'] * d['c4']))
    print(f"  n=7: excess/(c3*c4) has {len(ratios)} distinct values "
          f"(need 1 for linear formula)")
    if len(ratios) <= 5:
        print(f"    values: {sorted(ratios)}")

    # Maybe excess = a*c3*c4 + b*c3 + c*c4 + d for some (a,b,c,d)?
    # Use least squares to find coefficients
    import numpy as np
    X = np.array([[d['c3']*d['c4'], d['c3'], d['c4'], 1] for d in data])
    y = np.array([d['excess'] for d in data])

    try:
        coeffs, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
        max_err = max(abs(y - X @ coeffs))
        print(f"  Fit: excess = {coeffs[0]:.4f}*c3*c4 + {coeffs[1]:.4f}*c3 + "
              f"{coeffs[2]:.4f}*c4 + {coeffs[3]:.4f}")
        print(f"  Max error: {max_err:.4f}")
        if max_err < 0.1:
            print(f"  EXACT FIT!")
            # Round coefficients
            a, b, c_coeff, d_coeff = [round(x) for x in coeffs]
            print(f"  Rounded: excess = {a}*c3*c4 + {b}*c3 + {c_coeff}*c4 + {d_coeff}")
            # Verify
            exact_match = 0
            for d in data:
                pred = a*d['c3']*d['c4'] + b*d['c3'] + c_coeff*d['c4'] + d_coeff
                if pred == d['excess']:
                    exact_match += 1
            print(f"  Exact match: {exact_match}/{len(data)}")

        # Try more variables: c3^2, c4^2, c5, c3*c5, etc.
        X2 = np.array([[d['c3']*d['c4'], d['c3']**2, d['c4']**2, d['c5'],
                        d['c3']*d['c5'], d['c3'], d['c4'], 1] for d in data])
        coeffs2, _, _, _ = np.linalg.lstsq(X2, y, rcond=None)
        max_err2 = max(abs(y - X2 @ coeffs2))
        print(f"\n  Extended fit: excess = {coeffs2[0]:.4f}*c3*c4 + {coeffs2[1]:.4f}*c3^2 + "
              f"{coeffs2[2]:.4f}*c4^2 + {coeffs2[3]:.4f}*c5 + "
              f"{coeffs2[4]:.4f}*c3*c5 + {coeffs2[5]:.4f}*c3 + "
              f"{coeffs2[6]:.4f}*c4 + {coeffs2[7]:.4f}")
        print(f"  Max error: {max_err2:.4f}")
        if max_err2 < 0.1:
            rounded = [round(x) for x in coeffs2]
            print(f"  Rounded: {rounded}")
    except Exception as e:
        print(f"  numpy not available or error: {e}")


# ============================================================
# k=8: compound walks are (3,5), (4,4), (3,3,2=impossible)
# Actually: (3,5), (4,4), and (3,3,2) where 2 is impossible.
# Wait: length 8 = 3+5 or 4+4. Also 3+3+2 impossible (no 2-cycles).
# And 3+3+?? where the third sub-walk has length 2: impossible.
# So compound walks of length 8: (3,5) and (4,4) types.
# ============================================================
print("\n" + "=" * 70)
print("EXTENDING: tr(A^8) correction formula")
print("=" * 70)

for n in [8]:
    random.seed(42)
    m = n * (n - 1) // 2
    sample = [random.randint(0, (1 << m) - 1) for _ in range(200)]

    data = []
    for bits in sample:
        T = tournament_from_bits(n, bits)
        tr8 = matrix_power_trace(T, 8)
        c3 = c3_from_score(T)
        c4 = count_directed_k_cycles_dp(T, 4)
        c5 = count_directed_k_cycles_dp(T, 5)
        c8 = count_directed_k_cycles_dp(T, 8)
        excess8 = tr8 - 8 * c8
        data.append({'excess': excess8, 'c3': c3, 'c4': c4, 'c5': c5, 'c8': c8})

    # Try fitting: excess8 = a*c3*c5 + b*c4^2 + c*c3^2*c3 + ...
    try:
        X = np.array([[d['c3']*d['c5'], d['c4']**2, d['c3']**2, d['c3']*d['c4'],
                       d['c3'], d['c4'], d['c5'], 1] for d in data])
        y = np.array([d['excess'] for d in data])
        coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        max_err = max(abs(y - X @ coeffs))
        labels = ['c3*c5', 'c4^2', 'c3^2', 'c3*c4', 'c3', 'c4', 'c5', '1']
        print(f"  k=8 fit:")
        for label, coeff in zip(labels, coeffs):
            if abs(coeff) > 0.01:
                print(f"    {coeff:+.4f} * {label}")
        print(f"  Max error: {max_err:.4f}")
        if max_err < 0.1:
            rounded = [round(x) for x in coeffs]
            print(f"  Rounded: excess = ", end="")
            terms = []
            for label, coeff in zip(labels, rounded):
                if coeff != 0:
                    terms.append(f"{coeff}*{label}")
            print(" + ".join(terms))
    except Exception as e:
        print(f"  Error: {e}")


print("\nDone.")
