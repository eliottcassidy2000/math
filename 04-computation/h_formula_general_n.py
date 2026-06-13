#!/usr/bin/env python3
"""
Does H = tr(c_0) + 3*(#3-cycles) generalize beyond n=5?

At n=5:
  H = tr(c_0) + tr(c_2)/4 + tr(c_4)/16
  tr(c_4) = 5*24 = 120
  tr(c_2) = 12*t_3 - 30
  => H = tr(c_0) + (12*t_3 - 30)/4 + 120/16
       = tr(c_0) + 3*t_3 - 7.5 + 7.5
       = tr(c_0) + 3*t_3

At n=3:
  H = tr(c_0) + tr(c_2)/4
  We found c_2 = 2*I for all n=3 tournaments.
  tr(c_2) = 6.
  => H = tr(c_0) + 6/4 = tr(c_0) + 1.5
  Need to check: is H = tr(c_0) + 1.5 always true at n=3?

At n=7: M(r) = c_0 + c_2*r^2 + c_4*r^4 + c_6*r^6
  c_6 = 6!*I = 720*I (conjectured).
  Is tr(c_4) determined by #3-cycles? Or something else?

Also: does the c_2 off-diagonal pattern c_2[a,b] = 2*(score_b - score_a)?

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def count_paths_weighted(A, verts, r_val, start=None, end=None):
    total = 0.0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        w = 1.0
        for i in range(len(p)-1):
            w *= r_val + (A[p[i]][p[i+1]] - 0.5)
        total += w
    return total

def transfer_matrix_r(A, r_val):
    n = len(A)
    M = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = 0.0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    ea = count_paths_weighted(A, S_verts, r_val, end=a)
                    bb = count_paths_weighted(A, R_verts, r_val, start=b)
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

def ham_path_count(A):
    n = len(A)
    return sum(1 for p in permutations(range(n))
               if all(A[p[i]][p[i+1]] == 1 for i in range(n-1)))

def count_3cycles(A):
    n = len(A)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                count += (A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])
    return count

def extract_even_poly(A, max_even_deg):
    """Extract even-power polynomial coefficients using Vandermonde-like system."""
    n_mat = len(A)
    num_coeffs = max_even_deg // 2 + 1
    u_samples = np.linspace(0, 1, num_coeffs + 2)
    r_samples = np.sqrt(u_samples)

    M_samples = [transfer_matrix_r(A, rv) for rv in r_samples]

    coeffs = {}
    for k in range(num_coeffs):
        coeffs[2*k] = np.zeros((n_mat, n_mat))

    for a in range(n_mat):
        for b in range(n_mat):
            vals = np.array([M_samples[i][a][b] for i in range(len(r_samples))])
            poly_coeffs = np.polyfit(u_samples, vals, num_coeffs - 1)
            for k in range(num_coeffs):
                coeffs[2*k][a][b] = poly_coeffs[num_coeffs - 1 - k]

    return coeffs

def tournament_canonical(A):
    n = len(A)
    min_adj = None
    for perm in permutations(range(n)):
        adj = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if min_adj is None or adj < min_adj:
            min_adj = adj
    return min_adj

# =====================================================================
# n=3: verify H = tr(c_0) + 3/2
# =====================================================================
print("=" * 70)
print("n=3: H FORMULA")
print("=" * 70)

# All tournaments on 3 vertices
all_n3 = []
for bits in range(8):
    A = [[0]*3 for _ in range(3)]
    pairs = [(0,1), (0,2), (1,2)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    all_n3.append(A)

for A in all_n3:
    H = ham_path_count(A)
    M0 = transfer_matrix_r(A, 0.0)
    tr_c0 = np.trace(M0)
    t3 = count_3cycles(A)
    scores = tuple(sorted(sum(row) for row in A))

    # At n=3: M(r) = c_0 + c_2*r^2. c_2 = 2I always.
    # H = tr(c_0) + tr(c_2)/4 = tr(c_0) + 6/4 = tr(c_0) + 3/2
    predicted = tr_c0 + 1.5
    print(f"  scores={scores}, H={H}, t3={t3}, tr(c0)={tr_c0:.4f}, "
          f"predicted H={predicted:.4f}, match={abs(predicted-H)<1e-6}")

# Check 3/2 = 3*t3 + something?
# At n=3: max t3 = 1 (the 3-cycle), min t3 = 0 (transitive).
# H for 3-cycle = 3, H for transitive = 1.
# tr(c_0) for 3-cycle = 1.5, for transitive = -0.5.
# 1.5 + 1.5 = 3 ✓; -0.5 + 1.5 = 1 ✓.
# But 3*t_3: 3*1=3, 3*0=0. So H = tr(c_0) + 3*t_3 doesn't work at n=3.
# Instead H = tr(c_0) + 3/2 always works at n=3 (since c_2 is universal).

print("\n  At n=3: H = tr(c_0) + 3/2 (universal c_2 = 2I)")
print("  NOT H = tr(c_0) + 3*t_3 (would give wrong answer for t_3=0)")

# =====================================================================
# n=4: check
# =====================================================================
print()
print("=" * 70)
print("n=4: H FORMULA")
print("=" * 70)

# At n=4 (even): M(r) = c_0 + c_2*r^2 (no c_4 since degree 3 is odd, c_4 would be degree 4 > n-1=3)
# Actually at n=4: max even power for diagonal = n-1=3 (ODD, so excluded).
# Max even power = n-2 = 2. So M(r) = c_0 + c_2*r^2 only.

# tr(M(r)) = 0 at even n (since sum (-1)^{pos(a,P)} = 0).
# So tr(c_0) + tr(c_2)/4 = 0 always at even n.
# H is NOT given by tr(M) at even n!

# Let's just check the formulas.
def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = []
    for a in range(n):
        for b in range(a):
            if a - b >= 2:
                tiles.append((a, b))
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A, tiles

n = 4
_, tiles = tiling_to_tournament(0, n)
m = len(tiles)

iso_classes = defaultdict(list)
for bits in range(2**m):
    A, _ = tiling_to_tournament(bits, n)
    canon = tournament_canonical(A)
    iso_classes[canon].append(bits)

print(f"n=4: {2**m} tilings, {len(iso_classes)} classes")
print(f"  Note: tr(M(r)) = 0 at even n (cannot extract H from trace)")

for canon in sorted(iso_classes.keys()):
    bits = iso_classes[canon][0]
    A, _ = tiling_to_tournament(bits, n)
    H = ham_path_count(A)
    t3 = count_3cycles(A)
    M0 = transfer_matrix_r(A, 0.0)
    Mhalf = transfer_matrix_r(A, 0.5)
    tr0 = np.trace(M0)
    tr_half = np.trace(Mhalf)

    # At even n: tr(M) = 0. Let's verify.
    print(f"  H={H}, t3={t3}, tr(M(0))={tr0:.4f}, tr(M(0.5))={tr_half:.4f}")


# =====================================================================
# n=5: the formula H = tr(c_0) + 3*t_3
# =====================================================================
print()
print("=" * 70)
print("GENERALIZATION: FORMULA STRUCTURE")
print("=" * 70)

# At odd n, M(r) = sum_{k even} c_k * r^k (k = 0, 2, ..., n-1)
# tr(M(1/2)) = H = sum_k tr(c_k) / 4^{k/2}

# The top term: tr(c_{n-1}) = n * (n-1)! = n!.
# So the top contribution to H = n! / 2^{n-1}.

# At n=3: n!/2^2 = 6/4 = 3/2. ✓ (this is the universal part)
# At n=5: n!/2^4 = 120/16 = 7.5. ✓

# At n=3: H = tr(c_0) + 3/2 (only two terms: c_0 and c_2=2I)
# At n=5: H = tr(c_0) + tr(c_2)/4 + 7.5 = tr(c_0) + 3*t_3

# The "3*t_3" comes from: tr(c_2)/4 + 7.5 = (12*t_3 - 30)/4 + 7.5 = 3*t_3 - 7.5 + 7.5 = 3*t_3

# Is there a general formula for tr(c_2)?
# At n=3: tr(c_2) = 6 for ALL tournaments. And t_3 can be 0 or 1.
# So tr(c_2) ≠ f(t_3) at n=3 (it's constant 6 regardless of t_3).

# At n=5: tr(c_2) = 12*t_3 - 30.
# Let's check if this generalizes as tr(c_2) = a*t_3 + b for some a,b at each n.

# For this, we need n=7 data. Too slow for exhaustive, but let's check a few.

print("\nFormula structure by n:")
print("  n=3 (odd): H = tr(c_0) + n!/2^{n-1} = tr(c_0) + 3/2")
print("  n=5 (odd): H = tr(c_0) + 3*t_3")
print(f"           = tr(c_0) + [tr(c_2)/4 + 120/16]")
print(f"           = tr(c_0) + [(12*t_3 - 30)/4 + 7.5]")
print(f"           = tr(c_0) + 3*t_3")

# =====================================================================
# Off-diagonal c_2: is it 2*(score_b - score_a)?
# =====================================================================
print()
print("=" * 70)
print("OFF-DIAGONAL c_2: SCORE DIFFERENCE PATTERN?")
print("=" * 70)

n = 5
_, tiles = tiling_to_tournament(0, n)
m_tiles = len(tiles)

iso_classes = defaultdict(list)
for bits in range(2**m_tiles):
    A, _ = tiling_to_tournament(bits, n)
    canon = tournament_canonical(A)
    if canon not in iso_classes:
        iso_classes[canon] = bits

print("  Testing c_2^offdiag[a,b] = 2*(score_b - score_a):")
match_count = 0
total_count = 0

for canon, bits in iso_classes.items():
    A, _ = tiling_to_tournament(bits, n)
    M0 = transfer_matrix_r(A, 0.0)
    Mhalf = transfer_matrix_r(A, 0.5)
    c2 = (Mhalf - M0) / 0.25

    scores = [sum(A[i]) for i in range(n)]

    for a in range(n):
        for b in range(n):
            if a == b: continue
            total_count += 1
            predicted = 2 * (scores[b] - scores[a])
            actual = c2[a][b]
            if abs(predicted - actual) < 1e-6:
                match_count += 1

print(f"  Match: {match_count}/{total_count}")

# That probably doesn't work. Let me check what the pattern actually is.
print("\n  Actual off-diagonal c_2 values vs score differences:")
# Use Paley (class 9) and a non-Paley
for name, bits_val in [("Paley", None), ("class 8", None)]:
    if name == "Paley":
        A = [[0]*5 for _ in range(5)]
        for i in range(5):
            for j in range(5):
                if i != j and (j-i)%5 in [1,4]:
                    A[i][j] = 1
    else:
        # Find class 8 representative
        for canon, bits_val in iso_classes.items():
            A_test, _ = tiling_to_tournament(bits_val, n)
            if ham_path_count(A_test) == 11:
                A = A_test
                break

    M0 = transfer_matrix_r(A, 0.0)
    Mhalf = transfer_matrix_r(A, 0.5)
    c2 = (Mhalf - M0) / 0.25
    scores = [sum(A[i]) for i in range(n)]

    print(f"\n  {name}: scores={scores}")
    for a in range(n):
        for b in range(n):
            if a == b: continue
            sd = scores[b] - scores[a]
            edge = A[a][b]
            print(f"    c2[{a},{b}]={c2[a][b]:+.0f}, score_diff={sd:+d}, "
                  f"A[{a},{b}]={edge}")

# =====================================================================
# The off-diagonal c_2 formula
# =====================================================================
print()
print("=" * 70)
print("OFF-DIAGONAL c_2: FORMULA SEARCH (n=5)")
print("=" * 70)

# Try: c_2[a,b] = alpha * (score_a + score_b - (n-1)) + beta * A[a][b] + gamma
# Or simpler patterns

for canon, bits_val in iso_classes.items():
    A, _ = tiling_to_tournament(bits_val, n)
    H = ham_path_count(A)
    if H not in [9, 11, 15]:  # focus on interesting classes
        continue

    M0 = transfer_matrix_r(A, 0.0)
    Mhalf = transfer_matrix_r(A, 0.5)
    c2 = (Mhalf - M0) / 0.25

    scores = [sum(A[i]) for i in range(n)]
    A_np = np.array(A, dtype=float)

    # Check if c_2[a,b] = f(score_a, score_b, A[a,b])
    # For Paley (regular, all scores = 2): c_2[a,b] should be constant
    # Actually Paley has c_2 = 0... hmm.

    # Let me just check: c_2 = alpha * (out-degree matrix - in-degree matrix)
    out_deg = np.diag(scores)
    J = np.ones((n,n))

    # c_2[a,b] = f(s_a, s_b) for off-diagonal
    # From the data: c_2[a,b] appears to be 2*(s_b - s_a) minus some correction
    # Let's tabulate

    print(f"\n  H={H}: c_2 off-diagonal pattern")
    for a in range(n):
        row = []
        for b in range(n):
            if a == b:
                row.append("  .  ")
            else:
                val = c2[a][b]
                row.append(f"{val:+5.0f}")
        print(f"    [{', '.join(row)}]  score_a={scores[a]}")
    print(f"    scores: {scores}")

print()
print("=" * 70)
print("DONE")
print("=" * 70)
