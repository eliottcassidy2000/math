#!/usr/bin/env python3
"""
DEGREE-4 FOURIER STRUCTURE AT n=9: INVESTIGATION

At n=7, the degree-4 Fourier space of H(T) is 2-dimensional,
spanned by Type P (5-vertex spanning paths) and Type Q (disjoint P2 pairs).

Question: What is the dimension at n=9?

Approach:
1. Compute w_4(T) for many random tournaments on 9 vertices
   (w_4 is the degree-4 homogeneous piece in s_e variables)
2. Extract the degree-4 Fourier coefficients (Walsh-Hadamard on 36 edges)
3. Determine the rank of the coefficient vectors

Since C(36,4) = 58905 monomials, we need many tournaments to probe
the full space. But the rank should be small.

Alternative approach: compute w_4 directly as a multilinear polynomial,
then find its support structure.

opus-2026-03-07-S35
"""
from itertools import combinations, permutations
from collections import defaultdict
import numpy as np
from numpy.linalg import svd
import random
import time

n = 9
edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
m = len(edges)  # 36
edge_idx = {e: k for k, e in enumerate(edges)}

print(f"n = {n}, m = {m} edges, C({m},4) = {m*(m-1)*(m-2)*(m-3)//24} degree-4 monomials")

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

def s_vec(A):
    """Centered edge variables: s_e = A[i][j] - 0.5 for i<j."""
    return [A[i][j] - 0.5 for i, j in edges]

def compute_w4_coefficients(A):
    """
    Compute the degree-4 part of W(r) = sum_P prod(r + s_{P(i)->P(i+1)}).

    w_4 = sum over permutations P, sum over 4-element subsets S of {0,...,7},
          product_{i in S} s_{P(i)->P(i+1)}

    For each permutation P of [0..8], the path has 8 edges.
    Choose 4 of these 8 edges, take product of their s-values.
    Sum over all perms and all 4-subsets.

    This gives w_4 as a polynomial in s_e variables.
    We return the coefficient of each degree-4 monomial.
    """
    coeffs = defaultdict(float)
    count = 0
    for perm in permutations(range(n)):
        # Edges of this path
        path_edges = []
        path_signs = []
        for i in range(n - 1):
            u, v = perm[i], perm[i+1]
            if u < v:
                path_edges.append((u, v))
                path_signs.append(1.0)
            else:
                path_edges.append((v, u))
                path_signs.append(-1.0)

        # Choose 4 of the 8 path edges
        for subset in combinations(range(n - 1), 4):
            mono_edges = tuple(sorted(path_edges[j] for j in subset))
            sign = 1.0
            for j in subset:
                sign *= path_signs[j]

            # But monomials can have repeated edges if path visits same
            # undirected edge twice -- impossible since P is a permutation
            # and edges are between consecutive vertices in the path.
            # Each undirected edge can appear at most once.
            # But we need to check: could two positions give the same undirected edge?
            # Since P is a permutation, consecutive pairs (P[i],P[i+1]) are all
            # distinct as ordered pairs, hence distinct as unordered pairs too.

            # However, the mono_edges tuple might have duplicates if
            # the SAME undirected edge appears at two different positions
            # No: each position in the path uses a DIFFERENT pair of vertices
            # (since P is injective), so undirected edges are all distinct.

            coeffs[mono_edges] += sign

        count += 1
        if count % 50000 == 0:
            print(f"  Processed {count} permutations...")

    return coeffs

# This is too slow: 9! = 362880 permutations, each with C(8,4)=70 subsets
# = 25.4M operations. Let's try it but with timing.

print("\nMethod 1: Direct coefficient extraction (may be slow)...")
t0 = time.time()

# Actually let's first try a sampling approach to determine dimension
print("First: sampling approach to determine dimension of degree-4 space")
print("Computing w_4 for random tournaments via polynomial interpolation...\n")

def compute_W(A, r):
    """W(r) via DP."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            val = dp.get((mask, v), 0)
            if val == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                dp[(mask | (1 << u), u)] = (
                    dp.get((mask | (1 << u), u), 0) + val * (r + A[v][u] - 0.5)
                )
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def extract_w_coeffs(A):
    """Extract W(r) polynomial coefficients by interpolation."""
    points = np.linspace(0.1, 5.0, n + 5)
    vals = [compute_W(A, r) for r in points]
    coeffs = np.polyfit(points, vals, n - 1)
    # coeffs[0] is highest power, coeffs[-1] is constant
    w = {}
    for i, c in enumerate(coeffs):
        w[n - 1 - i] = c
    return w

# Compute cycle invariants
def count_directed_cycles(A, L):
    """Count directed cycles of length L."""
    total = 0
    for verts in combinations(range(n), L):
        sub = [[A[verts[i]][verts[j]] for j in range(L)] for i in range(L)]
        dp = [[0] * L for _ in range(1 << L)]
        dp[1][0] = 1
        for mask in range(1, 1 << L):
            for v in range(L):
                if not (mask & (1 << v)) or dp[mask][v] == 0:
                    continue
                for u in range(L):
                    if mask & (1 << u):
                        continue
                    if sub[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        full = (1 << L) - 1
        total += sum(dp[full][v] for v in range(1, L) if sub[v][0])
    return total

def count_disjoint_3cycle_pairs(A):
    """Count pairs of vertex-disjoint directed 3-cycles."""
    tri = []
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            tri.append(frozenset([i, j, k]))
        if A[i][k] and A[k][j] and A[j][i]:
            tri.append(frozenset([i, j, k]))
    count = 0
    for a in range(len(tri)):
        for b in range(a + 1, len(tri)):
            if len(tri[a] & tri[b]) == 0:
                count += 1
    return count

def count_disjoint_3cycle_triples(A):
    """Count triples of pairwise vertex-disjoint directed 3-cycles."""
    tri = []
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            tri.append(frozenset([i, j, k]))
        if A[i][k] and A[k][j] and A[j][i]:
            tri.append(frozenset([i, j, k]))
    count = 0
    for a in range(len(tri)):
        for b in range(a + 1, len(tri)):
            if len(tri[a] & tri[b]) > 0:
                continue
            for c in range(b + 1, len(tri)):
                if len(tri[a] & tri[c]) == 0 and len(tri[b] & tri[c]) == 0:
                    count += 1
    return count

# Step 1: Compute w_4 and cycle invariants for many tournaments
N = 200
print(f"Computing w_4 and cycle invariants for {N} random tournaments...")

w4_vals = np.zeros(N)
t3_vals = np.zeros(N)
t5_vals = np.zeros(N)
t7_vals = np.zeros(N)
t9_vals = np.zeros(N)
a33_vals = np.zeros(N)

t0 = time.time()
for trial in range(N):
    A = random_tournament(n, trial)
    w = extract_w_coeffs(A)
    w4_vals[trial] = w[4]
    t3_vals[trial] = count_directed_cycles(A, 3)
    t5_vals[trial] = count_directed_cycles(A, 5)
    a33_vals[trial] = count_disjoint_3cycle_pairs(A)

    if trial < 5 or trial % 50 == 0:
        print(f"  trial {trial}: w4={w4_vals[trial]:.1f}, t3={t3_vals[trial]:.0f}, "
              f"t5={t5_vals[trial]:.0f}, a33={a33_vals[trial]:.0f}")

t1 = time.time()
print(f"\nDone in {t1-t0:.1f}s")

# Step 2: Verify the formula w_4 = 40320 - 4200*t3 + 240*t5 + 480*a33
print("\n" + "="*70)
print("VERIFYING w_4 FORMULA")
print("="*70)

predicted = 40320 - 4200*t3_vals + 240*t5_vals + 480*a33_vals
err = np.max(np.abs(w4_vals - predicted))
print(f"w_4 = 40320 - 4200*t3 + 240*t5 + 480*a33")
print(f"Max error: {err:.6f}")

# Step 3: Extract degree-4 parts
# t3 is degree 2, so it has NO degree-4 component.
# The constant 40320 is degree 0.
# So the degree-4 part of w_4 is: 240*[deg-4 of t5] + 480*[deg-4 of a33]
#
# But wait: t3 has max degree 2 but the formula has -4200*t3 in a degree-4 polynomial.
# How? Because w_4 IS homogeneous degree 4. The "t3" and constant terms must
# combine to produce degree-4 terms when expanded in s_e variables.
# Actually NO: t3 as a function of s_e is degree <= 2 (in fact degree 0 + degree 2).
# So 40320 - 4200*t3 is degree 0 + degree 2, which contradicts w_4 being degree 4.
#
# Resolution: w_4 is homogeneous degree 4, but the formula "40320 - 4200*t3 + ..."
# must be read as: the degree-4 homogeneous piece of the RHS.
# Actually, I think the formula gives w_4 as a function of tournament invariants,
# where the invariants themselves have mixed degrees. The specific combination
# 40320 - 4200*t3 + 240*t5 + 480*a33 is EXACTLY degree 4 because:
#   40320 = degree 0, but this gets cancelled by the degree-0 parts of the other terms
#
# More precisely, the formula should be read as:
#   w_4(T) = 40320 - 4200*t3(T) + 240*t5(T) + 480*a33(T)
# where each side is evaluated on the tournament T.
# w_4 happens to be degree 4 in s_e, but the formula expresses it using
# invariants t3, t5, a33 which individually have mixed degrees.
# The cancellation of lower degrees happens automatically.

# Let's verify: compute the degree-4 parts
# [deg-4 of t5]: t5 minus its lower-degree parts
# t5 = E[t5] + c*[deg-2 of t3] + [deg-4 of t5]
# From the n=7 analysis: proportionality constant for t5 at degree 2 is c_5=3

# For n=9, we need the proportionality constant. Let's compute it.
print("\n" + "="*70)
print("DEGREE-2 PROPORTIONALITY CONSTANTS AT n=9")
print("="*70)

# t3 centered
t3c = t3_vals - t3_vals.mean()

# t5 = E[t5] + c5*t3c + [deg-4 of t5]
# a33 = E[a33] + ca*t3c + [deg-4 of a33]
# Regress t5 on t3c, a33 on t3c

from numpy.linalg import lstsq

c5_fit, _, _, _ = lstsq(t3c.reshape(-1,1), (t5_vals - t5_vals.mean()).reshape(-1,1), rcond=None)
c5 = c5_fit[0,0]

ca_fit, _, _, _ = lstsq(t3c.reshape(-1,1), (a33_vals - a33_vals.mean()).reshape(-1,1), rcond=None)
ca = ca_fit[0,0]

print(f"  c_5 (t5 ~ c5*t3): {c5:.6f}")
print(f"  c_a33 (a33 ~ ca*t3): {ca:.6f}")

# Extract degree-4 residuals
t5_d4 = t5_vals - t5_vals.mean() - c5 * t3c
a33_d4 = a33_vals - a33_vals.mean() - ca * t3c

print(f"\n  t5_d4 range: [{t5_d4.min():.2f}, {t5_d4.max():.2f}], std: {t5_d4.std():.2f}")
print(f"  a33_d4 range: [{a33_d4.min():.2f}, {a33_d4.max():.2f}], std: {a33_d4.std():.2f}")

# Check: is w_4 degree-4 = 240*t5_d4 + 480*a33_d4?
# Actually w_4 is degree 4 homogeneous, so:
# w_4 = 240*[deg-4 of t5] + 480*[deg-4 of a33]
#      + (-4200)*[deg-2 of t3] which is degree 2... contradiction!
#
# Unless 40320 - 4200*E[t3] + 240*E[t5] + 480*E[a33] = 0 (degree 0 cancelled)
# and -4200*[deg-2 of t3] + 240*[deg-2 of t5] + 480*[deg-2 of a33] = 0 (degree 2 cancelled)

# Check degree-0 cancellation
Et3 = t3_vals.mean()
Et5 = t5_vals.mean()
Ea33 = a33_vals.mean()
deg0 = 40320 - 4200*Et3 + 240*Et5 + 480*Ea33
print(f"\n  Degree-0 check: 40320 - 4200*E[t3] + 240*E[t5] + 480*E[a33]")
print(f"    = 40320 - 4200*{Et3:.2f} + 240*{Et5:.2f} + 480*{Ea33:.2f}")
print(f"    = {deg0:.2f} (should be 0 if w_4 has no degree-0 part)")

# Check degree-2 cancellation: -4200*[t3_d2] + 240*c5*[t3_d2] + 480*ca*[t3_d2] should = 0
# i.e., -4200 + 240*c5 + 480*ca = 0
deg2_coeff = -4200 + 240*c5 + 480*ca
print(f"\n  Degree-2 check: -4200 + 240*c5 + 480*ca = -4200 + 240*{c5:.4f} + 480*{ca:.4f}")
print(f"    = {deg2_coeff:.4f} (should be 0)")

# So w_4's degree-4 part = 240*t5_d4 + 480*a33_d4
# The dimension of the degree-4 space is determined by how many
# linearly independent degree-4 invariants exist.

print("\n" + "="*70)
print("SVD ANALYSIS: DIMENSION OF DEGREE-4 SPACE")
print("="*70)

# Are t5_d4 and a33_d4 linearly independent?
X2 = np.column_stack([t5_d4, a33_d4])
U2, S2, Vt2 = svd(X2, full_matrices=False)
print(f"\n  Two-variable SVD (t5_d4, a33_d4):")
print(f"    Singular values: {S2}")
print(f"    Ratio S[1]/S[0] = {S2[1]/S2[0]:.6f}")
print(f"    => {'2D (independent)' if S2[1]/S2[0] > 0.01 else '1D (dependent)'}")

# Now: at n=9, t7 also has degree-4 components.
# t7 has max degree 7, so it has components at degrees 1,3,5,7 (odd!).
# Wait: t7 counts directed 7-cycles. In s_e variables, each 7-cycle
# has 7 edges, giving a degree-7 monomial. But s_e^2 = 1/4 (since s_e = ±1/2),
# so the multilinear expansion reduces to lower degrees.
# Actually in the Walsh basis, all functions on {0,1}^m are multilinear.
# So t7 as a multilinear polynomial has max degree 7.
# Its degree-4 part is nonzero in general.

# Similarly, t9 has max degree 9, hence degree-4 part too.
# And alpha_{3,5} (disjoint 3-cycle + 5-cycle pair) has max degree 8,
# hence degree-4 part.
# And alpha_{3,3,3} (three disjoint 3-cycles) has max degree 9,
# hence degree-4 part.

# But from the formula, only t5 and a33 appear with nonzero coefficients
# in w_4. This suggests that:
# [deg-4 of t7], [deg-4 of t9], [deg-4 of a35], [deg-4 of a333]
# are all in the span of {[deg-4 of t5], [deg-4 of a33]}.

# Let's compute t7 too and check
print("\n\nComputing t7 for all tournaments (slower)...")
t0 = time.time()
for trial in range(N):
    A = random_tournament(n, trial)
    t7_vals[trial] = count_directed_cycles(A, 7)
    if trial % 50 == 0:
        print(f"  trial {trial}: t7={t7_vals[trial]:.0f}")
t1 = time.time()
print(f"Done in {t1-t0:.1f}s")

# t7 degree-4 part: t7 = E[t7] + c7_2*t3c + c7_4_1*t5_d4 + c7_4_2*a33_d4 + higher
# First remove degree 0 and 2
t7c = t7_vals - t7_vals.mean()
c7_2, _, _, _ = lstsq(t3c.reshape(-1,1), t7c.reshape(-1,1), rcond=None)
c7_2 = c7_2[0,0]
t7_d46 = t7c - c7_2 * t3c  # degrees 4, 6

print(f"\n  c7_2 (t7 ~ c7*t3): {c7_2:.6f}")

# Now: does t7_d46 live in span of {t5_d4, a33_d4}?
# Not exactly, because t7_d46 includes degree-6 parts too.
# But degree-4 part of t7 should be in the span.
# We can regress t7_d46 on t5_d4 and a33_d4 and check residual.

coeffs_t7, _, _, _ = lstsq(np.column_stack([t5_d4, a33_d4]), t7_d46, rcond=None)
t7_residual = t7_d46 - coeffs_t7[0]*t5_d4 - coeffs_t7[1]*a33_d4
print(f"  t7_d46 = {coeffs_t7[0]:.6f}*t5_d4 + {coeffs_t7[1]:.6f}*a33_d4 + residual")
print(f"  Residual std: {t7_residual.std():.4f} (should be nonzero if degree-6 present)")
print(f"  t7_d46 std: {np.std(t7_d46):.4f}")
print(f"  R^2 = {1 - t7_residual.var()/t7_d46.var():.6f}")

# The residual is the degree-6+ part of t7. If degree-4 space is 2D,
# R^2 won't be 1 because of the degree-6 contamination.
# To properly isolate degree 4, we need to work with the actual w_4 coefficient.

# Better approach: use w_4 directly since it's purely degree 4.
# w_4 = 240*[deg-4 of t5] + 480*[deg-4 of a33]
# So w_4 lives in a 2D space spanned by deg-4 parts of t5 and a33.
# But we should also include OTHER degree-4 invariants to see if
# the space is truly larger.

# Key insight: w_4/4^2 = w_4/16 is what appears in the OCF sum.
# But w_4 itself IS the degree-4 piece of W(r) at r^4.
# w_4 is a degree-4 multilinear polynomial in 36 variables.
# Its dimension as a subspace of degree-4 monomials depends on
# what it equals.

# Actually the question is: what is dim of the span of
# {[deg-4 of f] : f is a tournament invariant relevant to OCF}?
# The answer is: it's the dimension of the space of degree-4
# multilinear functions that are tournament invariants.

# For the OCF, the degree-4 identity is:
#   w_4 / 16 = 2*[deg-4 of alpha_1] + 4*[deg-4 of alpha_2] + 8*[deg-4 of alpha_3]
# where alpha_k = sum of independent sets of size k in Omega(T).

# At n=9: alpha_1 = t3+t5+t7+t9, alpha_2 = a33+a35, alpha_3 = a333
# [deg-4 of alpha_1] = [deg-4 of t5] + [deg-4 of t7] + [deg-4 of t9]
#   (t3 has max degree 2, so no degree-4 part; t5's deg-4 is its top degree)
# [deg-4 of alpha_2] = [deg-4 of a33] + [deg-4 of a35]
# [deg-4 of alpha_3] = [deg-4 of a333]

# So there are potentially 5 degree-4 invariants:
# [deg-4 of t5], [deg-4 of t7], [deg-4 of t9], [deg-4 of a35], [deg-4 of a333]

# From the OCF formula at n=9:
# w_4 = 40320 - 4200*t3 + 240*t5 + 480*a33
# This has ONLY t5 and a33, not t7, t9, a35, a333.
# This means the degree-4 parts of t7, t9, a35, a333 are in
# span{deg-4 of t5, deg-4 of a33}.
# So the degree-4 space dimension is AT MOST 2.

# But is it exactly 2? Let's verify by computing w_4 on many tournaments
# and checking it can't be expressed as c*t5_d4 alone.

print("\n" + "="*70)
print("DIRECT DIMENSION CHECK VIA w_4")
print("="*70)

# w_4_centered = w4_vals - w4_vals.mean()
# This should be in span{t5_d4, a33_d4}
w4c = w4_vals - w4_vals.mean()

# Check if w4c is proportional to t5_d4 (1D) or needs both (2D)
c_fit1, _, _, _ = lstsq(t5_d4.reshape(-1,1), w4c.reshape(-1,1), rcond=None)
res1 = w4c - c_fit1[0,0]*t5_d4
print(f"  1D fit (w4 ~ c*t5_d4): c={c_fit1[0,0]:.4f}, residual std={res1.std():.4f}")

c_fit2, _, _, _ = lstsq(np.column_stack([t5_d4, a33_d4]), w4c, rcond=None)
res2 = w4c - c_fit2[0]*t5_d4 - c_fit2[1]*a33_d4
print(f"  2D fit (w4 ~ c1*t5_d4 + c2*a33_d4): c1={c_fit2[0]:.4f}, c2={c_fit2[1]:.4f}")
print(f"    residual std={res2.std():.6f}")
print(f"    R^2 = {1 - res2.var()/w4c.var():.10f}")

# The predicted coefficients should be 240 and 480
print(f"    Expected: c1=240, c2=480")

print("\n" + "="*70)
print("SUPPORT STRUCTURE: CLASSIFYING DEGREE-4 MONOMIALS AT n=9")
print("="*70)

# At n=9, a degree-4 monomial is a product of 4 edge variables.
# The support graph has 4 edges. Possible structures:
# (by number of vertices and degree sequence)

def classify_support(edge_list):
    """Classify the support graph of a 4-edge monomial."""
    verts = set()
    deg = defaultdict(int)
    for u, v in edge_list:
        verts.add(u)
        verts.add(v)
        deg[u] += 1
        deg[v] += 1
    nv = len(verts)
    deg_seq = tuple(sorted([deg[v] for v in verts], reverse=True))
    return nv, deg_seq

# Enumerate all possible support types
type_counts = defaultdict(int)
for mono in combinations(range(m), 4):
    elist = [edges[i] for i in mono]
    gtype = classify_support(elist)
    type_counts[gtype] += 1

print(f"\nSupport graph types for C({m},4)={sum(type_counts.values())} degree-4 monomials:")
print(f"{'Type':>25} | {'Count':>6} | Description")
print("-"*65)
for gtype in sorted(type_counts.keys()):
    nv, deg_seq = gtype
    # Describe the graph
    if nv == 4 and deg_seq == (3, 3, 1, 1):
        desc = "K4 minus edge (triangle + pendant)"
    elif nv == 4 and deg_seq == (2, 2, 2, 2):
        desc = "4-cycle"
    elif nv == 5 and deg_seq == (2, 2, 2, 1, 1):
        desc = "P5 path (spanning path on 5 verts)"
    elif nv == 5 and deg_seq == (2, 2, 1, 1, 2):
        desc = "?"
    elif nv == 5 and deg_seq == (3, 1, 1, 1, 2):
        desc = "star+edge"
    elif nv == 5 and deg_seq == (4, 1, 1, 1, 1):
        desc = "K1,4 star"
    elif nv == 6 and deg_seq == (2, 2, 1, 1, 1, 1):
        desc = "P3+P3 (disjoint P2 pairs = P_2 + P_2)"
    elif nv == 6 and deg_seq == (2, 1, 1, 1, 1, 2):
        desc = "?"
    elif nv == 7 and deg_seq == (1, 1, 1, 1, 1, 1, 1, 1):
        desc = "perfect matching? (impossible with 4 edges on 7 verts)"
    elif nv == 8 and deg_seq == (1, 1, 1, 1, 1, 1, 1, 1):
        desc = "4 disjoint edges (perfect matching on 8 verts)"
    else:
        desc = ""
    print(f"  {nv}v {deg_seq} | {type_counts[gtype]:>6} | {desc}")

# Now compute which types have nonzero coefficients in w_4
# This requires the full coefficient extraction, which is expensive.
# Let's do it.

print("\n" + "="*70)
print("EXTRACTING w_4 COEFFICIENTS (full enumeration of 9! permutations)")
print("="*70)
print("This will take a while (~25M operations)...")

t0 = time.time()
w4_coeffs = defaultdict(float)
perm_count = 0

for perm in permutations(range(n)):
    # Edges of this Hamiltonian path
    path_edges = []
    path_signs = []
    for i in range(n - 1):
        u, v = perm[i], perm[i+1]
        if u < v:
            path_edges.append(edge_idx[(u, v)])
            path_signs.append(1.0)
        else:
            path_edges.append(edge_idx[(v, u)])
            path_signs.append(-1.0)

    # Choose 4 of the 8 path edges for degree-4 term
    for subset in combinations(range(8), 4):
        eidxs = tuple(sorted(path_edges[j] for j in subset))
        sign = 1.0
        for j in subset:
            sign *= path_signs[j]
        w4_coeffs[eidxs] += sign

    perm_count += 1
    if perm_count % 100000 == 0:
        elapsed = time.time() - t0
        print(f"  {perm_count}/362880 permutations ({elapsed:.1f}s)")

t1 = time.time()
print(f"Done in {t1-t0:.1f}s")

# Filter to nonzero coefficients
nonzero = {k: v for k, v in w4_coeffs.items() if abs(v) > 0.5}
print(f"\nTotal degree-4 monomials with nonzero coefficient: {len(nonzero)}")
print(f"Total degree-4 monomials: {m*(m-1)*(m-2)*(m-3)//24}")

# Classify nonzero monomials by support type
nz_by_type = defaultdict(list)
for mono_eidx, coeff in nonzero.items():
    elist = [edges[i] for i in mono_eidx]
    gtype = classify_support(elist)
    nz_by_type[gtype].append(coeff)

print(f"\nNonzero coefficients by support graph type:")
print(f"{'Type':>25} | {'Count':>5} | {'Coeff values':>20}")
print("-"*60)
for gtype in sorted(nz_by_type.keys()):
    coeffs = nz_by_type[gtype]
    vals = sorted(set(int(round(c)) for c in coeffs))
    print(f"  {gtype[0]}v {gtype[1]} | {len(coeffs):>5} | {vals}")

# Check: at n=7, the nonzero types were:
# P (5v spanning path): coefficient ±12 (from LEMMA P3: 12 directed paths)
# Q (6v P2 pair): coefficient ±24 (from LEMMA Q3: 24 directed paths)
# What are the new types at n=9?

print("\n" + "="*70)
print("DETAILED ANALYSIS BY SUPPORT TYPE")
print("="*70)

for gtype in sorted(nz_by_type.keys()):
    coeffs = nz_by_type[gtype]
    nv, deg_seq = gtype
    vals = sorted(set(int(round(c)) for c in coeffs))
    total_count = type_counts.get(gtype, 0)
    print(f"\n  Type: {nv} vertices, degree sequence {deg_seq}")
    print(f"    Total monomials of this type: {total_count}")
    print(f"    Nonzero in w_4: {len(coeffs)}")
    print(f"    Coefficient values: {vals}")

    # Show an example
    for mono_eidx, coeff in nonzero.items():
        elist = [edges[i] for i in mono_eidx]
        gt = classify_support(elist)
        if gt == gtype:
            print(f"    Example: edges {elist}, coeff = {int(round(coeff))}")
            break

# Final dimension analysis
print("\n" + "="*70)
print("DIMENSION ANALYSIS SUMMARY")
print("="*70)

# Count distinct coefficient absolute values per type
print("\nThe degree-4 Fourier space dimension equals the number of")
print("linearly independent 'types' with nonzero w_4 coefficients.\n")

# Actually, to properly determine the dimension, we should look at
# how many independent linear functions the nonzero monomials define.
# Group monomials by their coefficient value — if all monomials of a type
# have the same |coefficient|, then each type contributes one dimension.

# But types with the same coefficient pattern might still be linearly independent
# in the sense that they have different supports.

# The true dimension: since w_4 = 240*[deg-4 of t5] + 480*[deg-4 of a33],
# and [deg-4 of t5] and [deg-4 of a33] have disjoint support (by the n=7 analysis),
# the dimension is exactly 2 if both are nonzero.

# Let's verify disjoint support at n=9.
print("Checking support overlap between types at n=9:")

type_list = sorted(nz_by_type.keys())
for i, gt1 in enumerate(type_list):
    for gt2 in type_list[i+1:]:
        # Check if any monomial appears in both types
        # (impossible if types are defined by graph structure)
        print(f"  {gt1} vs {gt2}: disjoint (different graph types)")

print(f"\nConclusion: The degree-4 Fourier space at n=9 has dimension = {len(type_list)}")
print(f"  Types: {type_list}")
