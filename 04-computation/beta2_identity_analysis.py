#!/usr/bin/env python3
"""
beta2_identity_analysis.py - Search for the algebraic identity proving beta_2=0

ESTABLISHED FORMULA: dim(Omega_2) = |A_2| - J_2
where J_2 = #{(a,c) : c->a and exists b with a->b->c}

GOAL: Prove rk(d_2) + rk(d_3) = |A_2| - J_2 = dim(Omega_2) [equiv to beta_2=0]

We know: rk(d_2) = dim(Z_1) - beta_1 = C(n-1,2) - beta_1

So the identity becomes:
  rk(d_3) = |A_2| - J_2 - C(n-1,2) + beta_1

Can we find combinatorial formulas for each quantity?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def compute_J2(A, n):
    """J_2 = #{(a,c): c->a and k_{a,c} >= 1}"""
    J2 = 0
    for a in range(n):
        for c in range(n):
            if a == c or not A[c][a]:
                continue
            for b in range(n):
                if b != a and b != c and A[a][b] and A[b][c]:
                    J2 += 1
                    break
    return J2

def compute_J2_total(A, n):
    """Total k values: sum_{c->a} k_{a,c}"""
    total_k = 0
    for a in range(n):
        for c in range(n):
            if a == c or not A[c][a]:
                continue
            k = sum(1 for b in range(n) if b != a and b != c and A[a][b] and A[b][c])
            total_k += k
    return total_k

def compute_c3(A, n):
    """Number of directed 3-cycles (as vertex sets)."""
    c3 = 0
    for S in combinations(range(n), 3):
        a, b, c = S
        # Check all 2 cyclic orderings
        if (A[a][b] and A[b][c] and A[c][a]):
            c3 += 1
        elif (A[a][c] and A[c][b] and A[b][a]):
            c3 += 1
    return c3

def compute_full_data(A, n, max_p=4):
    """Compute full chain complex data."""
    paths = {}
    omega = {}
    for p in range(max_p+1):
        paths[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega[p] = np.eye(n)
        elif len(paths[p]) > 0 and len(paths[p-1]) > 0:
            omega[p] = compute_omega_basis(A, n, p, paths[p], paths[p-1])
        else:
            omega[p] = np.zeros((max(1, len(paths[p])), 0))

    dims = {p: (omega[p].shape[1] if omega[p].ndim == 2 else 0) for p in range(max_p+1)}
    rks = {}

    for p in range(1, max_p+1):
        if dims[p] > 0 and len(paths.get(p-1, [])) > 0:
            bd = build_full_boundary_matrix(paths[p], paths[p-1])
            bd_om = bd @ omega[p]
            if p >= 2 and dims[p-1] > 0:
                coords, _, _, _ = np.linalg.lstsq(omega[p-1], bd_om, rcond=None)
                rks[p] = np.linalg.matrix_rank(coords, tol=1e-8)
            else:
                rks[p] = np.linalg.matrix_rank(bd_om, tol=1e-8)
        else:
            rks[p] = 0

    betas = {}
    for p in range(max_p):
        z_p = dims[p] - rks.get(p, 0)
        b_p = z_p - rks.get(p+1, 0)
        betas[p] = b_p

    return paths, omega, dims, rks, betas


print("=" * 70)
print("IDENTITY ANALYSIS: rk(d_3) = |A_2| - J_2 - C(n-1,2) + beta_1")
print("=" * 70)

for n in [4, 5]:
    n_arcs = n*(n-1)//2
    total = 1 << n_arcs
    Z1 = (n-1)*(n-2)//2  # constant

    print(f"\nn = {n}, dim(Z_1) = {Z1}")

    errors = 0
    data_table = []

    for bits in range(total):
        A = build_adj(n, bits)
        paths, omega, dims, rks, betas = compute_full_data(A, n)

        A2 = len(paths[2])
        J2 = compute_J2(A, n)
        c3 = compute_c3(A, n)
        beta1 = betas[1]
        rk_d3 = rks[3]

        # Predicted rk_d3
        pred = A2 - J2 - Z1 + beta1
        if pred != rk_d3:
            errors += 1

        data_table.append({
            'bits': bits, 'A2': A2, 'J2': J2, 'c3': c3,
            'beta1': beta1, 'rk_d3': rk_d3, 'dim_Om2': dims[2],
            'dim_Om3': dims[3], 'rk_d2': rks[2]
        })

    print(f"  Identity verification errors: {errors}")

    # Explore correlations
    # Is J_2 related to c_3 (3-cycles)?
    j2_vs_c3 = defaultdict(set)
    for d in data_table:
        j2_vs_c3[d['c3']].add(d['J2'])

    print(f"\n  J_2 vs c_3:")
    for c3 in sorted(j2_vs_c3.keys()):
        vals = sorted(j2_vs_c3[c3])
        if len(vals) <= 5:
            print(f"    c3={c3}: J2 in {vals}")
        else:
            print(f"    c3={c3}: J2 in [{min(vals)}, {max(vals)}], {len(vals)} values")


# ANALYSIS 2: What ARE the J_2 pairs? How do they relate to directed cycles?
print(f"\n{'='*70}")
print("ANALYSIS 2: J_2 structure in terms of neighborhoods")
print("=" * 70)

n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

# For each tournament, compute J_2 using neighborhood structure
print("Computing J_2 via neighborhoods...")

# J_2 counts (a,c) with c->a and exists b: a->b, b->c
# This is: c->a and (out(a) cap in(c)) \ {a,c} is nonempty
# = c->a and a has a common out-neighbor with c's in-neighborhood
# = the number of "backed" pairs with at least one 2-step connection

# Let's express J_2 in terms of the adjacency matrix
# J_2 = #{(a,c) : A[c][a]=1 and (A^2)[a][c] > 0}
# where A^2[a][c] = sum_b A[a][b]*A[b][c] (excluding b=a,c)

print("Verifying J_2 = #{(a,c): A[c][a]=1 and A^2[a][c] > 0}...")
errors = 0
for bits in range(total):
    A = build_adj(n, bits)
    J2 = compute_J2(A, n)

    # Compute A^2 (excluding self and direct)
    A_mat = np.array(A, dtype=float)
    A2_mat = A_mat @ A_mat  # A^2[a][c] = sum_b A[a][b]*A[b][c]

    J2_mat = 0
    for a in range(n):
        for c in range(n):
            if a != c and A[c][a] and A2_mat[a][c] > 0.5:
                J2_mat += 1

    if J2 != J2_mat:
        errors += 1

print(f"  Errors: {errors}")
if errors == 0:
    print("  CONFIRMED: J_2 = #{(a,c): A[c][a]=1 and (A^2)[a][c] > 0}")


# ANALYSIS 3: Total NT paths = sum_{c->a} k_{a,c}
print(f"\n{'='*70}")
print("ANALYSIS 3: Total NT paths and relationship to |A_2|")
print("=" * 70)

# NT = |A_2| - TT. We know TT = #{(a,b,c): a->b->c, a->c}
# And NT = |A_2| - TT = #{(a,b,c): a->b->c, c->a}

# NT = sum_{c->a} k_{a,c}
# TT = sum_{a->c} k_{a,c} (where here k_{a,c} for a->c counts intermediates)

# Interesting: |A_2| = TT + NT = sum_{all ordered (a,c)} k_{a,c}
# And A^2[a][c] counts the number of a->b->c paths (for ANY b != a,c? or including?)
# Actually A^2[a][c] = sum_b A[a][b]*A[b][c] including b=a,c.
# But A[a][a] = 0 and A[c][c] = 0 (tournament), so b=a and b=c contribute 0.
# So A^2[a][c] = k_{a,c} for any a != c. Great!

# Therefore: TT = sum_{a->c} A^2[a][c] = sum_{a,c: A[a][c]=1} A^2[a][c]
# And NT = sum_{c->a} A^2[a][c] = sum_{a,c: A[c][a]=1} A^2[a][c]
# And |A_2| = sum_{a!=c} A^2[a][c] = tr(A^3)... wait no.
# sum_{a!=c} A^2[a][c] = sum of all entries of A^2 minus diagonal = sum(A^2) - tr(A^2)
# But tr(A^2) = sum_a sum_b A[a][b]*A[b][a] = 0 for tournaments (A[a][b]+A[b][a]=1 for a!=b,
# but A[a][b]*A[b][a] = 0 since exactly one is 1)
# So |A_2| = sum of all entries of A^2 = total(A^2)

# And J_2 = #{(a,c): A[c][a]=1 and A^2[a][c] > 0}
# This is the number of "reverse pairs with positive 2-step connection"

# In matrix terms: J_2 = |{(a,c): A^T[a][c] = 1 and A^2[a][c] > 0}|
# = number of entries where BOTH A^T and sgn(A^2) are 1

print("Verifying |A_2| = total(A^2)...")
errors = 0
for bits in range(total):
    A = build_adj(n, bits)
    A_mat = np.array(A, dtype=float)
    A2_mat = A_mat @ A_mat
    total_A2 = int(np.sum(A2_mat) + 0.5) - int(np.trace(A2_mat) + 0.5)
    A2_paths = len(enumerate_allowed_paths(A, n, 2))
    if total_A2 != A2_paths:
        errors += 1
print(f"  Errors: {errors}")

# Check if there's a simple relationship: J_2 = f(A, A^2)?
# J_2 = #{(a,c): A[c][a] * [A^2[a][c] > 0]}
# In terms of entry-wise: J_2 = sum_{a!=c} A^T[a][c] * sign(A^2[a][c])

# Is there a formula for J_2 in terms of tr(A^k) or similar spectral quantities?
# J_2 involves the SUPPORT of A^2, not just its values. This makes it non-polynomial.

# But maybe: J_2 = #{reverse pairs with 2-step path}
# = (total reverse pairs) - #{reverse pairs WITHOUT any 2-step path}
# Total reverse pairs: C(n,2) (each pair has one reverse direction)
# Wait: total ordered reverse pairs = n(n-1)/2 (the arcs c->a are just arcs in A^T,
# which is the same as A but with reversed orientation; there are n(n-1)/2 such pairs)

# #{(a,c): c->a, A^2[a][c]=0}: these are pairs where c->a but there's NO 2-step path a->b->c.
# This means: for every b != a,c, either a -/-> b or b -/-> c.
# I.e., out(a)\{c} and in(c)\{a} are DISJOINT sets.

# Since out(a)\{c} and in(c)\{a} are subsets of V\{a,c} (size n-2),
# and |out(a)\{c}| = d_a - [a->c] = d_a (since c->a, not a->c)... wait.
# d_a = out-degree of a. out(a) includes all vertices a points to.
# If c->a, then a does NOT point to c. So out(a)\{c} = out(a) (c not in it).
# |out(a)| = d_a.
# in(c)\{a}: all vertices pointing to c, except a.
# Since c->a, a is NOT in in(c) (because c->a means A[c][a]=1, not A[a][c]=1;
# wait, in(c) = {w: w->c, w!=c} = {w: A[w][c]=1}. Does a point to c? c->a means A[c][a]=1.
# But does a point to c? A[a][c] = 1 - A[c][a] = 0 (tournament). So a does NOT point to c.
# So a is NOT in in(c). in(c)\{a} = in(c) (since a not in it).
# |in(c)| = n-1 - d_c.

# For A^2[a][c] = 0 (disjoint): we need out(a) cap in(c) = empty (in V\{a,c}).
# out(a) subset V\{a}, |out(a)| = d_a.
# in(c) subset V\{c}, |in(c)| = n-1-d_c.
# Both are subsets of V\{a,c} (size n-2) if we exclude {a,c}.
# Actually: out(a) can include c or not. Since c->a, a -/-> c, so c NOT in out(a).
# And in(c) can include a or not. Since a -/-> c, a NOT in in(c).
# So both are subsets of V\{a,c}.

# Disjoint condition: |out(a)| + |in(c)| <= n-2
# i.e., d_a + (n-1-d_c) <= n-2
# i.e., d_a - d_c <= -1
# i.e., d_a < d_c

# Wait, this is a NECESSARY condition for A^2[a][c] = 0, but is it sufficient?
# Not quite: even if d_a + (n-1-d_c) <= n-2, the sets might overlap.
# d_a + (n-1-d_c) = d_a - d_c + n-1. For this to be <= n-2: d_a <= d_c - 1.
# If d_a <= d_c - 1, then |out(a)| + |in(c)| = d_a + (n-1-d_c) <= n-2.
# By pigeonhole, they CAN be disjoint, but don't HAVE to be.

# If d_a >= d_c, then |out(a)| + |in(c)| > n-2, so they MUST overlap.
# This means A^2[a][c] > 0 whenever d_a >= d_c (and c->a).
# So: if c->a and d_a >= d_c, then (a,c) is in J_2.

# Conversely, if d_a < d_c, the sets MIGHT or might not overlap.

print(f"\n{'='*70}")
print("ANALYSIS 4: When is A^2[a][c] = 0 for reverse pair c->a?")
print("=" * 70)

zero_by_scores = Counter()
nonzero_by_scores = Counter()

for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(row) for row in A]
    A_mat = np.array(A, dtype=float)
    A2_mat = A_mat @ A_mat

    for a in range(n):
        for c in range(n):
            if a == c or not A[c][a]:
                continue
            if A2_mat[a][c] < 0.5:
                zero_by_scores[(scores[a], scores[c])] += 1
            else:
                nonzero_by_scores[(scores[a], scores[c])] += 1

print("Reverse pairs (a,c) with A^2[a][c] = 0, by (d_a, d_c):")
all_score_pairs = set(zero_by_scores.keys()) | set(nonzero_by_scores.keys())
for sp in sorted(all_score_pairs):
    da, dc = sp
    z = zero_by_scores.get(sp, 0)
    nz = nonzero_by_scores.get(sp, 0)
    frac = z / (z + nz) if z + nz > 0 else 0
    print(f"  (d_a={da}, d_c={dc}): zero={z}, nonzero={nz}, frac_zero={frac:.3f}")

print("\nNote: d_a >= d_c should give frac_zero=0")
for sp in sorted(all_score_pairs):
    da, dc = sp
    z = zero_by_scores.get(sp, 0)
    if da >= dc and z > 0:
        print(f"  VIOLATION: d_a={da} >= d_c={dc} but {z} zeros found!")


# ANALYSIS 5: Formula J_2 = C(n,2) - #{reverse pairs with d_a < d_c and disjoint neighborhoods}
print(f"\n{'='*70}")
print("ANALYSIS 5: J_2 = C(n,2) - #{(a,c): c->a, A^2[a][c]=0}")
print("=" * 70)

errors_j2 = 0
for bits in range(total):
    A = build_adj(n, bits)
    j2 = compute_J2(A, n)
    A_mat = np.array(A, dtype=float)
    A2_mat = A_mat @ A_mat

    # Count reverse pairs with A^2=0
    zero_count = 0
    for a in range(n):
        for c in range(n):
            if a != c and A[c][a] and A2_mat[a][c] < 0.5:
                zero_count += 1

    # n(n-1)/2 total reverse pairs
    pred_j2 = n*(n-1)//2 - zero_count
    if j2 != pred_j2:
        errors_j2 += 1

print(f"  J_2 = C(n,2) - #{'{'}A^2[a][c]=0, c->a{'}'}: {errors_j2} errors")


# ANALYSIS 6: What IS the number of (a,c) with c->a and A^2[a][c]=0?
# These are pairs where c dominates a AND there's no 2-step path from a to c
# => out(a) and in(c) are disjoint (both subsets of V\{a,c})
print(f"\n{'='*70}")
print("ANALYSIS 6: Disconnected reverse pairs")
print("=" * 70)

# For tournament T with adjacency A:
# A reverse pair (a,c) with c->a is "disconnected" if A^2[a][c] = 0
# This means: every vertex w != a,c either (w->a) or (c->w) or both.
# In other words: no w has BOTH a->w AND w->c.
# = out(a) cap in(c) = empty (in V\{a,c})
# Since |out(a)| = d_a and |in(c)| = n-1-d_c, and both subsets of V\{a,c} (size n-2):
# disjoint requires d_a + (n-1-d_c) <= n-2, i.e., d_a <= d_c - 1

# Let's define: D(T) = #{(a,c): c->a, d_a <= d_c - 1, out(a) cap in(c) cap V\{a,c} = empty}
# Then J_2 = C(n,2) - D(T)

# Q: Is D(T) determined by score sequence?
D_by_scores = defaultdict(set)
for bits in range(total):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    A_mat = np.array(A, dtype=float)
    A2_mat = A_mat @ A_mat
    D = sum(1 for a in range(n) for c in range(n)
            if a != c and A[c][a] and A2_mat[a][c] < 0.5)
    D_by_scores[scores].add(D)

print("D(T) by score sequence:")
for scores in sorted(D_by_scores.keys()):
    vals = sorted(D_by_scores[scores])
    det = "DETERMINED" if len(vals) == 1 else "NOT determined"
    if len(vals) <= 5:
        print(f"  {scores}: D in {vals} -- {det}")
    else:
        print(f"  {scores}: D in [{min(vals)},{max(vals)}], {len(vals)} values -- {det}")


# ANALYSIS 7: Key relationship beta_1 and J_2
print(f"\n{'='*70}")
print("ANALYSIS 7: Relationship between J_2 and beta_1")
print("=" * 70)

j2_b1 = defaultdict(set)
b1_j2 = defaultdict(set)
for bits in range(total):
    A = build_adj(n, bits)
    j2 = compute_J2(A, n)
    paths, omega, dims, rks, betas = compute_full_data(A, n)
    b1 = betas[1]
    j2_b1[j2].add(b1)
    b1_j2[b1].add(j2)

print("J_2 -> beta_1 values:")
for j2 in sorted(j2_b1.keys()):
    print(f"  J_2={j2}: beta_1 in {sorted(j2_b1[j2])}")

print("\nbeta_1 -> J_2 values:")
for b1 in sorted(b1_j2.keys()):
    print(f"  beta_1={b1}: J_2 in {sorted(b1_j2[b1])}")


# ANALYSIS 8: The corrected identity at each tournament
print(f"\n{'='*70}")
print("ANALYSIS 8: rk(d_3) + J_2 + C(n-1,2) - beta_1 = |A_2| (should hold)")
print("=" * 70)

err = 0
for bits in range(total):
    A = build_adj(n, bits)
    j2 = compute_J2(A, n)
    paths, omega, dims, rks, betas = compute_full_data(A, n)
    b1 = betas[1]
    rk3 = rks[3]
    A2 = len(paths[2])
    Z1 = (n-1)*(n-2)//2

    lhs = rk3 + j2 + Z1 - b1
    if lhs != A2:
        err += 1
        if err <= 3:
            print(f"  FAIL: bits={bits}, rk3={rk3}, J2={j2}, Z1={Z1}, b1={b1}, A2={A2}, lhs={lhs}")

print(f"  Errors: {err}")

# ANALYSIS 9: Compute some ALTERNATIVE tournament invariants and correlate with J_2
print(f"\n{'='*70}")
print("ANALYSIS 9: J_2 in terms of classical tournament invariants")
print("=" * 70)

# Classical invariants: score sequence, c_3 (3-cycles), c_5 (5-cycles)
# Also: sum of d_i^2, sum d_i^3, etc.

# J_2 = C(n,2) - D where D = #{disconnected reverse pairs}
# = C(n,2) - #{(a,c): c->a, out(a) cap in(c) cap V\{a,c} = empty}

# Alternative: J_2 = #{(a,c): c->a, A^2[a][c] > 0}
# = sum_{a!=c} A^T[a][c] * (1 if A^2[a][c]>0 else 0)
# = C(n,2) - sum_{a!=c} A^T[a][c] * (1 if A^2[a][c]=0 else 0)

# Can we express this using tr(A^k)?
# sum_{a!=c} A^T[a][c] * A^2[a][c] = sum_{a!=c} A[c][a] * sum_b A[a][b]*A[b][c]
# = sum_{a,b,c distinct} A[c][a]*A[a][b]*A[b][c]
# = #{directed 3-cycles (c,a,b,c)}... sort of.
# Actually this counts ordered triples (a,b,c) where c->a, a->b, b->c, all distinct.
# This is 3 * c_3 (each 3-cycle counted 3 times for 3 starting vertices).

# Hmm wait: (c,a,b) with c->a->b->c is a 3-cycle starting at c.
# Differently: sum A^T * A^2 entrywise (excluding diagonal)
# = sum_{a,c: a!=c} A[c][a] * (A^2)[a][c]
# = sum_{a,c: a!=c} sum_b A[c][a]*A[a][b]*A[b][c]
# This overcounts because b can equal a or c... but A[a][a]=0, so b!=a guaranteed.
# And A[b][c]*A[c][a]*A[a][b] requires b!=c (otherwise A[c][c]=0).
# So this is sum over distinct (a,b,c) of A[c][a]*A[a][b]*A[b][c] = 3*c_3? No, wait.
# Each 3-cycle {a,b,c} with a->b->c->a contributes:
#   (a,b,c): A[c][a]*A[a][b]*A[b][c] = 1*1*1 = 1
# But the cycle a->b->c->a also appears as starting from b or c:
#   (b,c,a): A[a][b]*A[b][c]*A[c][a] = 1
#   (c,a,b): A[b][c]*A[c][a]*A[a][b] = 1
# So the total sum = 3 * c_3.

# But J_2 involves the INDICATOR A^2[a][c] > 0, not the VALUE A^2[a][c].
# These are different unless all nonzero A^2[a][c] equal 1.

# At n=5: A^2[a][c] can be 0, 1, 2, or 3.
# So J_2 != "something involving 3*c_3" in general.

# Let's check: sum_{c->a} A^2[a][c] vs J_2 vs c_3
inv_data = []
for bits in range(total):
    A = build_adj(n, bits)
    j2 = compute_J2(A, n)
    c3 = compute_c3(A, n)
    A_mat = np.array(A, dtype=float)
    A2_mat = A_mat @ A_mat

    NT = sum(int(A2_mat[a][c]) for a in range(n) for c in range(n)
             if a != c and A[c][a])
    sum_A2_rev = NT  # sum of A^2 over reverse pairs

    inv_data.append({'j2': j2, 'c3': c3, 'NT': NT, 'sum_A2_rev': sum_A2_rev})

    # Verify NT = 3*c3
    # NO! NT = total number of non-TT paths, not 3*c3.
    # 3*c3 = sum_{c->a} A^2[a][c] which is the SAME as NT? Let me check.

# Actually, sum_{c->a} A^2[a][c] = number of triples (a,b,c) with a->b, b->c, c->a
# = number of directed 3-cycles = 3*c_3? Not quite.
# A directed 3-cycle on vertex set {x,y,z} with x->y->z->x has 3 ordered triples:
# (x,y,z), (y,z,x), (z,x,y). Each gives a "reverse pair" contribution.
# For triple (x,y,z): a=x, b=y, c=z. c->a iff z->x (TRUE in cycle). A^2[x][z] includes y.
# So yes, each 3-cycle contributes 3 to the sum.
# But the sum also counts longer paths: (a,b,c) with a->b->c->a where {a,b,c} don't form
# a 3-cycle at the vertex-set level... wait, they do. If a->b, b->c, c->a, and a,b,c distinct,
# then {a,b,c} forms a directed 3-cycle. So sum_{c->a} A^2[a][c] = 3*c_3.

err_nt = 0
for d in inv_data:
    if d['NT'] != 3 * d['c3']:
        err_nt += 1
print(f"  NT = 3*c_3: {err_nt} errors")


# ANALYSIS 10: beta_1 in terms of c_3 or J_2?
print(f"\n{'='*70}")
print("ANALYSIS 10: beta_1 formula")
print("=" * 70)

# beta_1 = dim(Z_1) - rk(d_2) = C(n-1,2) - rk(d_2)
# rk(d_2) = rank of boundary map from Omega_2 to Omega_1
# Since Omega_1 = R_1 (all arcs), rk(d_2) = rank of d_2 restricted to Omega_2

# Dimension of Omega_2 = |A_2| - J_2 = TT + NT - J_2
# TT + NT = |A_2|
# J_2 = #{reverse pairs with at least one 2-step path}

# Can we express rk(d_2) combinatorially?
# d_2 maps Omega_2 -> R_1 (space of arcs, dim = n(n-1)/2)
# rk(d_2) is the rank of a certain matrix.

# At degree 1: Z_1 = ker(d_1) in Omega_1 = R_1
# d_1: R_1 -> R_0 (vertices): d_1(a,b) = b - a
# rk(d_1) = n-1 (always for tournaments, since they're strongly connected...
#            wait, some tournaments are NOT strongly connected)

# For a tournament that's NOT strongly connected: beta_0 > 1.
# For strongly connected tournaments: beta_0 = 1, rk(d_1) = n-1.
# dim(Z_1) = dim(R_1) - rk(d_1) = n(n-1)/2 - (n - beta_0)

# Hmm, actually for tournaments, rk(d_1) = n - beta_0. And the
# tournament is always strongly connected iff it has a Hamiltonian path (which all do).
# Wait, a tournament always has a Hamiltonian path, but might not be strongly connected.
# The transitive tournament on 3 vertices: 0->1->2, 0->2. Is this strongly connected?
# 2 cannot reach 0. So NO, the transitive tournament is NOT strongly connected.
# beta_0 of transitive tournament on 3 = 1 (there's always one connected component via
# the path 0->1->2... wait, path homology beta_0 counts connected components differently)

# Actually in path homology, beta_0 is always 1 for tournaments because every tournament
# has a Hamiltonian path, making it "path-connected" in the GLMY sense.

# Let me verify:
betas_dist = defaultdict(Counter)
for bits in range(total):
    A = build_adj(n, bits)
    _, _, _, _, betas = compute_full_data(A, n)
    for p in range(4):
        betas_dist[p][betas[p]] += 1

print("Betti number distributions at n=5:")
for p in sorted(betas_dist.keys()):
    print(f"  beta_{p}: {dict(sorted(betas_dist[p].items()))}")

# beta_1: what determines it?
b1_by_c3 = defaultdict(set)
for bits in range(total):
    A = build_adj(n, bits)
    c3 = compute_c3(A, n)
    _, _, _, _, betas = compute_full_data(A, n)
    b1_by_c3[c3].add(betas[1])

print("\nbeta_1 by c_3:")
for c3 in sorted(b1_by_c3.keys()):
    print(f"  c_3={c3}: beta_1 in {sorted(b1_by_c3[c3])}")

# c_3 ranges: c_3 = C(n,3) - sum C(d_i,2) / (n-2)... actually
# c_3 = C(n,3) - sum C(d_i, 2) for the score sequence.
# Wait: c_3 for tournaments is given by: c_3 = C(n,3) - sum C(d_i, 2)

err_c3 = 0
for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(row) for row in A]
    c3 = compute_c3(A, n)
    formula = n*(n-1)*(n-2)//6 - sum(s*(s-1)//2 for s in scores)
    if c3 != formula:
        err_c3 += 1
print(f"\nc_3 = C(n,3) - sum C(d_i,2): {err_c3} errors")


print("\nDone.")
