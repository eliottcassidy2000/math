#!/usr/bin/env python3
"""beta2_arcflip_proof_attempt.py - Try to prove arc-flip invariance algebraically

KEY FACTS:
1. beta2 = dim(Z2) - rk(d3) is arc-flip invariant (verified n<=8)
2. D|Om2| = p_{vu} - p_{uv} where p_{xy} = #{allowed 2-paths through x?y}
3. dim(Z1) = (n-1)(n-2)/2 is CONSTANT
4. D(rk d2) = -D(beta1)
5. We need: D(rk d3) = D(dim Z2) = D(dim Om2) - D(rk d2) = D|Om2| + D(beta1)

So the claim is: D(rk d3) = (p_{vu} - p_{uv}) + D(beta1).

Let's try to understand EACH term:
- D(dim Om2): counts how many TT triples + NT cancellation dimensions change
- D(rk d2): how the rank of d2: Om2 ? Om1 changes
- D(rk d3): how the rank of d3: Om3 ? Om2 changes

For the proof, we need to show the d3 rank changes exactly track Z2 changes.

NEW APPROACH: Can we decompose the change into "local" parts that each satisfy
the invariance? I.e., D(rk d3) and D(dim Z2) both decompose into contributions
from specific path types, and these contributions match term by term.

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved

random.seed(42)


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def compute_all_dims(A, n):
    """Compute dim(Om_k), dim(Z_k), rk(d_k), beta_k for k=1,2,3."""
    paths0 = [(i,) for i in range(n)]
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)

    # Om1
    if not paths1:
        return {'dim_O1': 0, 'dim_O2': 0, 'dim_O3': 0,
                'rk_d1': 0, 'rk_d2': 0, 'rk_d3': 0,
                'dim_Z1': 0, 'dim_Z2': 0, 'b1': 0, 'b2': 0}

    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    dim_O1 = omega1.shape[1] if omega1.ndim == 2 else 0

    D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
    if dim_O1 > 0:
        D1_om = D1 @ omega1
        sv1 = np.linalg.svd(D1_om, compute_uv=False)
        rk_d1 = int(sum(s > 1e-8 for s in sv1))
    else:
        rk_d1 = 0

    # Om2
    if paths2:
        omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    else:
        dim_O2 = 0
        omega2 = None

    if dim_O2 > 0:
        D2 = build_full_boundary_matrix([tuple(p) for p in paths2],
                                        [tuple(p) for p in paths1])
        D2_om = D2 @ omega2
        sv2 = np.linalg.svd(D2_om, compute_uv=False)
        rk_d2 = int(sum(s > 1e-8 for s in sv2))
    else:
        rk_d2 = 0

    # Om3
    if paths3:
        omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
        dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0
    else:
        dim_O3 = 0
        omega3 = None

    if dim_O3 > 0 and dim_O2 > 0:
        D3 = build_full_boundary_matrix([tuple(p) for p in paths3],
                                        [tuple(p) for p in paths2])
        D3_om = D3 @ omega3
        # Project into Om2 basis
        # Actually rk_d3 is the rank of d3 restricted to Om3 ? Om2
        # The image lives in Om2. We need to project D3_om into Om2.
        # D3_om gives coordinates in the paths2 basis. We need to check
        # how much of D3_om lands in Om2.
        # Actually, the image of d3 on Om3 automatically lands in Om2
        # because d3 maps Om3 ? Om2 (GRIGORYAN's definition).
        # But we need to compute rank in Om2 coordinates.
        # Method: project each column of D3_om onto Om2 space.
        # Actually, im(d3) is a subspace of Z2 <= Om2.
        # The rank of d3|_Om3 as a map to the full path space is what we want.
        # But we need to ensure the image is in Om2.
        # Simpler: compute rank of the composite Om3 ? full_paths_2 ? Om2.
        # Since d3(Om3) <= Om2 automatically (chain complex property), the rank
        # is just the rank of D3_om projected to Om2 coordinates.

        # Project D3_om into Om2 coordinates
        # D3_om has shape (|paths2|, dim_O3)
        # omega2 has shape (|paths2|, dim_O2)
        # We want: omega2^+ @ D3_om
        proj = np.linalg.lstsq(omega2, D3_om, rcond=None)[0]
        sv3 = np.linalg.svd(proj, compute_uv=False)
        rk_d3 = int(sum(s > 1e-8 for s in sv3))
    else:
        rk_d3 = 0

    dim_Z1 = dim_O1 - rk_d1
    dim_Z2 = dim_O2 - rk_d2
    b1 = dim_Z1 - rk_d2  # Wait, b1 = dim(ker d1 in O1) - rk(d2 in O2 -> O1)
    # Actually b1 = dim_O1 - rk_d1 - rk_d2
    b1 = dim_O1 - rk_d1 - rk_d2
    b2 = dim_Z2 - rk_d3

    return {
        'dim_O1': dim_O1, 'dim_O2': dim_O2, 'dim_O3': dim_O3,
        'rk_d1': rk_d1, 'rk_d2': rk_d2, 'rk_d3': rk_d3,
        'dim_Z1': dim_Z1, 'dim_Z2': dim_Z2,
        'b1': b1, 'b2': b2
    }


def flip_arc(A, n, u, v):
    """Flip arc u?v to v?u. Returns new adjacency matrix."""
    B = [row[:] for row in A]
    B[u][v] = 0
    B[v][u] = 1
    return B


# ============================================================
# Part 1: Verify D|A3| = (n-3) * D|A2| at n=6,7
# ============================================================
print("=" * 70)
print("PATH COUNTING FORMULA: D|A3| = (n-3) * D|A2|")
print("=" * 70)

for n in [5, 6, 7]:
    random.seed(42)
    trials = {5: 200, 6: 100, 7: 50}[n]
    ratio_dist = Counter()
    mismatches = 0
    tested = 0

    for trial in range(trials):
        A = random_tournament(n)
        # Try one random arc flip
        edges = [(i, j) for i in range(n) for j in range(n) if A[i][j]]
        u, v = edges[random.randint(0, len(edges)-1)]

        # Count allowed paths
        p2_before = len(enumerate_allowed_paths(A, n, 2))
        p3_before = len(enumerate_allowed_paths(A, n, 3))

        B = flip_arc(A, n, u, v)
        p2_after = len(enumerate_allowed_paths(B, n, 2))
        p3_after = len(enumerate_allowed_paths(B, n, 3))

        d2 = p2_after - p2_before
        d3 = p3_after - p3_before

        tested += 1
        if d2 != 0:
            ratio = d3 / d2
            ratio_dist[ratio] += 1
            if abs(ratio - (n - 3)) > 0.01:
                mismatches += 1
        else:
            if d3 != 0:
                mismatches += 1
                ratio_dist['d2=0,d3!=0'] += 1
            else:
                ratio_dist['both=0'] += 1

    print(f"\nn={n}: {tested} arc flips tested")
    print(f"  Mismatches (D3 != (n-3)*D2): {mismatches}")
    print(f"  Ratio distribution:")
    for r, cnt in sorted(ratio_dist.items(), key=lambda x: str(x)):
        print(f"    {r}: {cnt}")
    print(f"  Expected ratio: n-3 = {n-3}")


# ============================================================
# Part 2: Verify full dimension formula under arc flip
# ============================================================
print(f"\n{'=' * 70}")
print("FULL DIMENSION ANALYSIS UNDER ARC FLIP")
print("=" * 70)

for n in [5, 6]:
    random.seed(42)
    trials = {5: 100, 6: 50}[n]
    formula_ok = 0
    formula_fail = 0

    for trial in range(trials):
        A = random_tournament(n)
        edges = [(i, j) for i in range(n) for j in range(n) if A[i][j]]
        u, v = edges[random.randint(0, len(edges)-1)]
        B = flip_arc(A, n, u, v)

        dims_A = compute_all_dims(A, n)
        dims_B = compute_all_dims(B, n)

        d_O2 = dims_B['dim_O2'] - dims_A['dim_O2']
        d_O3 = dims_B['dim_O3'] - dims_A['dim_O3']
        d_Z2 = dims_B['dim_Z2'] - dims_A['dim_Z2']
        d_rk_d2 = dims_B['rk_d2'] - dims_A['rk_d2']
        d_rk_d3 = dims_B['rk_d3'] - dims_A['rk_d3']
        d_b1 = dims_B['b1'] - dims_A['b1']
        d_b2 = dims_B['b2'] - dims_A['b2']

        # Check: D(rk d3) = D(dim Z2)  iff  Dbeta2 = 0
        if d_b2 == 0:
            formula_ok += 1
        else:
            formula_fail += 1
            print(f"  FAIL at trial {trial}: Dbeta2={d_b2}")
            print(f"    dims_A: {dims_A}")
            print(f"    dims_B: {dims_B}")

        # Check individual relationships
        # D(dim Z2) = D(dim Om2) - D(rk d2) = d_O2 - d_rk_d2
        check_Z2 = (d_Z2 == d_O2 - d_rk_d2)
        if not check_Z2:
            print(f"  Z2 identity fail: d_Z2={d_Z2}, d_O2-d_rk_d2={d_O2-d_rk_d2}")

    print(f"\nn={n}: Dbeta2=0 verified: {formula_ok}/{formula_ok+formula_fail}")


# ============================================================
# Part 3: Decompose D(dim Om2) into explicit terms
# ============================================================
print(f"\n{'=' * 70}")
print("DECOMPOSING D(dim Om2)")
print("=" * 70)

# When we flip u?v to v?u:
# A TT triple (a,b,c) exists when a?c (not just through b).
# Om2 includes TT triples + NT cancellation pairs.
#
# TT triples involving the flipped edge:
# - u?v as first edge: (u,v,c) where u?c. TT iff u?c, which is independent of flip.
#   Wait, (u,v,c): allowed iff u?v and v?c. u?c makes it TT.
#   After flip: v?u, so (u,v,c) is no longer allowed (u?v gone).
#   But (v,u,c) becomes allowed if v?u and u?c.
#
# Let me think more carefully using the ABCD partition.
# V \ {u,v} partitioned into:
#   A = {w : u?w, w?v}  (u beats, v loses to)
#   B = {w : w?u, v?w}  (u loses to, v beats)
#   C = {w : u?w, v?w}  (both beat)
#   D = {w : w?u, w?v}  (both lose to)

n = 5
random.seed(42)
abcd_data = Counter()

for trial in range(500):
    A = random_tournament(n)
    edges = [(i, j) for i in range(n) for j in range(n) if A[i][j]]
    u, v = edges[random.randint(0, len(edges)-1)]

    # Classify vertices
    a_set = [w for w in range(n) if w != u and w != v and A[u][w] and A[w][v]]
    b_set = [w for w in range(n) if w != u and w != v and A[w][u] and A[v][w]]
    c_set = [w for w in range(n) if w != u and w != v and A[u][w] and A[v][w]]
    d_set = [w for w in range(n) if w != u and w != v and A[w][u] and A[w][v]]

    a, b, c, d = len(a_set), len(b_set), len(c_set), len(d_set)

    # Compute dimensions
    B_mat = flip_arc(A, n, u, v)
    p2_before = len(enumerate_allowed_paths(A, n, 2))
    p2_after = len(enumerate_allowed_paths(B_mat, n, 2))

    dims_A_full = compute_all_dims(A, n)
    dims_B_full = compute_all_dims(B_mat, n)

    d_O2 = dims_B_full['dim_O2'] - dims_A_full['dim_O2']
    d_b1 = dims_B_full['b1'] - dims_A_full['b1']
    d_Z2 = dims_B_full['dim_Z2'] - dims_A_full['dim_Z2']
    d_rk_d3 = dims_B_full['rk_d3'] - dims_A_full['rk_d3']

    abcd_data[(a, b, c, d, d_O2, d_b1, d_Z2, d_rk_d3)] += 1

print(f"\nn=5: ABCD partition vs dimension changes")
print(f"{'(a,b,c,d)':>12} {'DOm2':>5} {'Dbeta1':>5} {'DZ2':>5} {'Drk3':>6} {'count':>6}")
for (a, b, c, d, dO2, db1, dZ2, dr3), cnt in sorted(abcd_data.items()):
    print(f"  ({a},{b},{c},{d})    {dO2:>5} {db1:>5} {dZ2:>5} {dr3:>6}   {cnt:>5}")


# ============================================================
# Part 4: Check formula D(dim Om2) = b - a (conjecture)
# ============================================================
print(f"\n{'=' * 70}")
print("TESTING: D(dim Om2) = f(a,b,c,d)?")
print("=" * 70)

# From the ABCD decomposition:
# Paths through u?v: all (w,u,v), (u,v,w), (w,u,v,w') etc.
# After flip to v?u: these disappear, but new paths through v?u appear.
#
# 2-paths through u?v before flip: (w,u,v) where w?u, and (u,v,w) where v?w
#   # of (w,u,v) = |{w: w?u}| = |B|+|D| (since w?u for w in B and D)
#   # of (u,v,w) = |{w: v?w}| = |B|+|C| (since v?w for w in B and C)
#   Total 2-paths using edge u?v: |B|+|D| + |B|+|C| = 2|B|+|C|+|D|
#
# 2-paths through v?u after flip: (w,v,u) where w?v, and (v,u,w) where u?w
#   # of (w,v,u) = |{w: w?v}| = |A|+|D| (since w?v for w in A and D)
#   # of (v,u,w) = |{w: u?w}| = |A|+|C| (since u?w for w in A and C)
#   Total 2-paths using edge v?u: |A|+|D| + |A|+|C| = 2|A|+|C|+|D|
#
# So D(# allowed 2-paths) = (2|A|+|C|+|D|) - (2|B|+|C|+|D|) = 2(|A|-|B|) = 2(a-b)
#
# Note: a+b+c+d = n-2

print("Formula check: D(# 2-paths) = 2(a-b)?")
for trial in range(200):
    A = random_tournament(n)
    edges = [(i, j) for i in range(n) for j in range(n) if A[i][j]]
    u, v = edges[random.randint(0, len(edges)-1)]

    a_set = [w for w in range(n) if w != u and w != v and A[u][w] and A[w][v]]
    b_set = [w for w in range(n) if w != u and w != v and A[w][u] and A[v][w]]
    a, b = len(a_set), len(b_set)

    B_mat = flip_arc(A, n, u, v)
    d_p2 = len(enumerate_allowed_paths(B_mat, n, 2)) - len(enumerate_allowed_paths(A, n, 2))

    if d_p2 != 2 * (a - b):
        print(f"  MISMATCH: Dp2={d_p2}, 2(a-b)={2*(a-b)}")

print("  (No output means formula is correct)")


# ============================================================
# Part 5: Explore D(dim Om2) vs D(# allowed 2-paths)
# ============================================================
print(f"\n{'=' * 70}")
print("D(dim Om2) vs D(# allowed 2-paths)")
print("=" * 70)

# dim(Om2) is NOT the same as # allowed 2-paths.
# Om2 = ker(d0) at degree 2, which is the space of "allowed" combinations.
# For TT paths, each TT triple contributes 1 to dim(Om2).
# For NT paths, cancellation pairs contribute.

# Actually, Om_p = intersection of ker of all face maps except d_0 and d_p.
# At p=2: Om2 = ker(d1) in the 2-chain space.
# d1(a,b,c) = (a,c) - this must be an allowed 1-path (edge a?c) for (a,b,c) to be in Om2.
# So Om2 = {(a,b,c) : a?b?c is allowed AND a?c is an edge (TT)}
#   ? cancellation space from NT paths.

# Let me just directly compare dim(Om2) and #paths2 to see the relationship.
for trial in range(100):
    A = random_tournament(n)
    paths2 = enumerate_allowed_paths(A, n, 2)
    dims = compute_all_dims(A, n)
    if dims['dim_O2'] != len(paths2):
        print(f"  dim(Om2)={dims['dim_O2']}, #paths2={len(paths2)}")
        break
else:
    print(f"  dim(Om2) = #paths2 for all n={n} tested")

# Check at n=6
n = 6
for trial in range(50):
    A = random_tournament(n)
    paths2 = enumerate_allowed_paths(A, n, 2)
    dims = compute_all_dims(A, n)
    if dims['dim_O2'] != len(paths2):
        print(f"  n=6: dim(Om2)={dims['dim_O2']}, #paths2={len(paths2)}")
        break
else:
    print(f"  dim(Om2) = #paths2 for all n={n} tested")


# ============================================================
# Part 6: So D(dim Om2) = 2(a-b) and we need
# D(rk d3) = D(dim Om2) - D(rk d2) = 2(a-b) + Dbeta1
# ============================================================
print(f"\n{'=' * 70}")
print("KEY IDENTITY: D(rk d3) = 2(a-b) + Dbeta1")
print("=" * 70)

n = 5
random.seed(42)
formula_ok = 0
formula_fail = 0

for trial in range(500):
    A = random_tournament(n)
    edges = [(i, j) for i in range(n) for j in range(n) if A[i][j]]
    u, v = edges[random.randint(0, len(edges)-1)]

    a_set = [w for w in range(n) if w != u and w != v and A[u][w] and A[w][v]]
    b_set = [w for w in range(n) if w != u and w != v and A[w][u] and A[v][w]]
    a, b = len(a_set), len(b_set)

    B_mat = flip_arc(A, n, u, v)
    dims_A = compute_all_dims(A, n)
    dims_B = compute_all_dims(B_mat, n)

    d_rk_d3 = dims_B['rk_d3'] - dims_A['rk_d3']
    d_b1 = dims_B['b1'] - dims_A['b1']
    expected = 2 * (a - b) + d_b1

    if d_rk_d3 == expected:
        formula_ok += 1
    else:
        formula_fail += 1
        print(f"  FAIL: Drk3={d_rk_d3}, 2(a-b)+Dbeta1={expected}, a={a}, b={b}, Dbeta1={d_b1}")

print(f"\nn=5: D(rk d3) = 2(a-b) + Dbeta1: {formula_ok}/{formula_ok+formula_fail}")

# Now at n=6
n = 6
random.seed(42)
formula_ok = 0
formula_fail = 0

for trial in range(200):
    A = random_tournament(n)
    edges = [(i, j) for i in range(n) for j in range(n) if A[i][j]]
    u, v = edges[random.randint(0, len(edges)-1)]

    a_set = [w for w in range(n) if w != u and w != v and A[u][w] and A[w][v]]
    b_set = [w for w in range(n) if w != u and w != v and A[w][u] and A[v][w]]
    a, b = len(a_set), len(b_set)

    B_mat = flip_arc(A, n, u, v)
    dims_A = compute_all_dims(A, n)
    dims_B = compute_all_dims(B_mat, n)

    d_rk_d3 = dims_B['rk_d3'] - dims_A['rk_d3']
    d_b1 = dims_B['b1'] - dims_A['b1']
    expected = 2 * (a - b) + d_b1

    if d_rk_d3 == expected:
        formula_ok += 1
    else:
        formula_fail += 1
        if formula_fail <= 3:
            print(f"  FAIL: Drk3={d_rk_d3}, 2(a-b)+Dbeta1={expected}, a={a}, b={b}, Dbeta1={d_b1}")

print(f"\nn=6: D(rk d3) = 2(a-b) + Dbeta1: {formula_ok}/{formula_ok+formula_fail}")

# At n=7
n = 7
random.seed(42)
formula_ok = 0
formula_fail = 0

for trial in range(100):
    A = random_tournament(n)
    edges = [(i, j) for i in range(n) for j in range(n) if A[i][j]]
    u, v = edges[random.randint(0, len(edges)-1)]

    a_set = [w for w in range(n) if w != u and w != v and A[u][w] and A[w][v]]
    b_set = [w for w in range(n) if w != u and w != v and A[w][u] and A[v][w]]
    a, b = len(a_set), len(b_set)

    B_mat = flip_arc(A, n, u, v)
    dims_A = compute_all_dims(A, n)
    dims_B = compute_all_dims(B_mat, n)

    d_rk_d3 = dims_B['rk_d3'] - dims_A['rk_d3']
    d_b1 = dims_B['b1'] - dims_A['b1']
    expected = 2 * (a - b) + d_b1

    if d_rk_d3 == expected:
        formula_ok += 1
    else:
        formula_fail += 1
        if formula_fail <= 3:
            print(f"  FAIL: Drk3={d_rk_d3}, 2(a-b)+Dbeta1={expected}")

print(f"\nn=7: D(rk d3) = 2(a-b) + Dbeta1: {formula_ok}/{formula_ok+formula_fail}")


# ============================================================
# Part 7: If the formula holds, can we prove it?
# ============================================================
print(f"\n{'=' * 70}")
print("PROOF ATTEMPT SUMMARY")
print("=" * 70)
print("""
If D(rk d3) = 2(a-b) + Dbeta1 for all tournaments and all arc flips,
then beta2 is arc-flip invariant, hence beta2(T) = beta2(T_trans) = 0 for all T.

The formula says:
  D(rk d3) = D(dim Om2) + Dbeta1
            = D(dim Om2) - D(rk d2)

This is equivalent to:
  D(rk d3) + D(rk d2) = D(dim Om2)

I.e., the TOTAL rank change in d2 and d3 exactly matches the dimension change of Om2.
This means: new dimensions in Om2 are either captured by d2 (become boundaries)
or by d3 (get filled by Om3), with no new cycles left uncovered.

EQUIVALENT: beta2 = dim(Om2) - rk(d2) - rk(d3) is invariant under arc flip.
""")


print("\n\nDone.")
