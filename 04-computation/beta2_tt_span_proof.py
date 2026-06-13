#!/usr/bin/env python3
"""
beta2_tt_span_proof.py - Does Omega_2 = TT-span for ALL tournaments?

opus-S45 found that at n=5, Omega_2 is spanned entirely by TT (transitive
triple) paths. This is a powerful structural property.

If Omega_2 = TT-span for all n, then:
1. Every 2-cycle is a linear combination of TT paths only
2. The boundary d_2 of such cycles is "clean" (no junk faces)
3. This simplifies the beta_2 = 0 analysis enormously

Key: A TT path (a,b,c) has a->b, b->c, a->c (transitive).
     An NT path (a,b,c) has a->b, b->c, c->a (non-transitive, part of 3-cycle).

We test:
1. Does Omega_2 = TT-span at n=6? (exhaustive)
2. Does it hold at n=7,8? (sampled)
3. If so, what does this mean for dim(Omega_2)?
4. Can we prove it algebraically?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, random, time
import numpy as np
from collections import Counter
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


def check_tt_span(A, n):
    """Check if Omega_2 is spanned entirely by TT paths.
    Returns (is_tt_span, dim_Om2, n_TT, n_NT, max_nt_coeff).
    """
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap1 = enumerate_allowed_paths(A, n, 1)
    if not ap2:
        return True, 0, 0, 0, 0.0

    om2 = compute_omega_basis(A, n, 2, ap2, ap1)
    d2 = om2.shape[1] if om2.ndim == 2 else 0
    if d2 == 0:
        return True, 0, 0, 0, 0.0

    # Identify TT and NT paths
    tt_idx = [i for i, (a,b,c) in enumerate(ap2) if A[a][c]]
    nt_idx = [i for i, (a,b,c) in enumerate(ap2) if not A[a][c]]
    n_tt = len(tt_idx)
    n_nt = len(nt_idx)

    # Check if any Omega_2 basis vector has nonzero NT component
    max_nt = 0.0
    for col in range(d2):
        for i in nt_idx:
            v = abs(om2[i, col])
            if v > max_nt:
                max_nt = v

    is_tt = max_nt < 1e-8
    return is_tt, d2, n_tt, n_nt, max_nt


# ============================================================
# TEST 1: n=5 exhaustive (re-verify opus finding)
# ============================================================
print("=" * 70)
print("TEST 1: Omega_2 = TT-span at n=5 (exhaustive)")
print("=" * 70)

n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

violations_5 = 0
tt_dim_dist = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    is_tt, d2, n_tt, n_nt, max_nt = check_tt_span(A, n)
    if not is_tt:
        violations_5 += 1
    tt_dim_dist[(n_tt, d2)] += 1

print(f"Violations: {violations_5}/{total}")
if violations_5 == 0:
    print("  CONFIRMED: Omega_2 = TT-span for ALL n=5 tournaments")
print(f"\n(n_TT, dim Omega_2) distribution:")
for key in sorted(tt_dim_dist.keys()):
    eq = "=" if key[0] == key[1] else "!="
    print(f"  TT={key[0]}, Omega_2={key[1]} ({eq}): {tt_dim_dist[key]}")


# ============================================================
# TEST 2: n=6 exhaustive
# ============================================================
print(f"\n{'='*70}")
print("TEST 2: Omega_2 = TT-span at n=6 (exhaustive)")
print("=" * 70)

n = 6
n_arcs = n*(n-1)//2
total = 1 << n_arcs

violations_6 = 0
tt_dim_dist_6 = Counter()
max_nt_overall = 0.0
violation_examples = []

t0 = time.time()
for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        dt = time.time() - t0
        print(f"  ... {bits}/{total} ({dt:.0f}s), violations={violations_6}")

    A = build_adj(n, bits)
    is_tt, d2, n_tt, n_nt, max_nt = check_tt_span(A, n)
    if not is_tt:
        violations_6 += 1
        if max_nt > max_nt_overall:
            max_nt_overall = max_nt
        if len(violation_examples) < 5:
            scores = sorted([sum(row) for row in A])
            violation_examples.append((bits, scores, d2, n_tt, n_nt, max_nt))
    tt_dim_dist_6[(n_tt, d2)] += 1

dt = time.time() - t0
print(f"\nCompleted in {dt:.0f}s")
print(f"Violations: {violations_6}/{total}")

if violations_6 == 0:
    print("  CONFIRMED: Omega_2 = TT-span for ALL n=6 tournaments!")
else:
    print(f"  FAILED: {violations_6} violations, max NT coeff = {max_nt_overall:.6f}")
    for bits, scores, d2, n_tt, n_nt, max_nt in violation_examples:
        print(f"    bits={bits}, scores={scores}, Om2={d2}, TT={n_tt}, NT={n_nt}, max_nt={max_nt:.6f}")

# Print TT vs Om2 relationship
print(f"\nDistinct (TT, dim Omega_2) patterns: {len(tt_dim_dist_6)}")
eq_count = sum(v for k, v in tt_dim_dist_6.items() if k[0] == k[1])
neq_count = sum(v for k, v in tt_dim_dist_6.items() if k[0] != k[1])
print(f"  TT = dim Omega_2: {eq_count} ({100*eq_count/total:.1f}%)")
print(f"  TT != dim Omega_2: {neq_count} ({100*neq_count/total:.1f}%)")

for key in sorted(tt_dim_dist_6.keys()):
    eq = "=" if key[0] == key[1] else "!="
    print(f"  TT={key[0]}, Omega_2={key[1]} ({eq}): {tt_dim_dist_6[key]}")


# ============================================================
# TEST 3: n=7 sampled
# ============================================================
print(f"\n{'='*70}")
print("TEST 3: Omega_2 = TT-span at n=7 (sampled)")
print("=" * 70)

random.seed(42)
n = 7
samples = 5000
violations_7 = 0
max_nt_7 = 0.0

t0 = time.time()
for s in range(samples):
    if s % 1000 == 0 and s > 0:
        dt = time.time() - t0
        print(f"  {s}/{samples} ({dt:.0f}s), violations={violations_7}")
    bits = random.getrandbits(n*(n-1)//2)
    A = build_adj(n, bits)
    is_tt, d2, n_tt, n_nt, max_nt = check_tt_span(A, n)
    if not is_tt:
        violations_7 += 1
        if max_nt > max_nt_7:
            max_nt_7 = max_nt

dt = time.time() - t0
print(f"Done in {dt:.0f}s. Violations: {violations_7}/{samples}")
if violations_7 == 0:
    print("  CONFIRMED (sampled): Omega_2 = TT-span at n=7")


# ============================================================
# TEST 4: n=8 sampled
# ============================================================
print(f"\n{'='*70}")
print("TEST 4: Omega_2 = TT-span at n=8 (sampled)")
print("=" * 70)

n = 8
samples = 2000
violations_8 = 0

t0 = time.time()
for s in range(samples):
    if s % 500 == 0 and s > 0:
        dt = time.time() - t0
        print(f"  {s}/{samples} ({dt:.0f}s), violations={violations_8}")
    bits = random.getrandbits(n*(n-1)//2)
    A = build_adj(n, bits)
    is_tt, d2, n_tt, n_nt, max_nt = check_tt_span(A, n)
    if not is_tt:
        violations_8 += 1

dt = time.time() - t0
print(f"Done in {dt:.0f}s. Violations: {violations_8}/{samples}")
if violations_8 == 0:
    print("  CONFIRMED (sampled): Omega_2 = TT-span at n=8")


# ============================================================
# ANALYSIS: Why does TT-span hold?
# ============================================================
print(f"\n{'='*70}")
print("ANALYSIS: Why TT-span holds - constraint structure")
print("=" * 70)

# For an NT path (a,b,c) with c->a:
# The boundary face (a,c) = delta_1 is NOT an arc (c->a, not a->c).
# For (a,b,c) to be in Omega_2, we need d(a,b,c) in span(Omega_1).
# d(a,b,c) = (b,c) - (a,c) + (a,b)
# Since (a,c) is NOT in R_1 (allowed 1-paths), we need the coefficient of (a,c)
# to be cancelled by other NT paths sharing the same junk face.
#
# Paths sharing junk face (a,c): all NT paths (a,w,c) where a->w->c and c->a.
# The number of such paths = A^2[a,c] = #{w: a->w, w->c}.
#
# In Omega_2, the linear constraint is: sum of coefficients on all paths
# sharing junk face (a,c) must be 0. So the span of NT paths in Omega_2
# has dimension |NT| - J_2 (each junk pair gives one constraint).
#
# But this is NOT the same as Omega_2 dimension! The TT paths also
# have constraints among themselves (from the boundary map to Omega_1).
#
# Wait, do TT paths have NO constraints?
# d(TT path (a,b,c)) = (b,c) - (a,c) + (a,b)
# All three are valid 1-paths (arcs). So d(a,b,c) is in R_1.
# For d(a,b,c) to be in Omega_1, we need d(a,b,c) to be in im(compute_omega_basis(1,...))
# But Omega_1 = R_1 for tournaments! Because every allowed 1-path (i,j) has
# boundary (j) - (i) which is in R_0. So Omega_1 = {all 1-chains on arcs}.
# Since all arcs exist, Omega_1 = R^{C(n,2)}.
#
# Therefore d(a,b,c) is ALWAYS in Omega_1 for TT paths. No constraint!
# This means ALL TT paths are independently in Omega_2.
# And NT paths contribute to Omega_2 only through the subspace
# where junk cancellation occurs: dim = |NT| - J_2.
#
# BUT: do the NT elements in Omega_2 coincide with TT-span elements?
# The key question: is the TT subspace + NT-with-cancellation = Omega_2?
# And further: is Omega_2 = TT-span exactly (not just a subspace)?

# Let's verify: for n=5, check dim(Omega_2) = |TT| always.
print("\nn=5: checking dim(Omega_2) = |TT|:")
n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

for bits in range(total):
    A = build_adj(n, bits)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap1 = enumerate_allowed_paths(A, n, 1)
    if not ap2:
        continue
    om2 = compute_omega_basis(A, n, 2, ap2, ap1)
    d2 = om2.shape[1] if om2.ndim == 2 else 0
    n_tt = sum(1 for a,b,c in ap2 if A[a][c])
    if d2 != n_tt:
        scores = sorted([sum(row) for row in A])
        print(f"  MISMATCH: bits={bits}, scores={scores}, TT={n_tt}, Om2={d2}")
        break
else:
    print("  CONFIRMED: dim(Omega_2) = |TT| for ALL n=5 tournaments!")

# Same at n=6
print("\nn=6: checking dim(Omega_2) = |TT|:")
n = 6
n_arcs = n*(n-1)//2
total = 1 << n_arcs

mismatches = 0
mismatch_examples = []
t0 = time.time()
for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        dt = time.time() - t0
        print(f"  ... {bits}/{total} ({dt:.0f}s), mismatches={mismatches}")
    A = build_adj(n, bits)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap1 = enumerate_allowed_paths(A, n, 1)
    if not ap2:
        continue
    om2 = compute_omega_basis(A, n, 2, ap2, ap1)
    d2 = om2.shape[1] if om2.ndim == 2 else 0
    n_tt = sum(1 for a,b,c in ap2 if A[a][c])
    if d2 != n_tt:
        mismatches += 1
        if len(mismatch_examples) < 3:
            scores = sorted([sum(row) for row in A])
            mismatch_examples.append((bits, scores, n_tt, d2))

dt = time.time() - t0
print(f"\n  Completed in {dt:.0f}s. Mismatches: {mismatches}/{total}")

if mismatches == 0:
    print("  CONFIRMED: dim(Omega_2) = |TT| for ALL n=6 tournaments!")
    print("  This means: TT paths form a BASIS of Omega_2")
    print("  NT paths contribute ZERO additional dimension (all constrained away)")
else:
    print(f"  FOUND {mismatches} mismatches!")
    for bits, scores, n_tt, d2 in mismatch_examples:
        print(f"    bits={bits}, scores={scores}, TT={n_tt}, Om2={d2}, diff={d2-n_tt}")


# ============================================================
# KEY THEOREM CANDIDATE: dim(Omega_2) = |TT| = C(n,3) - 2*c_3
# ============================================================
print(f"\n{'='*70}")
print("KEY: dim(Omega_2) = |TT| = C(n,3) - 2*c_3")
print("=" * 70)

# |TT| = number of transitive triples = C(n,3) - c_3
# Wait: for each triple {a,b,c}, it's either:
#   - A 3-cycle (3 NT paths, contributes 0 TT)
#   - A transitive triple (2 TT paths and 1 NT path... no)
#
# Actually, each ordered triple (a,b,c) with a->b, b->c is an allowed 2-path.
# It's TT if a->c, NT if c->a.
# For a TRIPLE {x,y,z}, the tournament restricted to it is either transitive or a 3-cycle.
# - If transitive: say x->y->z->x is NOT possible, so it's x->y, y->z, x->z (for some ordering).
#   The allowed 2-paths from this triple: (x,y,z) with x->z [TT], (x,z,y) NOT allowed (z!->y).
#   Actually no: allowed means consecutive arcs. (x,y,z) needs x->y and y->z. TT needs x->z.
#   For a transitive triple x>y>z (meaning x->y, x->z, y->z):
#   Allowed 2-paths: (x,y,z) [x->y, y->z, and x->z: TT], (z,x,y) [z->x? NO: x->z]
#   Hmm, let me think. If x->y, x->z, y->z (transitive):
#   2-paths: need a->b, b->c.
#   (x,y,z): x->y YES, y->z YES. a->c = x->z YES. TT.
#   (y,x,z): y->x? NO (x->y). Not allowed.
#   (z,y,x): z->y? NO (y->z). Not allowed.
#   (x,z,y): x->z YES, z->y? NO. Not allowed.
#   (z,x,y): z->x? NO. Not allowed.
#   (y,z,x): y->z YES, z->x? NO. Not allowed.
#   So transitive triple gives exactly 1 TT path and 0 NT paths.
#
# - If 3-cycle: say x->y, y->z, z->x.
#   2-paths:
#   (x,y,z): x->y YES, y->z YES. x->z? No (z->x). NT.
#   (y,z,x): y->z YES, z->x YES. y->x? No (x->y). NT.
#   (z,x,y): z->x YES, x->y YES. z->y? No (y->z). NT.
#   So 3-cycle gives 3 NT paths and 0 TT paths.
#
# Total: |TT| = #{transitive triples} * 1 = C(n,3) - c_3
#        |NT| = c_3 * 3 = 3*c_3
#        |A_2| = |TT| + |NT| = C(n,3) - c_3 + 3*c_3 = C(n,3) + 2*c_3

print("Formula verification:")
print("  |TT| = C(n,3) - c_3")
print("  |NT| = 3 * c_3")
print("  |A_2| = C(n,3) + 2*c_3")
print()

from math import comb

# Verify at n=5
n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

for bits in range(min(total, 100)):
    A = build_adj(n, bits)
    ap2 = enumerate_allowed_paths(A, n, 2)
    n_tt = sum(1 for a,b,c in ap2 if A[a][c])
    n_nt = sum(1 for a,b,c in ap2 if not A[a][c])
    # Count c_3
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j]+A[j][i])*(A[j][k]+A[k][j])*(A[i][k]+A[k][i]) > 0:
                    # Check if it's a 3-cycle
                    s = A[i][j] + A[j][k] + A[i][k]
                    if s == 1 or s == 2:  # not all same direction = 3-cycle
                        # Actually: transitive iff one vertex beats both others
                        # 3-cycle iff none beats both
                        score_i = A[i][j] + A[i][k]
                        score_j = A[j][i] + A[j][k]
                        score_k = A[k][i] + A[k][j]
                        if max(score_i, score_j, score_k) < 2:
                            c3 += 1
    predicted_tt = comb(n, 3) - c3
    predicted_nt = 3 * c3
    if n_tt != predicted_tt or n_nt != predicted_nt:
        print(f"  FORMULA ERROR at bits={bits}: TT={n_tt} vs {predicted_tt}, NT={n_nt} vs {predicted_nt}")
        break
else:
    print(f"  Verified |TT| = C(n,3)-c_3 and |NT| = 3*c_3 for 100 n=5 tournaments")

# The KEY result: dim(Omega_2) = |TT| = C(n,3) - c_3
# Combined with |A_2| - J_2 formula: J_2 = |A_2| - |TT| = (C(n,3) + 2*c_3) - (C(n,3) - c_3) = 3*c_3
# Wait: J_2 = 3*c_3?? That would mean J_2 = |NT|!
print(f"\n  If dim(Omega_2) = |TT|, then |A_2| - J_2 = |TT| = C(n,3) - c_3")
print(f"  So J_2 = |A_2| - |TT| = (C(n,3) + 2*c_3) - (C(n,3) - c_3) = 3*c_3 = |NT|")
print(f"  This means EVERY NT path is individually constrained (A^2[a,c] = 1 always?)")
print(f"  OR: the constraints kill exactly as many dimensions as there are NT paths")

# Verify J_2 = 3*c_3 at n=5
def compute_J2(A, n):
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

print(f"\nVerifying J_2 = 3*c_3 at n=5:")
n = 5
total = 1 << (n*(n-1)//2)
j2_eq_3c3 = 0
for bits in range(total):
    A = build_adj(n, bits)
    j2 = compute_J2(A, n)
    # Count c3
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                si = A[i][j] + A[i][k]
                sj = A[j][i] + A[j][k]
                sk = A[k][i] + A[k][j]
                if max(si, sj, sk) < 2:
                    c3 += 1
    if j2 == 3 * c3:
        j2_eq_3c3 += 1
    else:
        print(f"  MISMATCH: bits={bits}, J2={j2}, 3*c3={3*c3}")
        break
else:
    print(f"  CONFIRMED: J_2 = 3*c_3 for ALL n=5 tournaments ({j2_eq_3c3}/{total})")

# Also check n=6
print(f"\nVerifying J_2 = 3*c_3 at n=6:")
n = 6
total = 1 << (n*(n-1)//2)
mismatches_j2 = 0
t0 = time.time()
for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        dt = time.time() - t0
        print(f"  ... {bits}/{total} ({dt:.0f}s), mismatches={mismatches_j2}")
    A = build_adj(n, bits)
    j2 = compute_J2(A, n)
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                si = A[i][j] + A[i][k]
                sj = A[j][i] + A[j][k]
                sk = A[k][i] + A[k][j]
                if max(si, sj, sk) < 2:
                    c3 += 1
    if j2 != 3 * c3:
        mismatches_j2 += 1
        if mismatches_j2 <= 3:
            scores = sorted([sum(row) for row in A])
            print(f"  MISMATCH: bits={bits}, scores={scores}, J2={j2}, 3*c3={3*c3}")

dt = time.time() - t0
print(f"\n  Completed in {dt:.0f}s. Mismatches: {mismatches_j2}/{total}")
if mismatches_j2 == 0:
    print("  CONFIRMED: J_2 = 3*c_3 for ALL n=6 tournaments!")
    print("  Combined with dim(Omega_2) = |TT|:")
    print("  dim(Omega_2) = C(n,3) - c_3 for ALL tournaments!")


# ============================================================
# PROOF SKETCH (if all confirmed)
# ============================================================
print(f"\n{'='*70}")
print("PROOF SKETCH: Why dim(Omega_2) = C(n,3) - c_3")
print("=" * 70)

print("""
Claim: For any tournament T on n vertices, dim(Omega_2(T)) = C(n,3) - c_3(T).

Proof:
1. The allowed 2-paths A_2 consist of:
   - TT paths: (a,b,c) with a->b, b->c, a->c (from transitive triples)
   - NT paths: (a,b,c) with a->b, b->c, c->a (from directed 3-cycles)

2. |TT| = C(n,3) - c_3, |NT| = 3*c_3, |A_2| = C(n,3) + 2*c_3.

3. For a TT path (a,b,c), all boundary faces (b,c), (a,c), (a,b) are arcs.
   So d(a,b,c) is in R_1 = Omega_1 automatically.
   Therefore every TT path is in Omega_2.

4. For an NT path (a,b,c) with c->a, the face (a,c) is NOT an arc.
   So d(a,b,c) has a component outside R_1.
   An NT path is in Omega_2 only if this "junk" component cancels.

5. The junk cancellation constraint: for each junk pair (a,c) with c->a,
   the sum of coefficients of all NT paths (a,w,c) through w must be 0.
   Each such path is from a different 3-cycle through {a,c}.

6. KEY CLAIM: J_2 = 3*c_3. This means:
   - Every NT path has its junk face (a,c) with A^2[a,c] >= 1
   - Each junk pair constraint removes exactly one dimension
   - Total constraints = J_2 = 3*c_3 = |NT|
   - So the NT paths contribute 0 additional dimensions to Omega_2

7. But wait: J_2 counts distinct junk PAIRS, not paths.
   If J_2 = 3*c_3 = |NT|, it means each NT path has a UNIQUE junk pair!
   Equivalently: no two NT paths share the same junk pair.
   Equivalently: for each reverse pair (a,c) with c->a and A^2[a,c]>0,
   there is exactly ONE intermediate b with a->b->c.
   This means A^2[a,c] = 1 for all junk pairs!

8. A^2[a,c] = 1 for all (a,c) with c->a and A^2[a,c]>0?
   This is a strong claim. Let me check...
""")

# Check A^2[a,c] = 1 for all junk pairs
print("Checking A^2[a,c] = 1 for all junk pairs at n=5:")
n = 5
total = 1 << (n*(n-1)//2)
a2_violations = 0
a2_dist = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    for a in range(n):
        for c in range(n):
            if a == c or not A[c][a]:
                continue
            a2 = sum(1 for b in range(n) if b != a and b != c and A[a][b] and A[b][c])
            if a2 > 0:
                a2_dist[a2] += 1
                if a2 != 1:
                    a2_violations += 1

print(f"  A^2 distribution for junk pairs: {dict(sorted(a2_dist.items()))}")
if a2_violations == 0:
    print("  CONFIRMED: A^2[a,c] = 1 for ALL junk pairs at n=5")
    print("  Each NT path has a unique junk pair!")
else:
    print(f"  VIOLATIONS: {a2_violations}")

print("\nChecking A^2[a,c] = 1 for all junk pairs at n=6:")
n = 6
total = 1 << (n*(n-1)//2)
a2_violations_6 = 0
a2_dist_6 = Counter()
t0 = time.time()
for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        dt = time.time() - t0
        if a2_violations_6 > 0:
            print(f"  ... {bits}/{total} ({dt:.0f}s), a2>1 found: {a2_violations_6}")
    A = build_adj(n, bits)
    for a in range(n):
        for c in range(n):
            if a == c or not A[c][a]:
                continue
            a2 = sum(1 for b in range(n) if b != a and b != c and A[a][b] and A[b][c])
            if a2 > 0:
                a2_dist_6[a2] += 1
                if a2 != 1:
                    a2_violations_6 += 1

dt = time.time() - t0
print(f"\n  Completed in {dt:.0f}s")
print(f"  A^2 distribution for junk pairs: {dict(sorted(a2_dist_6.items()))}")
if a2_violations_6 == 0:
    print("  CONFIRMED: A^2[a,c] = 1 for ALL junk pairs at n=6")
else:
    print(f"  A^2 > 1 occurs for {a2_violations_6} junk pairs at n=6")
    # This would mean J_2 < |NT| and dim(Omega_2) > |TT|
    # Let's check if this breaks the TT-span property


print("\nDone.")
