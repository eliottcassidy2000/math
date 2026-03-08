#!/usr/bin/env python3
"""
beta2_proof_strategy.py - Prove beta2=0 using tournament structure

Key insight: For tournaments, between any two vertices a,c with c->a
(backward pair), the number of 2-paths a->b->c equals A^2[a,c].
This is always >= 1 in a tournament (by pigeon).

Wait, actually A^2[a,c] counts walks, and we need at least one
vertex b with a->b and b->c. In a tournament with c->a:
  vertices other than a,c: n-2 total
  each is beaten by a or beats a (can't be both since tournament)
  if a->b (a beats b), then b->c or c->b
  we need a->b AND b->c

  a beats d_a vertices (out-degree of a)
  one of these is c if a->c, but here c->a
  so a beats d_a vertices among {v != a, c}
  of those d_a vertices, how many does c beat?
  c->a already used

  A^2[a,c] = sum_b A[a][b]*A[b][c] where b != a, b != c
  For a junk pair (c->a with A^2[a,c]>0): at least one b exists.

  But can A^2[a,c] = 0 when c->a?
  This means: no b with a->b->c.
  All vertices beaten by a also beat c... wait that means
  all out-neighbors of a are out-neighbors of c.

  a has d_a out-neighbors (not counting c since c->a)
  each of these is either beaten by c (b->c: NO) or beats c (c->b: YES for our search... wait)
  We need b->c (A[b][c]=1). So if a->b and b->c.
  A^2[a,c]=0 means: for all b with a->b, we have c->b (not b->c).

  So a's out-neighbors are all BEATEN by c.
  a beats d_a vertices among {all except c}
  c also beats all of these d_a vertices
  c beats d_c vertices total, and d_c >= d_a + 1 (c beats a plus all of a's targets)
  So d_c >= d_a + 1 (since c->a is one, plus a's d_a out-nbrs are also beaten by c)

  Wait: a has out-degree d_a. These d_a vertices are among vertices != a.
  c->a, so c is NOT in a's out-neighborhood.
  The d_a out-neighbors of a: all have c->b (from A^2[a,c]=0).
  c also has c->a.
  So c beats a + all d_a vertices = d_a + 1 vertices (at minimum).
  These d_a+1 are DISTINCT (a is one, the d_a out-neighbors of a are others).
  But d_a out-neighbors of a don't include a or c.

  So out-degree of c >= d_a + 1.
  But total out-degrees sum to C(n,2), so average is (n-1)/2.

  In a tournament: A^2[a,c] = 0 with c->a is POSSIBLE when d_c >> d_a.
  Example: transitive tournament, a=vertex with max score, c=vertex with min score.
  Then A^2[a,c] = |{b: a->b and b->c}|, but c has in-degree n-1 (beaten by everyone)
  and out-degree 0, so b->c = 0 for all b. A^2[a,c] = 0.

  Actually for transitive tournament: c->a is impossible since a has max score
  and c has min score (a->c).

OK let me just compute things directly.

Strategy:
For each "junk face" (a,c,d) arising from a 3-path (a,b,c,d) with c->a,
the constraint says: sum_{b: a->b->c} alpha(a,b,c,d) = 0
where b ranges over all intermediaries.

The MINIMUM number of such b is 1 (when A^2[a,c]=1).
When there are multiple b's, the constraint has multiple terms.

For beta2=0, we need: rk(bd3|_Om3) = dim(Z2).

Approach: show that the boundary map is surjective by explicit construction.
For each Z2 element, construct a preimage in Om3.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
import numpy as np
from collections import Counter, defaultdict
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


# ============================================================
# KEY TEST: Is beta2=0 a TOURNAMENT-SPECIFIC property?
# Test on NON-tournament digraphs
# ============================================================
print("=" * 70)
print("BETA2 ON NON-TOURNAMENTS")
print("=" * 70)

# Test directed cycles (known: beta1=1 always, beta2=0)
from path_homology_v2 import path_betti_numbers

for n_val in [3, 4, 5, 6, 7, 8]:
    A = [[0]*n_val for _ in range(n_val)]
    for i in range(n_val):
        A[i][(i+1) % n_val] = 1
    betti = path_betti_numbers(A, n_val, max_dim=min(n_val-1, 4))
    print(f"  C_{n_val} (directed cycle): beta = {betti}")

# Test bidirectional graphs (complete undirected = double edges)
print("\nBidirectional complete graphs:")
for n_val in [3, 4, 5]:
    A = [[1 if i != j else 0 for j in range(n_val)] for i in range(n_val)]
    betti = path_betti_numbers(A, n_val, max_dim=min(n_val-1, 4))
    print(f"  K_{n_val} (bidirectional): beta = {betti}")

# Test random digraphs (NOT tournaments)
print("\nRandom dense digraphs (edge prob=0.5):")
import random
random.seed(42)
for n_val in [5, 6, 7]:
    b2_nonzero = 0
    n_tests = 200
    for _ in range(n_tests):
        A = [[0]*n_val for _ in range(n_val)]
        for i in range(n_val):
            for j in range(n_val):
                if i != j and random.random() < 0.5:
                    A[i][j] = 1
        betti = path_betti_numbers(A, n_val, max_dim=2)
        if betti[2] > 0:
            b2_nonzero += 1
    print(f"  n={n_val}: beta2>0 in {b2_nonzero}/{n_tests} random digraphs")


# Test "almost-tournaments" (one missing edge)
print("\nAlmost-tournaments (remove one edge from tournament):")
n = 5
total = 1 << (n*(n-1)//2)
b2_nonzero = 0
b2_examples = []
for bits in range(total):
    A = build_adj(n, bits)
    # Remove each edge in turn
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                A_mod = [row[:] for row in A]
                A_mod[i][j] = 0
                betti = path_betti_numbers(A_mod, n, max_dim=2)
                if betti[2] > 0:
                    b2_nonzero += 1
                    if len(b2_examples) < 3:
                        b2_examples.append((bits, i, j, betti))
                A_mod[i][j] = 1  # restore

print(f"  n=5: beta2>0 in {b2_nonzero} out of {total * n*(n-1)} edge removals")
if b2_examples:
    for bits, i, j, betti in b2_examples:
        print(f"    bits={bits}, removed {i}->{j}: beta={betti}")


# ============================================================
# SPECIFIC PROPERTY: "completeness fills"
# ============================================================
print(f"\n{'='*70}")
print("COMPLETENESS PROPERTY")
print("=" * 70)

print("""
Key property of tournaments for beta2=0:

For any vertices a, c with c->a (backward pair),
and any vertex d with c->d or d->c:

If a->c does NOT hold (backward pair), the "junk face" (a,c,d) must cancel.

In a tournament, for every such pair (a,c) with c->a:
- A^2[a,c] = #{b : a->b->c} >= 0
- If A^2[a,c] >= 1, we have at least one intermediary b
- The cancellation requires: sum of coefficients of all (a,b,c,d) = 0

The COMPLETENESS ensures every vertex is either beaten by a or by c,
creating enough paths to cancel all junk faces.

Similar for d_2 faces (a,b,d) with d->b.
""")

# Verify: in tournaments, J2 counts exactly the junk pairs with A^2>0
n = 5
for bits in range(10):
    A = build_adj(n, bits)
    j2 = 0
    for a in range(n):
        for c in range(n):
            if a == c:
                continue
            if A[c][a]:  # backward pair
                a2 = sum(A[a][b] * A[b][c] for b in range(n) if b != a and b != c)
                if a2 > 0:
                    j2 += 1

    # In a tournament, is A^2[a,c] always > 0 when c->a?
    always_pos = True
    for a in range(n):
        for c in range(n):
            if a == c:
                continue
            if A[c][a]:
                a2 = sum(A[a][b] * A[b][c] for b in range(n) if b != a and b != c)
                if a2 == 0:
                    always_pos = False
                    # This means a has NO intermediate -> c is a "dead end" from a
                    # Possible in tournaments? Let me check.
                    da = sum(A[a])
                    dc = sum(A[c])
                    print(f"  A^2[{a},{c}]=0 in bits={bits}: d_a={da}, d_c={dc}")

    if not always_pos:
        print(f"  bits={bits}: NOT all A^2 positive")


# Check at n=6: how often is A^2[a,c]=0 for backward pairs?
print(f"\nA^2[a,c]=0 rate for backward pairs at n=6:")
n = 6
total = 1 << (n*(n-1)//2)
zero_count = 0
total_backward = 0

for bits in range(total):
    A = build_adj(n, bits)
    for a in range(n):
        for c in range(n):
            if a != c and A[c][a]:
                total_backward += 1
                a2 = sum(A[a][b] * A[b][c] for b in range(n) if b != a and b != c)
                if a2 == 0:
                    zero_count += 1

print(f"  Total backward pairs: {total_backward}")
print(f"  A^2=0 cases: {zero_count} ({100*zero_count/total_backward:.1f}%)")
print(f"  A^2>0 cases: {total_backward-zero_count}")
print(f"  Note: A^2=0 means NO intermediate b exists - this is the transitive closure case")


print("\n\nDone.")
