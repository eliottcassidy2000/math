"""
DEFINITIVE test of Grok's "compat-subspace S of dim 3" claim.

Grok says: when β₁(T)=0, H₁(T)=0 imposes exactly 3 global relations
on ⊕_v H₁(T\v), giving a 3D subspace S. Each bad vertex injects into S.

Most concrete interpretation: the "bad indicator vector"
  d(T) = (β₁(T\v₁), ..., β₁(T\vₙ)) ∈ {0,1}^n
always lies in a FIXED 3D subspace of R^n when β₁(T)=0.

Test: collect all bad indicator vectors across all β₁=0 tournaments.
Compute the span. If dim ≤ 3: Grok is right. If dim > 3: refuted.

opus-2026-03-09-S51n
"""
import numpy as np
from collections import Counter

def tournament_from_bits(n, bits):
    A = np.zeros((n,n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def get_tts(A, n):
    tts = []
    for a in range(n):
        for b in range(n):
            if b == a or A[a][b] == 0: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] == 1 and A[a][c] == 1:
                    tts.append((a,b,c))
    return tts

def boundary2_matrix(A, n, edges, tts):
    edge_idx = {e: i for i, e in enumerate(edges)}
    mat = np.zeros((len(edges), len(tts)), dtype=float)
    for j, (a,b,c) in enumerate(tts):
        if (b,c) in edge_idx: mat[edge_idx[(b,c)], j] += 1
        if (a,c) in edge_idx: mat[edge_idx[(a,c)], j] -= 1
        if (a,b) in edge_idx: mat[edge_idx[(a,b)], j] += 1
    return mat

def beta1(A, n):
    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    if not tts: return len(edges) - (n-1)
    mat = boundary2_matrix(A, n, edges, tts)
    r = np.linalg.matrix_rank(mat, tol=1e-8)
    return len(edges) - (n-1) - r


# ============================================================
# TEST 1: Span of bad indicator vectors in R^n
# ============================================================
print("=" * 60)
print("TEST 1: Span of bad indicator vectors d(T) ∈ R^n")
print("Grok claims dim ≤ 3 when β₁(T)=0")
print("=" * 60)

for n in [5, 6, 7]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")

    indicators = []
    indicator_set = set()
    tested = 0

    import random
    if n <= 5:
        bits_iter = range(1 << ne)
        limit = None
    elif n <= 6:
        bits_iter = range(1 << ne)
        limit = None  # exhaustive at n=6
    else:
        bits_iter = (random.randint(0, (1<<ne)-1) for _ in range(50000))
        limit = 5000

    for bits in bits_iter:
        if limit and tested >= limit:
            break

        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        d = np.zeros(n)
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1,n-1), dtype=int)
            for i2,i in enumerate(remaining):
                for j2,j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            if beta1(A_sub, n-1) == 1:
                d[v] = 1

        indicators.append(d)
        indicator_set.add(tuple(d.astype(int)))
        tested += 1

    if not indicators:
        print(f"  No β₁=0 tournaments found")
        continue

    # Compute span
    mat = np.array(indicators)
    rank = np.linalg.matrix_rank(mat, tol=1e-8)

    print(f"  Tested: {tested} tournaments with β₁=0")
    print(f"  Distinct indicator vectors: {len(indicator_set)}")
    print(f"  Span dimension: {rank}")
    print(f"  Grok predicts ≤ 3: {'✓ CONFIRMED' if rank <= 3 else '✗ REFUTED'}")

    # List distinct patterns
    if n <= 6:
        print(f"  Distinct patterns:")
        for pat in sorted(indicator_set):
            cnt = sum(1 for d in indicators if tuple(d.astype(int)) == pat)
            print(f"    {pat}: {cnt} tournaments (weight {sum(pat)})")


# ============================================================
# TEST 2: If span > 3, is there ANY fixed subspace containing all?
# ============================================================
print("\n" + "=" * 60)
print("TEST 2: Minimal subspace containing all bad indicators")
print("=" * 60)

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")

    indicators = []
    for bits in range(1 << ne):
        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        d = np.zeros(n)
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1,n-1), dtype=int)
            for i2,i in enumerate(remaining):
                for j2,j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            if beta1(A_sub, n-1) == 1:
                d[v] = 1
        if np.sum(d) > 0:  # only non-zero indicators
            indicators.append(d)

    mat = np.array(indicators)
    rank = np.linalg.matrix_rank(mat, tol=1e-8)

    # SVD to find the subspace
    U, S, Vt = np.linalg.svd(mat, full_matrices=False)
    basis = Vt[:rank]

    print(f"  Non-zero indicators: {len(indicators)}")
    print(f"  Span dimension: {rank}")
    print(f"  Basis vectors:")
    for i in range(rank):
        print(f"    v{i}: {np.round(basis[i], 4)}")

    # Check: does the all-ones vector lie in the span?
    ones = np.ones(n) / np.sqrt(n)
    proj = basis @ ones
    residual = np.linalg.norm(ones - basis.T @ proj)
    print(f"  All-ones in span: {'YES' if residual < 1e-8 else 'NO'} (residual={residual:.6f})")


# ============================================================
# TEST 3: The "3 global relations from H₁(T)=0"
# What exactly are the constraints that β₁(T)=0 imposes?
# ============================================================
print("\n" + "=" * 60)
print("TEST 3: What constraints does β₁(T)=0 impose on d(T)?")
print("=" * 60)

print("""
For β₁(T)=0: rank(∂₂|Ω₂) = C(n,2) - (n-1).
This gives dim(B₁) = C(n,2) - (n-1) = (n-1)(n-2)/2.

For β₁(T)=1: rank(∂₂|Ω₂) = C(n,2) - n.
dim(B₁) = C(n,2) - n.

The DIFFERENCE is 1: β₁=0 means one MORE boundary relation than β₁=1.
But this is about the chain complex of T, not about the indicator vector d(T).

The indicator d(T) = (β₁(T\\v₁), ..., β₁(T\\vₙ)) is a DERIVED quantity.
There's no obvious "3 global relations" connecting β₁(T)=0 to d(T).

In fact, the constraint β₁(T)=0 is a SINGLE equation (rank condition),
not "3 relations." The number 3 is the maximum of Σ d(v), not the codimension.
""")

# Verify: what's the relationship between β₁(T) and Σ β₁(T\v)?
n = 5
ne = n*(n-1)//2
print(f"\nn=5 exhaustive: β₁(T) vs Σ_v β₁(T\\v)")
cross = Counter()
for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    b1 = beta1(A, n)
    s = 0
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        if beta1(A_sub, n-1) == 1:
            s += 1
    cross[(b1, s)] += 1

for (b1, s), cnt in sorted(cross.items()):
    print(f"  β₁={b1}, Σβ₁(T\\v)={s}: {cnt}")


# ============================================================
# TEST 4: For β₁=0, which vertex SUBSETS can be bad?
# ============================================================
print("\n" + "=" * 60)
print("TEST 4: Which subsets of vertices can be bad when β₁=0?")
print("=" * 60)

n = 5
ne = n*(n-1)//2
bad_subsets = Counter()
for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue
    bad = tuple(v for v in range(n)
                if beta1(np.array([[A[r][c] for c in [x for x in range(n) if x!=v]]
                                   for r in [x for x in range(n) if x!=v]], dtype=int), n-1) == 1)
    bad_subsets[bad] += 1

print(f"n=5 β₁=0: bad subsets and their counts:")
by_size = {}
for subset, cnt in sorted(bad_subsets.items()):
    sz = len(subset)
    if sz not in by_size:
        by_size[sz] = []
    by_size[sz].append((subset, cnt))

for sz in sorted(by_size):
    print(f"\n  Size {sz} ({len(by_size[sz])} distinct subsets):")
    for subset, cnt in by_size[sz]:
        print(f"    {subset}: {cnt}")

# Check: which SIZE-3 subsets appear?
# Are they ALL C(5,3)=10 triples, or only some?
all_triples = set()
for subset, cnt in bad_subsets.items():
    if len(subset) == 3:
        all_triples.add(subset)
print(f"\n  3-bad subsets: {len(all_triples)} out of C(5,3)={5*4*3//6}")
print(f"  ALL triples appear: {len(all_triples) == 10}")

# Similarly for size 1 and 2
for sz in [1, 2]:
    subsets = [s for s in bad_subsets if len(s) == sz]
    from itertools import combinations
    all_possible = list(combinations(range(5), sz))
    print(f"  {sz}-bad subsets: {len(subsets)} out of C(5,{sz})={len(all_possible)}")
    print(f"  ALL appear: {len(subsets) == len(all_possible)}")
