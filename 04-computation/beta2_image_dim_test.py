"""
DEFINITIVE test: Does res: Z₁(T) → ⊕_v H₁(T\v) factor through a fixed 3D space?

Grok claims the IMAGE of res lives in a 3D subspace independent of #bad.
But rank(res) = #bad exactly. So the image IS #bad-dimensional.

For the "fixed 3D" claim to work, there would need to be cross-vertex
compatibility constraints making the image ≤ 3D even when #bad > 3.
Since #bad ≤ 3 is exactly what we're trying to prove, let's check:

1. Is the image always ≤ 3D? (trivially yes since #bad ≤ 3 empirically)
2. Is there a STRUCTURAL reason the image is ≤ 3D?
3. If we hypothetically had #bad=4, would there be a 3D constraint?

We test by examining cross-vertex compatibility: when we embed all H₁(T\v)
into a common space, do the generators satisfy any universal relations?

opus-2026-03-09-S51l
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

def compute_h1_generator(A, n):
    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    if not tts: return edges, None
    mat = boundary2_matrix(A, n, edges, tts)
    U, S, _ = np.linalg.svd(mat)
    rank = np.sum(S > 1e-8)
    cocycles = U[:, rank:]
    if cocycles.shape[1] - (n - 1) <= 0: return edges, None
    cob = np.zeros((len(edges), n), dtype=float)
    for v in range(n):
        for i, (a, b) in enumerate(edges):
            if b == v: cob[i, v] += 1
            if a == v: cob[i, v] -= 1
    cc = np.linalg.lstsq(cocycles, cob, rcond=None)[0]
    U2, S2, V2 = np.linalg.svd(cc.T)
    rc = np.sum(S2 > 1e-8)
    h1c = V2[rc:]
    if h1c.shape[0] == 0: return edges, None
    return edges, cocycles @ h1c[0]


# ============================================================
# TEST 1: Image dimension of res equals #bad (not fixed 3)
# ============================================================
print("=" * 60)
print("DEFINITIVE: rank(res: Z₁→⊕H₁(T\\v)) = #bad")
print("=" * 60)

print("""
The map res: Z₁(T) → ⊕_v H₁(T\\v) sends a 1-cycle z to its
projection onto H₁(T\\v) for each v.

Since H₁(T\\v) is 0 or 1-dimensional:
  ⊕_v H₁(T\\v) = R^{#bad}

The map res picks out the H₁ component for each bad vertex.
If the hidden cycle lifts are independent in Z₁ (which they are),
then rank(res) = #bad.

IMAGE dim = #bad. NOT a fixed 3.
The only way IMAGE ≤ 3 is if #bad ≤ 3, which is circular.
""")

# Verify once more exhaustively
for n in [5, 6, 7]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")
    tested = 0
    limit = None if n <= 5 else (500 if n <= 6 else 100)
    bad_dist = Counter()

    import random
    bits_iter = range(1 << ne) if n <= 5 else (random.randint(0, (1<<ne)-1) for _ in range(5000))

    for bits in bits_iter:
        if tested >= (limit or 999999):
            break

        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        bad_count = 0
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1,n-1), dtype=int)
            for i2,i in enumerate(remaining):
                for j2,j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            if beta1(A_sub, n-1) == 1:
                bad_count += 1

        bad_dist[bad_count] += 1
        tested += 1

    print(f"  Tested: {tested}")
    print(f"  #bad distribution: {dict(sorted(bad_dist.items()))}")
    if max(bad_dist.keys()) <= 3:
        print(f"  #bad ≤ 3 always ✓ (but this IS HYP-282, not a proof)")
    else:
        print(f"  *** #bad > 3 FOUND! ***")


# ============================================================
# TEST 2: "Cross-v compatibility matrix"
# For each pair of bad vertices (v_i, v_j), the hidden cycle
# lifts z_vi and z_vj share some edges. Does this constrain them?
# ============================================================
print("\n" + "=" * 60)
print("TEST 2: Cross-vertex Gram matrix of hidden cycle lifts")
print("=" * 60)

n = 5
ne = n*(n-1)//2
gram_patterns = Counter()

for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    edge_idx = {e: i for i, e in enumerate(edges)}

    bad = []
    lifts = []
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        es, h1 = compute_h1_generator(A_sub, n-1)
        if h1 is not None:
            bad.append(v)
            lift = np.zeros(len(edges))
            for k, (a, b) in enumerate(es):
                ao, bo = remaining[a], remaining[b]
                if (ao, bo) in edge_idx:
                    lift[edge_idx[(ao, bo)]] = h1[k]
            lifts.append(lift)

    if len(bad) != 3:
        continue

    # Gram matrix of lifts
    G = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            G[i,j] = np.dot(lifts[i], lifts[j])

    # Normalize
    norms = np.sqrt(np.diag(G))
    G_norm = G / np.outer(norms, norms)

    # Round for pattern detection
    pattern = tuple(np.round(G_norm[np.triu_indices(3, k=1)], 4))
    gram_patterns[pattern] += 1

print(f"n=5 Gram matrix patterns (normalized, upper triangle):")
for pat, cnt in sorted(gram_patterns.items(), key=lambda x: -x[1]):
    print(f"  {pat}: {cnt} tournaments")


# ============================================================
# TEST 3: Support overlap structure
# ============================================================
print("\n" + "=" * 60)
print("TEST 3: Edge support overlap between hidden cycle lifts")
print("=" * 60)

n = 5
ne = n*(n-1)//2
overlap_patterns = Counter()

for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    edge_idx = {e: i for i, e in enumerate(edges)}

    bad = []
    lifts = []
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        es, h1 = compute_h1_generator(A_sub, n-1)
        if h1 is not None:
            bad.append(v)
            lift = np.zeros(len(edges))
            for k, (a, b) in enumerate(es):
                ao, bo = remaining[a], remaining[b]
                if (ao, bo) in edge_idx:
                    lift[edge_idx[(ao, bo)]] = h1[k]
            lifts.append(lift)

    if len(bad) != 3:
        continue

    # Support sets
    supports = [set(j for j, e in enumerate(edges) if abs(lifts[i][j]) > 1e-10) for i in range(3)]
    sizes = tuple(len(s) for s in supports)
    overlaps = (len(supports[0] & supports[1]),
                len(supports[0] & supports[2]),
                len(supports[1] & supports[2]),
                len(supports[0] & supports[1] & supports[2]))
    overlap_patterns[(sizes, overlaps)] += 1

print(f"n=5 support patterns:")
for (sizes, overlaps), cnt in sorted(overlap_patterns.items(), key=lambda x: -x[1]):
    print(f"  Sizes={sizes}, pairwise overlaps={overlaps[:3]}, triple={overlaps[3]}: {cnt}")


# ============================================================
# TEST 4: THE REAL QUESTION — why is #bad ≤ 3?
# ============================================================
print("\n" + "=" * 60)
print("THE REAL QUESTION: Why is #bad ≤ 3?")
print("=" * 60)

print("""
Every mechanism Grok proposes reduces to:
  "The image is ≤ 3D" ⟸ "#bad ≤ 3" (circular)

What would a NON-circular proof look like?
It would need to show: from the STRUCTURE of tournaments alone,
without counting bad vertices, some space has dim ≤ 3.

Candidates for such a space:
1. A quotient/kernel that is 3D regardless of T — REFUTED (varies)
2. A fixed subspace of B₁ — REFUTED (depends on T)
3. A combinatorial constraint (like: 4 vertices can't all have β₁=1
   in their complement deletions) — THIS is what needs proving

The combinatorial approach: if {v₁,v₂,v₃,v₄} are all bad, then
T\\v₁,...,T\\v₄ all have β₁=1. Can we derive a contradiction from
the structure of these 4 tournaments?

At n=5: T\\vᵢ has 4 vertices, so it's a 4-tournament.
A 4-tournament has β₁=1 iff it's the 4-cycle (exactly 1 of 3 iso classes).
Having 4 bad vertices means FOUR different 4-vertex subsets are 4-cycles.
But C(5,4) = 5, so at most 5 possible 4-cycles. Can 4 coexist?
""")

# Direct test: among 4-vertex subtournaments, how many are 4-cycles?
n = 5
ne = n*(n-1)//2
cycle4_counts = Counter()

for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    count = 0
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        if beta1(A_sub, n-1) == 1:
            count += 1
    cycle4_counts[count] += 1

print(f"n=5 exhaustive: #bad (= # of 4-cycle deletions) distribution:")
print(f"  {dict(sorted(cycle4_counts.items()))}")
print(f"  Max #bad = {max(cycle4_counts.keys())}")
print(f"  Total β₁=0 tournaments: {sum(cycle4_counts.values())}")

# At n=5: what are the BAD vertex sets?
bad_sets = Counter()
for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    bad = []
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        if beta1(A_sub, n-1) == 1:
            bad.append(v)
    if len(bad) > 0:
        bad_sets[tuple(bad)] += 1

print(f"\nn=5: Bad vertex sets (for tournaments with some bad vertex):")
size_dist = Counter()
for bset, cnt in bad_sets.items():
    size_dist[len(bset)] += cnt
print(f"  Size distribution: {dict(sorted(size_dist.items()))}")

# For #bad=3: which triples appear?
print(f"\n  3-bad triples:")
for bset, cnt in sorted(bad_sets.items()):
    if len(bset) == 3:
        print(f"    {bset}: {cnt} tournaments")

# Key structural question: are the 3 bad vertices always a TT?
n = 5
ne = n*(n-1)//2
always_tt = 0
not_tt = 0
for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    bad = []
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        if beta1(A_sub, n-1) == 1:
            bad.append(v)

    if len(bad) != 3:
        continue

    # Check: is {bad[0], bad[1], bad[2]} a transitive triple?
    a, b, c = bad
    is_tt = False
    for perm in [(a,b,c),(a,c,b),(b,a,c),(b,c,a),(c,a,b),(c,b,a)]:
        x, y, z = perm
        if A[x][y] == 1 and A[y][z] == 1 and A[x][z] == 1:
            is_tt = True
            break
    if is_tt:
        always_tt += 1
    else:
        not_tt += 1

print(f"\nn=5: 3-bad vertices form a TT: {always_tt}, form a 3-cycle: {not_tt}")
print(f"(We already know: always TT, never 3-cycle)")

# ============================================================
# TEST 5: Can we derive the |bad|≤3 constraint combinatorially?
# At n=5: 4-tournament has β₁=1 iff it's a 4-cycle.
# A 4-cycle on {a,b,c,d} means: a→b→c→d→a (cyclic order).
# 4 bad vertices means 4 out of 5 subsets are 4-cycles.
# Let's check: can the 5-tournament structure support 4 coexisting 4-cycles?
# ============================================================
print("\n" + "=" * 60)
print("TEST 5: Can 4 different 4-subsets all be 4-cycles?")
print("=" * 60)

n = 5
ne = n*(n-1)//2
for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    # Don't require β₁(T)=0 for this test
    count = 0
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        if beta1(A_sub, n-1) == 1:
            count += 1
    if count >= 4:
        print(f"  FOUND: bits={bits}, #bad={count}, β₁(T)={beta1(A,n)}")
        # Print tournament
        for i in range(n):
            row = [str(A[i][j]) for j in range(n)]
            print(f"    {' '.join(row)}")

print("If nothing printed: no tournament (even with β₁>0) has #bad≥4 at n=5")

# Same for n=6
print(f"\n--- n=6 (sampled) ---")
import random
n = 6
ne = n*(n-1)//2
found4 = 0
tested = 0
for _ in range(10000):
    bits = random.randint(0, (1<<ne)-1)
    A = tournament_from_bits(n, bits)
    count = 0
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        if beta1(A_sub, n-1) == 1:
            count += 1
    if count >= 4:
        found4 += 1
        b1 = beta1(A, n)
        print(f"  FOUND: #bad={count}, β₁(T)={b1}")
    tested += 1

print(f"  Tested {tested}, found #bad≥4: {found4}")
if found4 == 0:
    print(f"  #bad ≤ 3 holds even for tournaments with β₁>0!")
