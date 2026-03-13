"""
vitali_sign_pattern_analysis.py -- kind-pasteur-2026-03-13-S61

The delta_sigma for each SX pair is EXACTLY ±1.
Key question: what determines the SIGN PATTERN?

For SX pair (s,x): delta_sigma = -sum_{w in S\{s}} sign(s->w)*sign(x->w)
This sum has 3 terms, each ±1, giving sum in {-3,-1,1,3}.
So delta_sigma in {-3,-1,1,3}.

But our data shows delta_sigma ALWAYS ±1 (never ±3).
This means the 3 terms sum to ±1 (two +1 and one -1, or vice versa).

Why? The (1,1,2,2) score constraint forces specific structure within S.

Analysis:
1. What determines whether delta_sigma = +1 or -1 for each SX pair?
2. Is the sign pattern related to dc7?
3. What is the "hidden geometry" of the sign pattern?
"""

import numpy as np
from itertools import combinations
from collections import Counter, defaultdict

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

def count_directed_k_cycles(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

n = 7
total_bits = n * (n-1) // 2

print("=" * 60)
print("SIGN PATTERN ANALYSIS OF DELTA_SIGMA")
print("=" * 60)

np.random.seed(42)
examples = []

for trial in range(5000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(L, lambda_graph(B, n)):
            continue

        S = set(subset)
        outside = sorted(k for k in range(n) if k not in S)

        # Compute sign pattern
        sign_pattern = {}
        raw_sums = {}
        for s in sorted(S):
            for x in outside:
                dsig = 0
                term_list = []
                for w in sorted(S):
                    if w == s:
                        continue
                    sign_sw = int(A[s][w]) - int(A[w][s])
                    sign_xw = int(A[x][w]) - int(A[w][x])
                    term = -sign_sw * sign_xw
                    dsig += term
                    term_list.append(term)
                sign_pattern[(s, x)] = dsig
                raw_sums[(s, x)] = tuple(term_list)

        # Compute dc7
        c7_A = count_directed_k_cycles(A, n, 7)
        c7_B = count_directed_k_cycles(B, n, 7)
        dc7 = c7_B - c7_A

        # Analyze the (1,1,2,2) tournament structure
        # Find which vertices have score 1 (losers) and score 2 (winners)
        s_list = sorted(S)
        internal_scores = {}
        for i, si in enumerate(s_list):
            sc = sum(A[si][sj] for sj in s_list if sj != si)
            internal_scores[si] = sc
        losers = [v for v in s_list if internal_scores[v] == 1]
        winners = [v for v in s_list if internal_scores[v] == 2]

        # The 3-cycle structure within S: (1,1,2,2) has exactly 1 directed 3-cycle
        # The 4th arc connects the two "middle" vertices
        # Actually: (1,1,2,2) means two losers (outdeg 1) and two winners (outdeg 2)

        examples.append({
            'bits': bits,
            'subset': subset,
            'sign_pattern': sign_pattern,
            'raw_sums': raw_sums,
            'dc7': dc7,
            'losers': losers,
            'winners': winners,
            'outside': outside,
            'A': A.copy(),
        })
        break

    if len(examples) >= 200:
        break

print(f"Collected {len(examples)} examples")

# Analysis 1: Sign pattern structure
print(f"\n--- Sign Pattern Structure ---")

# For each example, represent the sign pattern as a matrix |S| x |outside|
# S has 4 vertices, outside has n-4=3 vertices => 4x3 = 12 entries, each ±1
sign_matrices = []
for ex in examples:
    S_sorted = sorted(ex['subset'])
    mat = []
    for s in S_sorted:
        row = []
        for x in ex['outside']:
            row.append(ex['sign_pattern'][(s, x)])
        mat.append(row)
    sign_matrices.append(np.array(mat))

# What are the row sums and column sums?
row_sum_dist = Counter()
col_sum_dist = Counter()
for mat in sign_matrices:
    for r in range(4):
        row_sum_dist[int(np.sum(mat[r]))] += 1
    for c in range(3):
        col_sum_dist[int(np.sum(mat[:, c]))] += 1

print(f"Row sums (per s in S): {dict(sorted(row_sum_dist.items()))}")
print(f"Col sums (per x outside): {dict(sorted(col_sum_dist.items()))}")

# Total sum of sign pattern (= sum of all delta_sigma)
total_sums = Counter()
for mat in sign_matrices:
    total_sums[int(np.sum(mat))] += 1
print(f"Total sign sum: {dict(sorted(total_sums.items()))}")

# Matrix rank of sign pattern
rank_dist = Counter()
for mat in sign_matrices:
    r = np.linalg.matrix_rank(mat.astype(float))
    rank_dist[r] += 1
print(f"Sign pattern rank: {dict(sorted(rank_dist.items()))}")

# Analysis 2: Unique sign patterns
unique_patterns = Counter()
for mat in sign_matrices:
    unique_patterns[tuple(mat.flatten())] += 1
print(f"\nUnique sign patterns: {len(unique_patterns)}")
for pat, cnt in unique_patterns.most_common(10):
    mat = np.array(pat).reshape(4, 3)
    print(f"  count={cnt}: {mat.tolist()}")

# Analysis 3: Sign pattern vs dc7
print(f"\n--- Sign Pattern vs dc7 ---")
pattern_dc7 = defaultdict(list)
for ex, mat in zip(examples, sign_matrices):
    key = tuple(mat.flatten())
    pattern_dc7[key].append(ex['dc7'])

# Does the sign pattern determine dc7?
ambiguous = 0
for key, dc7s in pattern_dc7.items():
    if len(set(dc7s)) > 1:
        ambiguous += 1
print(f"Sign patterns: {len(pattern_dc7)}")
print(f"Ambiguous for dc7: {ambiguous}")

if ambiguous > 0:
    print("Sign pattern does NOT determine dc7")
else:
    print("Sign pattern DETERMINES dc7!")
    # Show the mapping
    for key, dc7s in sorted(pattern_dc7.items(), key=lambda x: x[1][0]):
        mat = np.array(key).reshape(4, 3)
        print(f"  dc7={dc7s[0]}: {mat.tolist()}")

# Analysis 4: Is dc7 a function of row sums or column sums?
print(f"\n--- dc7 vs Row/Col Sum Structure ---")
rc_dc7 = defaultdict(list)
for ex, mat in zip(examples, sign_matrices):
    row_sums = tuple(sorted(int(np.sum(mat[r])) for r in range(4)))
    col_sums = tuple(sorted(int(np.sum(mat[:, c])) for c in range(3)))
    rc_dc7[(row_sums, col_sums)].append(ex['dc7'])

for key, dc7s in sorted(rc_dc7.items()):
    dc7_dist = Counter(dc7s)
    print(f"  rows={key[0]}, cols={key[1]}: dc7={dict(sorted(dc7_dist.items()))} ({len(dc7s)} ex)")

# Analysis 5: Role of losers vs winners in sign pattern
print(f"\n--- Loser/Winner Role in Sign Pattern ---")
loser_row_sums = Counter()
winner_row_sums = Counter()
for ex, mat in zip(examples, sign_matrices):
    S_sorted = sorted(ex['subset'])
    for r, s in enumerate(S_sorted):
        rs = int(np.sum(mat[r]))
        if s in ex['losers']:
            loser_row_sums[rs] += 1
        else:
            winner_row_sums[rs] += 1

print(f"Loser row sums: {dict(sorted(loser_row_sums.items()))}")
print(f"Winner row sums: {dict(sorted(winner_row_sums.items()))}")

# Analysis 6: The raw term structure
print(f"\n--- Raw Term Analysis ---")
# Each SX pair (s,x) has 3 terms, each ±1.
# The (1,1,2,2) constraint on S means sign(s->w) is structured.
# For a loser s (score 1): s beats exactly 1 of the other 3 S-vertices
#   So sign(s->w) = +1 for 1 vertex, -1 for 2 vertices
# For a winner s (score 2): s beats exactly 2 of the other 3 S-vertices
#   So sign(s->w) = +1 for 2 vertices, -1 for 1 vertex

# This means: for loser s, sum_{w} sign(s->w) = 1-2 = -1
#              for winner s, sum_{w} sign(s->w) = 2-1 = +1

# So delta_sigma(s,x) = -sum_{w!=s} sign(s->w)*sign(x->w)
# For loser s: sign(s->w) = (+1, -1, -1) in some order
# For winner s: sign(s->w) = (+1, +1, -1) in some order

# The sign of delta_sigma depends on how sign(x->w) aligns with sign(s->w)

# Let's check: does the alignment between s-pattern and x-pattern
# have a nice interpretation?

# For each outside vertex x: sign(x->w) for w in S gives a vector in {-1,+1}^3
# The delta_sigma(s,x) = -<sign_s, sign_x> (dot product)

# Verify this interpretation
for ex in examples[:3]:
    S_sorted = sorted(ex['subset'])
    print(f"\n  Example: S={ex['subset']}, losers={ex['losers']}, winners={ex['winners']}")

    # sign vectors for S-vertices
    for s in S_sorted:
        others = [w for w in S_sorted if w != s]
        signs = [int(ex['A'][s][w]) - int(ex['A'][w][s]) for w in others]
        print(f"    s={s} ({'L' if s in ex['losers'] else 'W'}): sign_s = {signs}")

    # sign vectors for outside vertices
    for x in ex['outside']:
        signs_x = [int(ex['A'][x][w]) - int(ex['A'][w][x]) for w in S_sorted]
        print(f"    x={x}: sign_x(full) = {signs_x}")
        # For pair (s, x), we need sign_x restricted to S\{s}
        for s in S_sorted:
            others = [w for w in S_sorted if w != s]
            sign_s = [int(ex['A'][s][w]) - int(ex['A'][w][s]) for w in others]
            sign_x = [int(ex['A'][x][w]) - int(ex['A'][w][x]) for w in others]
            dot = sum(a*b for a, b in zip(sign_s, sign_x))
            dsig = -dot
            print(f"      (s={s},x={x}): sign_s={sign_s}, sign_x={sign_x}, dot={dot}, dsig={dsig}")

# Analysis 7: The "quadratic form" structure
print(f"\n{'='*60}")
print("QUADRATIC FORM INTERPRETATION")
print(f"{'='*60}")

# delta_sigma(s,x) = -sum_{w in S, w!=s} sign(s->w)*sign(x->w)
# This is a BILINEAR form in (sign_s, sign_x)!
#
# Define for each vertex v and each w in S:
#   e_v(w) = sign(v->w) = A[v][w] - A[w][v]
#
# Then delta_sigma(s,x) = -sum_{w in S, w!=s} e_s(w)*e_x(w)
#                        = -<e_s|_{S\{s}}, e_x|_{S\{s}}>
#
# The total dc7 must be some function of these bilinear products.
#
# Key insight: e_s(w) for w in S\{s} is the "internal signature" of s within S.
# e_x(w) for w in S is the "external projection" of x onto S.
#
# The sign pattern matrix M[s,x] = delta_sigma(s,x) encodes the
# "internal-external alignment" between S and the complement.

# Does sum of delta_sigma determine dc7?
print(f"\ndc7 vs features of sign matrix:")
features_dc7 = defaultdict(list)
for ex, mat in zip(examples, sign_matrices):
    # Features:
    total = int(np.sum(mat))
    abs_total = int(np.sum(np.abs(mat)))  # always 12
    pos_count = int(np.sum(mat > 0))
    neg_count = int(np.sum(mat < 0))

    # Structured features
    S_sorted = sorted(ex['subset'])
    loser_sum = sum(int(np.sum(mat[r])) for r, s in enumerate(S_sorted) if s in ex['losers'])
    winner_sum = sum(int(np.sum(mat[r])) for r, s in enumerate(S_sorted) if s in ex['winners'])

    key = (total, loser_sum, winner_sum)
    features_dc7[key].append(ex['dc7'])

for key, dc7s in sorted(features_dc7.items()):
    dc7_dist = Counter(dc7s)
    print(f"  total={key[0]:+d}, L_sum={key[1]:+d}, W_sum={key[2]:+d}: dc7={dict(sorted(dc7_dist.items()))}")

# Analysis 8: Compute dc7 directly from sign pattern using the known formula
print(f"\n{'='*60}")
print("DIRECT dc7 FORMULA SEARCH")
print(f"{'='*60}")

# We know dc7 != linear(delta_sigma). Try QUADRATIC features.
# The sign matrix M is 4x3 with entries ±1.
# Quadratic features: M⊗M has 144 entries, but by symmetry we can reduce.

# Actually, let's think about this differently.
# The sign matrix M[s,x] = delta_sigma(s,x).
# dc7 involves 7-cycles. Under arc reversal within S:
#   - All arcs within S (6 arcs) flip
#   - All arcs between S and outside (12 arcs) stay
#
# A 7-cycle must pass through all 7 vertices.
# Under the flip, the cycle changes when it traverses arcs within S.
#
# The number of S-internal arcs in a 7-cycle on 7 vertices passing through all 4 S-vertices:
# The 4 S-vertices partition the cycle into 4 "segments" separated by S-vertices.
# Each segment is a path through outside vertices.
# Since there are only 3 outside vertices, and the cycle visits all 7,
# the 4 segments connect S→outside paths.

# Let me think about this combinatorially...
# A Hamiltonian cycle visits all 7 vertices. The 4 S-vertices divide it into 4 arcs
# between consecutive S-vertices (wrapping around).
# Some of these arcs may go directly S->S (internal arc), others go S->X->...->S.

# Count arcs WITHIN S that the 7-cycle uses:
# The 7-cycle has exactly 7 arcs. Some connect S-S, some S-X, some X-S, some X-X.
# Since we visit 4 S-vertices and 3 X-vertices in a cycle of length 7:
#   The number of S-S arcs can be 0, 1, 2, 3, or 4.

# Under the flip, each S-S arc reverses direction.
# An arc s1->s2 that was part of the cycle becomes s2->s1.
# So a directed Hamiltonian cycle that used k arcs within S:
#   - If ALL k arcs are reversed, the cycle may or may not remain valid.

# This is getting complex. Let me just try to find the formula empirically.
# Use all 200 examples and try polynomial regression.

# Feature set: entries of the sign matrix + their products
X_rows = []
y_rows = []
for ex, mat in zip(examples, sign_matrices):
    features = list(mat.flatten())  # 12 entries
    # Add pairwise products of entries
    flat = mat.flatten()
    for i in range(len(flat)):
        for j in range(i, len(flat)):
            features.append(flat[i] * flat[j])
    X_rows.append(features)
    y_rows.append(ex['dc7'])

X = np.array(X_rows, dtype=float)
y = np.array(y_rows, dtype=float)

from numpy.linalg import lstsq, matrix_rank

print(f"Feature matrix: {X.shape}, rank: {matrix_rank(X)}")
coeffs, _, _, _ = lstsq(X, y, rcond=None)
pred = X @ coeffs
err = np.max(np.abs(pred - y))
print(f"Quadratic fit max error: {err:.6f}")
print(f"Perfect fit? {err < 0.001}")

# Try with just the 12 linear features
X_lin = np.array([mat.flatten().tolist() for mat in sign_matrices], dtype=float)
coeffs_lin, _, _, _ = lstsq(X_lin, y, rcond=None)
err_lin = np.max(np.abs(X_lin @ coeffs_lin - y))
print(f"\nLinear fit max error: {err_lin:.6f}")

# What about the OUTSIDE vertex's total sign w.r.t. S?
# For each x, define phi(x) = sum_{w in S} sign(x->w)
# phi(x) in {-4,-2,0,2,4} (sum of 4 terms each ±1)
# But actually we know the score of x in the full tournament, not just its S-arcs.

# Let's try: dc7 = f(phi_x1, phi_x2, phi_x3) where phi_xi = sum_{w in S} sign(xi->w)
print(f"\n--- Phi Analysis ---")
phi_dc7 = defaultdict(list)
for ex in examples:
    S_sorted = sorted(ex['subset'])
    phis = []
    for x in ex['outside']:
        phi = sum(int(ex['A'][x][w]) - int(ex['A'][w][x]) for w in S_sorted)
        phis.append(phi)
    phi_dc7[tuple(sorted(phis))].append(ex['dc7'])

for key, dc7s in sorted(phi_dc7.items()):
    dc7_dist = Counter(dc7s)
    print(f"  phi={key}: dc7={dict(sorted(dc7_dist.items()))} ({len(dc7s)} ex)")

# Analysis 9: The full external projection
# For each outside vertex x, its "S-projection" is (sign(x->w1), ..., sign(x->w4)) in {-1,+1}^4
# Combined: 3 vectors in {-1,+1}^4 = point in {-1,+1}^{12}
# This 12-bit binary code may determine dc7

print(f"\n--- Full External Projection ---")
proj_dc7 = defaultdict(list)
for ex in examples:
    S_sorted = sorted(ex['subset'])
    projs = []
    for x in ex['outside']:
        proj = tuple(int(ex['A'][x][w]) - int(ex['A'][w][x]) for w in S_sorted)
        projs.append(proj)
    # Sort by vertex to get canonical form (but this depends on x labeling)
    # Actually, keep x-ordering as is (sorted outside vertices)
    key = tuple(projs)
    proj_dc7[key].append(ex['dc7'])

print(f"Distinct projection patterns: {len(proj_dc7)}")
ambig = sum(1 for dc7s in proj_dc7.values() if len(set(dc7s)) > 1)
print(f"Ambiguous for dc7: {ambig}")

if ambig == 0:
    print("External projection DETERMINES dc7!")
    for key, dc7s in sorted(proj_dc7.items(), key=lambda x: x[1][0]):
        if len(dc7s) >= 3:
            print(f"  proj={key}: dc7={dc7s[0]} ({len(dc7s)} ex)")
else:
    # Show ambiguous cases
    for key, dc7s in sorted(proj_dc7.items()):
        if len(set(dc7s)) > 1:
            print(f"  AMBIG proj={key}: dc7={Counter(dc7s)}")

# Analysis 10: The INTERNAL structure of S combined with external projection
print(f"\n{'='*60}")
print("INTERNAL + EXTERNAL = FULL PICTURE")
print(f"{'='*60}")

# The (1,1,2,2) tournament on 4 vertices has specific structure.
# There are 4 tournaments on 4 vertices with score (1,1,2,2):
# Let's enumerate them.
for bits4 in range(1 << 6):
    A4 = bits_to_adj(bits4, 4)
    sc = tuple(sorted(int(sum(A4[i])) for i in range(4)))
    if sc == (1, 1, 2, 2):
        # Print the tournament
        arcs = []
        for i in range(4):
            for j in range(i+1, 4):
                if A4[i][j]:
                    arcs.append(f"{i}->{j}")
                else:
                    arcs.append(f"{j}->{i}")
        # Which 3-cycles exist?
        c3 = 0
        for a, b, c in [(0,1,2), (0,1,3), (0,2,3), (1,2,3)]:
            if (A4[a][b]*A4[b][c]*A4[c][a] + A4[c][b]*A4[b][a]*A4[a][c]) > 0:
                c3 += 1
        print(f"  bits={bits4:06b}: arcs={arcs}, c3={c3}")

# There should be exactly 4! / |Aut| = 12 or so labeled tournaments
# Let me count properly
count_1122 = 0
tournament_classes = []
for bits4 in range(1 << 6):
    A4 = bits_to_adj(bits4, 4)
    sc = tuple(sorted(int(sum(A4[i])) for i in range(4)))
    if sc == (1, 1, 2, 2):
        count_1122 += 1
        tournament_classes.append(bits4)

print(f"\n(1,1,2,2) tournaments on 4 vertices: {count_1122}")

# For each, what does the internal S-structure + external projection give?
# The dc7 formula must depend on:
# (a) which (1,1,2,2) tournament S induces, and
# (b) how the outside vertices project onto S

# Let's classify by internal tournament type + external phi
internal_ext_dc7 = defaultdict(list)
for ex in examples:
    S_sorted = sorted(ex['subset'])
    # Internal tournament encoding: arcs between S-vertices
    A = ex['A']
    internal = tuple(int(A[s1][s2]) for s1 in S_sorted for s2 in S_sorted if s1 != s2)

    # External: for each x, its score INTO S (how many S-vertices does x beat?)
    ext = []
    for x in ex['outside']:
        score_into_S = sum(int(A[x][s]) for s in S_sorted)
        ext.append(score_into_S)

    key = (internal, tuple(sorted(ext)))
    internal_ext_dc7[key].append(ex['dc7'])

print(f"\n(Internal, sorted-ext-scores) groups: {len(internal_ext_dc7)}")
ambig = sum(1 for dc7s in internal_ext_dc7.values() if len(set(dc7s)) > 1)
print(f"Ambiguous for dc7: {ambig}")

# Even more detail: external FULL projection (not just score)
internal_fullext_dc7 = defaultdict(list)
for ex in examples:
    S_sorted = sorted(ex['subset'])
    A = ex['A']
    internal = tuple(int(A[s1][s2]) for s1 in S_sorted for s2 in S_sorted if s1 != s2)

    ext = []
    for x in ex['outside']:
        proj = tuple(int(A[x][s]) for s in S_sorted)
        ext.append(proj)

    key = (internal, tuple(ext))
    internal_fullext_dc7[key].append(ex['dc7'])

print(f"\n(Internal, full-ext-projection) groups: {len(internal_fullext_dc7)}")
ambig = sum(1 for dc7s in internal_fullext_dc7.values() if len(set(dc7s)) > 1)
print(f"Ambiguous for dc7: {ambig}")

if ambig > 0:
    # Show ambiguous
    for key, dc7s in sorted(internal_fullext_dc7.items()):
        if len(set(dc7s)) > 1:
            print(f"  AMBIG: dc7={Counter(dc7s)}")
            break
else:
    print("(Internal, full ext projection) DETERMINES dc7!")

# What about the OUTSIDE-OUTSIDE arcs?
internal_fullext_xx_dc7 = defaultdict(list)
for ex in examples:
    S_sorted = sorted(ex['subset'])
    A = ex['A']
    internal = tuple(int(A[s1][s2]) for s1 in S_sorted for s2 in S_sorted if s1 != s2)

    ext = []
    for x in ex['outside']:
        proj = tuple(int(A[x][s]) for s in S_sorted)
        ext.append(proj)

    # Outside-outside arcs
    xx_arcs = tuple(int(A[ex['outside'][i]][ex['outside'][j]])
                    for i in range(3) for j in range(3) if i != j)

    key = (internal, tuple(ext), xx_arcs)
    internal_fullext_xx_dc7[key].append(ex['dc7'])

print(f"\n(Internal, full-ext, XX-arcs) groups: {len(internal_fullext_xx_dc7)}")
ambig = sum(1 for dc7s in internal_fullext_xx_dc7.values() if len(set(dc7s)) > 1)
print(f"Ambiguous for dc7: {ambig}")

if ambig == 0:
    print("Full local data DETERMINES dc7!")
    print("\nThis means: dc7 depends on the COMPLETE arc data around the Vitali atom.")
    print("The sigma information (Level 1.5) does NOT suffice.")

print("\nDone.")
