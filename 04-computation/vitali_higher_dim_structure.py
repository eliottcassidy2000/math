"""
vitali_higher_dim_structure.py -- kind-pasteur-2026-03-13-S61

DEEP EXPLORATION: The hidden higher-dimensional structure in tournaments.

The user asks about "the vitali set and how it relates to the hidden
higher dimensional structure buried in tournaments of n 2,1 and 0."

KEY STRUCTURAL INSIGHTS:
1. The lambda graph lives in a HIGHER-DIMENSIONAL space than the tournament.
   A tournament T on n vertices has C(n,2) binary degrees of freedom.
   The lambda graph L(T) has C(n,2) integer entries (each 0..n-2).
   The map T -> L(T) is many-to-one: multiple tournaments have the same lambda.

2. Vitali atoms are FIBERS of this projection:
   The set of tournaments with the same lambda graph forms an equivalence class.
   Vitali atoms move between tournaments in the same fiber.
   This is analogous to the Vitali set construction where rational translations
   preserve a "modular structure" (membership in Q-coset).

3. The {2,1,0} overlap weights define a SECONDARY projection:
   For each pair of odd cycles (C, C'), the overlap is |V(C) cap V(C')| = 2, 1, or 0.
   The overlap spectrum is a vector in Z^3 counting how many pairs have each overlap.

4. c3 and c5 are "LAMBDA-MEASURABLE" — they factor through the projection T -> L(T).
   c7+ are "LAMBDA-NON-MEASURABLE" — they depend on more than lambda.
   This is exactly like measurable vs non-measurable functions in the Vitali construction!

5. The WITNESS MATRIX W is the hidden higher-dimensional object:
   W is a n x C(n,2) binary matrix where W[k][(u,v)] = 1 iff k witnesses {u,v}.
   Row sums = delta(k) (3-cycles through vertex k).
   Column sums = lambda(u,v).
   The Vitali atom changes individual W entries but preserves column sums.

Let me explore this structure computationally.
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter

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

def witness_matrix(A, n):
    """Build the n x C(n,2) witness matrix."""
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    W = np.zeros((n, len(pairs)), dtype=int)
    for pi, (u, v) in enumerate(pairs):
        for k in range(n):
            if k == u or k == v:
                continue
            if (A[u][v]*A[v][k]*A[k][u] + A[v][u]*A[u][k]*A[k][v]):
                W[k][pi] = 1
    return W, pairs

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

# EXPLORATION 1: Witness matrix fiber structure
print("=" * 60)
print("WITNESS MATRIX FIBER ANALYSIS")
print("=" * 60)

n = 7
total_bits = n * (n-1) // 2
np.random.seed(42)

print(f"\nn={n}: C(n,2)={n*(n-1)//2} pairs, witness matrix is {n} x {n*(n-1)//2}")

# For a Vitali atom: how does the witness matrix change?
found_pairs = 0
w_changes_list = []

for trial in range(2000):
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

        # Found a Vitali pair!
        W_A, pairs = witness_matrix(A, n)
        W_B, _ = witness_matrix(B, n)

        # How does W change?
        diff = W_B - W_A  # entries in {-1, 0, 1}
        n_changed = np.count_nonzero(diff)
        n_plus = np.count_nonzero(diff == 1)
        n_minus = np.count_nonzero(diff == -1)

        # Column sums should be the same (lambda preservation)
        col_change = np.sum(diff, axis=0)
        assert np.all(col_change == 0), "Lambda changed!"

        # Row sum changes?
        row_change = np.sum(diff, axis=1)

        w_changes_list.append({
            'subset': subset,
            'n_changed': n_changed,
            'n_plus': n_plus,
            'n_minus': n_minus,
            'row_changes': tuple(row_change),
            'diff': diff.copy(),
        })
        found_pairs += 1
        break

    if found_pairs >= 100:
        break

print(f"Found {found_pairs} Vitali pairs")
print(f"\nWitness matrix change statistics:")
n_changed_counts = Counter(d['n_changed'] for d in w_changes_list)
print(f"  # of changed entries: {dict(sorted(n_changed_counts.items()))}")
print(f"  # of +1 entries: {Counter(d['n_plus'] for d in w_changes_list)}")
print(f"  # of -1 entries: {Counter(d['n_minus'] for d in w_changes_list)}")

# Row changes (delta(k) changes)
print(f"\nRow sum changes (delta_k change):")
for d in w_changes_list[:5]:
    print(f"  subset={d['subset']}: row_changes={d['row_changes']}")

# Are row changes always zero? delta(k) = # of 3-cycles through k.
# Since lambda is preserved, c3 is preserved (c3 = S1/3).
# But individual delta(k) may change!
delta_changes = [d['row_changes'] for d in w_changes_list]
all_zero = all(all(r == 0 for r in rc) for rc in delta_changes)
print(f"  delta(k) always preserved? {all_zero}")

if not all_zero:
    # Analyze which vertices change
    for d in w_changes_list[:3]:
        print(f"\n  Subset S = {d['subset']}:")
        print(f"    Row changes: {d['row_changes']}")
        for k in range(n):
            rc = d['row_changes'][k]
            if rc != 0:
                in_S = k in d['subset']
                print(f"    Vertex {k}: delta change = {rc:+d}, in S? {in_S}")

# EXPLORATION 2: The overlap weight decomposition in the witness space
print(f"\n{'='*60}")
print("OVERLAP WEIGHTS IN WITNESS SPACE")
print(f"{'='*60}")

# Two 3-cycles C1 = {a,b,c} and C2 = {d,e,f} have overlap weight:
# W=2 if |{a,b,c} cap {d,e,f}| = 2
# W=1 if |{a,b,c} cap {d,e,f}| = 1
# W=0 if they are disjoint

# In terms of the witness matrix:
# C1 = {a,b,c} is encoded as 3 witness entries: W[c][(a,b)]=1, W[b][(a,c)]=1, W[a][(b,c)]=1
# (each vertex of the 3-cycle witnesses the other two)

# Two 3-cycles share k vertices: this corresponds to sharing k witness-row indices.
# W=2: share 2 vertices, differ in 1. The shared pair is in both cycles -> lambda overlap.
# W=1: share 1 vertex.
# W=0: completely disjoint.

# The lambda graph COLUMN for pair (u,v) counts how many 3-cycles contain both u and v.
# This is the W=2 structure! lambda(u,v) = C(lambda(u,v), 2)... no.
# lambda(u,v) = #{w: {u,v,w} is a 3-cycle vertex set}
# = #{3-cycles containing both u and v}

# The number of 3-cycle pairs with W=2:
# P2 = sum_{(u,v)} C(lambda(u,v), 2)
# This IS lambda-determined (proved in THM-171).

# The number of W=1 pairs and W=0 pairs are also lambda-determined (THM-171).
# So the OVERLAP SPECTRUM is fully controlled by lambda.

# But the overlap spectrum of c5 cycles? Of c7 cycles?
# c5 overlap is more complex since 5-vertex sets have more overlap types.

# EXPLORATION 3: The c7 non-measurability
print(f"\n{'='*60}")
print("c7 NON-MEASURABILITY: WHAT EXTRA INFORMATION IS NEEDED?")
print(f"{'='*60}")

# c7 is NOT lambda-determined. So there exist two tournaments with the same
# lambda graph but different c7 counts. What distinguishes them?

# The witness matrix W has the same column sums (lambda) but different
# individual entries. Does the witness matrix determine c7?

# For any k-cycle count: the cycle visits k vertices and uses k arcs.
# A 7-cycle uses 7 arcs from the tournament adjacency A.
# lambda captures the 3-cycle information but not higher.

# The sigma function: sigma(u,v) = common successors + common predecessors.
# sigma is NOT lambda-determined (verified above).
# sigma captures information about (A^2)[u][v] beyond what lambda gives.

# Hypothesis: c7 = f(lambda, sigma)?
# Or: c7 = f(lambda, full A^2)?

# At n=7: check if sigma determines c7 within each lambda class.
print(f"Testing at n=7: does (lambda, sigma) determine c7?")

n = 7
total_bits = n * (n-1) // 2
np.random.seed(42)

groups = defaultdict(set)
for trial in range(50000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A2 = A @ A

    # Key for lambda + sigma
    lam_key = tuple(L[i][j] for i in range(n) for j in range(n))

    # sigma(u,v) for all u<v
    sigma_vals = []
    for u in range(n):
        for v in range(u+1, n):
            sigma = n - 2 - int(A2[u][v]) - int(A2[v][u])
            sigma_vals.append(sigma)
    sigma_key = tuple(sigma_vals)

    combined_key = lam_key + sigma_key

    # c7
    c7 = 0
    for combo in combinations(range(n), 7):
        for perm in permutations(combo[1:]):
            path = (combo[0],) + perm
            valid = True
            for i in range(7):
                if A[path[i]][path[(i+1) % 7]] != 1:
                    valid = False
                    break
            if valid:
                c7 += 1

    groups[combined_key].add(c7)

ambiguous = sum(1 for v in groups.values() if len(v) > 1)
print(f"  Groups: {len(groups)}, ambiguous: {ambiguous}")
print(f"  (lambda, sigma) determines c7? {ambiguous == 0}")

# Also check: does A^2 determine c7?
groups_a2 = defaultdict(set)
for trial in range(20000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    A2 = A @ A
    a2_key = tuple(int(A2[i][j]) for i in range(n) for j in range(n))

    c7 = 0
    for combo in combinations(range(n), 7):
        for perm in permutations(combo[1:]):
            path = (combo[0],) + perm
            valid = True
            for i in range(7):
                if A[path[i]][path[(i+1) % 7]] != 1:
                    valid = False
                    break
            if valid:
                c7 += 1

    groups_a2[a2_key].add(c7)

ambiguous_a2 = sum(1 for v in groups_a2.values() if len(v) > 1)
print(f"\n  A^2 groups: {len(groups_a2)}, ambiguous: {ambiguous_a2}")
print(f"  A^2 determines c7? {ambiguous_a2 == 0}")

# EXPLORATION 4: The hierarchical dimension of tournament invariants
print(f"\n{'='*60}")
print("DIMENSION HIERARCHY OF TOURNAMENT INVARIANTS")
print(f"{'='*60}")

# Level 0: Score sequence s = (s_0, ..., s_{n-1}) — n values, sum = C(n,2)
# Level 1: Lambda graph L[u][v] — C(n,2) values
# Level 2: Full adjacency A[u][v] — C(n,2) binary values
# Level 3: Witness matrix W[k][(u,v)] — n * C(n,2) binary values

# Score determines: c3 count (YES), c5 count (NO at n=5)
# Lambda determines: c3 (YES), c5 (YES), c7 (NO at n=7)
# A determines: all cycle counts (trivially)

# The KERNEL of each projection:
# ker(A -> score) = #{tournaments with same score sequence}
# ker(A -> lambda) = #{tournaments with same lambda graph}
# These kernels are the "fibers" of the projection.

# Vitali atoms move within the ker(A -> lambda) fiber.
# They preserve EVERYTHING that lambda determines:
# c3, c5, overlap spectrum of c3, ...

# At n=5:
n = 5
total_bits = n * (n-1) // 2

# Count fiber sizes
score_fibers = defaultdict(int)
lambda_fibers = defaultdict(int)
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    scores = tuple(sorted(int(sum(A[i])) for i in range(n)))
    L = lambda_graph(A, n)
    lam_key = tuple(L[i][j] for i in range(n) for j in range(n))
    score_fibers[scores] += 1
    lambda_fibers[lam_key] += 1

print(f"\nn={n}: {1 << total_bits} tournaments")
print(f"  Score classes: {len(score_fibers)}")
print(f"  Lambda classes: {len(lambda_fibers)}")
print(f"  Score fiber sizes: {sorted(Counter(score_fibers.values()).items())}")
print(f"  Lambda fiber sizes: {sorted(Counter(lambda_fibers.values()).items())}")

# At n=6:
n = 6
total_bits = n * (n-1) // 2

score_fibers = defaultdict(int)
lambda_fibers = defaultdict(int)
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    scores = tuple(sorted(int(sum(A[i])) for i in range(n)))
    L = lambda_graph(A, n)
    lam_key = tuple(L[i][j] for i in range(n) for j in range(n))
    score_fibers[scores] += 1
    lambda_fibers[lam_key] += 1

print(f"\nn={n}: {1 << total_bits} tournaments")
print(f"  Score classes: {len(score_fibers)}")
print(f"  Lambda classes: {len(lambda_fibers)}")
print(f"  Score fiber sizes: min={min(score_fibers.values())}, max={max(score_fibers.values())}, "
      f"mean={np.mean(list(score_fibers.values())):.1f}")
print(f"  Lambda fiber sizes: min={min(lambda_fibers.values())}, max={max(lambda_fibers.values())}, "
      f"mean={np.mean(list(lambda_fibers.values())):.1f}")

# Ratio: how much finer is lambda than score?
print(f"  Lambda/score refinement: {len(lambda_fibers) / len(score_fibers):.1f}x")

# EXPLORATION 5: The Vitali analogy in detail
print(f"\n{'='*60}")
print("THE VITALI ANALOGY")
print(f"{'='*60}")
print("""
VITALI SET IN R/Q              | TOURNAMENTS & LAMBDA
================================|================================
Real line R                     | Tournament space {0,1}^C(n,2)
Rational translations Q         | Vitali atoms (lambda-preserving)
Coset R/Q                       | Lambda fiber (same lambda graph)
Lebesgue measure on R           | c5_dir (or any lambda-measurable fn)
Non-measurable set V            | c7_dir (NOT lambda-measurable)
                                |
Key property: Q preserves       | Lambda atoms preserve c5
Lebesgue measure                |  (THM-172 + THM-173)
                                |
Key property: V has no           | c7 has no lambda-formula
well-defined measure under Q     |  (5 ambiguous groups at n=7)
                                |
Sigma-algebra of measurable sets | Lambda-measurable invariants:
                                |   c3, c5, overlap spectra,
                                |   independence #s i_k for c3, c5...
                                |
Non-measurable "residual"       | Lambda-non-measurable residual:
                                |   c7, c9, ..., and their overlap
                                |   contributions to i_k
""")

print("The {2,1,0} overlap weights connect to this as follows:")
print("  W=2 pairs: controlled by lambda (THM-171)")
print("  W=1 pairs: controlled by lambda + vertex degrees (THM-171)")
print("  W=0 pairs: controlled by lambda (THM-171)")
print()
print("For c3-c3 overlaps, ALL are lambda-measurable (THM-171).")
print("For c3-c5 overlaps: count controlled by both c3 and c5,")
print("  both lambda-measurable. So c3-c5 overlap count is measurable.")
print("For c5-c5 overlaps: count may be measurable (both c5 are).")
print("For c3-c7 overlaps: count involves c7, which is NOT measurable.")
print("  This is where the 'non-measurability' first enters the overlap structure.")

print("\nDone.")
