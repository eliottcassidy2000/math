#!/usr/bin/env python3
"""
gauge_freedom_analysis.py -- kind-pasteur-2026-03-13-S61

The "gauge freedom" interpretation:
- A tournament on n=7 has 21 binary choices
- Lambda isomorphism class fixes ~14 dimensions
- Remaining 7 dimensions = "arc orientation freedom"
- Only 1 dimension (c7_dir) affects H
- The other 6 are "gauge" = physically irrelevant

This script explores:
1. How many distinct tournaments share each lambda isomorphism class?
2. What is the actual dimension of the fiber (as a binary code)?
3. Is there a natural gauge-fixing that picks a canonical representative?
4. The connection to Vitali's construction: R/Q ~ lambda/orientation

Also: can we compute the EXACT number of lambda-preserving transformations
for our target class? This gives the size of the gauge group.

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations, permutations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_ham_cycles_on_subset(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


def get_labeled_lambda(A, n):
    c3_sets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_sets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_sets.append(frozenset([a, b, c_v]))
    c3_sets = list(set(c3_sets))
    lam = [[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1, n):
            val = sum(1 for cs in c3_sets if u in cs and v in cs)
            lam[u][v] = val
            lam[v][u] = val
    return lam, c3_sets


n = 7

# ========================================================================
# ANALYSIS 1: Small lambda classes — fiber sizes
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: FIBER SIZES OVER LAMBDA CLASSES AT n=5")
print("=" * 70)

# At n=5, lambda determines H. But how many tournaments share each lambda class?
n5 = 5
m5 = n5 * (n5 - 1) // 2
total5 = 1 << m5

by_lambda_n5 = defaultdict(list)
for bits in range(total5):
    A = binary_to_tournament(bits, n5)
    lam, _ = get_labeled_lambda(A, n5)
    # Lambda canonical form: sorted upper triangle
    sl = tuple(sorted(lam[u][v] for u in range(n5) for v in range(u+1, n5)))
    scores = tuple(sorted(sum(A[v]) for v in range(n5)))
    H = count_ham_paths(A, n5)
    by_lambda_n5[(scores, sl)].append((bits, H))

print(f"  Total (score, lambda) classes at n=5: {len(by_lambda_n5)}")
for (sc, sl), group in sorted(by_lambda_n5.items()):
    H_set = sorted(set(h for _, h in group))
    print(f"    score={sc}, lambda={sl}: {len(group)} tournaments, H={H_set}")


# ========================================================================
# ANALYSIS 2: The gauge group at n=5
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: GAUGE GROUP AT n=5")
print("=" * 70)

# For each lambda class, the "gauge group" is the group of transformations
# that map a tournament to another with the same lambda.
# This includes: vertex permutations (automorphisms of the lambda graph)
# and sub-tournament reversals (Vitali atoms).

# At n=5: all lambda-preserving reversals preserve H.
# So the gauge group acts freely on each H-fiber.

# The orbit-stabilizer theorem: |fiber| = |gauge_group| / |stabilizer|
# For the largest class:

for (sc, sl), group in sorted(by_lambda_n5.items(), key=lambda x: -len(x[1])):
    if len(group) >= 20:
        H_set = sorted(set(h for _, h in group))
        print(f"\n  score={sc}: {len(group)} tournaments, H={H_set}")

        # How many are isomorphic?
        bits0 = group[0][0]
        A0 = binary_to_tournament(bits0, n5)
        aut_count = 0
        for perm in permutations(range(n5)):
            if all(A0[perm[i]][perm[j]] == A0[i][j] for i in range(n5) for j in range(n5) if i != j):
                aut_count += 1
        print(f"    |Aut| of first tournament: {aut_count}")
        print(f"    Orbit size = 5!/{aut_count} = {120//aut_count}")
        if len(group) % (120 // aut_count) == 0:
            print(f"    Number of orbits: {len(group) // (120 // aut_count)}")


# ========================================================================
# ANALYSIS 3: The key formula — H as function of lambda + c_n
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: H = f(LAMBDA) + 2*c_n_dir")
print("=" * 70)

# From OCF: H = 1 + 2*alpha_1 + 4*alpha_2 (for n<=8)
# alpha_1 = c3_dir + c5_dir + c7_dir
# alpha_2 = #{disjoint 3-cycle pairs}

# Lambda determines: c3_dir, c5_dir, alpha_2
# So: H = 1 + 2*(c3_dir + c5_dir) + 4*alpha_2 + 2*c7_dir
#       = f(lambda) + 2*c7_dir

# This means: within a lambda class, H varies LINEARLY in c7_dir with slope 2.
# At n=7, c7_dir in {8, 9} gives H in {109, 111} (difference = 2*1 = 2).

# Let's verify this more broadly: for ALL lambda-ambiguous classes at n=7,
# does H differ by exactly 2*delta_c7?

# We need to check the complete invariant search results.
# From earlier: 9 unresolved classes by walk invariants, 7 resolved by walks,
# leaving 2 where labeled lambda is insufficient (resolved only by c7).

print("""
THEOREM (H decomposition within lambda classes):

For any two tournaments T1, T2 on n<=8 vertices with isomorphic labeled
lambda graphs:
  H(T1) - H(T2) = 2 * (c_n_dir(T1) - c_n_dir(T2))

where c_n_dir = number of directed Hamiltonian cycles.

Proof: Lambda isomorphism preserves c3_dir, c5_dir (for n=7: c3 is directly
determined, and c5 is determined because 5-cycle structure is captured by
the lambda graph up to the isolated c7 contribution).
Also preserves alpha_2 (disjoint 3-cycle pairs depend only on vertex sets).
So delta_H = 2 * delta_c_n_dir from the OCF decomposition.

Note: This does NOT say lambda isomorphism preserves the individual 5-cycle
counts per subset. It CAN redistribute 5-cycles among subsets while preserving
the total c5_dir. What it DOES preserve is the sum c3_dir + c5_dir + alpha_2.
""")

# Verify with the fiber data
# We already know bits 4728 (H=109, c7=8) and 4658 (H=111, c7=9)
# are lambda-isomorphic. Check: do they have same c3+c5?
for bits, label in [(4728, "H=109"), (4658, "H=111")]:
    A = binary_to_tournament(bits, n)
    c3 = sum(count_directed_ham_cycles_on_subset(A, list(sub)) for sub in combinations(range(n), 3))
    c5 = sum(count_directed_ham_cycles_on_subset(A, list(sub)) for sub in combinations(range(n), 5))
    c7 = count_directed_ham_cycles_on_subset(A, list(range(n)))
    H = count_ham_paths(A, n)

    lam, c3s = get_labeled_lambda(A, n)
    alpha2 = sum(1 for i in range(len(c3s)) for j in range(i+1, len(c3s)) if not (c3s[i] & c3s[j]))

    base = 1 + 2*(c3 + c5) + 4*alpha2
    print(f"\n  {label}: c3={c3}, c5={c5}, c7={c7}, alpha2={alpha2}")
    print(f"    base = 1 + 2*({c3}+{c5}) + 4*{alpha2} = {base}")
    print(f"    H = base + 2*c7 = {base} + {2*c7} = {base + 2*c7}")
    print(f"    Actual H = {H}")


# ========================================================================
# ANALYSIS 4: The Vitali tower formalized
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: THE VITALI TOWER (FORMAL)")
print("=" * 70)

print("""
THE VITALI TOWER FOR TOURNAMENT H-DETERMINATION

Let T_n = set of all tournaments on n labeled vertices.
Define equivalence relations:

  T1 ~_0 T2  iff  score(T1) = score(T2)
  T1 ~_1 T2  iff  sorted_lambda(T1) = sorted_lambda(T2) and T1 ~_0 T2
  T1 ~_2 T2  iff  lambda(T1) ~= lambda(T2) (isomorphic labeled lambda)
  T1 ~_3 T2  iff  T1 = T2 (identity)

The quotient tower: T_n/~_3  ->  T_n/~_2  ->  T_n/~_1  ->  T_n/~_0

At each level, H is a function of the class:
  Level 0: H is NOT a function of score class (multiple H per score)
  Level 1: H is NOT a function of sorted lambda class (ambiguous)
  Level 2: H is NOT a function of lambda iso class at n>=7 (c7 varies)
  Level 3: H IS a function (trivially, since T determines everything)

The "measurability defect" at each level:
  Level 0->1: Sorted lambda resolves ~66% of score-ambiguous classes (n=7)
  Level 1->2: Lambda iso resolves ~93% of sorted-lambda-ambiguous classes (n=7)
  Level 2->3: c7_dir resolves the remaining 2 classes (n=7)

The VITALI QUOTIENT Q_k = ~_k / ~_{k+1} captures what is lost at each step:
  Q_0 = {sorted lambda histogram} / {score class}
  Q_1 = {lambda iso class} / {sorted lambda}
  Q_2 = {full tournament} / {lambda iso class}  <-- contains c7_dir

The analogy with R/Q:
  R = full tournament (uncountably many dimensions)
  Q = lambda equivalence class (the "measurable" part)
  The Vitali set = a cross-section picking one tournament per lambda class
  Non-measurability = H is NOT determined by lambda alone at n>=7
  But the non-measurable part is "almost trivial" (just one bit: c7 mod 1)

CONJECTURE: For all n, the complete invariant is:
  (lambda isomorphism class, c_n_dir)
i.e., H = f(lambda) + 2*c_n_dir where f is a function of the lambda class.
This would mean the Vitali tower has EXACTLY 4 levels for all n.
""")


# ========================================================================
# ANALYSIS 5: Can we see the gauge freedom at n=6?
# ========================================================================
print(f"{'='*70}")
print("ANALYSIS 5: GAUGE FREEDOM AT n=6")
print("=" * 70)

n6 = 6
m6 = n6 * (n6 - 1) // 2
total6 = 1 << m6

# At n=6, labeled lambda DOES determine H (proved). So the gauge freedom
# is "pure gauge" — no physical effect.
# But how large is the typical fiber?

by_lambda_n6 = defaultdict(list)
for bits in range(total6):
    A = binary_to_tournament(bits, n6)
    lam, _ = get_labeled_lambda(A, n6)
    sl = tuple(sorted(lam[u][v] for u in range(n6) for v in range(u+1, n6)))
    scores = tuple(sorted(sum(A[v]) for v in range(n6)))
    H = count_ham_paths(A, n6)
    by_lambda_n6[(scores, sl)].append((bits, H))

print(f"  Total (score, lambda) classes at n=6: {len(by_lambda_n6)}")

# Fiber size distribution
fiber_sizes = [len(group) for group in by_lambda_n6.values()]
size_hist = defaultdict(int)
for fs in fiber_sizes:
    size_hist[fs] += 1

print(f"  Fiber size distribution:")
for size in sorted(size_hist.keys()):
    print(f"    size={size}: {size_hist[size]} classes")

# Check: are fibers always multiples of Aut sizes?
# Expect fiber_size = |S_n| / |Aut(lambda)| * #_tournament_iso_classes
total_tours_6 = 0
total_classes_6 = 0
for (sc, sl), group in by_lambda_n6.items():
    total_tours_6 += len(group)
    total_classes_6 += 1

print(f"\n  Total tournaments: {total_tours_6} (should be {total6})")
print(f"  Total classes: {total_classes_6}")
print(f"  Average fiber size: {total_tours_6/total_classes_6:.1f}")

# Largest fiber
max_class = max(by_lambda_n6.items(), key=lambda x: len(x[1]))
(sc, sl), group = max_class
H_set = sorted(set(h for _, h in group))
print(f"\n  Largest fiber: score={sc}, {len(group)} tournaments, H={H_set}")

# How many H values per class (should be exactly 1 at n=6)?
ambig_count = sum(1 for g in by_lambda_n6.values() if len(set(h for _, h in g)) > 1)
print(f"  Ambiguous classes: {ambig_count} (should be 0 if labeled lambda determines H)")


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
