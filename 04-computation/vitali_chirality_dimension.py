#!/usr/bin/env python3
"""
vitali_chirality_dimension.py -- kind-pasteur-2026-03-13-S61

The Vitali atom preserves lambda but changes c7. The 3-cycle orientation
signature changes dramatically (sum 5 vs sum 1). This suggests a
CHIRALITY structure.

Questions:
1. Is the 3-cycle orientation sum an invariant of the labeled lambda class?
   (No -- it changed! But what IS invariant about it?)
2. The 6-arc flip reverses a 4-vertex sub-tournament. In the complement,
   clockwise cycles become counter-clockwise. Is there a parity invariant?
3. The oriented lambda graph: each edge of lambda gets an orientation from
   the tournament. Is the "cycle space" of this oriented graph related to c7?
4. Hidden dimension: what is the dimension of the fiber over a lambda class?
   (How many non-isomorphic tournaments share a lambda graph?)

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


def lambda_canonical_form(lam, n):
    """Canonical form of lambda matrix under vertex permutation."""
    best = None
    for perm in permutations(range(n)):
        form = tuple(lam[perm[i]][perm[j]] for i in range(n) for j in range(i+1, n))
        if best is None or form < best:
            best = form
    return best


n = 7

# ========================================================================
# ANALYSIS 1: Fiber dimension over lambda classes
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: FIBER DIMENSION OVER LAMBDA CLASSES")
print("=" * 70)

# For the specific lambda class of our deepest ambiguity,
# how many non-isomorphic tournaments share this lambda graph?

# First: find ALL tournaments with this lambda canonical form.
# This is expensive at n=7 (2M tournaments). Let's focus on the specific score class.
# Score (2,2,3,3,3,3,5)

target_score = (2, 2, 3, 3, 3, 3, 5)
target_lam_canonical = None

A_target = binary_to_tournament(4658, n)
lam_target, _ = get_labeled_lambda(A_target, n)
# Compute canonical form just for this one
# (Full canonical is expensive; use sorted lambda histogram instead)
sorted_lam_target = tuple(sorted(lam_target[u][v] for u in range(n) for v in range(u+1, n)))

print(f"  Target lambda histogram: {sorted_lam_target}")

# Find all tournaments in this score class with this lambda histogram
matching = []
m = n * (n - 1) // 2
total = 1 << m

count = 0
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted(sum(A[v]) for v in range(n)))
    if scores != target_score:
        continue

    lam, _ = get_labeled_lambda(A, n)
    sl = tuple(sorted(lam[u][v] for u in range(n) for v in range(u+1, n)))
    if sl != sorted_lam_target:
        continue

    H = count_ham_paths(A, n)
    c7 = count_directed_ham_cycles_on_subset(A, list(range(n)))
    matching.append((bits, H, c7))
    count += 1

print(f"\n  Tournaments with score {target_score} and lambda histogram {sorted_lam_target}:")
print(f"  Total: {len(matching)}")

# Group by H
by_H = defaultdict(list)
for bits, H, c7 in matching:
    by_H[H].append((bits, c7))

for H in sorted(by_H.keys()):
    c7_set = set(c7 for _, c7 in by_H[H])
    print(f"    H={H}: {len(by_H[H])} tournaments, c7 in {sorted(c7_set)}")


# ========================================================================
# ANALYSIS 2: Lambda isomorphism classes within the fiber
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: LAMBDA ISOMORPHISM CLASSES IN FIBER")
print("=" * 70)

# Among the matching tournaments, compute lambda degree sequence
# as a cheaper proxy for lambda isomorphism

by_deg_seq = defaultdict(list)
for bits, H, c7 in matching:
    A = binary_to_tournament(bits, n)
    lam, _ = get_labeled_lambda(A, n)
    deg_seq = tuple(sorted(sum(lam[v]) for v in range(n)))
    by_deg_seq[deg_seq].append((bits, H, c7))

print(f"  Lambda degree sequence classes: {len(by_deg_seq)}")
for ds in sorted(by_deg_seq.keys()):
    group = by_deg_seq[ds]
    H_set = sorted(set(H for _, H, _ in group))
    c7_set = sorted(set(c7 for _, _, c7 in group))
    print(f"    deg_seq={ds}: {len(group)} tournaments, H in {H_set}, c7 in {c7_set}")


# ========================================================================
# ANALYSIS 3: The oriented lambda graph
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: THE ORIENTED LAMBDA GRAPH")
print("=" * 70)

# For each tournament, define the "oriented lambda graph":
# - Same vertices as the tournament (0..n-1)
# - For each pair (u,v) with lambda[u][v] > 0, put a directed edge
#   with weight lambda[u][v] in the tournament direction
# This is the lambda graph with tournament orientation.

# The key question: what properties of this oriented graph are lambda-invariant
# (same for all tournaments with same lambda graph) vs H-informative (determine c7)?

for bits, expected_H in [(4658, 111), (4728, 109)]:
    A = binary_to_tournament(bits, n)
    lam, c3s = get_labeled_lambda(A, n)
    H = count_ham_paths(A, n)

    print(f"\n  bits={bits}, H={H}:")

    # Oriented lambda: for each edge (u,v) with lam>0, record direction and weight
    oriented = []
    for u in range(n):
        for v in range(u+1, n):
            if lam[u][v] > 0:
                if A[u][v]:
                    oriented.append((u, v, lam[u][v], "->"))
                else:
                    oriented.append((u, v, lam[u][v], "<-"))

    # Count orientations by lambda value
    fwd_by_lam = defaultdict(int)
    bwd_by_lam = defaultdict(int)
    for u, v, w, d in oriented:
        if d == "->":
            fwd_by_lam[w] += 1
        else:
            bwd_by_lam[w] += 1

    for lv in sorted(set(w for _, _, w, _ in oriented)):
        total = fwd_by_lam.get(lv, 0) + bwd_by_lam.get(lv, 0)
        print(f"    lambda={lv}: {fwd_by_lam.get(lv,0)} fwd, {bwd_by_lam.get(lv,0)} bwd (total {total})")

    # For each vertex, compute "oriented lambda degree": sum of lambda*direction
    # where direction is +1 (out) or -1 (in)
    for v in range(n):
        out_sum = sum(lam[v][w] for w in range(n) if w != v and A[v][w] and lam[v][w] > 0)
        in_sum = sum(lam[v][w] for w in range(n) if w != v and A[w][v] and lam[v][w] > 0)
        # print(f"    v={v}: oriented lambda out={out_sum}, in={in_sum}, net={out_sum-in_sum}")

    # Net orientation per vertex
    net = []
    for v in range(n):
        out_s = sum(lam[v][w] for w in range(n) if w != v and A[v][w])
        in_s = sum(lam[v][w] for w in range(n) if w != v and A[w][v])
        net.append(out_s - in_s)
    print(f"    Net oriented lambda per vertex: {net}")
    print(f"    Sorted: {sorted(net)}")

    # Total "oriented lambda flow": sum of lambda[u][v] * sign(A[u][v] - 0.5)
    total_flow = sum(lam[u][v] * (1 if A[u][v] else -1) for u in range(n) for v in range(u+1, n) if lam[u][v] > 0)
    print(f"    Total oriented lambda flow: {total_flow}")


# ========================================================================
# ANALYSIS 4: Parity invariant of the 4-vertex reversal
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: CHIRALITY AND PARITY")
print("=" * 70)

# A tournament on n vertices has n*(n-1)/2 arcs. At n=4, this is 6 arcs.
# A 4-vertex tournament on {a,b,c,d} can be transitive (score 0,1,2,3)
# or have one 3-cycle (score 1,1,2,2).

# The 1,1,2,2 type has EXACTLY one 3-cycle and one "source-sink" pair.
# When we reverse the sub-tournament, the 3-cycle reverses.

# Key observation: the 4-vertex sub-tournament being reversed has
# scores (1,1,2,2) in BOTH cases. This means it has one 3-cycle.
# Reversing changes the 3-cycle direction but preserves that it HAS one.

# The 3-cycle orientation signature (sum=5 vs sum=1) is NOT preserved.
# But note: sum=5 + sum=1 = 6, and there are 11 3-cycle vertex sets total.
# Each vertex set has one 3-cycle vertex set, so 11 3-cycles.
# In the CW/CCW count, we saw:
# H=109: 5 clockwise + 6 counter-clockwise = 11
# H=111: 1 clockwise + 10 counter-clockwise = 11

# So the "chirality" = #CW - #CCW differs: 5-6=-1 vs 1-10=-9.
# The reversal changes chirality by... let me check.

# Actually, the 3-cycle orientation is defined relative to vertex ordering.
# Let me use a more canonical notion: the SIGN of the 3-cycle.
# For a 3-cycle on {a,b,c} with a<b<c, sign = +1 if a->b->c->a, -1 if a->c->b->a.

print("  3-cycle chirality (sign = +1 if a->b->c->a for a<b<c):")
for bits, label in [(4658, "H=111"), (4728, "H=109")]:
    A = binary_to_tournament(bits, n)
    lam, c3s = get_labeled_lambda(A, n)

    total_sign = 0
    for c3 in c3s:
        a, b, c = sorted(c3)
        if A[a][b] and A[b][c] and A[c][a]:
            total_sign += 1
        else:
            total_sign -= 1
    print(f"    {label} (bits={bits}): total chirality = {total_sign} ({len(c3s)} cycles)")


# ========================================================================
# ANALYSIS 5: The connection to c7 via winding number
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: WINDING NUMBER AND HAMILTONIAN CYCLES")
print("=" * 70)

# A Hamiltonian cycle visits all 7 vertices. It creates a directed 7-cycle.
# The "winding" of this cycle through the 3-cycle structure might connect
# to the chirality.

# For each directed Hamiltonian cycle, count how many 3-cycles it "agrees with"
# (i.e., the 3 consecutive vertices in the HC form a sub-arc agreeing with a 3-cycle direction)

for bits, label in [(4658, "H=111"), (4728, "H=109")]:
    A = binary_to_tournament(bits, n)
    c7 = count_directed_ham_cycles_on_subset(A, list(range(n)))
    print(f"\n  {label} (bits={bits}): c7={c7} directed Hamiltonian cycles")

    # Find all directed Hamiltonian cycles
    # A directed HC is a permutation (v0,v1,...,v_{n-1}) with
    # v0->v1->...->v_{n-1}->v0
    hc_count = 0
    for perm in permutations(range(n)):
        if perm[0] != 0:  # Fix first vertex to avoid counting rotations
            continue
        is_hc = all(A[perm[i]][perm[(i+1) % n]] for i in range(n))
        if is_hc:
            hc_count += 1
            # For this HC, how many consecutive triples are 3-cycle vertex sets?
            agrees = 0
            for i in range(n):
                triple = frozenset([perm[i], perm[(i+1)%n], perm[(i+2)%n]])
                a, b, c = sorted(triple)
                if A[a][b] and A[b][c] and A[c][a]:
                    agrees += 1
                elif A[a][c] and A[c][b] and A[b][a]:
                    agrees += 1
            print(f"    HC: {perm} -> {agrees}/7 consecutive triples are 3-cycles")

    # Count = c7/2 since we fix first vertex but not direction
    # Actually: fixing v0=0 and counting all HCs gives c7 * 2 / n... hmm.
    # A directed HC has n rotations, so fixing v0=0 gives c7 directed HCs / n
    # Wait no: a directed HC is a cyclic sequence, fixing first = dividing by n.
    # So we found c7/n * 2 = ... let me just verify
    print(f"    Found {hc_count} directed HCs (with v0=0)")
    print(f"    Total directed HCs = {hc_count} * {n} / 1 = ... hmm")
    # Actually fixing v0=0 in a DIRECTED cycle gives exactly c7_dir cycles
    # because each cycle visits vertex 0 exactly once, and we check all orderings starting from 0
    # Wait: no. Fixing v0=0, the remaining n-1 positions give (n-1)! permutations.
    # Each directed Hamiltonian cycle appears exactly 1 time (shifted to start at 0).
    print(f"    (c7_dir should be {c7}, found {hc_count})")


# ========================================================================
# ANALYSIS 6: Hidden dimension = oriented cycles in lambda
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 6: ORIENTED CYCLE SPACE OF LAMBDA GRAPH")
print("=" * 70)

# The lambda graph on n=7 has C(7,2) = 21 possible edges.
# With lambda values > 0, it's a weighted graph.
# The tournament orientation makes it a DIRECTED weighted graph.
# The cycle space of the undirected lambda graph has dimension m - n + c
# where m = #edges, n = #vertices, c = #components.

# For our target: all pairs have lambda >= 0, with some = 0.

for bits, label in [(4658, "H=111"), (4728, "H=109")]:
    A = binary_to_tournament(bits, n)
    lam, c3s = get_labeled_lambda(A, n)

    edges = sum(1 for u in range(n) for v in range(u+1, n) if lam[u][v] > 0)
    zero_pairs = sum(1 for u in range(n) for v in range(u+1, n) if lam[u][v] == 0)

    print(f"\n  {label} (bits={bits}):")
    print(f"    Lambda edges (>0): {edges}, zero pairs: {zero_pairs}")
    print(f"    Cycle space dimension: {edges} - {n} + 1 = {edges - n + 1}")
    # (Assuming lambda graph is connected)

    # Check connectivity
    visited = set()
    stack = [0]
    while stack:
        v = stack.pop()
        if v in visited:
            continue
        visited.add(v)
        for w in range(n):
            if w != v and lam[v][w] > 0 and w not in visited:
                stack.append(w)
    print(f"    Connected: {len(visited) == n}")


# ========================================================================
# ANALYSIS 7: The {2,1,0} Vitali structure at n=5
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 7: VITALI STRUCTURE AT n=5")
print("=" * 70)

# At n=5, lambda determines H completely. So the Vitali quotient is trivial.
# But does the 4-vertex reversal preserve lambda at n=5 too?
# At n=5, a 4-vertex sub-tournament reversal affects C(4,2)=6 arcs out of C(5,2)=10.
# That's a lot of the structure!

n5 = 5
m5 = n5 * (n5 - 1) // 2
total5 = 1 << m5

reversals_n5 = 0
preserving_n5 = 0
nonzero_delta_n5 = 0

for bits in range(total5):
    A = binary_to_tournament(bits, n5)
    lam, _ = get_labeled_lambda(A, n5)
    H = count_ham_paths(A, n5)

    for S in combinations(range(n5), 4):
        A_rev = [row[:] for row in A]
        for i in S:
            for j in S:
                if i != j:
                    A_rev[i][j] = 1 - A[i][j]
        lam_rev, _ = get_labeled_lambda(A_rev, n5)
        reversals_n5 += 1
        if all(lam_rev[i][j] == lam[i][j] for i in range(n5) for j in range(i+1, n5)):
            preserving_n5 += 1
            H_rev = count_ham_paths(A_rev, n5)
            if H_rev != H:
                nonzero_delta_n5 += 1

print(f"  n=5: {total5} tournaments x C(5,4)=5 subsets = {reversals_n5} total reversals")
print(f"  Lambda-preserving: {preserving_n5} ({100*preserving_n5/reversals_n5:.1f}%)")
print(f"  Lambda-preserving with delta_H != 0: {nonzero_delta_n5}")

if nonzero_delta_n5 == 0:
    print("  => At n=5, lambda-preserving reversals ALWAYS preserve H!")
    print("     Consistent with labeled lambda being complete at n=5.")


# Also check 3-vertex reversals at n=5
reversals_3_n5 = 0
preserving_3_n5 = 0
nonzero_delta_3_n5 = 0

for bits in range(total5):
    A = binary_to_tournament(bits, n5)
    lam, _ = get_labeled_lambda(A, n5)
    H = count_ham_paths(A, n5)

    for S in combinations(range(n5), 3):
        A_rev = [row[:] for row in A]
        for i in S:
            for j in S:
                if i != j:
                    A_rev[i][j] = 1 - A[i][j]
        lam_rev, _ = get_labeled_lambda(A_rev, n5)
        reversals_3_n5 += 1
        if all(lam_rev[i][j] == lam[i][j] for i in range(n5) for j in range(i+1, n5)):
            preserving_3_n5 += 1
            H_rev = count_ham_paths(A_rev, n5)
            if H_rev != H:
                nonzero_delta_3_n5 += 1

print(f"\n  3-vertex reversals at n=5:")
print(f"  Lambda-preserving: {preserving_3_n5}")
print(f"  Lambda-preserving with delta_H != 0: {nonzero_delta_3_n5}")


# ========================================================================
# ANALYSIS 8: The {2,1,0} at n=6
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 8: VITALI STRUCTURE AT n=6")
print("=" * 70)

n6 = 6
m6 = n6 * (n6 - 1) // 2
total6 = 1 << m6

reversals_n6 = 0
preserving_n6 = 0
nonzero_delta_n6 = 0
delta_H_dist_n6 = defaultdict(int)

for bits in range(total6):
    A = binary_to_tournament(bits, n6)
    lam, _ = get_labeled_lambda(A, n6)
    H = count_ham_paths(A, n6)

    for S in combinations(range(n6), 4):
        A_rev = [row[:] for row in A]
        for i in S:
            for j in S:
                if i != j:
                    A_rev[i][j] = 1 - A[i][j]
        lam_rev, _ = get_labeled_lambda(A_rev, n6)
        reversals_n6 += 1
        if all(lam_rev[i][j] == lam[i][j] for i in range(n6) for j in range(i+1, n6)):
            preserving_n6 += 1
            H_rev = count_ham_paths(A_rev, n6)
            delta = H_rev - H
            delta_H_dist_n6[delta] += 1
            if delta != 0:
                nonzero_delta_n6 += 1

print(f"  n=6: {total6} tournaments x C(6,4)=15 subsets")
print(f"  Lambda-preserving: {preserving_n6}")
print(f"  Lambda-preserving with delta_H != 0: {nonzero_delta_n6}")
print(f"  Delta H distribution:")
for dh in sorted(delta_H_dist_n6.keys()):
    print(f"    delta_H = {dh:+d}: {delta_H_dist_n6[dh]}")


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
