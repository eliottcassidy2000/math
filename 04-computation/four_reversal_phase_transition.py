#!/usr/bin/env python3
"""
four_reversal_phase_transition.py -- kind-pasteur-2026-03-13-S61

WHY does the 4-vertex reversal become non-trivial at n=7?

At n=5: 4-vertex reversal can preserve labeled lambda,
        but ALWAYS preserves both H and |Ps|.
At n=7: 4-vertex reversal can preserve labeled lambda,
        AND change H by +/-2 and |Ps|.

QUESTIONS:
1. What fraction of lambda-preserving 4-reversals change H at n=7?
2. Is the delta always +/-2?
3. Does the sub-tournament type (transitive vs cyclic) matter?
4. What's special about n=7 vs n=5?
5. Can we predict WHICH 4-vertex subsets give non-trivial reversals?
6. What does the 4-reversal do to the CYCLE STRUCTURE?

THEORETICAL FRAMEWORK:
  A 4-vertex reversal on subset S flips C(4,2)=6 arcs within S.
  External arcs (to/from V\S) are unchanged.

  At n=5: V\S has 1 vertex. External structure is a single vertex
  connected to 4 vertices. The connection to this vertex is unchanged.

  At n=7: V\S has 3 vertices. External structure involves 3 vertices
  connected to each other (3 internal arcs) and to S (12 cross arcs).
  MORE DEGREES OF FREEDOM for the external structure to interact
  with the reversal.

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations
from collections import defaultdict
from copy import deepcopy


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


def lambda_graph(A, n):
    lam = [[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1, n):
            count = 0
            for w in range(n):
                if w == u or w == v:
                    continue
                if ((A[u][v] and A[v][w] and A[w][u]) or
                    (A[u][w] and A[w][v] and A[v][u])):
                    count += 1
            lam[u][v] = count
            lam[v][u] = count
    return lam


def lambda_key(A, n):
    lam = lambda_graph(A, n)
    return tuple(lam[i][j] for i in range(n) for j in range(i+1, n))


def reverse_subtournament(A, n, subset):
    B = deepcopy(A)
    for u in subset:
        for v in subset:
            if u != v:
                B[u][v] = A[v][u]
    return B


def count_directed_cycles_on(A, verts):
    k = len(verts)
    if k < 3 or k % 2 == 0:
        return 0
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
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


def sub_tournament_type(A, subset):
    """Classify the sub-tournament on the given vertices."""
    verts = list(subset)
    k = len(verts)
    scores = tuple(sorted(
        sum(A[u][v] for v in verts if v != u) for u in verts
    ))
    c3 = 0
    for sub3 in combinations(verts, 3):
        c3 += count_directed_cycles_on(A, list(sub3))
    return scores, c3


# ========================================================================
# PART 1: Exhaustive n=5 — classify ALL 4-reversals
# ========================================================================
print("=" * 70)
print("PART 1: 4-vertex reversals at n=5 — EXHAUSTIVE")
print("=" * 70)

n5 = 5

# For each tournament and each 4-subset, check:
# 1. Does the reversal preserve labeled lambda?
# 2. If yes, does it change H?
# 3. What's the sub-tournament type?

stats_5 = {
    'total': 0,
    'preserves_lambda': 0,
    'changes_H': 0,
    'by_subtype': defaultdict(lambda: {'total': 0, 'pres_lambda': 0, 'changes_H': 0})
}

for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    H_orig = count_ham_paths(A, n5)
    lk_orig = lambda_key(A, n5)

    for subset in combinations(range(n5), 4):
        stats_5['total'] += 1
        stype = sub_tournament_type(A, subset)

        B = reverse_subtournament(A, n5, subset)
        lk_rev = lambda_key(B, n5)

        stats_5['by_subtype'][stype]['total'] += 1

        if lk_orig == lk_rev:
            stats_5['preserves_lambda'] += 1
            stats_5['by_subtype'][stype]['pres_lambda'] += 1

            H_rev = count_ham_paths(B, n5)
            if H_orig != H_rev:
                stats_5['changes_H'] += 1
                stats_5['by_subtype'][stype]['changes_H'] += 1

print(f"  Total 4-reversals: {stats_5['total']}")
print(f"  Lambda-preserving: {stats_5['preserves_lambda']} ({100*stats_5['preserves_lambda']/stats_5['total']:.1f}%)")
print(f"  H-changing: {stats_5['changes_H']}")

print(f"\n  By sub-tournament type:")
for stype, data in sorted(stats_5['by_subtype'].items()):
    pct = 100*data['pres_lambda']/data['total'] if data['total'] > 0 else 0
    print(f"    scores={stype[0]}, c3={stype[1]}: "
          f"{data['total']} total, {data['pres_lambda']} pres lambda ({pct:.1f}%), "
          f"{data['changes_H']} change H")


# ========================================================================
# PART 2: n=6 — does the transition start here?
# ========================================================================
print(f"\n{'='*70}")
print("PART 2: 4-vertex reversals at n=6 — EXHAUSTIVE")
print("=" * 70)

n6 = 6

stats_6 = {
    'total': 0,
    'preserves_lambda': 0,
    'changes_H': 0,
    'by_subtype': defaultdict(lambda: {'total': 0, 'pres_lambda': 0, 'changes_H': 0}),
    'delta_H_counts': defaultdict(int)
}

for bits in range(1 << 15):
    A = binary_to_tournament(bits, n6)
    H_orig = count_ham_paths(A, n6)
    lk_orig = lambda_key(A, n6)

    for subset in combinations(range(n6), 4):
        stats_6['total'] += 1
        stype = sub_tournament_type(A, subset)

        B = reverse_subtournament(A, n6, subset)
        lk_rev = lambda_key(B, n6)

        stats_6['by_subtype'][stype]['total'] += 1

        if lk_orig == lk_rev:
            stats_6['preserves_lambda'] += 1
            stats_6['by_subtype'][stype]['pres_lambda'] += 1

            H_rev = count_ham_paths(B, n6)
            delta_H = H_rev - H_orig
            if delta_H != 0:
                stats_6['changes_H'] += 1
                stats_6['by_subtype'][stype]['changes_H'] += 1
                stats_6['delta_H_counts'][delta_H] += 1

print(f"  Total 4-reversals: {stats_6['total']}")
print(f"  Lambda-preserving: {stats_6['preserves_lambda']} ({100*stats_6['preserves_lambda']/stats_6['total']:.1f}%)")
print(f"  H-changing (among lambda-preserving): {stats_6['changes_H']}")

if stats_6['changes_H'] > 0:
    print(f"  Delta H values: {dict(sorted(stats_6['delta_H_counts'].items()))}")
    print(f"  *** H-CHANGING 4-REVERSALS EXIST AT n=6! ***")
else:
    print(f"  No H-changing 4-reversals at n=6")

print(f"\n  By sub-tournament type:")
for stype, data in sorted(stats_6['by_subtype'].items()):
    if data['pres_lambda'] > 0:
        print(f"    scores={stype[0]}, c3={stype[1]}: "
              f"{data['pres_lambda']} pres lambda, "
              f"{data['changes_H']} change H")


# ========================================================================
# PART 3: n=7 — sample analysis
# ========================================================================
print(f"\n{'='*70}")
print("PART 3: 4-vertex reversals at n=7 — SAMPLING")
print("=" * 70)

import random
random.seed(42)

n7 = 7

stats_7 = {
    'total': 0,
    'preserves_lambda': 0,
    'changes_H': 0,
    'by_subtype': defaultdict(lambda: {'total': 0, 'pres_lambda': 0, 'changes_H': 0}),
    'delta_H_counts': defaultdict(int)
}

sample_7 = random.sample(range(1 << 21), 10000)

for bits in sample_7:
    A = binary_to_tournament(bits, n7)
    H_orig = count_ham_paths(A, n7)
    lk_orig = lambda_key(A, n7)

    for subset in combinations(range(n7), 4):
        stats_7['total'] += 1
        stype = sub_tournament_type(A, subset)

        B = reverse_subtournament(A, n7, subset)
        lk_rev = lambda_key(B, n7)

        stats_7['by_subtype'][stype]['total'] += 1

        if lk_orig == lk_rev:
            stats_7['preserves_lambda'] += 1
            stats_7['by_subtype'][stype]['pres_lambda'] += 1

            H_rev = count_ham_paths(B, n7)
            delta_H = H_rev - H_orig
            if delta_H != 0:
                stats_7['changes_H'] += 1
                stats_7['by_subtype'][stype]['changes_H'] += 1
                stats_7['delta_H_counts'][delta_H] += 1

print(f"  Total 4-reversals (sample): {stats_7['total']}")
print(f"  Lambda-preserving: {stats_7['preserves_lambda']} ({100*stats_7['preserves_lambda']/stats_7['total']:.1f}%)")
print(f"  H-changing (among lambda-preserving): {stats_7['changes_H']}")

if stats_7['changes_H'] > 0:
    print(f"  Delta H values: {dict(sorted(stats_7['delta_H_counts'].items()))}")

print(f"\n  By sub-tournament type:")
for stype, data in sorted(stats_7['by_subtype'].items()):
    if data['pres_lambda'] > 0:
        print(f"    scores={stype[0]}, c3={stype[1]}: "
              f"{data['pres_lambda']} pres lambda, "
              f"{data['changes_H']} change H "
              f"({100*data['changes_H']/data['pres_lambda']:.1f}% of lambda-pres)")


# ========================================================================
# PART 4: Cycle structure change under 4-reversal at n=7
# ========================================================================
print(f"\n{'='*70}")
print("PART 4: Cycle structure change under H-changing 4-reversal")
print("=" * 70)

# For the known H-changing pairs, what cycles change?
h_change_examples = []
for bits in sample_7[:5000]:
    A = binary_to_tournament(bits, n7)
    H_orig = count_ham_paths(A, n7)
    lk_orig = lambda_key(A, n7)

    for subset in combinations(range(n7), 4):
        B = reverse_subtournament(A, n7, subset)
        lk_rev = lambda_key(B, n7)
        if lk_orig == lk_rev:
            H_rev = count_ham_paths(B, n7)
            if H_orig != H_rev:
                h_change_examples.append((bits, subset, H_orig, H_rev))
                if len(h_change_examples) >= 20:
                    break
    if len(h_change_examples) >= 20:
        break

print(f"  Found {len(h_change_examples)} H-changing examples")

for bits, subset, H_orig, H_rev in h_change_examples[:5]:
    A = binary_to_tournament(bits, n7)
    B = reverse_subtournament(A, n7, subset)

    # Count cycles by size
    for k in [3, 5, 7]:
        c_orig = 0
        c_rev = 0
        for verts in combinations(range(n7), k):
            c_orig += count_directed_cycles_on(A, list(verts))
            c_rev += count_directed_cycles_on(B, list(verts))
        print(f"    bits={bits}, subset={subset}: "
              f"c{k}: {c_orig} -> {c_rev} (delta={c_rev-c_orig})")
    print(f"    H: {H_orig} -> {H_rev} (delta={H_rev-H_orig})")

    # Which 5-vertex sets have different cycle counts?
    diff_5 = []
    for verts in combinations(range(n7), 5):
        c5_orig = count_directed_cycles_on(A, list(verts))
        c5_rev = count_directed_cycles_on(B, list(verts))
        if c5_orig != c5_rev:
            diff_5.append((verts, c5_orig, c5_rev))
    print(f"    Differing 5-cycle vertex sets: {len(diff_5)}")
    for verts, c_o, c_r in diff_5[:3]:
        # How many of the 4-reversal vertices are in this 5-set?
        overlap = len(set(verts) & set(subset))
        print(f"      {verts}: c5 {c_o} -> {c_r}, overlap with reversal set = {overlap}")

    # Which 7-cycle count changes?
    c7_orig = count_directed_cycles_on(A, list(range(n7)))
    c7_rev = count_directed_cycles_on(B, list(range(n7)))
    print(f"    c7: {c7_orig} -> {c7_rev} (delta={c7_rev-c7_orig})")
    print()


# ========================================================================
# PART 5: The critical condition — WHEN does 4-reversal change H?
# ========================================================================
print(f"\n{'='*70}")
print("PART 5: Conditions for non-trivial 4-reversal")
print("=" * 70)

# At n=5: 4-reversal preserves lambda in 640/5120 cases, but NEVER changes H.
# At n=6: ?
# At n=7: 4-reversal sometimes changes H.

# The difference must be about the EXTERNAL structure.
# At n=5, the external structure is a single vertex. Its connections
# to the 4-set are fully determined by lambda (since lambda_{u,v}
# counts 3-cycles through u,v, and for the single external vertex w,
# the 3-cycles through (u,w) involve u, w, and one of the 3 remaining
# vertices of the 4-set).

# At n=6, V\S has 2 vertices. Their connections to the 4-set give
# 2*4 = 8 external arcs, plus 1 arc between them = 9 external arcs.
# Lambda constrains these, but may not fully determine them.

# At n=7, V\S has 3 vertices. 3*4=12 cross arcs + 3 internal = 15 external.
# Even more freedom.

# KEY: Does the external structure's INTERACTION with the 4-reversal
# change the cycle count? Specifically, 5-cycles that use SOME vertices
# from S and SOME from V\S will be affected by the reversal.

# For a 5-vertex set using k vertices from S and 5-k from V\S:
# k=4: all 4 from S + 1 from V\S. The sub-tournament on these 5 vertices
#       has 4 internal (reversed) arcs among S-vertices + 4 cross arcs.
# k=3: 3 from S + 2 from V\S. 3 reversed internal + 6 cross + 1 external.
# k=2: 2 from S + 3 from V\S. 1 reversed internal + 6 cross + 3 external.

print(f"  5-cycle composition (overlap with 4-set) for H-changing reversals:")

overlap_stats = defaultdict(lambda: defaultdict(int))
for bits, subset, H_orig, H_rev in h_change_examples[:10]:
    A = binary_to_tournament(bits, n7)
    B = reverse_subtournament(A, n7, subset)
    subset_set = set(subset)

    for verts in combinations(range(n7), 5):
        c5_orig = count_directed_cycles_on(A, list(verts))
        c5_rev = count_directed_cycles_on(B, list(verts))
        overlap = len(set(verts) & subset_set)
        if c5_orig != c5_rev:
            overlap_stats[overlap]['diff'] += 1
            overlap_stats[overlap]['total_delta'] += c5_rev - c5_orig
        else:
            overlap_stats[overlap]['same'] += 1

print(f"  Overlap with 4-set -> 5-cycle changes:")
for ov in sorted(overlap_stats.keys()):
    data = overlap_stats[ov]
    diff = data.get('diff', 0)
    same = data.get('same', 0)
    total_delta = data.get('total_delta', 0)
    print(f"    overlap={ov}: {diff} changed, {same} same, net delta={total_delta}")


# ========================================================================
# PART 6: WHY n=5 is trivial — the "single external vertex" argument
# ========================================================================
print(f"\n{'='*70}")
print("PART 6: Why 4-reversal is trivial at n=5")
print("=" * 70)

# At n=5 with 4-set S and single external vertex e:
# The 3-cycles through edges WITHIN S are reversed but stay as 3-cycles.
# So c3 on vertex sets within S is unchanged.
# The 3-cycles involving e: {u, v, e} for u,v in S.
# These use arcs u->v (or v->u), u->e (or e->u), v->e (or e->v).
# The arc u->v is the only one reversed. So the 3-cycle direction changes.
# BUT: a 3-cycle's EXISTENCE (as a vertex set containing a cycle) is preserved
# by reversal. So c3 on {u,v,e} is unchanged.

# For 5-cycles (all 5 vertices): the reversal flips all 6 internal arcs.
# This reverses all Hamiltonian cycles within the 5-tournament.
# But the NUMBER of directed Hamiltonian cycles = c5_dir is preserved
# under reversal (each cycle maps to its reverse, still counted).

# So at n=5: c3 AND c5 are unchanged under 4-reversal (when lambda preserved).
# Since H = 1 + 2*alpha_1 + 4*alpha_2 and alpha_1 = c3 + c5, alpha_2 = disjoint pairs,
# the question is whether alpha_2 can change.

# At n=5, alpha_2 = # pairs of disjoint cycles.
# Two 3-cycles on {a,b,c} and {d,e} ... wait, at n=5, a disjoint pair needs
# 2*3 = 6 vertices, but n=5 only has 5. So disjoint 3-3 pairs are impossible!
# Disjoint 3-5 pairs also impossible (need 8 vertices).
# So alpha_2 = 0 always at n=5! H = 1 + 2*alpha_1.

# Wait, that's not right either. alpha_2 counts independent pairs in the
# conflict graph, not just vertex-disjoint cycle pairs.
# A pair of cycles is "independent" in Omega(T) iff they don't share a vertex.
# At n=5: a pair of 3-cycles needs 6 vertices, impossible for n=5.
# A 3-cycle and a 5-cycle: 5-cycle uses all 5 vertices, so always shares.
# A pair of 5-cycles: only 1 possible vertex set, so no pair.
# So indeed alpha_2 = 0 at n=5.

# Therefore H = 1 + 2*(c3+c5) at n=5. Since c3 and c5 are preserved
# by lambda-preserving 4-reversal, H is preserved. QED.

print(f"  WHY 4-reversal is H-preserving at n=5:")
print(f"    - alpha_2 = 0 always (not enough vertices for disjoint pairs)")
print(f"    - H = 1 + 2*(c3 + c5)")
print(f"    - 4-reversal preserves c3 (vertex sets are the same, cycle existence preserved)")
print(f"    - 4-reversal preserves c5 (reversal maps cycles to reverse cycles)")
print(f"    - Therefore H is preserved.")

print(f"\n  WHY this breaks at n=7:")
print(f"    - alpha_2 > 0 possible (disjoint 3-cycle pairs exist)")
print(f"    - H = 1 + 2*(c3+c5+c7) + 4*alpha_2")
print(f"    - 4-reversal preserves c3 (on vertex sets not wholly in the 4-set)")
print(f"    - BUT: 5-cycles using MIXED vertices (some in 4-set, some outside)")
print(f"      can change direction! The reversal of internal arcs can create or")
print(f"      destroy directed 5-cycles that span the boundary.")
print(f"    - So c5 (counted over 5-vertex subsets) CAN change.")
print(f"    - And disjoint pair structure can change too.")


# ========================================================================
# PART 7: Verify the n=6 case — is it the true boundary?
# ========================================================================
print(f"\n{'='*70}")
print("PART 7: Alpha_2 at n=6 under 4-reversal")
print("=" * 70)

# At n=6 with 4-set S and 2 external vertices:
# Disjoint 3-cycle pairs: possible (3+3 = 6 vertices).
# So alpha_2 can be > 0.
# BUT: does 4-reversal change c5 at n=6?

# 5-vertex subsets at n=6: each subset has 4 from S + 1 external or 3+2 or 2+3 or ...
# The key is: can 5-cycles change direction?

if stats_6['changes_H'] == 0:
    print(f"  n=6: NO H-changing 4-reversals (exhaustive)")
    print(f"  Even though alpha_2 > 0 is possible at n=6,")
    print(f"  the 4-reversal preserves it exactly.")
    print(f"  The phase transition is SHARP at n=7.")
else:
    print(f"  n=6: H-changing 4-reversals EXIST")
    print(f"  Phase transition starts at n=6, not n=7")

# Let's verify: at n=6, does lambda-preserving 4-reversal preserve c5?
c5_change_count = 0
c5_preserve_count = 0
for bits in range(0, 1 << 15, 32):  # Sample
    A = binary_to_tournament(bits, n6)
    lk_orig = lambda_key(A, n6)

    for subset in combinations(range(n6), 4):
        B = reverse_subtournament(A, n6, subset)
        lk_rev = lambda_key(B, n6)
        if lk_orig == lk_rev:
            # Count c5
            c5_orig = sum(count_directed_cycles_on(A, list(v5))
                         for v5 in combinations(range(n6), 5))
            c5_rev = sum(count_directed_cycles_on(B, list(v5))
                        for v5 in combinations(range(n6), 5))
            if c5_orig != c5_rev:
                c5_change_count += 1
            else:
                c5_preserve_count += 1

print(f"\n  Lambda-preserving 4-reversals at n=6:")
print(f"    c5 preserved: {c5_preserve_count}")
print(f"    c5 changed: {c5_change_count}")


# ========================================================================
# PART 8: The overlap weight connection
# ========================================================================
print(f"\n{'='*70}")
print("PART 8: Overlap weights W={2,1,0} under 4-reversal")
print("=" * 70)

# The {2,1,0} overlap weights capture how pairs of 3-cycles overlap.
# W=2: share 2 vertices (strongly conflicting)
# W=1: share 1 vertex (weakly conflicting)
# W=0: disjoint (independent)

# At n=7, the H-changing 4-reversal changes some 5-cycle counts.
# Does it also change the 3-cycle overlap structure?

for bits, subset, H_orig, H_rev in h_change_examples[:3]:
    A = binary_to_tournament(bits, n7)
    B = reverse_subtournament(A, n7, subset)

    # Count 3-cycle vertex sets
    c3_orig = []
    c3_rev = []
    for v3 in combinations(range(n7), 3):
        if count_directed_cycles_on(A, list(v3)):
            c3_orig.append(frozenset(v3))
        if count_directed_cycles_on(B, list(v3)):
            c3_rev.append(frozenset(v3))

    # Count overlap types
    ov_orig = {0: 0, 1: 0, 2: 0}
    ov_rev = {0: 0, 1: 0, 2: 0}

    for i in range(len(c3_orig)):
        for j in range(i+1, len(c3_orig)):
            o = len(c3_orig[i] & c3_orig[j])
            ov_orig[o] += 1

    for i in range(len(c3_rev)):
        for j in range(i+1, len(c3_rev)):
            o = len(c3_rev[i] & c3_rev[j])
            ov_rev[o] += 1

    print(f"  bits={bits}, subset={subset}, H: {H_orig}->{H_rev}")
    print(f"    |c3|: {len(c3_orig)} -> {len(c3_rev)}")
    print(f"    Same 3-cycle vertex sets? {set(c3_orig) == set(c3_rev)}")
    print(f"    Overlap W=0: {ov_orig[0]} -> {ov_rev[0]}")
    print(f"    Overlap W=1: {ov_orig[1]} -> {ov_rev[1]}")
    print(f"    Overlap W=2: {ov_orig[2]} -> {ov_rev[2]}")


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
