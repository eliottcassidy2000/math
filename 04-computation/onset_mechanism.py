#!/usr/bin/env python3
"""
onset_mechanism.py -- kind-pasteur-2026-03-13-S61

WHY does the Vitali non-measurability onset at n=7?
At n=5,6: lambda-preserving 4-vertex reversals always preserve H.
At n=7: they CAN change c7 (and hence H).

Key question: what structural property makes c7 sensitive to
a 4-vertex reversal while c3, c5, alpha_2 are all preserved?

Hypotheses:
1. It's about the ratio: 6/21 arcs affected at n=7 vs 6/15 at n=6
2. It's about 7-cycles being able to "route around" the 4 reversed vertices
3. It's about the balance condition being satisfiable at n=7 but not n=5,6
4. It's related to the fact that 4+3=7 (partition of vertices into S and complement)

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


# ========================================================================
# ANALYSIS 1: The balance condition at different n
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: BALANCE CONDITION FREQUENCY BY n")
print("=" * 70)

# For a pair (u,v) and outside vertex w, the balance condition is:
# #{w: v->w and w->u} = #{w: u->w and w->v}
# This must hold for ALL 6 pairs in S and ALL outside vertices.

# But actually the condition is: for each pair (u,v) in S,
# #{w outside: A[v][w]*A[w][u] = 1} = #{w outside: A[u][w]*A[w][v] = 1}
# This is a sum over w outside S.

# At n=5: S has 4 vertices, outside has 1. One equation per pair, one outside vertex.
# For pair (u,v): A[v][w]*A[w][u] = A[u][w]*A[w][v] where w is the single outside vertex.
# This means: (v->w and w->u) iff (u->w and w->v).
# Expanding: A[v][w]*A[w][u] = A[u][w]*A[w][v]
# Let a = A[u][w], b = A[v][w] (so A[w][u]=1-a, A[w][v]=1-b).
# LHS: b*(1-a), RHS: a*(1-b). Equal iff b-ab = a-ab iff a = b.
# So: A[u][w] = A[v][w] for the single outside vertex w.
# This means u and w have same arc direction, and v and w have same arc direction.
# In other words: w treats u and v identically (beats both or loses to both).

# This must hold for ALL 6 pairs in S. So for each pair (u,v) in S:
# w must beat both u and v, or lose to both u and v.
# But S has 4 vertices, so there are C(4,2)=6 pairs.
# For w to beat all pairs: w beats all 4 vertices in S (out-degree to S = 4).
# For w to lose to all pairs: w loses to all 4 (out-degree to S = 0).
# BUT: for each pair, w can independently beat or lose. The condition is
# A[u][w] = A[v][w] for all u,v in S.
# This means A[u][w] is the SAME for all u in S. So w either beats ALL of S or loses to ALL.
# Out-degree of w to S is 0 or 4.

print("\n  n=5: Balance condition requires outside vertex w to beat ALL or NONE of S")
n5 = 5
total5 = 1 << (n5*(n5-1)//2)
balanced_count = 0
for bits in range(total5):
    A = binary_to_tournament(bits, n5)
    for S in combinations(range(n5), 4):
        outside = [v for v in range(n5) if v not in S]
        w = outside[0]
        # Check: does w beat all or none of S?
        out_deg_to_S = sum(A[w][s] for s in S)
        if out_deg_to_S == 0 or out_deg_to_S == 4:
            balanced_count += 1
print(f"    Total (tournament, S) pairs: {total5 * 5}")
print(f"    Balanced: {balanced_count} ({100*balanced_count/(total5*5):.1f}%)")

# At n=6: S has 4 vertices, outside has 2 (say w1, w2).
# For each pair (u,v) in S:
#   #{w in {w1,w2}: v->w and w->u} = #{w in {w1,w2}: u->w and w->v}
# This is sum over w in {w1,w2}. Each w contributes independently.
# The condition can be written: for each pair (u,v),
#   sum_w [A[v][w]*A[w][u] - A[u][w]*A[w][v]] = 0
# For each outside w and pair (u,v): contribution = A[v][w](1-A[u][w]) - A[u][w](1-A[v][w])
#   = A[v][w] - A[v][w]*A[u][w] - A[u][w] + A[u][w]*A[v][w]
#   = A[v][w] - A[u][w]
# So the condition is: for each pair (u,v) in S:
#   sum_{w outside} (A[v][w] - A[u][w]) = 0
# i.e., the out-degree of v to outside = out-degree of u to outside.

# Wait, that's for EACH pair. So ALL vertices in S must have the SAME
# out-degree to the outside complement!

print(f"\n  n=6: Balance condition requires ALL 4 vertices in S to have SAME out-degree to outside")
n6 = 6
total6 = 1 << (n6*(n6-1)//2)
balanced_count_6 = 0
for bits in range(total6):
    A = binary_to_tournament(bits, n6)
    for S in combinations(range(n6), 4):
        outside = [v for v in range(n6) if v not in S]
        out_degs = [sum(A[s][w] for w in outside) for s in S]
        if len(set(out_degs)) == 1:
            balanced_count_6 += 1
print(f"    Total (tournament, S) pairs: {total6 * 15}")
print(f"    Balanced: {balanced_count_6} ({100*balanced_count_6/(total6*15):.1f}%)")

# At n=7: S has 4 vertices, outside has 3.
# Same logic: for each pair (u,v) in S:
#   sum_{w outside} (A[v][w] - A[u][w]) = 0
# So all vertices in S have same out-degree to outside.
# But wait: |outside| = 3, so out-degree to outside is in {0,1,2,3}.
# All 4 vertices having same out-degree to 3 outside vertices.

print(f"\n  n=7: Balance requires ALL 4 vertices in S to have SAME out-degree to 3 outside vertices")

# Let's check: among the lambda-preserving reversals we found, is this condition met?
n = 7
A2 = binary_to_tournament(4658, n)
S = [1, 2, 3, 5]
outside = [v for v in range(n) if v not in S]
print(f"\n  Case 1 (bits=4658): S={S}, outside={outside}")
for s in S:
    od = sum(A2[s][w] for w in outside)
    print(f"    vertex {s}: out-degree to outside = {od}")


# Wait -- the simplified condition might not be right. Let me re-derive.
# For pair (u,v) in S and outside set W:
# condition: sum_{w in W} A[v][w]*(1-A[u][w]) = sum_{w in W} A[u][w]*(1-A[v][w])
# = sum_w A[v][w] - A[v][w]*A[u][w] = sum_w A[u][w] - A[u][w]*A[v][w]
# = sum_w A[v][w] - sum_w A[v][w]*A[u][w] = sum_w A[u][w] - sum_w A[u][w]*A[v][w]
# Since sum A[v][w]*A[u][w] = sum A[u][w]*A[v][w], both cross terms are equal.
# So: sum_w A[v][w] = sum_w A[u][w]
# YES! The condition simplifies to: out-degree of u to outside = out-degree of v to outside,
# for ALL pairs (u,v) in S.
# Equivalently: ALL vertices in S have the same out-degree to outside(S).

# So the necessary condition for lambda-preservation is very simple.
# But it's not sufficient (it's necessary from the pair-level balance, but
# we also need the lambda graph to actually be preserved, which involves
# the triangle count matching).

# Actually wait: the analysis showed that the balance condition for each pair (u,v) is
# #{w: v->w, w->u} = #{w: u->w, w->v}
# which we showed = sum_w A[v][w] = sum_w A[u][w] (same out-degree to outside).
# This is NECESSARY AND SUFFICIENT for lambda preservation under 4-vertex reversal!
# Because lambda[a][b] for a,b both in S changes only via third vertices in S (which
# get ALL arcs flipped, preserving the intra-S 3-cycle count by symmetry) or outside
# (which requires the balance condition). And for a in S, b outside, or both outside,
# the arcs between S and outside are unchanged.

# Wait, I need to be more careful. Lambda[a][b] for a in S, b outside:
# 3-cycles through {a,b,...}: third vertex c is either in S or outside.
# If c in S: arcs a<->c flip (both in S). But b<->c doesn't flip (c in S, b outside).
# Hmm, actually arc b->c or c->b doesn't change since b is outside and c is in S.
# But arc a->c flips. So the 3-cycle {a,b,c} might change!

# Oh, I think the analysis is more complex than I thought. Let me check directly.

print(f"\n\n  DETAILED: Lambda preservation for a in S, b outside S:")
# Check lambda[s][w] for s in S, w outside
for s in S:
    for w in outside:
        # Third vertex x: 3-cycles through {s,w,x}
        # x in S: arc s<->x flips. Does {s,w,x} status change?
        # x outside: nothing changes.
        changes = []
        for x in range(n):
            if x == s or x == w:
                continue
            triple = [s, w, x]
            cyc_old = count_directed_ham_cycles_on_subset(A2, triple) > 0

            # After reversal
            A_rev = [row[:] for row in A2]
            for i in S:
                for j in S:
                    if i != j:
                        A_rev[i][j] = 1 - A2[i][j]
            cyc_new = count_directed_ham_cycles_on_subset(A_rev, triple) > 0

            if cyc_old != cyc_new:
                changes.append((x, 'in_S' if x in S else 'outside'))
        if changes:
            print(f"    lambda[{s}][{w}]: changes from x={changes}")

# Lambda preservation verification
A_rev = [row[:] for row in A2]
for i in S:
    for j in S:
        if i != j:
            A_rev[i][j] = 1 - A2[i][j]
lam_old, _ = get_labeled_lambda(A2, n)
lam_new, _ = get_labeled_lambda(A_rev, n)

print(f"\n  Lambda comparison (old vs new):")
for i in range(n):
    for j in range(i+1, n):
        if lam_old[i][j] != lam_new[i][j]:
            loc = []
            if i in S: loc.append(f"{i}:S")
            else: loc.append(f"{i}:out")
            if j in S: loc.append(f"{j}:S")
            else: loc.append(f"{j}:out")
            print(f"    ({i},{j}) [{','.join(loc)}]: {lam_old[i][j]} -> {lam_new[i][j]}")

if all(lam_old[i][j] == lam_new[i][j] for i in range(n) for j in range(i+1, n)):
    print(f"    ALL lambda values preserved!")


# ========================================================================
# ANALYSIS 2: Why c7 changes
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: WHY c7 CHANGES UNDER THE REVERSAL")
print("=" * 70)

# A directed Hamiltonian cycle visits all 7 vertices.
# It crosses between S and outside multiple times.
# Each crossing uses an arc between S and outside (unchanged by reversal).
# The arcs WITHIN S are reversed.
# A HC that uses arcs within S will be affected.

# For a HC (v0,v1,...,v6) visiting all vertices:
# The arcs within S are those where consecutive vertices vi,vi+1 are both in S.
# After reversal, those arcs flip direction.

# A HC is preserved iff it uses zero arcs within S (impossible if |S|>=2 and HC visits all)
# OR it uses arcs within S that happen to form a valid HC segment even after reversal.

# Actually: |S|=4. A HC must visit all 4 S-vertices. Between consecutive S-vertices
# in the HC, there might be outside vertices. The S-vertices form a subsequence of the HC.
# The arcs used WITHIN S are the ones connecting consecutive S-vertices in the HC
# (when no outside vertex is between them).

# Case: HC visits S-vertices in order s1,s2,s3,s4 (with possible outside vertices between).
# If s_i and s_{i+1} are CONSECUTIVE in HC (no outside between): arc s_i->s_{i+1} flips.
# If there's an outside vertex between: no arc within S is used for this transition.

# So the effect depends on how many "S-consecutive pairs" the HC has.

for bits, label in [(4658, "H=111"), (4728, "H=109")]:
    A = binary_to_tournament(bits, n)
    if bits == 4658:
        S_set = set([1, 2, 3, 5])
    else:
        # Need to transform: the perm was (0,1,5,2,3,4,6) applied to A1
        # But let's just check which 4-vertex reversal preserves lambda for A1
        lam_orig, _ = get_labeled_lambda(A, n)
        found_S = None
        for S_cand in combinations(range(n), 4):
            A_test = [row[:] for row in A]
            for i in S_cand:
                for j in S_cand:
                    if i != j:
                        A_test[i][j] = 1 - A[i][j]
            lam_test, _ = get_labeled_lambda(A_test, n)
            if all(lam_test[i][j] == lam_orig[i][j] for i in range(n) for j in range(i+1, n)):
                found_S = list(S_cand)
                break
        if found_S:
            S_set = set(found_S)
            print(f"\n  {label}: lambda-preserving S = {found_S}")
        else:
            S_set = set()
            print(f"\n  {label}: NO lambda-preserving 4-reversal found")
            continue

    # Find all directed HCs
    hcs = []
    for perm in permutations(range(n)):
        if perm[0] != 0:
            continue
        if all(A[perm[i]][perm[(i+1) % n]] for i in range(n)):
            hcs.append(perm)

    print(f"    {len(hcs)} directed HCs")

    for hc in hcs:
        # Count consecutive pairs within S
        s_consec = 0
        s_arcs = []
        for i in range(n):
            u, v = hc[i], hc[(i+1) % n]
            if u in S_set and v in S_set:
                s_consec += 1
                s_arcs.append((u, v))
        print(f"    HC {hc}: {s_consec} S-internal arcs: {s_arcs}")

    # After reversal: which HCs survive?
    A_rev = [row[:] for row in A]
    for i in S_set:
        for j in S_set:
            if i != j:
                A_rev[i][j] = 1 - A[i][j]

    hcs_rev = []
    for perm in permutations(range(n)):
        if perm[0] != 0:
            continue
        if all(A_rev[perm[i]][perm[(i+1) % n]] for i in range(n)):
            hcs_rev.append(perm)

    print(f"    After reversal: {len(hcs_rev)} directed HCs")

    # Which HCs are shared?
    hc_set_old = set(hcs)
    hc_set_new = set(hcs_rev)
    shared = hc_set_old & hc_set_new
    destroyed = hc_set_old - hc_set_new
    created = hc_set_new - hc_set_old
    print(f"    Shared: {len(shared)}, Destroyed: {len(destroyed)}, Created: {len(created)}")
    print(f"    Net change: {len(created) - len(destroyed)}")

    for hc in sorted(destroyed):
        print(f"    DESTROYED: {hc}")
    for hc in sorted(created):
        print(f"    CREATED:   {hc}")


# ========================================================================
# ANALYSIS 3: The partition 4+3=7 and why it matters
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: THE PARTITION 4+3=7")
print("=" * 70)

# At n=5: 4+1=5. Outside has 1 vertex. Very constrained.
# At n=6: 4+2=6. Outside has 2 vertices.
# At n=7: 4+3=7. Outside has 3 vertices.

# The balance condition requires all S-vertices to have the same out-degree
# to outside. At n=5 with 1 outside vertex: all must beat w or all must lose.
# Only 2 options. Very rigid.

# At n=6 with 2 outside: out-degree to outside in {0,1,2}. All 4 must match.
# 3 options. Still rigid.

# At n=7 with 3 outside: out-degree to outside in {0,1,2,3}. All 4 must match.
# 4 options. More flexible.

# But the key is: does this flexibility ALLOW lambda preservation to happen
# while c7 changes?

# The c7 change happens because a 7-cycle visits ALL vertices. The 4-vertex
# reversal affects the internal routing. With 3 outside vertices, the HC
# can enter and exit S multiple times with different routing.

# At n=5: a 5-cycle visits all 5. With |S|=4, the single outside vertex
# creates a rigid structure. The 5-cycle must use exactly 3 internal S-arcs
# and 2 S-outside arcs.

# At n=7: a 7-cycle with |S|=4 uses between 0 and 3 internal S-arcs.
# The variety of routing options creates room for c7 to change.

# Let's count: how many internal S-arcs does each HC use?
print("\n  Distribution of S-internal arc counts in HCs:")
for bits, label in [(4658, "H=111")]:
    A = binary_to_tournament(bits, n)
    S_set = set([1, 2, 3, 5])

    arc_count_dist = defaultdict(int)
    for perm in permutations(range(n)):
        if perm[0] != 0:
            continue
        if all(A[perm[i]][perm[(i+1) % n]] for i in range(n)):
            count = sum(1 for i in range(n) if perm[i] in S_set and perm[(i+1)%n] in S_set)
            arc_count_dist[count] += 1

    total_hc = sum(arc_count_dist.values())
    for k in sorted(arc_count_dist.keys()):
        print(f"    {k} S-internal arcs: {arc_count_dist[k]} HCs ({100*arc_count_dist[k]/total_hc:.0f}%)")


# ========================================================================
# ANALYSIS 4: Can 5-vertex reversals work at n=8?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: PREDICTION FOR n=8")
print("=" * 70)

print("""
At n=7: 4-vertex reversal affects C(4,2)=6 out of C(7,2)=21 arcs (28.6%).
Balance requires all 4 S-vertices to have same out-degree to 3 outside vertices.

At n=8: 4-vertex reversal affects 6/28 = 21.4% of arcs. Outside has 4.
Same balance condition, but with more outside flexibility.

PREDICTION: Lambda-preserving 4-vertex reversals that change c_n exist for ALL n >= 7.
The onset is EXACTLY at n=7 because:
1. 4+3=7 is the first partition where outside has >= 3 vertices
2. 3 outside vertices = enough routing options for HCs
3. The balance condition (equal out-degree to outside) has multiple solutions

For n=5: outside=1, only two routing options. Too rigid.
For n=6: outside=2, three routing options. Still too rigid (no c6 to change anyway,
   since c6 is always 0 for tournaments and 6-cycles don't exist as they'd be even).
For n=7: outside=3, four routing options. Just enough freedom for c7 to vary.

WAIT: Actually the issue at n=6 is simpler: there ARE no Hamiltonian 6-cycles
in a tournament on 6 vertices (6 is even). So c6 doesn't exist.
The largest odd cycle at n=6 is c5. But c5 involves 5 vertices, leaving only 1 out.
A 4-vertex reversal on 5-cycles: S has 4 of the 5 cycle vertices,
plus 1 outside vertex = not enough freedom.

This gives the KEY INSIGHT: n=7 onset because it's the FIRST n where:
- n is odd (so c_n exists as directed Hamiltonian cycles)
- n - 4 = 3 >= 3 outside vertices (enough routing flexibility)
- The balance condition can be satisfied non-trivially
""")


# ========================================================================
# ANALYSIS 5: The 2:1 fiber ratio (H=111 twice as common)
# ========================================================================
print(f"{'='*70}")
print("ANALYSIS 5: THE 2:1 FIBER RATIO")
print("=" * 70)

# We found: 5040 tournaments with H=109 (c7=8) and 10080 with H=111 (c7=9).
# The 2:1 ratio. Why?

# Each tournament with c7=8 has exactly one lambda-preserving reversal
# (we found only 1 out of 35). The reversal takes it to c7=9.
# But each c7=9 tournament could have TWO different reversals back to c7=8.

# Or: the orbits under the reversal group have different sizes.
# If each c7=8 tournament pairs uniquely with one c7=9 tournament,
# but each c7=9 tournament can be reached from two different c7=8 tournaments...

# Let's check: for a few c7=9 tournaments, how many lambda-preserving
# 4-reversals do they have?

print(f"  Checking lambda-preserving reversals for H=111 tournaments:")
count_checked = 0
by_rev_count = defaultdict(int)
for bits, H, c7 in [(b, h, c) for b, h, c in
    [(4658, 111, 9)] +
    [(b, h, c) for b, h, c in
     sorted([(bits, H, c7) for bits, H, c7 in
             [(bits, count_ham_paths(binary_to_tournament(bits, n), n),
               count_directed_ham_cycles_on_subset(binary_to_tournament(bits, n), list(range(n))))
              for bits in range(0, 1 << 21, 2048)
              if tuple(sorted(sum(binary_to_tournament(bits, n)[v]) for v in range(n))) == (2, 2, 3, 3, 3, 3, 5)]
             if H == 111])[:5]]]:
    A = binary_to_tournament(bits, n)
    if count_ham_paths(A, n) != 111:
        continue
    lam, _ = get_labeled_lambda(A, n)
    rev_count = 0
    for S in combinations(range(n), 4):
        A_rev = [row[:] for row in A]
        for i in S:
            for j in S:
                if i != j:
                    A_rev[i][j] = 1 - A[i][j]
        lam_rev, _ = get_labeled_lambda(A_rev, n)
        if all(lam_rev[i][j] == lam[i][j] for i in range(n) for j in range(i+1, n)):
            H_rev = count_ham_paths(A_rev, n)
            print(f"    bits={bits}: rev at S={S}, H_rev={H_rev}")
            rev_count += 1
    by_rev_count[rev_count] += 1
    count_checked += 1
    if count_checked >= 5:
        break

# Also check some H=109 tournaments
print(f"\n  Checking lambda-preserving reversals for H=109 tournaments:")
for bits in [4728, 9388]:
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    if H != 109:
        continue
    lam, _ = get_labeled_lambda(A, n)
    rev_count = 0
    for S in combinations(range(n), 4):
        A_rev = [row[:] for row in A]
        for i in S:
            for j in S:
                if i != j:
                    A_rev[i][j] = 1 - A[i][j]
        lam_rev, _ = get_labeled_lambda(A_rev, n)
        if all(lam_rev[i][j] == lam[i][j] for i in range(n) for j in range(i+1, n)):
            H_rev = count_ham_paths(A_rev, n)
            print(f"    bits={bits}: rev at S={S}, H_rev={H_rev}")
            rev_count += 1
    print(f"    Total reversals: {rev_count}")


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
