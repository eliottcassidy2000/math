import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
c6_correction_formula.py
kind-pasteur-2026-03-07-S39b

Find the correction formula for tr(A^6) - 6*c_6.

At k=6, non-simple closed walks decompose into (3,3)-compound walks:
two directed 3-cycles sharing a junction vertex v.

Walk: v -> a -> b -> v -> c -> d -> v (length 6)

For each vertex v and each ORDERED pair of directed 3-cycles (C1, C2)
through v, we get ONE compound walk starting at v.

Question: what is the exact multiplicity?
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations
from collections import defaultdict


def count_directed_k_cycles_dp(T, k):
    n = len(T)
    if k > n or k < 3:
        return 0
    count = 0
    for verts in combinations(range(n), k):
        v = list(verts)
        dp = [[0] * k for _ in range(1 << k)]
        dp[1][0] = 1
        for mask in range(1, 1 << k):
            for last in range(k):
                if dp[mask][last] == 0 or not (mask & (1 << last)):
                    continue
                for nxt in range(1, k):
                    if mask & (1 << nxt):
                        continue
                    if T[v[last]][v[nxt]]:
                        dp[mask | (1 << nxt)][nxt] += dp[mask][last]
        full = (1 << k) - 1
        for last in range(1, k):
            if T[v[last]][v[0]]:
                count += dp[full][last]
    return count


def matrix_power_trace(T, k):
    n = len(T)
    Ak = [[int(i == j) for j in range(n)] for i in range(n)]
    for _ in range(k):
        new = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for l in range(n):
                    new[i][j] += Ak[i][l] * T[l][j]
        Ak = new
    return sum(Ak[i][i] for i in range(n))


print("=" * 70)
print("c_6 CORRECTION FORMULA")
print("=" * 70)

# For each vertex v, count the number of compound (3,3) walks starting at v.
# A compound walk at v: v -> a -> b -> v -> c -> d -> v
# This is determined by choosing:
#   - A directed 3-cycle through v: v->a->b->v (the 3-cycle is a->b->v->a, traversed as v->a->b->v)
#   - Another directed 3-cycle through v: v->c->d->v
# These can be the SAME cycle (walk goes around it twice) or different.
# The walk is valid iff v->a, a->b, b->v, v->c, c->d, d->v are all edges.
# Since v->a->b->v and v->c->d->v are directed 3-cycles through v,
# all required edges exist. So the count is t3(v)^2.

# But how many of these walks are captured as "non-simple walk starting at v"
# in tr(A^6)?

# A compound walk v->a->b->v->c->d->v has 6 rotations (length 6):
# Rot 0: v->a->b->v->c->d->v  (starts at v, first occurrence)
# Rot 1: a->b->v->c->d->v->a  (starts at a)
# Rot 2: b->v->c->d->v->a->b  (starts at b)
# Rot 3: v->c->d->v->a->b->v  (starts at v, second occurrence)
# Rot 4: c->d->v->a->b->v->c  (starts at c)
# Rot 5: d->v->a->b->v->c->d  (starts at d)

# Rotations 0 and 3 start at v (the junction vertex).
# Rotation 0 is the (C1-first, C2-second) arrangement.
# Rotation 3 is the (C2-first, C1-second) arrangement.

# So each UNORDERED pair {C1, C2} generates 2 compound walks starting at v:
#   (C1-first, C2-second) and (C2-first, C1-second).
# For the ORDERED pair (C1, C2): the walk starting at v is unique.
# For the ORDERED pair (C2, C1): a different walk starting at v.

# But wait: t3(v)^2 counts ALL ordered pairs (C1, C2), which gives
# both orderings naturally. So t3(v)^2 counts exactly the compound
# walks starting at v from rotations 0 and 3.

# The remaining 4 rotations (1,2,4,5) start at a,b,c,d.
# Are these also counted by some other vertex's t3^2?

# Rotation 1 (starts at a): a->b->v->c->d->v->a
# This walk repeats v (not a). It's a compound walk centered at v.
# It is NOT counted by t3(a)^2 because a appears only once.

# So the 4 non-junction rotations are NOT counted by any t3^2 sum.
# Total non-simple walk contributions to tr(A^6):
# For each ordered pair (C1, C2) through v: generates 1 walk, which
# has 6 rotations contributing to tr(A^6). Of these, exactly 2 are
# counted by the t3^2 approach (the two starting at v).

# Wait, no. t3(v)^2 counts ordered pairs. Each ordered pair gives
# ONE walk starting at v. So t3(v)^2 = number of walks starting at v.
# But the ordered pair (C1, C2) and (C2, C1) give rotation 0 and
# rotation 3 of the SAME walk. So t3(v)^2 counts BOTH rotations.

# So: sum_v t3(v)^2 counts EXACTLY the walks starting at junction
# vertices. Each walk has 6 rotations, of which 2 start at the
# junction. So sum_v t3(v)^2 = 2/6 * total_nonsimple_walks.
# Therefore: excess = 3 * sum_v t3(v)^2.

# BUT: the test showed this is WRONG. Why?

# Possibility: when C1 and C2 share vertices other than v, the
# walk has MORE than one junction vertex. If C1 and C2 share vertex w
# (besides v), then w appears at two positions in the walk. The walk
# then has TWO junction vertices: v and w. Rotation starting at w
# IS a compound walk centered at w, counted by t3(w)^2.

# In this case, the walk has 6 rotations, of which:
# - 2 start at junction v (counted by t3(v)^2)
# - 2 start at junction w (counted by t3(w)^2)
# - 2 start at non-junction vertices (not counted)
# Total counted = 4 (not 2). So factor would be 6/4 = 3/2, not 3.

# For walks where C1 = C2 (same cycle traversed twice), ALL 3 vertices
# of the cycle are junction vertices. So 6 out of 6 rotations are
# counted. Factor = 6/6 = 1.

# So the correction factor depends on the overlap structure!

# Let me compute this correctly by partitioning walks by overlap type.

for n in [5, 6]:
    m = n * (n - 1) // 2
    total_tours = 1 << m

    print(f"\nn={n} (exhaustive, {total_tours} tournaments):")

    # For each tournament, compute various candidate formulas
    formulas = defaultdict(int)

    for bits in range(total_tours):
        T = tournament_from_bits(n, bits)
        nn = len(T)

        tr6 = matrix_power_trace(T, 6)
        c6 = count_directed_k_cycles_dp(T, 6) if nn >= 6 else 0
        excess = tr6 - 6 * c6

        # Compute directed 3-cycles through each vertex
        # d3_through_v[v] = list of (a,b) such that v->a->b->v
        d3 = defaultdict(list)
        for v in range(nn):
            for a in range(nn):
                if a == v or not T[v][a]:
                    continue
                for b in range(nn):
                    if b == v or b == a or not T[a][b] or not T[b][v]:
                        continue
                    d3[v].append((a, b))

        t3v = {v: len(d3[v]) for v in range(nn)}

        # Compute exact non-simple walk count by categorizing ordered pairs
        # by overlap type
        walks_disjoint = 0  # C1, C2 share only v (cycles on disjoint vertex sets minus v)
        walks_overlap1 = 0  # C1, C2 share v and one other vertex
        walks_same = 0      # C1 = C2 (same cycle traversed twice)

        for v in range(nn):
            for c1 in d3[v]:
                for c2 in d3[v]:
                    s1 = {c1[0], c1[1]}
                    s2 = {c2[0], c2[1]}
                    overlap = s1 & s2  # overlap excluding v
                    if c1 == c2:
                        walks_same += 1
                    elif len(overlap) == 0:
                        walks_disjoint += 1
                    elif len(overlap) == 1:
                        walks_overlap1 += 1
                    # overlap == 2 means s1 == s2 but c1 != c2, which means
                    # same vertex set different ordering.
                    # In a tournament, on {v,a,b} there's exactly 1 directed 3-cycle.
                    # So s1 == s2 and c1 != c2 is impossible. Let me verify...
                    elif len(overlap) == 2:
                        # Same vertex set! But in tournament, unique directed 3-cycle.
                        # This shouldn't happen.
                        walks_overlap1 += 1  # shouldn't reach here

        sum_t3_sq = sum(t3v[v]**2 for v in range(nn))
        assert sum_t3_sq == walks_disjoint + walks_overlap1 + walks_same

        # Formula derivation:
        # - Disjoint pair: walk has 1 junction vertex (v). 2 out of 6 rotations
        #   are counted by t3^2. Contribution to excess = 6/2 * walks_disjoint = 3*walks_disjoint.
        # - Overlap-1 pair: walk has 2 junction vertices (v and the shared vertex).
        #   4 out of 6 rotations are counted. Contribution = 6/4 * walks_overlap1 = 1.5*walks_overlap1.
        #   But wait, this means each such walk is counted TWICE in sum_v t3^2:
        #   once at v, once at the other junction. But I'm counting per-v.
        #   Hmm, I need to think about this more carefully.

        # Let me just directly count walks through junction vertex analysis.
        # For each ordered pair (v, C1, C2), one walk is:
        # w = v -> C1[0] -> C1[1] -> v -> C2[0] -> C2[1] -> v
        # The junction vertices of this walk are the vertices that appear >= 2 times.
        # v always appears 3 times (positions 0, 3, 6).
        # C1[0] appears at position 1. C2[0] appears at position 4.
        # C1[1] appears at position 2. C2[1] appears at position 5.
        # Junction beyond v: any vertex in {C1[0],C1[1]} that also appears in {C2[0],C2[1]}.
        # So extra junctions = |{C1[0],C1[1]} & {C2[0],C2[1]}|.

        # For EACH junction vertex, 2 out of 6 rotations start there. But v appears
        # 3 times in the walk, not 2! Wait, v is at positions 0, 3, 6 (=0). So in the
        # circular walk of length 6 (positions 0-5), v is at positions 0 and 3.
        # An overlap vertex w (say C1[0]=C2[0]) is at positions 1 and 4.
        # So 2 rotations start at v, 2 at w, 2 at other vertices.

        # Total rotations: 6 per walk.
        # Rotations starting at junction vertices:
        # - disjoint: 2 (only v)
        # - overlap-1: 4 (v + 1 other)
        # - same cycle: 6 (v + both others)

        # Total contribution to sum_v t3(v)^2 from this walk:
        # Actually, sum_v t3(v)^2 counts ordered pairs per vertex.
        # Walk W = (v, C1, C2) is counted once (in t3(v)^2).
        # But if w is also a junction, the walk rotated to start at w
        # can be decomposed as (w, C1', C2') where C1' and C2' are
        # 3-cycles through w. This pair (w, C1', C2') is counted in t3(w)^2.

        # So for a walk with K junction vertices:
        # sum_v t3(v)^2 counts it 2K times (2 rotations per junction vertex,
        # each being an ordered pair at that junction).
        # Total excess from this walk = 6 (number of rotations).
        # So: 6 = (excess per walk) = (sum_t3_sq contribution) * (6 / (2K))
        # => excess per walk counted by sum_t3_sq = 2K times
        # => excess = sum (6 for each walk) = sum_walks 6
        # => sum_t3_sq = sum_walks 2K
        # => excess = sum_walks 6 = 3 * sum_walks 2K / K... no that's not right.

        # Let me think again. Let W be a distinct non-simple closed walk (as a
        # circular sequence). It has 6 rotations, each contributing 1 to tr(A^6).
        # So the total contribution of W to excess is 6.
        # Let K = number of junction vertices of W.
        # W has 2K rotations starting at junction vertices.
        # Each such rotation, when decomposed at its starting junction vertex,
        # gives an ordered pair of 3-cycles through that vertex.
        # So W contributes 2K to sum_v t3(v)^2.

        # Therefore: sum_v t3(v)^2 = sum_W 2K_W
        # excess = sum_W 6

        # If all walks had K=1: excess = 3 * sum_t3_sq. (Since 6/2 = 3.)
        # If all walks had K=2: excess = 6/(2*2) * sum_t3_sq = 3/2 * sum_t3_sq.
        # If all walks had K=3: excess = 6/(2*3) * sum_t3_sq = sum_t3_sq.

        # For the CORRECT formula, we need to weight by junction count:
        # excess = sum_v [disjoint_pairs(v) * 3 + overlap1_pairs(v) * ??? + same_pairs(v) * ???]
        # This doesn't factor nicely.

        # BUT: we can compute sum_v (disjoint + overlap + same) separately.
        # excess = 3*walks_disjoint/N_d + 3/2*walks_overlap1/N_o + 1*walks_same/N_s?
        # No, that's not right either. Let me just compute it directly.

        # Direct: count distinct walks and their rotations.
        # Each ordered pair (v, C1, C2) gives ONE walk starting at v.
        # Two ordered pairs give the SAME walk orbit iff one is a rotation of the other.

        # For (v, C1, C2):
        # The walk is: v, C1[0], C1[1], v, C2[0], C2[1]
        # Rotations:
        # Start at position 0: v, C1[0], C1[1], v, C2[0], C2[1] -> pair (v, C1, C2)
        # Start at position 1: C1[0], C1[1], v, C2[0], C2[1], v -> pair (C1[0], (C1[1],v), (C2[0],C2[1]))
        #   Wait, is this a valid decomposition? The walk starting at C1[0]:
        #   C1[0] -> C1[1] -> v -> C2[0] -> C2[1] -> v -> C1[0]
        #   Junction vertex is v (appears at positions 2 and 5).
        #   This decomposes as: first sub-walk C1[0]->C1[1]->v (length 2, not a closed walk!)
        #   NO. The junction is where the walk REPEATS a vertex. In this rotated walk:
        #   C1[0], C1[1], v, C2[0], C2[1], v
        #   v repeats at positions 2 and 5. So the decomposition at v gives:
        #   Sub-walk 1: C1[0] -> C1[1] -> v (open walk from C1[0] to v, length 2)
        #   That's NOT a closed walk. So this is NOT a compound walk in the sense I was thinking.

        # I was confused. A compound walk means: decompose at a REPEATED vertex.
        # The repeated vertex appears at positions i and j (0 <= i < j < 6).
        # Sub-walk 1: w_i, ..., w_j is a closed walk of length j-i.
        # Sub-walk 2: w_j, ..., w_6, w_0, ..., w_i is a closed walk of length 6-(j-i).
        # Both lengths >= 3.

        # For walk v, a, b, v, c, d (circular):
        # v at positions 0 and 3. Sub-walks: (0-3) length 3, (3-0) length 3. YES.
        # If a = c (overlap): a at positions 1 and 4. Sub-walks: (1-4) length 3, (4-1) length 3. YES.
        # If b = c: b at positions 2, c at position 4. That's b at 2 and c=b at 4.
        #   Sub-walks: (2-4) length 2, (4-2) length 4. Length 2 < 3. NO!
        #   So b=c does NOT create a new junction! The sub-walk of length 2 (b->v->b)
        #   is a 2-cycle which doesn't exist in tournaments.

        # Wait, but b=c means the walk is v->a->b->v->b->d->v. This requires edge v->b
        # AND b->v. But b->v is part of the first 3-cycle (a->b->v), so T[b][v]=1.
        # And v->b would be in the second 3-cycle (v->c->d->v with c=b), so T[v][b]=1.
        # But T[b][v]=1 and T[v][b]=1 contradicts tournament (no bidirectional edges).
        # So b=c is IMPOSSIBLE in a tournament!

        # Let me check: if the first 3-cycle is v->a->b->v (T[v][a]=T[a][b]=T[b][v]=1),
        # and the second is v->c->d->v (T[v][c]=T[c][d]=T[d][v]=1).
        # If b=c: T[b][v]=1 (from first cycle) and T[v][b]=T[v][c]=1 (from second cycle).
        # T[b][v]=1 and T[v][b]=1 => both directions, IMPOSSIBLE in tournament.

        # Similarly b=d is impossible: T[b][v]=1 (first cycle), T[d][v]=1 so T[v][d]=0
        #   but we need T[v][d]... wait d is in second cycle's last position, T[d][v]=1.
        #   If b=d: T[b][v]=1 and T[d][v]=T[b][v]=1, that's fine. No contradiction.
        #   Actually let me check edges more carefully:
        #   First cycle: v->a->b->v needs T[v][a]=1, T[a][b]=1, T[b][v]=1
        #   Second cycle: v->c->d->v needs T[v][c]=1, T[c][d]=1, T[d][v]=1
        #   If b=d: T[b][v]=1 (first) and T[d][v]=T[b][v]=1 (second). OK, no contradiction.
        #   But b appears at position 2, d at position 5. Distance 3. Junction at b:
        #   Sub-walks: (2-5) length 3: b,v,c,d=b -> walk b->v->c->b, which is closed if c->b.
        #   (5-2) length 3: d=b, v, a, b -> walk b->v->a->b, closed if a->b. Yes (T[a][b]=1).
        #   So b=d IS possible and creates a valid second junction.

        # What about a=c? T[v][a]=1 (first) and T[v][c]=T[v][a]=1 (second). OK.
        #   a at position 1, c=a at position 4. Distance 3.
        #   Sub-walks: (1-4) length 3: a,b,v,c=a -> walk a->b->v->a, closed (T[v][a]=1). YES.
        #   (4-1) length 3: c=a,d,v,a -> walk a->d->v->a, closed (T[v][a]=1). YES.
        #   So a=c IS possible.

        # What about a=d? a at 1, d at 5. Distance 4.
        #   Sub-walks: (1-5) length 4: a,b,v,c,d=a -> walk a->b->v->c->a. Length 4.
        #   (5-1) length 2: d=a, v, a -> walk a->v->a. Length 2. IMPOSSIBLE (bidirectional).
        #   Wait: T[a][v]... first cycle has v->a, so T[v][a]=1, hence T[a][v]=0.
        #   So a->v is NOT an edge! The walk d=a,v means T[d][v]=T[a][v]=0. But second
        #   cycle needs T[d][v]=1. So if a=d, T[a][v]=0 (from first cycle) and
        #   T[d][v]=T[a][v]=0 (need 1 for second cycle). CONTRADICTION.
        #   So a=d is IMPOSSIBLE.

        # Summary: the only possible overlaps between {a,b} and {c,d} in a tournament are:
        # - a=c (same first vertex of both cycles)
        # - b=d (same last vertex of both cycles)
        # - Both (a=c AND b=d, which means same cycle C1=C2)
        # None of the "cross" overlaps (b=c, a=d) are possible.

        # Let me verify this computationally.
        overlap_types = defaultdict(int)
        for v in range(nn):
            for c1 in d3[v]:
                for c2 in d3[v]:
                    s1 = set(c1)
                    s2 = set(c2)
                    ov = s1 & s2
                    # Check specific overlaps
                    a_eq_c = (c1[0] == c2[0])
                    a_eq_d = (c1[0] == c2[1])
                    b_eq_c = (c1[1] == c2[0])
                    b_eq_d = (c1[1] == c2[1])
                    key = (a_eq_c, a_eq_d, b_eq_c, b_eq_d)
                    overlap_types[key] += 1

        if bits == 0 or (excess > 0 and bits < 20):
            for key, cnt in sorted(overlap_types.items()):
                if cnt > 0:
                    print(f"  bits={bits}: overlap {key}: {cnt}")

    # Now compute exact formula
    # Walk types:
    # Type D (disjoint): no overlap between {a,b} and {c,d}. 1 junction (v). Factor = 6/2 = 3.
    # Type A (a=c): 2 junctions (v and a=c). Factor = 6/4 = 3/2.
    # Type B (b=d): 2 junctions (v and b=d). Factor = 6/4 = 3/2.
    # Type S (same, a=c AND b=d): 3 junctions (v, a=c, b=d). Factor = 6/6 = 1.

    # excess = sum_v [3*D_v + 3/2*A_v + 3/2*B_v + 1*S_v]
    # where D_v, A_v, B_v, S_v are ordered pairs at junction v of each type.

    # But wait: the A-type and B-type walks have 2 junctions, so they're
    # counted TWICE in sum_v. Let me be more careful.

    # Let W be a non-simple walk. W contributes 6 to excess.
    # W has some junction vertices. Each junction v contributes 2 to sum_v t3(v)^2
    # (two rotations starting at v).

    # For the FORMULA: excess = sum_v f(v) where f(v) = contribution of walks
    # whose junction vertex is v, weighted by 6/(number of junctions of walk).
    # This doesn't factor vertex-locally.

    # Better approach: compute excess directly from walk enumeration.
    # excess = 6 * (number of distinct non-simple walk orbits)
    # = (number of non-simple closed walks as sequences of length 6)

    # Since each ordered pair (v, C1, C2) gives one walk as a sequence:
    # excess = sum_v t3(v)^2
    # NO! Because different ordered pairs might give the SAME sequence.
    # E.g., (v, C1, C2) at junction v and the rotated version at junction w
    # give the SAME walk (different starting points).

    # Actually: each ordered pair (v, C1, C2) gives one walk STARTING AT v.
    # Different starting points = different entries in tr(A^6). So they're
    # counted separately in tr(A^6). The total is sum_v t3(v)^2.

    # WAIT. That can't be right because the test showed excess != sum_t3_sq.

    # Let me recount. tr(A^6) = sum of ALL closed walks of length 6 (as sequences
    # starting at each vertex). This equals 6*c_6 + (non-simple walks as sequences).
    # excess = non-simple walks as sequences.

    # Each ordered pair (v, C1, C2) gives ONE closed walk v->c1[0]->c1[1]->v->c2[0]->c2[1]->v.
    # This is a sequence of length 6 starting (and ending) at v.
    # It contributes 1 to tr(A^6)[v][v].

    # Is every non-simple closed walk of this form? YES, because any non-simple
    # walk of length 6 in a tournament decomposes into (3,3) at the first repeated vertex.
    # And the decomposition gives ordered pair (junction, C1, C2).

    # But WHICH vertex is the "first" junction? The walk visits vertices in order
    # w0, w1, ..., w5 with w6=w0. The first vertex that repeats is some w_i = w_j
    # where j is minimal. This gives the decomposition at w_i.

    # For walk v->a->b->v->c->d->v:
    # w0=v, w1=a, w2=b, w3=v. First repeat: w3=w0=v. Junction at v.
    # So the decomposition at the first repeat gives (v, (a,b), (c,d)).
    # Excellent!

    # For the ROTATED walk starting at a: a->b->v->c->d->v->a:
    # w0=a, w1=b, w2=v, w3=c, w4=d, w5=v. First repeat: w5=w2=v. Junction at v.
    # Decomposition: sub-walk w2..w5 = v->c->d->v (length 3), and
    # sub-walk w5..w6,w0..w2 = v->a->b->v (length 3, wrapping around).
    # But this gives the ordered pair (v, (c,d), (a,b)), which is the REVERSE order.
    # So this rotation is ALSO counted in t3(v)^2 as the pair (v, C2, C1).

    # So rotations of a compound walk starting at non-junction vertices
    # also have their first repeat at v (the junction), and decompose as
    # ordered pairs at v. But the ordered pair is DIFFERENT (swapped order).

    # Hmm wait, let me recheck. Walk starting at a: a,b,v,c,d,v.
    # When does the first repeat occur?
    # w0=a, w1=b: new
    # w2=v: new
    # w3=c: new (assuming c != a,b,v)
    # w4=d: new (assuming d != a,b,v,c)
    # w5=v: REPEAT! v first appeared at position 2.
    # So the first repeat is at position 5 where w5=v=w2.
    # Decomposition: sub-walk from position 2 to 5: v->c->d->v (3-cycle).
    # Remaining: position 5 to 2 (wrapping): v->a, a->b->v. Wait, that's
    # positions 5,0,1,2 in the circular walk: v->a->b->v.
    # So the ordered pair is (v, (c,d), (a,b)). YES, this is a DIFFERENT
    # ordered pair than (v, (a,b), (c,d)).

    # So BOTH rotations starting at a AND at the second v generate the
    # pair (v, C2, C1) (the reverse order). But wait, the rotation starting
    # at the second v (position 3):
    # w0=v, w1=c, w2=d, w3=v: REPEAT at position 3 where w3=w0=v.
    # Decomposition: (v, (c,d), ...) and the remaining wraps to (a,b).
    # So pair is (v, (c,d), (a,b)). Same as above!

    # What about rotation starting at b (position 2):
    # b, v, c, d, v, a.
    # w0=b, w1=v: new. w2=c: new. w3=d: new. w4=v: REPEAT at position 1.
    # Decomposition at v: sub-walk positions 1 to 4: v->c->d->v (3-cycle).
    # Remaining: positions 4,5,0,1: v->a->b->v.
    # Pair: (v, (c,d), (a,b)). Same as above!

    # Hmm, so ALL 6 rotations decompose to ordered pairs at v?
    # That seems wrong. Let me check rotation at c (position 4):
    # c, d, v, a, b, v.
    # w0=c, w1=d, w2=v, w3=a, w4=b, w5=v: REPEAT at position 5.
    # Decomposition at v: sub-walk positions 2 to 5: v->a->b->v (3-cycle).
    # Remaining: v->c->d->v.
    # Pair: (v, (a,b), (c,d)). The ORIGINAL order!

    # And rotation at d (position 5):
    # d, v, a, b, v, c.
    # w0=d, w1=v, w2=a, w3=b, w4=v: REPEAT at position 4.
    # Decomposition at v: sub-walk positions 1 to 4: v->a->b->v.
    # Remaining: v->c->d->v.
    # Wait, remaining is: positions 4,5,0,1: v->c->d->v. Wait...
    # position 4 = v, 5 = c. Then wrapping: 0 = d, 1 = v. So v->c->d->v.
    # Pair: (v, (a,b), (c,d)).

    # So the 6 rotations give:
    # Rot 0 (v): (v, (a,b), (c,d))
    # Rot 1 (a): (v, (c,d), (a,b))
    # Rot 2 (b): (v, (c,d), (a,b))
    # Rot 3 (v): (v, (c,d), (a,b))
    # Rot 4 (c): (v, (a,b), (c,d))
    # Rot 5 (d): (v, (a,b), (c,d))

    # So 3 rotations give (v, (a,b), (c,d)) and 3 give (v, (c,d), (a,b)).
    # Each ordered pair is counted 3 times in the walk rotations.
    # But in tr(A^6), each rotation contributes to a DIFFERENT diagonal entry.
    # So tr(A^6) = sum over all sequences = 6 per walk orbit = 3 per ordered pair.

    # THEREFORE: excess = 3 * sum_v t3(v)^2 / DOUBLE_COUNTING_FACTOR.

    # Wait, no. Each ordered pair (v, C1, C2) corresponds to 3 rotations of the
    # walk, each contributing 1 to tr(A^6) at different diagonal entries. But the
    # ordered pair itself is ONE entry in the sum sum_v t3(v)^2.
    # So sum_v t3(v)^2 counts the number of ordered pairs.
    # Each pair corresponds to 3 rotations contributing to tr(A^6).
    # excess = 3 * sum_v t3(v)^2.

    # BUT THE TEST SAID THIS FAILS! Let me recheck my computation of t3(v).

    # The issue: when there are overlap cases, some rotations map to ordered
    # pairs at a DIFFERENT junction vertex (not v). Let me re-examine with overlap.

    # Case: a=c (overlap). Walk: v->a->b->v->a->d->v (using c=a).
    # Vertices: {v, a, b, d} (only 4 distinct vertices).
    # Rotations:
    # Rot 0 (v): v,a,b,v,a,d. First repeat: w3=v=w0. Pair: (v, (a,b), (a,d)).
    # Rot 1 (a): a,b,v,a,d,v. First repeat: w3=a=w0. Pair: (a, (b,v), (d,v)).
    #   So the first repeat is at a, not v! This is decomposed at junction a.
    #   Sub-walk a->b->v->a (3-cycle through a with vertices b,v).
    #   Remaining: a->d->v->a (3-cycle through a with vertices d,v).
    #   Pair: (a, (b,v), (d,v)).
    #   This is counted in t3(a)^2, NOT in t3(v)^2!

    # AH HA! So with overlap, some rotations decompose at DIFFERENT junction vertices.
    # That's why excess != 3 * sum_v t3(v)^2 in general.

    # Actually, each rotation decomposes at SOME junction vertex and contributes
    # to that vertex's t3^2. So excess IS sum_v t3(v)^2 * (correction per vertex),
    # but the correction is uniform if we account for the fact that each rotation
    # contributes EXACTLY 1 to excess and is counted EXACTLY once in some vertex's
    # t3^2 contribution... wait, no. Each rotation is a closed walk. It may
    # have MULTIPLE first-repeat positions. We use the FIRST one.

    # Actually no, the "first repeat" is uniquely defined for each sequence.
    # So each rotation maps to exactly one ordered pair at one junction vertex.
    # Therefore, sum_v (# rotations decomposing at v) = total rotations = excess.
    # And for the rotations decomposing at v: they are SOME ordered pairs (v, C1, C2).
    # But not ALL rotations decomposing at v need be in t3(v)^2.
    # Wait, they DO! If a rotation decomposes at v, then v appears at positions
    # i and i+3 (and the two sub-walks are 3-cycles through v). So the ordered
    # pair IS an element of t3(v)^2.

    # So: excess = sum_v (# rotations whose first repeat is at v)
    #            = sum_v (something <= t3(v)^2 * 3 ???)

    # This is getting very complicated. Let me just compute the EXACT excess
    # by the definition and check against simpler formulas.

    # SIMPLEST CHECK: compute sum_v t3(v)^2 carefully and compare.

    break  # Enough analysis, let me just run a clean test.


# ============================================================
# Clean formula test
# ============================================================
print("\n" + "=" * 70)
print("CLEAN TEST: excess = ? * sum_v t3(v)^2")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    total_tours = 1 << m

    # Test excess = sum_v t3(v)^2
    m1 = 0
    m3 = 0
    for bits in range(total_tours):
        T = tournament_from_bits(n, bits)
        nn = len(T)
        tr6 = matrix_power_trace(T, 6)
        c6 = count_directed_k_cycles_dp(T, 6) if nn >= 6 else 0
        excess = tr6 - 6 * c6

        # t3(v): directed 3-cycles through v
        # A directed 3-cycle through v is an ordered pair (a,b) with v->a->b->v.
        t3v = [0] * nn
        for v in range(nn):
            for a in range(nn):
                if a == v or not T[v][a]:
                    continue
                for b in range(nn):
                    if b == v or b == a or not T[a][b] or not T[b][v]:
                        continue
                    t3v[v] += 1

        sum_t3_sq = sum(t3v[v] ** 2 for v in range(nn))

        if excess == sum_t3_sq:
            m1 += 1
        if excess == 3 * sum_t3_sq:
            m3 += 1

    print(f"  n={n}: excess = 1*sum_t3_sq: {m1}/{total_tours} "
          f"{'PASS' if m1 == total_tours else 'FAIL'}")
    print(f"  n={n}: excess = 3*sum_t3_sq: {m3}/{total_tours} "
          f"{'PASS' if m3 == total_tours else 'FAIL'}")

    if m1 != total_tours and m3 != total_tours:
        # Show some examples
        for bits in range(total_tours):
            T = tournament_from_bits(n, bits)
            nn = len(T)
            tr6 = matrix_power_trace(T, 6)
            c6 = count_directed_k_cycles_dp(T, 6) if nn >= 6 else 0
            excess = tr6 - 6 * c6
            if excess == 0:
                continue
            t3v = [0] * nn
            for v in range(nn):
                for a in range(nn):
                    if a == v or not T[v][a]:
                        continue
                    for b in range(nn):
                        if b == v or b == a or not T[a][b] or not T[b][v]:
                            continue
                        t3v[v] += 1
            sum_t3_sq = sum(t3v[v] ** 2 for v in range(nn))
            if excess != 0 and bits < 50:
                print(f"    bits={bits}: excess={excess}, sum_t3_sq={sum_t3_sq}, "
                      f"ratio={excess/sum_t3_sq:.4f}, t3v={t3v}")

        # Try: sum_v t3(v) * (t3(v) + something)
        # or: products of traces
        print(f"  Trying trace products...")
        m_tr = 0
        for bits in range(total_tours):
            T = tournament_from_bits(n, bits)
            nn = len(T)
            tr6 = matrix_power_trace(T, 6)
            c6 = count_directed_k_cycles_dp(T, 6) if nn >= 6 else 0
            excess = tr6 - 6 * c6
            tr3 = matrix_power_trace(T, 3)
            # tr(A^3)^2 = (3*c3)^2 = 9*c3^2
            # Does excess = tr3^2 / 9 or something?
            # excess involves "pairs of 3-cycles through same vertex"
            # while tr3^2 involves "pairs of 3-cycles (any vertices)"
            if excess * 9 == tr3 ** 2:  # excess = tr3^2/9?
                m_tr += 1
        print(f"    excess = tr3^2/9: {m_tr}/{total_tours}")

        # Actually: tr(A^3)^2 = (sum_v (A^3)[v][v])^2 = sum_{v,w} (A^3)[v][v] * (A^3)[w][w]
        # (A^3)[v][v] = number of closed walks of length 3 starting at v = 2*t3(v)
        # Wait, no. (A^3)[v][v] = sum_{a,b} T[v][a]*T[a][b]*T[b][v]. This counts
        # walks v->a->b->v, which equals t3(v) (directed 3-cycles through v, each
        # counted once as v->a->b->v).
        # So tr(A^3) = sum_v t3(v) = 3*c3 (each 3-cycle counted 3 times).
        # And tr(A^3)^2 = (sum_v t3(v))^2.
        # But sum_v t3(v)^2 != (sum_v t3(v))^2 in general (Cauchy-Schwarz).

        # tr(A^6) is NOT the same as tr(A^3)^2. tr(A^3)^2 = (sum diag entries of A^3)^2,
        # while tr(A^6) = sum diag entries of A^6. And A^6 != (A^3)^2... wait, yes it is!
        # A^6 = (A^3)^2. So tr(A^6) = tr((A^3)^2) = sum_{i,j} (A^3)[i][j]^2.
        # This is the Frobenius norm of A^3, not (tr(A^3))^2.

        # OK so: tr(A^6) = sum_{i,j} (A^3)[i][j]^2 = ||A^3||_F^2.
        # And (A^3)[v][v] = t3(v), so the diagonal part of ||A^3||_F^2 is sum_v t3(v)^2.
        # But tr(A^6) also includes off-diagonal entries: sum_{i!=j} (A^3)[i][j]^2.

        # So: tr(A^6) = sum_v (A^3)[v][v]^2 + sum_{i!=j} (A^3)[i][j]^2
        #             = sum_v t3(v)^2 + sum_{i!=j} (A^3)[i][j]^2

        # And excess = tr(A^6) - 6*c6.
        # Hmm, but I know 6*c6 = 6 * (number of simple 6-cycles), which has nothing
        # to do with t3(v)^2 directly.

        # Wait, I realize I can just check: is excess = sum_v t3(v)^2?
        # Let me recompute with correct t3(v).
        # Actually (A^3)[v][v] = sum_{a,b} T[v][a]*T[a][b]*T[b][v].
        # Each term corresponds to a walk v->a->b->v. If this forms a directed 3-cycle
        # (i.e., {v,a,b} are distinct with v->a->b->v directed cycle), it contributes 1.
        # So (A^3)[v][v] counts directed 3-cycles through v where v is the "start" = t3(v).
        # And t3(v) is exactly what I computed above.

        # So sum_v t3(v)^2 = sum_v (A^3)[v][v]^2.
        # And tr(A^6) = sum_v (A^6)[v][v] = sum_v row_v_of_A3 dot col_v_of_A3
        # = sum_v sum_j (A^3)[v][j] * (A^3)[j][v] = sum_v sum_j (A^3)[v][j]^2
        # Wait no: (A^6)[v][v] = sum_j (A^3)[v][j] * (A^3)[j][v].

        # Hmm, (A^3)[j][v] != (A^3)[v][j] in general (A is not symmetric).
        # So tr(A^6) = sum_v sum_j (A^3)[v][j] * (A^3)[j][v].

        # This is NOT ||A^3||_F^2 = sum_{v,j} (A^3)[v][j]^2.

        # So the relationship is:
        # tr(A^6) = sum_v sum_j (A^3)[v][j] * (A^3)[j][v]
        # excess = tr(A^6) - 6*c6

        # Let me test: excess = sum_v (A^3)[v][v]^2 = sum_v t3(v)^2?
        # That would mean the off-diagonal parts of A^3 contribute only 6*c6.
        # That seems unlikely but let me check.

        # Actually, I already tested this above: m1 counts excess = sum_t3_sq.
        # It failed. So excess != sum_v t3(v)^2.

        # Let me instead test: excess + 6*c6 = tr(A^6) = sum_{v,j} (A^3)[v][j]*(A^3)[j][v]
        # and see if we can decompose this.

        # Actually, let me look at small examples.
        print(f"\n  Detailed n={n} examples:")
        ct = 0
        for bits in range(total_tours):
            T = tournament_from_bits(n, bits)
            nn = len(T)
            tr6 = matrix_power_trace(T, 6)
            c6 = count_directed_k_cycles_dp(T, 6) if nn >= 6 else 0
            c3 = count_directed_k_cycles_dp(T, 3)
            excess = tr6 - 6 * c6
            t3v = [0] * nn
            for v in range(nn):
                for a in range(nn):
                    if a == v or not T[v][a]:
                        continue
                    for b in range(nn):
                        if b == v or b == a or not T[a][b] or not T[b][v]:
                            continue
                        t3v[v] += 1
            sum_t3_sq = sum(t3v[v] ** 2 for v in range(nn))
            if excess > 0 and ct < 5:
                print(f"    bits={bits}: c3={c3}, c6={c6}, excess={excess}, "
                      f"sum_t3_sq={sum_t3_sq}, t3v={t3v}")
                ct += 1

print("\nDone.")
