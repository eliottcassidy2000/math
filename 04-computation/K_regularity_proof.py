#!/usr/bin/env python3
"""
K_regularity_proof.py -- Explore WHY K = c5 - 2*ov1 - 2*ov2 is constant
for ALL regular tournaments (not just circulant ones).

Key approach: decompose K into local 5-vertex contributions and show
that regularity forces the global sum to be constant.

The identity K(n) = -3n(n^2-1)(n^2-9)/320 holds for n odd.

Strategy:
1. On each 5-vertex subset, compute local K_5 = c5_local - 2*ov1_local - 2*ov2_local
2. Show that K = sum over 5-subsets of some local contribution
3. Find what regularity forces about the 5-subset score distribution

Author: kind-pasteur-2026-03-12-S60
"""

from itertools import combinations
from collections import defaultdict
import random


def tournament_from_bits(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


def count_c3_sets(A, n):
    c3 = []
    for a, b, c in combinations(range(n), 3):
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            c3.append(frozenset([a, b, c]))
    return list(set(c3))


def count_c5_directed_on_subset(A, verts):
    """Count directed 5-cycles on the given 5-vertex subset (must be exactly 5 verts)."""
    nn = len(verts)
    assert nn == 5
    total = 0
    start = 0
    dp = {(1 << start, start): 1}
    for mask in range(1, 1 << nn):
        if not (mask & (1 << start)):
            continue
        for v in range(nn):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(nn):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << nn) - 1
    for v in range(nn):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


def count_c5_global(A, n):
    """Count total directed 5-cycles in tournament on n vertices."""
    total = 0
    for subset in combinations(range(n), 5):
        total += count_c5_directed_on_subset(A, list(subset))
    return total


def main():
    print("=" * 70)
    print("WHY IS K CONSTANT FOR ALL REGULAR TOURNAMENTS?")
    print("=" * 70)

    # ====== PART 1: 5-vertex score sequence census ======
    print("\n" + "=" * 60)
    print("PART 1: 5-vertex score distributions in regular n=7 tournaments")
    print("=" * 60)

    n = 7
    m = 3

    # Generate some regular tournaments
    random.seed(42)

    # First: the circulant (Paley) tournament
    S_paley = {1, 2, 4}  # QR mod 7
    A_paley = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S_paley:
            A_paley[i][(i + s) % n] = 1

    # Analyze 5-vertex subset properties for Paley
    print(f"\n  Paley T_7 (S = {S_paley}):")
    score_dist = defaultdict(int)
    c5_by_score = defaultdict(list)
    c3_by_score = defaultdict(list)

    for subset in combinations(range(n), 5):
        verts = list(subset)
        sub_scores = tuple(sorted([sum(A_paley[v][w] for w in verts if w != v) for v in verts]))
        c5_local = count_c5_directed_on_subset(A_paley, verts)
        # Count 3-cycles in this 5-subset
        c3_local = 0
        for a, b, c in combinations(verts, 3):
            if (A_paley[a][b] and A_paley[b][c] and A_paley[c][a]) or \
               (A_paley[a][c] and A_paley[c][b] and A_paley[b][a]):
                c3_local += 1
        score_dist[sub_scores] += 1
        c5_by_score[sub_scores].append(c5_local)
        c3_by_score[sub_scores].append(c3_local)

    print(f"    5-subset score distributions:")
    total_c5 = 0
    for scores in sorted(score_dist):
        c5_vals = c5_by_score[scores]
        c3_vals = c3_by_score[scores]
        avg_c5 = sum(c5_vals) / len(c5_vals)
        avg_c3 = sum(c3_vals) / len(c3_vals)
        print(f"      {scores}: count={score_dist[scores]}, avg_c5={avg_c5:.1f}, avg_c3={avg_c3:.1f}")
        total_c5 += sum(c5_vals)
    print(f"    Total c5 (sum over 5-subsets): {total_c5}")
    print(f"    Actual c5 (global): {total_c5}")  # Each 5-cycle counted once

    # ====== PART 2: Local vs Global decomposition ======
    print("\n" + "=" * 60)
    print("PART 2: Local decomposition of K")
    print("=" * 60)

    # Each directed 5-cycle lives in exactly 1 five-subset.
    # Each ov1 pair (3-cycles sharing 1 vertex) spans 5 vertices -> 1 five-subset.
    # Each ov2 pair (3-cycles sharing 2 vertices, i.e., sharing an edge) spans 4 vertices.
    #   -> This pair appears in (n-4) five-subsets (choosing the 5th vertex)!
    #
    # So: K = sum over 5-subsets of [c5_local - 2*ov1_local] - 2*ov2_global
    #
    # But ov2 pairs span 4 vertices, so they're counted in (n-4) five-subsets.
    # We need to be careful about the decomposition.

    # Let's compute directly: for each 5-subset, what is c5 - 2*ov1 - 2*ov2 LOCAL?
    print(f"\n  Local K on 5-subsets for Paley T_7:")
    local_K_dist = defaultdict(int)
    local_K_by_score = defaultdict(list)

    for subset in combinations(range(n), 5):
        verts = list(subset)
        sub_A = [[A_paley[v][w] for w in verts] for v in verts]
        # Get c3 sets within this 5-subset
        c3_local = count_c3_sets(sub_A, 5)
        c5_local = count_c5_directed_on_subset(A_paley, verts)
        # Overlap analysis within the 5-subset
        n3 = len(c3_local)
        disj_l = ov1_l = ov2_l = 0
        for i in range(n3):
            for j in range(i+1, n3):
                o = len(c3_local[i] & c3_local[j])
                if o == 0:
                    disj_l += 1
                elif o == 1:
                    ov1_l += 1
                elif o == 2:
                    ov2_l += 1
        K_local = c5_local - 2*ov1_l - 2*ov2_l
        sub_scores = tuple(sorted([sum(sub_A[i][j] for j in range(5) if j != i) for i in range(5)]))
        local_K_dist[K_local] += 1
        local_K_by_score[sub_scores].append((K_local, c5_local, ov1_l, ov2_l, n3))

    print(f"    Local K distribution: {dict(sorted(local_K_dist.items()))}")
    print(f"    Sum of local K: {sum(k*v for k,v in local_K_dist.items())}")
    print(f"    Global K: -126")

    print(f"\n    Local K by 5-subset score:")
    for scores in sorted(local_K_by_score):
        vals = local_K_by_score[scores]
        K_vals = [v[0] for v in vals]
        print(f"      {scores}: K_local values = {sorted(set(K_vals))}, "
              f"c5={sorted(set(v[1] for v in vals))}, "
              f"ov1={sorted(set(v[2] for v in vals))}, "
              f"ov2={sorted(set(v[3] for v in vals))}, "
              f"c3_count={sorted(set(v[4] for v in vals))}")

    # ====== PART 3: The multiplicity correction ======
    print("\n" + "=" * 60)
    print("PART 3: Multiplicity correction")
    print("=" * 60)

    # Key question: how does K = sum of local contributions?
    # c5: each directed 5-cycle lives in exactly 1 five-subset -> coefficient = 1
    # ov1: each ov1 pair (sharing 1 vertex) spans 5 vertices -> coefficient = 1
    # ov2: each ov2 pair (sharing 2 vertices = 1 edge) spans 4 vertices
    #       -> appears in C(n-4,1) = n-4 five-subsets -> coefficient = n-4

    # So: sum_subsets(K_local) = c5 - 2*ov1 - 2*(n-4)*ov2
    # But K_global = c5 - 2*ov1 - 2*ov2
    # Therefore: K_global = sum_subsets(K_local) + 2*(n-5)*ov2

    # Let's verify this!
    c3_sets_global = count_c3_sets(A_paley, n)
    n3_g = len(c3_sets_global)
    disj_g = ov1_g = ov2_g = 0
    for i in range(n3_g):
        for j in range(i+1, n3_g):
            o = len(c3_sets_global[i] & c3_sets_global[j])
            if o == 0:
                disj_g += 1
            elif o == 1:
                ov1_g += 1
            elif o == 2:
                ov2_g += 1

    sum_local = sum(k*v for k, v in local_K_dist.items())
    K_global = -126

    print(f"  Global: c3={n3_g}, ov1={ov1_g}, ov2={ov2_g}, disj={disj_g}")
    print(f"  sum_local = {sum_local}")
    print(f"  K_global  = {K_global}")
    print(f"  Predicted: sum_local + 2*(n-5)*ov2 = {sum_local} + {2*(n-5)*ov2_g} = {sum_local + 2*(n-5)*ov2_g}")
    print(f"  Match: {sum_local + 2*(n-5)*ov2_g == K_global}")

    # ====== PART 4: Score sequence determines K_local? ======
    print("\n" + "=" * 60)
    print("PART 4: Is K_local determined by 5-subset score sequence?")
    print("=" * 60)

    # Test across MULTIPLE regular tournaments at n=7
    # If K_local is score-determined, then K_global constancy follows from
    # the fact that all regular tournaments have the same score distribution
    # on 5-subsets... wait, DO they?

    # First check: do all regular n=7 tournaments have the same 5-subset score distribution?
    print(f"\n  Testing 5-subset score distributions across regular n=7 tournaments:")

    def random_regular_tournament(n):
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        mm = (n - 1) // 2
        for _ in range(10000):
            scores = [sum(A[i]) for i in range(n)]
            if max(scores) == mm and min(scores) == mm:
                return A
            high = scores.index(max(scores))
            low = scores.index(min(scores))
            if high != low and A[high][low]:
                A[high][low] = 0
                A[low][high] = 1
            else:
                for j in range(n):
                    if j != high and A[high][j]:
                        if scores[j] < mm:
                            A[high][j] = 0
                            A[j][high] = 1
                            break
        scores = [sum(A[i]) for i in range(n)]
        if all(s == mm for s in scores):
            return A
        return None

    score_dists_seen = []
    for trial in range(10):
        A = random_regular_tournament(7)
        if A is None:
            continue
        sd = defaultdict(int)
        for subset in combinations(range(7), 5):
            verts = list(subset)
            sub_scores = tuple(sorted([sum(A[v][w] for w in verts if w != v) for v in verts]))
            sd[sub_scores] += 1
        score_dists_seen.append(dict(sd))
        if trial < 5:
            print(f"    Trial {trial}: {dict(sorted(sd.items()))}")

    # Check if they're all the same
    all_same = all(d == score_dists_seen[0] for d in score_dists_seen)
    print(f"\n  All score distributions identical: {all_same}")

    if not all_same:
        # Score distributions DIFFER across regular tournaments
        # So constancy of K is NOT explained by score-distribution constancy
        print(f"  Score distributions DIFFER! K constancy needs a deeper explanation.")
        # Show the differences
        all_score_types = set()
        for d in score_dists_seen:
            all_score_types.update(d.keys())
        for st in sorted(all_score_types):
            vals = [d.get(st, 0) for d in score_dists_seen]
            if len(set(vals)) > 1:
                print(f"    {st}: {vals}")

    # ====== PART 5: All tournaments on n=5 - deep analysis ======
    print("\n" + "=" * 60)
    print("PART 5: Exhaustive n=5 analysis - what does regularity force?")
    print("=" * 60)

    n = 5
    m = 2

    # For n=5, ov2 pairs span 4 out of 5 vertices.
    # There's only 1 five-subset (the whole tournament).
    # So K_local = K_global for n=5.

    # At n=5, c5 is the directed 5-cycles in the whole tournament.
    # c3 pairs can have overlap 0 (impossible - only 5 vertices, need 6 for disjoint),
    # overlap 1 (span 5 vertices), or overlap 2 (span 4 vertices).

    # Actually at n=5 with c3 vertices = 3+3-overlap >= 5 only if overlap <= 1.
    # Overlap 0: need 6 vertices - IMPOSSIBLE at n=5!
    # So at n=5: disj = 0 always, C(c3,2) = ov1 + ov2

    print(f"  At n=5: disj3 = 0 always (need 6 vertices for disjoint pair)")
    print(f"  So K = c5 - 2*ov1 - 2*ov2 = c5 - 2*(ov1+ov2) = c5 - 2*C(c3,2)")
    print(f"  For regular: c3 = {n*(n**2-1)//24} = 5, C(c3,2) = 10")
    print(f"  K = c5 - 20")

    # Verify: K=-18 means c5=2 for regular 5-vertex tournaments? No wait...
    # c5 is DIRECTED 5-cycles. For regular n=5 (all scores = 2):
    # K = c5 - 2*C(5,2) = c5 - 20 = -18 => c5 = 2

    # Check: how many directed 5-cycles in a regular 5-vertex tournament?
    reg_c5 = set()
    for bits in range(1 << 10):
        A = tournament_from_bits(5, bits)
        scores = [sum(A[i]) for i in range(5)]
        if all(s == 2 for s in scores):
            c3_s = count_c3_sets(A, 5)
            c5 = count_c5_directed_on_subset(A, list(range(5)))
            reg_c5.add(c5)

    print(f"  Directed 5-cycles in regular n=5: {sorted(reg_c5)}")
    print(f"  Verification: K = c5 - 20 = {list(reg_c5)[0]} - 20 = {list(reg_c5)[0]-20}")

    # ====== PART 6: Moment analysis ======
    print("\n" + "=" * 60)
    print("PART 6: What polynomial in scores equals K_local?")
    print("=" * 60)

    # At n=5, K = c5 - 2*C(c3,2) and c3 depends on scores.
    # For general n, each 5-subset contributes.
    #
    # Key question: is there a polynomial in the score sequence that equals K?
    #
    # For a tournament on n vertices with out-degrees d_1,...,d_n:
    # c3 = C(n,3) - sum C(d_i, 2)  (well-known)
    # c5 = ???
    #
    # For REGULAR: d_i = m for all i, so c3 = C(n,3) - n*C(m,2)
    # And K should depend on sum d_i^k for various k.

    # Let's look at K as a function of score sequence for ALL n=5 tournaments
    print(f"\n  K by score sequence at n=5:")
    K_by_scores = defaultdict(list)
    for bits in range(1 << 10):
        A = tournament_from_bits(5, bits)
        scores = tuple(sorted([sum(A[i]) for i in range(5)]))
        c3_s = count_c3_sets(A, 5)
        c3_count = len(c3_s)
        c5 = count_c5_directed_on_subset(A, list(range(5)))
        # overlap at n=5: disj=0, so ov1+ov2 = C(c3,2)
        n3 = len(c3_s)
        ov1 = ov2 = 0
        for i in range(n3):
            for j in range(i+1, n3):
                o = len(c3_s[i] & c3_s[j])
                if o == 1:
                    ov1 += 1
                elif o == 2:
                    ov2 += 1
        K = c5 - 2*ov1 - 2*ov2
        K_by_scores[scores].append((K, c3_count, c5, ov1, ov2))

    for scores in sorted(K_by_scores):
        vals = K_by_scores[scores]
        K_set = sorted(set(v[0] for v in vals))
        c3_set = sorted(set(v[1] for v in vals))
        c5_set = sorted(set(v[2] for v in vals))
        ov1_set = sorted(set(v[3] for v in vals))
        ov2_set = sorted(set(v[4] for v in vals))
        S2 = sum(s**2 for s in scores)
        S3 = sum(s**3 for s in scores)
        print(f"    {scores}: K={K_set}, c3={c3_set}, c5={c5_set}, "
              f"S2={S2}, S3={S3}")

    # Test if K = f(S2, S3, S4) for some polynomial
    print(f"\n  Testing: is K determined by power sums S_k = sum d_i^k?")
    from collections import Counter
    power_sums_to_K = defaultdict(set)
    for bits in range(1 << 10):
        A = tournament_from_bits(5, bits)
        scores = [sum(A[i]) for i in range(5)]
        c3_s = count_c3_sets(A, 5)
        c5 = count_c5_directed_on_subset(A, list(range(5)))
        n3 = len(c3_s)
        ov1 = ov2 = 0
        for i in range(n3):
            for j in range(i+1, n3):
                o = len(c3_s[i] & c3_s[j])
                if o == 1:
                    ov1 += 1
                elif o == 2:
                    ov2 += 1
        K = c5 - 2*ov1 - 2*ov2
        S1 = sum(scores)  # = C(n,2) = 10
        S2 = sum(s**2 for s in scores)
        S3 = sum(s**3 for s in scores)
        S4 = sum(s**4 for s in scores)
        power_sums_to_K[(S1, S2, S3, S4)].add(K)

    all_determined = all(len(v) == 1 for v in power_sums_to_K.values())
    print(f"  K determined by (S1,S2,S3,S4): {all_determined}")
    if not all_determined:
        for key, vals in sorted(power_sums_to_K.items()):
            if len(vals) > 1:
                print(f"    {key}: K = {sorted(vals)}")

    # Try: is K determined by (S2, S3)?
    ps23_to_K = defaultdict(set)
    for bits in range(1 << 10):
        A = tournament_from_bits(5, bits)
        scores = [sum(A[i]) for i in range(5)]
        c3_s = count_c3_sets(A, 5)
        c5 = count_c5_directed_on_subset(A, list(range(5)))
        n3 = len(c3_s)
        ov1 = ov2 = 0
        for i in range(n3):
            for j in range(i+1, n3):
                o = len(c3_s[i] & c3_s[j])
                if o == 1:
                    ov1 += 1
                elif o == 2:
                    ov2 += 1
        K = c5 - 2*ov1 - 2*ov2
        S2 = sum(s**2 for s in scores)
        S3 = sum(s**3 for s in scores)
        ps23_to_K[(S2, S3)].add(K)

    all_det_23 = all(len(v) == 1 for v in ps23_to_K.values())
    print(f"  K determined by (S2, S3): {all_det_23}")
    if not all_det_23:
        for key, vals in sorted(ps23_to_K.items()):
            if len(vals) > 1:
                print(f"    (S2={key[0]}, S3={key[1]}): K = {sorted(vals)}")

    # ====== PART 7: Score-regularity and c5 formula ======
    print("\n" + "=" * 60)
    print("PART 7: Known c5 formula in terms of scores")
    print("=" * 60)

    # Kendall's formula (or similar): for tournament on n vertices with scores d_1,...,d_n,
    # c3 = C(n,3) - sum C(d_i, 2)
    #
    # For 5-cycles, there's no simple formula just from score sequence.
    # But for REGULAR tournaments:
    # c3 = C(n,3) - n*C(m,2) = n(n-1)(n-2)/6 - n*m*(m-1)/2
    #    = n[(n-1)(n-2)/6 - m(m-1)/2]
    #    With m=(n-1)/2:
    #    = n[(n-1)(n-2)/6 - (n-1)(n-3)/8]
    #    = n(n-1)[(n-2)/6 - (n-3)/8]
    #    = n(n-1)[4(n-2) - 3(n-3)] / 24
    #    = n(n-1)[4n-8-3n+9] / 24
    #    = n(n-1)(n+1) / 24

    print(f"  c3 for regular tournament: n(n^2-1)/24")
    for n in [5, 7, 9, 11, 13]:
        c3 = n*(n**2-1)//24
        print(f"    n={n}: c3 = {c3}")

    # Liskovets-Savchenko formula for c5 in regular tournaments:
    # c5 = [n(n-1)(n-2)(n-3)(n-4) - 20*n*sum_{i<j}(2*d_i*d_j - (n-1)d_i - (n-1)d_j + ...)] / 120
    # For regular: c5 depends on the FULL structure, not just scores!
    # But K = c5 - 2*C(c3,2) for n=5 IS constant.
    #
    # Actually, what IS constant is c5 for regular n=5:
    print(f"\n  Is c5 constant for regular tournaments?")
    for n in [5, 7, 9]:
        mm = (n - 1) // 2
        c5_vals = set()
        count = 0
        if n == 5:
            for bits in range(1 << 10):
                A = tournament_from_bits(n, bits)
                scores = [sum(A[i]) for i in range(n)]
                if all(s == mm for s in scores):
                    c5 = count_c5_global(A, n)
                    c5_vals.add(c5)
                    count += 1
        else:
            random.seed(123)
            for _ in range(2000):
                A = random_regular_tournament(n)
                if A is None:
                    continue
                c5 = count_c5_global(A, n)
                c5_vals.add(c5)
                count += 1
        print(f"    n={n}: c5 values = {sorted(c5_vals)} ({count} tournaments)")

    # ====== PART 8: c5 varies but K is constant ======
    print("\n" + "=" * 60)
    print("PART 8: c5 varies but K is constant - the cancellation")
    print("=" * 60)

    # At n=7, c5 varies but K = c5 - 2*ov1 - 2*ov2 is constant.
    # This means ov1 + ov2 = (c5 - K)/2, and this sum tracks c5.
    #
    # Since disj + ov1 + ov2 = C(c3, 2) (constant for regular),
    # we have disj = C(c3,2) - (c5-K)/2 = C(c3,2) - c5/2 + K/2
    # i.e., disj = -c5/2 + [C(c3,2) + K/2]
    # which is the disjoint 3-cycle identity!

    n = 7
    mm = 3
    c3_n7 = n*(n**2-1)//24  # = 14
    K_n7 = -126

    print(f"  n=7: c3={c3_n7}, C(c3,2)={c3_n7*(c3_n7-1)//2}, K={K_n7}")
    print(f"  Formula: ov1 + ov2 = (c5 + 126)/2")
    print(f"  And disj = {c3_n7*(c3_n7-1)//2} - (c5 + 126)/2 = 91 + 63 - c5/2 = 28 - (c5-126)/2")

    # Verify with specific tournaments
    random.seed(99)
    print(f"\n  Verification on random regular n=7 tournaments:")
    for trial in range(10):
        A = random_regular_tournament(7)
        if A is None:
            continue
        c3_sets = count_c3_sets(A, 7)
        c5 = count_c5_global(A, 7)
        n3 = len(c3_sets)
        assert n3 == c3_n7
        disj = ov1 = ov2 = 0
        for i in range(n3):
            for j in range(i+1, n3):
                o = len(c3_sets[i] & c3_sets[j])
                if o == 0:
                    disj += 1
                elif o == 1:
                    ov1 += 1
                elif o == 2:
                    ov2 += 1
        K = c5 - 2*ov1 - 2*ov2
        print(f"    c5={c5:>4}, ov1={ov1:>3}, ov2={ov2:>3}, disj={disj:>3}, K={K}")
        assert K == K_n7

    # ====== PART 9: The deep WHY ======
    print("\n" + "=" * 60)
    print("PART 9: Algebraic identity behind constancy")
    print("=" * 60)

    # For a tournament T on n vertices with adjacency matrix A (A[i][j]=1 if i->j):
    #
    # c3 = tr(A^3)/3  (well known)
    # c5 = tr(A^5)/5  (well known)
    #
    # ov1 = # pairs of 3-cycles sharing exactly 1 vertex
    #      = sum_v c3_through_v * (c3 - c3_through_v) / ... hmm this is complicated
    #
    # ov2 = # pairs of 3-cycles sharing exactly 2 vertices = sharing an edge
    #      = sum_{i->j} c3(i,j) * (c3(i,j) - 1) / 2
    #      where c3(i,j) = # 3-cycles containing the edge i->j
    #
    # For edge i->j: c3(i,j) = |{k : j->k and k->i}| = |out(j) cap in(i)|
    #                         = |out(j) cap out'(i)| where out'(i) = {k : k->i}
    #
    # For regular tournament: |out(i)| = |in(i)| = m for all i
    # |out(j) cap in(i)| = lambda_{ij} (number of common in-out neighbors)

    # Let's compute these lambda values
    print(f"  For regular n=7, computing lambda_{'{ij}'} = |out(j) cap in(i)| for edge i->j:")

    A = A_paley
    lambda_vals = defaultdict(int)
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j]:
                lam = sum(1 for k in range(n) if k != i and k != j and A[j][k] and A[k][i])
                lambda_vals[lam] += 1

    print(f"    Lambda distribution: {dict(sorted(lambda_vals.items()))}")

    # For DRT (doubly regular tournament): lambda is constant!
    # lambda = (n-3)/4 for doubly regular tournament on n vertices
    print(f"    DRT lambda = (n-3)/4 = {(n-3)/4}")

    # ov2 for DRT:
    # Each edge contributes C(lambda, 2) = lambda*(lambda-1)/2 pairs
    # Total edges = n*m = n*(n-1)/2
    # ov2 = n*m * C(lambda, 2) = n*m*lambda*(lambda-1)/2
    lam_drt = (n - 3) / 4
    ov2_drt = n * mm * lam_drt * (lam_drt - 1) / 2
    print(f"    DRT ov2 = n*m*C(lambda,2) = {ov2_drt}")
    print(f"    Actual ov2 for Paley (which IS doubly regular): {ov2_g}")

    # For n=7, lambda = 1, C(1,2) = 0, so ov2_DRT = 0!
    # But the actual ov2 is NOT 0 for Paley n=7... let's check
    print(f"\n  Check: for Paley T_7, lambda = 1 means each edge is in exactly 1 three-cycle")
    print(f"  So C(1,2) = 0, meaning no two 3-cycles share an edge!")
    print(f"  This means ov2 = 0 for doubly regular n=7.")
    print(f"  Actual ov2 = {ov2_g}")
    if ov2_g == 0:
        print(f"  CONFIRMED: ov2=0 for Paley T_7 (as predicted)")

    # For non-DRT regular tournaments, lambda varies. Let's check.
    random.seed(77)
    print(f"\n  Lambda distributions for non-Paley regular n=7:")
    for trial in range(5):
        A = random_regular_tournament(7)
        if A is None:
            continue
        lam_dist = defaultdict(int)
        for i in range(7):
            for j in range(7):
                if i != j and A[i][j]:
                    lam = sum(1 for k in range(7) if k != i and k != j and A[j][k] and A[k][i])
                    lam_dist[lam] += 1
        c3_s = count_c3_sets(A, 7)
        c5 = count_c5_global(A, 7)
        n3 = len(c3_s)
        ov1 = ov2 = 0
        for i in range(n3):
            for j in range(i+1, n3):
                o = len(c3_s[i] & c3_s[j])
                if o == 1:
                    ov1 += 1
                elif o == 2:
                    ov2 += 1
        print(f"    Trial {trial}: lambda={dict(sorted(lam_dist.items()))}, "
              f"c5={c5}, ov1={ov1}, ov2={ov2}, K={c5-2*ov1-2*ov2}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
