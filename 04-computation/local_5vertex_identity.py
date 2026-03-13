#!/usr/bin/env python3
"""
local_5vertex_identity.py -- Derive the 5-vertex identity behind disj3 = -c5/2 + const

KEY INSIGHT: On 5 vertices, disjoint 3-cycle pairs are IMPOSSIBLE (need 6 vertices).
So within each 5-subset, ALL 3-cycle pairs overlap (share 1 or 2 vertices).

For the global identity c5 - 2*ov1 - 2*ov2 = K(p):
  - Each global ov1 pair lives in exactly 1 five-subset
  - Each global ov2 pair lives in exactly (p-4) five-subsets
  - Each c5 cycle lives in exactly 1 five-subset

So: K = c5 - 2*ov1 - 2*ov2 = sum_V [c5(V) - 2*ov1(V)] - 2*ov2
    = sum_V [c5(V) - 2*ov1(V)] - 2/(p-4) * sum_V ov2(V)

For this to be orientation-independent, we need:
  sum_V [c5(V) - 2*ov1(V) - 2/(p-4)*ov2(V)] = K(p)

Since on V: C(c3(V),2) = ov1(V) + ov2(V):
  sum_V [c5(V) - 2*C(c3(V),2) + 2*(1-1/(p-4))*ov2(V)] = K(p)

This script:
1. For EVERY 5-vertex tournament type, tabulate c3, c5, ov1, ov2
2. Check if c5 - 2*ov1 = f(c3) for some universal function
3. If not, explore what residuals exist and how they cancel
4. Derive K(p) from the combinatorial structure

Author: kind-pasteur-2026-03-12-S60
"""

import math
from itertools import combinations, permutations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def all_5_tournaments():
    """Enumerate ALL non-isomorphic 5-vertex tournaments.
    Returns list of (adjacency, score_seq, count_in_class)."""
    n = 5
    seen_canonical = set()
    results = []

    for bits in range(1 << (n*(n-1)//2)):
        A = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        # Canonical form: sort by image under all permutations
        # Use score sequence + secondary invariants
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))

        # For canonical form, generate all relabelings and take min
        min_rep = None
        for perm in permutations(range(n)):
            rep = []
            for i in range(n):
                for j in range(i+1, n):
                    rep.append(A[perm[i]][perm[j]])
            rep = tuple(rep)
            if min_rep is None or rep < min_rep:
                min_rep = rep
        if min_rep in seen_canonical:
            continue
        seen_canonical.add(min_rep)

        # Count isomorphism class size
        # (actually we'll just count 1 for now)

        # Compute c3, c5, ov1, ov2
        c3_sets = []
        for a, b, c in combinations(range(n), 3):
            if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                c3_sets.append(frozenset([a, b, c]))
        c3 = len(c3_sets)

        start = 0
        dp = {(1 << start, start): 1}
        for mask in range(1, 1 << n):
            if not (mask & (1 << start)):
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                key = (mask, v)
                if key not in dp or dp[key] == 0:
                    continue
                cnt = dp[key]
                for w in range(n):
                    if mask & (1 << w):
                        continue
                    if A[v][w]:
                        nkey = (mask | (1 << w), w)
                        dp[nkey] = dp.get(nkey, 0) + cnt
        full = (1 << n) - 1
        c5 = 0
        for v in range(n):
            if v == start:
                continue
            key = (full, v)
            if key in dp and dp[key] > 0:
                if A[v][start]:
                    c5 += dp[key]

        ov1 = 0
        ov2 = 0
        for i in range(len(c3_sets)):
            for j in range(i+1, len(c3_sets)):
                overlap = len(c3_sets[i] & c3_sets[j])
                if overlap == 1:
                    ov1 += 1
                elif overlap == 2:
                    ov2 += 1

        results.append({
            'A': A, 'scores': scores, 'c3': c3, 'c5': c5,
            'ov1': ov1, 'ov2': ov2
        })

    return results


def all_orientations_circ(p):
    m = (p - 1) // 2
    orientations = []
    for bits in range(1 << m):
        S = []
        for i in range(m):
            chord = i + 1
            if (bits >> i) & 1:
                S.append(chord)
            else:
                S.append(p - chord)
        orientations.append((bits, S))
    return orientations


def count_c3_vertex_sets(A, p):
    c3_sets = []
    for a, b, c in combinations(range(p), 3):
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            c3_sets.append(frozenset([a, b, c]))
    return list(set(c3_sets))


def main():
    print("=" * 70)
    print("LOCAL 5-VERTEX IDENTITY ANALYSIS")
    print("=" * 70)

    # ====== PART 1: Complete 5-vertex tournament census ======
    print("\n" + "=" * 60)
    print("PART 1: Non-isomorphic 5-vertex tournaments")
    print("=" * 60)

    tourns = all_5_tournaments()
    print(f"\n  {len(tourns)} non-isomorphic 5-vertex tournaments")

    print(f"\n  {'scores':>18} {'c3':>4} {'c5':>4} {'ov1':>4} {'ov2':>4} "
          f"{'c5-2ov1':>8} {'c5-2(ov1+ov2)':>14}")

    for t in sorted(tourns, key=lambda x: (x['c3'], x['c5'])):
        combo1 = t['c5'] - 2 * t['ov1']
        combo2 = t['c5'] - 2 * (t['ov1'] + t['ov2'])
        print(f"  {str(t['scores']):>18} {t['c3']:>4} {t['c5']:>4} {t['ov1']:>4} {t['ov2']:>4} "
              f"{combo1:>8} {combo2:>14}")

    # Group by c3
    by_c3 = defaultdict(list)
    for t in tourns:
        by_c3[t['c3']].append(t)

    print(f"\n  Grouped by c3:")
    for c3 in sorted(by_c3):
        group = by_c3[c3]
        c5_vals = sorted(set(t['c5'] for t in group))
        ov1_vals = sorted(set(t['ov1'] for t in group))
        combo_vals = sorted(set(t['c5'] - 2*t['ov1'] for t in group))
        combo2_vals = sorted(set(t['c5'] - 2*(t['ov1']+t['ov2']) for t in group))
        print(f"    c3={c3}: c5={c5_vals}, c5-2*ov1={combo_vals}, c5-2*(ov1+ov2)={combo2_vals}")

    # ====== PART 2: What IS constant on 5 vertices? ======
    print("\n" + "=" * 60)
    print("PART 2: Searching for universal 5-vertex identity")
    print("=" * 60)

    # Try linear combinations: a*c5 + b*ov1 + c*ov2 = f(c3)
    # We need this to hold for ALL tournaments, not just circulant

    # For c3=3, we have two entries with different (c5, ov1, ov2):
    # (1, 0, 3) and (1, 1, 2)
    # For these: any a*1 + b*0 + c*3 = a*1 + b*1 + c*2
    # => b = c (from same c5=1, 3c-b = 2c+b => b = c)
    # With b=c: a + 3b = a + b + 2b = a + 3b (trivially same)

    # For c3=4, we have three entries: (1,2,4), (2,3,3), (3,3,3)
    # a + 2b + 4c = a + 2b + 4c (entry 1)
    # 2a + 3b + 3c (entry 2)
    # 3a + 3b + 3c (entry 3)
    # These must be equal:
    # a + 2b + 4c = 2a + 3b + 3c => -a - b + c = 0 => c = a + b
    # 2a + 3b + 3c = 3a + 3b + 3c => a = 0
    # So a = 0, c = b.
    # Then b*ov1 + b*ov2 = b*(ov1+ov2) = b*C(c3,2) must equal f(c3)
    # This is trivially true: ov1+ov2 = C(c3,2) on 5 vertices (no disjoint pairs)

    # So on 5 vertices, the only linear identity is ov1+ov2 = C(c3,2)
    # c5 is NOT linearly determined by c3 alone (or by ov1,ov2)

    print("  On 5 vertices, the ONLY linear identity is ov1 + ov2 = C(c3,2)")
    print("  c5 is NOT a function of c3 (or of ov1,ov2) alone")

    # So the GLOBAL identity c5 - 2*ov1 - 2*ov2 = K(p) does NOT reduce to a LOCAL 5-vertex identity!
    # It's a GLOBAL property of REGULAR tournaments.

    # ====== PART 3: What makes regularity special? ======
    print("\n" + "=" * 60)
    print("PART 3: Regularity and the 5-vertex distribution")
    print("=" * 60)

    # For circulant T on Z_p: each 5-subset V induces a sub-tournament
    # The DISTRIBUTION of 5-vertex tournament types depends on the orientation
    # But the weighted sum c5(V) - 2*ov1(V) turns out constant

    # The question is: what DISTRIBUTION property ensures constancy?
    # Key: for regular tournaments, the number of 5-subsets with each score sequence
    # is determined. But the sub-tournament type (which depends on arc directions,
    # not just scores) can vary.

    # Let's check: for two orientations at same p, are the score sequence distributions
    # of 5-subsets the same?

    for p in [7, 11]:
        m = (p - 1) // 2

        print(f"\n  p={p}:")
        dist_by_orient = {}

        for bits, S in all_orientations_circ(p):
            A = build_adj(p, S)

            # Score distribution of 5-subsets
            score_dist = defaultdict(int)
            # Type distribution
            type_dist = defaultdict(int)

            for subset in combinations(range(p), 5):
                verts = list(subset)
                # Sub-tournament scores
                sub_scores = []
                for v in verts:
                    d = sum(A[v][w] for w in verts if w != v)
                    sub_scores.append(d)
                sub_scores = tuple(sorted(sub_scores))
                score_dist[sub_scores] += 1

                # c3 and c5 of sub-tournament
                c3_local = 0
                c3_sets_local = []
                for a, b, c in combinations(verts, 3):
                    if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                        c3_local += 1
                        c3_sets_local.append(frozenset([a, b, c]))

                # c5 via DP
                nn = 5
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
                c5_local = 0
                for v in range(nn):
                    if v == start:
                        continue
                    key = (full, v)
                    if key in dp and dp[key] > 0:
                        if A[verts[v]][verts[start]]:
                            c5_local += dp[key]

                ov1_local = 0
                for i in range(len(c3_sets_local)):
                    for j in range(i+1, len(c3_sets_local)):
                        if len(c3_sets_local[i] & c3_sets_local[j]) == 1:
                            ov1_local += 1

                type_dist[(c3_local, c5_local, ov1_local)] += 1

            # Check if score distribution is orientation-independent
            sigma = tuple(1 if (bits >> i) & 1 else -1 for i in range(m))
            dist_by_orient[sigma] = (dict(score_dist), dict(type_dist))

        # Compare all orientations
        orientations = list(dist_by_orient.keys())
        ref_score = dist_by_orient[orientations[0]][0]
        ref_type = dist_by_orient[orientations[0]][1]

        score_same = all(dist_by_orient[s][0] == ref_score for s in orientations)
        type_same = all(dist_by_orient[s][1] == ref_type for s in orientations)

        print(f"    Score seq distribution same across orientations: {score_same}")
        print(f"    (c3,c5,ov1) type distribution same: {type_same}")

        if not score_same:
            # Find differences
            for s in orientations:
                if dist_by_orient[s][0] != ref_score:
                    print(f"    First different orientation: {s}")
                    for key in set(list(ref_score.keys()) + list(dist_by_orient[s][0].keys())):
                        r = ref_score.get(key, 0)
                        d = dist_by_orient[s][0].get(key, 0)
                        if r != d:
                            print(f"      scores {key}: ref={r}, this={d}")
                    break

        if not type_same:
            print(f"    (c3,c5,ov1) type distributions DIFFER but weighted sum is constant!")
            # Show the variation
            for s in orientations[:3]:
                td = dist_by_orient[s][1]
                weighted_sum = sum(cnt * (c5 - 2*ov1) for (c3, c5, ov1), cnt in td.items())
                print(f"      orient {s[:3]}...: sum(cnt*(c5-2*ov1)) = {weighted_sum}")

    # ====== PART 4: Prove K(p) formula ======
    print("\n" + "=" * 60)
    print("PART 4: K(p) = c5 - 2*ov1 - 2*ov2 formula derivation")
    print("=" * 60)

    # K(p) values:
    # p=7: -126, p=11: -1386, p=13: -3276
    # const(p) = C(c3,2) + K/2
    # p=7: 91 + (-126)/2 = 91 - 63 = 28 OK
    # p=11: 1485 + (-1386)/2 = 1485 - 693 = 792 OK

    print(f"\n  K(p) values and attempts at closed form:")
    K_vals = {}
    for p in [7, 11, 13, 17, 19, 23]:
        m = (p - 1) // 2
        c3 = math.comb(p, 3) - p * math.comb(m, 2)

        # Compute K from Interval tournament
        S = list(range(1, m + 1))
        A = build_adj(p, S)
        c5 = 0
        for subset in combinations(range(p), 5):
            verts = list(subset)
            nn = 5
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
                        c5 += dp[key]

        c3_sets = count_c3_vertex_sets(A, p)
        ov1_global = 0
        ov2_global = 0
        for i in range(len(c3_sets)):
            for j in range(i+1, len(c3_sets)):
                overlap = len(c3_sets[i] & c3_sets[j])
                if overlap == 1:
                    ov1_global += 1
                elif overlap == 2:
                    ov2_global += 1

        K = c5 - 2 * ov1_global - 2 * ov2_global
        K_vals[p] = K
        const_p = math.comb(c3, 2) + K // 2

        print(f"    p={p}: K={K}, c3={c3}, const={const_p}")

    # Factor K values
    print(f"\n  K/(p*(p-1)*(p-2)):")
    for p in sorted(K_vals):
        K = K_vals[p]
        base = p * (p-1) * (p-2)
        print(f"    p={p}: K={K}, p(p-1)(p-2)={base}, ratio={K/base:.6f}")

    print(f"\n  K/(p*(p-1)*(p-3)):")
    for p in sorted(K_vals):
        K = K_vals[p]
        base = p * (p-1) * (p-3)
        print(f"    p={p}: K={K}, base={base}, ratio={K/base:.6f}")

    # -K/(p*(p-1)):
    print(f"\n  -K/(p*(p-1)):")
    for p in sorted(K_vals):
        K = K_vals[p]
        base = p * (p-1)
        print(f"    p={p}: -K={-K}, base={base}, ratio={-K/base:.6f}")

    # -K/C(p,3):
    print(f"\n  -K/C(p,3):")
    for p in sorted(K_vals):
        K = K_vals[p]
        base = math.comb(p, 3)
        print(f"    p={p}: -K={-K}, C(p,3)={base}, ratio={-K/base:.6f}")

    # -K/C(p,4):
    print(f"\n  -K/C(p,4):")
    for p in sorted(K_vals):
        K = K_vals[p]
        base = math.comb(p, 4)
        print(f"    p={p}: -K={-K}, C(p,4)={base}, ratio={-K/base:.6f}")

    # Try: -K = C(p,3) * (something)
    print(f"\n  Factored -K:")
    for p in sorted(K_vals):
        K = K_vals[p]
        neg_K = -K
        # Factor
        factors = []
        n = neg_K
        for f in range(2, 100):
            while n % f == 0:
                factors.append(f)
                n //= f
        if n > 1:
            factors.append(n)
        print(f"    p={p}: -K={neg_K} = {'*'.join(map(str, factors))}")

    # ====== PART 5: Double counting attempt ======
    print("\n" + "=" * 60)
    print("PART 5: Double counting approach")
    print("=" * 60)

    # Identity: c5 = sum_V c5(V) (each 5-cycle in exactly 1 five-subset)
    # ov1 = sum_V ov1(V) (each ov1 pair in exactly 1 five-subset)
    # ov2 = (1/(p-4)) * sum_V ov2(V)
    #
    # So K = c5 - 2*ov1 - 2*ov2
    #      = sum_V c5(V) - 2*sum_V ov1(V) - 2/(p-4) * sum_V ov2(V)
    #      = sum_V [c5(V) - 2*ov1(V) - 2/(p-4)*ov2(V)]
    #
    # On V: ov1(V) + ov2(V) = C(c3(V), 2)
    # So ov1(V) = C(c3(V),2) - ov2(V)
    # => c5(V) - 2*ov1(V) - 2/(p-4)*ov2(V)
    #  = c5(V) - 2*(C(c3(V),2) - ov2(V)) - 2/(p-4)*ov2(V)
    #  = c5(V) - 2*C(c3(V),2) + ov2(V)*(2 - 2/(p-4))
    #  = c5(V) - 2*C(c3(V),2) + ov2(V)*2*(p-5)/(p-4)
    #
    # For p=7: 2*(p-5)/(p-4) = 2*2/3 = 4/3
    # So K = sum_V [c5(V) - 2*C(c3(V),2) + (4/3)*ov2(V)]

    # Let's verify this at p=7
    for p in [7, 11]:
        m = (p - 1) // 2
        S = list(range(1, m + 1))
        A = build_adj(p, S)
        c3_sets = count_c3_vertex_sets(A, p)

        coeff_ov2 = 2 * (p - 5) / (p - 4)
        total_K_local = 0

        for subset in combinations(range(p), 5):
            V = frozenset(subset)
            verts = list(subset)

            # c3 local
            c3_in = [c for c in c3_sets if c.issubset(V)]
            c3_local = len(c3_in)

            # c5 local
            nn = 5
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
            c5_local = 0
            for v in range(nn):
                if v == start:
                    continue
                key = (full, v)
                if key in dp and dp[key] > 0:
                    if A[verts[v]][verts[start]]:
                        c5_local += dp[key]

            # ov2 local
            ov2_local = 0
            for i in range(len(c3_in)):
                for j in range(i+1, len(c3_in)):
                    if len(c3_in[i] & c3_in[j]) == 2:
                        ov2_local += 1

            C_c3_2 = c3_local * (c3_local - 1) // 2
            local_contrib = c5_local - 2 * C_c3_2 + coeff_ov2 * ov2_local
            total_K_local += local_contrib

        print(f"\n  p={p}: sum of local contributions = {total_K_local:.4f}")
        print(f"    K(p) = {K_vals[p]}")
        print(f"    Match: {abs(total_K_local - K_vals[p]) < 0.01}")

    # ====== PART 6: Score sequence determines 5-subset type distribution ======
    print("\n" + "=" * 60)
    print("PART 6: Score-sequence constancy for 5-subsets of regular")
    print("=" * 60)

    # For a REGULAR tournament on p vertices (all scores = m = (p-1)/2):
    # What is the distribution of score sequences of 5-subsets?
    # Claim: it depends only on p (not on the specific regular tournament)
    # because the number of transitive triples through any vertex is determined
    # by the degree sequence alone.

    # This is related to the "Ramsey" structure of regular tournaments.
    # In a regular tournament, C(p, k) with score profile is determined by "switching classes"

    # Let's verify at p=7: all orientations give same 5-subset score distribution?
    p = 7
    m = 3
    print(f"\n  p={p}: checking 5-subset score distribution")

    all_score_dists = []
    for bits, S in all_orientations_circ(p):
        A = build_adj(p, S)
        score_dist = defaultdict(int)
        for subset in combinations(range(p), 5):
            verts = list(subset)
            sub_scores = tuple(sorted([sum(A[v][w] for w in verts if w != v) for v in verts]))
            score_dist[sub_scores] += 1
        all_score_dists.append(dict(score_dist))

    is_same = all(d == all_score_dists[0] for d in all_score_dists)
    print(f"    All orientations have same 5-subset score distribution: {is_same}")
    if is_same:
        print(f"    Distribution: {all_score_dists[0]}")
    else:
        # Show how many distinct distributions
        unique_dists = set(tuple(sorted(d.items())) for d in all_score_dists)
        print(f"    {len(unique_dists)} distinct score distributions among {len(all_score_dists)} orientations")
        for ud in sorted(unique_dists):
            count = sum(1 for d in all_score_dists if tuple(sorted(d.items())) == ud)
            print(f"      {count} orientations: {dict(ud)}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
