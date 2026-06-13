#!/usr/bin/env python3
"""
Why does the hereditary maximizer chain break at even n?

At n=7: all vertex deletions give H=45 (max at n=6). HEREDITARY.
At n=6: all vertex deletions give H=11 or H=13. NOT max at n=5 (max=15). NOT HEREDITARY.
At n=5: 24/64 maximizers have all deletions giving max at n=4.
At n=4: 0/24 maximizers have all deletions giving max at n=3.

Question: What are the score sequences of maximizers at each n?
And: What score sequences appear in vertex-deletions of maximizers?

Hypothesis: At even n, the maximizer has self-complementary score (s_i + s_{n-1-i} = n-1).
Deleting a vertex from a SC-score tournament at even n gives a non-SC score at odd n-1,
which cannot be the score class of the (n-1)-maximizer (which has regular score at odd n).

kind-pasteur-2026-03-06-S18f
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count

MAX_H = {2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189}

def score_seq(T):
    return tuple(sorted(sum(T[i]) for i in range(len(T))))

# ============================================================
# Phase 1: Score sequences of all maximizers
# ============================================================
print("=" * 70)
print("SCORE SEQUENCES OF MAXIMIZERS")
print("=" * 70)

maximizer_scores = {}

for n in range(3, 8):
    m = n * (n - 1) // 2
    best_h = 0
    score_groups = {}

    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        if h > best_h:
            best_h = h
            score_groups = {}
        if h == best_h:
            s = score_seq(T)
            score_groups[s] = score_groups.get(s, 0) + 1

    maximizer_scores[n] = score_groups
    print(f"\nn={n}: max H = {best_h}")
    for s, cnt in sorted(score_groups.items()):
        is_regular = all(x == (n-1)/2 for x in s) if n % 2 == 1 else False
        is_sc = all(s[i] + s[n-1-i] == n-1 for i in range(n))
        print(f"  score {s}: {cnt} tournaments (regular={is_regular}, self-comp={is_sc})")

# ============================================================
# Phase 2: Score sequences of vertex-deletions of maximizers
# ============================================================
print(f"\n{'='*70}")
print("DELETION SCORE SEQUENCES")
print("=" * 70)

for n in range(4, 8):
    m = n * (n - 1) // 2
    max_h = MAX_H[n]

    print(f"\nn={n} (max H={max_h}):")

    # Find first maximizer
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        if hamiltonian_path_count(T) == max_h:
            break

    parent_score = score_seq(T)
    print(f"  Parent score: {parent_score}")

    del_scores = {}
    for v in range(n):
        verts = [i for i in range(n) if i != v]
        sub = [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]
        h_sub = hamiltonian_path_count(sub)
        s = score_seq(sub)
        del_scores[s] = del_scores.get(s, (0, h_sub))
        del_scores[s] = (del_scores[s][0] + 1, h_sub)

    # What's the maximizer score at n-1?
    max_h_n1 = MAX_H[n-1]
    max_scores_n1 = set(maximizer_scores[n-1].keys())

    print(f"  Maximizer scores at n-1={n-1}: {max_scores_n1}")
    print(f"  Deletion scores:")
    for s, (cnt, h) in sorted(del_scores.items()):
        in_max = s in max_scores_n1
        is_sc = all(s[i] + s[n-2-i] == n-2 for i in range(n-1))
        is_regular = all(x == (n-2)/2 for x in s) if (n-1) % 2 == 1 else False
        print(f"    {s}: {cnt} vertices, H(T-v)={h}, in max scores={in_max}, "
              f"SC={is_sc}, regular={is_regular}")

# ============================================================
# Phase 3: The score obstruction explained
# At even n, the maximizer has self-comp score: s_i + s_{n-1-i} = n-1
# For n=6: score (2,2,2,3,3,3), deleting vertex with degree 2 or 3
# gives score at n=5 that's NOT regular ((2,2,2,2) from complement structure)
#
# At odd n, the maximizer has REGULAR score: all s_i = (n-1)/2
# Deleting vertex with degree (n-1)/2 gives score at n-1 that's...
# what? s_j -> s_j - T[j][v] for j != v
# If v has score (n-1)/2, half its arcs go out. The deletion changes
# scores non-uniformly.
# ============================================================
print(f"\n{'='*70}")
print("SCORE CHANGE UNDER VERTEX DELETION")
print("=" * 70)

# At n=7 (Paley T_7), all vertices have degree 3
# Deleting v: each j != v loses the arc j-v
# j->v: score of j decreases by 1
# v->j: score of j stays the same
# So 3 vertices lose 1 from their score (the successors of v)
# and 3 vertices keep their score (the predecessors of v)
# If original score is all 3's, new score is:
# 3 vertices with score 3 (predecessors of v: they lost their arc FROM v)
# Wait: v->j means j is a successor of v. When we delete v,
# j loses the incoming arc from v, so j's OUT-degree is unchanged.
# v->j: v has outgoing arc to j. For j's score: j's score = #{k: j->k}.
# v is not a target of j's outgoing arcs. Deleting v removes one
# potential target. If j->v, then j's score decreases by 1.
# If v->j, then j's score stays the same.

# So: delete v with out-degree d_v = (n-1)/2.
# Predecessors of v (vertices j with j->v): d_v = (n-1)/2 vertices, lose 1 score
# Successors of v (vertices j with v->j): n-1-d_v = (n-1)/2 vertices, keep score
# New score sequence (n-1 vertices):
# (n-1)/2 vertices with score (n-1)/2 - 1 = (n-3)/2
# (n-1)/2 vertices with score (n-1)/2

for n in [5, 7]:
    d = (n - 1) // 2
    print(f"\nn={n} regular tournament, d = {d}:")
    print(f"  Original score: all {d}")
    print(f"  After deleting v with out-degree {d}:")
    print(f"  - {d} predecessors of v: score {d} -> {d-1}")
    print(f"  - {d} successors of v: score {d} (unchanged)")
    del_score = tuple(sorted([d-1]*d + [d]*d))
    is_sc = all(del_score[i] + del_score[n-2-i] == n-2 for i in range(n-1))
    print(f"  Deletion score: {del_score}")
    print(f"  Self-complementary: {is_sc}")
    print(f"  Is this a maximizer score at n-1={n-1}? ", end="")
    print(f"{del_score in maximizer_scores.get(n-1, {})}")

for n in [4, 6]:
    print(f"\nn={n} SC tournament, score {tuple(sorted([(n-1)//2]*n + [n//2]*0)[:n])}:")
    # Actually get the real score
    for bits in range(1 << (n*(n-1)//2)):
        T = tournament_from_bits(n, bits)
        if hamiltonian_path_count(T) == MAX_H[n]:
            s = score_seq(T)
            print(f"  Score: {s}")
            for v in range(n):
                verts = [i for i in range(n) if i != v]
                sub = [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]
                ds = score_seq(sub)
                h_sub = hamiltonian_path_count(sub)
                print(f"  v={v}: deg(v)={sum(T[v])}, del score={ds}, H(T-v)={h_sub}")
            break

# ============================================================
# Phase 4: The n -> n-2 score analysis
# For even n -> even n-2: delete 2 vertices
# Score compatibility?
# ============================================================
print(f"\n{'='*70}")
print("n -> n-2 DELETION SCORE ANALYSIS")
print("=" * 70)

# n=7 -> n=5: pair deletions give H < 15 = max at n=5
# Check scores
for bits in range(1 << 21):
    T = tournament_from_bits(7, bits)
    if hamiltonian_path_count(T) == 189:
        print(f"\nn=7 maximizer (bits={bits}):")
        for v1 in range(7):
            for v2 in range(v1+1, 7):
                verts = [i for i in range(7) if i != v1 and i != v2]
                sub = [[T[verts[i]][verts[j]] for j in range(5)] for i in range(5)]
                ds = score_seq(sub)
                h_sub = hamiltonian_path_count(sub)
                if (v1, v2) in [(0,1), (0,2), (0,3)]:
                    print(f"  del ({v1},{v2}): score={ds}, H={h_sub}")
        # Count score types
        scores = {}
        for v1 in range(7):
            for v2 in range(v1+1, 7):
                verts = [i for i in range(7) if i != v1 and i != v2]
                sub = [[T[verts[i]][verts[j]] for j in range(5)] for i in range(5)]
                ds = score_seq(sub)
                h_sub = hamiltonian_path_count(sub)
                scores[(ds, h_sub)] = scores.get((ds, h_sub), 0) + 1
        print(f"  Pair deletion types:")
        for (ds, h), cnt in sorted(scores.items(), key=lambda x: -x[1]):
            is_reg = all(x == 2 for x in ds)
            print(f"    score={ds} H={h}: {cnt} pairs (regular={is_reg})")
        break

print("\nDone.")
