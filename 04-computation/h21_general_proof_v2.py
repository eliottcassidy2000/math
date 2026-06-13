#!/usr/bin/env python3
"""
GENERAL PROOF ATTEMPT for H=21 impossibility at ALL n.

THE PROOF STRUCTURE:
1. Base case: n <= 8, proved exhaustively (Part G).
2. Induction: For n >= 9, if vertex v not in any 3-cycle,
   remove v (Part J), reduce to n-1. By induction, H(T-v) != 21, so H(T) != 21.
3. Remaining case: every vertex in a 3-cycle, n >= 9.

For case 3, we need: alpha_1 > 10 (since H=21 needs alpha_1+2*alpha_2+...=10).

APPROACH A: Forced alpha_3 >= 1 (when 3 disjoint 3-cycles exist)
   Part C shows alpha_3 >= 1 => alpha_1+2*alpha_2+4*alpha_3 >= 13 > 10. Blocked.

APPROACH B: Even without alpha_3, if alpha_1 > 10 directly.

APPROACH A is cleaner. When does alpha_3 >= 1 hold?

CLAIM: If every vertex in a 3-cycle at n >= 9, then at least 3 pairwise-disjoint
3-cycles exist (as DIRECTED cycles in Omega).

PROOF OF CLAIM:
The 3-cycles cover all n vertices. Consider the 3-cycle VERTEX SETS (each a triple).
These triples cover [n]. We want 3 pairwise-disjoint triples among them.

This is related to the SUNFLOWER LEMMA / MATCHING THEORY for hypergraphs.

By the Erdos-Ko-Rado theorem: a pairwise-intersecting family of 3-element subsets
of [n] has at most C(n-1, 2) members (for n >= 7, the maximum family consists of
all triples through a fixed element).

The TURÁN-type question: the maximum number of 3-element subsets of [n] with
no 3 pairwise-disjoint members is characterized by the Hilton-Milner-Frankl theory.

But we don't have a MAXIMUM question — we have a COVERAGE question:
the t3 triples COVER all n vertices, and we want 3 pairwise-disjoint ones.

KEY LEMMA: If F is a family of 3-element subsets of [n] covering all of [n]
(union(F) = [n]), and n >= 9, then F contains 3 pairwise-disjoint members.

PROOF: Suppose for contradiction that F has no 3 pairwise-disjoint members.
Then the maximum matching in the 3-uniform hypergraph F has size at most 2.

By a theorem of Ray-Chaudhuri and Wilson (or direct argument):
If no matching of size 3 in F, then all triples in F are contained in the union
of 2 fixed triples, covering at most 6 vertices.

Wait, that's not quite right. Let me think more carefully.

If the max matching = 2: there exist A, B in F with A cap B = empty, and
every C in F intersects A or B. So every C in F has C cap (A union B) != empty.
Since |A union B| = 6, every triple in F uses at least one of 6 specific vertices.
Coverage: union(F) <= A union B union (extra from other triples) = at most [n].
But every triple in F hits {a1,...,a6} = A union B. The remaining vertices [n]\(A union B)
have at least n-6 >= 3 elements. Each must be covered by some triple, and that triple
has 2 elements in {a1,...,a6} and 1 in [n]\{a1,...,a6}. But that only covers 1 new vertex
per triple. Need n-6 such triples.

Wait, a triple C in F with C cap (A union B) != empty could have:
- 1 vertex in A union B, 2 outside. Then C is disjoint from some other triple that has
  1 vertex in A union B and 2 different outside vertices, potentially giving a 3rd
  disjoint triple.

Hmm, let me reconsider. We assumed max matching = 2. So C is disjoint from NEITHER A nor B.
But C could be disjoint from other triples. The condition is: no 3 PAIRWISE disjoint, meaning
no independent set of size 3 in the disjointness graph. This is weaker than max matching = 2.

Actually: 3 pairwise-disjoint triples = matching of size 3 = independent set of size 3
in the intersection graph. These are the same thing.

So: max matching M of disjoint triples in F has |M| <= 2.

By Konig-type reasoning for hypergraphs:
If |M| <= 2, then there exists a set S of at most 2*(3-1) = 4 vertices that
intersects every triple in F (a transversal).

Wait, that's the Konig/Erdos-Gallai theorem for matchings in hypergraphs.
For 3-uniform hypergraphs: if max matching <= 2, then there exists a cover of
size at most max(2, C) that intersects all triples. But this isn't tight.

Actually, the result I need is:

If F is a 3-uniform hypergraph with no matching of size 3, then F can be covered
by at most 6 vertices (the union of some 2 disjoint triples in F, or fewer).

Hmm, this is related to the Erdos matching conjecture (now proved by Frankl).

For 3-uniform, n >= 9:
Frankl's theorem: If F is 3-uniform on [n] with n >= 9 and no matching of size 3,
then |F| <= max(C(5, 3), ?) = ... well, the exact bound depends.

But I don't need the size bound. I need the COVERAGE bound:
If F covers [n] with n >= 9 and has no matching of size 3,
then... actually, this is impossible!

PROOF: If max matching = 2, take M = {A, B} with A cap B = empty, |A union B| = 6.
Every C in F has C cap (A union B) != empty (else C, A, B is a matching of size 3).
So every triple in F uses at least one vertex from A union B.
The vertices in [n] \ (A union B) have n - 6 >= 3 elements.
Each must be in some triple C in F. That triple has C cap (A union B) != empty
and C cap ([n]\(A union B)) != empty.
So C has at most 2 vertices in A union B and at least 1 outside.

But this only shows coverage is consistent. It doesn't give a contradiction.
The issue: n - 6 vertices outside can each be covered by a triple that also
intersects A union B.

If max matching = 1 or 0: ALL triples share a common vertex v.
Coverage requires: every other vertex is in some triple {v, x, y}.
So at most 1 + 2*t triples can cover 1 + 2t vertices. For n >= 9: need t >= 4 triples.

If max matching = 2 but NOT 3: we have {A, B} disjoint, covering 6 vertices.
The remaining n - 6 >= 3 vertices are each in some triple intersecting A union B.
Can we find 3 pairwise-disjoint triples?
We already have A and B (disjoint). Need C disjoint from both A and B.
C must have C cap (A union B) = empty (to be disjoint from both).
But then {A, B, C} is a matching of size 3. Contradiction with max matching = 2.

So: C must intersect A or B.
But we need C to contain some vertex outside A union B (to cover remaining vertices).
C has 3 vertices: at least 1 in A union B, and some outside.
But C is NOT disjoint from A or B. So C intersects at least one of {A, B}.
This doesn't mean C intersects BOTH.

If C intersects A but not B: then {B, C} might be disjoint? No, C intersects B too
(otherwise {A, B, C} would be a matching iff they're pairwise disjoint, but C intersects A).
Wait: matching of size 3 means THREE pairwise-disjoint triples.
{A, B, C}: A cap B = empty (by construction). Need A cap C = empty AND B cap C = empty.
If A cap C != empty, then {A, B, C} is NOT a matching of size 3 even if B cap C = empty.
So: it's possible for C to intersect A but not B, without forming a 3-matching.

But then: consider C and B. They're disjoint (C cap B = empty by assumption).
Now look for D disjoint from both C and A.
D must not intersect C (3 vertices) and not intersect A (3 vertices).
C and A share at least 1 vertex (C cap A != empty), so |C union A| <= 5.
D avoids C union A, so D subset of [n] \ (C union A), which has n - |C union A| >= 9 - 5 = 4 vertices.
Need D to be a triple in F from those 4+ vertices. But D must also be in F (a 3-cycle in the tournament).
This doesn't GUARANTEE D exists.

OK, the pure combinatorial argument about 3-uniform hypergraph matching doesn't directly work. The issue is that F (the family of 3-cycle vertex sets) might have max matching 2 while covering all 9 vertices, as long as every uncovered vertex is in some triple that intersects the matching.

BUT: the tournament structure adds extra constraints. The 3-cycles aren't arbitrary triples — they're constrained by the tournament's arc structure.

Let me try a DIFFERENT approach: count-based.

At n=9, no source/sink, every vertex in 3-cycle:
- Moon: t3 = 84 - sum C(s_v, 2).
- Scores in {1,...,7}.
- For t3 <= 10: sum C(s_v, 2) >= 74.

The minimum t3 we found was 3 (from score sequences where sum C(s_v,2) = 81, but actually
my score sequence generator showed min t3 = 2 at n=9 for no-source/sink. Can we check which
score sequences have t3 <= 10 AND allow every vertex in a 3-cycle?

For "every vertex in a 3-cycle": vertex v with score s has b_v >= 1 backward arcs.
The maximum of sum_v b_v = 3*t3. The minimum of b_v for vertex with score s:
b_v >= 1 always possible unless s = 0 or s = n-1. With s in {1,...,n-2}: b_v >= 1 is
always ACHIEVABLE (not forced). So the score sequence alone doesn't guarantee every vertex
is in a 3-cycle.

The computational check is needed.

Instance: kind-pasteur-2026-03-07-S33
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from math import comb

def main():
    # For n >= 9, the argument goes:
    # 1. Part J removes vertices not in 3-cycles
    # 2. After removal, if n' <= 8, base case applies
    # 3. If n' >= 9, every vertex in 3-cycle, t3 >= ceil(n'/3) >= 3
    #
    # CASE 3a: t3 = 3 at n' = 9.
    #   3 three-cycles covering 9 vertices => all disjoint => alpha_3 >= 1 => Part C
    #
    # CASE 3b: t3 >= 4 at n' >= 9.
    #   Need to show alpha_1 > 10.
    #   From the all-through-v analysis: if all 3-cycles share a vertex,
    #   then T-v is transitive, and 5-cycles through v grow as C(n-2, 3).
    #   At n=9: C(7,3) = 35 >> 10.
    #
    #   If NOT all share a vertex: there exist 2 disjoint 3-cycles {A, B}.
    #   These contribute 1 to alpha_2. If a 3rd disjoint cycle exists: alpha_3 >= 1 => Part C.
    #   If no 3rd disjoint: every other cycle intersects A or B.
    #   The 5-cycles on the remaining vertices still count.
    #
    # The issue: we need to handle the case "max matching = 2" carefully.
    # With 2 disjoint 3-cycles and t3 >= 4: at least 2 more 3-cycles.
    # Each intersects A or B.
    # In total: alpha_1 >= t3 + t5 + ... Many 5-cycles from inter-triple arcs.

    # Let me compute: at n=9, what score sequences have t3 <= 10, no src/sink,
    # AND could potentially have every vertex in a 3-cycle?

    # A necessary (not sufficient) condition: b_v >= 1 for each v, meaning
    # 3*t3 >= 9, t3 >= 3.

    # Score sequences with t3 in {3,...,10} and scores in {1,...,7}:
    n = 9
    cn3 = comb(n, 3)  # 84

    def gen_no_ss_scores():
        target = n * (n - 1) // 2  # 36

        def backtrack(idx, rem, prev, seq):
            if idx == n:
                if rem == 0:
                    yield tuple(seq)
                return
            lo = max(prev, 1)
            hi = n - 2  # 7
            for s in range(lo, hi + 1):
                if rem - s < 0:
                    break
                remaining = n - idx - 1
                if remaining > 0:
                    max_remain = remaining * hi
                    min_remain = remaining * max(s, 1)
                    if rem - s > max_remain or rem - s < min_remain:
                        continue
                new_seq = seq + [s]
                k = len(new_seq)
                if sum(new_seq) < comb(k, 2):
                    continue
                yield from backtrack(idx + 1, rem - s, s, new_seq)

        yield from backtrack(0, target, 1, [])

    print(f"=== n={n}: Score sequences with t3 in [3, 10] and no source/sink ===")

    viable = []
    for scores in gen_no_ss_scores():
        sum_cs2 = sum(comb(s, 2) for s in scores)
        t3 = cn3 - sum_cs2
        if 3 <= t3 <= 10:
            viable.append((scores, t3))

    print(f"Found {len(viable)} viable score sequences")
    for scores, t3 in sorted(viable, key=lambda x: x[1]):
        print(f"  {scores}: t3={t3}")

    # KEY OBSERVATION: At t3=3, the 3 three-cycles must partition [9].
    # This forces alpha_3 >= 1 => Part C blocks H=21.
    #
    # At t3=4,...,10: do tournaments with these scores AND every vertex in 3-cycle exist?
    # And if so, what is the minimum alpha_1?
    #
    # For the proof to work, we need EITHER:
    # (a) alpha_3 >= 1 (3 disjoint 3-cycles) => Part C, OR
    # (b) alpha_1 > 10 directly.
    #
    # At t3=3 at n=9: ALWAYS (a) — PROVED.
    # At t3 >= 4 at n=9: might not have (a), need (b).
    #
    # But with t3 >= 4 at n=9, the 5-cycle count is substantial.
    # Need to check computationally.

    print(f"\n=== At t3=3: alpha_3 >= 1 by partition argument ===")
    print(f"t3=3 scores: {[s for s, t in viable if t == 3]}")
    print("3*3 = 9 vertex-slots, exactly n=9. Each vertex in exactly 1 three-cycle.")
    print("=> 3 pairwise-disjoint three-cycles => alpha_3 >= 1")
    print("=> Part C: alpha_1+2*alpha_2+4*alpha_3 >= 3+6+4 = 13 > 10. H!=21.")
    print()

    print(f"=== Score sequences at t3=3 that might allow every vertex in 3-cycle ===")
    for scores, t3 in viable:
        if t3 == 3:
            # With t3=3: exactly 3 three-cycles. Each vertex in exactly 1 (since 3*3=9=n).
            # The backward arcs: sum b_v = 3*3 = 9. With 9 vertices, each b_v = 1.
            # This means: each vertex has exactly 1 backward arc.
            # b_v = arcs from N+(v) to N-(v).
            # For vertex with score s: b_v = 1. Total mixed arcs = s*(n-1-s).
            # So only 1 out of s*(n-1-s) cross arcs goes backward.
            print(f"  {scores}: b_v = 1 for each v")
            print(f"    Cross arcs: {[s*(8-s) for s in scores]}")
            print(f"    This structure: nearly transitive, 1 backward arc per vertex")


if __name__ == "__main__":
    main()
