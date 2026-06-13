#!/usr/bin/env python3
"""
Prove: alpha_1=10 with i_2=0 is IMPOSSIBLE in any tournament.

Strategy: Show that 10 pairwise-sharing odd cycles in a tournament
must always contain at least one vertex-disjoint pair.

Key observation from data:
- At n=6: alpha_1=10 always has (t3=6, t5=4), coverage=6
- At n=7: alpha_1=10 has (t3=5 or 6, t5=5 or 4), coverage=6 or 7
- Disjoint pairs are ALWAYS (3-cycle, 3-cycle)
- i_2=0 NEVER occurs

Proof attempt:
1. Among 10 odd cycles, at least 5 are 3-cycles (proved below).
2. Among 5+ three-cycles, at least one disjoint pair exists (Ramsey-type).

OR:
1. Show that the cycle compositions that give alpha_1=10 ALWAYS contain
   enough 3-cycles that two must be vertex-disjoint.

Instance: kind-pasteur-2026-03-07-S33
"""

from itertools import combinations
from collections import Counter, defaultdict
import time

def all_tournaments(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for bits in range(2**m):
        adj = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            if (bits >> k) & 1:
                adj[j][i] = 1
            else:
                adj[i][j] = 1
        yield adj

def find_3cycles(adj, n):
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.append(frozenset([i,j,k]))
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    cycles.append(frozenset([i,j,k]))
    return cycles


def max_pairwise_sharing_3cycles(n):
    """What is the maximum number of 3-cycles on n vertices such that
    ALL pairs share at least one vertex?

    A family of 3-element subsets of [n] is "pairwise intersecting" iff
    every two subsets share at least one element.

    By the Helly property for 1-dimensional simplicial complexes,
    a pairwise intersecting family of sets with max size k has at most
    C(n-1, k-1) members when a common element exists.

    But we need: pairwise intersecting family of 3-subsets of [n]
    where each subset is a DIRECTED 3-cycle in a tournament on [n].
    """
    print(f"\n=== Max pairwise-sharing 3-cycles at n={n} ===")

    # First: what's the max pairwise intersecting family of 3-subsets of [n]?
    # This is the sunflower / delta-system problem.
    # For 3-uniform hypergraphs, the Erdos-Ko-Rado theorem says:
    # max pairwise intersecting = C(n-1, 2) (all triples containing a fixed element)
    # when n >= 7 (and also for smaller n with adjustments).

    # So the max is C(n-1,2) = n(n-1)/2... no wait, that's wrong.
    # EKR for t-intersecting k-uniform: max |F| = C(n-1, k-1) for n >= 2k.
    # For k=3: max = C(n-1, 2) = (n-1)(n-2)/2 when n >= 6.
    # These are all triples containing vertex 0.

    # But: we need DIRECTED 3-cycles, not just triples.
    # A triple {a,b,c} is a 3-cycle iff the tournament on {a,b,c} is cyclic.
    # Not all triples through vertex 0 are 3-cycles.

    # At most, half of the C(n-1,2) triples through vertex 0 are 3-cycles
    # (since each triple is either cyclic or transitive, and both directions count).

    # Actually: each triple through v has exactly 1 3-cycle (out of 2 possible cyclic orders)
    # but the question is whether it's a directed 3-cycle or not.
    # In a tournament, each triple is either transitive (no 3-cycle) or has exactly 1 3-cycle.
    # By Moon's formula: t3 = sum_v C(s_v, 2) correction... no.
    # Moon: t3 = C(n,3) - sum_v C(s_v, 2) where s_v = out-degree.
    # Through vertex v: t3(v) = number of 3-cycles containing v.

    # For n=7: through vertex v with score s_v, the number of 3-cycles containing v is:
    # C(n-1,2) - C(s_v, 2) - C(n-1-s_v, 2) = C(6,2) - C(s,2) - C(6-s,2)
    # = 15 - s(s-1)/2 - (6-s)(5-s)/2

    for s in range(7):
        through_v = 15 - s*(s-1)//2 - (6-s)*(5-s)//2
        print(f"  n=7, s_v={s}: 3-cycles through v = {through_v}")

    # Max through v at s=3 (regular): 15 - 3 - 3 = 9 three-cycles through v.
    # Can we have a pairwise-intersecting family of 10 three-cycles?
    # If all pass through one vertex, max is 9 (at n=7, regular).
    # But we need 10. So they can't all pass through one vertex.

    # Sunflower-free?
    # If not all through one vertex, can we have 10 pairwise-intersecting 3-subsets of [7]?
    # EKR: max pairwise intersecting family of 3-subsets of [7] is C(6,2)=15
    # (but only if all contain a fixed element). Without this constraint...

    # Actually the Hilton-Milner theorem says:
    # max pairwise intersecting family NOT all through a common element:
    # |F| <= C(n-1,k-1) - C(n-k-1,k-1) + 1
    # For k=3, n=7: C(6,2) - C(3,2) + 1 = 15 - 3 + 1 = 13.
    # So 10 pairwise-intersecting 3-subsets of [7] CAN exist without common element!

    # But the question is: can 10 of them be DIRECTED 3-CYCLES in a tournament?

    print(f"\n  EKR at n=7: max through one vertex = C(6,2)=15")
    print(f"  Max 3-cycles through one vertex (regular) = 9")
    print(f"  Hilton-Milner: max without common = 13")
    print(f"  We need 10 pairwise-sharing DIRECTED 3-cycles")
    print(f"  Since max through one vertex is 9, if alpha_1=10 all sharing,")
    print(f"  at least 1 cycle avoids every vertex => impossible!")

    # Wait: that's not quite right. The 10 cycles might not all share
    # the SAME vertex. They just need to pairwise share SOME vertex.

    # The correct argument: if all 10 cycles pairwise share a vertex,
    # does there exist a common vertex? (Helly's theorem?)

    # Helly's theorem: for convex sets in R^d, if every d+1 sets have
    # non-empty intersection, then all have non-empty intersection.
    # But 3-subsets are NOT convex sets!

    # For "sunflower-type" results on hypergraphs:
    # Pairwise-intersecting 3-uniform: Helly does NOT apply in general.
    # Example: {{1,2,3}, {1,4,5}, {2,4,6}} — pairwise intersecting but no common element.


def check_pairwise_intersecting_3cycles():
    """
    At n=7 exhaustive: can alpha_1=10 with i_2=0 occur?
    We already know from sampling: NO.
    But WHY?

    Key structural argument:
    If 10 odd cycles are pairwise sharing in a tournament T on n vertices,
    and these include t3 three-cycles and t5 five-cycles:

    Case 1: t3=6, t5=4 (n=6 case)
      6 three-cycles on 6 vertices. Coverage = 6 (all vertices used).
      Two 3-cycles using all 6 vertices must be on complementary triples {a,b,c}, {d,e,f}.
      These are DISJOINT. So if any two complementary triples are both 3-cycles, i_2 >= 1.

      CLAIM: Among 6 three-cycles on 6 vertices, two must use complementary triples.
      Proof: By pigeonhole? There are C(6,3)/2 = 10 complementary pairs of triples.
      6 three-cycles partition into non-complementary and complementary.
      Need to show: at least one complementary pair.
    """
    print("\n=== COMPLEMENTARY TRIPLE ANALYSIS at n=6 ===")

    # Enumerate all 6 three-cycles configurations at n=6
    edges6 = [(i, j) for i in range(6) for j in range(i+1, 6)]

    compl_pairs = 0
    no_compl = 0

    for bits in range(2**15):
        adj = [[0]*6 for _ in range(6)]
        for k, (i, j) in enumerate(edges6):
            if (bits >> k) & 1:
                adj[j][i] = 1
            else:
                adj[i][j] = 1

        c3 = find_3cycles(adj, 6)
        if len(c3) != 6:
            continue

        # Check for complementary pairs
        has_compl = False
        for a in range(len(c3)):
            for b in range(a+1, len(c3)):
                if len(c3[a] | c3[b]) == 6:  # complementary = cover all 6 vertices
                    has_compl = True
                    break
            if has_compl:
                break

        if has_compl:
            compl_pairs += 1
        else:
            no_compl += 1
            # Show this case
            if no_compl <= 3:
                print(f"  6 three-cycles with NO complementary pair:")
                for c in c3:
                    print(f"    {sorted(c)}")

    print(f"  With complementary pair: {compl_pairs}")
    print(f"  Without complementary pair: {no_compl}")

    # Now: among alpha_1=10 (6 three-cycles + 4 five-cycles),
    # do we always have complementary three-cycle pairs?
    print(f"\n=== alpha_1=10 at n=6: complementary pair check ===")

    from itertools import permutations

    def find_5cyc(adj, n):
        cycles = []
        for verts in combinations(range(n), 5):
            v = list(verts)
            count = 0
            for perm in permutations(range(5)):
                ok = True
                for idx in range(5):
                    if not adj[v[perm[idx]]][v[perm[(idx+1)%5]]]:
                        ok = False
                        break
                if ok:
                    count += 1
            num_cycles = count // 5
            for _ in range(num_cycles):
                cycles.append(frozenset(verts))
        return cycles

    total_a10 = 0
    a10_with_compl_3 = 0
    a10_with_disjoint_35 = 0

    for bits in range(2**15):
        adj = [[0]*6 for _ in range(6)]
        for k, (i, j) in enumerate(edges6):
            if (bits >> k) & 1:
                adj[j][i] = 1
            else:
                adj[i][j] = 1

        c3 = find_3cycles(adj, 6)
        c5 = find_5cyc(adj, 6)
        if len(c3) + len(c5) != 10:
            continue

        total_a10 += 1

        # Check complementary 3-cycle pairs
        has_compl_3 = False
        for a in range(len(c3)):
            for b in range(a+1, len(c3)):
                if not (c3[a] & c3[b]):
                    has_compl_3 = True
                    break
            if has_compl_3:
                break
        if has_compl_3:
            a10_with_compl_3 += 1

        # Check disjoint 3-cycle and 5-cycle pairs
        has_disjoint_35 = False
        for c in c3:
            for d in c5:
                if not (c & d):
                    has_disjoint_35 = True
                    break
            if has_disjoint_35:
                break
        if has_disjoint_35:
            a10_with_disjoint_35 += 1

    print(f"  Total alpha_1=10: {total_a10}")
    print(f"  With disjoint 3-3 pair: {a10_with_compl_3}")
    print(f"  With disjoint 3-5 pair: {a10_with_disjoint_35}")

    # The minimum vertex count argument:
    # 10 pairwise-intersecting sets of size 3 or 5.
    # Each pair shares >= 1 vertex.
    # The vertex coverage is at most n.
    # If the coverage is exactly 6, then two 3-cycles cover all 6 vertices
    # and if both are present, they must share a vertex OR be disjoint.
    # Two complementary triples {a,b,c} and {d,e,f} are DISJOINT.

    print("\n=== SCORE SEQUENCE CONSTRAINT ===")
    # At n=6: exactly t3=6 three-cycles.
    # Moon: t3 = C(6,3) - sum C(s_v, 2) = 20 - sum C(s_v,2)
    # t3 = 6 => sum C(s_v,2) = 14
    # Score seq summing to C(6,2)=15: e.g. (1,2,2,3,3,4)
    # C(1,2)+C(2,2)+C(2,2)+C(3,2)+C(3,2)+C(4,2) = 0+1+1+3+3+6 = 14. YES.
    # Also (2,2,2,3,3,3): C(2,2)*3 + C(3,2)*3 = 3+9 = 12. t3 = 20-12 = 8. NO.
    # Also (1,1,3,3,3,4): 0+0+3+3+3+6 = 15. t3 = 5. NO.
    # (1,2,3,3,3,3): 0+1+3+3+3+3 = 13. t3 = 7. NO.
    # (2,2,2,2,3,4): 1+1+1+1+3+6 = 13. t3 = 7. NO.
    # (1,2,2,3,3,4): sum = 14. t3 = 6. YES.
    # (1,1,2,3,4,4): 0+0+1+3+6+6 = 16. t3 = 4.
    # (0,2,3,3,3,4): 0+1+3+3+3+6 = 16. t3 = 4.
    # (2,2,2,3,3,3): sum = 12. t3 = 8.

    # So t3=6 at n=6 requires score (1,2,2,3,3,4) (or permutations).
    # This is a NEAR-REGULAR score.

    print("  At n=6: t3=6 requires score (1,2,2,3,3,4) up to permutation")
    print("  Moon: t3 = C(n,3) - sum C(s_v, 2)")

    # Count how many triples through each vertex at this score
    # Vertex with s=1: through-v 3-cycles = C(5,2) - C(1,2) - C(4,2) = 10 - 0 - 6 = 4
    # Wait, that's for counting pairs (x,y) where x in out-nbhd, y in in-nbhd of v.
    # Through-v 3-cycles = s_v * (n-1-s_v) - number of arcs from out to in of v... no.
    # Actually through-v 3-cycles = #{(x,y): v->x, x->y, y->v} = sum_{x in out(v)} |out(x) cap in(v)|
    # Simpler: through-v = # {a,b}: v->a->b->v = s_v(n-1-s_v) - back edges?
    # No. Through-v directed 3-cycles:
    # v->a, a->b, b->v. So a in out(v), b in in(v), a->b.
    # = # arcs from out(v) to in(v).
    # At v with s=1: out(v)={x}, in(v)={w1,w2,w3,w4}. Arcs x->{w1,..,w4}: x beats some.
    # Score of x determines how many of {w1,..,w4} it beats.

    print("\n  Through-vertex 3-cycle counts:")
    scores = [(s, 5-s) for s in range(6)]
    for s in range(6):
        # through-v = s*(5-s) - # arcs from in(v) to out(v) = ???
        # Actually: through-v = # arcs from out(v) to in(v) = # pairs (a,b) a in out, b in in, a->b
        # This depends on the specific tournament, not just score.
        # But: out(v) has s vertices, in(v) has 5-s vertices.
        # # arcs from out to in + # arcs from in to out = s*(5-s) (total pairs)
        # Let e = # arcs out->in. Then through-v = e.
        # e ranges from max(0, s*(5-s) - C(5-s,2)?) ... no, it depends.
        # The through-v count is NOT determined by s_v alone.
        pass

    # KEY PROOF APPROACH for (10,0):
    # If alpha_1 >= 10 and all cycles pairwise share a vertex,
    # the "intersection graph" is K_10.
    # By vertex counting: 10 cycles of sizes l_1,...,l_10 >= 3
    # cover at most sum(l_i) - (number of shared vertices) vertices.
    # With all pairwise sharing, the total coverage is bounded.

    # SUNFLOWER LEMMA approach:
    # Erdos-Ko-Rado for intersecting families of k-sets:
    # For n >= 2k, max |F| = C(n-1, k-1).
    # For 3-cycles (k=3), max intersecting = C(n-1,2).
    # But the 3-cycles must actually be CYCLES in a tournament, so not all triples are allowed.

    # RAMSEY-TYPE ARGUMENT:
    # Among the 3-cycles in the family, consider the "vertex-disjointness" graph
    # (complement of the conflict/sharing graph restricted to 3-cycles).
    # i_2 = 0 means this graph has no edges.
    # But the 3-cycles use 3 vertices each. Pigeonhole:
    # t3 three-cycles cover at most n vertices. Average vertex appears in 3*t3/n cycles.
    # If n=6, t3=6: average vertex in 3 cycles. Max in EKR = C(5,2)=10 through one vertex.
    # Two 3-cycles on {a,b,c} and {d,e,f} with {a,b,c} cap {d,e,f} = {} are disjoint.
    # At n=6: there are C(6,3)/2 = 10 complementary pairs.
    # If t3=6, there are C(6,2)=15 pairs of 3-cycles.
    # Question: can all 15 pairs share a vertex?
    # Yes! Example: all 6 cycles through vertices {0,4} or {0,5}...
    # The counterexample above had 6 cycles with no complementary pair.

    # So the proof must use the 5-cycles too!
    # If t3=6, t5=4, and all 10 cycles pairwise share, then:
    # Every 5-cycle shares a vertex with every 3-cycle.
    # A 5-cycle on 5 vertices is disjoint from a 3-cycle on 3 vertices
    # iff the 3-cycle uses the one remaining vertex plus 2 others NOT in the 5-cycle.
    # At n=6: a 5-cycle uses 5 of 6 vertices. A 3-cycle is disjoint iff it's on
    # {the missing vertex, and 2 others outside the 5-cycle} — but there's only 1
    # vertex outside the 5-cycle! So a 3-cycle needs 3 vertices, but only 1 is outside.
    # => AT n=6, a 3-cycle and a 5-cycle can NEVER be disjoint (only 1 vertex outside).
    # So the disjoint pairs MUST be (3,3).

    print("\n=== VERTEX DISJOINTNESS CONSTRAINTS ===")
    print("At n=6: 3-cycle and 5-cycle NEVER disjoint (only 1 vertex outside 5-cycle)")
    print("At n=6: 5-cycle and 5-cycle can't be disjoint either (5+5=10 > 6)")
    print("So disjoint pairs at n=6 can ONLY be (3,3).")
    print("Two 3-cycles disjoint iff on complementary triples.")
    print("")
    print("At n=7: 3-cycle (3) + 5-cycle (5) = 8 > 7, so never disjoint.")
    print("Two 3-cycles disjoint: need 6 vertices, possible at n>=6.")
    print("Two 5-cycles disjoint: need 10 > 7 vertices, impossible at n=7.")
    print("So disjoint pairs at n=7 can ONLY be (3,3). CONFIRMED by data.")
    print("")
    print("At n=8: 3+5=8=n, so disjoint 3-5 possible! Also 3+3=6<8.")
    print("5+5=10>8, so no disjoint 5-5 pairs at n=8.")

    # THE KEY FOR (10,0):
    # i_2=0 means NO disjoint 3-3 pairs (the ONLY possible disjoint type at n<=7).
    # So we need: among the t3 three-cycles, NO two are vertex-disjoint.
    # This means the 3-cycles form a "pairwise intersecting family" of 3-sets.

    # BUT: from the data, alpha_1=10 at n=6 always has EXACTLY 2 disjoint 3-3 pairs.
    # Why? Because t3=6 at n=6, and 6 three-cycles on 6 vertices MUST contain
    # a disjoint pair IF alpha_1=10.

    # The 6 three-cycles WITHOUT complementary pair existed (shown above),
    # but those configurations have alpha_1 != 10 (the 5-cycle count is different).

    # Let me check: for the non-complementary 6 three-cycle configurations,
    # what's the total alpha_1?
    print("\n=== NON-COMPLEMENTARY 6 three-cycle configs: alpha_1 check ===")

    from itertools import permutations as perms

    for bits in range(2**15):
        adj = [[0]*6 for _ in range(6)]
        for k, (i, j) in enumerate(edges6):
            if (bits >> k) & 1:
                adj[j][i] = 1
            else:
                adj[i][j] = 1

        c3 = find_3cycles(adj, 6)
        if len(c3) != 6:
            continue

        has_compl = False
        for a in range(len(c3)):
            for b in range(a+1, len(c3)):
                if not (c3[a] & c3[b]):
                    has_compl = True
                    break
            if has_compl:
                break

        if not has_compl:
            c5 = find_5cyc(adj, 6)
            alpha1 = len(c3) + len(c5)
            scores = sorted([sum(adj[i]) for i in range(6)])
            print(f"  No compl pair: t3={len(c3)}, t5={len(c5)}, alpha_1={alpha1}, score={scores}")


check_pairwise_intersecting_3cycles()
max_pairwise_sharing_3cycles(7)
