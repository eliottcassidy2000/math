#!/usr/bin/env python3
"""
Why do 3 cyclic triples span exactly 5 vertices when n=6 and c3=3?

The key constraint: at n=6, if the 3 triples span 6 vertices, they must
use ALL vertices. With 3 triples on 6 vertices, each triple has 3 vertices,
union = 9 vertex-slots for 6 vertices. By inclusion-exclusion, some vertices
appear in multiple triples.

Pattern analysis: what configurations of 3 triples on 6 vertices are possible
such that EXACTLY these 3 are cyclic and no others?

At n=7 the disjoint configuration works: {a,b,c} disjoint from {d,e,f},
with {a,b,d} sharing. But 7 vertices needed.

At n=6 with 6 vertices used: 9 vertex-slots, 6 vertices -> 3 overlaps.
Possible intersection patterns:
- (0,0,2): two share 2 vertices, third disjoint. Span = 3+3+1=7? No: 3+3+3-2=7.
  Wait: T1={a,b,c}, T2={a,b,d}, T3={e,f,g}. Span=7. Can't fit in 6 vertices.
  With 6 vertices: T1={a,b,c}, T2={a,b,d}, T3={e,f,...}. Need T3 from {c,d,e,f}
  since only 6 vertices. If T3={c,d,e}, span=6. If T3={e,f,...}, need f in {a,..,f}.
  Actually with vertices {0,1,2,3,4,5}:
  T1={0,1,2}, T2={0,1,3}, T3={4,5,x} where x in {0,1,2,3}.
  If x=2: T3={2,4,5}. Span=6.

This is possible in principle. Does it create a 4th cyclic triple?
Let's check: with T1={0,1,2} cyclic and T2={0,1,3} cyclic and T3={2,4,5} cyclic,
we need to check that no other triple from C(6,3)=20 is cyclic.

The argument from the two-sharing analysis:
  T1={0,1,2} cyclic: say 0->1->2->0
  T2={0,1,3} cyclic: 0->1 already fixed. So 1->3->0 (to make cyclic with 0->1).
  Now: 2->0 and 3->0 (both go to 0). Also 1->2 and 1->3.
  Edge 2 vs 3: if 2->3, then {0,2,3} has 2->3->0->... and 2->0? No: 2->0.
    So 2->3, 3->0, 0->2? That's 0->2, but we have 2->0 from T1. Contradiction!
    Actually wait: T1 has 0->1->2->0. So 2->0. And 3->0 from T2.
    {0,2,3}: edges are 2->0, 3->0, and 2 vs 3.
    If 2->3: 2->3, 3->0, 0->... 2->0 means it's NOT 0->2. So {0,2,3} = 2->3->0, 2->0.
    That's vertex 2 beating both. Transitive. Not a cycle. GOOD.
    If 3->2: {0,2,3} = 3->2, 2->0, 0->... and 3->0.
    3 beats both 2 and 0. Transitive. GOOD.

  Now T3={2,4,5} cyclic: say 2->4->5->2.
  Edges involving 4 and 5 vs {0,1,3}:
  We need to ensure no other triple is cyclic.

  Consider {0,4,5}: need to assign 0 vs 4, 0 vs 5.
  Consider {1,4,5}: need 1 vs 4, 1 vs 5.
  Consider {3,4,5}: need 3 vs 4, 3 vs 5.
  Consider {0,2,4}: 2->4 from T3, 2->0 from T1. 0 vs 4?
  Consider {0,2,5}: 5->2 from T3, 2->0 from T1. 0 vs 5?
  Consider {1,2,4}: 1->2 from T1, 2->4 from T3. 1 vs 4?
  ... many constraints.

Let me just check this computationally.

kind-pasteur-2026-03-06
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations

def c3_from_scores(scores):
    n = len(scores)
    return (n*(n-1)*(n-2)//6) - sum(s*(s-1)//2 for s in scores)

def get_cyclic_triples(T, n):
    result = []
    for combo in combinations(range(n), 3):
        a, b, c = combo
        if T[a][b] and T[b][c] and T[c][a]:
            result.append(frozenset(combo))
        elif T[a][c] and T[c][b] and T[b][a]:
            result.append(frozenset(combo))
    return result

# At n=6: Can 3 triples span 6 vertices?
# We showed exhaustively they can't. But WHY?
# Let's try to construct such a tournament.
print("=== Trying to construct n=6 tournament with c3=3 spanning 6 vertices ===")
print()

# Configuration: T1={0,1,2}, T2={0,1,3}, T3={2,4,5}
# T1 cyclic: 0->1->2->0
# T2 cyclic with 0->1: must be 0->1->3->0, so 3->0 and 1->3.
# T3 cyclic: 2->4->5->2
# Fixed edges: 0->1, 1->2, 2->0, 1->3, 3->0, 2->4, 4->5, 5->2
# Free edges: 0 vs 4, 0 vs 5, 1 vs 4, 1 vs 5, 2 vs 3, 3 vs 4, 3 vs 5

# For each of 2^7 = 128 assignments, check c3
n = 6
free_edges = [(0,4), (0,5), (1,4), (1,5), (2,3), (3,4), (3,5)]
fixed = {
    (0,1): 1, (1,0): 0,
    (1,2): 1, (2,1): 0,
    (2,0): 1, (0,2): 0,
    (1,3): 1, (3,1): 0,
    (3,0): 1, (0,3): 0,
    (2,4): 1, (4,2): 0,
    (4,5): 1, (5,4): 0,
    (5,2): 1, (2,5): 0,
}

found = False
for mask in range(1 << len(free_edges)):
    T = [[0]*n for _ in range(n)]
    # Set fixed edges
    for (i,j), val in fixed.items():
        T[i][j] = val
    # Set free edges
    for idx, (i,j) in enumerate(free_edges):
        if mask & (1 << idx):
            T[i][j] = 1
            T[j][i] = 0
        else:
            T[j][i] = 1
            T[i][j] = 0

    triples = get_cyclic_triples(T, n)
    c3 = len(triples)
    if c3 == 3:
        all_verts = set()
        for t in triples:
            all_verts.update(t)
        span = len(all_verts)
        if span == 6:
            print(f"  FOUND span=6: triples={[sorted(t) for t in triples]}")
            found = True
            break

if not found:
    print("  Config T1={0,1,2},T2={0,1,3},T3={2,4,5}: cannot achieve c3=3 with span=6")
    print("  (always creates extra cyclic triples)")

# Try other configurations
print()
print("Trying all possible 3-triple configurations spanning 6 vertices at n=6...")
# Generate all sets of 3 triples from C(6,3)=20 triples that span all 6 vertices
from itertools import combinations as comb
all_triples_6 = list(comb(range(6), 3))
configs_tested = 0
configs_possible = 0

for t1, t2, t3 in comb(range(len(all_triples_6)), 3):
    trip1 = frozenset(all_triples_6[t1])
    trip2 = frozenset(all_triples_6[t2])
    trip3 = frozenset(all_triples_6[t3])
    union = trip1 | trip2 | trip3
    if len(union) != 6:
        continue  # Want span = 6
    configs_tested += 1

    # Try to build a tournament where exactly these 3 are cyclic
    # For each triple to be cyclic, try both orientations
    # Check if any assignment of all edges gives exactly c3=3

    # This is expensive but n=6 has only 15 edges, with some fixed by the triples
    # Let's just check: over all 2^15 tournaments, does any have exactly these 3 cyclic?
    # That's too slow for each config. Instead, use the constraint.

    # For each of 2^3=8 orientation choices for the 3 triples:
    possible = False
    for orient in range(8):
        # Orient each triple
        edges = {}  # (i,j) -> 1 means i beats j
        conflict = False

        for idx, trip in enumerate([sorted(trip1), sorted(trip2), sorted(trip3)]):
            a, b, c = trip
            if orient & (1 << idx):
                # a->b->c->a
                required = [(a,b,1),(b,c,1),(c,a,1)]
            else:
                # a->c->b->a
                required = [(a,c,1),(c,b,1),(b,a,1)]
            for i,j,v in required:
                key = (min(i,j), max(i,j))
                direction = 1 if i < j else -1
                if v == 1 and i < j:
                    val = 1  # i beats j
                elif v == 1 and i > j:
                    val = 0  # min beats max? No: j beats i means j<i, so (j,i) with j beating i -> val = 1 for (j,i)
                # Simpler: just record who beats whom
                if (i,j) in edges:
                    if edges[(i,j)] != 1:
                        conflict = True
                        break
                elif (j,i) in edges:
                    if edges[(j,i)] != 0:
                        conflict = True
                        break
                else:
                    edges[(i,j)] = 1
                    edges[(j,i)] = 0
            if conflict:
                break

        if conflict:
            continue

        # Determine fixed edges and free edges
        fixed_edges = {}
        for (i,j), v in edges.items():
            fixed_edges[(i,j)] = v

        # Count free edges
        all_pairs = [(i,j) for i in range(6) for j in range(i+1,6)]
        free = []
        for (i,j) in all_pairs:
            if (i,j) not in fixed_edges and (j,i) not in fixed_edges:
                free.append((i,j))

        # Try all 2^|free| assignments
        for fmask in range(1 << len(free)):
            T = [[0]*6 for _ in range(6)]
            for (i,j), v in fixed_edges.items():
                T[i][j] = v
            for fidx, (i,j) in enumerate(free):
                if fmask & (1 << fidx):
                    T[i][j] = 1
                else:
                    T[j][i] = 1

            # Verify tournament (each pair has exactly one edge)
            ok = True
            for i in range(6):
                for j in range(i+1,6):
                    if T[i][j] + T[j][i] != 1:
                        ok = False
                        break
                if not ok:
                    break
            if not ok:
                continue

            triples = get_cyclic_triples(T, 6)
            if len(triples) == 3:
                all_v = set()
                for t in triples:
                    all_v.update(t)
                if len(all_v) == 6:
                    configs_possible += 1
                    print(f"  POSSIBLE: triples={[sorted(t) for t in triples]}")
                    possible = True
                    break
        if possible:
            break

print(f"\nConfigurations spanning 6 vertices tested: {configs_tested}")
print(f"Configurations that can be realized with c3=3: {configs_possible}")

if configs_possible == 0:
    print("CONFIRMED: At n=6, no tournament with c3=3 has triples spanning 6 vertices.")
    print()
    print("This means: at n=6, c3=3 forces triples to span exactly 5 vertices,")
    print("so there is a 'spectator' vertex not in any triple.")

print("\nDone.")
