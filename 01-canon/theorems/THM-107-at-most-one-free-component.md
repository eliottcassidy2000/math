# THM-107: At Most One Free Component

**Status:** PROVED (algebraically)
**Filed by:** kind-pasteur-2026-03-08-S43
**Depends on:** THM-103 (TT span), THM-104 (cycle sum equality), THM-105 (dominant vertex forcing), THM-106 (free cycle bridge)

## Statement

In any tournament T on n >= 4 vertices, the 3-cycle adjacency graph (where two directed 3-cycles are adjacent iff they share a directed edge) has at most one connected component consisting entirely of free (non-dominated) cycles.

A cycle C = {a,b,c} is **dominated** if some external vertex d satisfies d->a, d->b, d->c (dominates all) or a->d, b->d, c->d (dominated by all). A cycle is **free** if it is not dominated.

## Corollary

**b_1(T) <= 1** for all tournaments T, where b_1 = dim H_1(T) in GLMY path homology.

This follows because:
- b_1 = #{free components} (by THM-103 + THM-104 + THM-105)
- #{free components} <= 1 (this theorem)

## Proof

Suppose for contradiction that C_1 and C_2 are free cycles in different connected components of the 3-cycle graph.

### Case 1: C_1 and C_2 share >= 2 vertices.
Two vertices determine a unique directed edge in a tournament. Both cycles contain this edge. So they share a directed edge and are in the same component. Contradiction.

### Case 2: C_1 and C_2 share exactly 1 vertex.
Let C_1 = {a,b,c} with a->b->c->a, and C_2 = {a,d,e} (sharing vertex a).

By THM-106 (free cycle bridge theorem), since C_1 is free, each external vertex creates a bridge cycle sharing an edge with C_1. In particular, vertices d and e each create bridges.

The bridge from v to C_1 shares one of the three edges {a->b, b->c, c->a} with C_1:
- Bridge via a->b: B_v = {v,a,b}. Contains a.
- Bridge via b->c: B_v = {v,b,c}. Does NOT contain a.
- Bridge via c->a: B_v = {v,c,a}. Contains a.

**Subcase 2a:** At least one of B_d, B_e contains a.
Say B_d = {d, a, w} for some w in {b,c}. Then B_d and C_2 share vertex pair {a,d}, hence share the directed edge between a and d. B_d is in C_1's component (shares edge with C_1). So C_2 and C_1's component share an edge. Contradiction.

**Subcase 2b:** Neither B_d nor B_e contains a.
Then both bridge via b->c: B_d = {d,b,c} and B_e = {e,b,c}.
The bridge cycle d->b->c->d requires d->b and c->d.
The bridge cycle e->b->c->e requires e->b and c->e.
From C_1: a->b.

So vertex b satisfies: a->b, d->b, e->b. That is, b is beaten by all three vertices of C_2 = {a,d,e}. Therefore C_2 is dominated (by vertex b, which loses to all of C_2). Contradiction with C_2 being free.

### Case 3: C_1 = {a,b,c} and C_2 = {d,e,f} are vertex-disjoint.

By THM-106, since C_1 is free, each of d,e,f creates a bridge to C_1. Each bridge B_v = {v, u, w} contains v and 2 vertices from {a,b,c}, sharing an edge with C_1.

Similarly, since C_2 is free, each of a,b,c creates a bridge to C_2. Each bridge B_u = {u, x, y} contains u and 2 vertices from {d,e,f}, sharing an edge with C_2.

Define bipartite "containment" graphs:
- G on {a,b,c} x {d,e,f}: edge (x,y) iff y in B_x (the bridge from x to C_2)
- H on {d,e,f} x {a,b,c}: edge (y,x) iff x in B_y (the bridge from y to C_1)

Each vertex in G has degree 2 (bridge uses 2 of 3 target vertices). So |E(G)| = 6.
Each vertex in H has degree 2. So |E(H)| = 6.

Consider H^T (transposed: edges (x,y) for (y,x) in H) as a graph on {a,b,c} x {d,e,f}.
|E(H^T)| = 6. The complement of G has |{a,b,c} x {d,e,f}| - |E(G)| = 9 - 6 = 3 edges.

By inclusion-exclusion: |G intersect H^T| >= |E(G)| + |E(H^T)| - 9 = 6 + 6 - 9 = 3.

So there exist at least 3 pairs (x,y) with x in C_1, y in C_2 such that:
- y in B_x (bridge from x to C_2 contains y)
- x in B_y (bridge from y to C_1 contains x)

Pick any such pair (x,y). Then:
- B_x = {x, y, ?_1} shares an edge with C_2 (by construction)
- B_y = {y, x, ?_2} is in C_1's component (shares edge with C_1)
- B_x and B_y share vertex pair {x,y}, hence the directed edge between x and y
- Therefore B_x is in C_1's component (connected to B_y via shared edge)
- Therefore C_2 is in C_1's component (connected to B_x via shared edge)

This contradicts C_1 and C_2 being in different components. QED.

## Verification

Exhaustive: n=4 (64/64), n=5 (1024/1024), n=6 (32768/32768).
Sampled: n=7 (2000/2000), n=8 (1000/1000), n=9 (500/500), n=10 (200/200).
Max free components always = 0 or 1.

## See Also
- THM-102 (beta_2 proof status)
- THM-103 (TT boundaries span im(d_2))
- THM-106 (free cycle bridge theorem)
- HYP-279 (b_1 <= 1, now PROVED as corollary)
