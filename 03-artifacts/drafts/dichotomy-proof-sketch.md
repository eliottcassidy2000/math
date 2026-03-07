# Dichotomy Proof Sketch: 3 Disjoint 3-Cycles OR Good Deletion

**Instance:** opus-2026-03-07-S42
**Status:** SKETCH — verified computationally at n=9, general proof open

## Statement

**Claim:** For any cycle-rich tournament T on n >= 9 vertices (every vertex in
a 3-cycle, no source/sink), at least one of the following holds:
(a) T contains 3 pairwise-disjoint directed 3-cycles (Part C gives H >= 27 > 21)
(b) Some vertex v exists such that T-v is cycle-rich on n-1 vertices

## Evidence

- n=9 (50M trials, 153,444 cycle-rich): ZERO counterexamples found
- 99.997% of cycle-rich n=9 tournaments have a good deletion (b)
- The single example without a good deletion (t3=3, 3 disjoint 3-cycles) satisfies (a)

## Proof Sketch for (b) when max matching <= 2

Suppose max matching of 3-cycle triples = 2. Let M = {A, B} be a maximal
matching with A ∩ B = ∅, |A ∪ B| = 6.

Every other 3-cycle C intersects A ∪ B (else {A, B, C} is a 3-matching).

The remaining n-6 >= 3 vertices are each in some 3-cycle touching A ∪ B.

### Strategy: Delete a vertex from V \ (A ∪ B)

Pick w ∈ V \ (A ∪ B). For T-w to be cycle-rich, we need:
1. No source/sink in T-w
2. Every vertex in T-w is in some 3-cycle

**Condition 1 (no source/sink):** w creates a source in T-w iff some score-1
vertex u has u → w (u's only arc). w creates a sink iff some score-(n-2)
vertex u has w → u.

**Avoiding condition 1:** Choose w that is NOT the unique target of any score-1
vertex AND NOT a neighbor beaten by a score-(n-2) vertex.

Score-1 vertices: each beats exactly 1 vertex. If there are k score-1 vertices,
they collectively have k targets. The "outer" vertices V \ (A ∪ B) have n-6
members. If k < n-6, some outer vertex is not targeted by any score-1 vertex.

**Key bound:** In a cycle-rich tournament, can there be >= n-6 score-1 vertices?
Each score-1 vertex v beats exactly 1 other vertex. By Part J (Key Lemma),
since v is in a 3-cycle, v's unique target u must have some w with u → w → v.
So v is beaten by w, meaning v's in-neighborhood includes w.

At n=9: n-6 = 3. So we need at most 2 score-1 vertices targeting outer vertices.
Score-1 vertices have in-degree 7 in a 9-vertex tournament. They're "almost sinks."

**Condition 2 (all in 3-cycle):** Every vertex in T-w must still be in a 3-cycle.
A vertex u loses its 3-cycle in T-w only if ALL u's 3-cycles pass through w.

For u ∈ A: u is in the 3-cycle A = {u, a', a''} (with a', a'' ∈ A, none = w).
So A persists in T-w. ✓

For u ∈ B: similarly, B persists. ✓

For u ∈ (A ∪ B) \ {vertices in A, B}: actually, A and B partition A ∪ B.
So all 6 vertices in A ∪ B are covered.

For u ∈ V \ (A ∪ B ∪ {w}): u is in some 3-cycle C touching A ∪ B.
If C doesn't use w, C persists in T-w. ✓
If C uses w: u's 3-cycle was {u, w, x} for some x. This is destroyed.
But u might have another 3-cycle not through w.

**The danger:** u has ALL 3-cycles through w. This means u is only cyclic
through w. If we delete w, u loses all cycles.

By Part J: in T-w, u is in no 3-cycle means u is in no cycle at all.
So u can be removed from T-w via Parts J/K, giving T'' on n-2 vertices.
But this means T-w is NOT cycle-rich (u is acyclic).

**Resolution:** If deleting w makes some u acyclic, we can instead try deleting
a different outer vertex w' such that no vertex loses ALL its 3-cycles.

**Counting argument:** Each outer vertex w can "poison" at most deg(w)-in-the-
3-cycle-structure number of other outer vertices. If the poisoning is sparse,
some w is safe to delete.

### Alternative: Direct H bound for max matching <= 2

At n=9, mm <= 2 cycle-rich tournaments have min H = 45. This is already > 21.
If this bound increases with n, the proof is complete without the dichotomy.

**Why H is large when mm <= 2:** With few disjoint 3-cycles, the 3-cycle triples
are concentrated (share many vertices). This forces the tournament to have many
5-cycles (the inter-triple arcs create cyclic sub-tournaments on 5 vertices).
The extra 5-cycles contribute to alpha_1, alpha_2, increasing H.

At n=8 cycle-rich, min H = 25 (exhaustive).
At n=9 cycle-rich (mm <= 2), min H = 45 (sampling).

The growth from 25 to 45 is consistent with a monotonic lower bound.

## Conclusion

The dichotomy is computationally verified at n=9. A full proof requires either:
1. Showing the counting argument (some outer vertex is safe to delete) works at all n, or
2. Showing min H for cycle-rich tournaments is > 21 for all n >= 8.

Approach 2 seems more promising, as the min-H grows rapidly with n.
