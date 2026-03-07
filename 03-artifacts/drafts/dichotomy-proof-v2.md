# Dichotomy Proof v2: Pigeonhole on Cyclic Neighborhoods

**Instance:** opus-2026-03-07-S43
**Status:** PROOF ATTEMPT — in progress

## Statement

**Claim (Dichotomy):** For any cycle-rich tournament T on n >= 9 vertices, at least one of:
(a) T has 3 pairwise-disjoint directed 3-cycles, or
(b) Some vertex v exists such that T-v is cycle-rich on n-1 vertices.

## Setup

T is cycle-rich: no source/sink (scores in [1, n-2]), every vertex in a 3-cycle.

**Definition:** For vertex v, let C(v) = set of 3-cycle vertex sets containing v.
For a 3-cycle {a,b,c}: it appears in C(a), C(b), and C(c).

**Definition:** Vertex v is a *bottleneck* for vertex u if ALL 3-cycles containing u
also contain v. Equivalently: u has no 3-cycle in T-v.

**Definition:** B(v) = {u : v is a bottleneck for u} = set of vertices that lose ALL
3-cycles when v is deleted.

## Key Lemma: Bottleneck Structure

**Lemma 1:** If v is a bottleneck for u, then in T-v, vertex u is in no directed cycle
of any length (by Key Lemma, Part J of THM-079).

**Lemma 2:** Each vertex u has at most one bottleneck v. (If u has 3-cycles through v
and through v', and these are different vertices, then neither is a bottleneck.)

**Proof of Lemma 2:** If v and v' are both bottlenecks for u, then ALL 3-cycles of u go
through v AND through v'. So every 3-cycle of u contains both v and v'.
But a 3-cycle has only 3 vertices: {u, v, v'}. So u's ONLY 3-cycle vertex set
is {u, v, v'}. This means u, v, v' form a directed 3-cycle with no other 3-cycles
involving u.

For {u, v, v'} to be u's only 3-cycle: every triple {u, a, b} with {a,b} != {v,v'}
must be transitive. This severely constrains u's neighborhood structure.

**Lemma 3:** B(v) is an independent set in the "3-cycle graph" restricted to V\{v}.
That is: no two vertices in B(v) share a 3-cycle not through v.
(Trivially true since both have ALL 3-cycles through v.)

**Stronger:** If u1, u2 ∈ B(v), then {u1, u2, a} is NOT a 3-cycle for any a ≠ v.
Because u1's only 3-cycles go through v, and {u1, u2, a} doesn't contain v (since a≠v).
Wait: {u1, u2, v} could be a 3-cycle, but {u1, u2, a} for a≠v cannot.

So T[B(v)] restricted to any triple that doesn't include v is transitive.
Since B(v) ⊂ V\{v}, any triple {u1, u2, u3} ⊂ B(v) is transitive.
This means T[B(v)] is TRANSITIVE (no 3-cycles among bottleneck vertices of v).

**Lemma 4:** |B(v)| <= something small.
If T[B(v)] is transitive and every u ∈ B(v) has ALL 3-cycles through v,
then u ∈ B(v) has 3-cycle {u, v, a} where a ∉ B(v) or a = v (impossible).
So each u ∈ B(v) requires a witness vertex a_u with {u, v, a_u} being a 3-cycle
and a_u ∉ B(v).

Can different u's share the same witness? If a is the witness for both u1 and u2:
{u1, v, a} and {u2, v, a} are both 3-cycles. So u1, u2, v, a are 4 vertices
containing 2 directed 3-cycles both using v and a. The sub-tournament on {u1,u2,v,a}
has at least 2 three-cycle vertex sets. Since T[{u1,u2}∪B(v)\{u1,u2}] is transitive
within B(v), we need u1 → u2 (or u2 → u1) without forming 3-cycles among B(v).

## Key Counting Argument

**Total bottleneck assignments:** sum_v |B(v)| <= n (by Lemma 2: each u has at most 1 bottleneck).

**Vertex deletion v fails (for cycle-richness) if:**
1. Some u in T-v becomes a source or sink [score obstruction], or
2. Some u in T-v is in no 3-cycle [cycle obstruction = u ∈ B(v)]

For (1): u becomes source iff score(u)=1 and u→v. u becomes sink iff score(u)=n-2 and v→u.
Let S1 = {u : score(u)=1} and S_{n-2} = {u : score(u)=n-2}.
|S1| vertices each block their unique target.
|S_{n-2}| vertices each block the unique vertex that beats them.

**Deletion v is "good" if:**
- No u ∈ S1 has u→v (v not the target of any score-1 vertex), AND
- No u ∈ S_{n-2} has v→u (v not beating any score-(n-2) vertex), AND
- B(v) = ∅ (no vertex loses all 3-cycles)

**Blocked vertices for deletion:**
- From (1): each score-1 vertex blocks 1 target, each score-(n-2) vertex blocks 1 other vertex.
  Total: |S1| + |S_{n-2}| vertices blocked by score obstruction.
- From (2): each v with B(v) ≠ ∅ is blocked.

**No good deletion means:** for ALL v ∈ V, at least one of the three conditions fails.

**Claim:** If ALL n vertices are blocked, then |S1| + |S_{n-2}| + #{v : B(v)≠∅} >= n.
But #{v : B(v)≠∅} <= n (trivially). And sum_v |B(v)| <= n.

The score obstructions block at most |S1| + |S_{n-2}| candidates (pigeonhole on targets).
The cycle obstructions block at most sum_v 1_{B(v)≠∅} candidates (which equals #{v : B(v)≠∅}).

**For the contrapositive:** If #{non-blocked} >= 1, then a good deletion exists.
#{non-blocked} >= n - |S1| - |S_{n-2}| - #{v : B(v)≠∅}.

So no good deletion implies: |S1| + |S_{n-2}| + #{v : B(v)≠∅} >= n.

**Key insight for n=9:**
By Landau's theorem, score sequence must satisfy partial sum conditions.
The maximum number of score-1 and score-7 vertices at n=9 is constrained.
Specifically: at most 3 vertices can have score 1 (since sum of 4 smallest >= C(4,2)=6,
but 4 scores of 1 sum to 4 < 6). Similarly at most 3 vertices with score 7.
So |S1| + |S_{n-2}| <= 6 at n=9.

This means #{v : B(v)≠∅} >= n - 6 = 3.
So at least 3 vertices have non-empty bottleneck sets.

sum_v |B(v)| <= n = 9. With at least 3 non-empty B(v)'s,
average |B(v)| <= 3 for the non-empty ones.

## Connection to 3-Disjoint Cycles

If no good deletion exists AND max matching <= 2:
Every vertex is "trapped" — deleting it either creates a source/sink or kills someone's cyclicity.

With at least 3 vertices having B(v) ≠ ∅, there are at least 3 "bottleneck centers."
Each center v has B(v) consisting of vertices whose ONLY 3-cycles go through v.
The bottleneck centers are each "hubs" of 3-cycle stars.

If the bottleneck centers are pairwise disjoint from each other's 3-cycles,
we get 3 disjoint 3-cycles (one from each center's star). This gives (a).

The challenge: prove the bottleneck centers can't all overlap.

## TO DO

1. Show that bottleneck centers v with |B(v)| >= 2 must be in "dense" 3-cycle stars
2. Show that 3 such centers with sufficient B(v) sizes force 3 disjoint 3-cycles
3. Handle the edge cases where bottleneck sets are small

This proof is incomplete but suggests the right framework.
The bottleneck function B(v) and Lemma 2 (unique bottleneck) are the key tools.
