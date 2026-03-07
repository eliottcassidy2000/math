# Good Deletion at n >= 10: Proof Attempt

**Instance:** opus-2026-03-07-S43
**Status:** IN PROGRESS — close to complete

## Statement

**Theorem:** Every cycle-rich tournament T on n >= 10 vertices has a vertex v
such that T-v is cycle-rich on n-1 vertices.

## Computational Evidence

- n=10: 0 no-good-deletion in 96M cycle-rich tournaments (100M samples)
- n=9: only 20 no-good-deletion in 93M cycle-rich (score (1,1,1,4,4,4,7,7,7), mm=3)
- n=9 failures have ALL vertices uniquely bottlenecked with sum|B(v)|=9

## Proof Framework

### Obstruction Analysis

A deletion of vertex w fails to preserve cycle-richness if:
(A) **Score obstruction:** Some u becomes source (score(u,T)=1, u→w) or sink (score(u,T)=n-2, w→u)
(B) **Cycle obstruction:** Some u ∈ B(w) (all u's 3-cycles go through w)

### Landau's Theorem Constraint

By Landau's theorem, a valid score sequence (s_1 <= ... <= s_n) must satisfy
sum_{i=1}^k s_i >= C(k,2) for all k.

For k vertices with score 1: sum = k >= C(k,2) = k(k-1)/2.
This gives k >= k(k-1)/2, so k-1 <= 2, so k <= 3.

**Similarly:** at most 3 vertices can have score n-2 (by the complementary Landau condition).

So |S_1| <= 3 and |S_{n-2}| <= 3, giving at most 6 score obstructions.

### Bottleneck Function Properties

**Lemma (Unique Bottleneck):** Each vertex u has at most one bottleneck vertex v.
If u has 3-cycles through both v and v', neither is a bottleneck.

**Corollary:** sum_w |B(w)| <= n.

**Lemma (Transitive Bottleneck Set):** B(w) induces a transitive subtournament.
No triple in B(w) forms a 3-cycle (else it would be a 3-cycle not through w).

### Counting Argument

Vertex w is a valid deletion candidate if:
1. No score-1 vertex targets w (at most 3 blocked by source creation)
2. No score-(n-2) vertex is beaten by w (at most 3 blocked by sink creation)
3. B(w) = ∅ (no vertex depends exclusively on w for cyclicity)

**Condition 1:** Score-1 vertices S_1 = {u : score(u)=1}. Each u ∈ S_1 has exactly
one out-neighbor. Call it target(u). Deleting target(u) creates a source.
So at most |S_1| <= 3 vertices are blocked by condition 1.

**Condition 2:** Score-(n-2) vertices S_{n-2} = {u : score(u)=n-2}. Each u ∈ S_{n-2}
is beaten by exactly one vertex. Call it beater(u). Deleting beater(u) creates a sink.
So at most |S_{n-2}| <= 3 vertices are blocked by condition 2.

**Condition 3:** Vertices w with B(w) ≠ ∅.
sum_w |B(w)| <= n, so at most n vertices have non-empty B(w).
But how many vertices w actually have B(w) ≠ ∅?

### Key Claim: |{w : B(w) ≠ ∅}| <= n-4 at n >= 10

Actually, we need the opposite: at most few w's have non-empty B(w), leaving room
for a good deletion.

**Alternative approach:** Consider the COMPLEMENT set: vertices w with B(w) = ∅.
Call these "safe" (from condition 3). A safe w is a good deletion candidate
if it also avoids conditions 1 and 2.

Number of safe vertices = n - |{w : B(w) ≠ ∅}|.

At n=9 with no-good-deletion: |{w : B(w) ≠ ∅}| = 6, safe = 3.
But all 3 safe vertices are blocked by score obstructions (they're the score-4 vertices
that are either targeted by score-1 or beating score-7).

At n >= 10: we need > 6 - |score_blocked| safe vertices.

### The n >= 10 Argument

**At n=9 (the critical case):**
- 3 score-1, 3 score-7 vertices block 6 candidates
- 6 bottleneck centers block... wait, the bottleneck centers ARE different vertices
  from the score-blocked ones. The score-1 and score-7 vertices have specific targets.

Actually, looking at the computational data: at n=9, the 3 safe vertices (B(w)=∅)
are all blocked by score obstructions. The 6 bottleneck centers + 3 safe = 9.
The 3 safe vertices each have either a score-1 vertex targeting them or they beat
a score-7 vertex.

**At n=10:** Adding one more vertex breaks this tight balance.

Suppose by contradiction n=10 has no good deletion. Then ALL 10 vertices are blocked.
At most 6 are blocked by score obstructions. So at least 4 must have B(w) ≠ ∅.
But sum|B(w)| <= 10, and the bottleneck vertices must be distinct.

Consider the "bottleneck graph": directed from u to bottleneck(u).
This is a partial function from V to V with at most n edges.
At n=9 no-good-del: this function is TOTAL (every vertex has a bottleneck) with
exactly 9 edges (a 9-vertex functional graph).

At n=10 with all 10 blocked: need at least 4 vertices as bottleneck centers.
sum|B(w)| >= 4 (at least 1 each). The remaining 10 - sum|B(w)| vertices are
bottleneck-free but score-blocked.

But score blocking requires: at most 3 from S_1 targets, at most 3 from S_{n-2} beaters.
Some targets/beaters might be in bottleneck centers too (double-blocked).

**The tight counting:**
- |blocked_by_score| <= 6 (from S_1 and S_{n-2})
- |blocked_by_cycle| >= 10 - 6 = 4 (need B(w) ≠ ∅ for these 4)
- sum|B(w)| >= 4 (at least 4 bottleneck assignments)
- But also: |{u : u has bottleneck}| = sum|B(w)| <= 10

The question is: can we have sum|B(w)| = 10 (all vertices bottlenecked)
AND |S_1| = 3, |S_{n-2}| = 3 at n=10?

**Landau constraint at n=10:** Can score sequence have 3 ones and 3 eights?
[1,1,1,s_4,...,s_7,8,8,8] with sum = 45. Used: 3+24=27. Remaining 4: sum=18, avg=4.5.
Check Landau:
k=3: 1+1+1=3 >= 3 ✓
k=4: 3+s_4 >= 6, so s_4 >= 3 ✓
...

Score [1,1,1,3,5,5,5,8,8,8]: sum = 3+3+15+24 = 45 ✓
Landau k=4: 1+1+1+3=6 >= 6 ✓
k=5: 6+5=11 >= 10 ✓
k=7: 21 >= 21 ✓

This is valid. But does a cycle-rich tournament with this score sequence exist?
And if so, can all vertices be bottlenecked?

**The key constraint:** If ALL 10 vertices are bottlenecked, then the bottleneck
graph is a total function V → V (every vertex mapped). This creates a very specific
functional graph structure. At n=9, this graph has 3 "trees" of depth 1 (center
with 2 children).

At n=10, can a 10-vertex functional graph support the required structure?
Each vertex u maps to its bottleneck v, meaning u's ONLY 3-cycles go through v.
This means {u, v, a} is a 3-cycle for some a, and ALL other triples through u are
transitive (except those also through v).

If bottleneck(u) = v and bottleneck(v) = u (mutual bottleneck), then u and v have
ONLY the 3-cycle {u, v, a} between them, and their ONLY 3-cycles all go through
each other. This means u and v are in exactly one 3-cycle (the one containing a),
and a is a "bridge" vertex.

But a must also be in a 3-cycle (cycle-rich). If a's 3-cycles all go through
some bottleneck w, then {a, w, b} is a 3-cycle for some b...

This mutual dependence creates an extremely rigid structure that seems impossible
at n >= 10 with the added vertex providing too many "escape routes" for cyclicity.

## Status

The proof is close but needs a clean formalization of why the functional graph
structure can't exist at n >= 10. The key is the interaction between Landau's
constraint on extreme scores and the bottleneck transitivity requirement.
