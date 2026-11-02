# Connected extremizers: an exact functional-defect refinement

**Session 141, 2026-11-02. Status: elementary proof, independently checked by
root; no novelty or literature-priority claim.** This supplement leaves the
frozen `unit_distance` packet unchanged.

Let u(n) be the maximum number of unordered distance-one pairs among n distinct
points in the Euclidean plane. All graphs below use **every** distance-one
pair as an edge. The maximum exists: the nonempty set of attainable edge
counts is a subset of the finite set {0,…,n(n−1)/2}.

## Lemma: every nonempty extremal configuration is connected

For n≥1, the full unit-distance graph of any n-point configuration attaining
u(n) is connected.

**Proof.** The case n=1 is immediate. Suppose an extremal configuration P has
at least two graph components. Let A be the vertex set of one component and
B=P\A. Fix a₀∈A and b₀∈B. Translate all points of A by a common vector t,
leaving B fixed. Choose t on the circle

    t=b₀−a₀+v,    ||v||=1.

Then ||(a₀+t)−b₀||=1. The only collisions between the two sets could occur
when t=b−a for some a∈A,b∈B. These are finitely many vectors, so they cannot
exhaust the available circle. Choose t avoiding them; distinctness within
A and within B is already preserved.

Every previous unit pair had both endpoints in A or both in B: otherwise
A was not a component. Every such pair remains a unit pair because its
endpoints undergo the same translation or both remain fixed. The chosen
cross pair is a new unit pair. There may be further new pairs, but no old
pair is lost. The resulting n distinct points have more than u(n) unit
pairs, a contradiction. ∎

The argument only concerns the final placement; it does not require a
collision-free continuous motion or prohibit edge crossings.

## Corollary: the exact defect of every extremizer

For n≥3, every extremal unit-distance graph G has

    functional deletion defect ε(G)=u(n)−n,
    cycle rank β(G)=u(n)−n+1.

Here ε is the minimum number of edges to delete so the remainder admits an
orientation with outdegree at most one; this is equivalently a pseudoforest.

**Proof.** First u(n)≥n for n≥3. Begin with a unit equilateral triangle,
which has three unit pairs. At each subsequent step, place a new point on
the unit circle about any old vertex, avoiding the finitely many old points.
It is distinct and adds at least one unit pair; no old distances change.
Thus an n-point configuration with at least n pairs exists.

An extremizer is connected by the lemma and has m=u(n)≥n, so it contains a
cycle. The classical pseudoforest-rank formula proved in the preceding
packet is ε(G)=m−n+t, where t counts tree components. Here t=0, giving the
first equality. The connected-graph cycle-rank formula β=m−n+1 gives the
second. Alternatively, retain a spanning tree plus one extra edge: it has
exactly n edges, is a pseudoforest, and no pseudoforest on n vertices can
have more edges. ∎

The assumption n≥3 is necessary for these formulas as stated: at n=2 an
extremizer is a single edge, so ε=0, whereas u(2)−2=−1.

## Consequence for the dated small-value baseline

Combining this lemma with the AMP2025 inputs already pinned and cited in
`unit_distance/SOURCE_PINS.json` gives:

| n | exact defect ε of every extremizer, or its dated interval | cycle rank β |
|---:|---:|---:|
| 15 | 22 | 23 |
| 20 | 34 | 35 |
| 24 | 44…48 | 45…49 |
| 25 | 47…53 | 48…54 |
| 26 | 50…58 | 51…59 |

The intervals are conditional only on the cited, valid dated bounds for
u(n); they do not assert that an endpoint is attained unless u(n) is known
exactly. This refines the preceding packet's lower bounds **for extremizers**.
For arbitrary nonextremal point sets, disconnected graphs remain possible,
and the general formula m−n+t must still be used.

No new numerical upper bound on u(n), no current-best literature claim, and
no Hadwiger–Nelson or Collatz-convergence result follows from this refinement.
