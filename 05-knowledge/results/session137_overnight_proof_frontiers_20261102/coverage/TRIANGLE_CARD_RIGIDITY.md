# Deleted-triangle cards have a canonical completion

**Status: PROVED, all orders.** This is a structural reduction of the zero-response question in `SYMMETRY_EDIT.md`, not its H-counting solution.

Let Q be any tournament on k>=1 vertices, and put `T=Q[C3,...,C3]`. Strongness of Q is not needed in this note. A homogeneous pair means two vertices having identical orientations to every vertex outside the pair; at order two the condition is vacuous.

## Theorem — a unique pair, and a unique family completion

Delete a vertex v of T, leaving B. Then:

1. B has exactly one homogeneous pair `{u,w}`: the two surviving vertices of v's triangle block.
2. Orient that pair as `u->w`. The deleted vertex's original beaten set is intrinsically recoverable from B as

       S0={u} union (N_B^+(u) minus {w}).                (1)

3. If any completion `B+x_S` is isomorphic to **any** tournament of the form `Q'[C3,...,C3]`, then `S=S0`.

### Proof

The two surviving vertices of v's original block have identical external orientations, so they form a homogeneous pair. Any other candidate pair has a vertex b in an intact triangle block. If its other vertex belongs to the same block, the third vertex distinguishes the pair. If its other vertex a lies outside that block, b's two blockmates have opposite relations to b but identical relations to a, so at least one distinguishes the pair. No other homogeneous pair exists.

To reconstruct v, complete the unique pair cyclically:

    x->u->w->x,

and give x the same external orientation as u and w. This is precisely (1).

Now suppose `B+x_S` belongs, up to isomorphism, to the C3-substitution family. Its new vertex x lies in some triangle block of that representation. Deleting x leaves that block's other two vertices as a homogeneous pair. By the already-proved uniqueness, they must be `{u,w}`. Cyclic orientation inside the block and homogeneity outside force exactly the attachment (1). Thus S=S0. QED.

## Consequences

* Every automorphism of B fixes the ordered pair `(u,w)` and S0. It extends uniquely to an automorphism of T fixing v. Hence

      Aut(B) = Aut(T)_v

  under restriction/extension.
* Two cards `T-v` and `T-w` are isomorphic if and only if v and w belong to the same automorphism orbit of T. Any card isomorphism extends by the canonical completion.
* Nonisomorphic C3-substitution tournaments have disjoint vertex-deletion-card supports. An isomorphism of cards would extend to an isomorphism of their canonical completions.
* A nonoriginal H-restoring cut in the zero-response census is necessarily a **genuinely nonisomorphic** successor and leaves the entire C3-substitution family. It cannot be dismissed as relabelling, a different choice of the deleted root, or another quotient presentation.

For completeness, the original triangle blocks are also intrinsically recognizable in T: they are exactly its homogeneous directed C3 triples. A purported such triple meeting a triangle block in one vertex is distinguished by one of that vertex's two unselected blockmates; a triple meeting a block in two vertices is distinguished by the third. Thus any homogeneous C3 triple must be a whole original block. Consequently Q is recovered by collapsing those blocks, and

    Aut(T) = (C3)^k semidirect Aut(Q).                  (2)

Equation (2) is a structural statement, not a claim that Q is determined by H(T).

## What this does and does not settle

For k>=4, the marked-boundary theorem already excludes H changes +/-2 under any one-vertex replacement. The present theorem says that any nonoriginal zero-response cut is a genuine exit from the structured family. Therefore **isolation of the isomorphism class** in the full one-vertex unit graph is equivalent to the precise scalar statement

    H(B+x_S)=H(T) implies S=S0.                         (3)

The converse direction in (3) is automatic. The difficult direction is OPEN in general, despite the finite censuses. A nonoriginal zero cut would refute class isolation; it would not by itself prove that a later unit path reaches a different H value.

The formal capacity perturbation in `SYMMETRY_EDIT.md` demonstrates why (3) needs actual path-kernel information: positivity, symmetry, parity, divisibility, pair-deletion lower bounds, and all cardinality-layer means still allow alternative zero cuts in a nonrealizable network.

## Exact controls

`card_rigidity_audit.py` checks every labelled quotient through order four, and every strong quotient isomorphism class at orders five and six. It enumerates homogeneous directed triangles in the full substitution, checks every deletion's unique homogeneous pair, and reconstructs the exact labelled original attachment from (1). It imports no inherited repository engine. Run:

    python3 -B card_rigidity_audit.py > card_rigidity_audit.out
    python3 -B -O card_rigidity_audit.py > card_rigidity_audit_optimized.out
    cmp card_rigidity_audit.out card_rigidity_audit_optimized.out
