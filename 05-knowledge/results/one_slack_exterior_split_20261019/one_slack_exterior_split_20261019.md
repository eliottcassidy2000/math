# One-slack exterior splitting resolves a concrete mixed (3,3,1) cone obstruction

Date: 2026-10-19. Status: root-reviewed and frozen. This is a theorem for one specified core and projected order, plus a general structural lemma. It is not a whole-morphology exclusion or a seven-thin theorem.

## 1. Ambient setting

Let P be a finite poset of maximum incomparability degree at most seven. Let x,y be a critical incomparable pair, oriented x before y, and let Q be the full root core: x,y and every vertex incomparable with either root. Every exterior vertex lies in A below both roots or B above both roots, and A<B. No core vertex can be below an A vertex or above a B vertex: that would make it comparable in the forbidden direction to both roots. Consequently all core/exterior pairs have the directions A<q, q<B, or incomparability.

For q in Q, its total exterior incomparabilities are at most 7−d_Q(q). In particular, if d_Q(q)=6, then

    either every a in A satisfies a<q,
    or every b in B satisfies q<b.

Indeed, failure of both alternatives produces distinct exterior vertices a in A and b in B, each incomparable with q, exceeding its one remaining incomparability. Both alternatives can hold. This is a disjunction about one actual ambient poset, not a permission to mix the two cases arbitrarily.

More generally, at incomparability bound k, every core vertex of degree k−1 satisfies the same disjunction. A degree-k vertex satisfies both alternatives.

## 2. Shared exterior letters and the boundary increment support

For a core extension o, let its prefix before the first root be pre(o), and its suffix after the last root be suf(o). Let L(u) be the number of ways to interleave A with the fixed word u respecting P, divided by e(A). Let R(v) be the analogous normalized B interleaving count. Then the number of ambient extensions restricting to o equals

    e(A)e(B)L(pre(o))R(suf(o)).

Every A is before the first root and every B after the last, so the other core/exterior constraints are automatically respected. Define U_empty=V_empty=1 and

    U_u=L(u)−L(u with its last letter deleted),
    V_v=R(v)−R(v with its first letter deleted).

Appending the last prefix letter, or prepending the first suffix letter, gives injections, so these increments are nonnegative.

A useful stronger support fact is that U_u>0 for a nonempty word u requires one shared exterior vertex a in A incomparable with EVERY letter of u. To see this, a merged extension counted by the increment is outside the append-letter injection. Its last core letter is therefore followed by some A letter a. Every letter of u precedes a in that merged extension. None can be below an A vertex in P, so each is incomparable with a. Dually, V_v>0 requires a shared b in B incomparable with every letter of v.

It follows immediately that if all A<q then every U_u with q in u is zero. If all q<B then every V_v with q in v is zero. If d_Q(q)=6, any pair of boundary words both containing q satisfies U_u V_v=0: otherwise the two shared letters just identified are distinct exterior incomparable neighbors of q.

The one-slack disjunction is stronger for cone certification than deleting only the overlapping support pairs: all generators containing q on ONE entire side vanish simultaneously, with the side fixed for the actual P. Thus the actual boundary data belong to a union of two nonnegative cones, rather than their convex hull.

With several one-slack vertices, one may branch on these alternatives separately. Forced A-comparability propagates upward in Q; forced B-comparability propagates downward. This gives a finite union-of-cones refinement. No sufficiency or ambient realizability claim is made for any resulting cone.

## 3. The concrete core and its two certificates

Use labels x=0,y=1 and predecessor masks

    (0,28,0,0,12,4,13,93,125).

The required projected strict majority order is

    (2,3,0,4,5,6,1,7,8).

The lower-exclusive set is {2,3,4}, the upper-exclusive set is {6,7,8}, and the common set is {5}. The core has exactly 540 linear extensions. The incomparable neighbors of 5 are {0,1,3,4,6,7}, so d_Q(5)=6. There is no degree-seven core vertex. The projected maximum incomparable span is four, attained by the root pair. The proof below does not need maximum-span or canonical tie-breaking assumptions; the specified projected order suffices.

For a pair of boundary words (u,v), let E_(u,v) count core extensions beginning with u and ending with v, with both roots between the words. Let N_(u,v)(a<b) count these extensions satisfying a<b. Expanding L and R along their word trees expresses ambient extension counts as the nonnegative sum of these columns, with coefficients e(A)e(B)U_u V_v.

Case A: every A<5. Delete every generator whose prefix contains 5. Exactly 112 generators remain. Direct exact enumeration verifies, for EACH such generator,

    2[3N(6<1)−2E] + [3N(4<5)−2E]
      + 2[3N(5<6)−2E] ≤ −2.

Case B: every 5<B. Delete every generator whose suffix contains 5. Again exactly 112 generators remain. For EACH such generator,

    2[3N(0<4)−2E] + 2[3N(4<5)−2E]
      + [3N(5<6)−2E] ≤ −2.

All listed pairs are incomparable and forward in the specified projected order. Multiplying the appropriate generator inequality by actual nonnegative boundary coefficients and summing proves that at least one listed comparison has ambient probability at most 2/3. This contradicts the requirement that every forward incomparable comparison in a strict majority order have probability strictly above 2/3.

The disjunction in Section 1 exhausts every allowed exterior, so no actual seven-thin ambient poset can have this full root core and this projected strict majority order. This resolves the specific relaxed-cone obstruction without assuming uniform core probabilities or realizing a formal generator mixture.

## 4. Independent fixed-certificate verification and scope

`04-computation/experiments/verify_one_slack_exterior_split_20261019.py` uses only the Python standard library. It reconstructs the 540 extensions, all boundary generators, both 112-column restrictions and every integer margin from the masks. In each case the integer margins range from −160 to −2. It checks transitivity, the common vertex, degree six, the projected order and its maximum span, and forward incomparability of all certificate pairs. Guards are explicit and remain active under optimization. Normal and optimized outputs are identical.

The certificate coefficients are the small fixed integers 2,1,2 and 2,2,1. No numerical solver is used by the verifier. A numerical discovery experiment is scratch only and unnecessary to verify the result. No broad core census, Git action, unproved packet139 exclusion, or secondary morphology transfer is used.

The elementary root-core and boundary-factorization setup agrees with packet126; the proof is reproduced here so the fixed-core result does not depend on that packet's complete morphology scan. Packet123 explains the canonical-root context but is not required for this exclusion. Packets141 and142 concern different antichain restrictions and are not used. No historical novelty claim is made.

Run `python 04-computation/experiments/verify_one_slack_exterior_split_20261019.py` from the repository root. The recorded outputs are `verification.json` and `verification_optimized.json` alongside this note.
