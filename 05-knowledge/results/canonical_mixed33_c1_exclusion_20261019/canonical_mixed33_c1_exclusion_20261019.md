# Canonical mixed (3,3,1) maximum-span root exclusion

Date: 2026-10-19. Completed exact exclusion; the fresh independent terminal replay passed. The six-thin theorem is an explicit external dependency.

## 1. Scoped theorem and complete core catalogue

Suppose a finite poset P has maximum incomparability degree at most seven and no balanced pair. In its strict >2/3 majority order, choose an incomparable root pair x before y with maximum index span, and among all such pairs maximize c=|I(x)∩I(y)|. This canonical rule is invariant under order duality. The conclusion is that this canonically selected root cannot have lower-exclusive count 3, upper-exclusive count 3 and common count 1; the dual statement follows. This is a canonical-root claim, not an exclusion of every maximum-span pair with this morphology.

Every maximum-span pair oriented x before y is critical. If z<x but z is not below y, then z is incomparable with y and precedes x in the majority order, giving a larger incomparable span; thus D(x)⊆D(y). The successor argument gives U(y)⊆U(x). Define L=D(y)\D(x), T=U(x)\U(y), and Z=I(x)∩I(y). The complete root core is Q={x,y}∪L∪T∪Z, containing every vertex incomparable with either root.

Label x=0,y=1, lower vertices 2,3,4, common vertex 5, and upper vertices 6,7,8. Mandatory relations are 2<y,3<y,4<y and x below all three upper vertices. The roots and common vertex are otherwise incomparable within these sets. Every optional core relation is a nonroot relation directed forward in the lower/common/upper block order. Reverse cross-block relations would compare the roots or compare the common vertex with a root. Internal lower and upper orders are arbitrary; natural labeling within each block loses no isomorphism type.

All 96,428 naturally labeled orders on the seven nonroots are generated. Only the common vertex can exceed degree seven in the nine-point core: it has degree eight exactly when isolated among the nonroots. Removing it leaves an arbitrary six-point naturally labeled poset, so exactly 4,824 cores are invalid. There are 91,604 degree-valid cores. Independent block-ideal generation checks equality of the complete core sets, and every analytic exclusion is independently replayed.

## 2. Exact exterior defect screens with forced-comparability closure

Every vertex outside Q={x,y}∪L∪T∪Z lies in A below both roots or B above both. We have A<B, and no core vertex can lie below A or above B. Any core vertex of degree seven is saturated and comparable with all exterior vertices, so A<s<B.

Let S consist of the two roots and all saturated core vertices. Define

    F_A={v : v=s or s<v for some s∈S},
    F_B={v : v=s or v<s for some s∈S}.

Then every A lies below F_A, and every B lies above F_B, by transitivity. In particular upper vertices belong to F_A and the lower vertices belong to F_B. Put s(v)=7−d_Q(v), α(v)=0 on F_A and s(v) otherwise, and β(v)=0 on F_B and s(v) otherwise. For any incomparable core u,v,

    δ_P(u,v)≤δ_Q(u,v)+α(u)+β(v)=H(u,v).

The extra predecessor difference is bounded by missing A-relations of u; the extra successor difference by missing B-relations of v. These are upper bounds valid for every actual exterior, even though some combined slack assignments are not realizable.

For incomparable u,v define the directed defect

    δ(u,v)=|D(v)\D(u)|+|U(u)\U(v)|.

The exchange lemma says δ(u,v)≤1 implies p(u<v)≤2/3. For completeness: when δ=0, swapping u and v injects extensions with u before v into those with v before u. If the only defect is an extra successor t of u, split the former extensions into u<v<t and u<t<v. Swapping u,v works on the first class. On the second, replace the three occupied labels (u,t,v) by (v,u,t). This is valid because D(v)⊆D(u) and every successor of t is a successor of v: any successor of t is a successor of u distinct from its sole extra successor t. Each class injects separately into v<u. Thus the forward class has at most twice the reverse count. A sole predecessor defect follows by order duality.

If H(u,v)≤1 and H(v,u)≤1, the pair is balanced and its entire core is excluded. Otherwise any projected majority order containing u before v with H(u,v)≤1 is excluded.

For completeness, an additional analytic core screen records the following necessary directed arcs: all comparabilities, x→y, x→z→y for common z, and v→u whenever H(u,v)≤1. Under no balance these all belong to the strict majority order. Any cycle excludes the core. The balanced-defect screen excludes 10,690 cores. The forced-cycle screen gives no additional exclusions in this family. For an upper vertex a, the correctly computed bound is H(a,y)=|L\D_Q(a)|+|U_Q(a)|. In particular, maximality of a alone does not make this bound at most one when multiple lower exclusive vertices remain incomparable with a.

## 3. Why the projected order test is necessary

With no balanced pair, every incomparable pair has exactly one direction of probability >2/3. Together with comparable pairs, these directions form a transitive tournament: a directed triangle would have three cyclic probabilities summing to more than two, whereas each linear extension satisfies at most two of the cyclic indicators. Its order M extends P.

All exterior A letters precede x and y, and all B letters follow both. Consequently no exterior letter is between the roots, and their span in M restricted to Q equals their full span. Projection cannot increase another core-pair span. Therefore the projected order must be a linear extension of Q with

    rank(y)>rank(x),
    |rank(v)−rank(u)|≤rank(y)−rank(x)

for every incomparable core pair u,v. All such projected orders are enumerated exactly. This is only a necessary condition, so retaining extra unrealizable orders is harmless.

### Explicit dependency: known six-thin theorem

If a candidate root span is at most three, it is also the full global maximum span, by the exterior argument above. Every incomparable neighbor of any vertex must then lie within three positions on either side of it in M, so the entire ambient poset has incomparability degree at most six.

Marcin Peczarski, *The Gold Partition Conjecture for 6-Thin Posets*, Order 25 (2008), 91–103, DOI https://doi.org/10.1007/s11083-008-9081-9, proves the Gold Partition Conjecture for every six-thin poset; that theorem implies the 1/3–2/3 balance conclusion. We therefore exclude all candidate spans≤3 using this known result. This is an explicit external theorem dependency, not a new six-thin proof. The publisher's primary abstract states both the six-thin theorem and the implication to balance.

The independent replay repeats the projected-span test and checks the number of excluded orders. It does not purport to reprove Peczarski's theorem. All larger-span cases below are handled by exact finite certificates.

### Additional canonical maximum-span screen

If an incomparable core pair s,t has projected span equal to the root span d, then its full ambient span is also d: projection only shrinks spans and d is the global maximum. No exterior vertex lies between s and t. Every common incomparable neighbor of a maximum-span pair must be between its endpoints, since otherwise it gives a longer incomparable pair with one endpoint. Therefore

    I_P(s)∩I_P(t)=I_Q(s)∩I_Q(t).

If its core common count exceeds the root count 1, this contradicts canonical selection. Such an order is discarded before cone certification. The independent replay repeats this test and counts its exclusions exactly. This is the [canonical tie lemma](../canonical_maxspan_transfer_20261019/canonical_maxspan_transfer_20261019.md), reproduced here for completeness.

No interval-transfer catalogue is used in this scan. In particular, this exclusion is not used as its own premise. The scope remains canonical(maximum span, maximum common count); no extra degree-sum tie-break is assumed.

## 4. Exact boundary-generator lift

For a core extension o, let pre(o) precede its first root and suf(o) follow its last root. Define L(u) as the number of interleavings of A with the fixed ordered word u, divided by e(A), and R(v) dually using B. Then the number of ambient extensions restricting to o is

    e(A)e(B)L(pre(o))R(suf(o)).

All A is before the first root, all B after the last. All other core/exterior relations are automatically respected outside these two boundary interleavings, proving this factorization.

Let U_empty=V_empty=1 and, for nonempty words,

    U_u=L(u)−L(u with last letter removed),
    V_v=R(v)−R(v with first letter removed).

These numbers are nonnegative: append the removed prefix letter, or prepend the removed suffix letter, to obtain an injection. Expanding along the two word trees gives a nonnegative linear combination of generators G_(u,v), where G_(u,v) counts the core extensions beginning with u and ending with v and having both roots between these boundary words.

If d_Q(s)=7, all exterior points are comparable with s, hence A<s<B. Every increment whose boundary word contains s vanishes, including increments later along that word-tree branch: after the first saturated letter in a prefix, all A has already appeared, and the normalized prefix count remains constant; the suffix statement is dual. Thus generators whose prefix or suffix includes a saturated vertex may be removed. No independence or realizability of the other generator coefficients is assumed.

For a retained generator g let E_g count its extensions, and N_g(a<b) count those having a before b. An exact certificate consisting of positive integer weights w_ab on forward incomparable pairs proves a contradiction whenever

    sum_ab w_ab [3 N_g(a<b)−2 E_g] ≤ 0

for every retained g. Multiplying by the actual nonnegative boundary coefficients and summing proves the same inequality for ambient extension counts. Yet every certificate pair is forward in M, so each ambient summand is strictly positive. This is impossible.

## 5. Shared slack-one exterior disjunction

If a core vertex q has incomparability degree six, at most one exterior vertex is incomparable with q. Consequently either every A vertex lies below q or every B vertex lies above q, possibly both. These two cases cover all actual exteriors. They must be treated as a union, rather than taking the convex hull of their coefficient cones.

In the first case all prefix increments whose word contains q vanish; the same holds for every vertex above q by transitivity. In the second case all suffix increments containing q or a vertex below q vanish. This is the [one-slack exterior split lemma](../one_slack_exterior_split_20261019/one_slack_exterior_split_20261019.md). It does not assume that every retained coefficient vector is realizable.

The continuation uses this rule only if the ordinary generator cone has no certificate. A complete binary tree records each selected degree-six vertex and both cases. Each leaf records the accumulated forced-A and forced-B closures, the exact retained generator indices, and a positive integer weighted comparison certificate valid on all those generators. Every leaf must be certified. Both children are required at every internal node; no unsupported case is discarded. Repeated splits only restrict the exterior possibilities, and induction up the finite tree establishes the contradiction for every actual exterior. A cached tree may be reused for another projected order only if every comparison appearing anywhere in its leaves is forward in that order.

The independent standard-library verifier reconstructs every extension and generator, checks the degree-six hypothesis and both children, recomputes both transitive closures and exact masks, and checks each integer margin. It also re-enumerates the complete set of projected orders and verifies that every retained order has a valid ordinary certificate or a complete split tree.

## 6. Why the refinement was needed

The unrefined scan stopped at core index 62484, with predecessor masks

    (0,28,0,0,12,4,13,93,125)

and projected order (2,3,0,4,5,6,1,7,8). This core has 540 extensions and degrees (5,5,2,3,3,6,3,3,2). Its root span is 4; the root is the only tied maximum-span incomparable pair, so earlier secondary-pair transfer gives no improvement. No core antichain satisfies the saturated-antichain pattern of packet 141.

A seven-generator rational probability mixture gives every forward incomparable comparison probability strictly greater than 2/3, with minimum signed margin 355/2515202. This is an obstruction to the unrestricted nonnegative generator relaxation, not an actual poset counterexample. Independent brute-force enumeration of all 9! permutations verifies the540 extensions, every conditional event count, and every rational margin.

Actual boundary extension weights factor as L(pre)R(suf). On the four realized boundary cells with prefix empty or[2] and suffix empty or[8], the mixture has nonzero determinant

    -22174681/56936169907236.

Therefore this particular extension-weight distribution is not realizable by an actual exterior. More decisively, vertex 5 has degree 6. In the A<5 branch the forward comparisons (6,1),(4,5),(5,6), with weights 2,1,2, give an integer generator margin at most−2 on all 112 retained generators. In the5<B branch the comparisons (0,4),(4,5),(5,6), with weights 2,2,1, do the same. Thus the actual exterior possibilities are excluded despite the failure of the convex relaxation. All 8 residual orders of this core are independently certified after the split refinement.

## 7. Three-phase provenance and scope

The original 54,355 successful core records remain immutable and are independently replayed by their original verifier. The second phase adds 2,876 successful cores using the slack-one fallback, then stops at the second relaxation obstacle. The third phase adds 23,683 cores using the ancestor-cap fallback and completes the catalogue. Each continuation has a separate record schema and checkpoint, hash-locks all prior streams and sources, and independently replays each successful record before committing it.

The terminal verifier regenerates the natural catalogue independently, invokes each phase's own immutable exact verifier, and checks that the three ordered core-index lists form the complete disjoint union of all 80,914 residual cores. It verifies each source lock, committed byte count, phase histogram and final record hash. The historical stopped checkpoints and two-phase prefix audit are retained as provenance; neither is presented as a complete morphology certificate. The fresh three-phase replay passed; its receipt is `complete_replay.json`.

This result is restricted to a canonical root selected by maximum span and then maximum common-neighbor count. It does not exclude every maximum-span root of this morphology and does not establish the full seven-thin balance conjecture.

## Attribution

The six-thin theorem is due to Peczarski, explicitly cited above. Majority-order, critical-pair, exchange-injection and boundary-extension methods belong to the established 1/3–2/3 and Gold-partition literature; no novelty claim is made for those general methods. The finite receipts certify the precise scoped statement and do not substitute for their mathematical hypotheses. The one-slack split is separately stated and proved in its linked dependency. No canonical(3,3,1) exclusion is used as its own premise.

## 8. Ancestor caps for actual boundary increments

The first continuation stopped at core index 65627, with predecessor masks (0,28,0,0,12,12,13,125,125), the same projected order (2,3,0,4,5,6,1,7,8), and 408 extensions. Its common vertex has degree four, so the slack-one disjunction is unavailable. Indeed, even a factorized boundary model has all forward core probabilities greater than 2/3, as separately documented by the rank-one obstruction result. Thus rank-one structure alone is insufficient here.

Actual interleaving increments satisfy stronger ancestor inequalities. For a prefix wq, deleting q from an interleaving counted by its positive increment leaves an interleaving counted by the preceding increment, and reinserting q has at most the number of A vertices incomparable with q possible positions. Consequently

    U_(wq) ≤ (7−d_Q(q)) U_w.

The same argument reversed gives V_(qv) ≤ (7−d_Q(q)) V_v. Empty ancestors have normalized increment one. These inequalities remain valid when the actual missing-A and missing-B bounds are separately enlarged to the full core slack. A detailed injection and the precise empty-word case are given in the [boundary increment cap lemma](../boundary_increment_caps_20261019/boundary_increment_caps.md).

For raw generator coefficients λ_(u,v)=U_u V_v, each supported child therefore satisfies

    λ_child ≤ s λ_parent,  s=7−d_Q(q),

where the prefix parent deletes its last letter, or the suffix parent deletes its first letter. Every supported child has a supported parent. Removal of saturated-boundary generators preserves this hereditary property. These inequalities concern raw λ, not the normalized mixture mass E_g λ_g. The search normalization is Σ_g E_g λ_g=1.

A cap certificate records positive integer weights on forward comparison rows W_g=Σ_ab w_ab(3N_g(a<b)−2E_g) and nonnegative integer multipliers t_e such that, coefficientwise,

    W + Σ_e t_e(s_e e_parent−e_child) ≤ 0.

Each cap expression is nonnegative on actual λ. Multiplying by actual λ and summing thus proves W·λ≤0, contradicting strict majority on every weighted comparison. This is a linear valid-inequality certificate; no factorization search or exhaustive exterior enumeration is needed.

The independent replay reconstructs the complete ordered generator list, checks each cap's side, parent, child and degree-derived slack, and evaluates every residual coefficient using Python integers. It also checks that every comparison weight is positive and every cap multiplier is positive, rejects duplicate cap rows, and independently checks the raw-coefficient convention. Cached certificates are reused only when all their comparison orientations are forward; cap inequalities are independent of the majority order.

The cap fallback precedes slack-one branching. If a cap-constrained LP suggests feasibility, its raw coefficients are rationally reconstructed and every cap and strict comparison margin is checked exactly before the object is called an obstruction. Numerical failure or failed rational reconstruction is an explicit error. Such a model would still be only a relaxation witness, not an actual poset counterexample.

For completeness, fix an extension a₁,…,aₘ of A. An extra interleaving counted by U_(wq) has an A letter after its final core letter q. Delete q. If w is nonempty, the image still has an A letter after all of w and is counted by U_w; if w is empty, it is the fixed A extension. If aⱼ is the last A predecessor of q (j=0 if absent), every later A letter is incomparable with q. Thus m−j≤s_A(q). Reinsertion before the final A letter and after every predecessor permits at most m−j A gaps; requiring q after all of w only reduces the options, with at most one insertion position per gap. Summing this bound over A extensions and dividing by e(A) proves the increment inequality. If A is empty, all nonempty increments vanish. Reversing the order proves the suffix inequality.

The third phase preserves the first two stopped streams, their source files, obstacle records and stopped-prefix audit by exact hashes. The complete terminal replay uses each phase's own immutable verifier and independently checks the full ordered disjoint union against the regenerated catalogue. The fresh terminal replay passed on all three phases in 1,440.7 seconds.

## 9. Complete coverage accounting

Independent block-ideal enumeration produces 96,428 natural nonroot orders and 91,604 degree-valid full cores. The defect screen excludes 10,690 whole cores. The remaining 80,914 cores are covered by three disjoint immutable certificate phases of 54,355, 2,876 and 23,683 cores.

Across those phases, 47,300,375 projected orders are excluded by directed-defect bounds, 28,176 by the explicit six-thin dependency, and 332,621 by the canonical common-count screen. The remaining 515,442 projected orders receive exact certificates: 440,021 single-row, 69,528 equal-two-row, 5,879 general weighted, 10 complete slack-one trees and 4 ancestor-cap certificates. These are method counts for certificates found, not claims about minimum support. The canonical screen entirely removes 15,209 otherwise nonempty residual cores.

The fast source/coverage audit verifies 71,579 distinct cached certificates and the exact phase union. Every successful record was independently replayed before commitment. A fresh complete replay separately regenerates the natural catalogue and all extensions/generators and rechecks every exact inequality. The final receipt `complete_replay.json` records `FULL_THREE_PHASE_INDEPENDENT_EXACT_REPLAY_PASSED` and agrees with every phase hash and count.
