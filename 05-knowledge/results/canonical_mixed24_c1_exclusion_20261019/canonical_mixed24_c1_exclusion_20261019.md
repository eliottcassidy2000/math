# Canonical mixed (2,4,1) maximum-span root exclusion

Date: 2026-10-19. Complete canonical-root exclusion, supported by exact certificates, independent catalogue verification and independent per-record replay. The six-thin theorem is an explicit external dependency.

## 1. Scoped theorem and complete core catalogue

Suppose a finite poset P has maximum incomparability degree at most seven and no balanced pair. In its strict >2/3 majority order, choose an incomparable root pair x before y with maximum index span, and among all such pairs maximize c=|I(x)∩I(y)|. This canonical rule is invariant under order duality. The conclusion is that this canonically selected root cannot have lower-exclusive count 2, upper-exclusive count 4 and common count 1; the dual statement follows. This is a canonical-root claim, not an exclusion of every maximum-span pair with this morphology.

Every maximum-span pair oriented x before y is critical. If z<x but z is not below y, then z is incomparable with y and precedes x in the majority order, giving a larger incomparable span; thus D(x)⊆D(y). The successor argument gives U(y)⊆U(x). Define L=D(y)\D(x), T=U(x)\U(y), and Z=I(x)∩I(y). The complete root core is Q={x,y}∪L∪T∪Z, containing every vertex incomparable with either root.

Label x=0,y=1, lower vertices 2,3, common vertex 4, and upper vertices 5,6,7,8. Mandatory relations are 2<y,3<y and x below all four upper vertices. The roots and common vertex are otherwise incomparable within these sets. Every optional core relation is a nonroot relation directed forward in the lower/common/upper block order. Reverse cross-block relations would compare the roots or compare the common vertex with a root. Internal lower and upper orders are arbitrary; natural labeling within each block loses no isomorphism type.

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

For completeness, an additional analytic core screen records the following necessary directed arcs: all comparabilities, x→y, x→z→y for common z, and v→u whenever H(u,v)≤1. Under no balance these all belong to the strict majority order. Any cycle excludes the core. The balanced-defect screen excludes 21,501 cores. The forced-cycle screen gives no additional exclusions in this family. For an upper vertex a, the correctly computed bound is H(a,y)=|L\D_Q(a)|+|U_Q(a)|. In particular, maximality of a alone does not make this bound at most one when multiple lower exclusive vertices remain incomparable with a.

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

## 5. Complete exact finite result

The independent catalogue contains 96,428 naturally labeled cores, of which 4,824 are over-degree and 91,604 are valid. The balanced-defect screen excludes 21,501 cores; no further core is excluded by a forced-majority cycle. All 70,103 residual cores are covered by the projected-order procedure.

The complete order counts are:

- Known-six-thin short-span exclusions: 12,716.
- Directed-defect exclusions: 37,396,742.
- Additional canonical common-count exclusions: 274,819.
- Remaining orders with exact positive certificates: 415,884.

Among the 70,103 residual cores, 17,760 already have no retained order before the canonical screen. The new screen entirely removes 15,711 additional cores that previously had at least one retained order. The remaining 36,632 cores have positive residual order counts and are handled by cone certificates. The canonical rule removes 274,819 of 690,703 orders surviving the earlier order screens; no whole-core coverage is inferred merely from a bounded sample.

Of the 415,884 certified orders, 333,056 use a stored single-row certificate, 77,386 use an equal-weight two-row certificate, and 5,442 use a positive integer certificate reconstructed from an LP suggestion. There are 55,622 distinct stored certificates. The largest coefficient is 299,659,835,262,193,117,992; products and sums use arbitrary-precision Python integers. The numerical LP is only a suggestion mechanism. Every coefficient is converted to an exact positive integer and every retained generator margin is checked exactly.

A found higher-support certificate does not imply a lower bound on the support of all alternative certificates. No support-minimality claim is made. No exact cone obstruction occurred.

## 6. Independent coverage and replay

The search adjoins maximal nonroot labels with predecessor ideals. The independent catalogue instead chooses the internal lower/common/upper orders separately, attaches predecessor ideals across the three blocks, and checks all transitivity restrictions. It verifies equality of the complete degree-valid core sets, not only their cardinalities. All 21,501 balanced-defect exclusions are independently checked.

For every residual core, a standard-library verifier independently reconstructs extensions, side bounds, maximum-span orders, the known-six-thin filter, and the canonical common-neighbor screen. It independently reconstructs the retained boundary generators and verifies every integer certificate inequality and deterministic order-assignment digest. Only after that replay succeeds is the record committed to the checkpointed gzip stream. Terminal union verification checks every core ID, record identity, source hash and count, with no skipped or duplicated residual core.

Run `python -O audit.py` for the fast complete coverage/provenance audit. Run `python -O replay_complete.py` for a fresh complete standard-library mathematical replay, including independent catalogue generation. The fresh replay is slower; its terminal receipt is `complete_replay.json`. Every mathematical record was already independently checked before commit, and the separate complete replay provides an additional end-to-end check. Explicit guards remain active under `-O`.

`run.py` is the checkpointed search entrypoint and needs NumPy/SciPy. The other search/replay source modules are helpers; use the documented entrypoints rather than legacy specialized entrypoints inside imported helper files. `progress.json` is terminal, `independent_cover.json` records the independent catalogue/analytic verification, `audit_receipt.json` records terminal provenance checks, and `certificates.jsonl.gz` contains the complete exact record stream. Source and record hashes are explicit in the receipts.

## 7. Dependencies, scope and attribution

The canonical tie lemma from [the prior structural packet](../canonical_maxspan_transfer_20261019/canonical_maxspan_transfer_20261019.md) is proved again in Section 3. Its dependency-aware interval-transfer procedure is not used here. No other morphology exclusion, and certainly not this result itself, is used to eliminate a secondary pair. There is no circular theorem dependency.

The only external mathematical theorem needed beyond the elementary reductions given here is the explicitly cited Peczarski six-thin theorem for short spans. Classical context includes G. Brightwell and C. D. Wright, *The 1/3–2/3 Conjecture for 5-Thin Posets* (1992), https://doi.org/10.1137/0405037, and E. J. Olson and B. E. Sagan, *On the 1/3–2/3 Conjecture*, https://arxiv.org/abs/1706.04985. General majority-order, exchange, extension-counting and positive-certificate methods are not claimed as original. These references provide context, not a claim that they contain this exact finite receipt. No novelty claim is made.

The theorem excludes only a root chosen by the rule (maximum span, maximum common count), in morphology (2,4,1) or its dual. Earlier any-root results remain stronger in their respective scopes and remain valid for a canonical root. A global future proof may choose a canonical root initially, but must still cover all its remaining morphologies. This is not an any-root theorem for (2,4,1), not a full mixed(2,4) branch theorem, and not a full seven-thin theorem. Other larger branches remain outside this certificate.
