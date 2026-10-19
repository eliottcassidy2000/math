# Mixed one-lower/four-upper maximum-span exclusion with two common neighbors

Date: 2026-10-19. Status: proved exact morphology theorem with arbitrary actual exterior. All residual certificates were independently replayed and the complete union audit passed. No full seven-thin theorem or novelty claim.

## 1. Theorem and complete core catalogue

Suppose a finite poset P of incomparability degree at most seven has no balanced pair. Let x precede y be a maximum-index-span incomparable pair in its strict >2/3 majority order. The following critical-pair morphology is impossible when c=2:

    D(x)⊆D(y), |D(y)\D(x)|=1,
    U(y)⊆U(x), |U(x)\U(y)|=4,
    c=|I(x)∩I(y)|.

The order-dual statement holds. Label the roots x=0,y=1, the lower vertex b=2, common vertices 3,...,c+2 and upper vertices c+3,...,c+6. The mandatory relations are b<y and x<a for each upper a. Apart from these, all core relations are among nonroots, directed forward in the block order lower/common/upper. Internal common and upper orders are arbitrary. Any reverse relation across the blocks would compare a common vertex with a root or compare x and y. Conversely every naturally labeled order on the 5+c nonroot vertices gives a valid core. Natural labeling within each block loses no isomorphism type.

There are 96,428 naturally labeled orders on the seven nonroot vertices. The only over-budget core vertices are isolated common vertices, whose incomparability degree is eight. With one specified common vertex isolated there are 4,824 remaining nonroot orders, and with both isolated there are 357. Hence inclusion-exclusion removes 2·4,824−357=9,291 cores, leaving 87,137 degree-valid cores.

The scan generates this catalogue by adjoining maximal labels with arbitrary predecessor ideals. The independent verifier generates common/upper orders, lower-to-common upsets, common-to-upper ideals and lower-to-upper upsets, enforcing transitivity. It checks equality of the complete core sets and all degree counts.

## 2. Exact exterior defect screens with forced-comparability closure

Every vertex outside Q={x,y,b}∪Z∪T lies in A below both roots or B above both. We have A<B, and no core vertex can lie below A or above B. Any core vertex of degree seven is saturated and comparable with all exterior vertices, so A<s<B.

Let S consist of the two roots and all saturated core vertices. Define

    F_A={v : v=s or s<v for some s∈S},
    F_B={v : v=s or v<s for some s∈S}.

Then every A lies below F_A, and every B lies above F_B, by transitivity. In particular upper vertices belong to F_A and the lower vertex belongs to F_B. Put s(v)=7−d_Q(v), α(v)=0 on F_A and s(v) otherwise, and β(v)=0 on F_B and s(v) otherwise. For any incomparable core u,v,

    δ_P(u,v)≤δ_Q(u,v)+α(u)+β(v)=H(u,v).

The extra predecessor difference is bounded by missing A-relations of u; the extra successor difference by missing B-relations of v. These are upper bounds valid for every actual exterior, even though some combined slack assignments are not realizable.

For incomparable u,v define the directed defect

    δ(u,v)=|D(v)\D(u)|+|U(u)\U(v)|.

The exchange lemma says δ(u,v)≤1 implies p(u<v)≤2/3. For completeness: when δ=0, swapping u and v injects extensions with u before v into those with v before u. If the only defect is an extra successor t of u, split the former extensions into u<v<t and u<t<v. Swapping u,v works on the first class. On the second, replace the three occupied labels (u,t,v) by (v,u,t). This is valid because D(v)⊆D(u) and every successor of t is a successor of v: any successor of t is a successor of u distinct from its sole extra successor t. Each class injects separately into v<u. Thus the forward class has at most twice the reverse count. A sole predecessor defect follows by order duality.

If H(u,v)≤1 and H(v,u)≤1, the pair is balanced and its entire core is excluded. Otherwise any projected majority order containing u before v with H(u,v)≤1 is excluded.

For completeness, an additional analytic core screen records the following necessary directed arcs: all comparabilities, x→y, x→z→y for common z, and v→u whenever H(u,v)≤1. Under no balance these all belong to the strict majority order. Any cycle excludes the core. Here the balanced-pair defect screen excludes 28,972 cores; the forced-majority cycle screen excludes an additional 3,255 cores. The remaining 54,910 cores are handled by the projected-order and boundary-generator procedure below. For an upper vertex a, H(a,y)=1_{b not<a}+|U_Q(a)|, because A<a and all B lie above y. Thus every maximal upper vertex is forced after y (H(a,y)≤1); an upper vertex with one upper successor is also forced after y when b<a. No such conclusion is asserted from this screen when b∥a and a has one upper successor.

## 3. Why the projected order test is necessary

With no balanced pair, every incomparable pair has exactly one direction of probability >2/3. Together with comparable pairs, these directions form a transitive tournament: a directed triangle would have three cyclic probabilities summing to more than two, whereas each linear extension satisfies at most two of the cyclic indicators. Its order M extends P.

All exterior A letters precede x and y, and all B letters follow both. Consequently no exterior letter is between the roots, and their span in M restricted to Q equals their full span. Projection cannot increase another core-pair span. Therefore the projected order must be a linear extension of Q with

    rank(y)>rank(x),
    |rank(v)−rank(u)|≤rank(y)−rank(x)

for every incomparable core pair u,v. All such projected orders are enumerated exactly. This is only a necessary condition, so retaining extra unrealizable orders is harmless.

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

## 5. Complete exact certificate

All 87,137 valid cores are covered:28,972 balanced-defect exclusions,3,255 forced-cycle exclusions and 54,910 residual cores. In the residual cores,20,923,177 necessary projected orders are excluded by directed-defect bounds. The remaining 519,617 orders have exact positive comparison certificates:392,383 single-row orders,122,562 two-row orders and 4,672 orders using integer certificates reconstructed from LP-dual suggestions. There are 52,059 distinct certificates.

Every numerical suggestion was converted into positive integer weights and checked on every retained generator in exact Python integer arithmetic. The largest weight is 19,367,567,609,030,429,021,940. The independent verifier uses only standard-library arbitrary-precision integer operations. No floating-point inequality establishes an accepted certificate.

## 6. Independent catalogue, replay and deterministic coverage

The independent verifier regenerated all 96,428 cores by a different common/upper ideal-upset construction, checked exact equality of the 87,137 valid-core set and verified all 32,227 analytic exclusions. All 54,910 residual core IDs were partitioned deterministically into four disjoint shards of 13,728,13,728,13,727,13,727. Each shard independently replayed each complete core record before committing it. All four have terminal independently-replayed success.

`audit_union.py` verifies the exact partition, every core ID and predecessor-mask record, terminal status, source and record hashes, and all numerical totals. There is an intentionally empty serial prefix; the prefix pause/partial labels describe that empty prefix, not the combined terminal result. `union_receipt.json` is the authoritative terminal combined receipt.

Run `python -O audit_union.py` for quick coverage/provenance verification. Run `python -O replay_full.py` for a fresh complete mathematical replay with no numerical dependencies. The latter regenerates the catalogue and repeats every certificate check, so it is slower. The recorded independent mathematical checks already cover every analytic exclusion and every shard record. Explicit guards remain enabled under -O.

The release contains the proof, complete catalogue/analytic witnesses, four certificate streams with replay-bearing receipts, independent verifier, reproducible search and hash manifest. The proof does not impose any exterior-size bound.

## 7. Full mixed(1,4) corollary and limits

Here the root y has incomparability degree 5+c, so c≤2. Together with the [c≤1 theorem](../mixed14_c01_exclusion_20261019/mixed14_c01_exclusion_20261019.md), this excludes every mixed one-lower/four-upper maximum-span critical-root morphology in a seven-thin balanced-pair-free poset, and the order-dual four-lower/one-upper morphology. The combined branch has 92,318 degree-valid cores and 536,923 residual orders with exact generator certificates.

Other two-sided multiplicities remain separate. This is not a full seven-thin balance theorem.

## 8. Attribution and priority

The classical setting is G. Brightwell and C. D. Wright, *The 1/3–2/3 Conjecture for 5-Thin Posets* (1992), https://doi.org/10.1137/0405037; see also E. J. Olson and B. E. Sagan, *On the 1/3–2/3 Conjecture*, https://arxiv.org/abs/1706.04985. General majority-order, exchange, counting and nonnegative-certificate methods are not claimed as original. These sources provide context rather than an attribution of this exact finite receipt. All necessary lemmas are proved above. No comprehensive literature comparison or novelty claim is made.
