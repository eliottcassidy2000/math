# Mixed two-lower/three-upper maximum-span exclusion with two common neighbors

Date: 2026-10-19. Complete exact certificate and independent replay. Known-six-thin external theorem dependency remains explicit.

## 1. Theorem and complete core catalogue

Suppose a finite poset P of incomparability degree at most seven has no balanced pair. Let x precede y be a maximum-index-span incomparable pair in its strict >2/3 majority order. The following critical-pair morphology is impossible when c=2:

    D(x)⊆D(y), |D(y)\D(x)|=2,
    U(y)⊆U(x), |U(x)\U(y)|=3,
    c=|I(x)∩I(y)|.

The order-dual statement holds. Label the roots x=0,y=1, the lower vertices b1=2,b2=3, common vertices 4,...,c+3 and upper vertices c+4,...,c+6. The mandatory relations are b1<y and b2<y and x<a for each upper a. Apart from these, all core relations are among nonroots, directed forward in the block order lower/common/upper. Internal lower, common and upper orders are arbitrary. Any reverse relation across the blocks would compare a common vertex with a root or compare x and y. Conversely every naturally labeled order on the 5+c nonroot vertices gives a valid core. Natural labeling within each block loses no isomorphism type.

There are 96,428 naturally labeled orders on the seven nonroots. The only over-budget vertices can be common vertices with no comparable core neighbor, giving degree eight. Specifying one common vertex isolated leaves 4,824 cores; specifying both leaves 357. Thus 2·4,824−357=9,291 cores are invalid, leaving 87,137 valid cores.

The search adjoins maximal labels with predecessor ideals. The independent catalogue instead chooses the internal lower/common/upper orders, lower predecessor ideals for each common vertex, and lower/common predecessor ideals for each upper vertex, imposing all isotonicity and transitivity conditions. Equality of the complete core sets and degree screens is checked exactly.

## 2. Exact exterior defect screens with forced-comparability closure

Every vertex outside Q={x,y,b1,b2}∪Z∪T lies in A below both roots or B above both. We have A<B, and no core vertex can lie below A or above B. Any core vertex of degree seven is saturated and comparable with all exterior vertices, so A<s<B.

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

For completeness, an additional analytic core screen records the following necessary directed arcs: all comparabilities, x→y, x→z→y for common z, and v→u whenever H(u,v)≤1. Under no balance these all belong to the strict majority order. Any cycle excludes the core. The balanced-defect screen excludes 27,833 cores, and a forced-majority cycle excludes 1,618 further cores. The remaining 57,686 cores require the projected-order and boundary-generator procedure. For an upper vertex a, H(a,y)=|{b1,b2}\D_Q(a)|+|U_Q(a)|. Only orientations satisfying the actual bound H≤1 are excluded by this screen; a maximal upper vertex need not satisfy it.

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

## 5. Complete finite result

The 87,137 degree-valid cores are partitioned into 27,833 balanced-defect exclusions, 1,618 forced-cycle exclusions and 57,686 residual cores. The projected-order procedure excludes 8,450 orders using the known six-thin theorem and 25,542,178 orders by directed defects. All 476,837 remaining orders have exact positive integer certificates: 397,286 single-row, 76,315 two-row and 3,236 reconstructed from LP-dual suggestions. There are 54,408 distinct stored certificates.

The largest coefficient is 1,941,651,185,962,341,358,890. Every product and sum in the certificate verification uses arbitrary-precision Python integers; fixed-width numerical arrays are not trusted for the final inequalities. Every completed core record was checked by the independent standard-library verifier before the shard committed it. The union audit checks source hashes, exact record hashes, deterministic core IDs and complete disjoint coverage of the residual catalogue. No cone obstruction occurred.

## 6. Verification and reproducibility

Run `python -O audit_union.py` for the fast union/provenance audit. Run `python -O replay_full.py` for a fresh full mathematical replay using only the Python standard library. The latter regenerates the independent catalogue, checks every analytic exclusion and independently verifies every residual record, projected-order assignment and exact generator inequality. Explicit guards remain active under `-O`.

`union_receipt.json` is the terminal combined receipt. The root `progress.json` and `independent_replay.json` preserve an intentionally empty serial prefix; their pause/partial status labels refer only to that prefix, not to the completed four-shard union. Each shard's terminal progress file authenticates its records and source dependencies. The complete independent mathematical replay was already performed per record before commit; the fast union audit does not replace those mathematical checks.

The arbitrary-exterior lift proves the stated theorem, with no bound on exterior size. This is not merely an enumeration of small ambient posets.

## 7. Combined scoped corollary and attribution

Together with the [c=0,1 theorem](../mixed23_c01_exact_exclusion_20261019/mixed23_c01_exact_exclusion_20261019.md), this excludes mixed (2,3) maximum-span critical roots with c≤2 and their order duals. The c=3 morphology is outside this theorem and receipt. This is not a full seven-thin theorem.

The classical context is G. Brightwell and C. D. Wright, *The 1/3–2/3 Conjecture for 5-Thin Posets* (1992), https://doi.org/10.1137/0405037; see also E. J. Olson and B. E. Sagan, *On the 1/3–2/3 Conjecture*, https://arxiv.org/abs/1706.04985. General majority-order, exchange, extension-counting and positive-certificate methods are not claimed as original. These references supply context, not a claim that they contain this exact finite receipt. Apart from the explicitly cited Peczarski six-thin theorem, all necessary lemmas are given above. No novelty or literature-priority claim is made.
