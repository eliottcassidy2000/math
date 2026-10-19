# Mixed one-lower/six-upper maximum-span exclusion with no common neighbors

Date: 2026-10-19. Status: proved exact morphology theorem with arbitrary actual exterior. All residual certificates were independently replayed and the complete union audit passed. No full seven-thin theorem or novelty claim.

## 1. Theorem and complete core catalogue

Suppose a finite poset P of incomparability degree at most seven has no balanced pair. Let x precede y be a maximum-index-span incomparable pair in its strict >2/3 majority order. The following critical-pair morphology is impossible when c=0:

    D(x)⊆D(y), |D(y)\D(x)|=1,
    U(y)⊆U(x), |U(x)\U(y)|=6,
    c=|I(x)∩I(y)|.

The order-dual statement holds. Label the roots x=0,y=1, the lower vertex b=2, common vertices 3,...,c+2 and upper vertices c+3,...,c+8. The mandatory relations are b<y and x<a for each upper a. Apart from these, all core relations are among nonroots, directed forward in the block order lower/common/upper. Internal common and upper orders are arbitrary. Any reverse relation across the blocks would compare a common vertex with a root or compare x and y. Conversely every naturally labeled order on the 7+c nonroot vertices gives a valid core. Natural labeling within each block loses no isomorphism type.

There are 96,428 naturally labeled orders on the seven nonroots. Every core has nine vertices, but every nonroot has a mandatory comparable root, so its core incomparability degree is at most seven. The root degrees are two for x and seven for y. Thus all 96,428 cores are degree-valid.

The scan generates naturally labeled nonroot orders by successively adjoining a maximal label with any predecessor ideal. The independent verifier instead chooses an arbitrary naturally labeled order on the six upper vertices and any upper upset as the vertices above the lower vertex b. It verifies equality of the complete core sets and all degree bounds.

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

For completeness, an additional analytic core screen records the following necessary directed arcs: all comparabilities, x→y, x→z→y for common z, and v→u whenever H(u,v)≤1. Under no balance these all belong to the strict majority order. Any cycle excludes the core. The balanced-pair screen excludes 32,026 cores, and there are no additional cycle-only exclusions. The remaining 64,402 cores require the projected-order and boundary-generator procedure below. For an upper vertex a, H(a,y)=1_{b not<a}+|U_Q(a)|, because A<a and all B lie above y. Thus every maximal upper vertex is forced after y (H(a,y)≤1); an upper vertex with one upper successor is also forced after y when b<a. No such conclusion is asserted from this screen when b∥a and a has one upper successor.

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

All 96,428 cores are covered:32,026 balanced-defect exclusions and 64,402 residual cores. Within the residual cores,26,878,012 necessary projected orders are excluded by directed-defect bounds. The remaining 696,230 orders have positive integer boundary certificates:573,566 single-row orders,107,920 two-row orders and 14,744 orders using integer certificates reconstructed from LP-dual suggestions. There are 86,529 distinct certificates. The largest integer weight is 70,952,808,454,000,004,685,213.

Every numerical suggestion was converted into positive integers and checked against every retained generator with arbitrary-precision integer arithmetic. Search acceptance explicitly converts generator entries to Python integers before multiplication. Independent replay uses only standard-library integer operations. No floating-point inequality establishes any accepted certificate.

## 6. Independent verification and coverage

The complete 96,428-core catalogue and all 32,026 analytic exclusions were independently regenerated and verified. The 64,402 residual IDs were partitioned deterministically into four disjoint shards of 16,101,16,101,16,100,16,100. Each shard invoked the independently implemented standard-library replay on every complete core record before committing it. All four finished with terminal independently-replayed success.

The final union audit checks the exact deterministic partition, core IDs and predecessor masks, terminal receipts, source and record hashes, and numerical totals. It verifies that nothing was silently skipped or duplicated. The serial prefix is intentionally empty; its pause/partial labels apply only to that empty prefix. `union_receipt.json` is the combined terminal receipt.

Run `python -O audit_union.py` for quick coverage/provenance verification. Run `python -O replay_full.py` to independently regenerate the catalogue and repeat every mathematical certificate check using only Python's standard library; the latter command is slower. The recorded independent mathematical checks already cover all analytic exclusions and all shard records. Explicit guards remain active under -O.

The release includes the full catalogue and analytic witnesses, exact certificate streams, replay-bearing receipts, independent verifier, reproducible search and hash manifest. No exterior-size cutoff is used.

## 7. Scope and combined reduction

For this morphology y has incomparability degree 7+c, so seven-thinness forces c=0. This theorem therefore excludes the entire one-lower/six-upper branch and its order-dual six-lower/one-upper branch.

Combined with the earlier zero- and one-exclusive-side results, the [two-sided two-exclusive reduction](two_sided_two_exclusive_reduction_20261019.md) shows that a maximum-span critical pair in a hypothetical seven-thin counterexample must have at least two lower exclusive neighbors and at least two upper exclusive neighbors. The corollary lists all dependencies and exact degree budgets.

The remaining genuinely two-sided morphologies are not settled here. No full seven-thin balance theorem is claimed.

## 8. Attribution and priority

The classical setting is G. Brightwell and P. Wright, *The 1/3–2/3 Conjecture for 5-Thin Posets* (1992), https://doi.org/10.1137/0405037; see also E. J. Olson and B. E. Sagan, *On the 1/3–2/3 Conjecture*, https://arxiv.org/abs/1706.04985. General majority-order, exchange-injection, extension-counting and nonnegative-certificate methods are not claimed as original. These references provide context rather than an attribution of this exact finite receipt. All necessary lemmas are proved above. No comprehensive literature comparison or novelty claim is made.
