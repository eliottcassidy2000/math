# Mixed one-lower/three-upper maximum-span exclusion with three common neighbors

Date: 2026-10-19. Status: proved exact morphology theorem with arbitrary actual exterior. Every residual certificate was independently replayed, and the deterministic coverage union passed. No full seven-thin balance theorem or novelty claim.

## 1. Theorem and complete core catalogue

Suppose a finite poset P of incomparability degree at most seven has no balanced pair. Let x precede y be a maximum-index-span incomparable pair in its strict >2/3 majority order. The following critical-pair morphology is impossible when c=3:

    D(x)⊆D(y), |D(y)\D(x)|=1,
    U(y)⊆U(x), |U(x)\U(y)|=3,
    c=|I(x)∩I(y)|.

The order-dual statement holds. Label the roots x=0,y=1, the lower vertex b=2, common vertices 3,...,c+2 and upper vertices c+3,...,c+5. The mandatory relations are b<y and x<a for each upper a. Apart from these, all core relations are among nonroots, directed forward in the block order lower/common/upper. Internal common and upper orders are arbitrary. Any reverse relation across the blocks would compare a common vertex with a root or compare x and y. Conversely every naturally labeled order on the 4+c nonroot vertices gives a valid core. Natural labeling within each block loses no isomorphism type.

For c=3 the nonroot poset has seven naturally labeled vertices. There are 96,428 such cores before the degree screen. The only over-budget vertices can be common vertices with no comparable core neighbor: they have degree eight. If one specified common vertex is isolated, the remaining nonroot catalogue has 4,824 orders; two isolated common vertices give 357, and three give 40. Inclusion-exclusion gives

    3·4,824−3·357+40=13,441

over-budget cores. Exactly 82,987 cores satisfy all-vertex degree≤7.

The scan generates naturally labeled nonroot orders by adjoining each new maximal label with any ideal as its predecessor set. The independent verifier instead enumerates common and upper orders, lower-to-common upsets, common-to-upper ideals and lower-to-upper upsets, enforcing transitivity. It checks equality of the complete resulting core sets and all degree counts.

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

For completeness, an additional analytic core screen records the following necessary directed arcs: all comparabilities, x→y, x→z→y for common z, and v→u whenever H(u,v)≤1. Under no balance these all belong to the strict majority order. Any cycle excludes the core. Here the analytic balanced-pair screen excludes 31,090 cores, and the cycle screen excludes an additional 8,485. The remaining 43,412 cores require the projected-order and boundary-generator analysis below. For an upper vertex a, H(a,y)=1_{b not<a}+|U_Q(a)|, because A<a and all B lie above y. Thus every maximal upper vertex is forced after y (H(a,y)≤1); an upper vertex with one upper successor is also forced after y when b<a. No such conclusion is asserted from this screen when b∥a and a has one upper successor.

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

## 5. Complete finite coverage and certificate totals

Of the 82,987 degree-valid cores:

- 31,090 are excluded by a two-direction defect bound proving a balanced pair;
- 8,485 additional cores are excluded by a forced-majority cycle;
- the remaining 43,412 cores are exhaustively handled by projected-order checks.

Within those residual cores,13,325,904 necessary projected orders are excluded directly by a forward directed-defect bound. Exactly 401,867 further orders require boundary-generator certificates. The positive certificates cover 303,794 orders by a single pair,92,811 by two rows and 5,262 by integer certificates reconstructed from numerical LP dual suggestions. There are 38,436 distinct certificates.

Every reconstructed suggestion is converted to positive integer weights and every generator margin is then checked exactly. The largest weight is 650,310,314,725,595,671,220,042,640. This is an ordinary Python arbitrary-precision integer, not a machine integer or approximate number. Search acceptance explicitly converts generator entries to Python integers before multiplication. The independent replay uses only standard-library integer arithmetic. No floating-point inequality establishes any accepted certificate.

The previous sections prove that each successful finite certificate contradicts the actual strict-majority order for every permitted exterior. This completes the c=3 theorem. The degree-invalid cores cannot occur in the stated seven-thin setting.

## 6. Independent replay, deterministic sharding, and reproduction

The catalogue and all 39,575 analytic exclusions were independently regenerated and verified. The first 2,031 residual cores were scanned serially and separately replayed. The remaining 41,381 core IDs were deterministically partitioned into four disjoint shards of 10,346,10,345,10,345,10,345 cores. Each shard invoked the independent standard-library replay on each entire core record before committing that record. All four shards finished with `ALL_SHARD_CERTIFICATES_REPLAYED`; the prefix receipt verifies all 2,031 prefix cores.

The final union audit checks that the prefix is exactly the initial residual-ID segment, that the four planned shards partition exactly the remaining IDs, that every record has its prescribed predecessor masks, that all source and record hashes match, and that every shard has terminal independently-replayed status. It recomputes all result counts. Thus no order or core is silently omitted through sharding.

Files and commands:

- `mixed13_c3_exclusion_20261019.md`: this self-contained proof.
- `analytic_cover.json.gz`: complete core catalogue and analytic witnesses.
- `residual_certificates.jsonl.gz`: serial prefix records.
- `shard_0` through `shard_3`: exact records and terminal replay-bearing receipts.
- `independent_replay.json`: independent catalogue, analytic-cover, and serial-prefix verification receipt.
- `union_receipt.json`: deterministic union/provenance audit.
- `scan_c3.py`, `cone_engine.py`, `run_shard.py`: scan implementation.
- `exact_replay_engine.py`, `replay_c3.py`, `replay_full.py`: independent standard-library implementation.

Run `python -O audit_union.py` for the quick coverage/provenance audit. This is not a fresh mathematical replay. Run `python -O replay_full.py` to independently regenerate the catalogue and repeat every mathematical certificate check; this longer command requires only Python's standard library. The recorded proof already includes those independent checks across the serial prefix and the four shards. Explicit guards remain active under -O. The manifest authenticates the complete release.

## 7. Full mixed(1,3) corollary and limits

For this root morphology, y has incomparability degree 4+c. Thus c≤3 in a seven-thin poset. Combining this theorem with the [c≤2 theorem](../mixed13_c012_exclusion_20261019/mixed13_c012_exclusion_20261019.md) excludes every mixed one-lower/three-upper maximum-span critical-root configuration. The order-dual three-lower/one-upper configuration is excluded as well. Across c=0,1,2,3 there are 88,208 degree-valid cores and 417,782 residual orders with exact generator certificates.

Other two-sided multiplicities remain separate. No conclusion that every seven-thin poset has a balanced pair is claimed.

## 8. Attribution and priority

The classical setting is G. Brightwell and C. D. Wright, *The 1/3–2/3 Conjecture for 5-Thin Posets* (1992), https://doi.org/10.1137/0405037; see also E. J. Olson and B. E. Sagan, *On the 1/3–2/3 Conjecture*, https://arxiv.org/abs/1706.04985. Majority-order, exchange-injection, linear-extension counting, and positive-certificate methods are not claimed as original. These sources provide context, not an attribution of this exact finite receipt to those authors. All lemmas needed here are proved above. A comprehensive literature comparison has not been performed, and no novelty or priority claim is made.
