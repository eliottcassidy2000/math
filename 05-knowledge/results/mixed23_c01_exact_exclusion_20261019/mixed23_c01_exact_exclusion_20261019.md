# Mixed two-lower/three-upper maximum-span exclusion for at most one common neighbor

Date: 2026-10-19. Status: finite morphology theorem relative to the known six-thin theorem with arbitrary actual exterior. This excludes c=0,1 only, with the explicitly stated known-six-thin theorem dependency; c=2,3 remain separate. No full seven-thin balance theorem or novelty claim.

## 1. Theorem and complete core catalogue

Suppose a finite poset P of incomparability degree at most seven has no balanced pair. Let x precede y be a maximum-index-span incomparable pair in its strict >2/3 majority order. The following critical-pair morphology is impossible when c≤1:

    D(x)⊆D(y), |D(y)\D(x)|=2,
    U(y)⊆U(x), |U(x)\U(y)|=3,
    c=|I(x)∩I(y)|.

The order-dual statement holds. Label the roots x=0,y=1, the lower vertices b1=2,b2=3, common vertices 4,...,c+3 and upper vertices c+4,...,c+6. The mandatory relations are b1<y and b2<y and x<a for each upper a. Apart from these, all core relations are among nonroots, directed forward in the block order lower/common/upper. Internal lower, common and upper orders are arbitrary. Any reverse relation across the blocks would compare a common vertex with a root or compare x and y. Conversely every naturally labeled order on the 5+c nonroot vertices gives a valid core. Natural labeling within each block loses no isomorphism type.

There are exactly 357 and 4,824 distinct labeled cores for c=0,1. Every core has at most eight vertices, hence automatically satisfies the core incomparability budget. The independent verifier generates all subsets of allowed nonroot edges, closes them transitively, deduplicates, and compares the resulting entire set with the certificate catalogue.

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

For completeness, an additional analytic core screen records the following necessary directed arcs: all comparabilities, x→y, x→z→y for common z, and v→u whenever H(u,v)≤1. Under no balance these all belong to the strict majority order. Any cycle excludes the core. In this c≤1 receipt, no additional cores require the cycle screen; its inclusion does not affect the exact counts below. For an upper vertex a, H(a,y)=|{b1,b2}\D_Q(a)|+|U_Q(a)|. Only orientations satisfying the actual bound H≤1 are excluded by this screen; a maximal upper vertex need not satisfy it.

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

## 5. Exact finite result

| c | Cores | Analytically balanced cores | Known-six-thin short-span orders | Directed-screen orders | Residual orders | Distinct certificates |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 357 | 0 | 1,452 | 20,870 | 865 | 250 |
| 1 | 4,824 | 142 | 4,442 | 810,230 | 22,332 | 4,445 |

All 23,197 remaining orders have positive integer boundary-generator certificates: 18,917 single-row orders, 3,974 two-row orders and 306 orders with integer certificates reconstructed from LP-dual suggestions. Each suggestion is converted to positive integers and checked on every generator exactly. The largest coefficient is 6,861,073,367,827,309,620; products and sums use Python arbitrary-precision integers, not fixed-width arrays.

### Why the known-theorem filter is substantive

Without the filter, the core with predecessor masks (0,12,0,4,1,29,5) and projected order (2,0,3,4,1,6,5) has a genuine relaxed boundary-cone obstruction. Its 73 core extensions support an exact positive mixture of four generators in which every forward incomparable comparison exceeds 2/3. The full rational mixture and an independent brute-force replay are preserved in `unfiltered_obstacle`.

The root span is three, so the known six-thin theorem excludes it globally. This example shows a limitation of the unfiltered core-only cone; it is not an actual counterexample poset or an asserted realizable exterior.

## 6. Independent verification and files

Run `python -O verify_mixed23.py`. The standard-library verifier independently regenerates every core, recomputes all forced-comparability and defect bounds, enumerates the necessary majority orders, reconstructs all retained boundary generators and verifies every positive weighted inequality exactly. It checks the complete order coverage and assignment digests. Explicit guards remain enabled under -O. `certificates.json.gz` contains all records, `summary.json` the search receipt and `replay.json` the independent receipt. `scan_mixed23.py` is the reproducible search; only suggestions come from numerical LP.

The arbitrary-exterior factorization and nonnegative boundary expansion above turn the complete finite certificate into the stated theorem. This is not an enumeration of small ambient posets, and no exterior-size cutoff is imposed.

## 7. Attribution and scope

The classical context is G. Brightwell and C. D. Wright, *The 1/3–2/3 Conjecture for 5-Thin Posets* (1992), https://doi.org/10.1137/0405037; see also E. J. Olson and B. E. Sagan, *On the 1/3–2/3 Conjecture*, https://arxiv.org/abs/1706.04985. General majority-order, exchange, extension-counting and positive-certificate methods are not claimed as original. These references provide context, not a claim that they contain this exact finite receipt. Apart from the explicitly cited known six-thin theorem, all needed lemmas are proved here, and no novelty or literature-priority claim is made. The c=2,3 cases are outside this theorem and this receipt.
