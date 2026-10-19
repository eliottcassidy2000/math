# Mixed (2,4) and (3,3) c=0 maximum-span critical-root exclusions

Date: 2026-10-19. Status: complete exact certificates and independent replay. The known six-thin theorem is an explicit external dependency for projected root spans at most three. No novelty claim is made.

## 1. Scoped theorem and complete catalogue

Let P be a finite poset of maximum incomparability degree at most seven with no balanced pair. Let x precede y be a maximum-index-span incomparable pair in its strict >2/3 majority order, with critical-pair conditions D(x)⊆D(y) and U(y)⊆U(x). Put

    L=D(y)\D(x), T=U(x)\U(y), Z=I(x)∩I(y).

Neither morphology (|L|,|T|,|Z|)=(2,4,0) nor (3,3,0) is possible. Their order-dual morphologies are excluded as well. No other morphology is asserted in this theorem.

Label x=0,y=1, the l lower exclusive vertices 2,...,l+1, and the u upper exclusive vertices l+2,...,l+u+1. Each lower vertex is below y and each upper vertex is above x. Apart from these mandatory relations, all core relations lie among nonroots and may run forward in the block order lower/upper; a reverse relation would compare the two roots. Internal orders in each block are arbitrary. Every such order can be naturally labeled within each block.

There are exactly 4,824 naturally labeled posets on the six nonroot vertices, and each gives one labeled core in either morphology. All cores have eight vertices and thus automatically satisfy the core degree bound at most seven. Exterior degree constraints still matter and are retained below.

Search generates natural orders by adjoining each maximal label with an arbitrary predecessor ideal. The independent verifier instead enumerates all 2^15 choices of forward relations and takes transitive closure, deduplicating the resulting predecessor tuples. It checks equality of the complete 4,824-core set with the supplied records, not just agreement of total counts.

## 2. Exact exterior defect screens with forced-comparability closure

Every vertex outside Q={x,y}∪L∪T lies in A below both roots or B above both. We have A<B, and no core vertex can lie below A or above B. Any core vertex of degree seven is saturated and comparable with all exterior vertices, so A<s<B.

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

For completeness, an additional analytic core screen records the following necessary directed arcs: all comparabilities, x→y, x→z→y for common z, and v→u whenever H(u,v)≤1. Under no balance these all belong to the strict majority order. Any cycle excludes the core. The balanced-defect screen excludes 223 cores for (l,u)=(2,4) and 202 for (3,3). The forced-cycle screen gives no additional exclusions in these two complete families. For an upper vertex a, the correctly computed bound is H(a,y)=|L\D_Q(a)|+|U_Q(a)|. In particular, maximality of a alone does not make this bound at most one when multiple lower exclusive vertices remain incomparable with a.

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

## 5. Complete exact results

| (l,u,c) | Valid cores | Balanced-defect cores | Short-span orders, known six-thin theorem | Directed-defect orders | Residual certified orders |
|---|---:|---:|---:|---:|---:|
| (2,4,0) | 4,824 | 223 | 5,958 | 845,820 | 25,408 |
| (3,3,0) | 4,824 | 202 | 9,109 | 904,460 | 30,857 |

There are no additional forced-cycle core exclusions. In the first family, the certified orders use 19,520 single-row, 5,540 equal-weight two-row and 348 reconstructed integer certificates. In the second, those counts are 25,646, 4,372 and 839. The numbers of distinct stored certificates are 5,450 and 6,506. The largest positive integer coefficients are 257,808,213,151,125 and 23,001,979,073,007,915, respectively. All generator products and sums use arbitrary-precision integers.

Every completed record was independently replayed before commit. A further full standard-library replay independently regenerated both catalogues and verified all analytic exclusions and all 56,265 residual orders. No cone obstruction occurred.

## 6. A shared certificate pattern, with its exact scope

Many stored certificates use two equally weighted forward comparisons: an internal lower comparison together with x before a lower vertex, or dually an internal upper comparison together with an upper vertex before y. Abstractly, if two forward comparisons C1,C2 satisfy

    3N_g(C1)+3N_g(C2)≤4E_g

for every retained boundary generator, then after the nonnegative exterior lift,

    p(C1)+p(C2)≤4/3.

They cannot both have probability strictly above 2/3. This supplies a shared two-comparison mechanism across the two morphologies, while allowing the low-probability witness to move from a root contact to an internal exclusive-vertex comparison. The argument is valid exactly when the generator inequalities have been verified; the comparison types alone do not imply the inequality for an arbitrary core.

Among the first family's distinct stored certificates, 720 have lower-internal plus x-to-lower types and 970 have upper-internal plus upper-to-y types. The corresponding second-family counts are 882 and 735. These are counts of stored certificates, not disjoint coverage counts of orders, and no universal two-row sufficiency claim follows. Higher-support integer certificates are also used. Their existence does not prove any lower bound on the support of an alternative certificate. The full weighted generator inequalities, rather than a conjectural uniform shortcut, establish the theorem.

## 7. Replay, files and attribution

Run `python -O replay_small.py` from this folder for a fresh full standard-library replay. It regenerates the two complete catalogues, checks record identity and coverage, reconstructs all projected orders and generators, and verifies every weighted inequality exactly. The scan entrypoint is `full_small.py`; `analytic_census.py` generates the two catalogues and analytic screens. `cone_engine.py` and `exact_replay_engine.py` are import-only helpers. Only the search needs NumPy/SciPy for numerical suggestions; the proof replay does not.

`full_progress_240.json` and `full_progress_330.json` record the completed search-plus-per-record-replay results. `independent_small_replay.json` records the further complete independent run. The compressed certificate files carry all records. Explicit guards remain enabled under `-O`.

The classical background includes G. Brightwell and C. D. Wright, *The 1/3–2/3 Conjecture for 5-Thin Posets* (1992), https://doi.org/10.1137/0405037; and E. J. Olson and B. E. Sagan, *On the 1/3–2/3 Conjecture*, https://arxiv.org/abs/1706.04985. General majority-order, exchange, extension-counting and positive-certificate methods are not claimed as original. Apart from the explicitly cited Peczarski six-thin theorem, the required mathematical reductions are given above. These references provide context, not a claim that they contain this exact finite receipt. No novelty claim is made.

The four larger branches (2,4,1),(3,3,1),(2,5,0),(3,4,0) had only preliminary analytic counts and bounded probes outside this packet. They are not covered by these certificates or this theorem. No full seven-thin conclusion is asserted.
