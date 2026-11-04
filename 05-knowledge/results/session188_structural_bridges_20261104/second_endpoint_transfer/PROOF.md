# A second endpoint test for arbitrary nontrivial child substitutions

## 1. A quotient invariant beyond endpoint parity

Let Q be any tournament on m fixed quotient names. Write S_i(Q),E_i(Q) for its Hamilton-path start and end counts and p_i(Q)=S_i+E_i modulo two.

For quotient names i,j, let R_ij(Q) count directed run words whose multiset contains every quotient name once and one additional copy of j, with first name i. Equal consecutive names are forbidden. The two copies of j are indistinguishable; there is no factor of two for labelling them. Every consecutive pair must be an edge of Q.

Define

    sigma_i(Q) = S_i(Q)-E_i(Q)
                  +2 sum_j [R_ij(Q)-R_ij(Q^op)] modulo4.             (1)

In particular sigma_i modulo two is p_i. The invariant can be computed entirely on words of lengths m and m+1, without constructing any inflated tournament.

**Second endpoint transfer.** Suppose p_i(Q)=0. For arbitrary nontrivial tournaments A_j, on fixed disjoint physical name sets, compare T=Q[A_j] and U=Q^op[A_j]. Profile-preserving target child replacements are also allowed. Use ports x->y and the physical block query z_i=(a,c)=(e_i,0): A_i points to x, x points to every other old name, and y points to every old name. Then

    Delta H(z_i)=2 sigma_i(Q) modulo8,
    Delta q_h(z_i)=0 modulo8 for every h>=2.                         (2)

Here Delta means source minus target and q_h counts ordered spanning covers by h nonempty paths. Therefore sigma_i=2 is a single-query separator, with contrast4 modulo8, through every fixed ambient root substitution. If sigma_i=0, the criterion is silent; it is not a response-collision assertion.

Every child must have at least two vertices. No quotient parity or total physical-order assumption is imposed.

## 2. Aggregate starts modulo four

Write S_(A_i)(T) for the sum of Hamilton-start counts over the names of A_i. A maximal-child-run decomposition gives this as a multilinear functional in child cover profiles. Let L=product_j H(A_j), which is odd. Modulo four, normalized child profiles have

    q1/H=1, q2/H=2, q3/H=0 or2, qk/H=0 for k>=4.                    (3)

The q2 assertion uses q2=2 modulo4 and H odd; multiplication of2 by an odd inverse leaves2. The q3 assertion needs only its evenness. The tail follows from k! divisibility. These facts are included in the independently established universal cover residues.

Any term with two repeated children vanishes modulo four. The possible one-child triple-return terms also vanish. Indeed a skeleton using j three times has the form

    prefix, j, B, j, C, j, suffix,

where B and C are nonempty words on disjoint sets of the other quotient names. Interchange B and C. All directed boundary conditions are unchanged, as are the first and last names of the whole word. The exchange is fixed-point-free because B and C have disjoint nonempty supports. Thus every aggregate-start triple-run coefficient is even, and its even q3 factor makes it zero modulo four. This argument concerns skeletons; it does not multiply independently chosen overlapping paths.

Only the all-one-run words and the single two-run words remain. Consequently

    S_(A_i)(T) = L [S_i(Q)+2 sum_j R_ij(Q)] modulo4.                  (4)

Apply the same equation to U. Every corresponding complete profile, and hence H, is preserved. The start-skew is therefore

    S_(A_i)(T)-S_(A_i)(U) = L sigma_i(Q) modulo4.                    (5)

When p_i=0, sigma_i is even. Multiplication by the odd L then fixes sigma_i modulo four. Notice that this simplification is valid only for the skew in the p_i=0 case; it is not a general claim that all normalized endpoint values equal their literal pair values.

For a transitive-pair child at each quotient name, q1=1,q2=2,q3=0. The same reasoning shows that sigma_i is exactly the aggregate start-skew of that pair inflation modulo four. This is a useful alternative computation, not a new hypothesis.

## 3. Recovering the physical anchor

The exact uniform-block start identity from `uniform_block_starts/PROOF.md`, Section1, states

    Delta(e_i,0)-2 Delta(0,e_i)
        =2[S_(A_i)(T)-S_(A_i)(U)].                                 (6)

Its operations are on complete physical responses. The complete coarse modulo-four law in `parity_profile_transfer/PROOF.md`, Section2, gives

    Delta(0,e_i)=2p_i(Q) modulo4.                                  (7)

It applies because all children are nontrivial. Under p_i=0, twice (7) vanishes modulo eight. Equations (5) and (6) therefore give Delta(e_i,0)=2sigma_i modulo eight. In particular the anchor equals the all-pair anchor modulo eight, independent of all nontrivial child sizes, orientations and unequal profiles.

This uses only exact physical-query identities and a complete residue reduction. No cubic or quartic term of the response is discarded.

## 4. The complete augmented profile and ambient lift

Let N be the common augmented-core order. For N even, the inherited modulo-eight cover formulas give Delta q2=2Delta H and Delta q3=0. Since Delta H is0 or4, both higher contrasts vanish modulo eight.

For N odd, let z(B) be the number of directed Hamilton cycles of the augmented tournament B, counted up to rotation. The same formulas give

    Delta q2=2Delta H-4Delta z,
    Delta q3=2Delta H-4Delta z modulo8.                            (8)

For the query z_i, the only old-or-port vertex entering y is x. Every directed Hamilton cycle must therefore contain x->y. Delete x and y from such a cycle. The result is an old Hamilton path whose last name belongs to A_i, because those and only those old names point to x. Conversely every old Hamilton path ending in A_i closes uniquely to such a cycle by appending x,y. This is a bijection on cycles up to rotation, with no extra multiplicity.

The old aggregate end-skew is p_i modulo two: only one run per child survives modulo two, and the product of child Hamilton counts is odd. Hence Delta z=p_i=0 modulo two. Equation (8), together with Delta H divisible by four, proves Delta q2=Delta q3=0 modulo eight. For h>=4, the factor h! makes each q_h separately zero modulo eight. This proves all of (2), including odd total core order.

In any fixed ambient, give both ports the original core's exterior incidences. The augmented core is a module and its ambient Hamilton contribution is sum_h c_h q_h, with nonnegative integer coefficients and c1 equal to the odd Hamilton count of the contracted ambient tournament. Therefore (2) lifts as

    Delta H_ambient = c1*2sigma_i =2sigma_i modulo8.                 (9)

The ambient statement rests on the full profile contrast. It does not assert a generic lift of arbitrary nonzero odd-core Hamilton residues.

## 5. Two useful prime quotients and a genuine boundary

For the prime-five quotient

    Q5=(20,1,2,7,14),

exact word counts give

    S=(1,1,1,3,3), E=(1,3,3,1,1), p=0,
    sigma=(0,2,2,2,2).

Thus every anchor i=1,2,3,4 supplies an all-order arbitrary-nontrivial-child separator with arbitrary fixed ambient lift. This is the alternating W5 class; the theorem does not rely on identifying it by name. Primality is separately checked on all proper subsets.

For the endpoint-balanced prime-eight quotient

    Q8=(238,4,240,118,195,146,34,74),

sigma=(2,0,0,0,0,2,2,2), so anchors0,5,6,7 work. The separate `prime_eight_uniform_anchor` packet independently proves the anchor0 full-profile identity by a complete three-coordinate, 129-run-vector certificate. The present argument explains why its arbitrary-child independence occurs.

The directed triangle and the regular prime-five quotient `(24,17,3,6,12)` have p=sigma=0 at every name. They remain blind to this sufficient test, though other two-port probes separate them. This is an explicit boundary, not evidence of a collision.

The nontrivial-child condition matters: for the eight-name quotient with singleton children, query(1,0) gives771 and2403, a contrast divisible by eight, whereas its all-pair contrast is4 modulo eight. The theorem does not silently include singleton children.

## 6. Original-module families and highest changed nodes

For a source core with prime quotient Q and a sigma_i=2 anchor, the pair {z_i,Jz_i}, where J(a,c)=(1-c,1-a), detects every original-module reversal family containing the core root oddly, provided all children are nontrivial. The two exact baseline contrasts are opposite under whole-augmented-profile duality, also in any fixed ambient. Choose the member with negative source-minus-baseline contrast. Every target child profile is coordinatewise no smaller by hereditary original-module monotonicity, so nonnegative run substitution preserves that strict negative sign. The original selected modules must all lie inside this core; changes outside are not being authorized or silently folded into this statement.

As a separate necessary-condition consequence in the existing scalar-equal highest-node reduction, a potential full named-response collision cannot have a highest changed prime node with nontrivial children and a sigma_i=2 anchor. The node's profile-preserving quotient converse is separated by (9), and other highest-node changes can be restored because the lifted ports see those nodes uniformly. This uses the inherited antichain/restoration theorem; it is not a reconstruction proof or a claim that the first two endpoint tests are complete.

## 7. Reproduction and scope of computation

The all-order proof is Sections1–4. The verifier retains complete quotient data through five names, checks the repeated-name recurrence against literal multiset words, checks the return-segment involution on all applicable words through six names, and uses literal physical cut-aware subset DP for complete q1,q2,q3 query controls. It also retains the exact five- and eight-name examples, mixed child types, odd/even physical orders and fixed ambient controls. Finite checks are not used to extrapolate the theorem.

Run `python3 verify.py > normal.json` and `python3 -O verify.py > optimized.json`, then compare. All checks remain active in optimized Python. The C++ physical counts are exact integers with an explicit maximum order of18; q_h for h<=3 is bounded by18!*binom(17,2)<2^64. No package is installed.

Frozen source allowlist: `PROOF.md`, `README.md`, `verify.py`, `physical.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`. The exploratory `explore.py` is excluded.
