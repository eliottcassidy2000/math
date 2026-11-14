# The first affordable complete phase companions

Session 328. **Proved; independent review passed.** The minimum is over all physical tournament sources, not only the four-switch palette from S327. It concerns two entire first 4-divisibility phase envelopes being complementary, and their exact target repairs. Target-specific partial repairs can occur at much smaller clocks; this minimum does not bound them. It is not full odd-spectrum coverage.

## 1. Which two first-phase envelopes can be complementary?

Fix an odd slope a≥3. For a positive odd H coprime to a, put

    E_H={cH mod4a : c=1,3,…,2a−1}.

This is the first necessary phase envelope for N=aQ+cH when H>1: every such physical source has 4|Q. For H=1 the same envelope applies to sources of order at least two, but the one-vertex Q=2 source is an exception; no fixed-H impossibility assertion below uses H=1. It has exactly half of the odd residue classes modulo4a.

**Theorem 1 (complete companion classification).** Two unit H values have complementary envelopes iff

    H₂≡−H₁ mod4a,  or  H₂≡H₁+2a mod4a.                 (1)

**Proof.** Let S be the odd half-interval1,3,…,2a−1 modulo4a. Its multiplicative stabilizer among the units is exactly{1,2a−1}. Indeed kS=S implies k∈S. Replacing k by2a−k preserves this property, since for odd x one has (2a−k)x≡2a−kx, and2a−S=S. Thus assume1≤k≤a. If k≥3, the open interval(2a/k,4a/k) contains an odd integer c: its length is greater than two when k<a; when k=a it is(2,4). Also c<2a. Then kc lies in the opposite half, a contradiction. Thus k=1; undoing the replacement gives k=1 or2a−1. Both work. Since −S is the complementary half, the multipliers interchanging halves are precisely−1 and2a+1. Multiplying by H₁ gives (1), using oddness of H₁.

This proves that the two S327 companion mechanisms exhaust all complete first-phase companions. It does not say their actual finite target sets at one Q cover all residues or all integers.

## 2. A bounded-clock census without a strong-source assumption

For every tournament T, vertex v and positive integer z,

    Q_T(z) ≥ z Q_(T−v)(z).                             (2)

For each old ordered z-cover, select one of its slots and insert v somewhere into that directed path, choosing a deterministic legal position. Such a position always exists: prepend or append if possible, otherwise cross the first change from vertices beating v to vertices beaten by v. Deleting v recovers the old cover and its selected slot, so this is an injection. Empty paths cause no exception.

Consequently Q_T(2)≥2^|T|. Moreover a nontransitive tournament contains a triangle. Starting from its Q2=12 and applying (2) to all other vertices gives

    Q_T(2)≥3·2^(|T|−1)  for nontransitive T.             (3)

These facts give a complete budget recursion. To enumerate order-n sources with Q≤B, take every order-(n−1) source with Q≤floor(B/2), add a vertex with every possible cut, retain the children satisfying the exact Q budget, and merge isomorphic graphs. If B<2^n there are none; if B<3·2^(n−1), the transitive tournament is the sole possible class. Canonicalizing by the sorted outdegree classes and then all permutations within equal-degree classes is complete: an isomorphism must preserve each degree class. No SCC hypothesis or finite-bank assumption is used.

At B=4096 the order cutoff is twelve. The executable census gives:

| order | isomorphism classes with Q≤4096 |
|---:|---:|
|1,2,3,4,5,6,7|1,1,2,4,12,56,456|
|8,9,10,11,12|2160,716,151,18,1|

The3578 classes yield919 different(H,Q) pairs. The C++ verifier has two independently implemented endpoint recurrences (outgoing updates and incoming summation); optimized outgoing and unoptimized incoming runs regenerate identical catalogues and receipts. Explicit runtime checks enforce the stated counts and companion list.

**Theorem 2 (global minimum clock).** Among all pairs of nonempty tournament sources with a common exact Q and H values coprime to99, the least Q admitting complementary first-phase envelopes is4096. At this first clock the only companion pairs are

    (H₁,H₂)=(19,217), (25,223).                         (4)

The complete catalogue has no sign or additive companions below4096. At4096 its entire H fibre is

    {1,5,19,25,217,223},

and (1) selects exactly (4). Thus the minimum is all-order: the size bound and recursion exhaust every possible source below this clock. It is not extrapolation from a sampled list of strong order-eight graphs.

This is much smaller than303104, the S327 minimum in its four-switch palette. The earlier scoped minimum is unchanged; adding different actual source cores improves the global budget.

## 3. Literal optimum sources

Here are explicit sources at the optimum clock. Ordinal joins mean every earlier block beats every later block.

- H19: the strong core(2,20,9,3,45,15), profile(19,256), followed by four singleton SCCs.
- H25: two copies of the strong core(2,4,9,3), each of profile(5,32), followed by two singleton SCCs.
- H217: the strong order-eight core(2,20,41,19,101,203,143,31), profile(217,4096).
- H223: the strong order-eight core(2,36,49,7,75,217,143,31), profile(223,4096).

The first two sources have order ten and the last two order eight. Each pair differs in H by exactly198, so their first-phase envelopes are complementary, while their actual clock is literally equal to4096.

## 4. Exact repairs, not just a common residue orbit

Let a=99 and suppose physical equal-clock sources have counts h and h+198d, where d is odd and both counts exceed one and are coprime to99. Let c be a published odd digit with phase r_c. For any k≥r_c, pad either base source by k−r_c singletons and apply its c-template. Both strong target families are physically realized:

    N_−=99q·2^k+c h,
    N_+=99q·2^k+c(h+198d).                              (5)

At target N_+, the other H value h has the **same forced digit c modulo198**, but its required effective clock is

    Q_req=q·2^k+2cd ≡2 mod4.                            (6)

At N_− the other H value requires q·2^k−2cd, again2 modulo4 (or is nonpositive). Thus neither target can be obtained by the **entire fixed-other-H fibre** of the99-digit formula, at any source order or nonnegative phase. The first source is not merely unavailable in a selected bank: its forced clock is physically impossible. Nevertheless the correct companion gives a literal strong realization at the paid clock q·2^k. This is an absolute target repair.

For the optimum pairs, q=4096 and d=1. At c=1,k=0 the four directly counted strong outputs are

| source H | exact target | strong output order |
|---:|---:|---:|
|19|405523|17|
|217|405721|15|
|25|405529|17|
|223|405727|15|

The target405529 is prime. Its H223 candidate is forced to clock4094 and fails. The phase-zero digit5 at source H25 produces

    405629=7·57947,

with prime cofactor57947; the H223 candidate clock is4086 and fails. Digit7 at H217 produces the prime407023, while the H19 candidate clock4110 fails. The replay directly counts all these constructed strong cores, including the order-nineteen405629 witness, and checks the finite primality claims.

The exact paid boundary is useful near the old30000 examples: no complete equal-clock first-phase companion can enter any positive99 chart below405523. At q=4096 the lowest companion source is H19 and c≥1; a larger q is at least4100 because physical clocks above two are divisible by four, making99q+1 already exceed405523. Equality is attained in the table. Thus this specific complete-companion strategy is genuinely too expensive to repair a target near30000. Such a repair needs a different clock, a partial phase alternative, or another cut.

## 5. What does not overlap, and what remains open

The new targets are repairs by changing the source count, not duplicate representations of the same target by these two sources at their common exact clock. For either optimum pair the two H values are coprime and the larger exceeds197, so equality cH₁=dH₂ with allowed digits is impossible. Independent dyadic retiming cannot fix this: if the exponents differ, the leading-term difference has magnitude at least99·4096=405504, whereas the absolute offset difference is less than197·223=43931. Hence the two optimum source target families are completely disjoint, even allowing independent padding and all99 digits. Their arithmetic phase envelopes complement each other, but their absolute outputs remain distinct.

The improvement is a sharp, globally certified first affordable phase companion and explicit strong repairs in both atom lanes. It is not an all-target algorithm. Higher phase levels, variable source grades and arbitrary exact budgets remain separate obligations. Chan–Pak remains open.

## Reproduce

Compile `bounded_catalog.cpp` with `-O3` for outgoing DP, and with `-O0 -DINCOMING` for the separate incoming DP. Both produce the same4096 catalogue. `replay_switch.py` checks all919 stored literal profiles, the companion classification on the finite unit group, the optimum source masks, the direct output counts, and the target-fibre rejections. Its normal and optimized Python receipts match. The finite negative theorem depends on the complete catalogue regeneration and its proved recursion, not merely checking the stored positive witnesses.
