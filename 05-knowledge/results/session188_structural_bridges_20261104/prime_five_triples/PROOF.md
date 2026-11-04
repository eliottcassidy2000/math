# Exact finite separation for prime-five quotients with triple children

**FINITE-EXACT:** every strong prime five-vertex quotient, inflated by arbitrary physically named three-vertex child tournaments, is separated from its quotient converse with all child interiors unchanged. This is a fixed fifteen-old-name classification, proved by an exact compressed certificate. It is not an all-order child-size or quotient theorem, and no collision is claimed outside the checked class.

## 1. Reduction to complete finite data

Let T=Q[A_0,...,A_4] and U=Q^op[A_0,...,A_4], where Q is strong and prime and every child has three names. Ports have fixed edge x->y and a_u=[u->x], c_u=[u->y]. Each query below is uniform on every child.

Every named tournament on three vertices is either transitive or a directed triangle. Its complete ordered nonempty path-cover profile is respectively

    (q_0,q_1,q_2,q_3)=(0,1,6,6) or (0,3,6,6).

The uniform port response depends on each child only through this profile: read a Hamilton path as maximal nonempty child runs separated by other child runs or ports. An ordered cover fills each child's runs in their order of appearance. This is an exact bijection. It includes all positions of x,y, including separated ports, so no cubic or quartic response contribution is discarded.

There are 2^5=32 possible profile patterns, covering all 8^5 choices of labelled child orientations. The checker independently classifies all eight three-vertex tournaments and verifies their profiles.

The checker enumerates all 2^10=1,024 named five-vertex tournaments. It tests every nontrivial proper subset for the module condition and tests strong connectivity. The 264 strong prime labelled quotients form exactly three classes under simultaneous relabelling, represented in outgoing-row notation by

    (20,1,2,7,14),
    (18,8,3,5,14),
    (24,17,3,6,12).

Canonicalization tries all 120 permutations. A quotient relabelling acts simultaneously on its children and the query masks; it never separately forgets names in the two compared responses. Thus it is enough to give a certificate for these three representatives and all 32 profile patterns.

## 2. The fixed physical query and retained certificates

For each listed representative, the same query works:

    a=all names of child0,  c=all names of child1.

With consecutive three-name blocks this is the literal mask pair (7,56). Every one of the 96 representative/profile comparisons has a nonzero integer difference. The minimum absolute gaps across the 32 patterns are respectively

    596632200, 1633888982, 24362496.

The sign need not be constant across patterns. No nonzero residue modulo four is asserted.

The receipt retains both exact query values for every profile pattern, not just a success flag. It also retains each response as a 32-coefficient multilinear polynomial in h_i=q_1(A_i), with q_2=q_3=6 fixed. This allows independent evaluation of all patterns and preserves the complete count at the chosen physical query. It is not advertised as the full two-port coefficient table at order15.

## 3. Exact run-word computation

The automaton alphabet is the five child colors and the two distinct ports. A child color may appear one, two or three times, representing that many nonempty child runs; a port appears exactly once. Consecutive old colors must follow the quotient edge and cannot be equal. Old/port transitions follow the actual query masks, and the only port-to-port transition is x->y.

A state records the five run counts in {0,1,2,3}, the two used-port bits, and the last symbol. Prefix multiplicities are exact integers. Every state with both ports used and every child count positive contributes the product of its five ordered-cover counts. The automaton then continues to account for longer words with more runs. No word exceeds 17 symbols, and every valid run word is counted once according to its successive states.

For a terminal count vector, each child visited once contributes h_i and each child visited two or three times contributes six. Hence its contribution is added to the corresponding multilinear h-monomial with integer weight 6 raised to the number of non-single-visit children. Evaluating h_i in {1,3} gives the exact actual query counts. Equality of the two complete named response functions would force equality at this particular query, so the retained nonzero differences certify separation throughout the stated finite class.

## 4. Independent physical-name controls

The separate C++ engine builds the actual seventeen-vertex augmented tournaments and counts Hamilton paths by a physical subset/last-vertex recurrence. It does not use the run-word or child-profile reduction. It checks all-transitive, all-cyclic and alternating child-type patterns for each representative, giving eighteen independent physical counts. They agree exactly with the compressed certificates.

All intermediate counts are bounded by 17!, which is below 2^64; the C++ unsigned 64-bit arithmetic is exact. Python run-word arithmetic uses unbounded integers. The checker rejects any control outside its declared order range. Explicit exceptions and nonzero process status retain every check in optimized mode.

## 5. Scope

This finite certificate reaches pair-free fifteen-name examples, including five directed-triangle children. It complements the unbounded marked-pair quotient results, since none of these children is a pair. It does not settle arbitrary child orders at least three, arbitrary larger prime quotients, arbitrary strong outer ambient lifting, or the full named-response reconstruction problem.

The proof consists of the exact run-word bijection, the complete quotient/profile classifications, and the retained finite integer certificates. The independent physical controls validate the implementation but are not substituted for the complete 96-pattern check.
