# The 28-name original-module boundary has no full-response collision

**Producer finite-exact theorem and source-relative integration, pending independent review.** No unrestricted named-response injectivity, arbitrary-child prime-nine theorem, or arbitrary-ambient lift of the new scalar residues is asserted.

## 1. Statement and imported boundary

For a named tournament T, write F_T(a,c) for the Hamilton-path count after adjoining the two actual ports x->y, where a_v=[v->x] and c_v=[v->y]. Original-module operations mean simultaneous XOR of edges in selected modules of the original source, never modules created by earlier operations.

The reviewed incoming results, pinned in DEPENDENCIES.json, establish that any nontrivial full-response original-module collision has every highest changed original strong quotient prime of order at least nine, and every child of that node is nontransitive. The reviewed minimum27 packet also excludes every source with at most 27 old names, using its complete prime-nine census and the all-C3 finite boundary. We import the certified quotient census and endpoint partition; we do not re-enumerate or newly certify them here.

**New finite theorem.** Let Q be any prime tournament on nine quotient names. Give eight children the directed triangle and the remaining child any nontransitive four-vertex tournament. Then Q[children] and its quotient-only converse Q^op[children], with all physical interiors unchanged, have different full two-port responses. A fixed battery of 14 child-uniform queries covers every endpoint-zero quotient and every position and type of the four-child, after one common transport to the certified quotient coordinates. Together with the nine endpoint anchors, this is a 22-query battery.

**Source-relative consequence.** A nonidentity original-module XOR target with the same full named response requires at least 29 old names. This theorem is the finite result plus the reviewed structural inputs and the one-exterior-vertex argument below, not an unrestricted reconstruction claim.

## 2. The only new four-child profiles

Every nontransitive four-vertex tournament contains a directed triangle. The fourth vertex either dominates that triangle, loses to it, or has mixed incidences. The first two are the source/sink extensions of C3; every mixed case is the unique strongly connected four-vertex type C3[TT2,1,1]. Rotations of the triangle and direct checking of the six mixed incidence patterns give that classification.

The complete ordered spanning-path-cover profiles are

- C3: (3,6,6);
- C3 with a uniform source or sink: (3,18,36,24);
- C3[TT2,1,1]: (5,22,36,24).

The retained four_child_records.json exhaustively records all 64 labeled four-name tournaments, every nonsingleton original module, every simultaneous original-module XOR family, and every resulting complete profile. There are 40 nontransitive labeled inputs and 160 family outputs. Every such output has its source profile. This is also immediate from the classification: the only proper nonsingleton module is the cyclic triple in a source/sink extension, or the transitive pair in the strong type. Toggling that module and/or the whole four-child preserves its isomorphism type up to converse, hence its cover profile. The cyclic children likewise retain (3,6,6) under every original-module family.

This profile equality is used only for child-uniform actual queries. It does not assert equality of named nonuniform child responses.

## 3. Exact normalized physical count modulo eight

A physical Hamilton path at a child-uniform query decomposes uniquely into maximal child runs, interspersed with the two singleton ports. Consecutive equal child symbols are forbidden, and transitions follow the augmented quotient. If child i occurs h times, its physical fillings are counted by q_h of that child. Thus summing weighted words is the exact physical count, including every named cubic and quartic contribution. It is not an endpoint-product or coefficient-truncated model.

Divide each child's cover profile by its odd q1 in Z/8Z. A C3 then has normalized weights (1,2,2). Both kinds of four-child have exactly the same normalized weights:

    (1,6,4,0) modulo eight.

For the source/sink type, multiply (3,18,36,24) by 3; for the strong type multiply (5,22,36,24) by 5. Visits four times have weight zero. Every repeated child contributes a factor two, so words with at least three repeated children vanish modulo eight.

The C++ recurrence retains the 94,208 states with nine occupancies 0,1,2,3, the two port-presence bits, and at most two repeated old children. It counts the baseline normalized C3 words using transition multipliers 2 on a child's second visit and 1 otherwise. For a complete state, let W be the sum of its endpoint counts. Replacing position j by a four-child changes that state's contribution by the following *integer* multiplier:

    t(1)=1, t(2)=3, t(3)=2.

Indeed baseline weights 1,2,2 become 1,6,4. There is no modular division by an even weight. Multiplying each complete-state sum by t at each of the nine possible positions computes all placements simultaneously. The physical factor is 3^8 q1, congruent to q1 modulo eight, so source/sink residues are three times the retained normalized residues and strong-four residues are five times them. Since both factors are odd, every retained normalized gap of four remains the physical gap four.

The recurrence omits no term that could be nonzero: a fourth visit to the special child vanishes independently; three repeated children have total 2-adic valuation at least three; all other child visits and both physical ports are counted. No restriction is placed on the total polynomial degree or on the number of selected physical names in a port mask.

## 4. Exact finite domain and battery

The approved incoming census has 69,445 prime representatives modulo isomorphism/converse. The separately proved all-order endpoint lemma covers the 67,097 representatives with some odd S_i(Q)+E_i(Q), using (a,c)=(e_i,0), with physical contrast two modulo four. That lemma applies to arbitrary nonempty corresponding profile-equal children and is not extrapolated from an order-eight census.

The remaining 2,348 prime representatives have all endpoint parities zero. We test all nine possible positions of the four-child, hence 21,132 quotient/position cases. The first twelve coarse queries are the inherited C3 battery:

    (0,1), (1,0), (1,2), (1,3), (3,1), (3,3),
    (3,5), (0,7), (0,18), (40,256), (4,3), (130,128).

The complete normalized bank for these queries is residues.tsv. Exactly 13 cases are not separated, retained in initial_misses.json. The two extra queries

    (65,9), (3,32)

separate every one of those cases. selected_residues.tsv retains both queries for every position on the 13 affected quotient nodes, including redundant entries; the result is not inferred from their earlier C3 witnesses. The 14-query fixed battery therefore covers the entire residual domain. Adding all nine endpoint anchors gives 22 distinct coarse queries because (1,0) occurs in both lists.

witnesses.tsv retains one chosen query and both normalized residues for every one of the 21,132 cases. physical_witnesses.jsonl.gz retains three fully expanded physical realizations per case, one for each four-child isomorphism type, with both complete 30-row augmented graphs, both physical masks, and both physical residues. In total it contains 63,396 physical witness pairs, or 126,792 augmented graphs. Compression is lossless and reproducible, with a fixed gzip timestamp.

Every physical graph uses old names 0,...,27. Child i starts with names 3i,3i+1,3i+2; the selected four-child at j also contains name 27. On that child's ordered labels [3j,3j+1,3j+2,27], the outgoing-row templates are

    cycle plus source: [2,4,1,7];
    cycle plus sink:   [10,12,9,0];
    strong four:       [6,4,8,3].

Every other child uses the cycle 0->1->2->0 in its own local coordinates. Ports are the fixed physical names x=28 and y=29 with x->y. Coarse mask m expands to sum_i (7*2^(3i)) over its selected children, together with 2^27 exactly when child j is selected. The target reverses only between-child edges.

For an arbitrary named quotient, use the same isomorphism and entire-child permutation on both sides and transport both physical masks with it. If its certified representative is the converse, source and target exchange roles. Internal isomorphisms also act simultaneously on both sides, and child-uniform masks are unchanged within a child. No independent source/target response relabeling is used.

## 5. Every original-module family and the exact size-28 reduction

Assume F_T=F_U with T!=U for an original-module XOR target U, and |T|<=28. Choose a highest changed original strong node M with k children. The reviewed structural restrictions give

    28 >= |T| >= |M| >= 3k >= 27.

Thus k=9, and |M| is 27 or 28. If |M|=27, every child is C3. If T has no exterior, the reviewed 27-name theorem applies directly. Otherwise T has exactly one exterior vertex v. Because M is an unchanged-boundary original module, v either dominates all of M or loses to all of M, on both sides. There is no second changed highest node: the only remaining old name is v, which has no edge to change internally.

Use the separating core query from the reviewed 27-name theorem, and set both ports' incidences to v to imitate any vertex of M. Then M together with both ports is uniformly dominated by v, or uniformly dominates v. Every augmented Hamilton path must put v first, or last, respectively. Deleting v is consequently a bijection with the augmented core Hamilton paths. The exact source and target counts, not merely a modular contrast, are unchanged. This excludes the one-vertex exterior case without claiming any general ambient lift.

We are therefore left with |T|=|M|=28. M is the root, has exactly nine children, and its child sizes are eight threes and one four. All are nontransitive. The three-children are C3 and the four-child is one of Section 2's types. There is no exterior ambient or other highest changed node.

A proper original module at a prime substitution root lies within one child. Selecting the root oddly reverses Q; the target interiors arise from original-module operations within each child, possibly combined with the root's internal converse. Section 2 shows every such interior has exactly the source's complete cover profile. Therefore, at every child-uniform physical query, the full family target has exactly the same count as the quotient-only target, on the same names and mask, by occupied-run substitution. The 22-query coarse battery separates every root-containing original family. Families without a changed prime root cannot be full-response collisions by the inherited highest-node restrictions.

Hence every nontrivial source-relative original-module full-response collision has at least 29 old names. The finite battery is family-uniform over root-containing families, not over families changing only child interiors. It is not a claim about arbitrary pairs of family outputs or later-created modules.

## 6. Verification and promotion boundary

Run python3 verify.py and python3 -O verify.py. The verifier uses explicit exceptions, validates all pinned dependencies, recomputes the complete twelve-query bank and the two-query residual bank, redoes all four-child classifications and family profiles, validates every retained physical graph and mask, and checks the complete named local-response controls. The latter retain every coefficient, including zero, cubic and quartic coefficients, for the three canonical four-children and their literal converses. They are small exact controls, not a full named 28-vertex coefficient bank.

The new finite computations are producer results until independently reviewed. The proof only imports the reviewed quotient census and inherited frontier; it does not treat the earlier 27-name certificate as an arbitrary-ambient or arbitrary-child theorem. RELEASE_FILES.json is the publication allowlist and excludes binaries, exploratory extra-query banks, and runtime logs.
