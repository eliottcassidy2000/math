# Independent uniform-block start tomography review

Status: PASS. No mathematical correction required. Producer proof a3e44ac3ef099401971e055e42186ebf2674328483f12a1c00b3513543be5345 and checker b7e944ad3278628b9470ed7347aaafd610bea988dbc01951926ccbcad70da0e5 are pinned in both receipts.

## Complete-word and functional audit

For query (1_I,0), no old name points into y. With x before y, the old segment between them must be empty; the prefix ends in I or is empty, and the suffix is arbitrary. With y before x, y must be first, the intervening old segment is nonempty and ends in I, and the suffix is empty or starts outside I. These give exactly (1+c)(1+h)+c(1+h-r).

For query (0,1_I), x is a universal source and must precede y. The middle old segment is empty or ends in I, and the suffix is empty or starts outside I. This gives exactly (1+c)(1+h-r). Both port orders and all empty-segment cases are covered. Child occurrences separated by ports are distinct ordered-cover runs, so they are correctly retained in products of series.

Quotient converse transposes the run matrix, fixing h and swapping r,c. Direct expansion yields the two contrasts 2(c-r)(1+h) and (c-r)(2+h); their stated linear combination is 2(r-c). The functional is linear but not multiplicative, and it is applied only after this formal cancellation. Equal complete profiles, not equal endpoint marginals or equal interiors, are the required hypothesis. A union of whole physical child sets is essential.

## Strictness and root-family audit

The cyclic aggregate-start formula has the correct orientation and indices. The Stirling recurrence log-concavity expansion is exact, including its wider minor; positive-support ratio comparison gives the stated sign, with zero-extension covering endpoints. When one transitive size is one, the k=2 minor is still 2^a-2^b and is nonzero. The arbitrary third child needs order at least two to supply positive q_2. No log-concavity condition is imposed on that child.

Equal transitive sizes, or a singleton queried third child, are not claimed to be covered by this argument. The retained all-three-start-skew-zero example still has an unequal actual response query; this is properly called blindness of this invariant, not a response collision.

For a root-containing family, all proper original modules lie inside one cyclic child. The root quotient is therefore conversed. Scalar strictness is detected by (0,0); hereditary equality supplies precisely the child profile equality needed above. Thus the fixed three-probe battery is valid even with crossing selected modules inside a child. It is original-versus-target detection at the whole root, not general all-pairs identification or an arbitrary strong ambient lift.

## Independent executable controls

The audit imports no producer code. It independently reconstructs all retained physical graphs and masks, computes starts by predecessor subset DP on the converse, checks internal cover profiles by literal permutation/cut enumeration, and computes the actual two-port counts with its separate C++ predecessor DP.

It rechecks all 1,165 retained start/query identities, all sixteen root families and their equality-branch profile restoration, and the strict minor signs through transitive order twenty using the independent inclusion-exclusion formula for ordered set partitions. Both producer modes replay exactly, and both audit modes are byte-identical. Full counts and source pins remain in the producer and audit receipts; finite checks validate the conventions rather than replace the all-order proof.
