# Independent prime-circulant cutoff review

Date: 2026-10-21 UTC. Verdict: PASS. This is a separate final refinement; the degree47 and24-coefficient checkpoints remain valid historical intermediate results, not silently rewritten.

## Primary source and analytic proof

The reviewer independently reopened Terence Tao, An uncertainty principle for cyclic groups of prime order, arXiv:math/0308286v6, printed page5: https://arxiv.org/pdf/math/0308286 . The displayed Cauchy–Davenport statement is the lower bound min(p, |A|+|B|-1) for nonempty subsets of Z/pZ. The later proof-ending MAX is an apparent typographical error and is not the statement invoked here.

For a separator S and a surviving nonempty SCC split X→Y, no Y-vertex has an outneighbor in X. Adding zero to the connection set therefore gives Y+(D union {0}) contained in Y union S. The latter has size p-|X|<p. Cauchy–Davenport consequently cannot saturate at p; its other term yields |S|>=|D|. Deleting a vertex's full outneighborhood attains that bound while leaving the vertex as a sink with another surviving vertex. Hence the prime circulant has exact vertex connectivity (p-1)/2.

The proof does not extend this additive estimate to composite orders. The specified Z/9Z connection set {1,3,4,7} is a genuine tournament connection set. Removing {1,4,7} leaves {2,5,8}→{0,3,6}. Its pertinent sumset has size6, not the prime-theorem lower bound7. Independently exhaustive connectivity computation confirms its minimum cut is3, strictly below its degree4.

For prime center order p>8, exact q<=8 and p=qm force q=1, so the tournament is circulant. At p>=37 connectivity is at least18, and deleting a directed16-cycle leaves a strong complement. Since p>=32, the half-support requirement is satisfied. An external vertex supplies at least eight entries or exits, yielding Delta>=14 via the previously reviewed lemma and self-converse entry dual. Thus prime center orders37,41,43,47 are excluded.

The final degree45 conclusion uses both exclusion of47 and the prior support restriction excluding even46. No composite cutoff is inferred from prime connectivity. Removing the four prime indices from the prior24-position support leaves20 potentially unknown coefficients, including composite39 and45. This is a necessary support list, not existence or evaluation of its coefficients. Low-order coefficients through9 retain their explicitly inherited status.

The numerator degree bound becomes45+12=57, so the inherited order18 recurrence is valid from n>=58. The independent prime lemma does not reprove the inherited rational identity or minimal eventual order, and does not claim earliest recurrence onset.

## Fresh scoped controls

Copied producer runs under normal Python and -O each pass17,159 always-active checks and agree on every nonvolatile receipt field. Actual clock fields differ.

A separate reviewer script uses cyclic bitmask unions for all49,961,16,129 nonempty subset pairs at primes3,5,7, respectively. It independently computes vertex cuts using transitive closure for all14 corresponding circulant orientations and the composite9 hostile. It also verifies the final support arithmetic. It passes17,160 always-active checks both normally and under -O. This is tiny scoped evidence, not a census of large prime-order centers or a replacement for Cauchy–Davenport.

## Files and publication scope

Compact current deliverables: REVIEW.md, independent.py, independent_results.json, source_pins.json. Original producer packet contains its frozen checker and receipts. Copied normal/optimized review directories need not be published.

Only session22-prime-review was created or edited for this audit; the earlier q2 review source pins were also refreshed after its authorized status-only editorial freeze. No core mathematical packet, Git state, or remote repository was changed by the reviewer.
