# Independent audit: bounded charge blocks and precedence

23 October 2026. **VERIFIED proof review / FINITE-EXACT replay.** Reviewed
`../algebra/BOUNDED_CIRCUITS_AND_PRECEDENCE.md` and `../algebra/check.py`.
Normal and optimized final replays both pass **210,035** checks and match their
frozen receipt. This audit did not edit the algebra files.

The paired zero-sum argument is sound. Equal margins permit pairing the
positive and negative changed slots; equal charge gives a zero-sum sequence
of paired differences. Extracting minimal zero-sum subsequences changes at
most2D(G) actual slots per block, preserves both integer margin and group
charge, and changes each differing slot once. The resource and shared-owner
conclusions genuinely use the all-mixture assumptions; they do not follow
from endpoint feasibility or the additive invariant alone.

The p-group proof is correct: the augmentation ideal in the truncated group
algebra has nilpotence degree1+sum(p^a_i-1), while a product with no nonempty
zero-sum subsequence would retain identity coefficient1. The displayed
minimal zero-sum sequence proves sharpness, including the separately handled
trivial group. The abstract2D(G) supply is sharp, but its occurrence in a
native tournament is not inferred. The original paired groups give8 and132
respectively, hence bounds16 and264 under the declared physical hypotheses.

The minimum-cardinality subset-sum DP maintains enough information: at each
processed prefix the best nonempty witness for each group element can be
extended using the next, previously unused index. Comparing against the
singleton next index handles empty-prefix starts. This gives an
inclusion-minimal witness without a false greedy single-deletion step.
At most m extractions, each of at most D(G) entries and |G| DP states, give
the stated O(m D(G)|G|) group-transition bound; witness storage and native
authentication remain separate costs.

The vacancy path has exactly suffix mixtures because its only collisions
are new p_(i+1) against old p_(i+1). Every proper nonempty suffix has nonzero
Z/2 charge. It therefore has exactly two guarded states and requires an
m-slot event inside its complete binary menus. No implication SCC is needed
for this obstruction. Five slots first exceed2D(Z/2)=4; three already expose
why separated charge blocks can violate retained precedence. The paper now
distinguishes those two claims.

Three requested precision repairs were accepted before this final audit:
both alternatives must satisfy the same prescribed original per-slot phase;
the extraction paragraph uses minimum-cardinality DP rather than ambiguous
deletion minimization; and the three-/five-slot thresholds are separated.
In particular, equality in the paired target quotient alone would not have
licensed equality of a finer per-slot phase. No remaining proof or scope
error was found. The196-state mod8 reference concerns Session45's different
sixteen-slot menu, not the earlier194-state menu; its provenance is explicit.

The subsequent balanced U/V extension also passes independent proof review.
For even m, alternating target tags give integer changes(-1)^i; the two
endpoint source labels contribute the only nonzero Z/2 changes. A proper
suffix of odd length violates the target margin, while one of even length
retains exactly the last endpoint charge. Empty and full are therefore the
only accepted states. The first/last pair and adjacent middle pairs give a
budget2 route when occupancy conflicts are removed. This is an exact
balanced-interface hostile from even m=6 onward against the4-slot bound,
still abstract rather than an asserted native harmonic realization.
