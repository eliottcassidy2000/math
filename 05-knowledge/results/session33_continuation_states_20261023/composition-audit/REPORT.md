# Session33 independent composition and realizability audit

Verdict: PASS, with inherited congruence wording clarified by the producer. No mathematical correction to the new clock/port theorems was needed. Tests and report live outside the repository. No Git mutation, attachment-related material, or publication was performed.

## Results first

The independent checker executes 31,209 explicit checks. Normal and optimized Python runs produce byte-identical receipts. It imports no producer functions and does not rerun a producer checker as evidence.

- Clock service: 128 binary-mask banks, 8,192 demands checked by enumerating directly attainable owner assignments, not matching/min-cut algorithms.
- Möbius inversion: all 16,384 rank tables on three labels with nonempty values 0..3 satisfy the algebraic inversion identity; 120 have nonnegative inverse coefficients. Also checked a fractional rational profile.
- Matching matrices: 2/4/6/8 ports have 1/3/15/105 distinct rows and rational ranks 1/3/15/105. Independently full rank modulo prime 1,000,003 proves these rational full-rank lower bounds; binary Gaussian elimination gives ranks 1/2/4/8. Every supplied row separator was checked using independently implemented connectivity.
- G25: independently constructed the square-sum graph, proved leaf18 forces one endpoint, enumerated every path from18, independently enumerated local edge/port states, and recovered all ten paths and exact solution-index tuples. Original state counts 5/101/149 compress to 3/5/5 deterministic classes (including dead), with rational ranks 2/4/4.
- All context universe sizes agree: 240,784 / 11,920 / 8,080. Full tensor size is 1,203,920. Sparse tensor completeness follows because accepted gluing is exactly a connected spanning degree-at-most-two graph with two ends, and conversely every full HP has a unique local restriction.
- Simultaneous quotient weights give10; every individual quotient lift was reconstructed, checked connected and degree-correct, and its endpoints counted. This independently verifies the endpoint polynomial z18(z2+z3+z4+z8+z9+z10+z11+z13+z22+z23).
- P4 versus K1,3: independent permutation HP counts1/0, exact Schur responses both one-half times [[1,-1],[-1,1]], interior determinants both2.

## Clock theorem proof audit

The Hall deficiency formula includes the empty set and is nonnegative. Its integral version follows from integral bipartite capacities. Neighborhood closure preserves rank and weakly increases demand, so testing closed sets suffices.

For arbitrary nonnegative real demands, feasible-region maximization recovers every r(S): allocate each donor touching S entirely to one eligible label of S. This uses no simultaneous duplication. Boolean inversion then uniquely recovers nonempty-mask totals. Conversely equal mask totals imply the same Hall inequalities. This proves universal continuation minimality only in the specified demand language.

For integral banks, the maximizing demand can be integral, hence all nonnegative integer demands suffice. This is NOT valid for arbitrary real capacities when only integer demands are allowed: one-label capacities1/5 and4/5 both admit only integer demand0 but real demand1/2 distinguishes them. The producer correctly restricts its integer-demand statement to integral banks. The bounded binary-demand observer is likewise not universally separating.

The realizability iff is correct for arbitrary real-valued r with r(empty)=0: the nonempty-mask Möbius coefficients must be nonnegative, and then one donor per positive mask constructs a realization. Integrality of all coefficients is necessary and sufficient for an integral-capacity realization; an integer-valued r has integer coefficients automatically. Empty-mask donors are intentionally unrecoverable and must never be included in a claimed minimal observable.

The rank min(|S|,2) shared-gate hostile has three pair-mask coefficients1 and triple coefficient-1. Its submodularity does not repair the failure of independent-bank realizability. The theorem concerns banks without extra edge capacities/shared gates, not arbitrary capacitated networks.

Finite time gluing is sound if retained cells resolve all changes and demand histories are cellwise in the stated continuation language. A final constant unbounded cell is harmless. Per-cell source allocation lifts via owner lists, and finite adjacent assignments can be concatenated. Cross-time donor identity correlations are irrelevant only because handoffs are free and there is no lifetime consumption/holding constraint. “No recycling” here means no same-cut double spending, not a prohibition on reusing the same live negative indicator at disjoint cuts. Physical source ancestry and semantic eligibility remain external prerequisites, not outputs of the compression.

## Port theorem proof audit

The deterministic quotient is equality of rows; linear minimal dimension is field-dependent matrix rank; tensor-coordinate equivalences may be performed together because each original slice agrees against all other original coordinates. These are distinct statements and are correctly separated. Mod-p compression preserves counts modp, not exact counts or existence. At four ports all three row-sums are positive2 but vanish mod2.

The universal pairing separator is constructive and valid: choose ab in P but not Q; arrange Q edges into an alternating a-to-b path with a near-perfect continuation. Under P, edge ab is isolated from the remaining terminals. k=1 is vacuous. Removing a closure edge from a union of alternating cycles produces one spanning path exactly when the original union was one cycle. Parallel-edge conventions are essential but admit the stated private-subdivision realization.

Fixed-context orbit merging requires row equality, not merely equivariance under a simultaneous action on states and contexts. Averaged invariant context weights legitimately weaken the observer. The G25 quotient is fixed-context only; universal or future-operation stability is not established or claimed.

The endpoint polynomial restoration is sound because every live class has the same selected port set, so core endpoint contributions are representative-independent. Local endpoint monomials retain identities in coefficients rather than discard them. A scalar quotient without this property would require polynomial-response refinement.

The electrical hostile shares both response and Gaussian determinant factor; the Schur determinant identity therefore gives equal boundary-modification determinants. It does not refute algorithms retaining extra variables, interior factors, or Hamiltonian-specific state.

## Inheritance/status

Directly read current THM-2176 and THM-1975 frontmatter and MISTAKE-249. THM-2176 has PROVED abstract statements plus CITED knot applications. The arbitrary-monoid two-sided congruence assertion was historically wrong: right continuations yield a right congruence; the corrected theorem uses commutativity (or one could use two-sided contexts). Producer amended its inheritance paragraph to say this explicitly after audit feedback. The new static row theorem does not need any monoid assumption.

THM-1975 is VERIFIED bounded evidence with a flagged kernel derivation, not an all-order undirected transfer theorem. Producer correctly does not use it as one. The inheritance report points to THM-3121 for the later directed all-order proof; no theorem/status is promoted by this audit.

No uniform153N, LRC closure, physical source construction, arbitrary-source Collatz entry, all-n square-sum path result, or new literature-priority claim follows from these results. Literature citations in the ports text are outside this independent local proof/test audit; the finite rank calculations here do not depend on them.

Artifacts: audit.py; receipt.json; receipt.normal.json; normal.out; optimized.out.
