# Independent first-descent and template-switch audit

**PASS.** No mathematical correction required. Reviewed RESULT.md, the author's verifier and its receipts, S179's actual parser and grammar proofs, and the named older run-block source. No claim of universal Collatz convergence or historical priority is certified.

## Proof review

1. The residue modulo 2Q is odd and unique because P and B are odd. Expanding the final carry after prefix j leaves precisely one term of valuation S_j; all later terms have greater valuation. The final odd-output congruence therefore forces the exact intermediate valuations, not merely divisibility. Every earlier positive state is strictly above the source because its slope is at least one and its carry is positive. The final strict threshold is exactly n>B/(Q−P). Removing finitely many progression members preserves the stated all-integer and odd-integer densities.
2. For W_q=(1,2)^q2^r, the slopes before correction and the first-crossing choice of r are correct. Equality of a positive power of3 with a positive power of2 is impossible. Both the closed carry and displayed examples agree with independent literal trajectories. The stated r=floor(alpha q)+1 follows, and all valuations remain in{1,2}.
3. The S179 mismatch statement is exact: after the first2q letters the next native valuation2 gives x=1 modulo8, so v2(x+5)=1. Since x+5=(9/8)^q(n+5), the source has v2(n+5)=3q+1. Thus the preceding template stops at exactly the claimed boundary.
4. A single final reset needs the stated multiplicative inequality; without it the positive affine carry prevents descent. The proposed bounded-valuation replacement charges extra steps instead. Irrationality of alpha follows by unique factorization. Density of its fractional multiples implies rho_q approaches1 from below. The finite fixed-margin obstruction is correctly restricted to representing these entire words, not eventual descent using another decomposition.
5. The non-context-free proof is sound. Inverse homomorphism followed by intersection with a^+b^+ recovers exactly a^q b^r(q), since the number of1 symbols fixes q uniquely. Any nonzero nonnegative integer period in a linear subset of its Parikh graph would have v=alpha u by bounded discrepancy. Irrationality excludes all such periods, so semilinearity would force finiteness. This uses the classical Parikh theorem already cited in S179; a direct citation is recommended for portability.
6. The rotation proof for arbitrary expanding words is valid: rotate after a minimum of the proper partial sums. Unwrapped prefixes are nonnegative; wrapped prefixes are positive after adding the positive total. Repetitions retain the property, and valuation-two correction has its first crossing at its last letter. No source coverage theorem results.

## Independent computational checks

The author script was executed in normal and optimized Python, with identical outputs matching 880 sources and135,132 native steps. The separate `independent.py` imports no author code. It finds compatible residues by explicit lifting and literal odd-map iteration, rather than modular inversion of affine triples. It checks every first-slope-crossing word through length6 with letters1..4 (49 words), three positive sources per word (147 total), both displayed numerical examples, and819 expanding words through length8 with letters1..3 for a legal rotation. Normal and optimized independent receipts are identical.

These finite controls support the proof and do not establish its universal quantifiers by enumeration. Neither they nor the proved infinite progressions provide a global entry selector, exclude all nontrivial cycles, or establish convergence.

## Presentation/source requests

Replace session-specific `/workspace/shared/session179/` references in the proof with repository-relative links before publication. Add the inherited Parikh source directly beside the formal-language proof. These are portability and attribution improvements, not mathematical blockers.
