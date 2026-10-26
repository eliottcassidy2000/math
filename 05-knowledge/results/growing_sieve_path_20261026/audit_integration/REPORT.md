# Independent review of the integrated synthesis

**PASS after one wording repair.** This review covers the packet's
`README.md` and `METHOD.md`, compared with the three frozen proof lanes and
their audits. It reviews the additional mathematical deductions rather than
rerunning the already audited exhaustive computations. Final reviewed hashes
are bound in `MANIFEST.json`.

## Exact consequences of tournament recovery

The additional identities are sound. Equivariance of the promised inverse
gives `Aut(U(T))=Aut(T)` on the same vertex set. Every Hamilton path of T
uses only arcs in U(T), and U(T) is a subgraph of T, so
`L(U(T))=L(T)`. The raw and automorphism-normalized path counts therefore
agree as stated. These facts also prove the isomorphism-class equivalence.
They supply a faithful carrier for tournaments on the image of U; they do
not characterize arbitrary digraphs in that image or produce an A000568
recurrence. The synthesis makes that distinction explicitly.

The fixed-path cube has `(n-1)(n-2)/2` free arcs. Its full named permutation
guard vector is its labelled Hamilton language, which the new theorem makes
injective. Thus the exact count `2^binom(n-1,2)` on X_n is justified at all
positive orders. The absorbing failure state is external to X_n; it is not
silently counted as another presentation. This closes the named-guard
question and leaves the different fresh successor-class observations open.

The five-vertex example and cut direction are correct: the last row changes
from 14 to 12, removing 4->1, and deleting 0 reveals the order
`4 -> 3 -> 2 -> 1`. Its remaining forward arcs do not change that order.
The strong-source restriction on the deletion lemma remains explicit.

## The field cover and its arithmetic boundaries

The determinant bracket and homogeneous-space statements agree with the
audited field proof. The root/pole conjugation is routed to the exact gauge
in that proof. At a fixed pole, the two sheets are converse Paley tournaments.
Their arc-transitive automorphism group carries any one used arc to every
arc, so the new deduction `U(T)=T` for this family is valid. The figures
`H=189`, `|Aut|=21`, and `tc=9` have the correct distinct meanings.

Consequently the raw 7/21 exclusions do transport to the *image* of U,
where H is unchanged. This says nothing about arbitrary partial digraphs
and does not re-prove the exclusion theorem from the order-21 group.
The through-eleven tc result keeps its finite scope.

The mod15 count 3072 is the complete reachable ambient observer, not an
assertion that the prime-filtered subset is transition-closed. The two
Euclid pairs, their twin radii, radicals, heights and shared modular state
agree with the frozen field receipts. The squared-Pythagorean abc packet is
correctly separated from the earlier centre identity; their heights and
radicals differ. Neither collision is used to refute abc or every possible
parity-dependent bound.

## The sieve's finite horizon and unbounded quantifiers

The strict-wheel count uses lower endpoint 3 and upper endpoint H. Its
prime-equality correction is exactly the actual twin pairs whose lower
prime is at most B. The recurrence and its specialization at
`B=floor(sqrt(H+2))` agree with W2--W3: the argument decreases for H>=3,
and the next prime is strictly above the relevant square root.

The square boundary, its delayed cases, and the cutoff-61 witness are
reported at their exact scope. In particular `4487=7*641` is removed,
`4757=67*71` survives the wheel, and 4759 is prime. The synthesis does not
promote this finite first-failure experiment to a universal square pattern.

The three limit assertions are different and are stated correctly:
arbitrarily large survivors at each finite strict stage; no positive
integer surviving all strict stages; and pointwise, bounded-height uniform
correctness of the exception-restored tests. The compatible class -1 does
not supply a bounded positive lift. No infinitude result follows merely
from the exact finite recurrence.

## Observation criteria and repair types

The factorization criterion is exact on the stated source and observation
image. The sidecar-cardinality converse is only set-theoretic; the text
explicitly denies an effectiveness or complexity conclusion. In ordinary
set-theoretic foundations, indexing the at-most-K target values separately
in each fibre gives the asserted sidecar. This is not a lower bound against
all primality algorithms.

The future-profile criterion retains operation names, legality and the
failure state, and therefore is the correct coarsest stable refinement.
The finite partition algorithm is not claimed to terminate at a uniform
depth on an unbounded class. Changing to unnamed or randomly resampled
operations is correctly treated as a different question.

Structural inference of missing tournament arcs, a composable one-bit
orientation extension, and height-dependent exhaustive arithmetic tests are
three distinct repairs. They do not establish equivalences among famous
open problems, solve arbitrary image realization, or supply unbounded
positive survivors.

## Repair applied during this review

The initial METHOD wording said augmentation recorded shape and arithmetic
evaluation recorded prime weights. Taken literally, these scalar outputs
need not recover either object. The owner replaced it with the precise
statements that augmentation counts proper-divisor monomials and evaluation
at primes computes their weighted sum, explicitly noting that neither scalar
generally recovers the shape or prime labels. The repaired text was reread
and passes. No other mathematical correction was requested.

This report does not audit new material outside the two named synthesis
files, nor does it claim a new independent replay of the producer tests.
The proof and computational audits linked from the synthesis remain the
authorities for their respective full universes.
