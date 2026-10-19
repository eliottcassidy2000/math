# Independent review: guarded parabolic compactification

19 October 2026. **Verdict: PASS, within the stated scope.**

The square-guard descent theorem, its unlimited-episode first-crossing bank
corollary, the repeated-source crossover exclusion, the exact K phase cycle,
and the source-sensitive compact-boundary equivalence are supported by the
written proofs and independent exact controls. No mathematical correction to
those proofs was required. The reviewed publication revision makes the
integer domains explicit and repairs references and the reproduction command.

This is a restricted, guarded result. It proves no arbitrary positive-integer
bank entry, universal rational coefficient stopping, or global Collatz
convergence statement.

## Symbolic audit

**Square guard.** A lawful word has a unique source class modulo its full
power-of-two modulus M. Therefore, if W occurs twice consecutively, its first
displacement belongs to M times the integers. For a contracting branch the
nonnegative displacement is strictly below its carry e, while the displayed
odd-position sum gives e <= (a-1)/4 < M. A positive displacement is impossible.
The exact-return case must remain: the two points of the trivial cycle supply
immediate examples. The shorter matching-bit statement uses exactly the same
argument with modulus 2^ell and includes ell=0.

**Unlimited episodes.** The inherited zero-displacement argument does not use
an episode cap. Before the first crossing, positive carry and noncontracting
boundary slopes put every proper macroblock endpoint above the source. The
expansion-prefix and reset-prefix comparisons cover every intermediate Terras
time. Hence a return has the claimed minimum and least period. Independently
lawful adjacent companion blocks compose to a companion return; its initial
two even steps give the smaller periodic point. The thirteen short exclusions
and surviving long source guards hold independently of the number of
episodes. Applying the new square guard therefore removes the episode cap
only under the additional self-matching hypothesis. No source-coverage step
is implicit in this deduction.

**Repetition and crossover.** Multiplying the repeated affine numerator by
the odd number a-M proves the stated divisibility equivalence. Endpoint
integrality authenticates the unique complete parity cylinder. The exact
valuation budget and displacement lattice follow. On the expanding side,
the displacement is a positive lattice multiple, yielding the nonlinear-ratio
bound. On the contracting side, a nonzero displacement must be negative when
the lattice spacing exceeds the carry. The fixed-point exception is correctly
retained. The near-neutral bank family has the stated bounded carry and
growing binary modulus, so the exclusion genuinely applies to its lawful
double repetitions.

The parabolic reciprocal map is a Möbius map. Its translation coordinate,
the cubic-coefficient invariant distinguishing it from the quadratic germ,
and the different parameter terms are correctly calculated. Neither the
real bottleneck clock nor formal iterability supplies missing parity guards.
The inherited rational three-cycle at -29/16 and algebraic parabolic cycle
at -7/4 are kept separate.

**K arithmetic phase.** The canonical source congruence gives the exact
integer phase recurrence, with range 1 through D-1 for L>=2. Coprimality
implies that the phase orbit has the order of Q modulo D. Its 146,768 distinct
phases yield exactly that many real subsequential mantissa limits, while all
unscaled residues approach the same rational 2-adic anchor. The source-height
bound and quantitative rational-anchor repair threshold follow without a
floating-point estimate.

**Boundary equivalence.** The simultaneous real/2-adic closure and continuous
extension are valid on the declared compact space. In Theorem C, contracting
prefix cylinders form an open cover of the actual unbounded-source boundary.
A finite subcover gives a common time bound and a finite maximum carry
threshold. Compactness forces every sufficiently large realized source into
that cover. Conversely, matching modulo 2^J transfers all J relevant parity
prefixes to a boundary point; noncontracting prefixes cannot descend below
the unshifted source when d is positive. This proves both directions,
including the empty-boundary case.

The source set S_d keeps the displacement equation. Each fixed contracting
word permits at most one positive source for that displacement. Diagonal
selection therefore produces an infinite noncrossing symbolic bank path from
an unbounded realized sequence. The reverse inclusion is appropriately not
asserted. In particular, the rational anchor must be tested for this
source-sensitive criterion only when it is actually in that boundary. The
endpoint filter does authenticate an entire completed word, but cannot by
itself prune unfinished infinite paths. The explicit all-odd boundary
examples correctly rule out a globally uniform finite repair deadline.

## Independent arithmetic

The supplied verifier produces its stored receipt byte-for-byte in both
normal Python and Python -O: **137,320 explicit checks**. All supplied
checksums pass.

The separate [independent checker](independent_controls.py) imports no supplied verifier. It
constructs affine carries by a closed odd-position sum, and reconstructs
source cylinders by one-bit lifting using actual parity rather than the
closed modular-inverse formula. Its normal and optimized outputs agree
byte-for-byte and perform **1,124,584 explicit checks**, including:

- All 8,190 binary words of lengths 1 through 12; exact repeat budgets and
  73,674 positive repeated lifts.
- 6,423 contracting-square inequalities at the least positive lawful source.
  Their negative affine displacement slopes certify all positive lifts.
- The four independent lower/companion interfaces, every K/V expansion
  prefix, every proper F/H reset-prefix comparison, the inherited exact
  episode carry constant, and all thirteen short all-lift exclusions.
- Enumeration of the complete 146,768-element modular phase orbit, rather
  than only the supplied multiplicative-order certificate.
- Exact near-neutral words of 336, 1,222, and 2,108 Terras bits and their
  genuine square cylinders.
- All 17 declared rational-anchor controls. In particular, d=4 stops in
  coefficient at 92 steps with sufficient K-depth 4; d=12852 stops at 135
  steps with sufficient K-depth 5.
- Endpoint-equation controls, positive-overshoot hostiles, exact returns,
  and all-odd boundary obstructions.

These computations are finite controls. The written symbolic and compactness
arguments, rather than the finite enumeration limits, supply the all-size
quantifiers.

The terminology citation to Rozier and Terracol was checked directly against
[Definition 1.2 and its surrounding discussion](https://arxiv.org/html/2502.00948v5).
It supports the distinction between coefficient stopping and actual descent;
the note does not attribute a rational stopping or pointwise entry theorem to
that source.

## Reviewed revision and reproduction

The original proof SHA-256 was
`981bb25c14ea474f84b800bf50622ff9711410b1f6caa0252b581833a3527416`.
The [reviewed publication proof](../GUARDED_PARABOLIC_COMPACTIFICATION.md) SHA-256 is
`563a91175935b7eb1bd8f5663f3166342848db976150c8f0f04d9f7bc3686114`.

The difference was inspected: it explicitly declares positive integer d in
Sections 4–6 and nonnegative integer ell in A2, adds repository-relative
links, and gives the full optimized reproduction command. It changes no
mathematical claim or proof. The supplied verifier and receipt are unchanged.

Supplied verifier SHA-256:
`d96a48bb79e5758330c4e9df64160569b353713e5eef260e3fe9a365c8ad5326`.

Supplied receipt SHA-256:
`d9dc7001a55cd0d076d2e8957f86fbed6d97e9199205d8130bcc3bfebc3c7b00`.

Independent checker SHA-256:
`d5a3a7e6881e357dc2f71e2ccfe8eca27e7440f0dccc50374c1d092cd8f37e08`.

From this review directory, run:

    python independent_controls.py
    python -O independent_controls.py

Both outputs reproduce [the normal receipt](independent-normal.json), which is byte-identical to
[the optimized receipt](independent-optimized.json). The checker needs only the Python standard
library and no external files or network access.

