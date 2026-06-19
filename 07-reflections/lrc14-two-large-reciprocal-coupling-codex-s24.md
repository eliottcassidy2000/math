# LRC14 Two-Large Reciprocal Coupling

**Source:** codex-2026-06-19-S24, HYP-2633 / T881.

HYP-2632 finished the finite phase job: the two-large repeated tail has a
small signed `chi_7`/affine/Q table, and the Dedekind packet says exactly why
frequency shells should be exposed before taking absolute values.  The new
session checked the next proof obligation: whether that finite sign table can
be read directly as the sign of the actual reciprocal relation-lattice tail.

It cannot, at least not before the lift lemma is proved.  The stored scout
`lrc14_two_large_reciprocal_coupling_codex_s24.py` compares HYP-2632 packet
weights with exact cumulative reciprocal sums through height `H=16`.  The
guardrail is sharp: `42_QR_a2` and `42_QR_a4` both have finite packet weight
`-25U`, but their bounded reciprocal lifts have opposite signs at `H=16`.
The finite-zero affine lane `(0,2)` also has a visibly nonzero bounded lift.

This does not weaken HYP-2632.  It tells us where the proof must spend its
next effort.  The two-large theorem needs the following chain:

```text
finite chi_7/affine/Q Dedekind packet
+ finite low-height wall deletion
+ residue-lift equidistribution on relation lattices
+ signed summation by parts inside additive frequency shells
=> two-large reciprocal tail bound using -108U+54U, not 162U.
```

The right analytic statement now looks like an Abel-summation lemma.  For each
additive-frequency shell, write the reciprocal tail as a Stieltjes sum over
relation-lattice height.  After deleting the finite low-height wall ledger,
prove that the cumulative discrepancy of residue lifts stays bounded or has
enough cancellation, then integrate it against the `1/prod(n_i)` denominator.

The assumption challenged is the tempting one that finite residue classes or
finite packet signs are enough vertices.  They are not.  Runners, raw supports,
and residue classes destroy the integer lift distribution; finite frequency
shells alone still destroy low-height resonance data.  The useful vertices for
this proof step are proof obligations: residue-lift equidistribution,
additive-frequency summation, finite character kernel, wall deletion,
exact-period packets, and boundary-face cancellation.

LRC(14) remains open.  The progress is that the remaining analytic target is
now a concrete lift-discrepancy theorem rather than an undifferentiated
cotangent/Dedekind wish.
