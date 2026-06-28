# LRC14 Owner-Cut Dual Certificate Synthesis

HYP-3410 made the owner-cut picture concrete, but the size-3 frontier still
left the next theorem shape too loose.  This pass turns the cuts into signed
owner-current certificates with a declared zero level.

The useful result is small and exact.  All three current HYP-3410 mixed fibers
have margin-1 labelled currents.  The first two have one-label selected cuts
or one-label island certificates; the `(72,20)` `10->20` frontier needs the
positive-debt cut `{2:g2, 11:g1, 13:g1}`.  That is the expected 14-factor split
in miniature: one even covering label plus two binding labels.

The surprising signal is the repeated `13` in the positive-debt Sophie audit.
For the old owner leak, core size `1` and cut size `2` give channels `13` and
`5`.  For the new frontier, core size `1` and cut size `3` give channels `25`
and `13`; the same frontier also has top finite-BDH variance `13:g1=49/50`.
This is not enough to assert algebraic necessity, but it is strong enough to
justify testing a bounded owner-current theorem rather than continuing to add
analogy names.

The Krasner lesson also became sharper.  The common owner core is exit-mixed
in every fiber, so local stability cannot mean "same common p-adic core."  It
must mean stability of the full owner/contact root packet, or the first
non-core owner current becomes named debt.

After rebasing over S257, I read this as a local certificate counterpart to the
Galois/GW residue-magnitude thread.  The frontier cut has one even-cover label
and two binding labels, matching the `14 = 2*7` split, but it does not replace
the global `q == 1 mod 3` GW-doubling criterion.

After S258, the guardrail is stricter: this is not the critical path unless it
feeds the decorrelation floor inequality.  Its best role is as a finite sidecar
for owner-mixed fibers.

After the HYP-3416 recursive quotient ladder landed, the placement is clearer:
this certificate is a labelled owner-current layer inside a larger quotient
stack, not an independent proof program.

Next work should extend the HYP-3406/HYP-3410 bank and ask whether every
surviving mixed residue/height fiber has one of: unit-island current,
positive-debt hitting current, exit-pure charal recursion, strict-open/AP-GW
terminal route, state-lift debt, or a new named residual.
