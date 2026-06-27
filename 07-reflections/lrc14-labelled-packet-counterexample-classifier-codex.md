# LRC14 Labelled-Packet Counterexample Audit

The classifier is the missing-picture bridge made executable.  In the shared
numbering it is HYP-2963: a bounded audit refinement of HYP-2961's conceptual
counterexample-family classifier, a companion to HYP-2962's fixed-margin packet
theorem and HYP-2956's labelled packet classification theorem, and a cross-check
against the HYP-2955 packet-migration gauntlet.

The important move was not adding another obstruction invariant.  It was
changing the classification question:

```text
Which rows are counterexamples?
```

became:

```text
Which labelled packet family would a counterexample have to inhabit?
```

On the bounded bank, the answer is clean:

```text
audited rows              21913
below threshold           0
tight rows                AP, GW 12->24
M<=2/27 low packets       7
unknown packets           0
```

The families are theorem-facing rather than aesthetic:

```text
q-witness                 direct denominator witness
AP/GW boundary            zero-open tight atoms
unit-petal                rank-0 C27 p=2 discharge
K33/state-lift            nonunit packet to HYP-2908/THM-572
covering-moment           exact-period boundary to gK8/L_y positivity
unknown                   zero-open, q>=14, unlabeled
```

The unknown count being zero matters, but it should not be oversold.  It is a
bounded AP-neighborhood and adversarial-bank statement, not a compactness proof.
The real theorem is still the global emitter: every primitive reduced residual
must produce the same packet fields.

## What Changed

S148 already had a strong gauntlet: no below-threshold named, single-swap, or
two-swap rows, and the `M<=2/27` layer was stable.  HYP-2953 and HYP-2954 said
the right object was a source-spectrum / missing-picture packet.  This session
made that object a classifier.

The key correction during the run was logical.  A row outside the S145 low
component atlas is not automatically unknown.  If it has positive strict
Haar/Baire interior, it is already discharged as loose.  Therefore "unknown"
must mean:

```text
zero-open + q>=14 + unlabeled packet
```

not merely:

```text
not in the low-frontier component dictionary.
```

That one distinction is the labelled-packet theorem in miniature.

## Swap-Chain Analogy

The arXiv:2606.22636 paper by Fu, Qin, and Wang was useful as a proof-shape
analogy.  Their binary fixed-margin swap-chain proof compares to heat-bath
moves, reduces to three rows, and separates scalar count sectors from Johnson
harmonic sectors.

For LRC14, the analogous split is:

```text
scalar sector       q_threshold, exact M, Farey excess
local reduction     AP/GW/petal/K33 packet atlas
harmonic sector     C27 owner, boundary owners, K33 incidence
comparison move     packet-preserving replacement
```

This is a good warning against premature scalarization.  Count data is useful
only after conditioning on the labelled fiber.

## New Pressure

The row

```text
drop(12,13)->add(26,36), M=3/37
```

is not a low `M<=2/27` atom, but it is a K33/state-lift family marker.  That
helps reframe "sporadic": AP/GW are isolated boundary atoms, while K33 and
petals are small labelled families with a theorem endpoint.  The classification
should therefore be family-plus-sporadic:

```text
families: q-witness, covering/open-Haar, K33, petal
sporadics: AP, GW, low named packet atoms
```

## Remaining Proof

The summit theorem can now be stated without poetry:

```text
For every primitive 13-speed residual S, construct P(S).
If P(S) is zero-open, q=14, and not positive-covering, then P(S) is
AP/GW, unit-petal, or K33/state-lift.
```

Everything else is safe by q-witness, positive open Haar/Baire mass, or
boundary-moment positivity.

That is still hard.  But it is a hard finite-labelled theorem, not a search for
one more scalar miracle.
