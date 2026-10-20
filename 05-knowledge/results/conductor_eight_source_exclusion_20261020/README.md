# Conductor-eight source exclusion

20 October 2026. The component proofs and exact-dependency synthesis are
**PROVED and independently audited PASS**.

This packet classifies every symmetric genus-four polynomial-plane degree
semigroup and every possible nonprincipal terminal degree pattern. All the
patterns are realized by explicit immersed curves. Actual constant-Jacobian
coefficient identities nevertheless exclude the quartic/sextic and
cubic/quintic sectors.

Combined with the separately audited all-quadratic-semigroup theorem and
the conductor-six exclusion, the
[synthesis](CONDUCTOR_TEN_SYNTHESIS.md) proves that a noninvertible planar
Keller map under the primitive monic source hypotheses has generic
source-conductor degree at least ten. The final synthesis and publication-copy
audit is recorded separately in `review/SYNTHESIS_PUBLICATION_AUDIT.md`.

## Proofs and controls

- [Complete genus-four classification and quartic/sextic exclusion](CONDUCTOR_EIGHT_CLASSIFICATION.md).
- [Direct cubic/quintic exclusion](CUBIC_QUINTIC_ADDENDUM.md).
- [Explicit immersed sextic/nonic realization](sextic_nonic/SEXTIC_NONIC_REALIZATION.md).
- [Classification and quartic/sextic audit](review/AUDIT.md).
- [Separate cubic/quintic audit](review/CUBIC_QUINTIC_AUDIT.md).
- [Exact synthesis dependencies](SYNTHESIS_DEPENDENCY_HASHES.json).
- [Publication-copy hash lineage](PUBLICATION_COPY_LINEAGE.json).

The two primary checkers have 57 and 31 exact controls; the independent
checkers have 42 and 32. All ordinary/optimized outputs agree. The nonic
realization contributes 11 more exact checks. Reproduction commands and
the distinction between universal arguments and finite controls are in
the proofs and audits.

The bound is not a function-field degree bound for the planar map, not a
uniform special-fiber bound, and not a resolution of the planar Jacobian
conjecture. Conductor ten and larger remain outside this packet.
