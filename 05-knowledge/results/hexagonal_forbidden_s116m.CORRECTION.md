# Correction to the historical centred-hexagon/Fano identification

**Current correction — 2026-10-25.** The original
[`hexagonal_forbidden_s116m.out`](hexagonal_forbidden_s116m.out) is preserved
byte-for-byte as computational history. Its assertion that the centred
hexagonal neighbourhood “IS” the Fano plane is not valid as an identification
of its existing Euclidean incidence or symmetry.

Seven labels do not specify seven Fano lines. The centred root hexagon has
three opposite-pair lines through its centre, and its central inversion has
cycle type `1·2³`, which preserves no Fano structure on seven labels.

The strongest repaired bridge is positive and explicit: complete the
Eisenstein lattice, choose one of its two norm-seven ideals, retain the
triangular two-cells, quotient to a torus, and choose one of its two face
classes. Each chosen class is a Fano design. The two ideals and two face
classes give exactly the four special free Campedelli codes of the marked
hourglass arrangement; inversion exchanges face classes rather than fixing
one. See the [audited construction](session61_hadwiger_six_seven_20261025/torus-fano/TORUS_FANO_CAMPEDELLI.md).

The missing sidecars are quotient lattice, triangle incidence, face class,
chirality and branch marking. No Euclidean plane chromatic lower bound lifts
backwards from the quotient K7. The original output's other claims and counts
are not revalidated by this correction.
