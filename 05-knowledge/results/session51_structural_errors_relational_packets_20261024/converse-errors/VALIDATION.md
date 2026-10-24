# Validation and source boundaries

All numerical claims are exact integer claims. No statistical inference or
full unrestricted order-ten enumeration is used.

- `check.py`: own refinement, exact directed isomorphism/automorphism search,
  160 ordinary positive representatives, unique converse partner for each,
  explicit anti-map and odd-power anti-involution witnesses, and actual stable
  cell action. Fifteen regular-nine representatives are rebuilt from pinned
  almost-regular-eight data. The 40 nonstrong ten-representatives are explicit.
- `matrix_check.py`: separate complete labelled matrix construction from
  permutations and pairs of disjoint permutations; translation quotient and
  every twisted fixed set for degrees one and two at p=5. No source imports.
- `prime_check.py`: independently memoized capacity recursions for ordinary
  and symmetric regular-matrix coefficients, and literal subset necklaces at
  p=3,5,7. It does not enumerate all p=7 matrices.
- `series_check.py`: all ordered words through total order16 in an explicit
  alphabet of actual strong cores, comparing direct reversed-dual fixed words
  against F=(1+Fs)A(z²,t²). The derivative factor2 is checked. The deliberately
  wrong t-marker substitution first disagrees at order14. Scope is bounded
  alphabet, not all tournaments.
- `port_check.py`: all Hamiltonian paths of each of the24 rigid self-converse
  positive-error ten-classes, checked by subset DP and literal path traversal.
  Blue paths are counted both from anti-involution half-systems and from a
  literal normalized-word staircase reflection test. Each node has both blue
  and black witnesses. These are ports, not edge pairings.

Each program is replayed under normal Python and `-O`; paired JSON outputs
are byte-identical. Checks use explicit RuntimeError, not Python assertions.
The main checker pins every inherited source before using it. A real-copy
relocation replays main/series/prime/matrix/port checks without access to the
original working directory. The manifest checks all non-cache packet files.

Inherited completeness is explicit: the strong-order-ten classification and
its representative receipt, the complete almost-regular-eight parent slice,
and the strong-nine structural exceptions. The all-prime coefficient formulas
are Session40 regular-block work; anti-involution existence is THM-647. Only
THM-644(a)'s proved grid-symmetry equivalence is used for blue incidence. The
rigid path-to-presentation bijection is reproved in this packet, not imported
from the older partly verified general fibre claim.

The exact order-ten positive-error merged extraction, its general weighted
SCC transfer, and the scoped prime-cell zero/positive-error split are the
consequence objects. No all-blue classification theorem, A000568(10) total,
blue/black adjacency, induced-edge count, or unqualified novelty is asserted.
