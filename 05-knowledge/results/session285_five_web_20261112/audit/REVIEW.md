# Session285 independent mathematical review

**Status: PASS WITH MINOR SCOPE/TERMINOLOGY REPAIRS.** All new mathematical lanes pass the scoped independent review below. No Chan–Pak completeness theorem is claimed.
This is not a re-audit of every inherited theorem cited in the historical web.

## Scope and method

Read the local stress derivation, discrete finite-difference closure, transport
connection, global period classification and their explicit hypotheses.
Independently derived the real equilibrium/complex equation correspondence,
the real-polyanalytic kernel, and the local/global parameter distinction.
`independent.py` solves real polynomial coefficient constraints rather than
regenerating the author's Hermitian basis; checks the symbolic CR flatness
condition for the five-state system; and computes residues by differentiating
a pole-cleared rational function rather than using the author's residue API.

## Initial findings

- The local five-dimensional result is sound for the stated smooth,
  body-force-free equilibrium and harmonic-trace hypotheses, with a prescribed
  nonzero holomorphic phase. It is not an arbitrary-traction or globally
  single-valued-displacement theorem.
- `q/h` real is the correct unordered principal-axis condition for the stated
  sign convention. The fifth pressure coordinate is physical, not Airy gauge.
- The Hermitian carrier and Lorentz determinant are correct. The coordinate
  action is deliberately not promoted to full elasticity covariance.
- The all-order real-polyanalytic proof and rectangular discrete closure are
  valid. Polynomial controls at m=1..4 give dimensions1,4,9,16. These controls
  support the proof and are not the source of its universal quantifier.
- The global dimension is exactly5−rank(period map), with at most three real
  constraints. All ranks0..3 are independently reproduced; the last period
  matrix has determinant72. Zeroes of h are excluded as required, and global
  injectivity of f is not assumed.
- Monodromy matrices have square-zero deviations. Suggested terminology repair:
  describe the group as abelian (class≤1), not ambiguously 'two-step nilpotent'.
  The formulas themselves were correct.

## Independent artifacts

- `independent.py` / `independent.json`: exact independent algebraic controls.
- Both normal and optimized Python executions agree.

## Further independent review

- The stellohedron marked-tournament realization retains its labels and
  signatures. Its all-order canonical-path proof correctly generates every
  loop from pentagons and append/swap squares. The stronger elementary-collapse
  proof gives an integral cycle-lattice basis, not merely real spanning.
- `loop_rank.py` independently builds cell boundary matrices at k=1..4; their
  ranks equal cycle ranks0,1,9,66. Normal and optimized outputs agree. The
  generic scalar-flat-state dimension remains v_k, not uniformly five.
- Direct independent subset-path DP (`tournament_counts.py`) verifies the
  rank-five determinant24, the C3/C5/Paley7/carousel7 fiber examples, and
  nonexistence of a differing-H central switch for every labelled deletion
  tournament of order1..4. This is the complete declared lower-order universe.
- Tournament real-fiber annihilator and restricted curvature calculations
  are algebraically correct. Determinants32,399,2646,8988 and negative
  witness−2 were independently recomputed. The midpoint wording needed a
  direction-specific qualifier, not a change to the theorem.
- Stress two-net superposition and the Pythagorean Hermitian slice are sound,
  with their explicitly retained degenerate-net and primitive-triple boundaries.
- The synthesis correctly labels cross-domain research prescriptions as
  proposals; no equality of numerals is used as a theorem.

## Final-source review

- The tournament affine-hull sharpening is valid: for every deletion order
  at least two some cross-cut weight is positive, so the unique multilinear
  Boolean polynomial for H has a genuinely quadratic coefficient. It cannot
  lie in the span of the four affine responses and the constant function.
  Thus the stated dimension is exactly rank(u,v,a,b)+1, not just an upper bound.
- Independent direct subset DP verifies the cyclic order-nine deletion core
  and cuts21,11 give H14207 and14209 with identical retained packet
  (h,Q,U,V,d,e)=(3357,50004,18906,34455,2238,1119).
  No target anchoring or repeatable interval-growth conclusion follows from
  this one pair. The independently replayed cyclic scan agrees with the
  author's exact declared presentation universe; it is not a tournament census.
- The final stellohedron collapse proof associates one elementary cell to
  each non-tree edge and eliminates them in descending top length. Its
  integral cycle-basis and all-order claims have an actual combinatorial
  argument, not extrapolation from the small rank controls.
- History reconciliation independently checks the pinned path/blob manifest,
  aggregate counts, omissions, catalog hashes and shallow-history scope.
  The34,864 read paths,7,509 seed-hit paths and10,069 expansion-only paths
  are lexical retrieval metrics. This review does not claim to have read
  them all or to have independently repeated the entire lexical extraction.
- Final synthesis and research proposals distinguish exact mechanisms,
  inherited results, conditional bridges and open questions. No mathematical
  map between domains is inferred solely from a shared cardinality of five.

## Repairs incorporated

1. Monodromy group is abelian (nilpotency class at most one); its deviation
   algebra has zero products. This avoids a misleading class-two label.
2. Midpoint stationarity in the tournament negative-curvature example is
   asserted along the displayed symmetric witness direction, not every
   response-fiber direction.
3. The concluding missing Chan–Pak step is a universal absolute-target
   section or interval-growth mechanism sufficient for coverage. Existing
   finite renewable interval constructions are not declared nonexistent.
4. History prose must agree with the frozen4,541-subject scope, rather than
   an availability-dependent fresh count. The lexical archive is not the
   complete historical corpus.

## Reproduction and limits

Each independent audit script has been run in normal and optimized Python
modes with identical JSON output. Producer stress, period, holonomy and
specified tournament computations were also replayed against their receipts.
The independent scripts use explicit failure checks, so optimization does not
remove them. Source hashes identify the reviewed snapshot separately.

This review does not re-audit all inherited corpus theorems, verify the original
1933 source line by line, or establish novel priority for every rediscovered
idea. The stress theorem has a self-contained proof under explicit hypotheses.
No continuum analogy supplies a discrete realization theorem. All universal
Chan–Pak, complete tournament-classification and target-selector obligations
remain open unless specifically proved in the displayed scoped theorem.

## Final wrapper review

`replay.py` invokes ten normal/optimized verifier pairs, rejects nonzero exits
and byte-level receipt differences, and verifies archived catalog hashes.
Its documentation separates the optional full lexical replay from the
wrapper, and records the pinned Git tree dependency of manifest reconciliation.
The wrapper itself is not an additional mathematical proof. Final execution
status is recorded by the parent in `REPLAY.json`.
