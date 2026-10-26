# Independent audit of the boundary lane

2026-10-26. Auditor: parallel `gale_arithmetic` agent, separately from the
boundary producer. **PASS** for C1-C3 of [the result](../boundary/RESULT.md),
subject to the stated fields and distinction of normal surface, divisor
normalization and surface resolution. No major or minor correction found.

* C1: direct projective normalization at each of the three triple ports
  gives exactly the displayed three square tests and eight lifts. A pair of
  passing triple ports forces all three qj square and gives a point of the
  smooth conic, hence infinitely many points over an infinite field. This
  validates the zero/eight/infinite trichotomy, including its exact finite
  criterion. The previously published whole-central dichotomy is false;
  its own fixed-twist example supplies the finite-eight counterexample.
* C2: independently substituting `E=delta*c^2/2`, `x=a^2-E`, `z=1-E`
  gives `b^2=delta*c^2-a^2` and exactly the three unit-root equations in
  the note. Their derivatives at all eight sign triples are nonzero.
  The completed local singularity and projectivized tangent conic follow.
  Every rational point of the total inverse-image divisor maps to one of
  the eight rational nodes; each exceptional conic is the full fiber.
  Thus its rational solubility gives the asserted zero/infinite answer.
  The central divisor's two local branches are conjugate over Q(i), so
  normalizing the divisor is a different operation and cannot supply a
  missed rational strict-transform point.
* C3: positivity and the even valuations at primes3 mod4 are necessary
  by clearing denominators and infinite descent. Gaussian Euclidean
  factorization supplies the converse, including negative valuations.
  The parametrization is multiplication by `(1+it)/(1-it)` and preserves
  the norm exactly. The delta2/delta3 contrast is a valid fixed-central-
  scheme hostile: varying d0 affects the transverse cone while disappearing
  from the central equations.
* The added local split-torification iff also passes: a soluble smooth
  conic is P1 and decomposes as Gm plus two rational points. An insoluble
  conic cannot contain even one nonempty split Q-torus stratum, whose
  identity would give a rational point. This says nothing about torifying
  the entire resolved surface or about other F1 conventions.

The producer checker was read and replayed (184,940 explicit gates).
Producer self-review subsequently forced Fraction conversion before each
fiber ratio, removing an incidental integer-to-float path in its draft;
I reread and replayed that version, with the same counts and gate total.
Its bounded search does not decide the all-square conic cases and says so.
An independent checker [audit_boundary.py](audit_boundary.py) imports no
producer code. It uses the two-bit squareclass group of +/-1,+/-2 and
recovers all4,096 twist rows /12,288 triple-port rows, including histogram
3888 zero,144 finite-eight,64 all-pair-square cases. It additionally
constructs Gaussian-prime norm witnesses by the Euclidean algorithm at
all primes1 mod4 below200, and independently checks the primitive modulo9
obstruction at delta3. These are finite controls; the written proof carries
the unbounded trichotomy and norm criterion.

Reproduce from the repository root:

    python3 -B 05-knowledge/results/normseven_arithmetic_f1_20261026/gale/audit_boundary.py
    python3 -B -O 05-knowledge/results/normseven_arithmetic_f1_20261026/gale/audit_boundary.py

The normal/optimized receipts are `boundary-audit.normal.json` and
`boundary-audit.optimized.json` in this directory. Parent-session manifest
records hashes. Required publication repair: replace the stale complete-
central dichotomy in THM-4773's title/status/source routing as well as its
source proof, preserving the correction lineage in MISTAKES.md.
