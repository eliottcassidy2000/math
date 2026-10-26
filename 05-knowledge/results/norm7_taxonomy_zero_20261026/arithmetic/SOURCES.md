# Source receipts for the quadratic exceptional-fibre atlas

Checked 2026-10-26. Primary text was opened and the specific statements below
were read, including hypotheses. The Springer chapter HTML is the published
2021 edition; no claim relies on an inaccessible search snippet. The larger
author-hosted PDF exceeded the web tool's size limit and was not treated as
read. No novelty is claimed for quaternion arithmetic or conic splitting.

## Primary classical inputs

- John Voight, *Quaternion Algebras*, GTM288, Springer2021,
  [Chapter5](https://link.springer.com/chapter/10.1007/978-3-030-56694-4_5),
  Theorem5.5.3 and the preceding pure-norm equation. Used only to identify
  rational points of the displayed conic with splitting of its quaternion
  algebra. The packet supplies the actual coordinate transformation.
- Same book,
  [Chapter12](https://link.springer.com/chapter/10.1007/978-3-030-56694-4_12),
  §12.4.12–16. Used for Hilbert symbols at odd primes and at two, including
  (-1,2)_2=1 and the odd-unit formula. The signed ramification list is the
  specialization to (-1,delta); direct primitive-congruence controls check
  it independently on a finite universe.
- Same book,
  [Chapter14](https://link.springer.com/chapter/10.1007/978-3-030-56694-4_14),
  Proposition14.6.7, Corollary14.6.5, and the discussion immediately following
  Proposition14.6.7. Used for global splitting and the criterion that no
  quaternion-ramified place split in the quadratic extension. This is the
  substantive local-global input, and is CITED rather than re-proved here.

## Inherited repository results

- `normseven_descent_moments_20261026/RESOLVED_FIBRES.md`, R1–R4:
  the exact family of rational quotient points, quadratic node fields,
  exceptional conics and the example degrees2/4. R5–R6 supply the separate
  good-prime zeta comparison. The current packet extends the degree decision
  to every rational pair (D,delta); it does not re-prove those cover equations.
- `normseven_arithmetic_f1_20261026/boundary/RESULT.md`, C2–C3:
  the rational-node local equation, Gaussian norm criterion over Q, and
  relative split-torus boundary. The present quadratic-field decision
  genuinely needs more than that Q-only prime-factor criterion.
- `normseven_arithmetic_f1_20261026/gale/RESULT.md`, TheoremGA:
  constant outer K3 twists and their Neron–Severi Galois character. Mentioned
  for comparison of retained arithmetic invariants; it is not a dependency
  of the conic theorem and does not give a nonconstant-twist classification.
- `01-canon/MISTAKES.md`, MISTAKE-506: the Hamilton algebra ramifies exactly
  at2 and infinity. Included as a hostile against treating denominator
  primes as quaternion ramification.

The numerical receipt and source hashes are recorded in `MANIFEST.json`.
All output is finite-exact. The all-rational-parameter conclusion uses the
written proof and the cited global theorem, not extrapolation from the grid.

## Conductor-support proof inputs, checked 2026-10-26

- Finn Knudsen, *The projectivity of the moduli space of stable curves, II*,
  Math. Scand.52 (1983), printed p179, Theorem2.7 and its proof:
  [author-hosted collected PDF](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1976b--DetDiv%2BKnudsenSequel-GT.pdf).
  This gives the stable-curve stack smooth proper over Spec Z with relative
  normal-crossing boundary. Its proof records the universal-curve induction
  used here to obtain schemes in genus zero. The exact printed page was
  visually read from PDF page index57 after local rendering, since web
  screenshot results returned text references without image pixels. The
  genus-zero labelled-tree description of intersections supplies the
  globally simple boundary required in the new deduction.
- [Stacks Proposition63.15.2](https://stacks.math.columbia.edu/tag/0GLG):
  finite locally constant higher direct images under a smooth proper map,
  for torsion coefficients of invertible order. Exact hypotheses read.
- [Stacks Lemma18.43.5](https://stacks.math.columbia.edu/tag/093U):
  finite locally constant module sheaves form a weak Serre subcategory.
  Used on the integral base, not on generic-point representations alone.
- [Stacks Theorem59.91.11 and Lemma59.91.12](https://stacks.math.columbia.edu/tag/095S):
  proper base change for torsion sheaves and bounded-below complexes.
  Applied to extension by zero on proper compactifications.

These sources do not identify the four proposed pure graph constituents or
their conductor exponents. The graph-cover construction and support deduction
are written explicitly in CONDUCTOR_SUPPORT.md and require their own audit.
