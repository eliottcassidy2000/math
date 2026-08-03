# Core papers sidecar: arrangements and tournament games

> Moved verbatim from `CORE-PAPERS.md` on 2026-08-03 (boxeph) to free shared
> startup bytes. Same maintenance rule as the main file: refresh primary
> link/version, name the first consumer, sharpen **does not prove**.

### De Concini--Procesi — *On the geometry of toric arrangements*

- **Primary / freshness:** [arXiv:math/0505351v4](https://arxiv.org/abs/math/0505351),
  revised 2005-06-08.
- **Imported role:** defines a toric arrangement as a finite family of
  codimension-one subtori or their cosets and computes complement cohomology.
  It supplies the standard object against which HYP-8830's LRC terminology is
  audited.
- **Repo consumer:** the
  [corrected phase-height reflection](../../07-reflections/orlik-solomon-is-a-repo-wide-pattern-toric-arrangements-are-the-lrc-lens-boxeph-S209.md).
- **Does not prove:** that the infinite Fourier annihilator
  `{k:k.v=0}` is a layer poset, that a thickened safe set is an ordinary
  complement, or that complement cohomology retains LRC height, wall side,
  owner, sign, or deletion data.

### Moci — *A Tutte polynomial for toric arrangements*

- **Primary / freshness:** [arXiv:0911.4823v5](https://arxiv.org/abs/0911.4823),
  revised 2010-11-09; final journal version in *Transactions of the AMS*.
- **Imported role:** introduces the multiplicity Tutte polynomial, proves
  deletion--restriction and positivity, and obtains characteristic and
  Poincare polynomials of a toric arrangement as specializations.
- **Repo consumer:** corrected HYP-8830 and MISTAKE-224.
- **Does not prove:** that the cutoff statistic `N_R`,
  [THM-1820's sinc-weighted Fourier series](../../01-canon/theorems/THM-1820-lrc-is-a-moment-nullcone-problem-relation-lattice-pairing.md),
  or an LRCMod count is an arithmetic Tutte or
  Mobius specialization. A finite character list and its actual layers must be
  stated first.

### Stanley — *Hyperplane Arrangements, Interval Orders and Trees*

- **Primary / freshness:** [author-hosted survey PDF](https://math.mit.edu/~rstan/papers/nas.pdf),
  version dated 1995-12-01.
- **Imported role:** records the braid characteristic polynomial and its `n!`
  real chambers, and the Shi formulas
  `chi(q)=q(q-n)^(n-1)` and `r=(n+1)^(n-1)`. These are the exact classical
  controls replayed in the corrected S209 computation.
- **Repo consumer:** corrected HYP-8830, MISTAKE-224, and the
  [S209 arrangement reflection](../../07-reflections/orlik-solomon-is-a-repo-wide-pattern-toric-arrangements-are-the-lrc-lens-boxeph-S209.md).
- **Does not prove:** that Shi walls are LRC safety walls, that braid
  cohomology is a per-tournament invariant, or that finite-field arrangement
  counts solve the LRC(14) AP-core extraction problem.

### Fisher--Ryan — *Tournament games and positive tournaments*

- **Primary / freshness:** [DOI 10.1002/jgt.3190190208](https://doi.org/10.1002/jgt.3190190208),
  *Journal of Graph Theory* **19** (1995), 217--236.
- **Imported role:** the tournament game has a unique optimal strategy with
  odd-cardinality support. The repo's elementary parity reading restores the
  essential second coordinate: an even principal tournament payoff block is
  nonsingular modulo two.
- **Repo consumer:** the
  [corrected antisymmetry reflection](../../07-reflections/antisymmetry-is-the-hinge-tori-odd-functions-saddles-and-tournaments-boxeph-S210.md).
- **Does not prove:** that pure optimum is equivalent to transitivity, that
  every intransitive tournament has recurrent replicator dynamics, or that a
  game saddle is a Morse saddle. Pure optimum means a Condorcet winner.
