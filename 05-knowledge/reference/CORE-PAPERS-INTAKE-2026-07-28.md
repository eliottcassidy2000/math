# External intake 2026-07-28 — owner puzzle bundle (klein-S691)

All four items below are **CITED-ABSTRACT only** (fetched 2026-07-28; proofs
not audited in-repo). Do not cite as proved dependencies; upgrade to CITED
only after a real read.

## Beke–Goh–Hatami–Jaffe–Naylor — *A characterization of idempotent Schur multipliers* (arXiv:2607.14316)

- **Result (abstract):** every idempotent Schur multiplier is a finite signed
  sum of contractive idempotent Schur multipliers (blocky matrices); boolean
  patterns of Schur norm ≤ γ decompose into `L ≤ 2^{Cγ⁶}` signed blocky pieces.
  Resolves Katavolos–Paulsen (2003); extends Cohen–Host beyond the
  translation-invariant case.
- **Consumer:** HYP-9046 (signed decomposition template for the LRC(14)
  cancellation wall / ONE-lattice pairing).
- **Guardrail:** usefulness hinges on a uniform-in-`V` Schur-norm bound for the
  LRC pattern — unverified; the constant `C` is not explicit for us.

## Nie–Wei — *On Feige's conjecture* (arXiv:2607.24528)

- **Result (abstract):** for independent nonnegative `X_i` with `E X_i = 1`,
  `P(Σ X_i < n+1) ≥ (n/(n+1))^n ≥ 1/e` (sharp; Feige 2006 had `1/13`). Proof
  reportedly developed with LLM assistance; builds on Vlassis–Thomas work on
  Gaffke's conjecture.
- **Consumer:** HYP-9047 (dimension-free budget floor on exactly-independent
  coprime CRT clock blocks; targeted at the THM-2588 clustered-tower residual).
- **Guardrail:** independence is exact only across coprime clocks; the bound
  is a mean-plus-one floor, never a loneliness (maximum) statement.

## Incudini–Mazzola — *Practical advantage beyond the quadratic speedup limit with fully-quantum walks* (arXiv:2607.22818)

- **Result (abstract):** fully-quantum Metropolis walks (quantum proposal and
  acceptance) with roughly sixth-degree total query speedup versus the best
  classical walks for low-temperature Ising Gibbs sampling; full
  fault-tolerant compilation.
- **Consumer:** none yet; methodological note for search-harness design
  (nonlocal proposals; the AK certificate search is a classical annealer).
- **Guardrail:** no in-repo quantum claims; treat as landscape.

## steven-le-thien/vandermonde-snp (GitHub) — *Powers of the Vandermonde determinant are eventually non-SNP* (Lean 4)

- **Result (README):** Lean 4 formalization with two independent verifiers —
  one conditional on Macdonald–Jack polynomial facts, one assumption-free via
  Ward recurrences. SNP here is a saturated-Newton-polytope-flavored
  positivity property of symmetric functions.
- **Consumer:** formalization-pattern exemplar (conditional plus
  assumption-free dual verifier) matching our Lean discipline; candidate
  reading for the owner's `μ₃`-fixed branch / `3+3+1=7` depth-two fragment
  ([decode note](../results/FRAGMENT-DECODE-mu3-depth2-tree-klein-S691.md)).
- **Guardrail:** repository claims unaudited; no dependency.
