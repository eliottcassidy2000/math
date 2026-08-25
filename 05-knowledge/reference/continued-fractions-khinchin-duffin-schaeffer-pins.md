# Continued fractions, Khinchin, and Duffin--Schaeffer primary pins

**Reference sidecar, 2026-08-24.** This file expands the bounded entry in
[CORE-PAPERS](CORE-PAPERS.md#continued-fractions-and-khinchin-content).

## Khinchin digit statistics

Cellarosi--Hensley--Miller--Wellens,
[arXiv:1402.0208v3](https://arxiv.org/abs/1402.0208), pins the classical
almost-everywhere digit-mean theorem.  The exact
[continuant probe](../../04-computation/jc_lrc_khinchin_continuant_sidecar_probe_20260823.py)
and Pell control
[THM-3744](../../01-canon/theorems/THM-3744-pell-prefix-loneliness-constant-carry-exact-formula.md)
show that this scalar recovers neither ordered words, target sidecars, nor
flatness.  It supplies no irrationality result for Khinchin's constant.

## Koukoulopoulos--Maynard — *On the Duffin--Schaeffer conjecture*

- **Primary / freshness:** [arXiv:1907.04593v3](https://arxiv.org/abs/1907.04593),
  published in [*Annals of Mathematics* 192 (2020), 251--307](https://annals.math.princeton.edu/2020/192-1/p05).
  **PUBLISHED / stable.**
- **Imported role:** for an arbitrary nonnegative approximation function
  `psi`, the set of real numbers admitting infinitely many reduced
  approximants `a/q` with `|alpha-a/q|<=psi(q)/q` has full Lebesgue measure
  when `sum_q psi(q)phi(q)/q` diverges.  The proof controls correlations
  between denominator layers through a second-moment argument and GCD-graph
  structure.
- **Repo consumer:**
  [THM-4056](../../01-canon/theorems/THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock.md)
  uses only the theorem's exact reduced-centre and first-moment input.  Its
  finite compiler identifies each divisor-complete family of primitive
  centres with one labelled LCM clock.
- **Does not prove:** THM-4056's finite cyclic bijection or Fourier formula,
  equality between formal multiplicity mass and finite union measure, any
  irrationality statement about a named constant, or LRC(14).  Conversely,
  THM-4056 does not reproduce the paper's infinite correlation estimate.
