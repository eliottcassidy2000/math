---
id: HYP-9047
title: "Feige's inequality (now a theorem) as a CRT-block anti-concentration sidecar for coverage budgets"
status: >
  WILDCARD HYPOTHESIS — external-import tool, no in-repo claim. External
  input is CITED-ABSTRACT (fetched 2026-07-28): Nie–Wei announce a proof of
  Feige's 2006 conjecture. Quantization sidecar identified as the load-bearing
  gap; decisive finite tests specified. Nothing here modifies any THM status.
source: klein-2026-07-28-S691
related:
  - THM-2051
  - THM-2100
  - THM-2126
  - THM-1290
  - HYP-9046-signed-idempotent-decomposition-for-the-lrc-wall
external:
  - "Z. Nie, J. Wei, *On Feige's conjecture*, arXiv:2607.24528 (2026): for independent nonnegative X_i with E X_i = 1, P(Σ X_i < n+1) ≥ (n/(n+1))^n ≥ 1/e; sharp at X_i = (n+1)·Bernoulli(1/(n+1))."
  - "U. Feige, *On sums of independent random variables with unbounded variance…* (2006) — original 1/13 bound and conjecture."
---

# HYP-9047 — Feige anti-concentration on coprime clock blocks

## The tool

Bonferroni-style certificates (the B5 layer, THM-2051's BONF5 closure) fight
the union bound's dimension loss with quintic corrections. Feige's
inequality is a *dimension-free* lower bound on staying below `mean + 1` for
sums of independent nonnegative variables — the loss is an absolute `1/e`,
independent of the number of summands. It requires genuine independence,
which runner phases do not have — **except across coprime clocks**: for `t`
uniform mod `N = lcm`, the residues `t mod q` over pairwise coprime `q` are
*exactly* independent (CRT). So group the coverage load by coprime clock
blocks and the independence hypothesis is exact, not asymptotic.

## Proposed mechanism

For a speed row with clock decomposition into pairwise coprime moduli
`q_1, …, q_B`, let `Y_b ≥ 0` be the (suitably normalized, mean-one) coverage
load contributed by block `b` at a uniform random time. Feige gives
`P(Σ_b c_b Y_b < Σ c_b + max c_b) ≥ 1/e`-type floors (after scaling each
block to mean one). If the row's *budget* step (the a-priori resolved-modulus
supply: divisor count + E3 budget + rigidity) currently pays a union-bound or
B-layer tax that grows with the number of blocks, replacing it by the Feige
floor caps that tax at an absolute constant.

**The honest gap (quantization sidecar):** Feige bounds
`P(load < mean + 1)`, not `P(load = 0)`; loneliness needs a *zero*-load (or
below-threshold) time. The sidecar that must be supplied per lane: an integer
quantization showing nonzero block loads are ≥ 1 after normalization, so
`load < mean + 1` with `mean` small forces the load profile into a finite
list that the existing exact machinery (clock closures, THM-2057-style
grids) can discharge. Without that, this is only a soft floor.

## Typed connection

- **Source:** the block-decomposed coverage load of a fixed frontier row.
- **Map:** CRT independence of coprime clocks (exact), then Feige on the
  block sums.
- **Preserved:** the mean budget (linearity of expectation); exactness of
  independence (no mixing loss).
- **Lost:** everything below `mean + 1` resolution; within-block dependence
  (blocks sharing a prime must be merged, growing block loads); the constant
  `1/e` is a probability floor, not a pointwise construction.
- **Sidecar:** the quantization step above, plus a merger bound for
  non-coprime clock families.
- **Cheapest decisive test (finite, exact):** pick two explicit small-height
  rows from the 165-row ledger with ≥ 3 pairwise coprime clocks; compute the
  exact distribution of block loads (they are finite exact objects mod lcm);
  check (a) whether the nonzero-load quantization holds, (b) whether
  `mean + 1` sits below the smallest nonzero total load — if yes, Feige is
  not even needed there (the exact computation closes it); if no, check
  whether the Feige floor beats the row's current B-layer tax. Either outcome
  is informative; record which.
- **Hostile control:** the tight AP at its own ruler (max-not-mean guardrail
  #7 — any averaging argument must be tested there first and is expected to
  fail; this hypothesis is only for *budget* steps, never for the good-period
  existence step itself).

## Guardrails acknowledged

Good-period existence is a MAX statement (guardrail 7; MISTAKE-127/129/130):
this tool is explicitly NOT proposed as a mean-value proof of loneliness.
Its only proposed home is the budget/supply layer of the covering-case
endgame, where means are the native currency. Independence is claimed ONLY
across exactly-coprime clocks (CRT), never heuristically.
