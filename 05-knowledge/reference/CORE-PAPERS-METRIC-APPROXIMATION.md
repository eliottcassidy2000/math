# Metric approximation and continued-fraction source sidecar

> **Freshness:** sources checked 2026-08-24. This is a routed extension of
> `CORE-PAPERS.md`, not an independent truth surface.

## Khinchin digit means

- **Primary:** Cellarosi--Hensley--Miller--Wellens,
  [arXiv:1402.0208v3](https://arxiv.org/abs/1402.0208), records the classical
  theorem that the geometric mean of continued-fraction digits tends to
  Khinchin's constant for Lebesgue-almost-every real, and records Euler's
  atypical continued fraction.
- **Consumers:** the exact
  [continuant probe](../../04-computation/jc_lrc_khinchin_continuant_sidecar_probe_20260823.py),
  [THM-3744](../../01-canon/theorems/THM-3744-pell-prefix-loneliness-constant-carry-exact-formula.md),
  and [THM-4056](../../01-canon/theorems/THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock.md).
- **Does not prove:** that a named number is Khinchin-typical; the arithmetic
  nature of Khinchin's constant; recovery of ordered words, target sidecars,
  recurrence carries, or flatness from a digit mean.

## Duffin--Schaeffer

- **Primary:** Koukoulopoulos--Maynard,
  [arXiv:1907.04593](https://arxiv.org/abs/1907.04593), published in *Annals of
  Mathematics* **192** (2020), proves that in the normalization
  `|x-a/q|<=psi(q)/q`, divergence of `sum phi(q)psi(q)/q` gives a full-measure
  limsup set of reduced approximants.
- **Quantitative primary:** Koukoulopoulos--Maynard--Yang,
  [arXiv:2404.14628](https://arxiv.org/abs/2404.14628), proves its stated
  quantitative counting theorem with finite mass
  `Psi(Q)=sum 2phi(q)psi(q)/q`.
- **Consumer:** THM-4056 recompiles that finite first moment as two exact LCM
  clocks and separately proves a pointwise irrationality firewall.
- **Does not prove:** interval independence by a clock-mean identity,
  membership of `e+pi` or Khinchin's constant, irrationality of either named
  constant, or a tournament/Berggren transfer.

## Irrationality exponent

- **Primary:** Sondow,
  [*Irrationality Measures, Irrationality Bases, and a Theorem of Jarnik*](https://arxiv.org/abs/math/0406300),
  records the continued-fraction formula
  `mu(x)=2+limsup log(a_(n+1))/log(q_n)`.
- **Consumer:** THM-4056 gives a short proof that existence of a finite
  infinite digit-geometric-mean limit forces `mu(x)=2` and records the exact
  Duffin--Schaeffer normalization consequence.
- **Does not prove:** that a single finite digit mean determines an exponent,
  or that a metric almost-everywhere statement applies to any named point.

## Reciprocal tree normalization

[THM-4057](../../01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md)
proves that reciprocal rationals share a projective Euclidean coefficient
word after deleting only a leading zero. The standard finite Khinchin word
always removes `a_0` and is not reciprocal-invariant; `3/5` versus `5/3` is
the exact hostile. The theorem also separates the Stern--Brocot mirror from
the Calkin--Wilf heap ordinal.
