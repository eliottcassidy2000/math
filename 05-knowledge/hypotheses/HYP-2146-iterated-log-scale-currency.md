---
id: HYP-2146
status: CONCEPTUAL PROOF-PROGRAM; S600 gives an inequality atlas, not a theorem
source: codex-2026-06-03-S600
related:
  - HYP-2064
  - HYP-2076
  - HYP-2145
  - HYP-2144
  - THM-369
  - THM-398
---

# HYP-2146: iterated-log inequalities are scale-currency ledgers

## Abstraction

The deeper pattern behind many `log log` and `log log log` inequalities is not
that those functions are magical.  They are accounting records for proof
currencies paid across nested scale ladders.

The common skeleton is

```text
residual_survival <= product_j (1 - saving_j)
                  <= exp(-sum_j saving_j).
```

Then the logarithm that appears is determined by the sum of the savings:

```text
sum 1/j                         ~ log J
sum 1/(j log j)                 ~ log log J
sum 1/(j log j log log j)       ~ log log log J
```

So an iterated log is a diagnostic.  It says how many times the proof has
compressed its scale space before finding averageable mass.  Ordinary dyadic
scales give `log`; prime harmonic mass, denominator tiers, or compressed scale
blocks often give `log log`; meta-blocks of those blocks give `log log log`.

## Relation To HYP-2145

Incoming mainline work in HYP-2145 frames iterated logs as the inverse
hyperoperation tower: `log`, `log log`, `log log log`, and `log*` measure how
many nested exponential/product hierarchies have been peeled away. This
hypothesis is the complementary proof ledger. It asks, once the relevant tower
level is identified, whether that level is being paid as a tax, earned as an
independent dividend, or routed through a local certificate family.

In short: HYP-2145 locates the height of the hierarchy; HYP-2146 tracks the
currency flow along that hierarchy.

## Relation To Tao-Style LRC Bounds

In the LRC lower-bound thread, the trivial `1/(2k)` bound is a first-moment
statement.  Tao's improvement has the shape

```text
1/(2k) + c log k / (k^2 (log log k)^2).
```

In this language, the surplus is an overlap dividend: a second-moment or
Fourier/Riesz-type argument finds overlap mass, but that mass is spread across
a compressed scale space and pays two `log log` taxes.  The repo's sieve work
does not beat this as a universal theorem; it often exits the Tao regime
entirely by finding a full `1/(k+1)` witness in a structured case.

## S600 Atlas

`04-computation/iterated_log_scale_laws_s600.py` records the scale laws and a
small numerical comparison.  The table includes:

```text
Tao_shape_unit      = log k / (k^2 (loglog k)^2)
third_log_tax       = Tao_shape / logloglog k
third_log_dividend  = Tao_shape * logloglog k
rank_tax_shape      = log k / (k^2 (loglog k + r logloglog k)^2)
```

The values are only orientation constants.  The important object is the
currency ledger: which scale classes are being averaged over, which are
unresolved, and whether the proof earns or pays one more meta-scale factor.

## New Inequality Templates

These are proposed proof templates, not established LRC theorems.

### 1. Scale-Harmonic Product Lemma

If a residual process has one scale-saving per level,

```text
R_{j+1} <= R_j (1 - alpha/j),
```

then

```text
R_J <= R_m (m/J)^{c alpha}
```

for an absolute `c` once the savings are small enough.  This is the ordinary
log-scale product inequality.

### 2. Compressed-Scale Product Lemma

If the proof only earns savings after one compression,

```text
R_{j+1} <= R_j (1 - alpha/(j log j)),
```

then

```text
R_J <= R_m (log m / log J)^{c alpha}.
```

When `J` itself is a logarithmic cutoff, this becomes a `log log` inequality.
This is the abstraction behind "the proof found mass in harmonic scale blocks,
not in every raw scale."

### 3. Meta-Compressed Product Lemma

If the proof has to average over scale-blocks of scale-blocks,

```text
R_{j+1} <= R_j (1 - alpha/(j log j log log j)),
```

then

```text
R_J <= R_m (log log m / log log J)^{c alpha}.
```

With `J ~ log N`, this is where `log log log N` enters.

### 4. Meta-Scale Dividend Template

If a Tao-style second-moment proof finds one profitable block in each of
`~log log log k` independent meta-scale classes, try for a surplus of the form

```text
Delta(k) >= c log k logloglog k / (k^2 (loglog k)^2).
```

This is a dividend version: the third log appears upstairs because meta-scale
classes produce independent overlap certificates rather than extra losses.

### 5. Rank-Tax Template

If the residual has `r` orthogonal scale channels that must all be paid, replace
the Tao denominator by

```text
(loglog k + r logloglog k)^2
```

and aim for

```text
Delta_r(k) >= c log k / (k^2 (loglog k + r logloglog k)^2).
```

This reads `log log log` as a rank tax: each unresolved channel consumes one
meta-scale layer.  It is a natural way to describe determinant/CRT residuals
where several independent prime-power lanes must align.

### 6. Helly-Scale Template

S599/HYP-2144 suggests a different place for iterated logs to appear.  If every
live bounded automaton state of size `M` has a determinant witness among the
first `H(M)` rows, then a global CRT search can be replaced by a sum over
minimal Helly witnesses.  If the witness search has harmonic mass

```text
sum_{h <= H(M)} 1/h ~ log H(M),
```

then the log factor comes from certificate size, not from time or denominator
height.  If `H(M)` itself is logarithmic in the ambient scale, the proof pays a
`log log` or `log log log` tax automatically.

S601/HYP-2152 makes the accounting explicit.  For `M` live determinant
component languages and Helly cutoff `H`, the local certificate entropy is

```text
Lambda_H(M) = log sum_{h<=H} binom(M,h).
```

For bounded `H`, this is `H log M + O_H(1)`.  Thus a Helly logarithm is paid in
the component-count variable, not automatically in denominator height.  The
S601 sample found only singleton and pair certificates in the S599 two-block
determinant residual, so the present n=14 branch is in the Helly-dividend
regime rather than a high-rank Helly-tax regime.

## LRC Payoff

For the LRC n=14 program, this says:

1. The bulk exits by hard certificates such as THM-369, cap overloads, endpoint
   owners, and determinant Helly witnesses.  These are not Tao-regime
   inequalities; they are full-bound exits.
2. The Tao/Bedert style lower-bound surplus matters only on the residual core
   that survives those exits.
3. The right question on that core is not "can we make logs prettier?" but
   "which scale currency is still unpaid?"

For the current determinant branch, the likely unpaid currencies are:

```text
prime-power CRT lanes,
minimal determinant Helly size,
owner-rank of the residual,
and meta-scale blocks of denominator witnesses.
```

## Tournament Analysis

Vertices should be scale currencies rather than runners:

```text
raw scale, compressed scale, meta-compressed scale,
prime harmonic block, determinant Helly size, residual rank, proof obligation
```

Pair observable:

```text
(saving exponent, compression depth, unresolved rank, certificate locality)
```

Switch/gauge: a lower compression depth beats a higher one; ties prefer larger
saving exponent and smaller residual rank.  The expected fingerprint is
transitive until two currencies trade off, at which point directed cycles mark
genuine proof-routing ambiguity.

## Assumption Challenge

Candidate vertices considered: integers, primes, denominator scales, prime
blocks, proof obligations, determinant rows, residual packets, Fourier modes,
and entropy channels.  The chosen quotient is scale currencies.

Preserved predicate:

```text
where the proof pays or earns logarithmic factors.
```

Destroyed information: the original labels of runners, components, and
denominators.  Challenged assumption: an iterated log is a nuisance term to be
simplified.  In this view it is a coordinate telling us which scale ladder the
proof is actually using.

## Honest Status

This is a conceptual proof-program.  The product inequalities are standard
scale calculus; the LRC-specific dividend/rank/Helly templates are proposals.
They should be tested against generated residuals from HYP-2076 and HYP-2144
before being promoted to a theorem.

**See:** `04-computation/iterated_log_scale_laws_s600.py`
(+ `05-knowledge/results/iterated_log_scale_laws_s600.out`),
`07-reflections/iterated-log-scale-currency-s600.md`, HYP-2064, HYP-2076,
HYP-2144, HYP-2145, THM-369, THM-398.
