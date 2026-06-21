---
id: HYP-2722
title: LRC14 miss-zeta word compatibility filters out cheap abstract q0-hiding atom moves
status: OPEN; exact frontier scout
source: codex-2026-06-21-S71
depends_on:
  - HYP-2721
  - HYP-2720
  - HYP-2719
  - HYP-2698
  - HYP-2702
  - THM-558
  - THM-561
related:
  - HYP-2697
  - HYP-2696
  - THM-557
  - THM-556
  - OPEN-Q-108
---

# HYP-2722 - Miss-Zeta Word Compatibility Filter

## Claim

The cheap abstract atom-cone moves from HYP-2721 are not proof-relevant until
they pass a generated miss-zeta word compatibility test.

At fixed slow coordinate `x`, HYP-2698 says an actual decorrelated context is
not an arbitrary residual vector.  It is a product word in miss-zeta
coordinates:

```text
Z_context,x(A) = product_i Z_i,x(A),
Z_x(A) = Pr(A subset residual mask).
```

HYP-2722 proposes the next filter:

```text
abstract atom move
  -> generated miss-zeta product word
  -> missed-count atom delta
  -> q0 boundary evaluation.
```

The abstract LP in HYP-2721 found low-moment-hidden moves with small
non-origin tax.  The HYP-2698/HYP-2702 generated sparse-tail frontier does not
produce those cheap hidden directions.  Generated words force either positive
`q0` margin, visible low-factorial leakage, or positive Bonferroni4 readout.

## Exact Scout

Script:

```text
04-computation/lrc14_miss_zeta_word_compatibility_codex_s71.py
```

Stored output:

```text
05-knowledge/results/lrc14_miss_zeta_word_compatibility_codex_s71.out
```

The scout takes every sparse-coordinate challenger shape from the
HYP-2702 frontier and every coherent generated context of complementary total
size.  For consecutive block `K` and challenger `C`, it computes the final
missed-count atom move

```text
q_t = Pr_K(T=t) - Pr_C(T=t)
```

after context and candidate cluster have been OR-convolved at the same slow
coordinate `x` and integrated exactly.

It then normalizes by `q_0` and compares the generated move to the cheap
abstract LP directions from HYP-2721.

## Evidence

The full frontier audit has:

```text
tests                 = 318
q0 failures           = 0
global min q0         = 20/16807
global min atom tax   = 136551/61661 ~= 2.214544039
global min W1 leak    = 0
global min W1+W2 leak = 144154/63487 ~= 2.270606581
global min W1+W2+W3   = 284898/63487 ~= 4.487501378
U4/q0 nonpositive     = 0
tail45/q0 negative    = 0
```

So `W_1` alone can be silent in a generated row, but `W_1+W_2` is not silent
on this frontier.  This matters because the cheapest abstract `r=1` and `r=2`
LP directions rely on low-factorial invisibility.  The generated word keeps
enough low-order leakage to expose them.

The closest cheap `r=1` direction is still not close:

```text
size=5, shape=(0,1,2,3,5), context=[2]
q0                  = 61661/2469600
nonorigin tax       = 136551/61661 ~= 2.2145
W1+W2 leak          = 163599/61661 ~= 2.6532
U4/q0               = 69930/61661 ~= 1.1341
L1 distance to LP r=1 direction = 918112/308305 ~= 2.9779
```

The closest cheap `r=2` direction is likewise exposed:

```text
size=4, shape=(0,1,2,4), context=[3]
q0                  = 2479/98784
nonorigin tax       = 7533/2479 ~= 3.0387
W1+W2 leak          = 10274/2479 ~= 4.1444
U4/q0               = 2865/2479 ~= 1.1557
L1 distance to LP r=2 direction = 9682/2479 ~= 3.9056
```

The smallest q0 margin is the already-known HYP-2698 singleton frontier:

```text
size=3, shape=(0,1,3), context=[1+1+1+1]
q0 = 20/16807.
```

Its normalized atom profile is far from cheap:

```text
(1, 4, -3/2, -109/30, 1/40, 3/28, 1/840).
```

## Proof Program

The compatibility proof should not try to characterize all signed atom moves.
It should prove a generated-word exclusion principle:

```text
If a q0 move comes from an HYP-2698 miss-zeta product word, then either
  (a) q0 has positive generated-context margin,
  (b) W1+W2 leakage is quantitatively visible,
  (c) Bonferroni4/tail45 is positive enough for THM-558/HYP-2696, or
  (d) the packet is low-support/low-height and goes to the finite ledger.
```

Concrete next lemmas:

1. **Singleton compatibility lemma.**  Prove the all-singleton generated word
   excludes the cheap LP directions by a signed death-chain inequality in the
   HYP-2702 kernel `g_r(t)`.
2. **Context merge monotonicity.**  Prove merging singleton context carriers
   into coherent blocks cannot reduce the relevant `W1+W2` leakage or q0
   margin below the singleton frontier.
3. **Low-support relation handoff.**  Combine with HYP-2719: support-3 Schur
   packets that survive the compatibility filter are finite-ledger objects,
   not high-tail analytic objects.
4. **Bonferroni split.**  Use THM-558 for the `U4`-visible cheap moves, but do
   not rely on `U4` alone because HYP-2721 found cheap abstract moves with
   `U4_delta=0`.

## Tournament Analysis

Vertices: coherent generated context partitions.

Pairwise observable:

```text
larger minimum W1+W2 leakage,
then larger minimum non-origin atom tax,
then larger minimum q0 margin.
```

Switch/gauge:

```text
arbitrary residual coordinate -> generated miss-zeta product word
  -> missed-count atom delta -> low-factorial leakage.
```

Fingerprint:

```text
vertices = 11
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1}
directed_3cycles = 0
```

Pressure path:

```text
[1+1+1+1]
> [1]
> [1+1+1]
> [3+1]
> [2+1+1]
> [2+2]
> [4]
> [1+1]
> [3]
> [2+1]
> [2].
```

This ranks compatibility barrier strength, not danger.  The weakest barrier in
the frontier is context `[2]`, but it still has `W1+W2` leakage
`144154/63487`.

## Assumption Challenge

This session considered vertices as residual masks, sparse-tail challenger
shapes, context partitions, missed-count atoms, low-factorial observers,
Bonferroni4 readouts, LP rays, and proof obligations.

The chosen quotient is the generated miss-zeta product word.  It preserves the
fixed-`x` compatibility structure HYP-2698 needs, then evaluates the atom move
only after OR/deletion composition.  It destroys arbitrary cone directions;
that is intentional, because HYP-2721 showed the arbitrary cone is too loose.

Challenged assumption: a cheap signed atom move in the abstract six-sector
cone corresponds to a real LRC packet.  The S71 frontier says no: at least on
the generated sparse-tail testbed, cheap moves are exposed by low factorial
leakage and positive `U4/q0`.

## Status

This is not an LRC14 proof.  It is a compatibility filter and a proof program
for the HYP-2721 q0/Vitali cone route.  The next proof step is to turn the
observed `W1+W2` leakage floor for generated words into a symbolic
singleton-context inequality, then prove coherent context merging.
