---
id: HYP-2721
title: LRC14 q0 origin atom is the Vitali atom-cone boundary functional
status: OPEN; exact scout plus proof program
source: codex-2026-06-21-S70
depends_on:
  - HYP-2720
  - HYP-2719
  - HYP-2718
  - THM-561
  - THM-170
  - THM-169
related:
  - HYP-2157
  - HYP-2104
  - THM-406
  - THM-558
  - HYP-2698
  - OPEN-Q-108
---

# HYP-2721 - Q0 Origin / Vitali Atom Cone

## Claim

The `Q_0` origin atom in the LRC14 carrier-product gap is the LRC analogue of
the tournament Vitali atom only after changing quotient:

```text
raw cover/profile data
  -> miss-zeta factorial moments W_j
  -> missed-count atom discrepancies Q_t
  -> relation-support/generated-word ledger.
```

The analogy is not "Vitali means small measure" and not "odd support is a
signed cone."  The shared structure is:

```text
low observers can be preserved,
while a hidden packet channel still moves the target functional.
```

In tournaments, THM-169/170 shows that a lambda-preserving 4-reversal can keep
the visible low cycle channels fixed while changing `H` through `c7` or
disjoint-pair packets.  In LRC14, THM-561 shows that a high factorial packet
can keep lower factorial moments unchanged while moving the origin atom

```text
Q_0 = ProductCover - ActualCover = sum_j (-1)^j W_j.
```

So the useful "Vitali atom cone" is the finite-difference cone spanned by the
unit factorial moves

```text
B_j(t) = (-1)^(j-t) binom(j,t), 0 <= t <= j.
```

A unit `W_j` creates atom profile `B_j`, preserves every lower factorial moment
`W_i` for `i<j`, and moves the origin atom by `(-1)^j`.  This is the abstract
LRC substitute for a Vitali move.

## Dictionary

```text
tournament Vitali 4-reversal
  <-> factorial finite-difference atom move B_j

lambda/c3/c5 visible data
  <-> low factorial observers W_1,...,W_r and cheap cut-space shadows

c7 / i2 hidden channels
  <-> high W_j tail, Q_t non-origin atom cone, and carrier relation support

endpoint completion matrix
  <-> generated miss-zeta word / sector-ownership context

H-change
  <-> q0 boundary evaluation ProductCover - ActualCover
```

This makes HYP-2719 and HYP-2720 sequential rather than competing.  HYP-2719
selects the Fourier relation-support packets, especially support-3
Schur/additive triangles.  HYP-2720/HYP-2721 then taxes those packets in the
factorial-origin basis before evaluating `Q_0`.

## Exact Scout

Script:

```text
04-computation/lrc14_q0_vitali_atom_cone_codex_s70.py
```

Stored output:

```text
05-knowledge/results/lrc14_q0_vitali_atom_cone_codex_s70.out
```

The scout reuses the exact S69 miss-zeta row bank and reports:

```text
origin_share = |Q_0| / sum_t |Q_t|
tax_ratio    = sum_{t>0}|Q_t| / |Q_0|
tail2_L1     = sum_{j>=3}|W_j|
```

Aggregate signal over the ten-row bank:

```text
sum |Q0|                    = 0.095623745...
sum atom_L1                 = 0.620915350...
aggregate origin share      = 0.154004479...
aggregate nonorigin tax     = 5.493317626...
sum tail2_L1 / sum |Q0|     = 15.926622148...
sum tail3_L1 / sum |Q0|     = 7.188663423...
cone_margin positive rows   = 10/10
```

Thus `Q_0` is a boundary functional on a much larger atom cone.  Trying to
dominate the whole missed-count law is too strong; the origin atom should be
bounded after relation support and generated-word compatibility have already
discarded fake cone directions.

The largest row-level Vitali-blindness pressure is the `3+3+2` split:

```text
tail2_L1 / |Q0| = 129.394025...
nonorigin tax   = 30.895795...
```

The biggest `|Q0|` row is the five-2-block row:

```text
Q0              = -0.026288364...
tail2_L1/|Q0|  = 9.159569...
nonorigin tax  = 5.962277...
```

So there are two different pressures: hidden-cone blindness and cap-risk size.
They should not be collapsed into one scalar.

## Abstract Cone Tax LP

Follow-up script:

```text
04-computation/lrc14_vitali_cone_tax_lp_codex_s70.py
05-knowledge/results/lrc14_vitali_cone_tax_lp_codex_s70.out
```

This solves the exact abstract LP:

```text
normalize q_0=1,
preserve W_1,...,W_r,
minimize sum_{t>0}|q_t|.
```

The local finite-difference move `B_{r+1}` is canonical, but it is not always
the cheapest move in the full six-sector atom cone:

```text
r=0: min tax 1,   B_1 tax 1
r=1: min tax 7/5, B_2 tax 3
r=2: min tax 3,   B_3 tax 7
r=3: min tax 7,   B_4 tax 15
r=4: min tax 19,  B_5 tax 31
r=5: min tax 63,  B_6 tax 63
```

The cheap abstract optimizers use high missed-count states, for example

```text
r=1: q_1=-6/5, q_6=1/5
r=2: q_1=-9/5, q_3=1, q_6=-1/5
```

This is a guardrail.  The atom-cone basis identifies the correct boundary
functional, but the abstract cone is too permissive.  A proof must retain the
generated LRC compatibility constraints, and it must keep the high-tail
Bonferroni/transfer ledgers from THM-558/HYP-2696 in view.  Otherwise high
`t` atoms can fake low-moment preservation too cheaply.

The LP also prints the Bonferroni4 readout

```text
U4_delta = q_0 + q_5 + 5q_6.
```

The cheapest abstract moves at `r=1` and `r=3` are visible to `U4`:

```text
r=1: U4_delta=2
r=3: U4_delta=5/2
```

But the cheapest moves at `r=2,4,5` have `U4_delta=0`.  Therefore THM-558 is
a necessary high-tail ledger but not a complete replacement for generated-word
compatibility and relation-support selection.

## Proof Route

The proof target should be a `Q_0` Vitali-cone lemma:

```text
For every HYP-2714 moderate-balanced split row E,
after removing the finite low-support relation packets,
the generated high-factorial atom cone satisfies
|Q_0(E)| <= cap_k - ProductCover(E).
```

Proposed split:

1. Use HYP-2719 to select low-support carrier packets by relation support,
   especially support-3 additive triangles.
2. Route low-height support-3 and small-denominator packets to the HYP-2714
   finite ledger.
3. For the remaining high-height/high-denominator tail, use the HYP-2720
   odd-support `L1` envelope in factorial coordinates.
4. Preserve generated miss-zeta words and sector ownership until the
   compatibility cone is known to be genuine.
5. Only then evaluate the boundary atom
   `Q_0=sum_j(-1)^j W_j` and compare with `cap_k-ProductCover`.

The Vitali lesson is that a low-observer quotient is allowed only when the
hidden channel has a separate completion ledger.  For tournaments that ledger
is endpoint completion of 7-cycles/disjoint pairs.  For LRC14 it should be
relation support plus generated miss-zeta words.

## Assumption Challenge

This session explicitly considered vertices as:

```text
runners,
carrier blocks,
residual masks,
low factorial moments,
finite-difference moves B_j,
missed-count atoms Q_t,
relation-support packets,
support-3 Schur triples,
generated miss-zeta words,
proof obligations.
```

The productive quotient for this hypothesis is neither runners nor arcs.  It is
the atom-cone boundary basis `B_j -> Q_0`, decorated by relation support.  This
quotient preserves the exact cover correction `Q_0` and the low-observer
blindness mechanism.  It destroys sector ownership and generated-word
compatibility, so those must be restored before using the cone as a proof
object.

Challenged assumption: a positive-cone or stochastic-dominance statement for
the whole missed-count law should prove the cap.  The data says that is
overkill.  The origin atom is only about `15.4%` of aggregate atom `L1` in the
tested bank.

## Tournament Analysis

Vertices: tested split rows.

Pairwise observable:

```text
tail2_L1/|Q0|, then non-origin atom tax, then support-3 Schur proxy.
```

Switch/gauge:

```text
low-factorial observer -> full atom cone -> relation-support proxy.
```

Fingerprint:

```text
score histogram = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed 3-cycles = 0
```

Hamiltonian pressure path:

```text
3+3+2 split
> two 4-blocks, moderate gap
> two 4-blocks, positive Q0 high
> two 4-blocks, positive Q0
> two 4-blocks, ratio 2:1
> seven singleton carriers
> two 4-blocks, high relation phase
> 5+3 split
> five 2-blocks
> two 4-blocks, wider gap.
```

The tournament is transitive in this scout, but the path ranks blindness
pressure, not danger.  The five-2-block row remains the largest absolute `Q_0`
in the bank.

## Status

This is not a proof of LRC(14).  It sharpens the current target: the remaining
moderate-balanced carrier gap should be proved as a Vitali-cone boundary
estimate for `Q_0`, after relation-support selection and before any scalar
absolute-value bound is taken.
