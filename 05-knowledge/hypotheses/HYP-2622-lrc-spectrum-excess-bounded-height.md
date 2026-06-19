---
id: HYP-2622
title: LRC spectrum excess ledger and bounded-height lower-bound filter
status: PARTIALLY-CONFIRMED
source: codex-2026-06-19-S17
depends_on:
  - HYP-2621
  - THM-524
related:
  - HYP-2052
  - OPEN-Q-108
---

# HYP-2622 - LRC Spectrum Excess Ledger and Bounded-Height Lower-Bound Filter

## Claim

Let

```text
M(S) = max_t min_{v in S} ||v t||
```

for a primitive `k`-speed set, and let the AP floor be `1/(k+1)`.
If `M(S)=p/q` in lowest terms, define the AP-floor excess

```text
e(k,S) = p(k+1) - q.
```

Then

```text
M(S) - 1/(k+1) = e(k,S) / (q(k+1)).
```

The dangerous second-spectrum rows are therefore not merely rows below the
doubled-top value `2/(2k+1)`.  They are rows with both:

```text
small excess e(k,S), especially e=1,
large denominator q relative to k.
```

The KPS S9 denominator lemma gives the complementary filter: exact `M(S)=p/q`
is attained at a binding-pair denominator with

```text
q | (v_i +/- v_j), hence q <= 2 max(S).
```

Thus every non-tight row satisfies

```text
M(S) - 1/(k+1) >= 1 / (2 max(S)(k+1)).
```

Consequently an `o(1/k^2)` dip in the true second point `sigma_2(k)` can only
come from extremizers with `max(S)/k -> infinity`.  Every bounded-height or
bounded-multiplier family has a built-in `Theta(1/k^2)` lower floor.

## Computation

Script:

- `04-computation/lrc_spectrum_ap_defect_excess_codex_s17.py`
- output: `05-knowledge/results/lrc_spectrum_ap_defect_excess_codex_s17.out`

This continues HYP-2621 on the one-defect family

```text
A_{k,r} = {1,...,k} \ {k-1} union {r(k-1)}.
```

It records:

```text
excess e = p(k+1)-q,
reciprocal depth = q(k+1)/e,
normalized depth = q/(e(k+1)).
```

## Bounded Multiplier Audit

The exact grid `k<=36`, `r<=12` found state counts:

```text
unit-other: 33
positive: 308
unit-target: 16
AP-tight: 6
```

The top normalized reciprocal depths are fixed-multiplier unit-excess rows:

```text
k=31, r=4, M=4/127, e=1, normalized depth=127/32
k=25, r=3, M=3/77,  e=1, normalized depth=77/26
k=19, r=3, M=3/59,  e=1, normalized depth=59/20
k=13, r=3, M=3/41,  e=1, normalized depth=41/14
k=7,  r=3, M=3/23,  e=1, normalized depth=23/8
```

The doubled-top branch `r=2` supplies its known unit-excess rows, but with
smaller normalized constants than the `r=3` and `r=4` branches on their
residue classes.

No `r>=5` row in the stored grid improves the normalized depth constant.
In the high-`r` probe on the `k == 1 mod 30` branch where `r=3` is AP-tight:

```text
k=31: r=4 has e=1, but r=5,6 have e=3,5; r=7 gives M=1/31; r>=8 has larger excess.
k=61: r=4 has e=1, but r=5,6 have e=3,5; r=7 gives M=1/61; r>=8 has larger excess.
```

So this continuation gives no evidence for a growing multiplier constant
inside the AP-defect ladder.

## Symbolic Witness Seeds

For the `r=3` branch, HYP-2621 observed

```text
M(A_{k,3}) = 3/(3k+2)
```

on tested residue classes `k == 7,13,19,25 mod 30`.  S17 extracts explicit
witness numerators.  With `q=3k+2`, the following `t=a(k)/q` clears every
speed in `A_{k,3}` to distance at least `3/q` in the tested range:

```text
k == 7  mod 30: a(k) = (3k - 1)/5
k == 13 mod 30: a(k) = (6k + 7)/5
k == 19 mod 30: a(k) = (6k + 1)/5
k == 25 mod 30: a(k) = (3k + 5)/5
```

This proves the witness half of the third-mediant branch symbolically.  The
remaining proof half is the upper certificate `M(A_{k,3}) <= 3/(3k+2)`, i.e. a
modular cover or no-better-crossing argument.

## Integration With KPS S9

The concurrent KPS S9 scripts add two critical context points:

- Exhaustive primitive boxes through `k=7` show the true second-spectrum
  picture is more sporadic than the doubled-top mediant.  In particular,
  `k=6` has `sigma_2=5/33 < 2/13`, and `k=7` has `sigma_2=3/23 < 2/15`.
- The denominator bound `q <= 2 max(S)` means that any unboundedly smaller
  lower-bound scale must force `max(S)/k -> infinity`.

Combined readout:

```text
fixed AP-defect ladders improve constants,
sporadic sets can beat the doubled-top mediant at small k,
but neither phenomenon threatens Theta(1/k^2) unless height grows with k.
```

Thus the live lower-bound problem should be restated as a height-escape
problem:

```text
Can true sigma_2 extremizers with max(S)/k -> infinity keep excess small enough
to make g(k) = sigma_2(k)-1/(k+1) = o(1/k^2)?
```

## Tournament Analysis

The session explicitly rejects raw runners as the default tournament vertices.
The quotient preserving the LRC-spectrum predicate is the family/excess ledger:

```text
numerator excess,
denominator depth,
height ratio max(S)/k,
residue class,
witness formula,
upper-cover obligation.
```

The proof-route tournament is transitive:

```text
infinite_lower_bound_for_g(k)
> excess_one_classifier
> r3_symbolic_witness
> high_r_excess_audit
> bounded_multiplier_grid
> doubled_top_packet
> raw_runner_vertices
```

Score histogram: `{0:1,1:1,2:1,3:1,4:1,5:1,6:1}`.  Directed `3`-cycles: `0`.

The challenged assumption: tournament vertices need not be runners, arcs, or
pair crossings.  Here the useful vertices are proof obligations/families,
because they preserve the gap scale and the exact lower-bound threat while
discarding irrelevant witness-time geometry.

## Status

Partially confirmed by exact rational computation and by integration with the
KPS denominator lemma.  HYP-2622 does not prove the global `Theta(1/k^2)` lower
bound for `g(k)`.  It reduces the live obstruction to height escape:

```text
bounded height => Theta(1/k^2) lower floor;
o(1/k^2) requires max(S)/k -> infinity plus small excess.
```

Next steps:

```text
1. Prove the upper cover for the r=3 residue classes.
2. Build a search organized by excess and height ratio, not by raw speed boxes.
3. Hunt specifically for small-excess rows with max(S)/k growing.
4. Relate those height-escape rows, if any, to the LRC(14) support-six coimage tails.
```
