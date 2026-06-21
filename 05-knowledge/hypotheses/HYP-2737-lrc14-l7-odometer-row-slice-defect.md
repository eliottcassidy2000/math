---
id: HYP-2737
title: LRC14 L7 integer-grid defect collapses to a single two-odometer row-slice inequality
status: SUPERSEDED by HYP-2739 residue closed form; retained as Lean-facing route
source: codex-2026-06-21-S73
depends_on:
  - HYP-2736
  - HYP-2733
  - HYP-2730
  - THM-562
related:
  - HYP-2734
  - HYP-2731
  - HYP-2722
  - HYP-2698
  - THM-546
  - OPEN-Q-108
---

# HYP-2737: L7 Odometer Row-Slice Defect

## Claim

HYP-2736 expresses the L7 torus-line discrepancy as an integer defect

```text
D_{p,q} = sum_ij |49*c_ij - 7*p*q| / (343*p*q).
```

The formal proof target can be made smaller.  The same cell counts are the
period counts of the two-clock word

```text
t -> (floor(t/p), floor(t/q)) mod 7,   0 <= t < 7*p*q.
```

For a fixed first-coordinate row, define

```text
rowdef_i(p,q) = sum_j |7*c_ij - p*q|.
```

The exact scout supports the stronger row-slice statement

```text
rowdef_i(p,q) is independent of i,
rowdef_i(p,q) <= 12*p.
```

Then HYP-2736 follows immediately:

```text
sum_ij |49*c_ij - 7*p*q|
  = 7 * sum_i rowdef_i(p,q)
  = 49 * rowdef_0(p,q)
  <= 588*p.
```

Thus the sharp observed tail `D_{p,q} <= 12/(7q)` is reduced to one
seven-bin row-balance inequality for a generated two-odometer word.

## Exact Evidence

Script:

```text
04-computation/lrc14_l7_odometer_row_slice_codex_s73.py
```

Stored output:

```text
05-knowledge/results/lrc14_l7_odometer_row_slice_codex_s73.out
```

The script cross-checks the odometer counter against the HYP-2736
common-breakpoint counter on the `q<=12` bounded-ratio atlas, then scans all
primitive `1 < p/q <= 43/20` with `q<=300`.

Observed facts:

```text
rowdef is identical across all 7 rows: failures=0
HYP-2733 row-rotation / zero-law checks: failures=0
defect_sum == 49*rowdef: failures=0
rowdef <= 12*p violations: 0
max rowdef/p = 12 at p/q=3/2
max q>=5 rowdef/p = 4 at p/q=11/10
```

The row rotation is the concrete HYP-2733 shift structure.  When `7` does not
divide `q`, the rows are cyclic shifts of row zero by

```text
s = p*q^{-1} mod 7.
```

When `7|p` or `7|q`, the rows are uniform, giving the apex-prime zero law.

## Proof Target

**Update after pull:** HYP-2739 proves a stronger exact residue-only closed
form for this row defect and closes the `rowdef_0 <= 12p` target
combinatorially.  HYP-2737 remains useful as a Lean-facing decomposition of
the same proof, but it is no longer the open analytic target.

Prove the single-row inequality:

```text
For coprime 1 < p/q <= 43/20,
rowdef_0(p,q) = sum_j |7*c_0j - p*q| <= 12*p.
```

The row can be viewed as `q` length-`p` windows sampled along the periodic
seven-bin word of the `q`-clock.  This is a one-dimensional balance problem,
not a 49-cell discrepancy problem.  The equality case `3/2` suggests that the
sharp constant is controlled by the shortest nonuniform Christoffel/Sturmian
row word, while the bounded-ratio condition makes larger denominators rapidly
slack.

## Relation To Earlier Routes

This is a refinement of the already usable L7 closure, not a prerequisite for
it.  Incoming HYP-2730 already supplies `D_{p,q} <= 14/p`, sufficient after
the finite atlas.  HYP-2737 gives a cleaner formalization route for the sharper
`12/(7q)` constant:

```text
continuous torus line
  -> common integer grid
  -> two-odometer generated word
  -> row-slice balance
  -> HYP-2736 integer defect inequality.
```

This also explains why the HYP-2722/HYP-2698 "generated word before scalar"
guardrail was right.  The L7 generated language is much simpler than the
64-mask miss-zeta context cone: it is just a two-clock odometer, and scalar
`D_{p,q}` should be evaluated only after the row-word balance is used.

## Tournament Analysis

Vertices are proof obligations, not runners or arcs:

```text
rowdef<=12p
> row-rotation collapse
> odometer identity
> small-q atlas
> r>=3 pairwise lift
> Delsarte packet handoff
> generated tail45 strip
> lonely-transfer guardrail.
```

Pairwise observable: which obligation is closest to a sorry-free L7 closure
interface, with risk of falsehood as tie-break.  Switch/gauge: retain the
generated odometer or atom word before scalarizing to `p0`.

Fingerprint in the script is transitive: no directed 3-cycles, singleton SCCs,
and one Hamiltonian path.

## Assumption Challenge

The useful vertices are not runners, ratio channels, cells, or even the full
49 cell counts.  The quotient to a single row preserves the exact HYP-2736
integer defect after the HYP-2733 rotation law, but it destroys the visible
two-dimensional geometry.  That loss is acceptable only because the rotation
law is itself exact; without it, a row proof would not control the full cell
discrepancy.
