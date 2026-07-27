---
id: HYP-9040
title: "Cancellation-completeness synthesis: parity, maximal rank, and Euler-characteristic selection"
status: >
  PARTLY REFUTED / SUPERSEDED by HYP-9045. The mod-two 91-line lane is
  structurally blind and the rank-twelve box is not an independent CRT-finite
  channel. The Euler lane remains OPEN after its first fast stratification was
  refuted at a degenerate stratum. This file is not a proved dependency.
source: mac-mini-2026-07-27-S142
related:
  - LEM-020-redei-involution-parity-layers
  - THM-1340
  - THM-2059
  - HYP-8871
---

# HYP-9040 -- cancellation-completeness synthesis

> **CURRENT-TRUTH REPAIR (2026-07-27):** HYP-9045 gives exact negative
> verdicts for Sections 1--2 and refutes the first fast fibre stratification
> proposed in Section 3. The sections below preserve the original questions,
> not current positive evidence.

This synthesis compares three possible ways to turn a difficult
noncancellation problem into a rigid finite obstruction.  The comparison is
heuristic: the lanes have different native carriers, and no reduction among
them is claimed.

## 1. Parity-force the 91-line

Apply the repository's Redei/palindrome mechanism--an involution with an odd
fixed layer forces nonvanishing--to the rank-eleven residual's CRT clock
discharge.  The exact candidate carrier is THM-2059's reduction-histogram dot
product modulo two.  The target would be a parity certificate for the mixed
`7 x 13` character rather than an analytic no-cancellation estimate.

**Cheapest decisive test:** enumerate the bounded rank-eleven packet with all
labels retained, construct the proposed involution, and determine whether its
fixed stratum has odd cardinality in every live histogram class.  A single
even fixed stratum is a hostile witness; an unlabelled count is insufficient.

## 2. Test the rank-twelve box first

The motivating pattern is that known collision obstructions occur at maximal
ambient rank (the dimension-three Jacobian counterexample and the rank-four
GMC core).  The conjectural LRC analogue is that a counterexample, if any,
should occur first in the rank-twelve maximal-minor box rather than the
rank-eleven star.  Exhausting that finite geometric channel would leave the
spectral rank-eleven channel isolated, but would not itself prove its rigidity.

**Cheapest decisive test:** discharge the exact finite rank-twelve box of
HYP-8871 with the inherited row, clock, and phase labels intact, and compare
its surviving obstructions with the rank-eleven cancellation kernel.  The
maximal-rank analogy alone has no logical force.

## 3. Euler-characteristic selection

THM-1340's dimension-three Keller engine contains a cuspidal rational curve
of Euler characteristic one.  The proposed selection law is that genuine
escape engines must contain a `chi=1` stratum, while the planar Keller ledger
cannot balance such a stratum in degree at least three.

For the dimension-three model, the first concrete audit is the stratified
Euler ledger obtained by viewing `K` as quadratic in `c`, hence as a double
cover of the `(a,b)`-plane:

```text
chi(K) = 2 - chi(branch curve),
1 = 3(1-chi(K)) + 2(chi(K)-chi(Sigma)) + |F_Sigma| chi(Sigma).
```

**Cheapest decisive test:** compute the branch curve, `Sigma`, and the fibre
cardinalities exactly, including points at infinity and singular corrections.
Until that audit is complete, the displayed balance is a proposed ledger, not
an established identity or a planar Jacobian obstruction.

## Guardrail

Parity, maximal rank, and Euler characteristic are three representations of
different objects.  A future bridge must name its map, preserved predicate,
lost labels, restoration sidecar, and hostile control before transferring a
conclusion between lanes.
