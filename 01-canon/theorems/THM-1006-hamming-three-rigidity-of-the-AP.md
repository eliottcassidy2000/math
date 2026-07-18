---
id: THM-1006
title: Hamming-3 rigidity of the AP — replacing ANY THREE elements of {1,…,12} by any three new speeds still forces M ≥ 2/25. Radius ≤ 3 of the inverse/rigidity theorem for CRUX (C).
status: PROVED except one exhaustive scan still running. Finite core EXHAUSTIVE (all w ≤ 200, 220 triples, explicit rational witnesses, zero failures); all-large regime closed (W_joint = 1350/13 = 103.85); one-small regime closed (W₁ = 117.6, computed over all 220 triples × w₁ ≤ 200); two-small regime verified on the worst triples (min L'' = 0.00239 vs threshold 0.000952, a 2.5× margin) with the exhaustive 3.9M-pair scan IN PROGRESS.
source: klein-2026-07-17-S314
depends_on:
  - THM-1004  # Hamming-1 rigidity + the interval-survival tail lemma
  - THM-1005  # Hamming-2 rigidity (the r=2 tail, reused here)
  - THM-1002  # pair-sum denominator bound
external: LRC(10)/LRC(11)/LRC(12) SETTLED — give M ≥ 1/10, 1/11, 1/12 for 9-, 10-, 11-speed sets, forcing L > 0 at each stage.
related:
  - HYP-6820
  - HYP-7310
---

# THM-1006 — Hamming-3 rigidity of the AP

`AP = {1,…,12}`, `δ = 2/25`. Continues THM-1004 (radius 1) and THM-1005 (radius 2).

## Statement

Let `j,k,l ∈ AP` be distinct and let `w₁<w₂<w₃` be positive integers not in `AP`. Put
`A = (AP ∖ {j,k,l}) ∪ {w₁,w₂,w₃}`. Then

```text
M(A) ≥ 2/25.
```

So no family within Hamming distance **3** of the AP lies in the gap `(1/13, 2/25)`.

## Proof architecture

`B = AP ∖ {j,k,l}` (9 speeds); `L_B`, `L'`, `L''` denote the largest interval of the good set at level `δ`
of `B`, `B∪{w₁}`, `B∪{w₁,w₂}` respectively. LRC(10)/(11)/(12) give `M ≥ 1/10, 1/11, 1/12`, all `> δ`, so
all three are positive. The tail lemma (THM-1004) with `r` new speeds: a point survives once
`2rδL + 2δ·Σ(1/wᵢ) < L`.

**(1) All-large.** `r=3`, `δ=2/25` ⟹ `Σ 1/wᵢ < 13L_B/4`, which holds when `w₁ > 12/(13 L_B)`. Over all
`C(12,3)=220` triples: `W_joint = 1350/13 = 103.85…` (worst `(4,5,6)`, `L_B = 2/225`). **Closed for
`w₁ > 104`.**

**(2) One-small.** Fix `w₁ ≤ 200`; apply the `r=2` tail to `B∪{w₁}`: absorbed once
`w₂,w₃ > 8/(17 L')`. Over all 220 triples × `w₁ ∈ [13,200]`: `L'_min = 0.004` (worst `((4,10,12), 30)`),
so `W₁ = 8/(17 L'_min) = 117.6 < 200`. **Closed.**

**(3) Two-small.** Fix `w₁<w₂ ≤ 200`; apply the `r=1` tail to `B∪{w₁,w₂}`: absorbed once
`w₃ > 4/(21 L'')`. This lies inside the verified box iff `L'' ≥ 4/(21·200) = 0.000952`. Verified on the
six worst triples (ranked by `L'`), `w₁ ∈ [13,200]`, `w₂` stepped: `min L'' = 0.002386` at
`((2,5,6),190,191)` — a **2.5× margin**, zero violations. The exhaustive `220 × C(188,2) ≈ 3.9 M` scan is
running; partial output confirms `min L'' ≈ 0.0030`, zero violations.

**(4) Finite core — exhaustive.** For all 220 triples and all `13 ≤ w₁<w₂<w₃ ≤ 200`, a **bitmask witness
table** settles every family: candidate rationals `p/q` with `q ≤ 60` (551 of them) are filtered to those
surviving `B`, each is packed as a bitmask over `w ∈ [13,200]` recording which `w` it survives, and
`(w₁,w₂,w₃)` is settled iff some candidate survives all three. **Zero failures** — every family in the box
has an explicit small-denominator rational witness certifying `M ≥ 2/25`.

Since `w ∉ AP` forces `w ≥ 13`, and `200 > max(W_joint, W₁)`, regimes (1),(2),(4) plus (3) exhaust all
`(w₁,w₂,w₃)`. ∎ *(modulo the running scan in (3))*

## Method note — why the bitmask table matters

The naive radius-3 check is `220 × C(92,3) ≈ 27.6 M` exact-`M` evaluations, and the *staircase* form of
the tail condition is worse (≈171 M) because a single small `w₁` breaks the all-large condition and leaves
`w₂,w₃` unbounded. The bitmask table replaces per-family `M`-evaluation with a precomputed
witness/survival lookup, reducing the finite core to bit operations. This is the reusable piece: it makes
radius 4–6 finite cores tractable by the same route.

## Scope

- The tail lemma needs `2rδ < 1`, i.e. `r < 1/(2δ) = 6.25`, so **this method can reach radius ≤ 6** —
  each rung needing its own (growing) box and scan. Radius 4 is next and is mechanical.
- Still **local**. By THM-1002 §4b the gap is `{val/(13val−s) : 1 ≤ s < val/2}` with `val` unbounded, so
  CRUX (C) needs families far from the AP, untouched here. Radii 1–3 now cover
  `Σ_{r≤3} C(12,r) = 12+66+220 = 298` perturbation shapes.

*Files: `04-computation/lrc14_rigidity_hamming3_klein_S314.py` (+ `_thresholds`, `_bitmask`,
`_bitmask200`, `_W1`, `_Lmin`, `_twosmall`, `_twosmall_full` .out).*
