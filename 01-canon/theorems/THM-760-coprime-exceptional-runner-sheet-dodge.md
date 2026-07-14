---
id: THM-760
title: Coprime exceptional-runner sheet dodge — a scaled core plus one coprime runner preserves the core margin up to the sheet clearance (the 12-core case has M at least 1/13)
status: PROVED (elementary; uses only the settled LRC(<=13) floor in the LRC14 corollary)
source: codex-2026-07-14-S2
depends_on:
  - LRC(<=13)  # settled external input accepted by project policy
related: [THM-522, THM-531, THM-636, THM-724, THM-755, THM-757, HYP-6780, HYP-6785]
---

# THM-760 — Coprime exceptional-runner sheet dodge

## Statement

Let `P` be a finite set of positive integer speeds, let `c >= 2`, and let
`w > 0` satisfy `gcd(c,w)=1`.  Write `cP={cp:p in P}`.  If `t0` is any time
at which the core has margin

`min_{p in P} ||p t0|| >= delta`,

then one of the `c` lifted times

`t_k = (t0+k)/c`, `0 <= k < c`,

satisfies

`min_{v in cP union {w}} ||v t_k|| >= min(delta, 1/2-1/(2c))`.

Consequently,

`M(cP union {w}) >= min(M(P), 1/2-1/(2c))`.

In particular, if `3 <= |P| <= 12`, the settled LRC up to 13 gives
`M(P) >= 1/(|P|+1)`, while `1/2-1/(2c) >= 1/4 >= 1/(|P|+1)`; hence

> **`M(cP union {w}) >= 1/(|P|+1)`, and for a 12-speed core
> `M(cP union {w}) >= 1/13 > 1/14`.**

For arbitrary `|P| <= 12`, the unconditional consequence is
`M(cP union {w}) >= min(1/(|P|+1),1/4)`.

Equivalently: every primitive 13-speed family in which twelve speeds share a
common divisor `c>=2` is already strictly lonely at the stronger `1/13`
threshold.  No AP, interval-core, covering, diameter, or upper-speed
assumption is needed.

## Proof

For every `p in P` and every integer `k`,

`||(cp)t_k|| = ||p(t0+k)|| = ||p t0||`,

so all `c` lifted sheets preserve the entire core margin exactly.

For the exceptional runner, the phases are

`w t_k = w t0/c + wk/c (mod 1)`.

Because `gcd(c,w)=1`, multiplication by `w` permutes the residue classes
modulo `c`.  Thus these phases are a translate of the complete `c`-grid
`{0,1/c,...,(c-1)/c}`.  Some grid point lies within circular distance
`1/(2c)` of `1/2`; its distance from the nearest integer is therefore at
least `1/2-1/(2c)`.  Choose the corresponding `k`.  At `t_k` the core has
margin at least `delta` and the exceptional runner has margin at least
`1/2-1/(2c)`, proving the stated minimum.  Taking `t0` to maximize the core
proves the minimax inequality.  The LRC(<=13) corollary follows immediately.
∎

## What this closes and what it does not

- It strictly generalizes THM-724's AP-specific shallow-witness lemma
  `c*{1,...,12} union {w}`: the core can be any set of at most 12 speeds.
- It closes the whole codimension-one scaled-core lane exposed by HYP-6780,
  including the unbounded near-dilate rays of THM-757, without a raw-speed
  cutoff or finite enumeration.
- In the endogenous pair-sum blocker complex (HYP-6785), the AP-core special
  case lives on the sixfold ruler `q=13c`; the proof above is more general and
  works even when the core witness uses another pair-sum ruler.
- It does **not** close scale-quotient families with two or more exceptional
  residue classes, nor families in which no 12-speed subfamily has a common
  divisor.  Those are the genuine residual after this peel.

## Assumption challenge and Tournament Analysis

The proof's natural vertices are the `c` **witness sheets**, not runners or
arcs.  The sheet quotient preserves every core constraint and exposes only the
exceptional-runner phase; it destroys endpoint-owner and pair-sum labels, which
must be restored when more than one exceptional runner is present.  The
companion HYP-6785 computation therefore keeps the exact proof-obligation
hypergraph and uses a runner tournament only as a lossy diagnostic quotient.

*Verification:* `04-computation/lrc14_endogenous_blocker_complex_codex_S2.py`
and `05-knowledge/results/lrc14_endogenous_blocker_complex_codex_S2.out`.
