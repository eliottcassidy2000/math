---
id: THM-645
title: The offset tent law — pair masses of danger windows at arbitrary rational offsets equal θ² + (s·Λ(ψ) − r1r2)/(s²q1q2), where ψ = frac(α2q1 − α1q2) is the offset phase and Λ is an explicit arc-overlap tent with peak min(r1,r2) and valley (r1+r2−s)₊; THM-638's same-sign and mixed-sign laws are exactly the tent's peak and valley, and half-shifted window pairs have ψ ∈ {0, ½} by the parity of q — the even/odd duality as tent phase
status: PROVED (same two floor identities as THM-638, with the offset phase carried through; ~half page). Verified exhaustively: θ ∈ {1/7, 2/7} × five offset pairs × all coprime q ≤ 40 (0 violations; lrc14_offset_tent_law_klein_S162.out). Half-shift parity corollary verified.
source: klein-2026-07-07-S162 (HYP-4861)
depends_on:
  - THM-638   # the ψ = peak/valley special cases
related:
  - HYP-4841  # the 2-anchor tail / double cover (consumer)
  - HYP-4801  # boxeph's 2-anchor discharge (consumer: joint anchor-window masses)
  - THM-643   # mac-mini's parity structure (the sibling involution frame on the tournament side)
external: none (elementary).
---

# THM-645 — the offset tent law

## Statement

Fix `θ = c/s ∈ (0,1)` in lowest terms. For `q ≥ 1` and a rational offset `α ∈ [0,1)` let
`A_q^{(α)} := {x ∈ [0,1) : frac(qx) ∈ (α, α+θ]}` (meas = θ always; `α = 0` is THM-638's
same-sign window, `α = 1−θ` its mixed-sign window).

> **Theorem.** For coprime `q1, q2 ≥ 1`, offsets `α1, α2`, with `r_i = c·q_i mod s` and the
> **offset phase** `ψ := frac(α2·q1 − α1·q2)` (convention pinned by the exhaustive test — the mirrored convention with swapped r-roles is equivalent):
> `meas(A_{q1}^{(α1)} ∩ A_{q2}^{(α2)}) = θ² + (s·Λ(ψ) − r1r2)/(s² q1 q2)`,
> where `Λ(ψ) = s·|arc(ψ, ψ + r1/s] ∩ arc(0, r2/s]|` (arcs mod 1) — a piecewise-linear
> TENT in ψ with maximum `min(r1, r2)` (at aligned phase) and minimum `(r1 + r2 − s)₊`.

**Special cases.** `ψ = 0` gives `Λ = min(r1,r2)`: THM-638(i). `ψ = frac(−θ(q1−q2))`-type
anti-alignment gives the valley `(r1+r2−s)₊`: THM-638(ii). So the signed pair-mass law is
the peak and valley of one tent, and the general correction ranges over the whole tent:
`−min(r1r2, (s−r1)(s−r2)) ≤ s·Λ − r1r2 ≤ min(r1,r2)(s − max(r1,r2))`, always `O(1/(4q1q2))`.

**Half-shift corollary (the even/odd duality as phase).** For `α1 = 0, α2 = ½`:
`ψ = frac(−q1/2) = frac(q1/2) = 0` if `q1` is even (peak: positively correlated) and `½` if `q1` is odd —
where, whenever `r1/s, r2/s ≤ ½` and `r1 + r2 ≤ s`, `Λ(½) = 0` and the mass is
`θ² − r1r2/(s²q1q2) ≤ θ²`: **half-shifted windows on odd directions are negatively
correlated,** exactly. (This is the mechanism behind the S160 parity dichotomy of the
2-anchor tail, now quantified at the pair level.)

## Proof

As in THM-638. `A_{q1}^{(α1)} = ⊔_{j∈Z_{q1}} ((j+α1)/q1, (j+α1)/q1 + θ/q1]` and similarly
for `q2`. Left-endpoint differences are `(jq2 − lq1)/N + (α1q2 − α2q1)/N`, `N = q1q2`; by
Bézout `jq2 − lq1` covers `Z_N`, so the offsets form the shifted grid `(Z_N + ψN')/N`
(`ψ` as defined, lifted). Hence `m = Σ_{δ∈Z_N} f((δ+ψ')/N)` — the sampled cross-correlation
on a shifted grid — and exchanging sum and integral,
`m = ∫₀^{w1} c_ψ(t) dt`, `c_ψ(t) = #{δ : (δ+ψ')/N ∈ [t−w2, t)} = a1 + 1[frac(Nt − ψ') < r1/s]`
by the floor identity (with `Nw2 = q1c/s = a1 + r1/s`). Integrating the indicator over
`(0, w1]` with `Nw1 = a2 + r2/s`:
`m = a1w1 + (1/N)( a2·r1/s + λ )`, `λ = meas{y ∈ (0, r2/s] : frac(y − ψ') < r1/s}` —
exactly the arc overlap `Λ(ψ)/s`. Expanding `a_i = (q_i c − r_i)/s` as before gives the
stated form. ∎

## Applications

- **Joint 2-anchor window masses** (boxeph's `PA₂` program): the events "arc `(0,θ)` empty"
  and "arc `(½, ½+θ)` empty" are intersections of complements of `A_d^{(0)}` and
  `A_d^{(½·…)}`-type windows; every pairwise term in a Hunter/cherry/LP treatment of the
  UNION `W₀ ∪ W_½` is an instance of this law. The parity split (even vs odd `q`) is the
  exact skeleton of the S160 dichotomy.
- **Shifted-anchor Chung–Erdős / CE pair sums** (monad's program at anchors other than the
  aligned 14-grid): all cross-anchor pair terms are tent values.
- **Conditional ledgers**: G_P-intersected pair terms decompose over P's arcs into offset
  windows — the ledger's pair layer at any anchor.

## Verification

`04-computation/lrc14_offset_tent_law_klein_S162.py` (+ .out): exhaustive law-vs-engine
agreement over θ ∈ {1/7, 2/7} × offset pairs {(0,½), (0,⅓), (¼,½), (0,2/7), (⅓,¼)} ×
all coprime pairs q ≤ 40 (0 violations), plus the half-shift parity corollary cases.
