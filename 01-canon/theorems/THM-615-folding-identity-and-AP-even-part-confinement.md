---
id: THM-615
title: The folding identity M(2U ∪ {w₁,w₂}) = max_t min(g_E(t), Ψ(t)) for even 2U + two odd tighteners (Ψ = max(min a,b; ½−max a,b)), and — via a finite mod-24 check — M(2·{1..11} ∪ {w₁,w₂}) = 1/12 for ALL odd w₁,w₂, proving the m=2,f=2 confinement (gap = 1/84 > 0) for the EXTREMAL AP even parts c·{1..11}. Reduces inf_U gap(U) > 0 to the crisp conjecture M(2U ∪ 2odd) ≥ 1/12 (verified min = 1/12); soft measure bounds are vacuous (λ(U) ≪ 2/7), so the general case needs the argmax arithmetic.
status: PARTIAL — the FOLDING IDENTITY and the AP-even-part confinement (M(2·{1..11} ∪ 2odd)=1/12 ∀ odd, hence for all dilations c·{1..11}) are PROVED (elementary + finite mod-24 check; verified exactly). The general uniform gap inf_U gap(U) ≥ 1/84 is a CONJECTURE (M(2U ∪ 2odd) ≥ 1/12, verified min 1/12) — the confinement core, unproven.
source: opus-2026-07-03-S64
depends_on:
  - THM-612   # Lemma D (mac-mini): the m=2,f=2 confinement setup; this is a global (folding) complement to its R-covering view
  - LRCUpTo13 # M(U) ≥ 1/12 for 11-runner U (used in the reduction)
related:
  - HYP-4068  # opus-S62 REFUTED (MISTAKE-101): u_max NOT bounded; confinement is a uniform gap — this makes that precise
  - HYP-4062  # kps: tight-locus AP rigidity — the AP even part here is that AP, dilated
  - THM-613   # slope idea; MISTAKE-101 (the sampling artifact that reframed this)
results:
  - 04-computation/lrc14_folding_identity_AP_confinement_opus_S64.py
  - 05-knowledge/results/lrc14_folding_identity_AP_confinement_opus_S64.out
---

# THM-615 — the folding identity and AP-even-part confinement

**Setup.** `S = E ∪ {w₁,w₂}`, `E = 2U` (11 even runners, `U` any 11 positive integers), `w₁,w₂` odd.
`g_E(t) = min_{u∈U} ‖2ut‖`, `M(S) = max_t min_{v∈S} ‖vt‖`. Confinement (m=2,f=2) = "no such `S` is tight
at `1/14`", i.e. `gap(U) := min_{odd w₁,w₂}(M(S) − 1/14) > 0`. Since `M(S) ≤ M(E) = M(U)` and (LRC≤13)
`M(U) ≥ 1/12`, the clean target is `M(S) ≥ 1/12` (⟹ `gap ≥ 1/12 − 1/14 = 1/84 > 0`).

## Lemma 1 (folding identity, PROVED)
For odd `w`, `‖w(t+½)‖ = ½ − ‖wt‖`; and `g_E(t+½) = g_E(t)` (each `‖2u(t+½)‖ = ‖2ut + u‖ = ‖2ut‖`,
`u∈ℤ`), so `g_E` is `(+½)`-periodic. Pairing `t` with `t+½` and using
`max(min(g,X), min(g,Y)) = min(g, max(X,Y))`:
> **`M(S) = max_t min( g_E(t), Ψ(t) )`,  `Ψ(t) = max( min(a,b),\ ½ − max(a,b) )`,  `a=‖w₁t‖, b=‖w₂t‖`.**

`Ψ(t) ≥ 1/12 ⟺` NOT extremity (not "one `‖wᵢt‖ < 1/12` and the other `> 5/12`"). *(Verified exactly.)*
So `M(S) ≥ 1/12` **iff** some point with `g_E ≥ 1/12` is non-extremity for `(w₁,w₂)`.

## Lemma 2 (AP even part, PROVED via finite mod-24 check)
For `E = 2·{1..11}`, `g_E` attains its max `M = 1/12` at exactly **8 points, all of denominator 24**. Hence
`M(S) = 1/12` (its max possible value) **iff** some argmax point has both `‖wᵢ·(a/24)‖ ≥ 1/12` — a condition
on `(w₁ mod 24, w₂ mod 24)`. Checking all 78 odd residue-pairs mod 24: **every one gives `M(S) = 1/12`**
(distinct residues verified directly; equal residues give `Ψ ≥ 1/4` automatically). Therefore
> **`M(2·{1..11} ∪ {w₁,w₂}) = 1/12` for ALL odd `w₁,w₂`**,  so `gap = 1/84 > 0`.
By scale-invariance (`M(2c·{1..11} ∪ {w₁,w₂})`, the tighteners re-chosen) the dilations `c·{1..11}` — the
**min-gap extremizers** (MISTAKE-101) — are covered too. **Confinement is proved for the AP even parts.**

## The reduction, and why the general case is hard
`inf_U gap(U) > 0` follows from the conjecture **`M(2U ∪ {w₁,w₂}) ≥ 1/12`** for all 11-runner `U` and odd
`w₁,w₂`. Verified: an M-minimizing descent (commensurate seeds included) bottoms at **exactly `1/12`**,
extremized by the AP even parts. By Lemma 1 it is equivalent to: *for every `U`, some point with
`g_E ≥ 1/12` avoids extremity for every odd pair.* Soft/measure proofs of this are **vacuous**: the bound
`meas{Ψ < b} ≤ 4b` needs `λ(U) = meas(lonely U) > 2/7`, but 11-runner `U` have `λ(U) ≈ 0.05–0.09 ≪ 2/7`
(the danger sets genuinely nearly cover). So the general case needs the **argmax/commensurability
arithmetic** — precisely the confinement core (mac-mini's Lemma D, THM-612), here reduced to a sharp,
scale-invariant target with the extremal case dispatched.

## Lemma 3 (large tightener, PROVED — the loose end of the general case)
For `U` with `M(U) > 1/12`: if `max(w₁,w₂) > u_max / (6(M(U) − 1/12))`, then `M(S) ≥ 1/12`.
**Proof.** Let `I₀` be the component of `{g_E ≥ 1/12}` around the global argmax `t₀` (`g_E(t₀)=M(U)`).
`g_E` is `2u_max`-Lipschitz, so `|I₀| ≥ (M(U)−1/12)/u_max =: L`. WLOG `w₁ > 1/(6L)`. On `I₀` (length `≥ L`),
`‖w₁t‖` cannot stay entirely in `[0,1/12)` nor entirely in `(5/12,1/2]` (each such arc has length
`1/(6w₁) < L`); being continuous on the interval, it must take a value in `[1/12,5/12]` at some `t*∈I₀`.
There `Ψ(t*) ≥ 1/12` (if `‖w₂t*‖≥1/12` then `min ≥ 1/12`; else `max ≤ 5/12` so `½−max ≥ 1/12`), and
`g_E(t*) ≥ 1/12`, so `min(g_E,Ψ)(t*) ≥ 1/12`, whence `M(S) ≥ 1/12`. ∎ *(Verified: 360 large-tightener
families, 0 violations.)*

This is a **discrepancy statement**: a large tightener's orbit `{w₁t}` is dense enough on `I₀` to hit the
"moderate" band, so it cannot stay extremal. It disposes of the *loose* end; the residual is `M(U) > 1/12`
with **both** tighteners `≤ u_max/(6(M(U)−1/12))` — a low-frequency regime that shrinks to the AP (`M(U)→1/12`,
where the bound blows up) handled by Lemma 2. The hard core is exactly the small-tightener × near-AP corner.

## Status
- **Proved:** Lemma 1 (folding identity); Lemma 2 (`M(2·{1..11} ∪ 2odd) = 1/12` ∀ odd → confinement for AP
  even parts `c·{1..11}`, the extremizers); Lemma 3 (large tightener ⟹ `M ≥ 1/12`).
- **Conjecture (verified min = 1/12):** `M(2U ∪ 2odd) ≥ 1/12` ⟹ `inf_U gap(U) ≥ 1/84 > 0` = m=2,f=2
  confinement as a uniform gap.
- **Not claimed:** the general `U` (`M(U) > 1/12`) case — the argmax non-extremity, unproven. This is the
  confinement hard core; measure bounds cannot reach it.
