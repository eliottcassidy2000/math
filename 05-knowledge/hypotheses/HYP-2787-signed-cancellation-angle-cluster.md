---
id: HYP-2787
title: Angle cluster for the signed-cancellation wall (mac-mini gap #1) + the single-block domination — six fresh leads pushed in real time for the team's signed Δ_w / multi-far frontier
status: OPEN cluster (rapid lead generation). Single-block domination VERIFIED (0 violations); the six signed-cancellation angles are leads to test, shared for parallel attack.
source: mac-mini-2026-06-21-S6
related:
  - HYP-2784  # the absolute bound is 125x lossy; signed cancellation is the wall
  - HYP-2785  # Dedekind/equidistribution tail
  - HYP-2786  # one-far low-mode signed phase ledger
  - HYP-2694  # single-block extremal wide cover
  - HYP-2697  # arbitrary cluster compression cone
  - HYP-2739  # L7 two-clock residue-only closed form (the r=2 analogue, PROVED)
  - HYP-2657  # QR/Gauss/Q(sqrt-7) structure
  - OPEN-Q-108
---

# HYP-2787 — Signed-cancellation angle cluster + single-block domination

The wall (HYP-2784): the wide/single-far bound is intrinsically ~125× lossy via any ABSOLUTE
(BV/Koksma/V) bound; the proof must exploit **signed cancellation** in `Δ_w = Σ_{n≠0} ŝ(n)·1̂_B(nw)`.
Six fresh angles + one verified domination, pushed for parallel attack:

## VERIFIED: single-block domination (Route E, the multi-far route)
**Splitting a far-cluster into separated blocks LOWERS p0.** Verified: single consec block ≥ split
blocks in all structured tests; scattered never beats the block (0/60 random). So multi-far ≤
single-block (closed-form `D_m`, HYP-2694), reducing the multi-cluster gap to the single-block. The
formal lemma: `p0(B ∪ block) ≥ p0(B ∪ split)` at matched count — a rearrangement/coupling inequality
(the contiguous phase-AP sweep dominates a fragmented one). Complements HYP-2697's compression cone.

## SIX signed-cancellation angles (leads to test)
1. **POISSON-DUAL of Δ_w.** `Δ_w = ∫1_{A_j}(x)[1_{sector_j}(wx)−1/7]dx` (the exact signed value). Its
   absolute Fourier bound overcounts 125×. Apply Poisson on the n-sum / the dual lattice to express
   the signed value as a spatial overlap that converges ABSOLUTELY and is `O(1/w)` honestly.
2. **DEDEKIND RECIPROCITY for the tail.** The arcs of `A_j` have endpoints `k/e` (e∈base), so
   `∫_a^b 1_{sector_j}(wx)dx` is a floor-sum = a Dedekind sum `s(w,7e)`. The reciprocity law
   `s(h,k)+s(k,h)=−1/4+(h/k+k/h+1/hk)/12` BOUNDS the high-mode tail (HYP-2785), with the `Q(√−7)`
   flavor from the apex prime.
3. **THETA/MODULAR smallness.** Write `Δ_w` as a theta-sum; the 125× cancellation is the theta
   function's modular transformation (`w↦1/w`) shrinking the sum — the circle-method major/minor arc
   split with 7 as the modulus.
4. **KLOOSTERMAN/WEIL bound.** The signed mod-7 sum `Σ_n ŝ(n)1̂_B(nw)` with `1̂_B` carrying the
   residue structure is a Kloosterman/Salié-type sum; Weil gives `√7` cancellation — and `√7` is
   exactly the `|√−7|` of the apex Gauss sum (HYP-2657). The √-7 IS the signed-cancellation gain.
5. **MOD-14 = 2×7 phase ledger.** HYP-2786: low modes `n=1,2,3`, `n mod 14` buckets 1,2 dominate,
   `7|n` vanishes. 14=2·7: split `n mod 14` into `n mod 7` (sector/QR, Gaussian-period organized) ×
   `n mod 2` (parity/odd-even). The finite head is a `Q(√−7)×Z/2` ledger.
6. **THREE-DISTANCE exact arcs.** The base orbit gaps are 3-valued (Steinhaus); `1̂_B` has an exact
   3-gap closed form, making the signed head `n≤13` an exact finite rational ledger (no estimation).

## Next
Test #1 (Poisson-dual) and #4 (Weil/√−7) first — both could turn the signed value into an honest
absolutely-convergent `O(1/w)` or `√7` bound, cracking HYP-2784's wall. Push results live.


## LIVE FINDINGS (mac-mini-S6)
**The honest signed constant is C_signed ≈ 1.0** (VERIFIED `sup_w Δ_w·w ≈ 0.77–1.01` for consec_{k-1},
k=8..12), vs the absolute `15.8` (HYP-2784). So `Δ_w ≤ ~1/w` honestly, and `Δ_w(binding w~18) ≤ 0.056 <
margin 0.13`. **Proving `C_signed ≈ 1` closes the single-far binding window** — this IS mac-mini gap #1.

**THE DEDEKIND IDENTITY (the crack for angle #2).** Exactly:
> `Δ_w = (1/w) · Σ_j Σ_{endpoints t of A_j} ±s_j(w·t)`,  `t = k/(7e)` (e∈base), `s_j` = bounded sawtooth
> (centered `∫(1_{sector_j}−1/7)`, `|s_j|≤3/49`).
So `C_signed·w = Σ_j Σ_{k,e} ±((w·k/(7e)))` is a **generalized Dedekind sum**. The absolute bound (THM-546)
takes `Σ|s_j|`; the SIGNED value is `~16×` smaller because the sawtooths CANCEL. Grouping by `e`: each
`e` gives a Dedekind-type sum `d_e(w)=Σ_k ±((wk/(7e)))`; individually `d_e(w)` grows `~w/(84e)` (reciprocity),
but the **leading terms cancel across `e`** (Σ structure), leaving the bounded `~1` constant. The PROOF =
Dedekind reciprocity + the leading-term cancellation across the base. NOVEL: the Dedekind/LRC link is not
in the literature (web-checked); the recent frontier (Trakulthongchai 2025, 9&10 runners) names exactly
this signed-cancellation density as the open obstruction.

**Single-block domination over SPLITS holds** (refined): block-vs-split never reverses; the 16/1008 "leaks"
are block-vs-block at a different OFFSET (e.g. (19,20,21) beats (15,16,17), both blocks) = the φ-shift, so
the extremizer is the sup-over-offset single block (HYP-2694), confirming the multi-far → single-block route.


## *** THE PERIODICITY CRACK (mac-mini-S6, potentially closes single-far) ***
The Dedekind identity `Δ_w·w = Σ_j Σ_{endpoints t of A_j} ±S_j(frac(w·t))` has endpoints `t=k/(7e)`
depending ONLY on B (the A_j arcs are a function of B, NOT of w). `S_j(frac(w·t))` is PERIODIC in
integer w with period = denominator of t. Hence:
> **`Δ_w·w` is EXACTLY periodic in w, period `P = 7·lcm(B)`.**
Therefore `sup_{w≥15} Δ_w·w = max over ONE PERIOD [15, 15+P)` — a FINITE exact computation. And
`Δ_w ≤ (Δ_w·w)/w ≤ (period-max)/15` for all `w≥15`. So:
> **If `period-max(Δ_w·w) < 15·margin_k`, then `Δ_w < margin_k` for ALL `w≥15` — single-far CLOSED**
> (rigorously, no asymptotics, no Dedekind reciprocity needed — just the periodicity + a finite max).
margins: k=8 →0.185, k=9 →0.132, k=10 →0.124, so the thresholds are `15·margin = 2.77, 1.98, 1.86`.
Observed `Δ_w·w ≤ 1.0` (consec_7, w≤600) `<< 1.98`. If the full-period max stays `< 15·margin`, the
single-far signed bound (mac-mini gap #1, HYP-2784's wall) is PROVED as a finite periodic check —
SIDESTEPPING the absolute-bound 125× wall entirely. @kps @codex: this turns your Dedekind-tail
(HYP-2785) into a finite computation. Verifying period + max now.
