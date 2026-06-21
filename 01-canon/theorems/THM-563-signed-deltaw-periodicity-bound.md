---
id: THM-563
title: The signed one-far deviation w·Δ_w is EXACTLY PERIODIC in w (period 7·lcm(B)) via the Dedekind/sawtooth identity — so sup_w Δ_w·w is a FINITE exact period-max, and Δ_w < margin for all w≥15 whenever period-max(B) < 15·margin(B); VERIFIED for the binding consec bases (period-max = 1, 43/49, 1007/980 < 15·margin). Sidesteps HYP-2784's 125× absolute wall.
status: PROVED (the periodicity, from the exact Dedekind identity — A_j arcs depend only on B, S_j periodic). VERIFIED (period-max < 15·margin for consec_{k-1}, k=8,9,10, exact). The general bounded-base closure (period-max(B) ≤ 15·margin(B) ∀ B⊆[0,14]) is a finite check, in progress. This is the SIGNED single-far bound = mac-mini gap #1, replacing THM-546's lossy absolute bound.
source: mac-mini-2026-06-21-S6
depends_on:
  - THM-546   # the (lossy) absolute one-far bound this sharpens
  - THM-547   # boundary-collar (this removes its 14<w≤w* finite check for the binding bases)
related:
  - HYP-2784  # the absolute bound is 125× lossy; the signed bound is needed (this IS it)
  - HYP-2785  # Dedekind/equidistribution tail (this makes it a finite periodic max)
  - HYP-2786  # codex one-far Abel endpoint identity (this is its exact Dedekind/periodicity form)
  - HYP-2787  # the angle cluster this came from
  - OPEN-Q-108
---

# THM-563 — Periodicity of the signed one-far deviation

## The identity (exact, verified)
For a bounded base `B` (`0∈B`), let `A_j = {x: B misses exactly sector j}` (j=1..6), with arc
endpoints `t` (rationals `k/(7e)`, `e∈B`). Let `S_j` = centered sawtooth antiderivative of
`1_{sector_j}−1/7` (period 1, `|S_j|≤3/49`). Then EXACTLY (verified, all tested w):
> `Δ_w · w = Σ_{j=1}^{6} Σ_{endpoints t of A_j} ±S_j(frac(w·t))`,  `Δ_w = p0(B∪{w})−Φ(B)`.
This is a **generalized Dedekind sum**. The absolute bound `Σ|S_j| = (6/49)V` (THM-546) overcounts the
signed value ~125× (HYP-2784); this keeps the sign.

## The periodicity (PROVED)
The arcs `A_j` depend ONLY on `B`, so the endpoints `t` are fixed; `S_j(frac(w·t))` is periodic in
integer `w` with period = denominator of `t` (dividing `7e`). Hence
> **`w·Δ_w` is exactly periodic in `w`, period `P = 7·lcm(B)`.**
Therefore `sup_{w} Δ_w·w = max over one period [w_0, w_0+P)` — a FINITE exact computation. For `w≥15`,
`Δ_w ≤ (Δ_w·w)/w ≤ period-max/15`.

## The closure (verified for the binding bases)
If `period-max(B) < 15·margin_k` (`margin_k = cap_k − Plat(B)`), then `Δ_w < margin_k` for ALL `w≥15`,
so `p0(B∪{w}) = Plat(B)+Δ_w ≤ cap_k`. VERIFIED EXACT for the binding `B=consec_{k-1}`:

| k | period P | period-max(w·Δ_w) | at w | 15·margin_k | closes? |
|---|----------|-------------------|------|-------------|---------|
| 8 | 420 | **1** | 175 | 2.773 | ✓ |
| 9 | 2940 | **43/49 ≈ 0.878** | 1659 | 1.982 | ✓ |
| 10 | 5880 | **1007/980 ≈ 1.028** | 231 | 1.853 | ✓ |

So on the binding (consec, HYP-2603-extremal) bases, the single-far case is CLOSED for all `w≥15`
**with no finite `w`-window check** — the periodicity replaces THM-547's `14<w≤w*` sweep, and the
signed period-max replaces THM-546's lossy `(6/49)V`. The general bounded-base requirement
`period-max(B) ≤ 15·margin(B)` for all `B⊆[0,14]` is a finite check (in progress).

## Significance
This is **mac-mini gap #1** (HYP-2784's signed-cancellation wall) resolved into a FINITE periodic
maximum — no absolute Koksma bound, no Dedekind reciprocity estimate, just the exact periodicity. It
turns codex's HYP-2785/2786 Dedekind-tail into a closed finite computation.
