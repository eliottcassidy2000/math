---
id: THM-504
title: The |T|=3 conditional-convergence barrier of the LRC singular series — why the Paley–Zygmund within-level cancellation cannot improve the almost-Sidon route, and where the real cancellation lives
status: PROVED-in-scope — the half-period positivity lemma (A) is proved; the conditional-convergence structure (B) is proved heuristically + verified numerically; the ruling-out of within-level √-cancellation (C) is established computationally with mechanism. The full inf L>0 (C'(14)) remains open.
source: kind-pasteur-2026-06-14-S2
depends_on:
  - THM-501   # the singular series L(S) (kind-pasteur S6)
  - THM-503   # structure: 7-vanishing, almost-Sidon loose class (mac-mini S1)
related:
  - THM-446   # the additive-relation / Sidon ladder
  - HYP-2510  # the Paley–Zygmund signed-short-vector hope — REFUTED here (within-level)
---

# THM-504 — the |T|=3 conditional-convergence barrier

Builds on THM-501/503: `L(S) = (6/7)^13 + Σ_{exact relations Σ t_v v = 0}
(6/7)^{13−|T|}(−1)^{|T|} ∏_{v∈T} s(t_v)`, `s(t) = sin(πt/7)/(πt)`. THM-503 carved out
the **almost-Sidon** loose class (`|T|=2` relations only) via the absolutely-convergent
pair bound `|P(a,b)| ≤ g²/(3 v_a v_b)`. The open prize is `inf L>0` over the dilated
cores `d·{1,…,12}∪{r}`, which carry abundant `|T|=3` relations. The dispatch: improve
the **Paley–Zygmund signed-short-vector** route — exploit the `(−1)^{|T|}` signs for
cancellation where the absolute bound fails. This file determines exactly what
cancellation is and isn't available at `|T|=3`.

## A. Half-period positivity (PROVED)

> **For `|t| ≤ 6`, `s(t) > 0`.** Proof: `s` is even (`s(−t)=s(t)`), and for
> `1 ≤ t ≤ 6`, `πt/7 ∈ (0, π)`, so `sin(πt/7) > 0` and `s(t) = sin(πt/7)/(πt) > 0`. ∎

Consequence: **every `|T|=3` relation with all `|t_v| ≤ 6` has `∏ s(t_v) > 0`** — the
dominant (small-coefficient) relations contribute with a *single* sign. Sign variation
in `s` requires a coefficient `|t| ≥ 8` (`t mod 14 ∈ {8,…,13}`), where `|s(t)| ≤
1/(π|t|) ≤ 1/(8π)` is already small. So the negatives are both *rare* (need large
coefficients) and *small* (the `1/|t|` decay) — they cannot cancel the small-coefficient
positives. This is the structural reason the next two parts hold.

## B. The |T|=3 correction is CONDITIONALLY convergent (proved heuristically + verified)

Fix the support to `|T|=3`. The signed sum `Σ_3(S) = Σ_{exact 3-term relations} ∏ s(t_v)`
**converges** to a definite limit, while the absolute sum `A_3(S) = Σ |∏ s(t_v)|`
**diverges** (logarithmically).

*Divergence of `A_3` (heuristic, confirmed):* 3-term exact relations with
`max|t_v| ≤ T` number `∼ T²` (two free coefficients on each rank-1 relation lattice),
and a relation with coefficients of size `∼T` has `|∏ s| ∼ 1/|t_a t_b t_c| ∼ 1/T³`. So
each dyadic coefficient-scale contributes `∼ T²·T^{−3} = 1/T`, and `A_3 ∼ Σ 1/T = ∞`
(log-divergent). *Convergence of `Σ_3`:* the sinc oscillation `sign s(t) = sign
sin(πt/7)` (period 14) makes the signed sum a conditionally convergent oscillatory
lattice sum.

*Verified* (`04-computation/lrc14_pz_signed_shortvector_kps2.py`, coefficient cutoff
`Tmax = 13, 20, 27`): the **signed** sum is stable (`|Σ_3| ≈ 0.205` generic-coprime,
`≈ 0.95` for `7·{1,…,12}∪{r}`, moving `< 3%` across cutoffs) while the **absolute** sum
grows (`A_3 = 0.314 → 0.342 → 0.356` generic; `1.41 → 1.57 → 1.66` core) and the ratio
`|Σ_3|/A_3` falls `0.65 → 0.60 → 0.58` — the signature of conditional convergence.

> **Corollary.** THM-503's absolute/triangle-inequality method (the engine of the
> almost-Sidon class) **cannot reach `|T| ≥ 3`**: not because the bound is loose, but
> because `A_3 = ∞`. The almost-Sidon class (`|T|=2` only, where the pair sum is
> absolutely convergent) is *exactly* the absolute-convergence boundary of the route.

## C. The Paley–Zygmund within-level √-cancellation is RULED OUT (with mechanism)

The natural Paley–Zygmund improvement (HYP-2510) hoped the signed `|T|=3` sum enjoys
*square-root* cancellation `|Σ_3| ∼ √N₃ · (typical |∏s|) ≪ A_3` (signs quasi-random).
The data refute this: the cancellation is **constant-factor** (`ratio ≈ 0.58–0.65`),
not `∼ 1/√N₃ ≈ 0.03` (`N₃ ∼ 10³–10⁴`). The signs are *positively correlated*, by part A:
the dominant terms are uniformly positive, so the numerator `|Σ_3|` *converges to a
nonzero limit* rather than decaying like `√N₃`. The ratio drops only because the
*denominator* `A_3` diverges, not because the signed sum shrinks.

> **So within-level cancellation does not help.** The `|T|=3` correction is a definite,
> non-cancelling (in the PZ sense) negative quantity `−(6/7)^{10} Σ_3(S)`; for the dense
> cores `Σ_3 ≈ 0.95`, contributing `≈ −0.20` to `L` (verified consistent with
> `L(core) ≈ 0.03` once the `|T|=2` positive and `|T|≥4` terms are included).

## D. Where the real cancellation lives — the reframe (the genuine improvement direction)

The cancellation that keeps `L > 0` is **across support levels**, in the `(−1)^{|T|}`
alternation, not within a level:
`L = (6/7)^13 + (6/7)^{11} C_2 − (6/7)^{10} Σ_3 + (6/7)^9 Σ_4 − ⋯`,
a **conditionally convergent alternating-in-`|T|` series** with each level itself
conditionally convergent (B) and the level masses growing (the Vitali wall of S603 /
THM-501's caveat). The correct tool is therefore **Abel/summation-by-parts across the
support filtration combined with a Pólya–Vinogradov/Weil bound on each level's signed
sinc-lattice sum** — exactly the "archimedean singular-integral positivity" THM-503 (4)
identified, now localized: the object to bound is the *convergent* signed sum `Σ_k(S)`
per level, and the cross-level alternation. The Paley–Zygmund second-moment instinct is
right that absolute/first-moment bounds fail (B makes this sharp: `A_3=∞`); it is wrong
to expect the saving to be within-level √-cancellation (C). This redirects the route's
open core from "bound `Σ|∏s|`" (impossible) to "bound the conditionally-convergent
signed `Σ_k` and their alternating sum."

## Honesty / scope

- A is fully proved. B's divergence of `A_3` is a (standard) counting heuristic confirmed
  numerically; the convergence of `Σ_3` is numerical (stable to `<3%` across cutoffs),
  consistent with conditional convergence but not a closed-form proof. C is computational
  with a proved mechanism (A). None of this proves `inf L>0`; it determines the *shape* of
  any proof: cross-level Abel summation over conditionally-convergent signed sinc sums.
- Credits: the singular series is THM-501 (kp S6); the almost-Sidon class and the
  "Weil/Pólya–Vinogradov" flag are THM-503 (mac-mini S1). This file resolves HYP-2510
  (within-level PZ cancellation: REFUTED) and sharpens THM-503's open core.
