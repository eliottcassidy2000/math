---
id: LEM-010
title: Deterministic good-period existence for the finite-Vmax glue (THM-527-A) — for a k-element cluster E of co-offsets (k≤13) and ruler Vmax, a NONZERO good period j (maxgap{frac(e_i j/Vmax)}>1/7, i.e. the observer can sit ≥1/14 from all cluster runners in period j) is produced EXPLICITLY, with no equidistribution: (i) if spread(E) < 6·Vmax/7 then j=1 works (the wraparound gap 1−spread/Vmax > 1/7); (ii) if Vmax > 3^{k-1} then Dirichlet pigeonhole gives some j ≤ 3^{k-1} with all phases in a 2/3-arc, hence an empty arc ≥ 1/3 > 1/7. Since k≤13, only the bounded region {Vmax ≤ 3^12 AND spread ≥ 6Vmax/7} needs a finite check
status: PROVED (both (i) and (ii) are elementary; verified). (i) machine-verified 0 failures over 7669 clusters with spread<6Vmax/7; (ii) Dirichlet pigeonhole (standard), and the hard j=1-fails clusters empirically need only j*≤4 (≪ 3^{k-1}), never absent (0 no-good-j over the swept clusters). This replaces the soft "#arcs < ρ*·Vmax" discrepancy route (VACUOUS for structured 2-block/AP clusters, where #arcs≈1-2.6·spread) and the quantitative-equidistribution residual with an elementary closure plus a bounded finite check. Reduces THM-527-A to {Vmax ≤ 3^12=531441 AND spread ≥ 6Vmax/7}, extending kps-S30's exact M(S) check (V₀≤1001)
source: mac-mini-2026-07-08-S58
depends_on:
  - THM-527   # lonely-density reformulation; good period ⟺ maxgap{frac(e_i x)}>1/7 at x≈j/Vmax
  - THM-663   # the covering-case assembly (this is the advance on its remaining item (1))
related:
  - THM-530   # ρ*_{1/7} ≥ m_P (the density-floor floor; gives good periods have positive DENSITY, LEM-010 gives EXISTENCE)
  - THM-518   # Weyl decorrelation (the route this ELEMENTARY argument replaces for large spread)
---

# LEM-010 — Deterministic good-period existence (finite-Vmax glue)

## Setup

From the lonely-density reformulation (THM-527): a covering-saturated `S = P ∪ {Vmax − e_i}` is
lonely (`M(S) ≥ 1/14`) as soon as some **good period** `j ∈ {1,…,Vmax−1}` exists — one where the
`k` cluster phases `{frac(e_i · j/Vmax)}` leave a circular gap `> 1/7` (so the observer, runner
`Vmax`, can sit `> 1/14` from every cluster runner). Here `E = {e_i}` are the co-offsets,
`0 = e_0 < e_1 < … `, `spread = max e_i`, `k = |E| ≤ 13`, all `e_i < Vmax`.

The density floor (THM-530/THM-661) gives the good set positive **measure** `ρ* ≥ m_P`; LEM-010
gives an explicit good **period**, closing the finite-`Vmax` discretization THM-527-A — WITHOUT the
soft discrepancy bound `#{good j} ≥ ρ*·Vmax − #arcs` (which is vacuous when `#arcs > ρ*·Vmax`, as it
is for structured 2-block/AP clusters where `#arcs ≈ 1–2.6·spread`) and WITHOUT quantitative
equidistribution.

## (i) The wraparound lemma: `spread < 6·Vmax/7 ⟹ j = 1` is good

At `j = 1` the phases are `frac(e_i/Vmax) = e_i/Vmax` (as `0 ≤ e_i < Vmax`), all lying in the arc
`[0, spread/Vmax]`. The circular gap from the top phase `spread/Vmax` forward to `0` has length
`1 − spread/Vmax`. If `spread/Vmax < 6/7` this wraparound gap exceeds `1/7`, so
`maxgap{frac(e_i/Vmax)} ≥ 1 − spread/Vmax > 1/7` — `j = 1` is good. ∎

(This is the observer's own free arc: below `6/7` compression the cluster leaves a `>1/7` window
open at the far side of the circle. Verified: **0 failures over 7669** clusters with `spread<6Vmax/7`,
`k = 8,11,13`.)

## (ii) The Dirichlet lemma: `Vmax > 3^{k−1} ⟹` some `j ≤ 3^{k−1}` is good

Map each `j ∈ {0,1,…,3^{k−1}}` to `(⌊3·frac(e_i j/Vmax)⌋)_{i=1}^{k−1} ∈ {0,1,2}^{k−1}` (`3^{k−1}`
cells; `e_0 = 0` is omitted). With `3^{k−1}+1` values of `j`, two collide: `j_1 > j_2` with, for
every `i`, `frac(e_i j_1/Vmax)` and `frac(e_i j_2/Vmax)` in the same third ⟹ `‖e_i(j_1−j_2)/Vmax‖ <
1/3`. Set `j* = j_1 − j_2 ∈ {1,…,3^{k−1}}`. Then every phase `frac(e_i j*/Vmax)` lies within `1/3`
of `0` on the circle, i.e. in `[0,1/3) ∪ (2/3,1)` — an arc of length `2/3`. The complementary arc
`(1/3, 2/3)` of length `1/3 > 1/7` is empty, so `maxgap ≥ 1/3 > 1/7` — `j*` is good. ∎

`j* ≤ 3^{k−1} < Vmax`, so `j*` is a valid nonzero period. (The `2/3`-arc / `1/3`-box is the crude
choice that suffices; `k ≤ 13 ⟹ 3^{k−1} ≤ 3^{12} = 531441`. Verified: the hard clusters where `j=1`
fails need only `j* ≤ 4` in practice — the `3^{k−1}` guarantee is enormously loose, and a good period
was found for EVERY swept cluster.)

## Consequence — THM-527-A reduced to a bounded finite check

Combining (i) and (ii): a good period exists whenever `spread < 6Vmax/7` **or** `Vmax > 3^{k−1}`.
The complement — `spread ≥ 6Vmax/7` **and** `Vmax ≤ 3^{k−1}` — forces `spread ≤ Vmax ≤ 3^{12} =
531441`, a **bounded** configuration region. There the finite-`Vmax` glue is a finite check
(`ρ* ≥ m_P` and, empirically, `ρ_K ≥ 0.27 > 0` throughout, with kps-S30's exact `M(S) ≥ 1/14`
already covering `Vmax ≤ 1001`). So THM-527-A — the last analytic item of the covering case
(THM-663) — is reduced from "an integer-vs-real / equidistribution estimate at `Vmax ≤ 91^12`" to
"two elementary lemmas + a bounded finite check on `Vmax ≤ 531441`".

## (iii) The sharp bound: `j* = O(k)` — AP case PROVED (mac-mini-S59)

The smallest good period `j*` is empirically `≈ k`, not `3^{k−1}`: over `>90k` adversarial
spread-dense clusters, **`max j* = 2 / 11 / 13` at `k = 8 / 11 / 13`**, never absent. The worst
cases are **arithmetic progressions**, and for them `j* = O(k)` is PROVED:

> **AP lemma.** For `E = {0, d, 2d, …, (k−1)d}` (co-offsets), `j* ≤ ⌈7(k−1)/6⌉`. Proof: by
> Dirichlet, `∃ j ≤ ⌈7(k−1)/6⌉` with `‖jd/Vmax‖ < 6/(7(k−1))`; then the `k` phases
> `{i·frac(jd/Vmax)}` form a `k`-term AP of step `< 6/(7(k−1))`, spanning `< 6/7`, so they leave a
> gap `> 1/7`. ∎ (Verified: `max j* ≤ ⌈7(k−1)/6⌉ = 9/12/14` at `k=8/11/13`, tight-ish.)

If `j* = O(k)` holds for ALL clusters (AP case done; general strongly verified, max `j* = k`, no
exceptions), the bounded finite check collapses to `Vmax ≤ O(k)` — **trivially inside kps-S30's exact
`M(S)` sweep** — and THM-527-A is **fully closed**. The general `j* = O(k)` is the single remaining
analytic lemma of the whole covering case; its extremal case (APs) is proved and it has zero
counterexamples over 90k+ adversarial clusters. File:
`04-computation/lrc14_maxjstar_search_macmini_S59.{py,out}`.

## Why this is the right tool

The good period is the observer's own loneliness window, and (i)–(ii) find it by **clustering the
cluster**: compress all `k` phases into a `< 6/7` arc so a `> 1/7` window opens. (i) does it for free
when the cluster is not too spread; (ii) does it by an integer dilation `j*` that simultaneously
pulls every `e_i j*/Vmax` toward `0`. Neither needs the good set's measure, its arc-count, or its
equidistribution — the three quantities the soft route stalled on.

**Open (a strong, well-supported conjecture): `j* = O(k)`, not `3^{k−1}`.** The smallest good period
was measured over the `j=1`-fails clusters INCLUDING the adversarial APs `{0,d,…,(k−1)d}` with
`d ∈ [Vmax/14, Vmax/7]` (built precisely to defeat `j=1`): `max j* = 3 (k=11) / 7 (k=13)` over
`≈5.5k` such clusters, and `max j* = 3` over dense/perforated hard clusters — **never absent**. So
empirically `j* ≤ 7` for `k ≤ 13`, astronomically below the `3^{k−1}` guarantee. Proving `j* ≤ c·k`
(a successive-minima / three-distance bound on the dilation that clusters an AP: `j ≈ Vmax/d ≈ k`)
would shrink the finite check to `Vmax ≤ O(k)` and make THM-527-A **fully elementary**. Files:
`lrc14_deterministic_goodperiod_macmini_S58.{py,out}` (+ the adversarial-AP `j*` sweep).
