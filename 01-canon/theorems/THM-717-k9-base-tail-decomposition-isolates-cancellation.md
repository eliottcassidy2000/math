---
id: THM-717
title: The k=9 base tail-decomposition isolates the signed cancellation — via the exact Abel identity J = 6T₁+4T₂+2T₃−2T₅−4T₆ (Tⱼ = P(N≥j)), the base row J ≥ 432/91 SEPARATES into a cancellation-free covering bound (POS) 6T₁+4T₂+2T₃ ≥ 4717/882 and a diameter-controlled bunching bound (BUNCH) p₅+3p₆ ≤ 1/7, with 4717/882 − 2/7 = 4465/882 = J(consec) ≥ 432/91; every prior absolute/Bonferroni bound on J failed on the −2T₅−4T₆ cancellation, which this decomposition confines to the single small (BUNCH) term (equality at consec) — leaving the bulk (POS) as a monotone-covering functional the density-floor machinery can attack without the wall
status: IDENTITY PROVED (Abel summation + N ≤ 6 pointwise); SEPARATION VERIFIED EXACT (both piece-bounds universal over 92 377 primitive 9-cores in [1..19], zero violations, both extremal at consec); the two piece-bounds are SHARP CONJECTURES (equality at consec). STRENGTHENS THM-711/716: consec = the EXACT global J-argmin over 167 950 primitive 9-cores in [1..20] (J = 4465/882, decorrelation budget J_iid − floor = +3.71, J_iid = 6963918/823543). Does NOT prove the base outright — it REDUCES it to [(POS) coupled covering extremality, cancellation-free] + [(BUNCH) p₅+3p₆ ≤ 1/7 bunching], each cleaner than raw J.
source: klein-2026-07-11-S254 (HYP-6030)
depends_on:
  - THM-711   # mac-mini: the k=9 base J = E[N(7-N)] >= 432/91 (this decomposes its proof)
  - THM-716   # mac-mini: the (mu,Var) finite-dimensional reduction (complementary decomposition)
  - THM-661   # the density floor nu >= bar (the T_1 lower bound; POS generalizes it to higher tails)
related:
  - THM-710   # eigen-transfer: k>=10 inherit, so k=9 (+k=8) is the whole base
  - THM-714   # the k=8 cubic base (the deg-3 analog awaits its own tail decomposition)
external: Abel summation; Bonferroni inequalities; the "AP is the extremal coverer" principle (Freiman-adjacent).
---

# THM-717 — the k=9 base tail-decomposition isolates the cancellation

## Setup

`N(x)` = number of empty sectors among the seven arcs `[s/7,(s+1)/7)` for the phases
`{frac(e x): e ∈ E}`; `pⱼ = P(N = j)`, `Tⱼ = P(N ≥ j)`. The k=9 base row (THM-711) is
`J := E[N(7−N)] = 6m₁ − m₂ ≥ 432/91` for every 9-core.

## The exact Abel identity

`N(7−N)` takes values `[0, 6, 10, 12, 12, 10, 6]` at `N = 0,…,6` (and `N ≤ 6` always, since
at least one sector is occupied). Abel summation `E[v(N)] = Σⱼ Tⱼ (vⱼ − vⱼ₋₁)` with
`v = [0,6,10,12,12,10,6]` gives the **exact identity**

> **`J = 6T₁ + 4T₂ + 2T₃ − 2T₅ − 4T₆`**   (the `T₄` coefficient `v₄−v₃ = 0` drops out).

Equivalently `J = (6T₁+4T₂+2T₃) − 2(p₅+3p₆)`, since `2T₅+4T₆ = 2p₅+6p₆ = 2(p₅+3p₆)`.

## The separation (isolating the cancellation)

Write **POS** `:= 6T₁+4T₂+2T₃` and **BUNCH** `:= 2T₅+4T₆ = 2(p₅+3p₆)`. Then `J = POS − BUNCH`,
and the base follows from two bounds that are each **extremal at consec** (verified universal,
92 377 primitive 9-cores in [1..19], zero violations):

- **(POS) — cancellation-free covering bound:** `6T₁+4T₂+2T₃ ≥ 4717/882 ≈ 5.3481`.
  A nonnegative-weighted sum of the *monotone* good-set tails `Tⱼ = P(≥ j sectors empty) =
  P(≥ j runner-free 1/7-arcs)`. It is `E[w(N)]` for the nondecreasing weight
  `w = [0,6,10,12,12,12,12]` (= `N(7−N)` with its high-N dip filled flat). No subtraction —
  the absolute/Bonferroni methods that DIE on `J` (nine documented cancellation failures) have
  no cancellation to fight here.
- **(BUNCH) — diameter-controlled bunching bound:** `p₅ + 3p₆ ≤ 1/7`, equality at consec.
  `p₆ = P(all 9 phases in one 1/7-arc)` and `p₅ = P(in exactly two)` are near-origin/near-rational
  events whose measure shrinks with the diameter; consec (smallest spread) maximizes them,
  at `p₅+3p₆ = 1/21 + 6/63 = 1/7` exactly.

**Assembly:** `J ≥ POS − BUNCH ≥ 4717/882 − 2/7 = 4465/882 = J(consec) ≥ 432/91` (margin
`4465/882 − 432/91 = +0.3151`). The separation is TIGHT (both parts extremal at the same core,
consec), so it loses nothing — yet it confines the signed part to (BUNCH).

## Why this matters

The base extremality resisted every *absolute* bound because `J`'s `−2T₅−4T₆` term cancels
against `6T₁+…` (MISTAKES / the standing "no order-blind bound" law, 9 confirmations). This
decomposition **quarantines that cancellation** in (BUNCH) — a single term that is (a) small
(≤ 2/7), (b) sharp only at consec, and (c) diameter-controlled (bunching → 0 as spread → ∞).
Everything else, (POS), is a monotone covering functional in reach of the density-floor toolkit
(THM-661 already gives `T₁ ≥ bar`; the coupled `6T₁+4T₂+2T₃` floor is the remaining covering
content — note the *individual* tail minima are NOT aligned (`min T₂, min T₃` occur off consec:
`6·minT₁+4·minT₂+2·minT₃ = 4.876 < 5.033`), so (POS) is a genuine *coupled* covering extremality,
not a sum of separate floors). Complementary to THM-716's `J = μ(7−μ) − Var` frontier: that
reduces to a 1-parameter minimization; this isolates the cancellation.

## Status & handoff

- PROVED: the Abel identity; `J ≥ 6·P(N≥1) = 6ν` (weak form, insufficient alone since
  `6·bar₉ = 3.37 < floor`); the exact rational values at consec.
- CONJECTURED (verified universal, sharp at consec): (POS) `≥ 4717/882`, (BUNCH) `p₅+3p₆ ≤ 1/7`.
- NEXT: (BUNCH) is the crux and is self-contained — a sharp bunching bound (three-gap governs
  consec's phases `{jx}`, so its LHS is a three-gap computation). (POS) is the coupled covering
  floor — attack with the moment-LP / Bonferroni machinery, now safe from cancellation. The k=8
  deg-3 base awaits its own tail decomposition (its requirement function is cubic, not `N(7−N)`).

## Files

`04-computation/lrc14_J_landscape_twopole_klein_S254.py`,
`lrc14_J_decorrelation_full_klein_S254.py` (+ `.out`s).
