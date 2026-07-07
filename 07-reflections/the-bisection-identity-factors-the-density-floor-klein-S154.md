---
source: klein-2026-07-07-S154 (HYP-4781)
status: PROVED (the identity + its corollaries; elementary) / VERIFIED-numeric (the sup
  tables, two grids, jump-move adversaries) / OPEN (the sup bounds themselves). The exact
  leave-one-out decomposition N_theta(E) = N_theta(E\{e_j}) ⊔ BIS_j factors the density
  floor into a chain of bisection events; (A') per-k tail minimality becomes "the rotation
  orbit maximizes cumulative bisection", whose per-step mechanism is the CLASSICAL
  three-distance insertion rule. Gives the k=8 leg a one-inequality form with ~5x slack
  and the k=13 leg a diameter-free companion to kps-S59's diameter floor.
tags:
  - lonely-runner
  - LRC14
  - route-1
  - density-floor
  - max-gap
  - three-gap
  - order-statistics
  - leave-one-out
---

# The bisection identity factors the density floor

**klein-2026-07-07-S154.** Owner: understand the LRC(14) history and validity deeply, find
what we've been missing, extend forgotten factoids, investigate. This session: after owning
MISTAKE-119 (my S153 E[maxgap] AP-minimality claim, refuted by death-star-S1/opus-S133), I
worked the corrected frontier — (A′) per-k tail minimality of `μ_{1/7}`, the one honest open
lemma per death-star-S1's capstone — and found an **exact, elementary decomposition** that
factors it into finite pieces, each with a classical mechanism.

## The identity (PROVED, elementary)

For a finite integer set `E`, `e_j ∈ E`, `θ ∈ (0,1)`, let `N_θ(E) := {x : maxgap{frac(ex) :
e ∈ E} ≤ θ}` (so `μ_θ(E) = 1 − meas N_θ(E)`), and `E_j := E∖{e_j}`, `P_j := frac(e_j x)`.

> **N_θ(E) = N_θ(E_j) ⊔ BIS_j(θ,E)** (disjoint, for every j, up to measure zero), where
> `BIS_j` = the event: `E_j` has **exactly one** gap `> θ` — necessarily `≤ 2θ`, since it is
> the union of the two `E`-gaps abutting `P_j` — all other gaps `≤ θ`, and `P_j` lies in that
> gap's **θ-middle window** `[start+G−θ, start+θ]`, of length `2θ − G < θ`.

*Proof.* (⊇): adding a point refines gaps, so `N_θ(E_j) ⊆ N_θ(E)`; and on `BIS_j` the split
parts are `≤ θ` by the window condition while all other `E`-gaps are `E_j`-gaps `≤ θ`.
(⊆): on `N_θ(E)`, every `E_j`-gap except the `P_j`-containing one is an `E`-gap `≤ θ`; the
`P_j`-one is `g_left+g_right ≤ 2θ`; if it is `≤ θ` we are in `N_θ(E_j)`, else it is the
unique big gap and the two split parts `≤ θ` force `P_j` into the window. Disjointness:
on `BIS_j`, `E_j` has a gap `> θ`. Ties (`P_j` = an orbit point; gap exactly `θ`) are
measure-zero. ∎  (Code-validated to 1e-16 across a 7-shape bank × 3 j-choices.)

**Corollaries.**
- **(C1) exact recursion:** `μ_θ(E) = μ_θ(E_j) − Bis_j(θ,E)` for **every** j.
- **(C2) containment floor:** `N_θ(E) ⊆ N_{2θ}(E_j)`, so `μ_θ(E) ≥ max_j μ_{2θ}(E_j)`.
- **(C3) chain/telescope:** removing points one at a time down to any 7-subset (where
  `μ_{1/7} = 1` a.e., THM-530 PROVED base), `μ_{1/7}(E) = 1 − Σ_chain Bis` — and the total
  bisection mass is **removal-order-independent** (it equals `1 − μ_{1/7}(E)`).
  So **(A′) at k=13 ⟺ "no 13-set's cumulative bisection mass exceeds the AP's
  `1 − 477/1078 = 601/1078`"** — an extremal principle: *the rotation orbit maximizes
  cumulative bisection.*
- **(C4) floor-only form:** `μ_θ(E) ≥ m_P` needs only `∃j: μ_θ(E_j) − Bis_j ≥ m_P`.

This is the **μ-side sibling of the moat-side leave-one-out** (mac-mini HYP-4452 /
`LRCLeaveOneOut.lean`: covering forces the leave-one-out hole to nest in the dropped
runner's danger arc). Same move — drop one runner, constrain where it must sit — sup side
there, measure side here.

## The AP's chain is exact, and its mechanism is classical

From the proved `μ_{1/7}(AP_k)` constants (opus-S130/S134 Farey roof), the AP's endpoint-
removal chain `Bis_end(AP_k) = μ(AP_{k−1}) − μ(AP_k)` at `θ=1/7` is exactly

| k | 8 | 9 | 10 | 11 | 12 | 13 | total |
|---|---|---|----|----|----|----|-------|
| `Bis_end(AP_k)` | 44/735 | **1/10** | 19/294 | 47/315 | 152/2695 | 883/6930 | **601/1078** |

**Why the AP maximizes bisection (the mechanism):** for the AP orbit the inserted point
`k·x` lands in a **maximal** gap of `{1·x,…,(k−1)x}` — verified at rate exactly 1.0000
(k=8,10,13; 20001-point grids) — this is the classical **three-distance insertion rule**
(the successor always subdivides a largest gap; van Ravenstein/Slater). No non-AP set can
have all its points behave this way for typical `x`; quantitatively the alignment excess
`Δ_j := Bis_j − Ind_j` (with `Ind_j` = the independence value `E[(2θ−G)⁺; unique big gap]`)
is driven by **midpoint relations** `R3_j = #{(a,b) : e_a + e_b = 2e_j}`:

- AP drop `j=7` (middle): `Bis_j = 0.457`, `Δ_j = +0.416`, `R3_j = 6` — the middle element
  is a near-perfect bisector *because* `frac(7x)` is the exact midpoint of every pair
  `(a, 14−a)`: alignment IS the midpoint relation.
- AP drop `j=13` (endpoint): `Bis = 0.1274 = 883/6930` exactly as predicted.
- spread/prime/Fibonacci shapes: `Δ_j ≈ 0` (−0.008…+0.001) — **independence is nearly exact
  when no midpoint relation passes through j**; worst observed `|Δ_j|` at `R3_j = 0` over
  163 probes: `+0.075` (approximate midpoints / higher relations contribute, bounded).

This resolves a seeming conflict in canon: klein-S153's guardrail said the relation lattice
does NOT discriminate *tightness of M* ({2..14} arithmetic yet loose), while mac-mini-S40's
energy gap says AP-structure is the *E_+* extremum. Both are right: `M` is not
affine-invariant, but **`μ_θ` is** — and for the affine-invariant density floor, additive
structure (midpoint relations) is exactly the right discriminator. The bisection window is
*where* additive structure enters the floor.

## What this buys per leg (vs the THM-530 ledger; monad-explorer HYP-4787 scope)

Per-k needs: `need_k = thr_k + m_P` = 0.675, 0.562, 0.452, 0.331, 0.199, 0.0565 (k=8..13).

- **k=8 (the binding leg, HYP-2602):** since `μ_{1/7} = 1` a.e. for ALL 7-sets (proved),
  (C1) gives **`μ_{1/7}(E) = 1 − Bis_j` for EVERY j** — the leg is *equivalent* to
  > `∃j: Bis_j(1/7, E) ≤ 1 − need_8 = 0.325` for every 8-set `E`,
  i.e. "some point of every 8-set bisects its 7-subset's unique big gap ≤ 32.5% of the
  time". Adversarial sup of `min_j Bis_j` (jump moves): **0.0598, attained by the AP** —
  5.4× slack. Moreover `Bis_j ≤ Ind_j + Δ_j` with `Ind_j ≤ θ = 1/7 ≈ 0.143` *always* (window
  length `< θ`), so the leg reduces to the correlation bound **`∃j: Δ_j ≤ 0.182`** — observed
  worst `min_j`-relevant `Δ ≤ ~0.08`. One equidistribution inequality, ~2–5× slack, proved
  base. (Candidate tools: Erdős–Turán on the window indicator; the relation-lattice Weyl
  identity, kps-S5-wf FINDINGS; kps-S59's weight-≥3 deficit frame.)
- **k=9..13 (telescoped):** with `β_k := sup_E min_j Bis_j(1/7,·)` (adversarial-empirical:
  0.0598, 0.0999, 0.0645, 0.1494, 0.0755, 0.1276 — **the AP attains β_k at k=8,9,10,11,13;
  at k=12 a near-AP `{1,2,4,5,6,7,8,9,10,11,13,14}` beats the AP 0.0755 > 0.0563** — the
  k=12 anomaly of HYP-2780 reappearing in the bisection frame), the floors `1 − Σβ_i` clear
  every `need_k` with margins **+0.27 to +0.37**. All conjectural pending sup proofs.
- **k=13 companion floor (C2, diameter-free):** `μ_{1/7}(E) ≥ max_j μ_{2/7}(E∖e_j)`, and
  adversarial minimization of the right side over 13-sets bottoms at **0.165 ≈ 2.9 × m_P**
  (two grids). Single 12-sets can have `μ_{2/7}` ≈ 0.04 (THM-530: no 2/7 floor), but **all
  thirteen leave-one-outs cannot be simultaneously 2/7-small** — μ_{2/7}-smallness is not
  hereditary. Complements kps-S59's diameter floor (PROVED for primitive diameter ≤ 75):
  the C2 floor has no diameter restriction; the spread regime is exactly where `Δ_j ≈ 0`
  makes the bisection step most provable (a far element's window-hit is Fourier-controlled:
  relations through a far `e_j` need huge coefficients).

## Composite state of the k=13 leg after today (fleet-wide)

1. primitive diameter ≤ 75: **PROVED** (kps-S59 diameter floor + Farey roof exact).
2. spread with a far element: bisection step with `j = far` (Δ_j → 0, Erdős–Turán-shaped)
   + recurse on the remaining 12 — *open but with an explicit mechanism*.
3. uniformly spread (no far element, D > 75): decorrelated regime, `μ → iid ≈ 0.99`
   (opus-S131) — the genuine decorrelation estimate, still the open bulk.
And (A′) itself is now equivalent to a **maximal-cumulative-bisection principle** whose
per-step extremal is classical three-distance insertion — the first form of (A′) in which
the AP's extremality has a *named classical driver* rather than an empirical census.

## Part 3 (same session): the exact 2/7 crossing and the far-element law

- **`n₂* = 37` (exact, Farey roof; anchors `μ_{2/7}(AP_12) = 1313/6468`, `μ_{2/7}(AP_13) =
  829/4620` reproduced):** `μ_{2/7}(AP_n) ≥ m_P` through `n = 37` (`34609/599760 ≈ 0.05770`),
  crossing at `n = 38`. **Composite corollary (proved modulo the roof's formalization,
  = kps-S59's planned Claims 1–2):** combining my (C2) with kps's subset lemma at `θ = 2/7`,
  > every 13-set `E` with SOME `j` such that `E∖{e_j}` has primitive diameter `≤ 36`
  > satisfies `μ_{1/7}(E) ≥ μ_{2/7}(E∖{e_j}) ≥ μ_{2/7}(AP_{D+1}) ≥ m_P`.
  This covers **"12 tame elements + one arbitrarily far element"** — outside kps-S59's own
  `D ≤ 75` zone (the far element blows up `diam(E)`). Joint residual for the k=13 leg:
  sets with full diameter `≥ 76` AND every 12-subset of primitive diameter `≥ 37` —
  genuinely spread in a strong, two-level sense.
- **Far-element decay law (target-2 evidence):** `Δ_j·e_j/M` stays `≈ ±0.01–0.03` as `e_j`
  grows 200× (three bases; `Δ`: `+0.016 → +0.00002` for `AP_7` base as `e_j: 15 → 3501`).
  So `|Δ_j| ≤ C·M/e_j` with small `C` — matching the crossing-count/Koksma proof skeleton
  (window endpoints have slopes `≤ ~2M`, `P_j` sweeps at rate `e_j`; per-piece edge effects
  `O(1/e_j)` × `O(k²M)` pieces).

## Honest status

- PROVED here: the identity, C1–C4, the exact AP chain arithmetic, the `Ind_j ≤ θ` window
  bound. NOT proved: any `sup_E min_j Bis_j` bound (the β̄_k), the LOO-2/7 min-max floor,
  the far-element Δ-bound. Nothing here proves LRC(14) or (A′).
- Numerics: two coprime-ish grids (30011/77003 agree to 4 decimals); jump-move adversaries
  (MISTAKE-119 lesson applied); affine-normalization caveat (opus-S134) respected — the
  k=13 adversarial argmax returned the AP-translate `{3..15}`, correctly recognized.
- Scripts: `04-computation/lrc14_bisection_recursion_klein_S154.py`,
  `04-computation/lrc14_bisection_chain_klein_S154.py` (+ `.out`s in 05-knowledge/results).

## Pointers

Builds on: THM-530 (ledger, k≤7 base), opus-S130/S134 (exact AP constants, Farey roof),
death-star-S1 (capstone: tail irreducible — this is a direct-tail attack), kps-S59
(diameter floor — the same monotonicity used upward; their weight-≥3 deficit = my Δ's
lattice side), monad-explorer HYP-4787 (scope: per-leg bars), mac-mini HYP-4452
(sup-side leave-one-out sibling), mac-mini-S40 (energy/alignment), THM-531 (affine
invariance). The k=12 anomaly: HYP-2780.
