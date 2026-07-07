---
source: kind-pasteur-2026-07-06-S52
status: creative synthesis — the coarse/scale reduction, re-aimed at the SUPREMUM (Route 1)
  after opus-S130 retired Route 2; it survives, grounds the multi-scale families in the
  settled LRC(<=13), and sharpens the open core to SINGLE-SCALE decorrelation.
tags:
  - lonely-runner
  - LRC14
  - coarse-reduction
  - scale
  - decorrelation
  - route-1
  - synthesis
---

# The coarse reduction is a *sup* fact — it survives the Route-2 collapse and isolates the single-scale core

The owner asked to **prioritize creative reasoning and synthesis over formalization, and
work the scale/decorrelation argument or any better idea.** The same hour, opus-S130 +
mac-mini-S37 retired **Route 2** (the J-K → rank-2 → 1-D `(C)` thread my recent covering
work lived in): its top link bounds the spectrum's *accumulation points*, not the
*supremum* the LRC needs (MISTAKE-117), and its bottom `(C)` is not a finite covering
(MISTAKE-116). So closing `(C)` would **not** close LRC(14). The correctly-aimed route is
**Route 1: bound `M(v) ≥ 1/14` directly — the supremum, the right object.**

This is the best thing that could have happened to the scale idea, because the coarse
reduction was never really about `(C)`. It is a **general fact about the sup `M`**, so it
detaches cleanly from Route 2 and applies **directly to LRC(14) via the settled
LRC(≤13)** — with **no J-K, no accumulation points, no covering system.**

## The coarse bound (a sup fact, rigorous)

For nonzero speeds `v = (v₁,…,v₁₃)`, `M(v) = sup_t min_i ‖vᵢ t‖`. Suppose the speeds carry
a **two-scale structure** at some scale `L`: `vᵢ = aᵢ + L·kᵢ` with all `kᵢ ≥ 1` and
`|aᵢ| ≤ A`. Let `K = {distinct kᵢ}` (the **coarse family** — which cluster each speed sits
in). Take `K`'s optimal witness `t_K` (`min_j ‖kⱼ t_K‖ = M(K)`, `t_K ∈ (0,1]`) and set
`t = t_K / L`. Then

    ‖vᵢ t‖ = ‖aᵢ·t_K/L + kᵢ·t_K‖ ≥ ‖kᵢ t_K‖ − |aᵢ|/L ≥ M(K) − A/L,

so **`M(v) ≥ M(K) − A/L`.** One line, kernel-clean; the only hypothesis is `kᵢ ≥ 1` (no
ground speeds below the scale). Verified on 13-speed families: **0 violations.**

## The dichotomy this forces on LRC(14)

Let `r = #distinct kᵢ` = number of clusters at scale `L`.

- **`r ≤ 12` (two speeds share a cluster = a close pair at scale `L`).** Then `K` has
  `≤ 12` speeds, so **LRC(≤13) (SETTLED) gives `M(K) ≥ 1/13`**, and
  `M(v) ≥ 1/13 − A/L > 1/14` for `L > 182·A`. **LONELY — grounded entirely in the
  settled cases, no new analysis, no circularity.** (Verified: every `r ≤ 12` sample had
  `M(v) > 1/14`.)

- **`r = 13` (all clusters distinct).** `K` is again a 13-speed family, but of **strictly
  smaller height** (`height(K) = max kᵢ ≈ height(v)/L`). A counterexample `v` (`M < 1/14`)
  forces `M(K) ≤ M(v) + A/L < 1/14 + A/L` — a smaller **near**-counterexample. **DESCEND.**
  Each step divides the height by the scale gap, so it terminates at either a
  bounded-height family (a *finite* verifiable base) or a **single-scale** family.

- **Single-scale / compressed families (no scale gap).** Here `kᵢ = vᵢ`, `aᵢ = 0`, `L = 1`,
  so `K = v` and the bound is vacuous — the coarse reduction says nothing. **This is the
  residue.** The extremal single-scale family is the **AP `{1,…,13}`**, and `M = 1/14`
  *exactly* (verified) — the tight LRC(14) locus lives precisely here.

## What this buys, honestly

It does **not** close LRC(14). What it does:

1. **Re-grounds the compression / peeling reduction on the correct object.** The fleet's
   "peel the far element, reduce to compressed families" (THM-620/608) was framed inside
   Route 2 / `(C)`. The coarse bound shows it is really a **statement about the sup `M`** —
   so it *survives* opus-S130's Route-2 retirement intact, and it is now a clean
   *quantitative* mechanism (`M(v) ≥ M(K) − A/L`) rather than a structural hand-wave.

2. **Grounds the multi-scale branch in the settled cases.** Every family that clusters
   into `≤ 12` groups at any scale is lonely *by LRC(≤13)* — a genuine reduction of an
   infinite subclass to the settled frontier, with an explicit witness `t = t_K/L`.

3. **Sharpens the open core.** The genuine analytic difficulty — the **decorrelation /
   witness-density floor** that opus-S130 and mac-mini-S37 name as the load-bearing open
   piece — is now needed **only for single-scale (compressed) families**, not for all
   `Fin 13 → ℤ`. The multi-scale families are handled structurally. This shrinks the
   domain the hard Fourier estimate must cover to exactly the bounded-ratio families —
   which is also the domain of **Tao's finite reduction (2018)** (LRC reduces to speeds
   bounded by a function of `n`). The coarse reduction is a **self-contained instance** of
   Tao's "large speeds reduce," proved directly from LRC(≤13) with an explicit constant
   (`L > 182·A`).

## The synthesis: structured ⊕ generic

LRC(14) splits, on the sup, into two complementary halves:

| branch | families | tool | status |
|---|---|---|---|
| **structured** | multi-scale (a scale gap; `≤ 12` clusters or a descent) | coarse bound `M(v) ≥ M(K) − A/L` + **LRC(≤13)** | **grounded** (this session) |
| **generic** | single-scale / compressed (no scale gap) | **decorrelation** (Riesz/Fourier witness-density floor, Route 1) | **open** — the honest core |

The two halves are exhaustive and disjoint: a family either has a scale gap (structured,
reduces) or does not (single-scale, generic). The coarse reduction is exactly the tool
that clears the structured half and hands the generic half — sharpened to bounded-ratio
families — to Route 1's Fourier floor. That is where the "entirely new way of looking at
things" (Trakulthongchai, per Quanta 2026) is genuinely required; the multi-scale families
never needed it.

## Ledger

- **New (rigorous, verified):** the coarse bound `M(v) ≥ M(K) − A/L` (a sup fact, J-K
  free); `r ≤ 12 ⟹ M(v) > 1/14` via LRC(≤13); the height-descent for `r = 13`; the
  residue = single-scale families (AP `{1..13}` extremal, `M = 1/14` exact).
- **Reframed:** the compression/peeling reduction is a **Route-1 (sup)** fact, so it
  survives opus-S130; the analytic core (decorrelation) is **restricted to single-scale
  families**.
- **Does NOT claim:** a proof of LRC(14). The single-scale core is open (decorrelation),
  and the `r = 13` descent's finite base needs the small-family census (partly in hand).
- **Honestly absorbs:** opus-S130 (Route 2 retired, MISTAKE-116/117), mac-mini-S37
  (analytic core is decorrelation). My own Route-2 `(C)`/covering work (S41–S51) stands as
  *correct conditional math about the spectrum gap*, now known not to close LRC(14) via
  J-K — the coarse reduction is the piece of it that re-aims at the sup.

- **Files:** `lrc_coarse_reduction_direct_kps_S52.py`, output
  `lrc_coarse_reduction_direct_kps_S52.out`. HYP-4697.
- **Pointers:** opus-S130 (`the-route-2-audit-two-broken-links`); mac-mini-S37
  (`the-honest-state-…-decorrelation-is-the-core`); LRC14-PROOF-MAP (Route 1); Tao 2018
  finite reduction; Bedert arXiv:2511.16636 (Riesz products, the generic-branch technique).
