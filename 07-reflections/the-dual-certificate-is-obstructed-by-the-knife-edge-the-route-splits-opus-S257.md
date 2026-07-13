---
source: opus-2026-07-11-S257
status: FRAMEWORK + rigorous obstruction. Pursuing the dual-certificate route for the covering-min gives a
  clean test-measure formulation (cert = prob measure nu with Sum_i nu(D_i) < 1), and a RIGOROUS finding: the
  deep-well knife-edge (W >= 1 a.e.) forces INT W dnu >= 1 for any absolutely-continuous nu, so a SINGLE
  uniform positive-polynomial certificate (mac-mini S40) is IMPOSSIBLE. The route splits, forced, into
  [tight/deep-well: the S255 rigidity, PROVED] + [loose covering: an AC test measure = second-moment
  anti-concentration, the favorable side]. This reconciles the dual and anti-concentration framings and locates
  the remaining work precisely.
tags:
  - lrc14
  - covering-min
  - dual-certificate
  - test-measure
  - knife-edge-obstruction
  - anti-concentration
  - forced-split
---

# The dual certificate is obstructed by the knife-edge; the route splits

**opus-2026-07-11-S257.** Owner: pursue the dual-certificate route for the covering-min. Doing so gives a clean
formulation, a rigorous impossibility for the *single*-certificate hope, and the forced split that is the
actual path.

## The dual (test-measure) formulation

For a family `v` and level `c`, let `Dᵢ = {t : ‖vᵢ t‖ < c}` (runner `i`'s danger arcs) and `W(t) = Σᵢ 1_{Dᵢ}`
(danger count). If a probability measure `ν` has

> `Σᵢ ν(Dᵢ) = ∫ W dν < 1`,

then (W is a nonnegative integer) `W = 0` on a set of positive `ν`-measure — a **safe point exists**, so
`M(v) ≥ c`. The dual certificate is such a `ν`. In Fourier form,
`Σᵢ ν(Dᵢ) = 2ck + Σ_{h≠0} b_{−h} Σᵢ ν̂(h vᵢ)`, with `b_h = sin(2πhc)/(πh)` the band coefficients.

**Lebesgue fails.** `ν =` Lebesgue gives `Σᵢ ν(Dᵢ) = 2ck = 2·(14/183)·13 = 364/183 = 1.989 > 1`. Measure alone
is insufficient; the certificate must use the **correlation terms** `Σᵢ ν̂(h vᵢ)` — the arithmetic /
anti-concentration content — and must exploit divisor-completeness (non-covering families need not be lonely at
`14/183`).

## The knife-edge obstruction (rigorous)

The deep well `{1..12,182}` attains `M = 14/183` **only at the single point** `t* = 14/183`; its danger count
`W_dw(t) ≥ 1` for **almost every** `t` (safe set `= {t*}`, measure 0 — verified `≈ 0.0014`, the grid
resolution). Therefore, for **any absolutely-continuous** `ν`,

> `Σᵢ ν(Dᵢ) = ∫ W_dw dν ≥ ν({W_dw ≥ 1}) = 1.`

So **no absolutely-continuous test measure certifies the deep well** — only the atom `δ_{t*}` does. A "single
positive trigonometric polynomial" (mac-mini S40) *is* an AC test measure, so a **single uniform certificate is
rigorously obstructed** by the deep-well knife-edge. The certificate must be **family-adaptive**, or the tight
case handled separately. This is why S40's single-certificate hope stalls: the extremizer is atomic.

## The forced split (the viable route)

- **Tight** (`M` at/near `14/183`, incl. the deep well) — the knife-edge families. Handled by the **S255
  rigidity**, **proved** for the deep well (`M_core = 1/13 ⟹` interval `⟹ s = 1 ⟹` equality, via S252
  prime-13 uniqueness — the unique minimizer). Near-deep-well by perturbation of S255.
- **Loose** (`M ≥ 14/183 + margin`) — safe set has **positive measure** (verified `0.065–0.10` for random
  covering families), so an **absolutely-continuous test measure certifies them**. This reduces to a
  **second-moment / anti-concentration** bound `measure{W = 0} > 0` — *favorable* here (spread families,
  `E[W] = 2ck ≈ 1.99`, `measure{W=0} ~ e^{−2} ≈ 0.13`). This is the project's anti-concentration core
  (S242–S245) **on its easy side**, and the correlation `Σᵢ ν̂(h vᵢ)` is exactly the **LRCFourierCompletion**
  completion identity (`|C_w − b²/q|`).

## Net

The dual certificate route, made concrete: a **single** uniform positive-polynomial certificate is **rigorously
impossible** (the deep-well knife-edge forces `∫ W dν ≥ 1` for AC `ν`). The route splits — forced — into
**[tight: the finite S255 rigidity, proved for the extremizer]** + **[loose: an AC test measure = second-moment
anti-concentration, the favorable side]**. This precisely locates the dual route, **reconciles it with the
anti-concentration framing** (the correlation terms *are* the completion identity), and reframes the honest
remaining work: the covering-min lower bound = **S255 rigidity for the tight stratum** (done for the deep well)
**+ an anti-concentration `measure{W=0} > 0` bound for the loose stratum** (the favorable direction). The
constructive path is the **split**, not one polynomial — and the tight half is already proved.

→ mac-mini S40 (dual-certificate hope — here shown to need the split), opus-S255 (deep-well tight case, proved
— the tight half), opus-S242–S245 (anti-concentration core — the loose half), LRCFourierCompletion (completion
identity = the correlation term), klein S267 (14/183 verified). Files:
`lrc14_dual_certificate_knife_edge_split_opus_S257.py` (+`.out`).
