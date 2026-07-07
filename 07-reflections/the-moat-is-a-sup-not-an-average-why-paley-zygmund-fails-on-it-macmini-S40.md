---
source: mac-mini-2026-07-06-S40
status: structural / negative-with-mechanism — Paley-Zygmund and additive-energy (unsigned second-moment) methods reframe but CANNOT close the moat; the moat is a signed-cancellation SUP phenomenon, not an average. Explains the fleet's PZ division (density floor closes; moat resists).
tags:
  - lonely-runner
  - LRC14
  - moat
  - paley-zygmund
  - second-moment
  - additive-energy
  - signed-cancellation
  - three-gap
---

# The moat is a sup, not an average: why Paley–Zygmund fails on it

Owner: *work the moat, think Paley–Zygmund.* I did — and the honest result is a
**negative one with a precise mechanism**, which sharply positions the fleet's
(successful) Paley–Zygmund work on the density floor.

## Setup

The moat: single-scale non-AP 13-family ⟹ `M(V) ≥ 1/13`. To show `M ≥ β` it suffices
that `safe(V,β) = |G_β| > 0`, `G_β = {t : ‖vᵢt‖ ≥ β ∀i}`. The second-moment
(Paley–Zygmund) route: `Z(t) = ∏ᵢ g(vᵢt)` with `g ≥ 0` supported on `{‖x‖ ≥ β}`, so
`{Z>0} ⊆ G_β` and

    |G_β| ≥ E[Z]² / E[Z²]     (Paley–Zygmund).

Both moments expand over the **relation lattice** `L(V) = {k ∈ ℤ¹³ : Σ kᵢvᵢ = 0}`:
`E[Z] = Σ_{k∈L(V)} ∏ĝ(kᵢ)`, and `E[Z²]` likewise with `g²`. This is my S24 frame.

## The finding: the moat is a signed-cancellation effect

With the triangular bump `g` on `[β, 1−β]`, `ĝ(m) = (−1)ᵐ h·sinc²(mh)`, `h = ½−β`,
so `E[Z] = h¹³·(1 + Corr(V))`, `Corr(V) = Σ_{k≠0}(−1)^{Σkᵢ} ∏sinc²(kᵢh)`. Two facts,
computed exactly:

1. **The unsigned bound is useless.** The total *absolute* correction
   `AbsCorr(V) = Σ_{k≠0} ∏sinc²(kᵢh) ≫ 1` for **every** family (min ≈ 6.9 over a
   broad scan; AP ≈ 41.4), while the main term is `1`. So a rigorous PZ-via-absolute-value
   bound (`AbsCorr < 1 ⟹ E[Z] > 0`) **never** fires. The smallness of `E[Z]` — and the
   AP's exact `E[Z] = 0` at `β = 1/13` — comes **entirely from signed cancellation** of the
   oscillatory lattice sum, not from small magnitude. The AP's tightness is *cancellation*,
   not smallness.

2. **`E[Z] > 0` is circular.** For an open bump, `E[Z] > 0 ⟺ |G_β| > 0 ⟺ M > β`. So the
   PZ numerator *presupposes* the moat. Paley–Zygmund cannot bootstrap it.

**The additive-energy byproduct (exact, clean, but unsigned).** The AP `{1,…,13}`
**uniquely maximizes additive energy** `E₊ = Σ_s r(s)²` among single-scale 13-families:
`E₊(AP) = 1469`, next value `≤ 1425` — an exact **gap of 44** (Schur-triple count ties at
42, so weight-4 energy, not weight-3, is the discriminator). This is the Freiman /
minimal-sumset signature (AP = minimal sumset = max energy). But `E₊` is the *second
moment* `E[Z²]`, which is **unsigned**: it only weakly correlates with `M` (`r ≈ −0.25`),
and at the moat boundary `β = 1/13` the AP is **sub-threshold** (`M = 1/14 < β`), so its
moments *vanish* — the energy↔margin link breaks exactly where the moat lives. Additive
energy is a genuine but **insufficient** proxy.

## Why this matters: the average-vs-sup division

`M(V) = maxₜ minᵢ ‖vᵢt‖` is an **L^∞ sup** (worst case over `t`). Moment methods bound
**L¹/L² averages** (`∫`). Passing from an average to the sup loses precisely the signed
cancellation — which is where the moat's entire content lives. This **explains the fleet's
Paley–Zygmund division**:

- **Density floor (an AVERAGE): PZ works.** `μ_{1/7} = P(maxgap > 1/7)` reduces by
  reverse-Markov to `inf_E E[maxgap] > 1/7` (kps-S57, a comfortable `+0.06` margin;
  opus-S131 `E[U]`; opus-S133). These are `L¹` order-statistic means — moment methods
  have full purchase, and kps-S57 notes "the density floor is NOT tight where the raw
  loneliness is," the same point from the other side.
- **Moat (a SUP): PZ resists.** The worst-case margin is the razor-thin `1/14 → 1/13`
  boundary, carried by signed cancellation. No second-moment estimate reaches it.

So the fleet's PZ successes are all — and can only be — on averaged quantities. The moat
is the residual **sup**, and my finding is the reason it must be attacked with a
**signed / location** tool, not a moment.

## The right tool (constructive complement)

The signed/sup nature points back to the **three-gap / continued-fraction frame**
(mac-mini-S15: *"the density floor IS the quantitative three-gap rigidity; detuning the
AP raises `g` and jumps `M` to the next rung"*) and the **denominator sieve** (oracle-S18).
The moat is: *the AP is the unique single-scale family pinned to `1/14` with no better `t`;
every perturbation opens a better-than-`1/13` `t`.* A **resonance-destruction** mechanism
is visible — perturbing the AP breaks its exact `q = 14` resonance (a moved element hits
residue `0` or collides mod 14) — but the transversal-mod-14 residue proxy is
**necessary, not sufficient**: `{1,…,12,27}` has transversal residues yet `M = 1/13`
(a strictly better `t` exists — the *values*, not just residues, decide). So the location
content is genuinely finer than a residue count; it is the three-gap rigidity itself.

## Where it stands

- **Ruled out for the moat:** Paley–Zygmund / second-moment / additive-energy (any unsigned
  estimate). Mechanism: `AbsCorr ≫ 1`; the moat is signed cancellation of the sup. Don't
  spend more fleet effort here.
- **Proper scope of PZ (confirmed):** the density floor (kps-S57/S58, opus-S131/133) —
  averaged quantities with comfortable margins.
- **Exact byproduct:** AP uniquely maximizes additive energy among single-scale (gap 44) —
  a Freiman certificate that a single-scale family *is* the AP, possibly useful to the sieve.
- **The moat's real tool:** three-gap location (S15) / denominator sieve (oracle-S18) —
  signed, not second-moment.

→ HYP-4767; cites mac-mini-S15 (three-gap frame), oracle-S18 (sieve); complements
kps-S57/S58 + opus-S131/S133 (density-floor PZ, the average side) and opus-S132
(consolidation); continues mac-mini-S39/S40 (the moat, additive energy).
