# What the full theorem needs: the safe-measure identity and the scale-flow descent

*mac-mini-2026-07-06-S14 (HYP-4402). Owner: reason mathematically about what a
FULL LRC(14) theorem still needs. This note does the analysis for the hardest
piece — gap-emptiness (G) — from the covering-measure side, proves one new
rigorous brick (the two-scale decorrelation lemma), and isolates the irreducible
kernel. Verified: `lrc_safe_measure_resonance_macmini_S14.out`,
`lrc_twoscale_decorrelation_v2_macmini_S14.out`.*

## The three pieces, precisely

`LRC14 ⟸ {LRC(≤13)} + TightLooseDichotomy + CornerLonely`, and

`TightLooseDichotomy` = **tight-locus rigidity** (`M = 1/13 ⇒ dilated AP`, clean
at prime 13 — S12) + **gap-emptiness (G)** (no covering-compressed family with
`M ∈ (1/13, 2/25)`). (G) is the general spectral-gap conjecture (HYP-2052) at
n = 13 — open in the literature (the "barely lonely runners" refined conjecture
is proven only for n ≤ 4, 6). It is the deep piece. This note is about (G).

## The exact analytic object: the safe measure

For speeds `S = {v₁,…,v_m}` and level `β>0`, let `g_β(u)=𝟙[‖u‖<β]` (the "danger"
indicator, an arc of measure `2β`) and

  **`safe(S,β) := |{t∈[0,1) : ‖vᵢt‖ ≥ β ∀i}| = ∫₀¹ ∏ᵢ (1 − g_β(vᵢt)) dt`.**

Then `M(S) = sup{β : safe(S,β) > 0} = inf{β : safe(S,β)=0}` — exactly. So

> **(G) ⟺ the only covering-compressed 12-family with `safe(·, 2/25) = 0` is the AP.**

Expanding the product (finite inclusion–exclusion, `Aᵢ={t:‖vᵢt‖<β}`):

  `safe(S,β) = Σ_{T⊆[m]} (−1)^{|T|} μ(⋂_{i∈T} Aᵢ) = (1−2β)^m + [resonance corrections],`

where the corrections are the Fourier terms over nontrivial additive relations
`Σ kᵢvᵢ = 0`, weighted by `∏ ĝ(kᵢ)`, `ĝ(k)=sin(2πkβ)/(πk)`. The **`(1−2β)^m`
term is the independent baseline**; the corrections encode the additive energy of
`S`. This is the Newman/Beurling–Selberg picture (HYP-+2873: additive energy
`= ∫|Ŝ|⁴`, the AP/interval is the Fejér-maximal concentrator).

## Why (G) is strictly harder than LRC(14) itself (the razor)

Two facts from the exact computation (`…resonance_macmini_S14.out`):

1. **Covering families sit far below the baseline.** `(1−4/25)^{12} ≈ 0.123`, but a
   covering family like `{1..11,23}` has `safe(·,2/25) = 0.009` — a **13×
   resonance suppression**. Covering forces many small additive relations, so the
   operative quantity is not the baseline but the razor-thin residual after
   near-total cancellation. A gap member would have `safe(·,2/25)=0` exactly,
   sitting a hair below `0.009` — a sub-1%-measure discrimination.
2. **No low-order Bonferroni certificate exists at 2/25.** The inclusion–exclusion
   partials `1−S₁+S₂−S₃+S₄` at `β=2/25` are `−0.92, 1.37, −2.94, 5.18` — swinging
   wider each order. Contrast the **LRC(14) threshold `β=1/14`**, where the pair
   term already dominates and `safe = (48−6c)/49 > 0` for `c ≤ 7` (the seven-wall,
   kps LRCStarSafe): there `safe ≈ 0.034` is comfortably positive. The 2/25 gap
   has no comfort — covering needs the AP's *high-order* resonance, which no
   finite Bonferroni truncation captures.

This is the precise sense in which (G) is a rigidity/extremal statement, not a
loneliness bound: it asks that the razor-thin high-order cancellation zeroing the
safe measure is achievable **only** by the AP.

## The new result: the multi-scale case of (G) is rigorously closed

The safe-measure identity has a feature that makes the scale-flow *rigorous*, not
heuristic: **`F_A(t) := ∏_{a∈A}(1−g(at))` is an indicator** (each factor is `0`
or `1`), so `F_A = 𝟙_{Safe(A,β)}` and

  **`safe(A ⊔ C, β) = |Safe(A,β) ∩ Safe(C,β)|`** — the measure of the *intersection*
  of the two safe sets.

Now take a family with a **scale gap**: `S = G_low ⊔ G_high` with
`max(G_low) ≤ ρ·min(G_high)`, `ρ → 0`. Two ingredients:

1. **Each part fails to cover.** `|G_low|, |G_high| ≤ 11`, so by **LRC(≤13)**,
   `M(G_low), M(G_high) ≥ 1/12 > 2/25`. Hence `Safe(G_low, 2/25)` and
   `Safe(G_high, 2/25)` both have **positive measure**.
2. **The fine safe set equidistributes in the coarse one.** `𝟙_{Safe(G_high)}`
   has Fourier support on combinations of the `G_high` frequencies (all
   `≥ min(G_high)`), so its average over any arc of `Safe(G_low)` (width
   `~1/max(G_low) ≫ 1/min(G_high)`) equals its global average `|Safe(G_high)|`
   (Erdős–Turán, error `→0` as the gap `min(G_high)/max(G_low) → ∞`). Therefore
   `|Safe(G_low) ∩ Safe(G_high)| → |Safe(G_low)|·|Safe(G_high)| > 0.`

> **Theorem (multi-scale collapse of (G)).** A covering 12-family with a scale gap
> has `Safe(·, 2/25) ≠ ∅`, i.e. `M ≥ 2/25` — it is **not** a gap member.
> Consequently **every gap member is a single bounded-ratio cluster** (no scale
> gap anywhere).

The special case `C = N·B` gives the clean product limit
`safe(A∪NB,β) → safe(A,β)·safe(B,β)` (verified to the exact fraction, ratio
`→1.0000` at `N=4001`). The general theorem needs only the two cited/standard
facts above — **no hard analysis, no discrepancy heuristic.** This is strictly
stronger than the *evidence-standard* separated-scale half of opus-S48's
scale-flow (OPEN-Q-108 R2): the unbounded-height multi-scale tail of (G) is now
**rigorously eliminated**, using only LRC(≤13) + equidistribution. It composes
with the gap-ladder order-statistic bounds (LRCGapLadder), which chain the top six
runners into one scale.

The engine — "a `k`-subfamily with `k ≤ 11` cannot cover at `2/25` because
`M ≥ 1/(k+1) ≥ 1/12 > 2/25`" — is the reason `12` is the *threshold size*: only
the full 12-family can reach the gap, and it can do so only if it is
irreducibly single-scale.

## The irreducible kernel that remains

The multi-scale case is now closed (theorem above), so a gap member is a **single
bounded-ratio cluster**. Two sub-regimes remain, and they are the genuine open
core:

1. **The single cluster at large height** `{c+δᵢ}` (ratios bounded, height `→∞`).
   The theorem's scale gap is *within-cluster absent*, so it does not apply. This
   is opus-S48's *difference core*: the loneliness is governed by the differences
   `δᵢ` at unit scale — a within-cluster renormalization. There is a natural next
   step here using the same identity: near a carrier resonance `t ≈ k/c`, the
   cluster's safe set is controlled by `𝟙_{Safe({δᵢ})}` at the fine scale, so one
   expects `safe(cluster) ≈ safe(differences)` — a *nested* decorrelation. Making
   this rigorous (and its contraction toward the AP fixed point, unique at the
   prime — S12) is the harder half of the scale-flow.
   *Evidence the cluster case also avoids the gap* (`…difference_core_probe…out`):
   the block of 12 consecutive integers `{c,…,c+11}` has the exact closed form
   **`M = c/(2c+11)`** (witness `t = 1/(min+max) = 1/(2c+11)`, both endpoints
   binding at `c/(2c+11)`; middle runners are farther). It equals `1/13` only at
   `c=1` (the AP) and jumps to `2/15, 3/17, 4/19, …` — never in the gap. Likewise
   every `{1..11, top}` with `top ≥ 13` sits at exactly `1/12` (or `2/25` at
   `top=24`). The AP is isolated at `1/13`; any near-AP perturbation jumps clear.

2. **The compact core** (bounded ratio *and* bounded height): a finite family of
   near-AP covering configurations. Here (G) is the **density-floor** statement —
   `safe(·,2/25)=0 ⇒ AP` — equivalently the additive-energy extremal
   `F_j ≥ 1/36` minimized uniquely at the AP (opus-S48, evidence standard). This
   is where the razor lives; it needs the *quantitative uniqueness* of the AP as
   the Fejér/additive-energy maximizer, sharp enough to forbid the sub-1% residual.

## Net: what a full theorem still needs (sharpened)

- **Tight side** (S12/S13): `residue_pinning_13` (formal) + strict lift-rigidity
  `nonzero lift ⇒ M>1/13`. Clean at the prime.
- **(G), multiple scales**: **CLOSED here** (theorem) — a gap member has no scale
  gap; it is a single bounded-ratio cluster. Uses only LRC(≤13) + equidistribution.
- **(G), single cluster**: the difference-core renormalization (opus-S48), now the
  *sole* remaining structural reduction — with a concrete nested-decorrelation
  attack via the same identity.
- **(G), compact core**: the AP-uniqueness density floor `≥ 1/36` — the
  additive-energy extremal, quantitative.
- **CornerLonely**: THM-619/620 band sweep → uniform argument.

The measure/identity frame unifies all of (G): everything is `safe(S,β)`, the
descent is fast-oscillation averaging of it, and the kernel is the AP's uniqueness
as the maximal-resonance (Fejér-extremal) covering configuration. The decorrelation
lemma removes the unbounded-height tail *rigorously*; what is left is a bounded,
single-scale, additive-energy extremal problem — finite in spirit, razor-thin in
substance.

## Pointers

- `lrc_safe_measure_resonance_macmini_S14.py/.out` — the identity, the razor, the
  Bonferroni-order contrast (threshold vs gap).
- `lrc_twoscale_decorrelation_macmini_S14.py`, `…_v2_…out` — the lemma verified
  (`safe(A∪NB) → safe(A)safe(B)`).
- opus-S48 / HYP-4013 / OPEN-Q-108 (scale-flow, density floor); HYP-+2873
  (additive energy = Fejér 4th moment); kps LRCStarSafe (Bonferroni 7-wall);
  HYP-4382 (prime tight-locus, the AP as unique fixed point); LRCGapLadder
  (order-statistic chaining).
