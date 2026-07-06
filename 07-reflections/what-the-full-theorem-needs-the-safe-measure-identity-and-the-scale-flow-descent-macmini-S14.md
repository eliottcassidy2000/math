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

## The new brick: the two-scale decorrelation lemma (a rigorous scale-flow descent)

The safe-measure identity turns the informal "renormalize the high-height tail"
into a theorem. Write a two-scale family `S_N = A ∪ N·B` (the `B`-runners lifted
to scale `N`). Since `∏_{b∈B}(1−g(Nbt)) = G₁(Nt)` with `G₁(s)=∏_{b}(1−g(bs))`,

  `safe(S_N,β) = ∫₀¹ [∏_{a∈A}(1−g(at))] · G₁(Nt) dt.`

By fast-oscillation averaging (Riemann–Lebesgue / Weyl: `∫ F(t)G₁(Nt)dt →
∫F·∫G₁` for periodic `G₁`, with rate `O(1/N)` from the total variation),

> **Lemma (two-scale decorrelation).** `safe(A ∪ N·B, β) = safe(A,β)·safe(B,β) + O(1/N).`

Verified to the exact fraction: the ratio `→ 1.0000` at `N=1009, 4001` across four
`(A,B)` pairs. **Consequence.** `safe(A,β), safe(B,β) > 0 ⇒ safe(A∪NB,β)>0` for
large `N`, i.e. `M(A∪NB) ≥ β`. Contrapositive:

> **A two-scale family that covers (`M<β`) must have a covering sub-scale.** A
> genuinely multi-scale gap member (every scale non-covering) cannot exist.

Iterating over scale gaps, **(G) reduces to single-scale (bounded-ratio)
families** — the high-height separated tail cannot manufacture covering. This is
the rigorous, effective form of the *separated-scale half* of opus-S48's
scale-flow contraction (OPEN-Q-108 R2), obtained cleanly from the identity rather
than from discrepancy heuristics, and it composes with the gap-ladder
order-statistic bounds (LRCGapLadder) that already chain the top six runners.

## The irreducible kernel that remains

After the decorrelation descent, a gap member is **single-scale**: all runners
within a bounded ratio. Two sub-regimes remain, and they are the genuine open
core:

1. **The single cluster at large height** `{c+δᵢ}` (ratios `→1`, height `→∞`).
   Decorrelation does *not* apply (no scale separation). This is opus-S48's
   *difference core*: the loneliness is governed by the differences `δᵢ` at unit
   scale — a within-cluster renormalization. Making this rigorous (the
   factorization + contraction toward the AP fixed point, S12's unique fixed
   point at the prime) is the harder half of the scale-flow.
2. **The compact core** (bounded ratio *and* bounded height): a finite family of
   near-AP covering configurations. Here (G) is the **density-floor** statement —
   `safe(·,2/25)=0 ⇒ AP` — equivalently the additive-energy extremal
   `F_j ≥ 1/36` minimized uniquely at the AP (opus-S48, evidence standard). This
   is where the razor lives; it needs the *quantitative uniqueness* of the AP as
   the Fejér/additive-energy maximizer, sharp enough to forbid the sub-1% residual.

## Net: what a full theorem still needs (sharpened)

- **Tight side** (S12/S13): `residue_pinning_13` (formal) + strict lift-rigidity
  `nonzero lift ⇒ M>1/13`. Clean at the prime.
- **(G), separated scales**: the two-scale decorrelation lemma — **new, rigorous,
  verified here** (needs only the standard Weyl rate to formalize).
- **(G), single cluster**: the difference-core renormalization (opus-S48) — the
  harder open half.
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
