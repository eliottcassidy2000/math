# THM-609 — Base Good-Region Floor (the one remaining lemma of the far-peel)

**Status:** VERIFIED (elementary proof below). This is **step 1** of opus's single-step far-peel decomposition of `CoveringFarLonely 22` — the sole remaining mathematics of the LRC(14) endgame per opus-S49. Handed to mac-mini; solved here as the continuity step of [[THM-608]] specialised to the LRC(≤13) margin.
**Source:** mac-mini-2026-07-03-S24
**Serves:** `LRCFarPeelCore.lean` / `LRCGoodPipeline.exists_lonely_of_goodRegion_pos` (opus/kps). Provides the hypothesis `0 < length (goodRegion B (1/14))` for any ≤12-speed base `B`.

Notation: `‖x‖ = dist(x, ℤ)`; the safe set at band `r` is `Safe_r(B) = { t ∈ ℝ : ‖b t‖ ≥ r  ∀ b ∈ B }`; `goodRegion B (1/14)` is (up to measure) `Safe_{1/14}(B) ∩ [0,1)`.

## Statement

Let `B` be a finite list of `m ≤ 12` nonzero integers, `V = max_{b∈B} |b|`, and assume the citation node
`LRCUpTo13`. Then

```
    length ( Safe_{1/14}(B) ∩ [0,1) )  ≥  1 / (91 · V)  >  0.
```

Concretely, `LRCUpTo13` yields a `t₀` with `‖b t₀‖ ≥ 1/13` for all `b ∈ B`, and the whole interval
`[t₀ − 1/(182 V), t₀ + 1/(182 V)]` lies in `Safe_{1/14}(B)`.

## Proof

`LRCUpTo13` applied to the `m ≤ 12` speeds gives `t₀ ∈ ℝ` with, for every `b ∈ B`,
```
    ‖b t₀‖  ≥  1/(m+1)  ≥  1/13.
```
The **margin** over the target band is `1/13 − 1/14 = 1/182`. Set `ρ = 1/(182 V)`. For any `t` with
`|t − t₀| ≤ ρ` and any `b ∈ B`, since `‖·‖` is `1`-Lipschitz,
```
    ‖b t‖  ≥  ‖b t₀‖ − |b|·|t − t₀|  ≥  1/13 − V·ρ  =  1/13 − 1/182  =  13/182  =  1/14.
```
Hence `[t₀ − ρ, t₀ + ρ] ⊆ Safe_{1/14}(B)`. This interval has length `2ρ = 1/(91 V)`. Intersecting with a unit
period (loneliness is `1`-periodic in `t`, so we may take `t₀ ∈ [0,1)` and the interval inside `[0,1)` after a
harmless unit shift) gives `length(Safe_{1/14}(B) ∩ [0,1)) ≥ 1/(91 V) > 0`. ∎

## Why it is exactly step 1

opus-S49's roadmap: `CoveringFarLonely 22` = one peel = **[step 1] base good-region floor** + [step 3, done]
`far_peel_length_pos` + [step 4, done] `exists_lonely_of_goodRegion_pos` + [data] steps 2, 5. Step 1 asks to
"upgrade the LRC(≤13) POINT to `length(goodRegion B (1/14)) > 0`, floor `≥ 1/(91·max B)`." That is exactly this
theorem, and the floor matches on the nose (`1/(91 V)`). The slow-runner-vs-wide-far tension that blocked the
*elementary* base closure (THM-608 scope note) does **not** arise here: the base is closed **wholesale by the
LRC(≤13) citation**, and only its continuity inflation is needed — which is unconditional.

## Lean formalization plan (for `TournamentH7`)

Target signature (in the `goodRegion` framework of `LRCGoodPipeline`):
```lean
theorem base_goodRegion_floor (cite : LRCUpTo13) (B : List ℤ)
    (hB : ∀ b ∈ B, b ≠ 0) (hlen : B.length ≤ 12) :
    0 < length (goodRegion B (1/14))
```
Steps:
1. `cite B.length (by omega) (fun i => B.get i) (…)` ⟹ `∃ t₀ : ℝ, Lonely (B.length+1) … t₀`, i.e.
   `‖b t₀‖ ≥ 1/(B.length+1) ≥ 1/13` for all `b`.
2. 1-Lipschitz margin (mathlib `abs_sub_round` / the corpus's `‖·‖` bound): `[t₀−ρ, t₀+ρ] ⊆ Safe`.
3. Bridge to `length (goodRegion …) > 0`: use the real-region calculus (`RRegion`,
   `exists_pos_interval_of_rlength_pos`, kps `LRCRealRegions`) — a rational sub-interval of `[t₀−ρ, t₀+ρ]`
   avoids every danger comb (each comb is closed and misses `t₀` by `≥ 1/182`), so it is `mem` of
   `goodRegion`, and a nondegenerate `mem`-interval forces `length > 0`. The margin machinery
   `FarElementRate` / `LRCFarElementRate` is the intended tool (opus-S49).

The only non-mechanical piece is step 3's ℝ→ℚ-region bridge, which lives in kps's/opus's Region API. The math
(steps 1–2) is complete and unconditional above.

## Consequence

With THM-609 (step 1) + opus `far_peel_length_pos` (step 3) + kps `exists_lonely_of_goodRegion_pos` (step 4) +
the DATA steps 2 (wrapped-comb length identity) and 5 (finite window `22 < w < ~78·max B`), the single-step
far-peel closes `CoveringFarLonely 22`, and `lrc14_of_covering_far_22` becomes unconditional modulo the
`LRCUpTo13` citation — i.e. LRC(14) is proved on the citation policy. THM-609 removes the last *genuine* lemma;
what remains (steps 2, 5) is bookkeeping and a finite computation.
