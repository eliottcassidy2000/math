# The far-element discrepancy is a Riesz–Dedekind ladder — but it is *dual* to the plateau, not the adversary's weapon

*kind-pasteur-2026-06-30. A reframing ("attacker = defender: the extremizers are Riesz products") put to a clean test on the verified `w·Δ_w` engine (`lrc14_uniform_C_growth_kps.py`), which confirmed the structural half and **corrected** the adversarial half. The correction is the interesting part.*

## The claim under test

The hard multiscale bases that make `sup_w |w·Δ_w|` large — e.g. `{0,1,2,30,31,32,60,61}` (published `≈3.91`) — are **coherent towers** `B ⊕ σ·B`, whose indicator Fourier transform factorizes as a Riesz product `∏_j \hat 1_B(σ^j ξ)`. The proof-side test measure (mac-mini, `the-key-is-a-riesz-product`) is *also* a Riesz product. So the reframing was: attacker and defender live in the same space of self-similar measures, `inf_S L(S)` is a minimax over Riesz products, and the adversary's optimal far element is the one that makes the base *more* self-similar.

## What the experiment confirmed

Reusing the verified engine (`G0` tent kernel, `cells_of`, `wD`, `supw`; benchmarks `0.73/2.54/3.91` reproduced exactly):

- **P1 — depth ladder (confirmed, and it is a *ladder*).** `sup|w·Δ_w|` grows **linearly** in the number of coherent blocks: `0 (degenerate) → 1.49 → 3.29 → 4.61` at scale 30 (slope ≈ 1.56), `0 → 2.23 → 5.13` at scale 50 (slope ≈ 2.5). Linear-in-depth is the signature of an **additive-over-scales** object — one Dedekind rung per scale (THM-563 / the Dedekind-ladder reflection). The per-block increment is *not* scale-free: it grows with scale separation (1.5 → 2.5), consistent with each rung approaching its frozen limit `C_sat` from below as the far block decorrelates.
- **P3 — exact Riesz factorization (confirmed, structural).** The clean tower `E = B ⊕ 30·B` satisfies `\hat 1_E(ξ) = \hat 1_B(ξ)·\hat 1_B(30ξ)` to `5·10⁻¹⁴`. This is true but near-definitional (a sumset always factorizes); it *names* the object, it doesn't do work.
- **P2 — "peak extends the tower" (partial).** Block-extending far elements (`w mod σ ∈ {0,1,2}`) are **strongly overrepresented** among the top peaks — the published extremizer's global peak `w=62 = 2·30+2` literally completes the third block — but they are not always the global argmax (the clean 3-block tower peaks at `w mod 30 = 8,9`). "The adversary extends the tower" is a tendency, not a law.

## What the experiment corrected — the load-bearing finding

A genuine **Sidon control** (Mian–Chowla `[0,1,3,7,12,20,30,44,65,80,96,122]`) broke the adversarial claim:

| card-9 base | span | plateau `Φ` | `sup\|w·Δ\|` | regime |
|---|---|---|---|---|
| consec `{0..8}` | 8 | **0.448** (cap `= 0.494`) | 0.93 | **binding** |
| tower `B⊕30B` | 62 | 0.261 | 3.29 | slack |
| Sidon | 65 | **0.093** | 3.65 | slack |

Two corrections fall out:

1. **The driver is spread / number of distinct scales, not coherence.** A wide Sidon set scores *as high as or higher than* the coherent tower (`9.74 > 4.61` at card 12). Both reach many scales — the tower by coherent repetition, the Sidon set by spreading. The tower is **not** the unique adversary. (My first "incoherent" control — random integers in a small box — was worthless: by pigeonhole such sets are dense in differences, i.e. *coherent*, so they scored high and looked like they refuted P1. Real incoherence must be Sidon.)

2. **The discrepancy is dual to the plateau; neither alone is a counterexample.** `sup|w·Δ|` is large *exactly where* `Φ` is small. The full coverage is `p_0 = Φ + Δ_w`, and it stays below `cap` in every regime: the binding case is the **narrow** consec (large `Φ`, tiny `Δ`), while wide/deep bases are automatically **slack** (tiny `Φ`, large but harmless `Δ`, since `Δ_w ~ (w·Δ_w)/w` is divided by a large `w`). This is the repo's `Φ`–`Δ` trade-off (HYP-2779, kps-S23) re-derived from a clean three-family cross-section.

## The reframing, sharpened

"Attacker = defender" survives, but not as "the extremizer is a Riesz product." The correct statement is a **conjugate-pair duality**:

> The plateau `Φ` (the certificate side, large for narrow/coherent-dense bases) and the far-element discrepancy `Δ` (the Riesz–Dedekind ladder, large for wide/deep bases) are a **conjugate pair on the near-cover**. They trade off: `Φ + scale·Δ`-type sums are what the cap bounds, and increasing one shrinks the other. The conjecture lives in their **sum**, and the binding locus is where they **balance** — the narrow, near-covering consec block — not at either extreme.

So the object to watch is not the discrepancy `Δ` by itself (unbounded, a red herring for a counterexample — MISTAKE-078's harmonic divergence is the same red herring) and not the plateau `Φ` by itself, but the **Legendre/Fourier-conjugate structure** linking them. That is why every "control the discrepancy absolutely" route dies and every route that *pairs* the two (the plateau recursion HYP-2675, the period-max certificate THM-563) makes progress: they respect the duality instead of fighting one half of it.

## Net

- **Not a new proof.** It tests a reframing, confirms its structural half (Riesz factorization; linear Dedekind depth-ladder), and **refutes its adversarial half** (spread, not coherence, drives the discrepancy; the tower is not uniquely worst).
- **The keeper:** the `Φ`–`Δ` conjugate duality is the thing to center. The discrepancy's unboundedness is not a threat because it is exactly compensated by plateau collapse; the binding case is the balanced narrow near-cover. Future "wide bound" work should target `Φ + Δ` jointly (2-D ET–Koksma / plateau recursion), never `Δ` alone.
- Artifacts: `04-computation/lrc14_riesz_depth_ladder_kps.py`, `05-knowledge/results/lrc14_riesz_depth_ladder_kps.out`.

## Addendum (same session): two follow-ups — a null and a clean two-scale picture

**The "spectral flatness beats additive energy" reframing is REFUTED.** Predictor shoot-out
over all 1716 sets `{0}∪(7-subset of {1..13})`, Kendall-τ vs `p_0`: additive energy `+0.393`
(27.5% discordant) *beats* spectral flatness `−0.341` (33%), diff-set size `−0.240`, width
`−0.184`. Additive energy remains the best *scalar* proxy — but none is clean (27–46%
inversions), confirming the repo's "consec-max is irreducibly aggregate" (the even-band
moments M2/M4, not any one scalar). Consec `{0..7}` is verified the `p_0`-argmax and the
AE-argmax. (`lrc14_spectral_predictor_and_dft_peaks_kps.py`, part A.)

**The `w·Δ_w` peak structure is a two-scale (arc × tower) resonance, and the DFT predicts the
envelope.** For the extremizer `E'={0,1,2,30,31,32,60,61}`: (B1) exact defected-Riesz identity
`\hat 1_{E'}(ξ) = \hat 1_B(ξ)\hat 1_B(30ξ) − e(62ξ)` (error 5.5e-14) — it is a 2-fold Riesz
product *minus one element*. (B2) `|\hat 1_{E'}|²` peaks exactly at `ξ = k/30` (the tower
scale). (B3) `|w·Δ_w|` peaks cluster at `w ≡ 0 (mod 30)` (the scale grid), while the DFT of the
sequence `w ↦ w·Δ_w` is dominated by **period 7** (the 1/7 danger arc). So `Δ` is a fast
period-7 carrier under a slow scale-30 envelope; the DFT predicts *where peaks can live* (the
envelope), not a single argmax. This is the honest form of "the extremizer's Fourier structure
organizes the discrepancy."

— Related: [[lrc14-thread]], `lrc-the-dedekind-ladder-of-far-coherence.md`, `lrc14-is-the-lonely-measure-and-the-key-is-a-riesz-product.md`, THM-563 (single-far period-max), HYP-2779 / kps-S23 (the trade-off), HYP-2675 (plateau recursion), OPEN-Q-108. Artifacts: `04-computation/lrc14_spectral_predictor_and_dft_peaks_kps.py`, `05-knowledge/results/lrc14_spectral_predictor_and_dft_peaks_kps.out`.
