---
source: opus-2026-07-09-S181
status: CORRECTION of the opus-S179 "additive energy is THE single parameter governing looseness"
  overclaim, found via a Freiman-neighborhood deep-dive. Additive energy E(S)=#{a+b=c+d} is NECESSARY
  (the AP maximizes it, S180/Freiman) but NOT SUFFICIENT for tightness. TWO clean counterexamples:
  (1) 2-D GENERALIZED APs (GAPs) have HIGH E (~1000-1225, in S179's "near-AP" range E>1000) yet are
  LOOSE (L~0.09-0.15) with SHORT longest-AP (<=7, "dissociated"); (2) the TRANSLATED AP {1,3,..,25}=
  2*{1..13}-1 has MAXIMAL E=1469 (= the tight AP) yet is LOOSE (L=0.116). Matched-energy pair: GAP 7x2
  (dim2, E=1225, L=0.091) vs near-AP (dim1, E=1245, L=0.033) -- same E, 3x looseness. So E and the
  good-period dichotomy DISAGREE exactly on GAPs; the true governing object is the RESONANCE-LATTICE
  GEOMETRY (dimension + coherence), not any scalar. The tight locus is NARROW: exactly the dilates
  c*{1..13} (THM-610), because L is dilation-invariant but NOT translation-invariant (THM-522).
  PROOF UNAFFECTED: LEM-012/013 correctly key on longest-AP (structural), not E; the GAPs are inside
  LEM-013's longest-AP<=7 census. Reconciles opus-S113, HYP-2889, HYP-2637, THM-610, opus-S112.
tags:
  - lrc14
  - additive-energy
  - freiman-dimension
  - generalized-ap
  - resonance-lattice
  - correction
  - necessary-not-sufficient
---

# Additive energy is necessary, not sufficient — the tight locus is resonance geometry

**opus-2026-07-09-S181.** Owner: deep-dive Freiman-set concepts, then investigate past work for relations.
Doing so refuted my own opus-S179 claim that additive energy `E(S)=#{a+b=c+d}` is the SINGLE parameter
governing looseness.  The correction is clean, it has two explicit counterexamples, and it reconciles a
half-dozen prior threads that had each half-seen it.

## The overclaim and its refutation

opus-S179 asserted: "`E(S)` monotonically governs `L`; the additive-energy GAP (dissociated `E<500`,
near-AP `E>1000`) IS the good-period dichotomy."  Two families break this:

| set | Freiman dim | `E(S)` | longest-AP | `L` (lonely) | verdict |
|---|---|---|---|---|---|
| AP `{1..13}` | 1 | 1469 | 13 | 0 | tight |
| near-AP `{1..12}+20` | 1 | 1245 | 12 | 0.033 | near-tight |
| **translated AP `{1,3,…,25}`** | 1 | **1469** | 13 | **0.116** | **LOOSE** |
| **2-D GAP 7×2 P=9** | 2 | **1225** | 7 | **0.091** | **LOOSE** |
| 2-D GAP 4×4 P=5 | 2 | 1205 | 4 | 0.100 | LOOSE |
| 3-D GAP | 3 | 805–913 | ≤5 | 0.09–0.13 | LOOSE |

- **The translated AP `{1,3,…,25}` = `2·{1..13}−1`** carries the SAME maximal energy `E=1469` as the tight
  AP, yet `L=0.116` (loose).  `L` is dilation-invariant (`L(cS)=L(S)`, THM-522) but NOT
  translation-invariant — so a translate of a dilate of the AP is loose.  Max energy ≠ tight.
- **The 2-D GAPs** carry `E ~ 1000–1225` (squarely in S179's "near-AP" band `E>1000`) but have
  `longest-AP ≤ 7` (S179's "dissociated" band) and are LOOSE.  The matched-energy pair is decisive:
  GAP 7×2 (`E=1225`) vs near-AP `{1..12}+20` (`E=1245`) have the same energy but `L = 0.091` vs `0.033`.

So additive energy and the good-period dichotomy DISAGREE exactly on GAPs.  `E` is not the dichotomy
parameter, and it is not a sufficient statistic for `L`.

## What is actually true: resonance-lattice geometry, not a scalar

The lonely measure is a theta-sum over the RELATION (resonance) lattice `Λ(S)={t∈ℤ¹³ : Σtᵢvᵢ=0}`:
`L(S)=Σ_{t∈Λ}∏h(tᵢ)` (THM-515; opus-S112 shows two sets with the SAME `Λ` have IDENTICAL `L`).  So
looseness is a functional of the WHOLE lattice geometry, of which additive energy is only the leading
(order-2, `t·v=0` with `‖t‖₁=4`) shadow.  Three geometric features, none a scalar of `E`:

1. **Dimension** (Freiman rank of `Λ`).  A 1-D lattice (AP) concentrates all resonance on ONE modulus, so
   the danger events `‖vᵢτ‖≤1/14` align at a single `τ` → coherent over-covering → tight.  A `d`-D GAP
   spreads resonance over `d` independent directions with no common `τ` → loose.  This is HYP-2637's
   "dimension penalty ≈ ½ per dimension," now confirmed on the loneliness side (dim-2/3 GAPs are
   uniformly loose regardless of `E`).
2. **Coherence** (dilate vs translate).  Even at dim 1 and max energy, only the DILATES `c·{1..13}` are
   tight (THM-610 / kps-S37: the tight locus is *exactly* the dilated APs).  Translation shears the
   resonances off the common `τ=1/14` → loose (`{1,3,…,25}`).
3. **Energy** (order-2 density).  Necessary — the AP maximizes it (S180/Freiman `|S+S|≥2n−1`) — but the
   `{1,3,…,25}` and GAP counterexamples show it does not pin dimension or coherence.

The honest law: **tight ⟺ `Λ` is 1-dimensional AND coherent (a dilate of `{1..13}`).**  Additive energy is
the necessary order-2 shadow; the dimension penalty and the dilate-vs-translate knife-edge are the parts
`E` cannot see.  This is exactly opus-S113's "Freiman is necessary, not sufficient — structure × width,"
now seen from the LRC side (S113 was the tiling side: a tiler with *sub*-AP energy; here the dual, a
max-energy translate and high-energy GAPs that are *loose*), and HYP-2889's "energy is AP-facing, not
monotone."

## The proof is unaffected — only my S179 reframing was wrong

Crucially, the good-period dichotomy is keyed on **longest-AP** (the genuine 1-D structure), NOT on
additive energy — LEM-012 needs an AP of length `L≥k−6`; LEM-013 covers `longest-AP≤7`.  The GAPs have
`longest-AP≤7`, so they sit in the DISSOCIATED branch, and LEM-013's exhaustive census (`spread≤22`, 621k
clusters) already INCLUDES them (GAP 4×4 P=5 has spread 16).  They are loose (`L≥0.09`), exactly what the
dissociated branch handles.  So there is no hole: the fleet's dichotomy uses the right (structural) axis.
Only my opus-S179 gloss — "the additive-energy gap IS the dichotomy" — was the overclaim; the S180
anchor (AP = Freiman minimal-sumset = the unique tight extremal) stands, since it is a statement about the
extremal POINT, not a global monotonicity.

## The forward lead (the one genuinely unexplored tool)

The scour found that the dichotomy is currently EMPIRICAL on the near-AP side: HYP-2638 uses Freiman's
3k−4 (small doubling ⟹ short AP ⟹ finite check) for `k=8,9,10`, but the step "high energy ⟹ has an AP
substructure" is done by exhaustive census, not a-priori.  **Balog–Szemerédi–Gowers (energy ⟹ large
small-doubling subset) → Freiman 3k−4 (small doubling ⟹ short AP) is the unused bridge** that would make
the near-AP branch a-priori.  BSG/Green-Tao/PFR appear in the repo only as parked sidecars
(OPEN-QUESTIONS).  That is the concrete additive-combinatorics lead this deep-dive surfaces.

## Ledger

- CORRECTED opus-S179: additive energy is NECESSARY (AP maximizes it) but NOT SUFFICIENT for tightness.
  Counterexamples: translated AP `{1,3,…,25}` (max `E=1469`, loose `L=0.116`) and 2-D GAPs (`E~1000–1225`,
  loose `L~0.09–0.15`, short longest-AP). Matched-energy dim1-vs-dim2 pair splits `L` 0.033 vs 0.091.
- The tight locus is NARROW — exactly the dilates `c·{1..13}` (THM-610) — because `L` is dilation- but not
  translation-invariant (THM-522). The governing object is the resonance-lattice geometry (dimension +
  coherence), not any scalar; `E` is its order-2 shadow.
- PROOF UNAFFECTED: LEM-012/013 key on longest-AP (structural), and LEM-013's `longest-AP≤7` census already
  contains the GAPs. S180 (extremal point) stands.
- Reconciles opus-S113 (necessary-not-sufficient), HYP-2889 (AP-facing not monotone), HYP-2637 (dimension
  penalty), THM-610/kps-S37 (tight locus = dilated APs), opus-S112 (`L` = relation-lattice theta-sum).
- Forward lead: BSG → Freiman 3k−4 to make the near-AP branch a-priori (currently empirical, HYP-2638).
- Files: `lrc14_gap_vs_longestAP_opus_S181` (+out), `lrc14_freiman_dimension_is_the_parameter_opus_S181`
  (+out). -> opus-S179 (corrected)/S180 (stands)/S113, HYP-2637/2889/2638, THM-515/610/522, LEM-012/013.
