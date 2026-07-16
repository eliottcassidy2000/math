---
id: THM-886
title: THE SELECTOR RE-COORDINATIZATION VERDICT — the two-sheet rows through the Verblunsky α-map (THM-883): (i) THE TWINS SURVIVE THE MAP AT DEPTH 10: THM-789's fingerprint pair U₀+(13,9) / U₀+(17,13) — identical raw signed-component fingerprints — have NEARLY identical α-vectors too (HUG 0.544 vs 0.547; ALT 0.78 both; profiles matching to ~0.01–0.06): the α-differences are nonzero but small, so the raw-safe-measure spectral coordinate does not crack what fingerprints could not; the twins' true difference (return-thickened component incidence, per THM-789) names the next coordinate: THE RETURN-THICKENED α (OPUC of the E+R measure); (ii) HUG READS ATOMICITY, ORTHOGONALLY TO M: the THM-824 diagnostic row U₀+(13,5) is the most boundary-hugging (HUG 0.811: α₉ ≈ −0.81, α₃ ≈ −0.56 — consistent with its all-fail radius-budget geometry) while the tightest-M control (AP-core, M = 1/11) is the most diffuse (HUG 0.25, μ(G) = 0.0998 over 28 arcs): the α-map coordinatizes safe-set GEOMETRY (atomic vs diffuse) as an axis independent of the maximin; (iii) the alternation statistic (ALT ≈ 0.78) is shared across deep rows and controls — period-2-ness alone does not select
status: COMPUTED VERDICT (6-row battery, exact M + arcs + α to depth 10; the hypothesis "deep ⟹ period-2 boundary-hugging" is REFINED: hugging = atomicity (true for the 824 row, moderate for the twins), alternation non-selective, and the twin-separation requires the return-thickened measure — an honest partial with the next coordinate identified)
source: opus-2026-07-16-S323 (owner: run the selector re-coordinatization through the alpha-map)
depends_on:
  - THM-883 (the α-map)
  - THM-789 (U₀, the trapped liar, the fingerprint twins, the return set (−1/858, 1/858))
related: [THM-824 (the (13,5) diagnostic), THM-797/803/817 (the selector arc this coordinatizes), THM-840 (the Markov question — same lesson: raw compressed states lose the twins)]
verification: 05-knowledge/results/selector_alpha_recoordinatization_opus_S323.out
---

# THM-886 — the selector re-coordinatization verdict

Battery: S = 2U ∪ {x, y} for U₀ = {1,2,3,5,7,8,9,10,11,12} with the canon
pairs (13,9), (17,13), (13,5), plus controls (odd-pair, perturbed-core,
AP-core). For each: exact maximin, safe arcs at δ = 1/14, α₀..α₉,
HUG = max_{n≥2}|α_n|, ALT = period-2 sign-alternation score.

| row | M | arcs | μ(G) | HUG | ALT |
|---|---|---|---|---|---|
| U₀+(13,9) trapped liar | 2/19 | 20 | 0.044 | 0.544 | 0.78 |
| U₀+(17,13) fingerprint twin | 2/19 | 16 | 0.038 | 0.547 | 0.78 |
| U₀+(13,5) THM-824 pair | 2/17 | 16 | 0.047 | **0.811** | 0.11 |
| U₀+(15,9) control | 2/17 | 28 | 0.061 | 0.453 | 0.56 |
| perturbed-core control | 1/9 | 24 | 0.057 | 0.501 | 0.78 |
| AP-core control | 1/11 | 28 | 0.100 | 0.250 | 0.78 |

**Findings.**
1. **The twins survive.** The (13,9)/(17,13) pair — built to defeat raw
   fingerprints — also nearly defeats the α-map: profiles agree to
   0.01–0.06 through depth 10. The difference is nonzero (the map is not
   blind) but far below selection threshold. Their true separation lives in
   the return-thickened incidence (THM-789), so the next coordinate is the
   **return-thickened α**: OPUC of the measure on E + closure(R) rather
   than the raw safe set. (Same lesson as THM-840's endpoint-Markov
   negative: compressed states that ignore the return structure lose
   exactly these rows.)
2. **HUG is an atomicity axis, orthogonal to M.** The most boundary-hugging
   row is the THM-824 diagnostic (0.811) whose 214-core replay all-failed
   the radius budget — its safe geometry is nearly-atomic; the tightest-M
   control is the most diffuse. The α-map contributes a genuinely new
   coordinate (atomic ↔ diffuse) rather than re-measuring the maximin.
3. **Alternation alone does not select** (0.78 across deep rows and
   controls alike): the two-band signature needs the hug AND the context.

**Named next.** The return-thickened α (E + R measures for the THM-789/817
cell structures); depth-20 tails for the twins (does the separation grow?);
the α-map along the THM-817 satellite ray B_n = 506 + 360360n (does N_R =
Θ(B) show up as α-tail period growth = torus dimension?).
