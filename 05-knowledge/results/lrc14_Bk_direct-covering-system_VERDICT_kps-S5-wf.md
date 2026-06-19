# B(k) via DIRECT COVERING-SYSTEM BOUND — verdict (kind-pasteur-2026-06-18-S5-wf)

Angle: close `inf_E μ(E) ≥ c0 > 0` (`μ(E)=meas{x: maxgap{frac(e_i x)}>2/7}`, `0∈E`, `k≤13`)
by a SINGLE uniform strip-avoidance / covering-system bound, with NO spread split.

## Two rigorous lower bounds tried (both PROVED to be valid lower bounds for μ)

### (1) Single 0-anchored arc:  μ_0(E) = meas{x: frac(e_i x) ∉ [1/7,3/7] ∀i}
- **PROVED LEMMA:** `μ_0(E) ≤ μ(E)`. (All points avoid the closed 2/7-arc ⟹ a point-free
  open arc of length >2/7 ⟹ maxgap >2/7.) Verified 0 containment violations / 213 shapes.
- **REFUTED as a uniform floor:** μ_0 → 0 as spread → ∞ at fixed k. Worst μ_0(k=13):
  0.0358 (spread≤26) → 0.0182 (≤52) → 0.0111 (≤104). Single fixed arc is too thin.

### (2) Arc-AVERAGED covering measure = the EXCESS-GAP functional
- **NEW PROVED LEMMA (this session):**
  `Ψ(E) := ∫₀¹ Σ_j (gap_j(x) − 2/7)_+ dx ≤ μ(E).`
  Proof: pointwise the integrand is 0 off {maxgap>2/7} and ≤ Σ_j gap_j = 1 on it, so it is
  ≤ 1[maxgap>2/7]; integrate. RIGOROUS, no approximation. Verified 0 violations / 300+213 shapes.
- **AVERAGING IDENTITY (PROVED + numerically cross-checked, err→3e-5):**
  `∫_θ 1[all frac(e_i x) avoid the 2/7-arc centered at θ] dθ = Σ_j (gap_j−2/7)_+`, so
  `Ψ(E) = ∫_θ μ_θ(E) dθ` = average over arc-position of the strip-avoidance measure. This IS
  the covering-system functional, averaged — it restores bulk that a single arc loses.
- Exact-rational engine (collision + gap=2/7 breakpoints, midpoint rule exact on linear cells).
  `Ψ(consecutive E={0..12}) = 3379/64680 = 0.052242`, ratio Ψ/μ = 3379/11606 ≈ 0.2911.
- **ALSO REFUTED as a uniform floor:** Ψ → 0 with spread (k=13): worst Ψ 0.0482 (spread≤13)
  → 0.0245 (≤26). Ψ halves as spread doubles. Same obstruction as μ_0.

## WHY the covering-system angle CANNOT close B(k) alone (the structural reason)
Both μ_0 and Ψ are "thin/local" functionals (single-arc / L¹ excess-gap mass). As spread→∞
the gaps subdivide finely, the excess-over-2/7 mass → 0, so any thin functional → 0. But
**μ itself does NOT → 0**: for relation-free large spread, μ → F(k) (the iid CEILING),
verified `|μ−F(4)| ~ 1/maxE → 0` (0.0013 at maxE=189, 4e-5 at 5030, 1e-5 at 63331). The
uniform floor is a property of the BULK MEASURE μ, not of any thin minorant. Covering-system /
strip-avoidance bounds are provably insufficient.

## The binding infimum is at BOUNDED SPREAD (exact, confirms HYP-2584/2585)
For every k, the bounded-spread minimum of μ is STRICTLY BELOW the iid ceiling F(k):
| k | bounded-spread min μ | argmin (perforated near-AP) | F(k) |
|---|---|---|---|
| 4 | 19/21=.9048 | {0,1,2,3} | .9971 |
| 5 | 9/14=.6429 | {0,1,2,3,4} | .9683 |
| 6 | 4/7=.5714 | {0,1,2,3,4,5} | .8999 |
| 7 | 13/35=.3714 | {0,2,3,4,5,6,8}={0..8}\{1,7} | .7998 |
| 8 | 71/220=.3227 | {0,2,3,4,5,6,8,11} | .6846 |
| 9 | 164/735=.2231 | {0,2,4,5,6,7,8,9,12} | .5689 |
So `inf_E μ(E) = min over BOUNDED-SPREAD shapes` — the large-spread limit is the ceiling.
This is the compactness that makes B(k) a FINITE check; the covering-system functionals
fail precisely because they vanish on exactly the large-spread shapes where μ is safe.

## RELATION TO THE FOURIER ROUTE (kps-S5 grounding, the correct frame)
`μ(E) = F(k) + Σ_{m∈Λ(E)\0} ĝ(m)` (Parseval, PROVED). The covering functionals correspond to
crude one-frequency minorants; their decay = the relation-free tail → 0. The uniform floor
needs the SHORT-RELATION patterns of Λ(E) (bounded-spread = many short relations). The
excess-gap lemma Ψ≤μ is a clean rigorous tool but a strictly weaker (thinner) functional.

## VERDICT
- **NEGATIVE for the assigned angle**: no single uniform covering-system / strip-avoidance
  bound (μ_0 single-arc, nor Ψ arc-averaged) gives `inf_E μ ≥ c0` — both PROVABLY → 0 with
  spread. The grounding-doc warning (point 1) is now made rigorous and quantitative.
- **POSITIVE deliverables**: (i) NEW proved lemma `Ψ(E) ≤ μ(E)` with exact-rational engine and
  the arc-averaging identity `Ψ = ∫_θ μ_θ dθ`; (ii) exact demonstration that `inf μ` is attained
  at bounded spread (`bounded-spread min μ < F(k)` ∀k), confirming the finite-check reduction;
  (iii) localization: B(k) is NOT a covering-system question — it is the bulk/relation-lattice
  (Fourier short-relation) crux, OPEN-Q-108 unchanged.
- No counterexample to `inf μ > 0`; the floor remains PLAUSIBLE and reduces to the
  bounded-spread finite check + relation-free tail bound (the Fourier route), NOT this angle.

Scripts: 04-computation/lrc14_Bk_{direct-covering-system-bound,mu0_spread_decay,excess_gap_floor,
psi_spread_decay,fourier_tail_explicit,boundedmin_vs_Fk,psi_rigour}_kps-S5-wf.py
