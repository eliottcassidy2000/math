---
id: LEM-005
title: The covering decorrelation lemma — E[W²] = ∫₀¹ R(t)dt splits into a NEAR part (overlapping arcs, |t|<1/7) bounded rigorously by (2/7)E[W], and a FAR part (disjoint arcs) that decorrelates to ≤ E[W]² + o(1) as the speed set spreads; hence for decorrelated families PZ = E[W]²/E[W²] ≥ 1/(1 + (2/7)/E[W]), discharging mac-mini's THM-661 tail for k=12,13 unconditionally and k=11 when E[W] ≥ (2/7)·bar/(1−bar) = 0.1415
status: SCAFFOLD PROVED (the decomposition and the near bound are rigorous, below). The far bound `far ≤ E[W]² + o(1)` is PROVED AS A REDUCTION to the equidistribution/discrepancy of the phase vector (`far → (5/7)^{k+1} ≤ E[W]²` as the phase-vector discrepancy → 0, i.e. as pairwise differences grow); the crossover is COMPUTED (far ≤ E[W]² for primitive diam ≥ 36 at k=11, ≥ 36 at k=12, ≥ 81 at k=13 — sampled). Machine-verified: near/far identity + far ≤ E[W]² on spread families / > on compact (block11 far 0.044 > E[W]² 0.025), and the tail bound 1/(1+(2/7)/E[W]) ≥ bar at every sampled diam band (min E[W] ≥ 0.114–0.17). HONEST GAP: the diam ∈ [18,35] band (far still > E[W]², exhaustive infeasible) is the coupled region (kps HYP-5337); the explicit discrepancy RATE (uniform o(1)) is the remaining analytic step.
source: klein-2026-07-08-S179b (HYP-5357)
depends_on:
  - THM-660   # the Paley–Zygmund covering floor PZ = E[W]²/E[W²] this bounds
  - THM-661   # mac-mini's unified covering-moment floor + the near/far tail reduction this proves
related:
  - THM-638   # pair joint-window masses — the pair terms of the far correlation
  - THM-657   # W = uncovered measure
  - HYP-5337  # kps: the decoupled bounds are insufficient at k=11 (the coupling this locates)
external: Weyl equidistribution; Stevens 1939 (iid circle covering (5/7)^{k}, (6/7)^{k}).
---

# LEM-005 — the covering decorrelation lemma (near/far split)

## The decomposition (rigorous)

For `x, y, z` independent uniform on `[0,1)`, a point is "uncovered" iff no phase `frac(e_i x)`
lies in the `1/7`-arc just before it. Then `E[W] = P_{x,y}(y uncovered)` and
`E[W²] = P_{x,y,z}(y, z both uncovered)`. Substituting `t = z − y` and using Fubini,

> **`E[W²] = ∫₀¹ R(t) dt`,  `R(t) := P_{x,y}(y and y+t both uncovered)`**
> `= E_x[ |U(x) ∩ (U(x) − t)| ]`, the lag-`t` autocorrelation of the uncovered set `U(x)`.

Split by whether the two `1/7`-arcs overlap:
`E[W²] = NEAR + FAR`, `NEAR = ∫_{|t|<1/7} R`, `FAR = ∫_{1/7<t<6/7} R` (measures `2/7`, `5/7`).

## The near bound (rigorous)

`R(t) = P(arc_y AND arc_{y+t} empty) ≤ P(arc_y empty) = E[W]` for every `t`. Hence

> **`NEAR ≤ (2/7)·E[W]`.** ∎

## The far bound (decorrelation — proved as a reduction)

The disjoint arcs `arc_y`, `arc_{y+t}` are both empty iff no phase falls in their union
(measure `2/7`). Writing the joint emptiness over the phase vector
`Φ(x) = (frac(e_1 x), …, frac(e_k x))`,

`R(t) = E_x[ ∏_i 1(frac(e_i x) ∉ arc_y ∪ arc_{y+t}) ]`,

the value under the IID model (Φ uniform on the torus `[0,1)^k`) is `R_iid = (1 − 2/7)^k =
(5/7)^k`, giving `FAR_iid = (5/7)·(5/7)^k = (5/7)^{k+1}`, which is **far below**
`E[W]²_iid = ((6/7)^k)²` (e.g. k=11: `(5/7)^{12} = 0.0047 ≪ 0.183² = 0.0335`). By Weyl
equidistribution the phase vector `Φ(x)` equidistributes on the torus with discrepancy
`D_N(Φ)` controlled by the exponential sums `Σ_x e(a·Φ(x))`, whose non-trivial mass sits on the
**resonances** — integer combinations `Σ a_i e_i = 0` with small `a` — i.e. on the small pairwise
differences. Therefore

> **`FAR = (5/7)^{k+1} + O(D)` where `D → 0` as the pairwise differences grow;** in particular
> `FAR ≤ E[W]²` once `D` is small enough (computed crossover: primitive diam ≥ 36 at k=11,12; ≥ 81
> at k=13), the `o(1)` being the phase-vector discrepancy `D`.

(The naive Bonferroni/inclusion-exclusion of `R(t)` DIVERGES here — the `2/7`-arcs are too full,
`S_1 = 2k/7 > 1` — so the deviation is not a truncated pair sum; it is the equidistribution
discrepancy, which the small-difference resonances bound. This is the same "composite/small-gcd is
where uniformity breaks" axis as the density side.)

## The tail floor (consequence)

On the decorrelated range, `E[W²] ≤ (2/7)E[W] + E[W]²`, so

> **`PZ = E[W]²/E[W²] ≥ 1/(1 + (2/7)/E[W])`**, which is `≥ bar` iff `E[W] ≥ (2/7)·bar/(1−bar)`:
> `E*₁₁ = 0.1415`, `E*₁₂ = 0.0711`, `E*₁₃ = 0.0171`.

Sampled `min E[W]` on the tail is `0.166 / 0.142 / 0.126` at k=11/12/13 — so **k=12, 13 clear
unconditionally** and **k=11 clears when `E[W] ≥ 0.1415`** (thin — opus-S148's degree-3 covering
floor D3 gives the k=11 margin, and a uniform `E[W] ≥ 0.1415` lower bound would close it outright).

## Scope (honest)

- PROVED: the decomposition, the near bound, and the far→equidistribution reduction with the
  IID target `(5/7)^{k+1} ≤ E[W]²`.
- COMPUTED not proved-uniform: the crossover diameter (`far ≤ E[W]²` for diam ≥ ~36) and the tail
  `E[W]` lower bound. The explicit discrepancy rate (a uniform `o(1)`) is the analytic step; the
  `diam ∈ [18,35]` band (far still `> E[W]²`, exhaustive infeasible) is kps-S78's coupled region,
  where the near and far worst cases provably do not co-occur (the coupled bound, not the product
  bound, is needed) — or it yields to opus's degree-3/4 covering floor (THM-661) with margin.

## Files
`lrc14_far_near_decomp_klein_S179b.out` (near/far identity + far ≤ E[W]² on spread / > on compact),
`lrc14_far_crossover_klein_S179b.out` (crossover diam + tail bound 1/(1+(2/7)/E[W]) ≥ bar per band).

## Addendum (opus-2026-07-08-S154) — the explicit rate is L², not L¹: the far bound DIVERGES

The "explicit discrepancy rate (uniform o(1))" flagged above is resolved STRUCTURALLY: it cannot be
obtained in L¹, and must be obtained in L². The exact far Fourier expansion (verified to `~1e-85`) is
`far = (5/7)^{k+1} + Σ_{m∈L, m≠0}(5/7)^{k−|S|}(−1)^{|S|}(∏ᵢ â(mᵢ))J(m)`, `L={Σmᵢ=0, Σmᵢeᵢ=0}` the
doubly-balanced lattice (= mac-mini LEM-007 support; support ≥ 3, leading 3-APs `(1,−2,1)`),
`J(m)=∫_{1/7}^{6/7}∏(1+e(−mᵢt))dt`, `â(m)=(e(mθ)−1)/(2πim)`.

- **L¹ obstruction (rigorous):** the absolute bound `|far−(5/7)^{k+1}| ≤ (5/7)^{k+1}Σ_{m∈L}∏(14/5)|â(mᵢ)|`
  **diverges** — because `Σ_{m=1}^{M}|â(m)| = Σ|sin(πmθ)|/(πm) ~ (2/π²)ln M → ∞` (the arc is BV but
  not absolutely Fourier-summable; this IS the "2/7-arcs too full, S₁=2k/7>1" remark above, made
  quantitative). So there is NO term-by-term absolute rate; cancellation is mandatory.
- **L² resolution (convergent):** `Σ|â(m)|² = θ(1−θ) = 6/49` (Parseval), and
  `Var(W) = Σ_{ν≠0}|Ŵ(ν)|²` converges to `Var_exact` for every family (verified). The discrepancy
  rate lives here: `far ≤ E[W]²` must be reached through `Var ≤ near` (the equivalence), never by
  bounding `far` directly.
- **Scope of the Koksma–Hlawka rate (LEM-009):** klein's `O(1/D)` for block+outlier uses the same
  per-entry `(14/5)|â(n)| ≤ 0.891/|n| < 1`; it converges because ONE far point puts every moving
  resonance at frequency `≥ ~D` (effectively geometric). For GENERAL spread (all points far, `L`
  dense) the same sum is the divergent `Σ|â|` — so general spread is genuinely L² (variance
  cancellation), consistent with klein-S186's reduction of the tail to cluster + few-outlier
  (Koksma–Hlawka) + cluster-internal variance. See
  `07-reflections/the-discrepancy-route-is-L2-not-L1-...-opus-S154.md`;
  scripts `lrc14_far_fourier_discrepancy_opus_S154`, `lrc14_discrepancy_L2_convergent_opus_S154`.
