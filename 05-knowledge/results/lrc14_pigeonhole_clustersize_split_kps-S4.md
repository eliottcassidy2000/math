# The pigeonhole cluster-size split of the LRC(14) residual (kind-pasteur-2026-06-18-S4)

**Convergence note:** mac-mini-S1 (THM-527 stub) and kind-pasteur-S2/S4 independently reached the
SAME reformulation of the residual: the via-max criterion `W(S\{Vmax}) > 1/(7 Vmax)` is governed,
in the carry-phase limit `V0 → ∞`, by whether the `c = |L|−1` cluster phase-points
`{frac(e_i x)}` (`e_i = Vmax − u_i`, `u_i ∈ L\{Vmax}`) leave a **circular gap > 2/7** at some
`x ∈ G_P` (= kind-pasteur's "offset-fit: circ_width < 5/7"). Density of good `x` = `ρ*(Δ,P)`.

## The pigeonhole bound (clean)

The `c` phase-points divide the circle into `c` gaps summing to 1, so the **max gap ≥ 1/c**
(pigeonhole), giving `circ_width ≤ 1 − 1/c` and the limit via-max margin

> **margin → 7·w_θ ,  w_θ = 6/7 − circ_width ≥ 1/c − 1/7 ,  so  margin ≥ 7/c − 1.**

| `c = |L|−1` | `|L|` | `7/c − 1` | verdict |
|---|---|---|---|
| 1 | 2 | 6 | automatic |
| 2 | 3 | 5/2 | automatic |
| 3 | 4 | 4/3 | **automatic (pigeonhole, ANY offsets)** |
| 4 | 5 | 3/4 | **pigeonhole fails → needs ρ\*>0** |
| 5 | 6 | 2/5 | needs ρ\* |
| ≥6 | ≥7 | ≤1/6 | needs ρ\* |

## The cluster-size split of the residual

- **|L| ≤ 2** (k ≤ 2 large speeds): **PROVED** — S1 (k≤1) + the k=2 slice (HYP-2581a).
- **|L| ∈ {3,4}**: cluster' = `c ≤ 3` points ⟹ pigeonhole forces limit margin ≥ 4/3 > 1 for ANY
  offsets. Closes **exactly like the k=2 slice**: via-max criterion for large Vmax (the limit
  regime) + a **finite core** for small Vmax (criterion-failures live at Vmax ≤ ~22, large speeds
  barely above 13). The remaining rigor is the carry-phase finite-V0 correction (mac-mini THM-527
  machinery), NOT a new Diophantine obstruction.
- **|L| ≥ 5** (cluster' `c ≥ 4`): pigeonhole insufficient (`7/c−1 ≤ 3/4 < 1`). This is the GENUINE
  hard case = the uniform `ρ*(Δ,P) ≥ c0 > 0` floor (mac-mini THM-527 / HYP-2581d / OPEN-Q-108).

## Supporting facts (VERIFIED, kps-S4, exact)

- **k≥3 floor:** min M over ~2100 sampled k≥3 covering S3 sets (Vmax≤130) = **2/21 ≈ 0.0952**
  (M·14 ≈ 1.33), LOOSER than the global covering floor 2/23 (k≤2). Tightness is ALWAYS driven by
  an adjacent SMALL pair `{a,a+1}` (e.g. {10,11}, {8,13}) at a small crossing `k/(a+b)`; the
  cluster never binds tightly. So the genuinely-tight covering sets are all k≤2 (PROVED).
- **Criterion C (free v) holds on 0 failures / ~8000 covering S3 sets with |L|≥3** — only |L|=2
  has C-failures (S*, MISTAKE-076). For |L|=3,4 the via-max criterion holds for Vmax≥~60; fails at
  small Vmax (finite core).
- **7-adic mechanism:** at τ=k/7 (gcd(k,7)=1), a runner is dangerous iff 7|v; covering forces a
  multiple of 7, so τ=k/7 is never a global witness — the covering-forced mult-of-7 cluster runner
  is the binding large runner whenever a large speed binds (mod7=0 in every observed case), but it
  clears 1/14 comfortably (M ≈ 0.10+). Ties to codex's private-q=7 obligation (HYP-2579).

Scripts: `04-computation/lrc14_{slowfast_offsetfit,k3_minM,L34_criterion_threshold}_kps-S4*.py`.
→ THM-526, THM-527 (mac-mini, pending), HYP-2581d, OPEN-Q-108.
