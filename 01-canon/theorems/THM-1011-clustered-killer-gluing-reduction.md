---
id: THM-1011  # renumbered from 1010 (codex sharp-descent-recursion first-pushed 1010 at 08:22 vs 09:43)
title: THE CLUSTERED-KILLER REDUCTION — the last stratum of "covering ⟹ M > 1/14" reduced to an explicit computable criterion, with the crude form REFUTED, the sharp form largely closing, and the true obstruction identified (and one hypothesis honestly refuted). SETTING: S = P ∪ K, core P slow (|P| = 13−j), killer block K fast (min K > 13·max P) but internally CLUSTERED, so THM-1007's balance iteration fails (its telescoping needs each killer 13× the previous). REDUCTION: apply THM-933's two-block gluing at ρ = 1/14. (I) THE CRITERION: with d_P ≤ δ(S_ρ(P)), d_K ≤ δ(S_ρ(K)), q(K) the block discrepancy and K_P the number of circular components of the core safe set, the sharp form (BG-K) gives μ(safe) ≥ d_P·d_K − q(K)·K_P, so the clustered case closes exactly when **q(K)·K_P < d_P·d_K** — fully computable in exact rationals; (II) THE CRUDE FORM IS HOPELESS: (BG) with the universal cap K_P → M(P) = Σ(core speeds) ≈ 66 closes 0/14 sampled clustered configs, the implied gap requirement being ~800× rather than the actual 13× — the component-count refinement is NOT optional; (III) THE SHARP FORM LARGELY CLOSES: exact component counts give K_P = 14–18 (4.5× smaller than the crude cap) and (BG-K) closes **9/14 configs at the actual gap 13×**, rising to 12/12 at gap 60× (scan: 9/12 at 13, 9/12 at 20, 11/12 at 26, 11/12 at 40, 12/12 at 60 — NON-monotone, so the obstruction is not scale alone); (IV) THE OBSTRUCTION IDENTIFIED: failures are exactly the blocks with large q(K), and q(K) is large precisely when the killer speeds are NEARLY EQUAL to one another (q = 1.45×10⁻² at K = [257,258,263], spread 2.3%; q = 6.6×10⁻⁴ at K = [417,544], spread 30%) — near-equal killers have nearly coincident bad sets, so the block safe set has long runs and large H-oscillation, while well-separated killers interleave and equidistribute; (V) HONEST REFUTATION: the hypothesis "large q(K) = internal arithmetic resonance (small reduced ratio), routing to the structured strata" is FALSE — median smallest-ratio-sum is 524 in the low-q half vs 515 in the high-q half, statistically identical, gcd 1 in both. The obstruction is multiplicative tightness, not resonance. STATUS: the clustered stratum is NOT closed; it is reduced to a sharp computable criterion, with the residual confined to tightly-clustered killer blocks at gaps below ~60×
status: (I) reduction exact (THM-933 (BG-K) instantiated at m = 2, all inputs exact rationals); (II)(III) measured (14-config battery + 5-point gap scan, exact δ, q, component counts); (IV) obstruction characterized empirically; (V) the resonance hypothesis REFUTED by its own test (the script's concluding lines were written before the data and are corrected in the .out — recorded as a hypothesis-refutation, not a finding). The clustered case remains OPEN; nothing here is claimed as its closure
source: kind-pasteur-2026-07-18-S128 (cont.53; owner: now close the clustered killer case)
depends_on:
  - THM-933   # the two-block gluing (BG)/(BG-K) this instantiates
  - THM-1007  # the single-killer + lacunary-chain closure whose gap this stratum is
  - LRC(≤13)  # core and block LRC floors
related:
  - THM-995 (XI)(XII) — the weak-target localization; THM-726 — the sharp-target multi-killer half
  - THM-935 — the relation-content axis (the resonance route tested and refuted here)
scripts: 04-computation/clustered_killer_gluing_kps_S128c53.py, clustered_bgk_kps_S128c53.py, clustered_gap_threshold_kps_S128c53.py (+ .out, incl. the corrected resonance test)
---

# THM-1011 — the clustered-killer reduction (renumbered from 1010)

## The criterion

THM-933 (BG-K) at m = 2 with blocks P (core) and K (killers):

> μ(S_ρ(P) ∩ S_ρ(K)) ≥ d_P·d_K − q(K)·K_P,  ρ = 1/14,

so **the clustered case closes iff q(K)·K_P < d_P·d_K** — every input exactly computable
(δ by interval sweep, q by the piecewise-linear H-oscillation, K_P by component count).

## What the battery says

| form | cap used | closes at gap 13× |
|---|---|---|
| (BG) crude | M(P) = Σ core ≈ 66 | **0/14** (needs ~800× gap) |
| (BG-K) sharp | K_P = 14–18 exact | **9/14** |

Gap scan for (BG-K): 9/12 at 13×, 9/12 at 20×, 11/12 at 26×, 11/12 at 40×, **12/12 at 60×**
— non-monotone, so scale alone is not the obstruction.

## The obstruction (and a refuted hypothesis)

Failures are exactly the large-q(K) blocks, and q(K) is large when the killers are NEARLY
EQUAL — [257,258,263] (2.3% spread) gives q = 1.45×10⁻², versus [417,544] (30% spread) at
6.6×10⁻⁴. Near-equal killers have nearly coincident bad sets ⟹ long runs in the block safe
set ⟹ large H-oscillation.

**Refuted:** the natural guess that large q(K) means internal arithmetic resonance (a small
reduced ratio) — hence routing to THM-935's structured strata — is FALSE: median
smallest-ratio-sum 524 (low-q) vs 515 (high-q), gcd 1 in both. Multiplicative tightness,
not resonance.


## (VI) THE BEAT LAW — and the REFUTATION of the effective-LRC(13) collapse (cont.54)

THM-1011's named-next proposed that near-equal killers "contribute essentially ONE
constraint", suggesting an effective-LRC(13) collapse. **The data refutes that as a MEASURE
statement, and replaces it with a sharper discrepancy law.**

**Kernel fact (trivial, verified).** If ‖v₁t‖ < λ and ‖v₂t‖ < λ then ‖dt‖ ≤ ‖v₁t‖+‖v₂t‖ < 2λ,
d = v₂−v₁ — two near-equal killers can fail together only when the BEAT dt is near an
integer. (Containment checked exactly on grids for d = 1,3,7,15,40.)

**(a) NO measure collapse — the effective-LRC(13) route is FALSE.** For K = {v, v+d} the
safe measure is δ(K) = 0.7347 for EVERY d tested (1…144) — exactly 5/7 + 1/49. The pair
always costs 2λ·2 − 1/49 of measure no matter how close the killers are; the overlap is the
generic 4λ² = 1/49 and does not grow as d → 1. A near-equal pair does NOT become one
constraint in measure, so no LRC(13) collapse is available.

**(b) THE BEAT LAW — the collapse is a DISCREPANCY phenomenon.** What d controls is q(K):

| d | 1 | 2 | 3 | 5 | 8 | 13 | 21 | 55 | 144 |
|---|---|---|---|---|---|---|---|---|---|
| q(K), v=200 | 1.60e-2 | 8.46e-3 | 5.96e-3 | 3.96e-3 | 2.79e-3 | 2.13e-3 | 1.62e-3 | 1.22e-3 | 9.9e-4 |

and q is essentially INDEPENDENT of v: at d = 1, q = 1.597e-2 (v = 200) vs 1.548e-2
(v = 400) while q·v changes by 2×. Empirically **q(K) ≈ c/d_min with c ≈ 0.016**, d_min the
smallest internal difference. The beat, not the speed scale, sets the block's discrepancy.

**(c) THE SHARPENED CRITERION AND A THIN RESIDUAL.** Substituting into (I):

> (BG-K) closes as soon as **d_min > c·K_P/(d_P·d_K)**, which for the observed
> K_P ≈ 15, d_P·d_K ≈ 0.05, c ≈ 0.016 gives **d_min ≳ 5**.

This matches the battery: the (BG-K) failures all have inferred d_min ≈ 2–5, the successes
d_min ≈ 6–12. So the clustered residual is not "clustered killers" broadly — it is exactly
**killer blocks containing two speeds within about 5 of each other**, a very thin stratum
(and there the beat speed d ≤ 5 is itself core-sized, which is the structure a finite/
covering argument should exploit).

**Status change:** the named-next "effective-LRC(13) collapse" is WITHDRAWN (refuted in (a)).
The replacement route is the d_min ≲ 5 stratum with its core-sized beat speed.

## Named next
- ~~The near-equal/effective-LRC(13) route~~ — REFUTED in (VI)(a): no measure collapse.
- **The d_min ≲ 5 stratum:** two killers within ~5 of each other, where the beat speed d is
  itself core-sized — a bounded-d finite analysis or a covering/counting exclusion is the
  live route to closing exactly what (BG-K) misses.
- Sharper q(K) bounds for tight blocks, or a better core component bound K_P.
