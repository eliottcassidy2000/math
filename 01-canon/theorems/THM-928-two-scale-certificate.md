---
id: THM-928
title: THE TWO-SCALE CERTIFICATE THEOREM — (A) the lacunary cascade: consecutive ratios ≥ 15 give uncovered ≥ (6/7)^13 − 2/R > 0 (LRC(14) for all 15-lacunary packets, 15-line proof; uniform version R ≥ 19 proves LRC(n) for EVERY n; fixed-point witnesses for dilate families); (B) the Bézout-rotation orbit law + Denjoy–Koksma: an INDEPENDENT SECOND PROOF of T2 (primary: mac-mini THM-926) — W's component edges march on the two cofactor rotations B*/A and B̃/B, the deviation = difference of two DK sums, V₁ = 44/(169C), verified 108/108; (C) THE FIRST ACCESSIBLE-SCALE FULL-PROBLEM LEVEL-5 CERTIFICATE: the corrected weighted-capped filter finds {300, 406, 511, 652, 862, 963, 1074, 1357, 1459, 1571, 1776, 1991, 2208} with BONF5 ≥ +0.039131 > 0 (exact rational), actual uncovered 0.117744
status: (A) PROVED (proof in-file; ingredients verified exactly, see script note); (B) VERIFIED 108/108 (identity exact everywhere, DK every family, bound covers everywhere) — complementary to THM-926, does not override it; (C) FOUND & CERTIFIED exact-rationally — the S329 existence question resolves POSITIVE behind the honest filter
source: opus-2026-07-16-S332 (owner: work the unified target, pull frequently, many concrete steps)
depends_on: [THM-926 (T2 primary, mac-mini), THM-864 (beat localization; the weight (q+p−1)), THM-927 (the third blocker), MISTAKE-151 (the two filter loopholes this corrects), THM-924 (trichotomy: resonant = exactly periodic)]
scripts: 04-computation/cascade_theorem_opus_S332.py, t2_structure_recon_opus_S332.py, t2_denjoy_koksma_opus_S332.py, bothclean_v2_weighted_opus_S332.py → 05-knowledge/results/*_opus_S332.out
---

# THM-928 — the two-scale certificate theorem

The loose-family certificate program splits by scale. This file proves the
exponential half outright, gives a second, structurally different proof of
the certificate half's key lemma (T2), and delivers the first coercive
accessible-scale certificate behind the corrected filter.

## (A) THE CASCADE THEOREM (exponential half) — PROVED

**Theorem A.** Let x₁ < x₂ < … < x₁₃ be distinct positive integers with
x_{k+1}/x_k ≥ R for all k. With D_x = {t ∈ ℝ/ℤ : ‖xt‖ < 1/14},

> μ( ℝ/ℤ ∖ ∪ₖ D_{x_k} ) ≥ (6/7)¹³ − 2/R.

Since (6/7)¹³ = 13060694016/96889010407 ≈ 0.134801, **any R ≥ 15 gives
positive uncovered measure, hence LRC(14) holds for every 15-lacunary
packet** — with ≥ 0.03 to spare at R ≥ 20. (Verified: per-step lemma and
component growth exact on 60 random configs, zero violations; κ_k ≤ 2x_k
along real cascades; theorem bound respected by 5-step exact prefixes at
R = 15/20/30 and by FULL exact cascades at n = 5, 6, 7.)

*Proof.* W₀ = ℝ/ℤ; W_k = W_{k−1} ∖ D_{x_k}; κ_k = #components(W_k);
μ_k = μ(W_k).

(L1) *Per-step loss.* D_x is x teeth of width 2/(14x) = 1/(7x), one per
cell [j/x, (j+1)/x). A component of length ℓ meets ≤ ⌊xℓ⌋ + 1 ≤ xℓ + 1
cells, each holding one tooth of measure 1/(7x), so
μ(W ∩ D_x) ≤ (1/7)μ(W) + κ(W)/(7x), i.e.
μ_k ≥ (6/7)μ_{k−1} − κ_{k−1}/(7x_k).

(L2) *Component growth.* Removing one tooth splits ≤ 1 component; the
number of teeth meeting W is ≤ x·μ(W) + κ(W). So
κ_k ≤ 2κ_{k−1} + x_k·μ_{k−1} ≤ 2κ_{k−1} + x_k.

(L3) *Induction* κ_k ≤ 2x_k for R ≥ 4: base κ₁ = x₁ (x₁ teeth cut the
circle into x₁ arcs); step κ_k ≤ 4x_{k−1} + x_k ≤ (4/R + 1)x_k ≤ 2x_k.

(L4) *Unroll.* Step 1 is exact (μ₁ = 6/7). For k ≥ 2,
κ_{k−1}/x_k ≤ 2x_{k−1}/x_k ≤ 2/R, so
μ₁₃ ≥ (6/7)¹³ − (2/(7R))·Σ_{j≥0}(6/7)ʲ = (6/7)¹³ − 2/R. ∎

**(A′) THE NESTED-GAP SHARPENING (S336): R ≥ 7/3 suffices — and it is tight.**
A far better constant by a far simpler proof. Maintain a closed interval
[a, b] inside a SAFE GAP of every processed comb (‖w t‖ ≥ 1/14 on all of
[a, b]). One step: if w·(b − a) ≥ 2, the period [⌈wa⌉, ⌈wa⌉ + 1] fits inside
[wa, wb], and its gap pulls back to a subinterval of length (6/7)/w. Since
2/(6/7) = 7/3, consecutive ratios ≥ 7/3 sustain the recursion forever; any
point of the final interval is a witness. No measure theory, no components —
pure floor arithmetic, and the witness is RATIONAL with explicit denominator.
TIGHT: boundary chains at ratio exactly 7/3 land min distance EXACTLY 1/14
(verified exactly, nested_gap_verify_opus_S336; 200 random ratio-3 chains
pass). Uniform-n version: gap 1/n, step length (1 − 2/n)/w, constant
2/(1 − 2/n) = 2n/(n − 2) → 2: ratio ≥ 3 proves LRC(n) for ALL n ≥ 5
simultaneously (n = 3: ratio 6; n = 4: ratio 4). LEAN: LRCLacunaryNest.lean
(lonely_of_pos_lacunary — the first universal quantifier-level branch,
reusing arcSafe/norm_ge_of_arcSafe). The measure cascade (A) retains value
where the INTERVAL structure is unavailable (measure-only block data); for
pure lacunarity (A′) supersedes it.

**Uniform-n corollary.** The same algebra at gap 1/n (tooth width 2/(nx),
per-component loss 2/(nx), geometric sum ≤ n/2) gives, for any n:
uncovered ≥ (1 − 2/n)^{n−1} − 2/R. The prefactor's minimum over n is 1/9
(n = 3), so **R ≥ 19 proves the Lonely Runner Conjecture for every
R-lacunary speed set at every n simultaneously.** (Related: Dubickas 2011
gets ratio 1 + O(log n/n), asymptotically far stronger; the value here is
the uniform explicit constant, the 15-line self-contained proof, and exact
constants that feed the certificate program.)

**Witness lemma (dilate families).** For {c·qʲ}, j = 0..12, integer q ≥ 2:
t = a/(qˢ − 1) has multiply-by-q orbit of length ≤ s in ℤ/(qˢ − 1), so all
13 distances take ≤ s values; s ≤ 2 suffices for every q ≤ 20 (verified
exactly; q = 14: t = 6/13, all distances 6/13; q = 2: t = 1/3, distances
1/3). THM-927's dilate escape family is lonely by inspection — the
parallel-class walk has a fixed-point witness.

## (B) THE BÉZOUT-ROTATION ORBIT LAW — second proof of T2

Primary: **THM-926** (mac-mini, Fourier resonance identity, refereed
10/10). This section is an independent route with new structure, verified
exactly (108 triples; the identity asserted EXACT on every one).

**Structure law (Layer 1, verified).** For coprime A < B, the components
of W = D_A ∩ D_B are indexed by the beat offset m = iB − jA,
|m| < (A+B)/13, exactly ONE component per admissible m; widths follow the
sawtooth min(2A, 2B, (A+B) − 13|m|)⁺/(13AB) (codex's T-law; μ(W) matches
exactly on all test pairs; κ_W = #admissible m). With BB* − AB̃ = 1
(Bézout), the LEFT endpoints march in m as an exact AP with step B*/A on
one sawtooth branch and B̃/B on the other — and since the width slope
∓1/(AB) converts one step into the other (B*/A − 1/(AB) = B̃/B), the RIGHT
endpoints march with the OTHER cofactor. ≤ ~4 regimes (branch change at
m ≈ 0, cap at |m| ≈ (B−A)/13) + junction singletons.

**Engine (Layer 2, verified).** Let ψ(t) = μ([0,t] ∩ D_C) − (2/13)t
(continuous, period 1/C, per-period variation V₁ = 44/(169C)). Then

> err := μ(W ∩ D_C) − (2/13)μ(W) = Σ_components [ψ(right) − ψ(left)]

**exactly**, so err is a difference of two orbit sums of mean-zero ψ̃ along
the rational rotations β = {C·B*/A} and {C·B̃/B}. Denjoy–Koksma over any
continued-fraction denominator q of β gives |Σ| ≤ (⌊N/q⌋ + 1)·V₁ per
family/edge; assembled with junctions it covers the true err on 108/108
triples (worst bound-tightness ratio 0.19; DK inequality holds for every
family of every triple).

**The resonance dictionary.** Exact triple resonance (qC = rA + sB at
small height; e.g. C = A + B) forces q* = 1 and errors up to 1.4×10⁻²;
clean triples sit ≤ 3.6×10⁻⁴. **THM-927's 3-AP blocker IS the q = 2,
r = s = 1 exact triple resonance** — the third blocker and T2's resonance
parameter are one object. Coercivity regime of this CF-form bound:
leading term ≈ 0.16/q* at equal scale, additive floor ≈ 10·V₁ → 0 with C;
below the crossover everything is finitely checkable — the two-scale shape
appearing inside the lemma itself.

## (C) THE FIRST ACCESSIBLE-SCALE FULL-PROBLEM LEVEL-5 CERTIFICATE

The corrected pair condition (post-MISTAKE-151), in THM-864's own shape:

> **weighted-capped cleanliness:** min_{1≤q,p≤13} |q·b − p·a|·(q+p−1) ≥ Θ.

Both coordinates capped (no beat-Dirichlet vacuity); q = p = 1 demands
|b − a| ≥ Θ (kills gcd-quantization); the weight is THM-864's error
functional. Plus Sidon, no-3AP (= no q ≤ 2 triple resonance), and
line-cleanliness at height 3 (no exact qC = rA + sB with q ≤ 3,
|r|,|s| ≤ 3 — justified by THM-926's 1/|k₁k₂k₃| line decay).

**DFS at Θ = 100 on [300, 40000] finds, lexicographically first:**

> **X = {300, 406, 511, 652, 862, 963, 1074, 1357, 1459, 1571, 1776,
> 1991, 2208}**  (max ratio 7.36 — genuinely one band, same-scale)
>
> S₂ = 1.846144 (equid 1.846154), S₃ = 1.039406 (equid 1.041420),
> S₄ = 0.432927 (equid 0.400546), S₅ = 0.200534 (equid 0.110920)
>
> **BONF5 = 1 − S₁ + S₂ − S₃ + S₄ − S₅ ≥ +0.039131 > 0 (exact rational
> arithmetic), actual uncovered = 0.117744.**

The filter did not optimize for BONF5 — the certificate is out-of-sample:
honest admissibility ⟹ coercive, first instance. The S329 emptiness
question, behind the corrected filter, resolves **POSITIVE**: accessible-
scale full-problem level-5 certificates EXIST. (S₂ lands within 10⁻⁵ of
equidistribution — the weighted filter tames the pair terms almost
perfectly; the residual inflation is S₅'s, and the alternating sum
absorbs it.)

## The two-scale certificate, assembled

- **Exponential scale:** Theorem A — every 15-lacunary packet is lonely;
  uniform R ≥ 19 for all n; dilate families have fixed-point witnesses
  (the parallel-class escape is tame).
- **Accessible scale:** the level-5 wall is complete (W, T1 — opus S327/8;
  T2 — THM-926 primary, (B) an independent second engine); c₀ = 17/84 > 0
  (THM-925/926); and the certificate machinery is now DEMONSTRABLY
  COERCIVE at accessible scale — (C).
- **The bridge is resonance:** exact triple resonance = mac-mini's
  q-periodic finite branch (THM-924) = the third blocker (THM-927) =
  q* = 1 in the DK engine (B) = the excluded lines in the filter (C).
  One object, five appearances.
