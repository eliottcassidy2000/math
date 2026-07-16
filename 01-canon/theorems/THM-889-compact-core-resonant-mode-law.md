---
id: THM-889
title: THE COMPACT-CORE RESONANT-MODE LAW — "compact-core flatness" resolved. (I) The independent-limit miss-measures are exact rationals A* = 360/7⁵, B* = 120/7⁵ (A* = 3B*), giving m̂*_s(a) = −(A*+B*)e(as/7) − A* ≠ 0 with max_a |1−ω^a||m̂*(a)| = 0.0495. (II) Incoherent compact cores obey the owner-resonant law |S(peak near ea)| ≈ e·|1−ω^a|·|m̂*(a)| — validated to ratios 1.03/1.09 on the largest clean owners — so absolute-constant flatness is FALSE asymptotically with EXPLICIT small slope ε ≈ 0.05–0.09·Vmax. (III) Internally-coherent compact cores (consecutive = near-dilates of c·(1,…,1)) instead follow the dilate/difference-flow mode: |S| ≈ 0.98·Vmax at n ≈ 6·mean — the steep species. (IV) PROGRAM-RANGE FLATNESS (the usable statement): for Vmax ≤ 500, max|S|²/M ≲ 12 (incoherent) resp. bounded by the coherent species' explicit ~0.27c only for near-dilate cores, which are the sheet tiles' territory; every needed instance is decidable (THM-881 w-freeness)
status: (I) PROVED-EXACT (full Z₇⁵ enumeration, 360/120 counts); (II)(III) LAW VALIDATED (probes to c = 2000 consecutive / c = 300+ generic; clean-owner ratios 1.03, 1.09; slope in the predicted band); the asymptotic law statement is census-grade, its equidistribution proof (m̂^{(e)} → m̂* with rate for incoherent cores) is the named residual; (IV) explicit-constant statement on the program range, per-instance decidable
source: death-star-2026-07-16-S19 (owner directive: "prove compact-core flatness"); extends klein's THM-883-resonant-mode (cont.4, scale-separation route) into the balanced regime where separation is unavailable
depends_on:
  - THM-883 (klein resonant-mode, cont.4)   # the far-bank law this extends
  - THM-881/880/729                          # descent frame, bilinear form, application
  - HYP-7017 (S18)                           # the refutation + census this completes
related: [HYP-3901 (difference flow — the coherent species' slow system), THM-757/759 (sheet tiles owning near-dilates), HYP-7016 (compact-core seam)]
verification: 04-computation/compact_core_resonant_mode_law_deathstar_S19.py -> 05-knowledge/results/compact_core_resonant_mode_law_deathstar_S19.out
---

# THM-889 — the compact-core resonant-mode law

**The directive was "prove compact-core flatness." The honest resolution: flatness with an
absolute constant is false even for balanced clusters — what is true, and provable, is the
resonant-mode LAW with exact constants, whose slope is small enough that flatness holds with
explicit constants on the entire range the LRC(14) assembly needs.**

## (I) The independent-limit miss-measures (exact)

For a stationary runner pinning section 0 and five independent uniform sections (the
incoherent-core limit of the "others" system): with s ∉ {0}, c ∉ {0, s},

- `A* := P[occupied = Z₇∖{s,c}] = 360/16807` (inclusion–exclusion: 5⁵−4·4⁵+6·3⁵−4·2⁵+1),
- `B* := P[occupied = Z₇∖{s}] = 120/16807` (6⁵−5·5⁵+10·4⁵−10·3⁵+5·2⁵−1), `A* = 3B*`,
- `A_0 = 0` (the pinned section cannot be missed — the stationary runner breaks c-symmetry).

Hence the miss-measure Fourier coefficient (klein's m̂ in the balanced limit):
> **m̂*_s(a) = Σ_{c≠s}(A_c + B*)e(ac/7) = −(A*+B*)·e(as/7) − A*,  |m̂*| ≥ B* > 0,**
> **max_a |1−e(a/7)|·|m̂*_s(a)| = 0.0495.**

Verified by full Z₇⁵ enumeration (counts exactly 360 per admissible c, 120; script P1).

## (II) The incoherent species: the law, validated

For compact cores with generic internal ratios, each owner e carries resonant modes near
the frequencies n = ea whose height follows the law `e·|1−e(a/7)|·|m̂*(a)|`. Probes
(c = 300 core [0,300,412,573,850,1239,1762], s = 3): owner 1239, a = 2: measured 89.82 vs
predicted 87.45 (**ratio 1.03**); owner 1762, a = 5: 135.25 vs 124.36 (**1.09**); at c = 1000
(Vmax = 5871): owner 2831, a = 2: 211.60 vs 199.81 (**1.06**), owner 4130, a = 2: 306.01 vs
291.49 (**1.05**) — with the peak locations landing at n = 2e EXACTLY. Owners
with residual internal near-relations deviate upward (flagged in the .out) — the deviation
is itself a coherence meter. Consequently max|S| ≈ ε·Vmax with ε ≈ 0.05–0.08 measured (0.0777 at Vmax = 5871)
(predicted band [0.014, 0.099]), i.e. **max|S|²/M grows ~ ε²·Vmax²/M — absolute flatness is
false, but with a slope ~200× shallower than the trivial bound**.

## (III) The coherent species: consecutive cores are near-dilates

`[0, c..c+5] = c·(0,1,…,1) + (0,0,1,2,3,4,5)`-type: ONE coherent block with slow system =
the differences (the HYP-3901 difference flow appears exactly here). The spectrum follows
the dilate mode: |S| ≈ 0.98·c at n ≈ 6·(mean speed) (measured 0.982c at c = 1000, 0.981c
at c = 2000; ratio |S|²/M ≈ 0.27c — the steep species). Exact scaling face: **dilate
covariance** S_{cE}(cn) = c·S_E(n), M_{cE} = c·M_E. My S18 "compact cores are flat"
census (c ≤ 30) was pre-asymptotic for exactly the same reason klein's C = 14 was at
t ≤ 50 — logged as the instructive pattern (three pre-asymptotic traps in 24h).

## (IV) Program-range flatness (the usable theorem) + division of labor

- Incoherent compact cores, Vmax ≤ 500: resonant modes ≤ 0.05·500 = 25 ⟹
  max|S|²/M ≤ flat floor + 625/M ≲ 12 (measured interpolation: 8.05 → 12.0 → 20.6 at
  Vmax = 177/588/1762). Every specific instance the assembly needs is decidable by ONE
  finite scan (THM-881 w-freeness) — and the resonant frequencies to check are now
  EXPLICIT (n ≈ ea), so certified per-instance constants come at probe cost, not
  full-scan cost.
- Coherent/near-dilate compact cores: the steep species — but these are precisely the
  structured families the sheet/coherent-ray tiles (THM-757/759) close; their resonant
  w-classes (w ≈ ℓ⁻¹·6·mean etc.) feed klein's peel-audit handoff.
- The named analytic residual (honest): prove m̂^{(e)} → m̂* with an explicit rate for
  incoherent cores (a finite-orbit equidistribution statement on the rational 6-torus;
  Erdős–Turán over the exact lattice, coherence parameter = min small-combination
  |Σkᵢeᵢ|), upgrading (II) from census-grade to proved. The law's algebraic content —
  (I) and the peak locations — is exact already.

## One-line summary

**Compact-core flatness holds where the program lives (explicit constants, decidable
instances), fails asymptotically by an exact small-slope law (m̂* = exact rationals), and
the steep exceptions are exactly the coherent families other tiles already own.**
