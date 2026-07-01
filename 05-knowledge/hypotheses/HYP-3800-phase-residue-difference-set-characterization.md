---
id: HYP-3800
title: THE PHASE-RESIDUE / DIFFERENCE-SET CHARACTERIZATION of the far-speed interaction with L_C -- new definitions + general-n structural facts + a proof scaffold for covering-min = n/Phi6. A far speed w interacts with the lonely set L_C through (i) its PHASE-RESIDUE p(w)=n*w mod Phi6 in Z/Phi6 (the DIRECTION of coupling; resonant iff p near 0, anti iff p near Phi6/2) and (ii) a coupling AMPLITUDE that decays with arc width (single-far O(1/w); pairwise-via-difference O(1/δ), SCALE-INVARIANT). The whole multi-far structure is governed by the DIFFERENCE-SET {w_i-w_j} and its phase-residues p(δ)=nδ mod Phi6 (translation-invariant): corr2(a,a+δ) depends ONLY on δ (verified: ~0.10 for a=300..2000 at δ=13), large&POSITIVE(redundant) iff δ∈(n-1)Z small. GENERAL-n structural facts (all proved-elementary): CF(t*)=[0; n-1, n] for ALL n (since Phi6/n=n-1+1/n); the KILLER IDENTITY n(n-1)=Phi6-1≡-1 mod Phi6 makes n-1 the UNIT phase-step (p(w+(n-1))=p(w)-1) = why the resonance lattice is (n-1)Z (the dropped speed = 1st CF convergent denominator); p(w) is a BIJECTION on Z/Phi6 (gcd(n,Phi6)=1); Phi6(n) factors over EISENSTEIN primes (≡1 mod 6, +3 if n≡2 mod3), so the finite phase problem Z/Phi6 CRT-factors over them, mod the antipode Z/2 (prime 2, iota: p<->-p). PROOF SCAFFOLD for covmin=n/Phi6: [bounded speeds<=n(n-1): lazy-cut ILP HYP-3779 DONE] + [single & <=6 far: Fourier+TV mac-mini HYP-3787 DONE] + [multi-far>=7: survival_r=P_{L_C}(all safe)>0, governed by the difference-phase histogram; REDUNDANCY(+)=helps survival, SPREADING(-)=the only threat & weak; close via a correlation inequality (Janson/FKG/kps-S4 signed-CS) bounding the negative part by the off-lattice (equidistributed) differences]
status: SYNTHESIS (new defs + general-n facts proved-elementary + scaffold). VERIFIED (grid N=6e5, n=14): p(w) predicts corr1 sign/direction (11/15; the 4 misses are large-w with |corr1|<0.04 = magnitude decayed, phase!=amplitude); corr2(a,a+δ) scale-invariant in a (δ=13: 0.099/0.103/0.104/0.087 for a=1000/2000/500/300) and depends on p(δ); killer identity, CF, bijection, Eisenstein-CRT all exact. The GENERAL-n facts are proved (elementary mod-Phi6 arithmetic). The PROOF SCAFFOLD is a route, NOT a proof: the multi-far>=7 step (bounding the spreading/negative part uniformly) is OPEN-Q-108, unclosed. Phase governs direction, amplitude is the width -- so the reduction is NOT purely finite (amplitude decay is essential); the finite object is the difference-phase HISTOGRAM.
source: klein-2026-07-01-S68
depends_on:
  - HYP-3791   # S67: 13-lattice self-similarity + resonance=>redundancy (renumbered from 3789; the pairwise structure this generalizes)
  - HYP-3790   # S66: single-far signed correction O(1/w) (the amplitude decay; = mac-mini HYP-3787)
related:
  - HYP-3787   # mac-mini: single & <=6 far RIGOROUS (Fourier+TV) -- scaffold step 3
  - HYP-3788   # mac-mini: |H|>=2 multi-patch equidistribution -- scaffold step 4 baseline
  - HYP-3792   # mac-mini-S77: safe-band residue frame -- STRONG convergence (same residue frame + killer identity k(n-1)->-k + CF[0;n-1,n])
  - THM-515    # additive energy governs L; the difference-set additive energy IS the redundancy magnitude
  - HYP-2873   # additive energy = spectral 4th moment (the difference-set mass on (n-1)Z)
  - THM-503    # 7-vanishing (prime n/2 gating); L singular integral not Euler product
  - HYP-3762   # three-gap / continued fraction (CF=[0;n-1,n] source)
  - HYP-3715   # t*=n/Phi6 hexagonal/Eisenstein point
  - OPEN-Q-108 # the multi-far>=7 uniform bound (scaffold's open step)
  - HYP-3132   # kind-pasteur multi-far crux; kps-S5 moment relaxation (inf meas(L_C)>6^{-r}) = complementary
results:
  - 04-computation/phase_residue_reduction_klein.py
  - 05-knowledge/results/phase_residue_reduction_klein.out
---

# HYP-3800 — the phase-residue / difference-set characterization + a proof scaffold

## New definitions (the key angles to track)
Let `n` be the runner count, `Phi6 = Phi6(n) = n^2-n+1`, `t* = n/Phi6` the binding (hexagonal) point,
`L_C = {t : ‖v t‖ > r ∀ v ∈ core}` the lonely set, `core = {1,…,n-2}`.

- **Coupling phase** of a speed `w`:  `φ(w) = w·t* mod 1 = nw/Phi6 mod 1`.
- **Phase-residue**:  `p(w) = nw mod Phi6 ∈ Z/Phi6`  (integer form; `φ = p/Phi6`). This is the **direction**
  of `w`'s coupling to `L_C`: **resonant** (danger aligns with the binding atom) iff `p(w)` near `0`;
  **anti** iff `p(w)` near `Phi6/2` (the antipode).
- **Resonance period**  `ρ = n-1`: the first continued-fraction convergent denominator of `t*`, and the
  speed the construction DROPS. `p(w+ρ) = p(w) - 1` (see killer identity) — the **unit phase-step**.
- **Coupling amplitude**: `p(w)` gives direction, but the SIZE decays with arc width — **single-far
  `O(1/w)`** (needle), **pairwise-via-difference `O(1/δ)`**, and the pairwise part is **scale-invariant**
  (`corr2(a,a+δ)` depends only on `δ`, not `a`). So the reduction is **not** purely finite: the amplitude
  decay is what makes far speeds impotent.
- **Difference-phase histogram** of a far-speed set `H`: the multiset `{ p(w_i - w_j) : i≠j } ⊂ Z/Phi6`.
  This translation-invariant finite object governs the multi-far correction; its mass near `0` (= additive
  energy of `H` on `ρZ`) is the **redundancy**.
- **Redundancy–spreading dichotomy**: a pairwise correlation is **redundant** (`corr2 > 0`, dangers lock,
  double-cover the same part of `L_C`, HELPS survival) iff its difference `δ ∈ ρZ` is small (`p(δ)≈0`);
  otherwise **spreading** (`corr2 < 0`, covers different parts, the only threat to survival) but WEAK.

## General-n structural facts (elementary, proved)
1. **Continued fraction** `t* = n/Phi6 = [0; n-1, n]` for all `n`  (`Phi6/n = n-1 + 1/n`).
2. **Killer identity** `n(n-1) = Phi6 - 1 ≡ -1 (mod Phi6)`. Hence `p(w+(n-1)) = p(w) + n(n-1) = p(w) - 1`:
   `n-1` is the unit phase-step ⟹ the resonance lattice is `(n-1)Z` (S67's `13Z`). Also the construction's
   killer `n(n-1)` is literally `≡ -1 mod Phi6` — the "`-1`" that anchors the phase lattice.
3. **Bijection**: `gcd(n, Phi6) = gcd(n, n^2-n+1) = gcd(n,1) = 1`, so `p(w)=nw` is a bijection on `Z/Phi6`.
4. **Eisenstein-prime CRT**: `Phi6(n)` is a product of primes `≡ 1 (mod 6)` (plus `3` once, when `n ≡ 2
   mod 3`) — the Eisenstein/hexagonal primes (S56/S57). `Z/Phi6` CRT-factors over them; the antipode
   `ι: p ↔ -p` (the prime-`2`/`Z2` sign grading, from `t ↔ 1-t`) acts on the whole thing. (n=14:
   `Phi6=183=3·61`.)

## The proof scaffold for covering-min `= n/Phi6`
No covering set `S` beats `M_C = n/Phi6`. Split `S = core ∪ far` (small `≤ n(n-1)` vs large).
- **[Step 1 — bounded]** speeds `≤ n(n-1)`: **lazy-cut ILP** (HYP-3779). DONE (exact).
- **[Step 2 — single & ≤6 far]** `M(core ∪ {w}) ≥ r` via `covered(w) = 2r|L_C| + signed correction`,
  `|correction| ≤ N/(3w)` (arc-count TV bound), union bound to `≤6` far speeds: **mac-mini HYP-3787
  (= klein HYP-3790)**. DONE (rigorous, explicit thresholds).
- **[Step 3 — multi-far `≥7`]** survival `= P_{L_C}(all far safe) > 0`? THE OPEN CRUX (OPEN-Q-108). The
  synthesis pins what must be shown:
  `survival_r = (1-2r)^r  +  Σ (signed correlations)`, and the correlations are governed by the
  **difference-phase histogram**. **Redundant (positive)** correlations only INCREASE survival (combs
  double-cover) — favorable. The threat is the **spreading (negative)** part, which S67 found is WEAK and
  lives OFF the resonance lattice (equidistributed differences). So the crux is a single analytic
  inequality: **bound Σ(negative correlations) < (1-2r)^r**, i.e., the spreading power of far combs cannot
  overcome independence. Tools now aligned to this exact statement: **kind-pasteur-S4 signed
  Cauchy–Schwarz** (finite `|corr| ≤ √meas·‖g‖₂`), **kind-pasteur-S5 moment relaxation** (`inf meas(L_C) >
  6^{-r}`), **mac-mini per-speed TV decay**, and **the difference-set additive-energy control** (THM-515/
  HYP-2873: energy on `(n-1)Z` = redundancy = positive = harmless). A correlation inequality
  (Janson/FKG-type, or the signed-CS made uniform) over the difference-phase histogram is the last mile.

## Why this is progress toward a proof (not just a picture)
- It **collapses the scale**: the pairwise (and by extension `r`-fold) interaction is translation-invariant
  — a function of the difference-set, not the unbounded speeds. The "search over huge speeds" becomes a
  finite question about a phase-histogram on `Z/Phi6`.
- It **fixes the sign**: the strong correlations are provably-structurally redundant (locking), which is
  the FAVORABLE direction. A proof no longer needs to fight the resonances; it needs only to bound the
  weak, off-lattice, equidistributed spreading part.
- It **names the finite arithmetic**: everything is governed by `Z/Phi6`, the Eisenstein-prime CRT, and
  the antipode `Z/2` — the same hexagonal/Eisenstein structure the covering-min value `n/Phi6` came from.
  The prime `2` is the sign (antipode); `n/2` gates (7-vanishing); the primes of `Phi6` carry the phases.

## Honest status
General-n facts (1)–(4) are proved (elementary). The characterization (phase = direction, amplitude =
width, pairwise = difference-governed & scale-invariant) is grid-verified. The **scaffold's Step 3 is
NOT closed** — bounding the spreading part uniformly over all far configurations is OPEN-Q-108. This is a
consolidation + reduction, converting the open crux into one analytic inequality on a finite phase-object.
