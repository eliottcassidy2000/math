---
id: HYP-3787
title: BOUNDING THE SIGNED RESONANCE CORRECTION on the lonely set L_C via the RIESZ/THETA form. The S65 coverage frac_w = meas(L_C ∩ {||wt||<r'})/meas(L_C) has an EXACT Fourier identity for its deviation from the equidistribution mean 2r': correction_w := frac_w - 2r' = (2/L) Σ_{j>=1} hat1(jw) sin(2π j r')/(π j), where hat1(k)=\hat{1_{L_C}}(k) is REAL & EVEN (t->1-t symmetry), hat1(0)=L=meas(L_C). This is the OFF-DIAGONAL (k=jw) extension of THM-515's ON-diagonal (k=0) theta form L=Σ_{t∈Λ}∏h. THREE results: (1) EXACT identity verified (FFT reconstruction = direct to <2e-4, all w); (2) TWO-ATOM SIGN LAW hat1(k) ≈ L cos(2π k t*) E(k), t*=n/Phi6 the binding hexagonal point, E(k)<=1 a decaying interval-width envelope -- so the SIGN of the correction is cos(2π k t*): RESONANT(+) at w near harmonics of 1/t*=Phi6/n=13.07 (peaks k=13,26,39,52,65, hat1(13)~+L), ANTI(-) at half-harmonics; (3) RIGOROUS FAR-ELEMENT BOUND: 1_{L_C}=finite union of I(r) intervals => |hat1(k)|<=I/(πk) (TV/jump bound, 0 pointwise violations) => |correction_w| <= (2 I(r)/(π L w)) S(r') = O(1/w), S(r')=Σ_j(1/j)min(2r',1/(πj))<∞ (=185/w at n=14,r=0.07). This makes S65 far-element impotence QUANTITATIVE: a far speed covers L_C only to O(1/w) beyond equidistribution => cannot patch the lonely set. Riesz(Bedert-2025)/theta = the right tool for THIS tail (THM-515)
status: MIXED (exact identity + geometric law + rigorous rate, on a fixed grid). VERIFIED (FFT N=6e5, n=14, r=r'=0.07, L=0.0324): (1) correction_w = (2/L)Σ hat1(jw)sin(2πjr')/(πj) matches direct frac_w-2r' to <2e-4 for w=13,26,39,6,7,12,61,182,183,91,5000,5001. (2) hat1/L vs cos(2πkt*): 0.968/0.999, 0.908/0.998, 0.819/0.995 (k=13,26,39) -- two-atom + decaying envelope. (3) TV bound |hat1(k)|<=I/(πk) (I=28): 0 violations over k=1..400; |correction_w|<=185/w holds for all tested w (200..10000), incl residual-resonant large w. RIGOROUS given I(r) finite. NOT a full theorem: I(r)->∞ as r->M_C, general n not covered, constant loose. Complements HYP-3786 (mechanism) with a rate.
source: klein-2026-07-01-S66
depends_on:
  - HYP-3786   # S65: equidistribution on L_C, far-element impotent (the MECHANISM this quantifies)
  - THM-515    # L_C measure = singular series = theta-form; Riesz-product is the right tool for inf L>0
related:
  - HYP-2873   # additive energy = spectral 4th moment (the L^2 of the correction; AP-core = max energy)
  - THM-527    # Riesz "wrong tool" for anti-lacunary AP -- but RIGHT tool for the far element (dissociated harmonics)
  - THM-501    # the singular series / lonely measure (= hat1(0) = L)
  - HYP-3715   # t* = n/Phi6 = zeta_6 hexagonal point (the two-atom center)
  - HYP-3763   # large-multiples-forced (complementary M-raising mechanism)
  - HYP-3745   # dodge != patch (this bounds the "patch" side quantitatively)
  - OPEN-Q-108 # the r>=2 multi-far uniform bound (this closes r=1 with a rate; r>=2 open)
  - HYP-3132   # the r-far uniform crux (kind-pasteur-S3: same object, multi-far)
results:
  - 04-computation/signed_resonance_correction_klein.py
  - 05-knowledge/results/signed_resonance_correction_klein.out
---

# HYP-3787 — bounding the signed resonance correction on L_C (Riesz/theta form)

## The object (from S65/HYP-3786)
For the covering-min construction `C = {1..n-2, n(n-1)}` (n=14: `{1..12,182}`), the lonely set at level `r`
is `L_C = {t : min_{v∈C}||vt|| > r}`, a Cantor set concentrated at the binding `t* = n/Phi6` (hexagonal
point). A beater must **cover** `L_C`; a test speed `w`'s coverage is `frac_w = meas(L_C ∩ {||wt||<r'})/L`,
`L = meas(L_C)`. S65 found `frac_w = 2r'` on average (Weyl) with high-variance **resonances**. The
**signed resonance correction** is the deviation

> `correction_w := frac_w − 2r'`.

## Result 1 — the EXACT Fourier identity
`1_{L_C}` is real and even (`||v(1−t)||=||vt||`), so `hat1(k) := \hat{1_{L_C}}(k)` is real & even,
`hat1(0)=L`. The danger band `D_w=1_{||wt||<r'}` has `\hat{D_w}(jw)=sin(2πjr')/(πj)` (`j≠0`), `=2r'`
(`j=0`), zero off multiples of `w`. Parseval gives

> **`correction_w = (2/L) Σ_{j≥1} hat1(jw) · sin(2πjr')/(πj)`**  — a SIGNED harmonic sum.

VERIFIED: FFT reconstruction equals the direct grid value to `<2×10⁻⁴` for every tested `w` (resonant,
anti, far). The correction is *literally* the Fourier coefficients of the lonely-set indicator, sampled at
the harmonics `jw` of the test speed, weighted by the danger-band sinc kernel.

## Result 2 — the TWO-ATOM SIGN LAW (why the sign is what it is)
`L_C` concentrates at `t*` and `1−t*` (the binding pair). Two symmetric atoms give

> **`hat1(k) ≈ L · cos(2π k t*) · E(k)`**, `E(k)≤1` a decaying envelope (the FT of the finite interval
> widths of `L_C`).

VERIFIED: `hat1(k)/L` vs `cos(2πk t*)` = `0.968/0.999, 0.908/0.998, 0.819/0.995, 0.697/0.991, 0.589/0.985`
at `k=13,26,39,52,65` (harmonics of `1/t*=Phi6/n=183/14=13.07`); the peaks of `|hat1|` are exactly these.
So the **SIGN of the correction is `cos(2π k t*)`**: `w` is **RESONANT** (`correction>0`, covers `L_C`)
when `w` is near a harmonic of `1/t* = Phi6/n`, and **ANTI-resonant** (`correction<0`, avoids `L_C`) at the
half-harmonics. The construction core `{6..12}` and apex `61=Phi6/3` sit at `cos<0` (anti). The signed
correction is the **phase of `k` against the binding frequency `Phi6/n`** — the E2/Eisenstein hexagonal
frequency (S56/S57) made spectral.

## Result 3 — the RIGOROUS far-element bound (the Riesz/TV certificate)
`1_{L_C}` is a **finite union of `I(r)` intervals** (Cantor set at finite level; `I(0.07)=28`). The FT of
an indicator of `I` intervals decays **pointwise** like `1/k` (each endpoint is a jump):

> **`|hat1(k)| ≤ I(r)/(π k)`**  (0 pointwise violations over `k=1..400`).

Feeding this into Result 1 with `|sin(2πjr')/(πj)| ≤ min(2r', 1/(πj))`:

> **`|correction_w| ≤ (2 I(r)/(π L w)) · S(r')`, `S(r') = Σ_{j≥1} (1/j) min(2r', 1/(πj)) < ∞`  = `O(1/w)`.**

At `n=14, r=r'=0.07`: `S=0.336`, so `|correction_w| ≤ 185/w`. VERIFIED to hold for all tested `w`
(`200..10000`), including large `w` that still hit a residual resonance (e.g. `w=10000 ≈ 765·1/t*`, which
resonates but is still capped by `185/w`). This upgrades S65's *qualitative* "far element impotent" to a
**quantitative rate**: a far speed's coverage of `L_C` exceeds the equidistribution mean `2r'` by only
`O(1/w)`, so it **cannot patch** the lonely set. The `1/w` is the Riesz/dissociation rate: the far
element's harmonics `jw` are dissociated from the short relation lattice of `C`, so they miss `hat1`'s
support (which lives on the short lattice / low frequencies).

## Convergence with kind-pasteur-2026-07-01-S3 (same day, independent)
kind-pasteur-S3 reframed the unbounded LRC(14) case as `survival(r) = (6/7)^r·meas(L_C) − [signed
resonance correction]` (= the singular series localized to a fixed core's `L_C`), and asked precisely for
**"an analytic bound on a signed sum over a fixed set (Riesz product, decorrelation, Erdős–Turán),"**
conjecturing the far-element decorrelation rate `Δ_W = O(1/W)`. This session supplies exactly that
machinery, and closes the single-far case:
- Their **"signed sum over a fixed set"** = my **exact Fourier identity** (Result 1), `correction_w`
  written literally as `Σ_j hat1(jw)·(danger kernel)` over the fixed set `L_C`.
- Their **"ι-odd resonance obstruction (Borsuk–Ulam), vanishes for one comb"** = my **two-atom sign law**
  (Result 2): the sign is `cos(2π k t*)`, an antipodal (`t↔1−t`, ι) phase against `Phi6/n`.
- Their conjectured rate **`Δ_W = O(1/W)`** = my **Result 3**, now a RIGOROUS TV/Riesz certificate
  `|correction_w| ≤ (2I/(πL))·S/w = O(1/w)` for a single far speed — **single-far closed with a rate.**
- **Still open (OPEN-Q-108):** the `r≥2` multi-far uniform bound (signed correction `<` main term over ALL
  bounded cores and far combs). My Fourier machinery extends to the `r`-fold product `∏_i D_{w_i}` (its
  coefficients live on `Σ_i j_i w_i`), but the uniform `r≥2` estimate is not done here.
NOTE (consistency): kind-pasteur's "single-far has no resonance" is about a single **far/large** comb `W`
— exactly my `O(1/w)` regime. My strong single-`w` resonances (`w=13`) are **small core harmonics**
(`≈ n−1 = 1/t*`), part of the fixed structure that DEFINES `L_C`, forbidden as beater additions
(HYP-3778) — not "far combs." No contradiction.

## Why Riesz products are the RIGHT tool here (repo synthesis)
- **THM-515**: `L = meas(L_C) = Σ_{t∈Λ}∏h(t_i)` (theta form, relation lattice `Λ`); the Riesz-product
  method (Bedert 2025, arXiv:2511.16636, *Riesz products and the Lonely Runner Conjecture*) is the state
  of the art for `inf L>0`. My `correction_w` is the **off-diagonal (`k=jw`) extension** of THM-515's
  **on-diagonal (`k=0`)** theta sum. Same object, sampled off the relation lattice.
- **THM-527 / "Riesz wrong tool"**: Riesz products need a **dissociated / lacunary** frequency set; the
  AP-core `{1..n−2}` is **anti-lacunary** (max additive energy, HYP-2873), so Riesz is the wrong tool for
  the *core*. But the **far element `w`** IS dissociated from the core (its harmonics `jw` are lacunary
  against the short lattice) — so Riesz is exactly the **right** tool for the far-element **tail**, which
  is the part S64 left open. The correction splits: near-core resonances (arithmetic, two-atom, NOT
  Riesz-small) vs far-element tail (dissociated, Riesz-controlled, `O(1/w)`).
- **HYP-2873**: `Σ_w correction_w²` is governed by the additive energy `= Σ_k |hat1(k)|²` (spectral 4th
  moment) — the AP-core's maximal energy is exactly why the near resonances are strong and concentrated.

## Net
The signed resonance correction on `L_C` is an EXACT signed harmonic sum of the lonely-set Fourier
coefficients (Result 1); its SIGN is `cos(2π k t*)`, the phase against the hexagonal binding frequency
`Phi6/n` (Result 2, two-atom); and it is RIGOROUSLY `O(1/w)` for the far element (Result 3, TV/Riesz
certificate), quantifying S65's impotence. Riesz products are the right tool for the far-element tail
(dissociated), not the anti-lacunary core (arithmetic). A grid-verified quantitative certificate, not yet
a uniform-in-`n` theorem (the interval count `I(r)` and its `r→M_C` blow-up are the next obligation).
