# Spectral flatness vs additive energy as an LRC monotone: the winner is set by the CAP's SCALE — additive energy (= spectral 4th moment) wins for the COARSE mod-7 coverage proxy measS7, but the Wiener entropy of the FINE danger cover wins for the TRUE LRC gap M — decisively disambiguated on ONE family (winner flips with the cap, not the family); yet NEITHER scalar certifies consec-max (HYP-2738 impossibility: the AP does not even minimize the entropy, so averaged/entropy lenses stay blind to the pointwise tight locus, and certification needs a SIGNED certificate)

*opus-2026-06-30. Owner's hint: "the right order parameter is spectral flatness / Wiener entropy of the
danger cover, not additive energy; additive energy is a scalar shadow that throws away multi-scale phase
info; concrete falsifiable test: does spectral flatness have fewer inversions than additive energy against
p0?" I ran the exact test both ways. The hint is RIGHT for the right target and the split is clean.*

## The test, run faithfully both ways
The repo's additive-energy monotone analysis (`lrc_q108_threadB_addenergy_monotone_kpswf4.py`, HYP-2738)
uses cap `p0 = measS7(E)` = measure of `t` where all 7 sectors mod 7 are hit — a COARSE (mod-7) coverage
proxy — over the family of 8-subsets of `{1..15}` grouped by span; additive energy `A(E)=#{a+b=c+d}` had
1368 within-span inversions at span 12. I added the owner's order parameter (Wiener entropy of the danger
cover, `W = -<log m(t)>`, and spectral flatness `GM/AM` of its spectrum) and counted inversions against
BOTH caps: the coarse `measS7` and the TRUE fine LRC gap `M(S)=max_t min_v ||vt||`.

| family | cap | additive energy A | Wiener entropy W | winner |
|---|---|---|---|---|
| near-AP single-swaps (n=8,9,10) | fine `M` | inv-frac 0.39–0.45 | **inv-frac 0.19–0.23** | **W** |
| span-12 (8-subsets of 1..15) | coarse `measS7` | **235 181** | 411 288 | **A** |
| span-12 (SAME family) | fine `M` | 369 437 | **233 735** | **W** |

Pearson on the span-12 family: `|A vs measS7| = 0.517`, `|W vs M| = 0.606`. On near-AP, `|W vs M| = 0.74`
vs `|A vs M| = 0.09`.

## The decisive disambiguation: the winner flips with the CAP's SCALE, not the family
On the SINGLE span-12 family, holding the sets fixed and changing only the cap, the winner flips:
**coarse `measS7` -> additive energy; fine `M` -> Wiener entropy.** So the effect is NOT a family artifact
— it is a **scale-matching law**:
> The order parameter must match the target's scale. Additive energy IS the spectral 4th moment
> (HYP-+2873, `A = int|T_S|^4`), a COARSE concentration statistic — it tracks the coarse mod-7 coverage
> `measS7` (both maximized by the interval, Fejer-optimality). The **Wiener entropy of the fine danger
> cover** is a FINE (log-weighted, all-multiplicity) statistic — it tracks the fine `L^inf` gap `M`, which
> lives at the low-multiplicity points additive energy averages over.

**The owner's instinct is vindicated for the right target.** "Additive energy throws away the multi-scale
phase information" is exactly right *for the true LRC gap*: the entropy retains the fine-scale structure and
beats additive energy as an `M`-monotone on two independent families. The repo's additive-energy framing was
measuring the COARSE proxy `measS7` (where `A` legitimately wins); against the gap we actually care about,
the winner is the entropy.

## But NEITHER scalar certifies consec-max (the honest ceiling — HYP-2738)
The win does NOT make `W` "the missing monotone" that proves consec-max, for a proven reason:
- **The AP does not minimize `W`.** Over the near-AP family, `min` swap-`W` `< ` AP-`W` at n=8,9,10
  (`wiener.py`): a better fine-`M`-*correlate* is still an AVERAGE, so it stays blind to the pointwise
  tight LOCUS (`M=1/n` exactly). This is the same meta-insight as
  `the-tight-skeleton-is-COMBINATORIAL-COVERING`: tightness is `min_t m(t) >= 1` (`L^inf`), and every
  averaged/entropy functional washes out the single widest hole that decides it.
- **HYP-2738's impossibility theorem is formal.** Any nonneg test function `phi` with `phi(0)=1` gives
  `p0 <= E[phi(N)]`; a `phi` that is itself consec-MAX yields the LOOSEST bound at consec, so it cannot
  certify consec is the argmax. Certification forces a SIGNED (Bonferroni / inclusion-exclusion)
  combination — the repo's `L_y = 1 - E[N] + sum_a w_a C_a` with negative weights. No single nonneg scalar
  (A, `W`, or flatness) can do it.

So the correct reading of the owner's hint: it identifies a **strictly better PROXY** for the fine gap
(entropy over energy — real, quantified above), but the barrier to consec-max is not "find a better scalar"
— it is that consec-max is a saturated LP (Delsarte), which no scalar order parameter can certify. The
missing object is a **signed certificate**, and the true open crux is the CRT-linkage residual (HYP-3749),
not a monotone.

## Why this still MOVES the program (the redirection it sharpens)
The scale-matching law reinforces the standing redirection "attack the lowness POINTWISE/FINE, not
coarse-averaged":
- It explains WHY the coarse lenses (additive energy, `int m^2`, mod-7 `measS7`) plateaued — they are the
  wrong SCALE for the fine `L^inf` gap, not merely the wrong functional.
- It aligns with **klein's live route** (S51 / HYP-3762): "tight => support >= n-5" via the wide-hole
  inequality (HYP-3749) is a FINE / pointwise support bound — the right-scale object. The entropy result
  says: if you want a continuous surrogate for `M`, use the log-multiplicity (fine), and it will track the
  gap; but to CERTIFY the pointwise extremum, use the wide-hole/support bound (klein), not a scalar.
- The fine surrogate that DID track `M` is `W = -<log m(t)>` — a candidate soft-max stand-in for the widest
  hole; worth trying as the smoothing in a Beurling-Selberg / Fourier-positivity certificate (the
  `L^inf` route, HYP-2974), where a log-barrier is exactly the kind of signed penalty HYP-2738 demands.

## Status
- **Confirmed (opus, two families):** against the true fine LRC gap `M`, the Wiener entropy of the fine
  danger cover beats additive energy as a monotone (near-AP inv-frac 0.19–0.23 vs 0.39–0.45, Pearson 0.74
  vs 0.09; span-12 233 735 vs 369 437 inversions, Pearson 0.606). The owner's hypothesis HOLDS for the
  right target.
- **Falsified (opus, exact repo setup):** against the coarse mod-7 proxy `measS7`, additive energy beats
  the entropy (235 181 vs 411 288) — `A` is the right scale for `measS7`, the wrong scale for `M`.
- **Decisive (opus):** the winner flips with the CAP on ONE fixed family — a scale-matching law
  (coarse->A, fine->W), not a family artifact.
- **Ceiling (HYP-2738, confirmed here):** no scalar monotone certifies consec-max — the AP does not even
  minimize `W`; averaged/entropy lenses are blind to the pointwise tight locus; certification needs a
  signed certificate. The real crux stays the CRT residual (HYP-3749).

Related: `the-tight-skeleton-is-COMBINATORIAL-COVERING` (averaged lenses blind, pointwise), HYP-2738
(consec-max irreducibly aggregate / impossibility), HYP-+2873 (additive energy = spectral 4th moment),
HYP-2974 (Fourier-Toeplitz PSD, the `L^inf`/signed route), HYP-3749 (wide-hole / CRT residual), HYP-3762 +
klein-S51 (tight => support >= n-5, the fine/right-scale route), HYP-3763 (this); OPEN-Q-108.
Scripts: `04-computation/lrc_flatness_vs_energy_scale_opus_20260630.py` (nearap+wiener+exacttest+disambig).
