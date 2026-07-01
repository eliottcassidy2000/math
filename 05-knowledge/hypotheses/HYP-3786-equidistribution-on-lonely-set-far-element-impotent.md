---
id: HYP-3786
title: EQUIDISTRIBUTION ON THE FIXED LONELY SET L_C -- the far element is IMPOTENT on L_C (narrow arcs => equidistributes => can't patch the lonely set), extending/closing the S64 far-element tail. For the covering-min construction C, the lonely set L_C(r)={t: min_{v in C}||vt||>r} is a CANTOR set (84,48,28,24,...,2 intervals as r: 0.05->M_C) CONCENTRATED at the binding t*=n/Phi6 (lonely measure 0.024 at the floor 1/n). A beater must COVER L_C to reduce M. KEY: (Weyl) the average coverage frac_w = meas(L_C∩{||wt||<r'})/meas(L_C) over w EQUALS 2r' EXACTLY (equidistribution on L_C on average), with HIGH variance from RESONANCE. Resonance law: w covers L_C iff (i) ||w t*|| < r' (w near a harmonic of 1/t*=Phi6/n=183/14) AND (ii) w is SMALL (arc width 2r'/w wide enough to cover L_C's intervals). A HUGE w has arc width 2r'/w -> 0, so it covers only ~2r'*meas(L_C) of L_C REGARDLESS of tuning: verified, harmonic-tuned huge w (~k*Phi6/n, w~5000) gives frac~0.149~2r', NOT concentrated. So the FAR ELEMENT CANNOT PATCH L_C -- it equidistributes -- closing the S64 residual (speeds>n(n-1)) by equidistribution, complementary to HYP-3763 (large multiples raise M). Only SMALL speeds (the core {1..n-2}, which DEFINE L_C by being far there, + the missing n-1) can cover L_C, and the covering structure (ILP, HYP-3778/3779) forbids a beater. Reconciles HYP-3745: a huge speed can DODGE small-modulus witnesses (be a valid covering member) but cannot actively COVER L_C (beat)
status: MIXED (empirical grid + clean mechanism). VERIFIED (grid N=4e5, n=14): L_C(r) Cantor structure + measures; mean frac_w = 0.1530 = 2r' exactly (Weyl on L_C); resonance law frac decreasing in ||w t*|| (w=13:frac 1.0, ||wt*||=0.0055; w=182: frac 0.016; w=7,61: frac 0); HUGE w (random AND harmonic-tuned) -> frac~0.149-0.153 (equidistributes, arc width 2r'/w->0). MECHANISM (not a theorem): far element impotent on L_C by narrow-arc equidistribution; closes the S64 tail heuristically. Complementary to HYP-3763 (M-raising) and HYP-3745 (dodge!=patch). Grid-approximate.
source: klein-2026-07-01-S65
depends_on:
  - HYP-3784   # S64: the Delsarte dual defeated by far-element equidistribution (this explains why)
  - HYP-3763   # large-multiples-forced (the complementary M-raising mechanism)
related:
  - HYP-3745   # CRT-escape: dodge (satisfy witnesses) vs patch (cover L_C) -- reconciled here
  - HYP-3779   # lazy-cut (no beater <= n(n-1)); this addresses > n(n-1) via equidistribution
  - HYP-3762   # three-gap (the arithmetic structure of L_C)
  - THM-501    # the lonely measure (meas L_C)
  - HYP-3715   # t* = n/Phi6 = zeta_6 hexagonal point (L_C concentrates there)
results:
  - 04-computation/lonely_set_equidistribution_klein.py
  - 05-knowledge/results/lonely_set_equidistribution_klein.out
  - 05-knowledge/results/lonely_set_harmonic_klein.out
---

# HYP-3786 — equidistribution on the fixed lonely set; the far element is impotent

## The fixed lonely set L_C
For the covering-min construction `C = {1..n-2, n(n-1)}`, the lonely set at level `r` is
`L_C(r) = { t : min_{v in C} ||v t|| > r }`. It is a **Cantor set** concentrated at the binding
`t* = n/Phi6` (the `zeta_6` hexagonal point, HYP-3715): as `r` rises from `0.05` to `M_C = 14/183`, the
measure shrinks `0.232 -> 0.116 -> 0.032 -> 0.024 -> 0.00001` and the interval count `84 -> 48 -> 28 ->
24 -> 2` (collapsing to the binding pair `{t*, 1-t*}`). At the LRC floor `1/n`, `meas(L_C) = 0.024` > 0 --
the construction's lonely measure (its looseness; `M_C > 1/n`), the object of THM-501.

## Equidistribution on L_C: Weyl on average, resonance in the variance
A beater must **cover** `L_C` (kill the loneliness) to reduce `M` below `M_C`. Measure the coverage
`frac_w = meas(L_C ∩ {t: ||w t|| < r'}) / meas(L_C)` (`r' = M_C`). Then:
- **(Weyl) the AVERAGE over `w` is `2r'` EXACTLY** (`mean frac_w = 0.1530 = 2*0.0765`): a generic speed's
  danger equidistributes on `L_C`, covering the same proportion `2r'` it covers of the whole circle.
- **HIGH variance (std 0.11) = RESONANCE**: some `w` cover ~all of `L_C`, some ~none.

## The resonance law: harmonics of `1/t* = Phi6/n`, AND small
`w` covers `L_C` (concentrated at `t*`) iff **both**:
- (i) `||w t*|| < r'` -- `w` is near a **harmonic** of `1/t* = Phi6/n = 183/14 = 13.07`
  (`frac` decreases with `||w t*||`: `w=13 -> 1.00` (`||w t*||=0.0055`), `w=39 -> 0.87`, `w=65 -> 0.63`,
  `w=182 -> 0.016`); and
- (ii) `w` is **SMALL** -- the danger arc width `2r'/w` must be wide enough to cover `L_C`'s intervals.
The core `{6..12}` and the apex `7`, `61` (`Phi6=183=3*61`) are **anti-resonant** (`frac = 0`): they
DEFINE `L_C` (it is where they are far).

## The far element is IMPOTENT on L_C (the extension / S64-tail closure)
A **huge** speed `w` has danger arc width `2r'/w -> 0`, so even if `||w t*||` is small it covers only a
vanishing sliver of `L_C`'s interval -- it **equidistributes** (`frac -> 2r'`). Verified: harmonic-tuned
huge `w ~ k*Phi6/n` (`k~400`, `w~5000`) gives `frac ~ 0.149 ~ 2r'`, **NOT** concentrated -- the same as
random huge `w`. So:
> **A huge speed cannot patch the lonely set `L_C` -- it equidistributes, regardless of tuning.**
This closes the S64 far-element residual (speeds `> n(n-1)`) *by equidistribution*: a beater cannot use a
huge speed to cover `L_C` (narrow arcs), and cannot use a small speed either (the core defines `L_C`; the
only small resonant speed is the missing `n-1 = 1/t*`-harmonic, whose use changes the covering and raises
`M` elsewhere, HYP-3763/ILP HYP-3778). Complementary mechanisms: **HYP-3763** (large multiples raise `M`)
and **equidistribution** (large speeds can't cover `L_C`) both forbid the far-element beater.

## Reconciling the CRT-escape (HYP-3745)
HYP-3745's CRT-escape: a huge tuned speed can DODGE every small-modulus witness (be a valid covering
member). This is consistent: dodging witnesses (`||w t|| < 1` at small moduli, tunable by CRT) is NOT the
same as **covering `L_C`** (needing wide arcs at the binding). A huge speed can be *safe* but cannot
*patch* the lonely set. **Dodge != patch.** The lonely set `L_C` -- concentrated at `t*` with narrow
intervals -- is exactly what a huge speed's narrow arcs cannot cover.

## Extensions (as far as this goes)
- **Fourier / Salem.** The variance of `frac_w` = the Fourier coefficients of the `L_C`-measure at the
  harmonics of `Phi6/n`. `L_C` is NOT Salem (arithmetic peaks at the harmonics) -- the arithmetic
  structure IS the resonance. A Salem (random) lonely set would have no resonances; the arithmetic one
  does, but only small speeds exploit them.
- **Dynamical.** `t -> {v t}` on `L_C` is the runner flow; equidistribution = unique ergodicity on
  average; the harmonics `k*Phi6/n` are the "eigenfrequencies" of `L_C`.
- **The binding as a `zeta_6` atom.** `L_C` concentrates at `t* = n/Phi6` (the hexagonal/Eisenstein
  point, S57); its harmonics are multiples of `Phi6/n`. The lonely set is an Eisenstein-point atom, and
  the far element (a huge speed) is blind to it.

## Net
The lonely set `L_C` is a Cantor set at the binding `t* = n/Phi6`; speeds equidistribute on it on average
(`frac = 2r'`, Weyl) with resonances at the small harmonics of `Phi6/n`. The **far element is impotent**:
a huge speed's danger arcs are too narrow to cover `L_C`, so it equidistributes and cannot patch the
lonely set -- closing the S64 far-element tail by equidistribution (complementary to HYP-3763). A beater
would need a small resonant speed, which the covering structure forbids (HYP-3778). Dodge (HYP-3745) is
not patch: the huge speed is safe but blind to the lonely set. Empirical/mechanistic, not yet a theorem.
