# The LRC(14) infimum is a quantized gap off the tight locus — the extremizer is a sporadic-tight perturbation, the repo's 0.0052 was a restricted-search artifact, and inf L>0 reduces to "extremizers have bounded lcm"

**Source:** kind-pasteur-2026-06-16-S7. Dispatch: spend a long session on new angles toward an
LRC(14) proof, taking every concept as a path. I worked the lonely-measure form
`L(S)=meas{τ:||vτ||>1/14 ∀v∈S}` (THM-515; `inf L>0 ⟹ C'(14) ⟹ LRC(14)`) with **exact** rational
arc-sweep arithmetic (not the grid the prior sessions used), and two clean levers fell out plus a
genuine correction of the inf.

## Two exact levers (THM-522)

1. **Scale-invariance `L(cS)=L(S)`** (proof: `τ↦cτ` is measure-preserving on the circle).
   `L` lives on the **primitive projective shape**; the "stranger scale" and overall scale are
   one degree of freedom, not two.
2. **Quantization `L(S) ∈ (1/(14·lcm S))·ℤ`** (every danger-arc endpoint `(14k±1)/(14v)` is a
   multiple of `1/(14·lcm)`), so **`L>0 ⟹ L ≥ 1/(14·lcm S)`**. The lonely measure is a *rational
   with controlled denominator* — never an irrational infinitesimal.

These give a **compactness reduction** of the prize: `inf L>0 ⟺ the L-minimizing configs have
bounded lcm`. The bounded-lcm part is free (quantization); scale-invariance kills the dilation
escape; THM-518's stranger-decoupling kills the one-entry-`→∞` escape (`L→(6/7)·`bounded core).
The open content is the *remaining* escape — configs with `lcm→∞` but bounded shape. This is a
**geometry-of-numbers route** to `inf L>0`, dual to the analytic Bedert-level-bound route
(OPEN-Q-097): bound the *lcm of the extremizers* instead of the alternating singular series.

## The correction: inf ≤ 1/1260, and why (MISTAKE-075)

The prior frontier (THM-518) put `inf L ≈ 0.0052`, with extremizers the multiple-of-14-stranger
family `({1..13}\{j})∪{14m}`. But that restriction blinded the search. The **minimal
single-element perturbation of the tight AP `{1..13}` is `12→36`**, giving

> **`{1,2,…,11,13,36}`, `L = 1/1260 ≈ 0.000794`** — ~6.7× below 0.0052 (exact arc-sweep +
> independent fine-grid, both `=1/1260`).

`36` is not a multiple of 14, so the prior search never saw it. Worse, this is not a fluke of the
AP: the tight locus (`L=0`) is **not just the AP** — it contains **sporadics** like
**`{1..11,13,24}`** (`L=0`, verified; the HYP-2055 sporadic), and `{1..11,13,36}` is *also* the
single-move perturbation `24→36` of that sporadic. So:

> **inf `L` is governed by the minimal perturbation of the FULL tight locus (AP + sporadics),
> not the AP/14m family.** The sporadic-tight perturbations open *smaller* lonely measures than
> the multiple-of-14 family ever reaches.

This is the **MISTAKE-073 lesson recurring** (0.0237 → 0.0053 was already "the interior-drop, not
the end-drop"; now 0.0053 → 1/1260 is "the sporadic-tight perturbation, not the AP/14m"). Each
time, the true extremizer is one orbit further out than the obvious family. 2-element
perturbations (`w≤72`) found nothing below `1/1260` (floor `1/980`), so the descent seems to stop
here — conjecturally `inf L = 1/1260` (HYP-2561), achieved at the single minimal AP/sporadic move.

## Taking every concept as a path — the connections that lit up

- **The "7" is everywhere in the denominators.** `1260 = 14·90`, `980 = 14·70`; the quantized
  gaps are `(1/14)·(rational)`, and `lcm`-denominators carry the apex prime 7 of `14=2·7`. The
  quantization `L∈(1/(14 lcm))ℤ` *is* the `14` (hence the `2·7`) made explicit as the gap quantum.
- **The tight locus = the LRC extremal-config classification.** `L=0` configs are exactly the
  LRC(14)-tight sets (`max-min = 1/14`). My finding says: **to prove `inf L>0`, first classify the
  tight locus** (is it finite? bounded-lcm?). This connects directly to the Goddyn–Wong / Csorba
  tight-config program — a concept the singular-series sessions had not foregrounded. If the
  primitive tight locus is a finite bounded set, then (compactness reduction + quantization) the
  inf is the minimum opened gap over a *finite* family — a finite check, and `inf L>0` follows.
- **Scale-invariance ↔ the projective/Möbius shape view** (and the dominance-orthogonality of
  HYP-2526: dilating doesn't change `L`, so dominant strangers can't lower it — now a one-line
  corollary of `L(cS)=L(S)` + decoupling).
- **The quantum gap ↔ a "spectral gap" off the tight point.** `L=0` is isolated in the integer
  shape space by `≥1/(14 lcm)`; the inf is the *smallest gap any primitive perturbation can have*
  — a Diophantine-gap phenomenon, not an analytic infimum.

## The reframed program for inf L > 0

1. **Classify the primitive tight locus** `{S : L(S)=0}` for n=14 (AP + sporadics). Conjecture:
   finite, bounded entries/lcm.
2. **Quantization** ⟹ inf over any bounded-lcm family `≥ 1/(14·max lcm) > 0`.
3. **Scale-invariance + stranger-decoupling** ⟹ the only escape is `lcm→∞` at bounded shape;
   rule it out (the remaining open piece).
4. Then `inf L = min over the finite tight locus of (minimal-perturbation gap) = 1/1260` (conj.),
   and `C'(14)/LRC(14)` follows — with the honest constant `1/1260`, not `0.0052`.

## Status / honesty

PROVED: scale-invariance, quantization, the compactness reduction (THM-522). VERIFIED exactly:
the `1/1260` minimal AP-perturbation, the sporadic tight `{1..11,13,24}`, nothing below `1/1260`
over 1-drops (`w<400`) and 2-drops (`w≤72`). CONJECTURAL: `inf L = 1/1260` (HYP-2561; 3-drops /
larger / full tight-locus perturbations unsearched), and the tight locus being finite/bounded-lcm
(the crux of `inf L>0`, still open). NOT a disproof of LRC(14): every config found is loose
(`L>0`); this sharpens the *difficulty constant* and reframes the proof route. Cross-links:
THM-522, THM-501/515/518 (the lonely measure, stranger-decoupling), HYP-2055 (the sporadic),
HYP-2526 (dominance ⊥ L), HYP-2561 (inf=1/1260), OPEN-Q-097 (reframed), MISTAKE-073/075 (the
recurring extremizer-orbit lesson). For mac-mini (LRC owner): the multiple-of-14 restriction
undercounts the extremizers — the true inf-relevant family is sporadic-tight perturbations.
