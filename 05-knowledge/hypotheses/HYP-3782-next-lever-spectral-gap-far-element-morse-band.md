---
id: HYP-3782
title: THE NEXT LEVER (SDP/Delsarte) HAS A SPECTRAL GAP -- loneliness is pointwise, so averaged/PSD certificates are blind (S54); the tight covering-min lower bound is COMBINATORIAL (the lazy-cut witness packing, HYP-3779). Working the HYP-3783 menu through E2 + F7(apex-7) + far-element resonance + Morse/band barrier. (1) FAR-ELEMENT RESONANCE: the construction killer n(n-1) = Phi6(n)-1 == -1 mod Phi6 for ALL n (verified 8,10,12,14) -- the far element sits at the "-1 slot" (the ceiling of the Stern-Brocot ray, the zeta_6/E2-hexagonal direction; Dedekind s(n,Phi6)->-1/12, S56). (2) MORSE/BAND LANDSCAPE: G(t)=min_v||v t|| (loneliness) for the construction is a POINTWISE SPIKE at the binding t*=n/Phi6 (global Morse max = M) atop a LADDER of lower local maxima (n=14: 0.0765, 0.0756, 0.0755, 0.0752,... = the radius-1 band-barrier critical points); the max sits at 1-n/Phi6 too (iota-symmetric, S55). (3) F7/SPECTRAL GAP: the Fejer-weighted AVERAGE of loneliness is 2.6-11x BELOW M (F_7: 0.0295, F_14: 0.0189, F_61: 0.0070 vs M=0.0765) and gets WORSE at higher degree -- the averaged/PSD lens is BLIND to the spike (S54). So a Delsarte/Fejer SDP has a spectral gap; the E2/Eisenstein bulk is regularizable (->-1/12, S56) but the apex-7/F7 cusp-form residual is the un-relaxable core the PSD relaxation cannot see. CONCLUSION: the tight lower-bound certificate is the COMBINATORIAL witness packing (lazy-cut), not the SDP -- the next lever is gap-limited by the pointwise nature of loneliness
status: MIXED (structural facts + honest assessment of the next lever). VERIFIED: far-element killer == -1 mod Phi6 (exact, all n); the Morse local-max ladder of the construction's loneliness landscape (n=14 grid); the Fejer-average spectral gap (2.6-11x below M, worsening with degree). ASSESSMENT (S54-backed, not a theorem): averaged/PSD (Delsarte/Fejer/theta) certificates are blind to the pointwise loneliness spike, so the SDP 'next lever' is expected gap-limited; the lazy-cut combinatorial witness packing (HYP-3779) remains the tight certificate. cvxpy unavailable (no full SDP tested); the Fejer-average is a naive spectral proxy (a proper Delsarte SDP could be tighter but, per S54, likely still gapped).
source: klein-2026-07-01-S63
depends_on:
  - HYP-3783   # the menu (this works its Lovasz/SDP item + Farey/E2/F7)
  - HYP-3766   # the witness is pointwise, the average is blind (S54) -- the spectral-gap principle
related:
  - HYP-3779   # lazy-cut combinatorial certificate (the tight one)
  - HYP-3768   # E2/Eisenstein bulk = regularizable, residual = cusp form (S56)
  - HYP-3767   # iota-symmetry (the landscape max at t* and 1-t*)
  - HYP-3763   # Steinhaus scaling (the far-element / band ladder)
  - HYP-3715   # Phi6 = zeta_6 hexagonal (the far-element -1 slot)
results:
  - 04-computation/covering_min_morse_band_e2_f7_klein.py
  - 05-knowledge/results/covering_min_morse_band_e2_f7_klein.out
---

# HYP-3782 — the next lever has a spectral gap; the far-element, the Morse band, and F7

## The next lever, assessed
HYP-3783 flagged Lovasz-theta/SDP (a tighter, PSD, multi-modulus relaxation) as the next lever after the
LP (too weak) and the Farey binding-first (too weak). This session assesses it through the requested
lenses -- E2, F7(apex-7), far-element resonance, Morse/band barrier -- and finds a **fundamental
spectral gap**: averaged/PSD certificates are blind to loneliness, which is pointwise (S54).

## (1) Far-element resonance: the killer sits at -1
The construction's killer `n(n-1)` equals `Phi6(n) - 1`, so **`n(n-1) == -1 (mod Phi6)`** for all `n`
(verified `8,10,12,14`). At the binding modulus `Phi6` the far element is at the `-1` slot -- the
ceiling of the Stern-Brocot ray `1/(n-1)`, the `zeta_6` / E2-hexagonal direction (HYP-3715), and the
Dedekind sum `s(n, Phi6) -> -1/12` (the E2 anomaly, S56). The far element is not "far" arithmetically:
it is the reflection anchor `-1` at the hexagonal modulus.

## (2) The Morse / band landscape
The loneliness `G(t) = min_v ||v t||` for the construction (`n=14`) is a **pointwise spike** at the
binding `t* = 14/183` (the global Morse maximum, `= M`), sitting on a **ladder of lower local maxima**:
`0.0765, 0.0760, 0.0756, 0.0755, 0.0752, 0.0747, ...` -- the radius-1 **band-barrier** critical points
(the moduli `n < D < 2n`). By `iota`-symmetry (S55) the global max also appears at `1 - 14/183 = 0.9235`.
So the landscape is: one sharp global spike (the binding) + a decaying ladder of band critical points.
The Morse structure is the loneliness function's, and its top value is the covering-min.

## (3) F7 / apex-7 and the spectral gap
Test whether a low-degree positive-definite (Fejer) minorant -- the natural Delsarte/SDP certificate,
with `F_7` the apex-7 kernel -- can certify `M >= n/Phi6`. The Fejer-weighted AVERAGE of the loneliness
is far below `M`, and worsens with degree:
```
  F_7  avg 0.0295   F_14 avg 0.0189   F_61 avg 0.0070    vs    M = 14/183 = 0.0765
  gap  2.6x                                       10.9x
```
Higher-degree Fejer concentrates the average and sees the spike LESS. This is the **spectral gap**: an
averaged / PSD lens is blind to the pointwise loneliness spike (S54, "the average is blind"). A full
Delsarte/theta SDP optimizes the kernel and could beat the naive average, but the same principle
(loneliness is pointwise, not spectral) predicts it remains gapped -- and cvxpy is unavailable to test
a full SDP here.

## Synthesis (E2 bulk vs F7 residual; combinatorial vs spectral)
- **E2 / Eisenstein = the regularizable BULK**: the Dedekind sum `s(n,Phi6) -> -1/12 = zeta(-1)` (S56),
  the far-element `-1` slot, the hexagonal `zeta_6` -- all the "average" structure, which a PSD/spectral
  certificate CAN see.
- **F7 / apex-7 cusp-form = the un-relaxable RESIDUAL**: the pointwise spike, the genus-1 obstruction
  `f_14` (S56/HYP-3768) -- what the averaged/PSD relaxation CANNOT see.
- So the tight lower-bound certificate is the **COMBINATORIAL witness packing** (the lazy-cut, HYP-3779:
  a finite list of ~thousands of lonely-witness cuts), NOT a Delsarte/Fejer/theta SDP. The next lever is
  gap-limited by the pointwise nature of loneliness.

## Net
Working the SDP lever through E2/F7/far-element/Morse-band: the far element sits at the `-1` slot of the
hexagonal binding modulus; the loneliness landscape is a pointwise spike (global Morse max = the
covering-min) on a band-barrier ladder; and the Fejer/PSD average is 2.6-11x below `M` -- a spectral gap
confirming (S54) that averaged/PSD certificates are blind to the spike. The regularizable E2/Eisenstein
bulk is spectrally visible, but the apex-7/F7 cusp-form residual is not -- so the tight certificate is
the combinatorial witness packing (lazy-cut), and the SDP next lever is gap-limited. The honest lesson:
loneliness is pointwise; the certificate must be a list of points (witnesses), not a spectrum.
