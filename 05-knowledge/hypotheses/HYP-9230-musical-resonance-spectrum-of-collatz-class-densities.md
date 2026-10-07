---
id: HYP-9230
title: "Musical resonance law: the density of any Collatz class (or end-saturated class) along log2 x has a log-periodic Fourier spectrum supported near frequencies f (cycles per octave) with small ||f log2 3||: the convergent and semiconvergent denominators of log2 3 (5, 12, 41, 53, 94, 147, 200, 253, 306, 359, 665, ...) and harmonics such as 106 = 2*53; mode f decays like exp(Re w(theta_f) log2 x) with drift Im w/2 pi, where w is the root (continued from 0) of exp(2 pi i theta + w log2(3/2)) + exp(-w) = 2, theta = f log2 3 - round(f log2 3); so the approach to the uniformity of THM-4590 runs on the tuning spectrum of the fifth 3/2"
status: >
  OPEN (the existence of the law and of nonzero amplitudes). NUMERICAL, exact class membership for all n <= 2^28, both signs;
  windowed DFT with 2048 windows per octave, plus an independent windowless log-measure DFT up to f = 1100 (audit A2).
  Spectrum for k >= 20: the -1 basin's largest mode is f = 106 (0.075 at k = 28), then 53, 12, 253, 265, 147, 94, 306, 41. For entry via 85
  at k = 28 the order is 53, 106, 147, 159, 200, 306, 253, 41, 94, 12. A_665 is nonzero (about 0.009 at k = 24).
  The DECAY LAW is the exact characteristic root of the transfer psi(u) = 1/2 psi(u + log2(3/2)) + 1/2 psi(u - 1) (audit A2). It matches the
  measured decay per octave for every measurable mode (f = 12: 0.105 against 0.105/0.106; 41: 0.085 against 0.088/0.084; 65: 0.124 against
  0.118/0.125; 147: 0.044 against 0.046/0.045; 265: 0.075 against 0.073/0.082) and the phase drift (f = 12: 0.0327 against 0.0324
  cycles/octave). The earlier "factor-2 discrepancy" was the failure of the quadratic (Gaussian) approximation: RESOLVED.
  HEURISTIC: the mechanism is the walk of log2 x by fifths up and octaves down (the irrational rotation in THM-4590 (3)).
  PRIOR ART: a log-periodic correction is known in Wirsching's predecessor-density program, Berg-Krueppel 1998 (an infinite product) and
  Tavares arXiv:2608.27617 (a Fourier series, proved non-constant); also Benford for 3x+1 (Kontorovich-Miller; Lagarias-Soundararajan).
  No earlier spectrum of CLASS densities was found.
source: mac-mini-2026-10-07-golden, 05-knowledge/results/golden_collatz_resonance_20261007.md
related:
  - 01-canon/theorems/THM-4590-collatz-classes-are-equidistributed-residues-windows-slowly-varying-density.md (asymptotic uniformity; this hypothesis is its rate)
  - 01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md (Ellison's bound; see the common-cause note below)
  - Berg-Krueppel (1998); Tavares, arXiv:2608.27617; Kontorovich-Miller, Acta Arith. 120 (2005); Lagarias-Soundararajan, J. London Math. Soc. 74 (2006)
scripts:
  - 04-computation/experiments/golden_20261007_lead/ (negfourier*.c, posentry85*.c + .out)
  - 04-computation/experiments/golden_20261007_audits/A2/ (spec.c, roots2.py, decayfit.py + .out)
---

# HYP-9230 — the musical resonance spectrum of Collatz classes

**Statement.**
* Let `B` be a Collatz class or an end-saturated class (THM-4590 (5)). Let `β(s, φ)` be its density at log-position `s + φ`, where `φ` is the fractional part of `log2 n`, at scale `2^s`.
* Its Fourier coefficients satisfy

      β̂(s, f) = A_f · exp( w(θ_f)·s ) · (1 + o(1)),   θ_f = f log2 3 − round(f log2 3).

* Here `w(θ)` is the root, continued from `w(0) = 0`, of

      exp(2πiθ + w·log2(3/2)) + exp(−w) = 2 .

  `Re w < 0` is the decay per octave and `Im w/2π` is the phase drift.
* For small `θ`: `Re w = −2π²·27.97·θ² + O(θ⁴)`, the variance of odd steps per octave being 27.97. The effective constant `−Re w/θ²` falls from 533 (`f = 53`) through 273 (`f = 12`) to 97 (`f = 17`).
* The amplitudes `A_f` are nonzero exactly near small `θ_f`: the convergent and semiconvergent denominators of `log2 3` and their harmonics.

**The spectrum, with its musical names.**

| `f` | `θ_f` | approximation | tuning |
|---|---|---|---|
| 12 | 0.0196 | `2^19 ≈ 3^12` (Pythagorean comma) | 12-TET |
| 41 | −0.0165 | `2^65 ≈ 3^41` | 41-TET |
| 53 | 0.0030 | `2^84 ≈ 3^53` (Mercator's comma) | 53-TET |

* Further strong modes:
  * 94, 147, 200, 253 and 359 (semiconvergents between 65/41 and 485/306);
  * 306 (`485/306`) and 665 (`1054/665`, `θ = 6.3·10^-5`);
  * harmonics 106 = 2·53, 159, 265 and 318;
  * 65 = 53 + 12.
* Modes 5 and 17 are tiny at large `k` (≤ 0.0003 at `k = 28`).

**Evidence (NUMERICAL).**

| class, scale | top modes (amplitude) |
|---|---|
| −1 basin, `k = 18` | 12 (.101), 106 (.092), 53 (.061), 265 (.060), 94, 171, 65, 147, 253 |
| −1 basin, `k = 28` | 106 (.075), 53 (.057), 12 (.035), 253, 265, 147, 94, 306, 41 |
| entry via 85, `k = 28` | 53 (.031), 106 (.025), 147, 159, 200, 306, 253, 41, 94, 12 |

* The entry-via-85 class has mean density 0.0236, and its 53-mode amplitude exceeds that mean: about 53 bands per octave.

**Decay per octave over `k = 16..28`:**

| `f` | −1 basin | entry 85 | exact root | Gaussian |
|---|---|---|---|---|
| 12 | .105 | .106 | .105 | .211 |
| 41 | .088 | .084 | .085 | .151 |
| 65 | .118 | .125 | .124 | .281 |
| 147 | .046 | .045 | .044 | .061 |
| 265 | .073 | .082 | .075 | .125 |

**Ellison's exceptions: a common cause, not evidence.**
* Ellison's (16,10), (19,12) and (27,17) and the low modes 10 = 2·5, 12 and 17 are the same near-coincidences `2^x ≈ 3^y`: both lists are small values of `‖y log2 3‖`.
* The dominant modes at large scale (53, 106, 253, 306) are not Ellison exceptions.
* The identification is correct arithmetic, but it is a shared continued fraction, not an independent confirmation.

**What would follow.**
* The approach to THM-4590's uniformity is non-uniform in frequency. Its slowest modes sit at the convergent denominators `q_n` of `log2 3`, with time constant `≈ 1/(−Re w(θ_(q_n)))` octaves, between `4·10^5` and `10^6` octaves for 1054/665.
* So any quantitative form of THM-4590 must depend on the irrationality measure of `log2 3`.

**Open.**
* A proof of the transfer law for class densities: the walk model is heuristic for actual classes.
* A characterization of which `A_f` vanish.
* The relation to the Berg–Krüppel / Tavares log-periodic predecessor-density correction.
