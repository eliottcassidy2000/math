---
id: HYP-9230
title: "Musical resonance law: the density of any Collatz class (or end-saturated class) along log2 x has a log-periodic Fourier spectrum supported near frequencies f (cycles per octave) with small ||f log2 3||, namely the denominators of the convergents and semiconvergents of log2 3 (5, 12, 41, 53, 306, 665, 15601, ...; continued fraction [1;1,1,2,2,3,1,5,2,23,2,2,1,1,55,...]: the 5-, 12-, 41-, 53-TET temperaments) and their sums; mode f decays like exp(-C ||f log2 3||^2 log2 x) with C between about 280 and 550 (2 pi^2 x the odd-step variance 27.9 per octave gives 551); so the approach to the uniformity of THM-4590 is controlled by the Diophantine approximation of log2 3 (Ellison's exceptions (16,10), (19,12), (27,17) are the modes 5, 12, 17)"
status: >
  OPEN. NUMERICAL, exact class membership for all n <= 2^28, both signs, with Fourier analysis in log2 position over 2048 windows per octave.
  The dominant modes are the same for the negative-cycle basins and the positive-integer entry classes:
  53 (theta = 0.0030: amplitude 0.0595 -> 0.0570 for the -1 basin, 0.0332 -> 0.0306 for entry via 85, over k = 16..28),
  12 (theta = 0.0196: 0.125 -> 0.035 and 0.026 -> 0.007), 41 (theta = -0.0165), 65 = 53 + 12, 29 = 17 + 12, 24, 36, 17, 5.
  The decay rates are ordered by theta^2. The Gaussian model with the measured odd-step variance 27.9 per octave predicts
  C = 2 pi^2 x 27.9 = 551. It fits f = 53 (predicted 0.94, observed 0.92 over 12 octaves) but overestimates f = 12 and f = 41 by a factor of about 2: OPEN.
  HEURISTIC mechanism: the landing phase after k odd steps is k log2 3 mod 1, i.e. the irrational rotation used in THM-4590 (3);
  this is the Benford mechanism of Kontorovich-Miller and Lagarias-Soundararajan.
source: mac-mini-2026-10-07-golden, 05-knowledge/results/golden_collatz_resonance_20261007.md
related:
  - 01-canon/theorems/THM-4590-collatz-classes-are-equidistributed-residues-windows-slowly-varying-density.md (asymptotic uniformity; this hypothesis is its rate)
  - 01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md (Ellison's bound and exceptions 13, 14, 16, 19, 27)
  - Kontorovich-Miller, Acta Arith. 120 (2005); Lagarias-Soundararajan, J. London Math. Soc. 74 (2006) (Benford's law for 3x+1)
scripts:
  - 04-computation/experiments/golden_20261007_lead/ (negfourier.c, negfourier70.c, posentry85.c, posentry85_70.c + .out)
---

# HYP-9230 — the musical resonance spectrum of Collatz classes

**Statement.** Let `B` be a Collatz class or an end-saturated class (THM-4590 (5)). Let `β(s, φ)` be its density in the window at log-position `s + φ`, with `φ ∈ [0,1)` the fractional part of `log2 n`, at scale `2^s`. Then the Fourier coefficients satisfy

    | β̂(s, f) | = A_f · exp( −(C + o(1)) ‖f log2 3‖² s ) ,

with `A_f` depending on the class and on small-scale data. The persistent modes are the best rational approximations `p/f` of `log2 3`.

**The spectrum, with its musical names.**

| `f` | `‖f log2 3‖` | approximation | tuning | Ellison exception |
|---|---|---|---|---|
| 5 | 0.0752 | `2^8 ≈ 3^5` | 5-TET | `(16, 10) = 2·(8, 5)` |
| 12 | 0.0196 | `2^19 ≈ 3^12` (Pythagorean comma) | 12-TET | `(19, 12)` |
| 17 | 0.0556 | `2^27 ≈ 3^17` | 17-TET | `(27, 17)` |
| 29 | 0.0361 | `2^46 ≈ 3^29` | 29-TET | — |
| 41 | 0.0165 | `2^65 ≈ 3^41` | 41-TET | — |
| 53 | 0.0030 | `2^84 ≈ 3^53` (Mercator's comma) | 53-TET | — |
| 306, 665, 15601 | 1.5·10^-3, 6.3·10^-5, 2.6·10^-5 | `485/306`, `1054/665`, `24727/15601` | — | — |

**Evidence (NUMERICAL, `n ≤ 2^28`).**

| class | scale `k` | 53 | 12 | 41 |
|---|---|---|---|---|
| −1 basin (mean 0.327) | 16 | 0.0595 | 0.1248 | 0.0404 |
| −1 basin | 28 | 0.0570 | 0.0354 | 0.0141 |
| entry via 85 (mean 0.0236) | 16 | 0.0332 | 0.0262 | 0.0326 |
| entry via 85 | 28 | 0.0306 | 0.0073 | 0.0118 |

* The entry-via-85 class is concentrated in about 53 bands per octave: its 53-mode amplitude exceeds its mean.
* Per-octave decay factors are about 0.993 (53), 0.917 (41), 0.90 (12), 0.885 (65), 0.82 (24, 29), 0.78 (17) and 0.72 (36), ordered by `θ²`.

**What would follow.**
* The approach to THM-4590's uniformity is non-uniform in frequency.
* Its worst modes sit at the convergent denominators `q_n` of `log2 3`, with time constant `≈ 1/(C ‖q_n log2 3‖²)` octaves. For the convergent 1054/665 (partial quotient 23, `θ = 6.3·10^-5`) that is between `4·10^5` and `10^6` octaves.
* So any quantitative form of THM-4590 must depend on the irrationality measure of `log2 3`. This is the place where Baker-type bounds and Ellison's exceptions enter the geometry of Collatz classes.

**Open points.**
* The factor-2 discrepancy in `C` at `f = 12` and `f = 41`.
* A proof of the decay law for one mode.
* Whether `A_f` is nonzero for `f = 665`.
