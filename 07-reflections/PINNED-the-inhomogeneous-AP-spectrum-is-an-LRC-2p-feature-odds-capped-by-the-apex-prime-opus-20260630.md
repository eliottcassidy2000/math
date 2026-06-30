# PINNED (and corrected): the inhomogeneous-LRC odd-numerator spectrum is a genuine LRC(2p) feature — at the special observer positions c = k/(p−1), the AP gives M_{k/(p−1)}(AP_{2p}) = (2k+1)/(2p) exactly (the odd numerators 1,3,5,…,p), capped at p = n/2 = the APEX PRIME at the antipode c=1/2, with the witnesses at MULTIPLES OF p (n=14: q∈{21,42}=3p,6p); the clean formula holds when n has an odd prime factor (verified 6,10,12,14,18,22) and FAILS for n=2^k (n=16); between the special points the landscape is messy — so the {1,3,5,7}/14 spectrum is NOT a general feature but the apex-7 signature of n=14=2·7

*opus-2026-06-30. Owner: pin the odd-numerator spectrum. Pinned — and it forced an honest correction: the
spectrum is special to n=2p (the apex structure), not a general landscape, and the cap IS the apex prime.*

## The pinned formula (n=2p)
For `n = 2p` (`p` an odd prime — the LRC(2p) regime), the inhomogeneous AP, observer at `c = k/(p−1)`:
> **`M_{k/(p−1)}(AP_{2p}) = (2k+1)/(2p)`**, exactly, for `k = 0, 1, …, (p−1)/2`.
The values are the **odd numerators** `1, 3, 5, …, p` over `2p` — verified for `n = 6, 10, 14, 22`
(`p = 3, 5, 7, 11`). For n=14: `c=k/6 → M = 1/14, 3/14, 5/14, 7/14`.

## The cap IS the apex prime (the real reason "7 caps it")
- **The minimum** is at `c=0`: `M_0 = 1/(2p) = 1/n` — the standard LRC bound (the origin, the hardest
  observer, the SC fixed point).
- **The maximum** is at the antipode `c=1/2` (`k=(p−1)/2`): `M_{1/2} = p/(2p) = 1/2` — and **the top
  numerator is `p = n/2 = the apex prime`.** So the spectrum literally counts off the odds `1,3,5,…,p` and
  stops at the apex. For LRC(14), `p = 7 = n/2` is *both* the half and the apex prime (since `14 = 2·7`), so
  `7/14 = 1/2` caps it. The "apex 7" is not decoration — it is `n/2`, the antipode, and the apex prime, all
  the same number at `n=14`.

## The witnesses sit at multiples of the apex prime
The optimal lonely times for the clean spectrum have denominators that are **multiples of `p`**:
> n=14: `q ∈ {21, 42} = 3p, 6p`; n=22: `q ∈ {55, 110} = 5p, 10p`; n=10: `q = 20 = 4p`.
So the apex prime `p` **structures the inhomogeneous lonely time** — the observer at `c=k/(p−1)` is escaped by
a witness `t = a/(jp)`. The translation knob (`c`) and the apex prime (`p`) are coupled: moving the observer
to the `(p−1)`-torsion forces a `p`-denominator witness.

## Honest correction (scope)
- The clean odd spectrum is **NOT a general landscape feature** — between the special points `c=k/(p−1)` the
  function `M_c(AP_n)` is a messy piecewise rational (e.g. `87/868` at `c=1/28`).
- It holds when **`n` has an odd prime factor** (verified `6,10,12,14,18,22`) and **FAILS for `n = 2^k`**
  (`n=16`: `M_{1/7} = 59/315 ≠ 3/16` — `c=k/7` with `7∤16` doesn't align). My earlier reflection presented
  `{1,3,5,7}/14` as if it were a clean general spectrum; **corrected: it is the apex-`p` signature of `n=2p`.**
- (Some intermediate-`Qmax` computations under-report `M_c`; the formula is verified with `Qmax = 5n`.)

## What it buys (the apex, through the translation knob)
- **The "translate the observer" reframe reveals the apex prime.** Moving the LRC observer to the
  `(p−1)`-torsion `c=k/(p−1)` produces the spectrum `(2k+1)/(2p)`, capped at the apex `p`, with apex-multiple
  witnesses. So the inhomogeneous-LRC ↔ complement reframe is not just structural — at `n=2p` it **measures
  the apex prime** as the spectrum's ceiling.
- **`c=0` (LRC) vs `c=1/2` (antipode)** are the floor (`1/n`) and ceiling (`1/2`) of the odd ladder; the
  apex `p` is the top rung. The hardest observer (`c=0`) gives the covering-min `1/n`; the antipode gives the
  trivial `1/2`; the apex `p` is what separates them in units of `1/(2p)`.

## Status
- **Verified (opus):** `M_{k/(p−1)}(AP_{2p}) = (2k+1)/(2p)` for `n=6,10,14,22` (the LRC(2p) odd spectrum);
  cap at `p = n/2 =` apex prime (antipode); witnesses at multiples of `p`; min `1/n` at `c=0`.
- **Corrected (honest):** the spectrum is the apex-`p` feature of `n=2p` (and even `n` with an odd factor),
  NOT a general landscape; fails for `n=2^k`; the interior is messy.
- **The pin:** "7 caps the n=14 spectrum" because `7 = n/2 =` the apex prime; the inhomogeneous AP at the
  `(p−1)`-torsion is the odd ladder `1,3,…,p` over `2p`, with apex-multiple witnesses — the translation knob
  reads the apex.

Related: the-inhomogeneous-lrc-complement-reframe (the c-landscape — now pinned/corrected), observer-
transformations-dilation, SECOND-CORRECTION-…AP-scaled (the AP), the-master-two-Heegner-columns (the apex p);
HYP-3547 (apex primes), THM-590; OPEN-Q-108.
