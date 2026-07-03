---
id: HYP-3877
title: THE ARITHMETIC BAND ROUTE closes the NEAR-EQUAL far case (regime C c=7 AND c>=8 near-equal) by a small FINITE denominator band {15..33} -- NOT a floor. At t=a/q (q in 15..27) the danger residues are EXACTLY {0,1,q-1}; a covering family is 1/14-lonely at a/q iff va !in {0,+-1} mod q for every speed v. Near-equal far blocks are LOOSE (M>=0.13 >> 1/14, never in the tight {AP,GW} locus), so a small-q witness always exists (worst band-q=33 over thousands of families, 0 failures, nfar=7..13). The Hunter/drift-floor FAILS here not because the family is tight but because it is an ANALYTIC (measure) method applied to an ARITHMETIC (small-q) case; the right tool is the small-q witness directly. Boundary: the tight DILATED AP {d,2d,..,13d} (SPREAD, M=14/183) needs a LARGE q (28,42,70,196 for d=2,3,5,14) -- that is regime B, handled by the pair-floor.
status: PARTIALLY-CONFIRMED (strong numeric: ~11000 adversarial covering families across regime C + c>=8 near-equal + all far scales to 50000, worst band-q=33, ZERO failures; dilated-AP control needs large q as predicted). A clean structural ROUTE + proof mechanism (CRT bounds far-blocking), NOT yet a theorem. Contradicts kps-S25 MISTAKE-072 'no floor closes regime C' -- correct, because it is NOT a floor.
source: mac-mini-2026-07-03-S21
related:
  - HYP-3737   # covering-min forced; the band (14, 2n-2]=(14,26] is exactly this denominator band
  - HYP-3551   # covering-min 14/183 (the TIGHT dilated-AP case the band does NOT close -> pair-floor)
  - HYP-3982   # kps star-safe (the c<=7 measure route); this is the ARITHMETIC complement for near-equal
  - HYP-3981   # kps c=7 trichotomy: regime C = the open near-equal-small case THIS closes
  - MISTAKE-072 # kps 'no floor closes regime C' -- true; the band is not a floor, it's the arithmetic witness
results:
  - 04-computation/regimeC_covering_families_macmini_20260703.py
  - 04-computation/regimeC_band_argument_macmini_20260703.py
  - 04-computation/regimeC_band_adversarial_macmini_20260703.py
  - 04-computation/band_route_boundary_macmini_20260703.py
  - 05-knowledge/results/band_route_boundary_macmini_20260703.out
---

# HYP-3877 -- the arithmetic band route for the near-equal far case

The `hge7` leg (>=7 far runners, compressed) splits at the far-cluster GEOMETRY, not just the count:
- **NEAR-EQUAL far** (small drift): this file -- the arithmetic band `{15..33}`.
- **SPREAD far** (dilated-AP-like, tight): regime B -- the pair-floor (kps/klein).
The near-equal case includes BOTH regime C (`c=7` near-equal) AND `c>=8` near-equal -- two of the fleet's
open cruxes -- and both close by the SAME small-denominator witness.

## The mechanism (exact)
At `t = a/q`, `gcd(a,q)=1`, a runner `v` is `1/14`-lonely iff `min(va mod q, q - va mod q) >= q/14`. For
`q in {15,...,27}` this is EXACTLY `va mod q not in {0, 1, q-1}` (the danger set is 3 residues, since
`ceil(q/14)-1 = 1`). So a covering family is `1/14`-lonely at `a/q` iff
```
    v*a  ≢  0, +1, -1   (mod q)   for every speed v.
```
`va ≡ 0` iff `q | v` (`gcd(a,q)=1`). Covering only forces coverage of `q <= 14`, so a band modulus
`q in {15..27}` may divide NO speed -- then only the `+-1` conditions remain, easily met by a generic `a`.

## Why near-equal closes and dilated-AP does not
- **Near-equal far blocks are LOOSE.** Numerically `M >= 0.128 >> 1/14` for every near-equal family
  (nfar `7..13`, far scale to `50000`); they are never in the tight locus `{AP, GW}` (which has ratio `~13`,
  not near-equal). A loose family's covering-min is achieved at a SMALL denominator -> the band closes it.
  Worst `band-q = 33` over `~11000` adversarial covering families (highly-composite + strided `w1`), `0`
  failures.
- **The dilated AP `{d,2d,...,13d}` is TIGHT and SPREAD.** `M = 14/183` (dilation-invariant), achieved at
  `q ~ 89`. The band `{15..33}` does NOT close it: smallest band-`q` = `28, 42, 70, 196` for `d = 2,3,5,14`
  (grows with `d`). This is regime B, handled by the pair-floor / covering-min -- NOT the band.

## The proof mechanism (sketch)
To FAIL the band, a family must have EVERY `q in {15..33}` block it (`q | some speed`, or `+-1` unavoidable).
The `<=6` near runners (`<=22`) block `<=6` band moduli (their own values `15..22`). The near-equal far
cluster (span `<= D`, `D` small) blocks `q` only if `q | (w1 + j)` for some `j` in a window of width `D`;
different band PRIMES `23, 29, 31, ...` need distinct multiples inside that small window, and by CRT
covering many of them forces `w1` far beyond the near-equal-small range. So a near-equal cluster cannot
block the whole band -> some `q` is free. (The dilated AP evades this because its SPREAD far runners
`d,2d,...,13d` supply many divisors across the band at once -- exactly the tight case the band excludes.)

## Impact
- **Regime C (c=7 near-equal) CLOSES** by the finite band `{15..33}` -- resolving the case kps-S25 found no
  floor could close. The Hunter/floor methods failed because they are analytic on an arithmetic case.
- **c>=8 NEAR-EQUAL CLOSES** too (same band, nfar up to 13) -- the c>=8 crux for near-equal blocks.
- Remaining: SPREAD far (dilated-AP / GW, tight) via the pair-floor (regime B, kps/klein); and the exact
  near-equal/spread transition (drift threshold where the band hands off to the pair-floor).

## Lean feasibility
The band condition is DECIDABLE per family (`exists q in {15..33}, a: forall v, va mod q not in {0,1,q-1}`).
Regime C is a FINITE family set (`w1 <= 7392`); the near-equal case reduces to residue classes mod
`lcm(15..33)`. A census-style `native_decide` over the (finite) near-equal residue patterns could formalize
it -- the same shape as the `hwindow` census.

## Honest scope
Strong numeric route (`~11000` families, `0` failures, clean dilated-AP boundary), with a proof mechanism
(CRT far-blocking bound). NOT yet a theorem. The claim is scoped to NEAR-EQUAL far; the spread/tight case is
explicitly OUT (regime B). The exact drift threshold and the rigorous CRT bound are the next steps.
