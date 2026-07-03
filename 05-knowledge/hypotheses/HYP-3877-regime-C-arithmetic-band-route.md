---
id: HYP-3877
title: THE ARITHMETIC BAND ROUTE -- the NEAR-EQUAL far case (regime C c=7 AND c>=8 near-equal) closes by a small FINITE denominator band via the small-q witness t=a/q, NOT a floor. [CORRECTED S22: honest band-below is {15..~50} not {15..33} -- see MISTAKE-095; the {15..33} below held only for generic/unaligned drifts.] At t=a/q (q in 15..27) the danger residues are EXACTLY {0,1,q-1}; a covering family is 1/14-lonely at a/q iff va !in {0,+-1} mod q for every speed v. NEAR-EQUAL far blocks are LOOSE and the far cluster is a SHORT AP mod q (span too small to block the band): worst band-q=33 over thousands of families, 0 failures, nfar=7..13, span up to w1. The Hunter/drift-floor FAILS here not because the family is tight but because it is an ANALYTIC (measure) method applied to an ARITHMETIC (small-q) case; the right tool is the small-q witness directly. IMPORTANT LIMITATION (S21, corrected): the band is NOT uniform over ALL hge7 -- COMPOSITE band-blocker families (far runners divisible by many small moduli, e.g. 17.19.23=7429, 29.31=899) block {15..33} by DIVISION and need q=34..53, GROWING ~log(max speed) (max q = 38,46,47,53,53 at mag 10^3..10^7). So a uniform finite band exists ONLY with an a-priori speed-magnitude bound (regime C has one: w1<7392 => band {15..~46}). Band-blockers are SPREAD/composite (ratio ~4), NOT near-equal, so the near-equal claim survives.
status: CORRECTED (see MISTAKE-095). The `{15..33}` figure holds only for GENERIC (unaligned) near-equal drifts; ALIGNED near-equal drifts (each runner `= q*round(N/q)`, span-ratio ~1.00, `≡0 mod q`) are band-blockers with witness `q` up to 47-49. The HONEST band-below is `{15..~50}` for `max-speed < 22638` (0 failures over ~164k adversarial families incl. aligned near-equal + spread). NO uniform band even for near-equal: witness `q ~ log(max-speed)` ([[HYP-4040]] proves `q(S_X) > X`). So: PARTIALLY-CONFIRMED that regime C / bounded-magnitude hge7 closes by a FINITE arithmetic band (size `~50`, was over-stated as `33`); REFUTED as a fixed/uniform band. The mechanism (small-q witness, danger set `2*ceil(q/14)-1` residues) and the magnitude-split architecture SURVIVE. Contradicts kps-S25 MISTAKE-072 'no floor closes regime C' (correct -- it's not a floor, it's the arithmetic witness), but the honest closure is magnitude-split, not a free lunch.
source: mac-mini-2026-07-03-S21
related:
  - HYP-+2876  # lcm families {1..11,13,lcm(2..X)} = the band-blockers; unbounded witness denom (REFUTED uniform band). THIS is the residual, handled by magnitude split.
  - HYP-3737   # covering-min forced; the band (14, 2n-2]=(14,26] is exactly this denominator band
  - HYP-3983   # kps-S27 mutual-independence closer; near-equal = where independence FAILS = where THIS arithmetic route succeeds (complementary)
  - HYP-3982   # kps star-safe (the c<=7 measure route); this is the ARITHMETIC complement for near-equal
  - HYP-3981   # kps c=7 trichotomy: regime C = the open near-equal-small case THIS closes
  - HYP-3901   # deep-cluster renormalization (n->8); the large-magnitude route for the >22638 side
  - MISTAKE-072 # kps 'no floor closes regime C' -- true; the band is not a floor, it's the arithmetic witness
results:
  - 04-computation/regimeC_band_adversarial_macmini_20260703.py     # regime C worst q=33, 0 fail
  - 04-computation/band_route_boundary_macmini_20260703.py          # near-equal vs dilated-AP boundary
  - 04-computation/band_route_fullspan_macmini_20260703.py          # full regime-C span (span up to w1)
  - 04-computation/hge7_gcd1_band_closure_macmini_20260703.py       # gcd=1 hge7 random: worst q=32
  - 04-computation/band_blocker_construction_macmini_20260703.py    # CONSTRUCTED band-blockers need q=34-36
  - 04-computation/max_lonely_q_hge7_macmini_20260703.py            # max q grows ~log mag: 38,46,47,53,53
  - 05-knowledge/results/max_lonely_q_hge7_macmini_20260703.out
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

## The band-blocker limitation (S21, the honest boundary of the route)
The band `{15..33}` is NOT uniform over all `hge7`. A **band-blocker** is a covering `gcd=1` family whose
speeds are divisible by many small moduli: taking speeds `= 17.19.23 = 7429`, `29.31 = 899`, `25.32 = 800`,
`27.16 = 432`, ... makes EVERY `q in {15..33}` divide some speed (`va ≡ 0`, danger), so no band `q <= 33`
is lonely -- the smallest lonely `q` is then `34, 35, 36`. Searching to maximize the smallest lonely `q`
over covering `gcd=1 hge7` families:
```
    max speed <= 10^3   10^4   10^5   10^6   10^7
    smallest lonely q =  38     46     47     53     53
```
So the smallest lonely `q` GROWS (~`log` of the magnitude), with NO universal finite bound. To need a large
`q`, a family must block `{2..q-1}` by division -> its `<=13` speeds must be divisible by all primes up to
`q-1`; `13` speeds carry `O(13 * (log mag / log p))` prime factors, so `q` grows only logarithmically -- but
it does grow. Band-blockers are SPREAD/composite (ratio `~4`), NOT near-equal, so the regime-C claim stands.

## Impact (corrected)
- **Regime C (c=7 near-equal) CLOSES** by `{15..33}` -- resolving the case kps-S25 found no floor could
  close, because the failure was analytic-method-on-arithmetic-case, not tightness. Regime C has an a-priori
  magnitude bound (`w1 < 7392`), so its band is finite/decidable.
- **c>=8 NEAR-EQUAL CLOSES** too (same band, nfar up to 13).
- **NOT a uniform closure of hge7.** The SPREAD/composite band-blockers need `q` up to `~53` (growing with
  magnitude); those are regime B, still the pair-floor's job -- unless an a-priori speed-magnitude bound is
  established (then the band `{15..Q(B)}` is finite and closes them by `native_decide`). This magnitude bound
  is the real missing piece (see [[HYP-3737]] -- covering-min forced at bounded denominator).

## Lean feasibility
The band condition is DECIDABLE per family (`exists q in {15..Q}, a: forall v, va mod q not in {0,1,q-1}`).
For a magnitude-bounded slice it reduces to residue classes mod `lcm(15..Q)` -- a `native_decide` census like
`hwindow`. Regime C (near-equal, `w1<7392`) is such a slice with `Q ~ 46`.

## The magnitude-split synthesis (S21, the constructive way this closes hge7)
The band-blockers are EXACTLY the lcm families of [[HYP-+2876]] (`{1,..,11,13,lcm(2..X)}`) whose witness
denominator is unbounded in `X`. HYP-+2876 REFUTED a UNIFORM finite band -- correct, and my growth table
(`q -> 53` at mag `10^7`) is its quantitative face. But those families have max speed `= lcm(2..X) -> infinity`.
So split hge7 by MAGNITUDE, not with one band:
- **max speed `< 22638`** (below kps's singles-bound threshold): the band `{15..50}` closes EVERY covering
  `gcd=1 hge7` family -- `189128` adversarial families (incl. lcm/band-blockers up to the threshold), `0`
  failures, worst smallest-lonely-`q = 50`. Finite because `lcm(2..11) = 27720 > 22638` bounds the lcm
  parameter `X <= 10`, so the witness denominator is bounded there.
- **max speed `> 22638`**: kps's SINGLES bound (regime B) / `cite_cluster7` (regime A). The unbounded-witness
  lcm families live here.
- **NEAR-EQUAL far, ANY magnitude**: band `{15..33}` regardless of `w1` -- a `k`-consecutive cluster is `k`
  consecutive residues mod `q`, and `k <= 13 < q-3` fits in the safe set `{2,..,q-2}` for `q >= 17`. So
  near-equal never needs the magnitude split at all.

Net structure: `{max speed < 22638}` -> band `{15..50}`;  `{max speed > 22638}` -> {near-equal: band `{15..33}`;
spread: singles/cluster}. Every covering `gcd=1 hge7` family is in one bucket -> no gap. This is a CLOSURE
STRUCTURE (magnitude split), consistent with HYP-+2876 (no *uniform* bound). Rigor still needed: the band-below
as a finite `native_decide` (the residue space mod `lcm(15..50)` is large -- needs a smarter finite argument or
the CRT far-blocking count), and kps's singles-above (joint-measure treatment).

## Honest scope
SOLID for NEAR-EQUAL (regime C + c>=8 near-equal): band `{15..33}`, `~11000` families, `0` failures, clean
dilated-AP boundary, CRT far-blocking mechanism. REFUTED as a uniform closure of all `hge7`: composite
band-blockers grow the smallest lonely `q` (`->53` at mag `10^7`). The clean full closure of `hge7` via the
band requires an a-priori speed-magnitude bound; without it, the SPREAD/composite part stays with the
pair-floor. Net: a genuine NEW tool that dissolves the specific stuck cases (regime C, c>=8 near-equal), plus
a sharp identification of the residual obstruction (band-blockers) and what would remove it (magnitude bound).
