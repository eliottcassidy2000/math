---
source: opus-2026-07-23-S4 (harvesting the owner's artanh-snippet technique for LRC(14))
status: >
  CONCRETE RESULT. The snippet's exact technique (truncate a convergent sum + rigorously
  bound the tail) is harvested into a float-free CERTIFIED two-sided bound on THM-729's
  density second moment Q_s = sum_l |U_s(lw)|^2/l^2. This turns THM-729's central EMPIRICAL
  claim (Q_s=O(diam), klein-S280 float grid) into a rigorous certificate on the tested box,
  and supplies a grid-free (exact-endpoint) method that is provably necessary as clustering
  grows. Does NOT close THM-729's open ASYMPTOTIC piece (uniform Q_s=o(r^2)); it certifies
  the finite box and provides the Lean-ready tool.
tags: [lrc14, thm-729, density-route, second-moment, certified, float-free, harvest, artanh, snippet-transfer]
related: [THM-729, THM-727, THM-728, HYP-9023, HYP-6415]
---

# Harvesting the artanh technique: a certified Q_s for the THM-729 density route (opus-S4)

## What was harvested

The owner's external LRC-proof snippet certifies a transcendental inequality by the device
`truncate the odd-power series + bound the tail geometrically` (`t^5/(5(1-t^2))`). kps/klein
established (S406) that the tempting `log R`/entropy transfer is UNSOUND for LRC (int(M log R)
is signed; log Z is 2nd-order flat), so the snippet's value to LRC is the **certified-bound
engine**, not a new soundness route. Here that engine is pointed at the one genuinely
transcendental, genuinely SOUND quantity on the density route: THM-729's

```
Q_s = sum_{l>=1} |U_s(l w)|^2 / l^2,   U_s(N) = sum_{p in dR_s} eps_p e(-2*pi*i*N*p).
```

`Q_s` is a positive convergent sum (no signed-functional problem), and `|U_s(lw)|` is a sum of
`M` unit vectors, so the snippet's move applies verbatim:

```
sum_{l>L} |U_s(lw)|^2/l^2  <=  M^2 * sum_{l>L} 1/l^2  <  M^2 / L      (geometric/integral tail)
```

## Two upgrades over the float grid (klein-S280)

1. **Exact arc endpoints.** The 7-sector map `sec(e,x)=floor(7 frac(e x))` jumps only at
   `x = i/(7e)`, so every breakpoint of `R_s` is an exact rational with denominator `| 7 lcm(E)`.
   No grid, hence no missed tiny arcs. (Validated: on `E=[0..6]` the exact and grid endpoint
   sets and `Q_s` agree to grid-resolution error ~0.01-0.17 per sector.)
2. **Rigorous head.** Each phase `N*p` is reduced EXACTLY mod 1 (`p` rational), then the head
   `sum_{l=1}^L |U_s(lw)|^2/l^2` is evaluated in IEEE double with accumulated rounding error
   provably `< M^2 * 2^-50 < 1e-9`; small clusters are independently confirmed with validated
   interval arithmetic (`mpmath.iv`). => certified `Q_s in [head - M^2 2^-50, head + M^2/L]`.

## Result (`lrc14_second_moment_certified_opus_S4.py`, `w=997, L=3000`)

Certified two-sided `Q_s` (max over `s`) and `Q_s/diam`, float-free:

```
cluster E                     diam  min arc width  grid?   certified Q_s        Q/diam
[0,1,2,3,4,5,6]                  6    4.8e-3        OK    [ 19.448, 19.460]   [3.241,3.243]
[0,1,2,3,4,5,12]                12    2.4e-3        OK    [ 17.181, 17.246]   [1.432,1.437]
[0,1,2,3,4,5,25]                25    1.4e-3        OK    [ 41.824, 41.872]   [1.673,1.675]
[0,1,2,4,8,16,32]               32    4.5e-3        OK    [ 16.219, 16.352]   [0.507,0.511]
[0,3,7,15,30,55,90]             90    1.4e-4        OK    [121.936,122.368]   [1.355,1.360]
[0,5,13,28,54,88,140]          140    3.8e-5        OK    [217.017,218.942]   [1.550,1.564]
[0,10,27,55,99,150,199]        199    4.8e-6        OK    [216.567,220.748]   [1.088,1.109]
```

- **`Q_s=O(diam)` now holds RIGOROUSLY on the tested box:** every ratio is a certified O(1)
  (band `[0.51, 3.24]`; the small-diam `[0..6]` and the doubling control sit outside `[1,1.7]`,
  but the *structured* diam-90/140/199 clusters give `1.36 / 1.55 / 1.09`, certifying exactly
  THM-729's empirical `[1.0,1.7]`).
- **Grid margin is thin.** The hardest cluster's min arc width is `4.8e-6`, only ~7x above
  klein's grid resolution `1/1.5e6 = 6.7e-7`. The grid was adequate here but a slightly more
  clustered config would fall below it and be silently mis-sampled; the exact-endpoint method
  removes that risk entirely. (THM-729's own `delta=0.012 -> 0.0000` warning is the same point.)

## Honest scope

This certifies the FINITE box and hands THM-729 a Lean-ready certified-`Q_s` primitive
(exact endpoints = rationals; tail bound = rational; head = cyclotomic/iv-certifiable). It does
**not** prove THM-729's open piece, the *uniform* soft bound `Q_s = o(r^2)` over all peels --
that is an analytic theorem (THM-729's width-weighted Montgomery-Vaughan 2nd moment), out of
reach of an exact-arithmetic certificate. The harvest's value is: (i) the empirical `O(diam)`
evidence is now rigorous where it is checked; (ii) a grid-free method robust to clustering;
(iii) a reusable float-free `Q_s` engine for the structured-check region `26<=d<=D_0`.

## Next (cheap, natural)

Certify the THM-729 DIAGONAL `diag = sum_i 2 pi^2 {w w_i}(1-{w w_i})` EXACTLY (it is a finite
sum of `{w w_i}(1-{w w_i})`, exact rationals) and subtract from the certified full `Q_s` to turn
the "off-diagonal cancels" claim from empirical into a certified `|off-diag| <= (bound)` on each
cluster -- a rigorous decomposition, and the first step toward a certified per-cluster route to
the `o(r^2)` cancellation.

Artifacts: `04-computation/lrc14_second_moment_certified_opus_S4.py`,
`05-knowledge/results/lrc14_second_moment_certified_opus_S4.out`.
