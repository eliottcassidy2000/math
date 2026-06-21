---
id: HYP-2757
title: LRC(14) L7 resonant tail = a FINITE atlas of 1-D torus-curve coverages (commensurable far pairs trace closed geodesics)
status: VALIDATED (curve limit converges, c-independent, <cap); finite-atlas check + non-resonant bound remain
source: kind-pasteur-2026-06-21-S24
depends_on:
  - OPEN-Q-108
related:
  - HYP-2708   # codex two-far live-depth (the r=2 instance)
  - HYP-2750   # mac-mini L7-tail
  - HYP-2745   # mac-mini G_P/Markoff per-cell width closed form (measS7=sum W_a)
  - THM-546    # single-far comb bound (r=1)
  - THM-557    # single-block
---

# HYP-2757 - L7 resonant tail = finite atlas of 1-D torus-curve coverages

## The setup (where L7 lives)
LRC(14) is reduced (kps-S23 ledger) to L7: the BALANCED multi-cluster wide bound, narrowed to the
**bounded ratio window** `f_2/f_1 in (1, 2.15)` for r=2 (above 2.15 the geometric comb chain L6 closes;
near 1 the elements merge into one cluster, THM-557). L7 = [resonant atlas] + [non-resonant 2D bound].

## The claim (the resonant atlas is FINITE and 1-DIMENSIONAL)
For a COMMENSURABLE far pair `f_2/f_1 = p/q` (small denom q, gcd(p,q)=1) in the window, write
`{f_1,f_2} = c*{q,p}`, `c=f_1/q`. The joint carrier orbit is
```text
(frac(f_1 x), frac(f_2 x)) = (frac(q*(cx)), frac(p*(cx))) = (frac(q y), frac(p y)),  y := c x.
```
As `c -> infinity`, `y` equidistributes on [0,1), so the joint orbit traces the **closed 1-D CURVE**
(a (q,p)-geodesic) `{(frac(qy), frac(py)) : y in [0,1)}` on the 2-torus -- it does NOT fill the torus.
Therefore
> **p0(B u c*{q,p}) -> p0_curve(p,q,B) := meas_y{ B-colors u {sector(qy), sector(py)} hit all 6 inner },**
a **1-DIMENSIONAL** integral over `y`, INDEPENDENT of `c` (for large c, by Weyl on `y`). The base `B` is
bounded (<=14); its colors at slow time `x=y/c` are frozen on each `y`-cell... [precise coupling: B-colors
vary on scale 1, the curve on scale 1/c; the limit decouples B's color-distribution from the (q,p)-curve].

CONSEQUENCE: the resonant part of L7 is a **FINITE ATLAS** -- one 1-D curve-coverage `p0_curve(p,q,B)` per
ratio `p/q in (1, 2.15)` with `q <= Q0` (finitely many: 5/4,4/3,7/5,3/2,8/5,5/3,7/4,9/5,2/1,...) and per
bounded base B (finitely many). Check `p0_curve(p,q,B) + finite-c correction <= cap_k` for all of them. The
finite-c correction is the single-far comb error (THM-546, PROVED). The NON-resonant ratios equidistribute
on the full 2-torus => p0 = the decorrelated carrier product (small) + a bounded-ratio 2D ET-Koksma error.

## Why this should close L7
The curve-reduction turns the "infinite two-free-scale" obstruction into: (i) a FINITE list of 1-D
coverage integrals (the geodesics, indexed by the small-denom ratios in the bounded window), each computable
exactly, plus (ii) the generic 2-torus case where the ratio is bounded so the relevant ET-Koksma frequencies
are bounded and a lossy constant suffices (margin 0.21). This is the SAME architecture as the single-far
closure (finite window + comb bound), lifted one dimension via the geodesic decomposition.

## VALIDATION (kps-S24): the curve limit is c-independent and < cap
p0(B u {cq,cp}), B=[0,2,..,12], as c grows (c=7,13,23,41,71): CONVERGES to a c-INDEPENDENT limit ~0.24-0.25 for
ratios 2/1, 7/4, 4/3, 5/3 (3/2 noisier, settling ~0.24-0.26) -- all < cap_9=0.494 (margin ~0.25). Small-c values
(0.36-0.40 at c=7) are the finite-c / near-merge regime (the comb correction, decaying). So the resonant-ratio
curve coverages are FINITE and SAFE. Generic (non-resonant) ratios decorrelate LOWER (full 2-torus, ~0.05) at
large c. So L7's worst LARGE-scale cases are the commensurable (curve) ones, all comfortable.

## Tests / next
- Compute `p0_curve(p,q,B)` for all `p/q in (1,2.15)`, `q<=8`, and bounded bases B; verify `<= cap_k` with
  margin (the resonant atlas). Connect to mac-mini's `measS7 = sum_a W_a` (HYP-2745): on the curve each W_a
  becomes a 1-D survival width.
- Bound the finite-c correction by THM-546 and the decoupling of B from the curve.
- The non-resonant 2D ET-Koksma constant on the bounded ratio window.
-> OPEN-Q-108, HYP-2708 (codex r=2), HYP-2750 (mac-mini L7-tail), HYP-2745 (W_a).
