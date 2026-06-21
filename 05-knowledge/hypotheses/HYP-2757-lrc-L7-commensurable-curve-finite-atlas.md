---
id: HYP-2757
title: LRC(14) L7 resonant tail = a FINITE atlas of 1-D torus-curve coverages (commensurable far pairs trace closed geodesics)
status: ATLAS RUN -- finite (24 ratios q<=8) & SAFE (max p0_curve 0.247/0.401 k=9/10, 0 violations, margin>=0.20); non-resonant tail + rigor remain
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

## UNIFICATION with codex's cone (kps-S24): the curve atlas IS the |R|<=2 cone
In the curve limit the far PAIR (2 elements) sweeps 2 sectors at a time, so it can only cover the missed-set S
if |S|<=2. Hence
> p0_curve(p,q,B) = sum_{|S|<=2} meas_x(B misses exactly S) * P_curve(S),  P_curve(S)=meas_y{ (sector(qy),sector(py)) covers S }.
So the resonant L7 atlas is EXACTLY codex's small-|R| context cone (HYP-2697/2698), with |R|<=2 (= the hard
regime codex/mac-mini localized), computed through the (q,p)-torus curve: P_curve(S) is a 1-D survival width
(connecting to mac-mini's W_a, HYP-2745). The cone's "context-weighted" structure = meas_x(B misses exactly S)
as the weight, P_curve(S) as the capacity. The atlas (24 ratios, 0 violations, margin>=0.20) is the cone's
top-of-lattice check for the far pair. So HYP-2757 (curve) + HYP-2697/2698 (cone) + HYP-2700 (Z/7-coloring) are
ONE object: the |R|<=2 coverage of the far pair, finite over the commensurable atlas + decorrelated otherwise.

## VALIDATION (kps-S24): the curve limit is c-independent and < cap
p0(B u {cq,cp}), B=[0,2,..,12], as c grows (c=7,13,23,41,71): CONVERGES to a c-INDEPENDENT limit ~0.24-0.25 for
ratios 2/1, 7/4, 4/3, 5/3 (3/2 noisier, settling ~0.24-0.26) -- all < cap_9=0.494 (margin ~0.25). Small-c values
(0.36-0.40 at c=7) are the finite-c / near-merge regime (the comb correction, decaying). So the resonant-ratio
curve coverages are FINITE and SAFE. Generic (non-resonant) ratios decorrelate LOWER (full 2-torus, ~0.05) at
large c. So L7's worst LARGE-scale cases are the commensurable (curve) ones, all comfortable.

## ATLAS RESULT (kps-S24, lrc_L7_resonant_atlas_kps.out): FINITE & SAFE
24 commensurable ratios p/q in (1,2.15] with q<=8; curve coverage p0(B u {701q,701p}) over 60+ bounded bases:
- k=9: 1488 points, **max p0_curve = 0.2468** (cap 0.4943, margin 0.2475), **0 violations**, worst 5/3.
- k=10: 1488 points, **max p0_curve = 0.4009** (cap 0.6044, margin 0.2035), **0 violations**, worst 2/1.
So the resonant atlas is FINITE and SAFE. For q>8 the ratio is "more irrational" => the 2D discrepancy is
smaller => p0 decorrelates LOWER (toward the full 2-torus value ~0.05). So the WORST L7 cases are the low-q
resonances (in the atlas), all < cap. L7 closure = [finite atlas q<=Q0~8, SAFE, VERIFIED] + [q>Q0 decorrelated,
2D ET-Koksma tail bound] -- mirroring the single-far [finite window + comb bound]. 

## Tests / next
- Compute `p0_curve(p,q,B)` for all `p/q in (1,2.15)`, `q<=8`, and bounded bases B; verify `<= cap_k` with
  margin (the resonant atlas). Connect to mac-mini's `measS7 = sum_a W_a` (HYP-2745): on the curve each W_a
  becomes a 1-D survival width.
- Bound the finite-c correction by THM-546 and the decoupling of B from the curve.
- The non-resonant 2D ET-Koksma constant on the bounded ratio window.
-> OPEN-Q-108, HYP-2708 (codex r=2), HYP-2750 (mac-mini L7-tail), HYP-2745 (W_a).
