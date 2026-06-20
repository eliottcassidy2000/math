# The true-wide branch decorrelates to exactly the boundary-collar plateau

*kps-2026-06-20-Sx-wf, while proving the HYP-2675 true-wide branch.*

## The coincidence that is not a coincidence

The LRC(14) sector crux splits into three regions:
- **finite half** (`max(E) ≤ 14`): proved exhaustively (kps-S19).
- **boundary collar** (`2nd-largest ≤ 14`, one far element): closed by THM-547, where
  the single far runner peels off and `p0(E) = Plat(E') + Δ_w ≤ Qb(k-1) + (6/49)V/w`.
- **true-wide** (`2nd-largest > 14`, ≥ 2 far elements): the last open region.

I attacked true-wide by a *different* mechanism — multi-cluster Weyl decorrelation. As the
inter-cluster gaps go to infinity, `p0(E)` converges to a decorrelated limit `p0_inf` in
which all clusters share the common dilation `x` and the inter-cluster offsets become
independent uniform rotations on the circle.

The worst case of this decorrelated limit, over **every** ≥2-cluster shape-multiset with
`k` points, turned out to be `[k-1 consecutive] + [1 singleton]`, with value EXACTLY

  `Qb(k-1) = 0.19660, 0.36210, 0.44789, 0.53125, 0.60224`  (k = 8..12)

— the **same** plateau `Qb(k-1)` that THM-547 uses for the boundary collar.

## Why it had to be this

Two ostensibly different physical pictures — *one* far runner slowly decorrelating
(THM-547) versus *several* far clusters decorrelating against each other (true-wide) —
land on the identical extremal value. The reason is that the decorrelated limit only ever
adds ONE free rotation phase to the optimization that matters: concentrating all the
"extra" points into a single tight consecutive cluster (which alone can cover the most
sectors) and leaving the rest as the thinnest possible far contribution (a singleton, one
random sector) is the configuration that maximizes "all six sectors covered." More
clusters means more points spread thinner across more independent phases, which *lowers*
the joint-coverage probability, not raises it. So the multi-cluster optimum collapses back
to the single-far peel.

This is the **same multi-scale residue collapse** the C(k) workflow predicted from the
other side: the spreading that destroys the `Δ_w` cancellation is exactly the spreading
that thins sector coverage. Here it shows up as: *the more genuinely wide the set, the more
its `p0` is pinned below the one-far plateau.* True-wide is never worse than the boundary
collar — it is the boundary collar's comfortable margin, inherited through decorrelation.

## The rigor lesson (a near-miss)

My first engine used an **independent-`x` convolution** of per-cluster coverage
distributions. It gave `17/343 = 0.0496` for two `consec_4` clusters — but the true
`M→∞` limit is `23/196 = 0.1173`. The independent-`x` model *under*-counts because it
forgets that all clusters see the *same* dilation `x`; only the offsets are independent.
An under-counting model is not an upper bound, so the whole worst-case search would have
been unsound. The fix — the **shared-`x`** engine (integrate over the common `x` and the
independent rotation phases) — both matches the limit and reproduces the THM-547 plateau,
which is the cross-check that confirms it is the right object.

The mathematics keeps insisting that the far structure of LRC(14) is governed by one
plateau constant `Qb(k-1)`, reachable from at least three independent directions
(single-far peel, multi-cluster decorrelation, and the finite-check argmax). When the same
rational constant arrives by three roads, it is the invariant, not the route.
