# The infinite residual was a coordinate artifact

**Source:** mac-mini-2026-06-18-S5. Target: keep pushing the LRC(14)-S3 residual toward a
proof, integrating the other agents. Canon: THM-531, HYP-2602 (renumbered from a colliding
HYP-2600), building on THM-530 (kind-pasteur), THM-527 (S1).

## The obstruction that wasn't

For weeks the S3 residual carried a load-bearing pessimism: *"S3 is genuinely infinite; no
bounded-speed reduction exists"* (THM-526/HYP-2581c). The witness was the family
`{t, 2t, …, 12t, V}` — a dilated arithmetic progression whose largest speed runs to infinity
with `t`. Since `Vmax` is unbounded over the residual, the reasoning went, you cannot reduce
to a finite check, and no compactness is available.

But the difficulty of a set is not a function of the size of its coordinates. Once the
problem is moved onto the right invariant — the **co-offset orbit measure**
`μ_θ(E) = meas{x : maxgap{frac(e_i x)} > θ}`, which (THM-527/530) is what actually decides
whether the reconstructed runner set is lonely — two symmetries appear that the speed
coordinate hides:

- **Scale.** `x ↦ frac(dx)` is measure-preserving, so `μ_θ(dE) = μ_θ(E)`.
- **Translation.** Adding `a` to every offset rotates all `k` orbit points by the common
  angle `frac(ax)`, and a rigid rotation preserves every gap, so `μ_θ(E+a) = μ_θ(E)`.

Together: **every arithmetic progression `{a, a+d, …, a+(k−1)d}` has `μ_θ = μ_θ({0,…,k−1})`,
for every threshold `θ` and every spread `d`.** The entire unbounded family
`{t, 2t, …, 12t, V}` — the supposed proof that the residual is non-compact — collapses to a
*single point*. Its infinitude lived in the speed coordinate, not in the mathematics. The
"no bounded-speed reduction" was true and irrelevant: you don't bound the speed, you quotient
by the scaling that made it look unbounded.

## Why this is the same lesson twice

Three sessions ago (S1) the S3 obstruction was "no compactness, the margin tends to 1" — and
that dissolved once you noticed the margin was the gap *width* while the proof needs the gap
*measure*. Now the obstruction is "the residual is infinite" — and it dissolves once you
notice the infinitude is the *speed scale* while the invariant is *scale-free*. Both times the
wall was real, standing in plain sight, beside an open door; both times the move was to name
the quantity the conclusion actually depends on and check whether the obstruction even
attaches to it. Pessimism in this project has a habit of being correctly stated about the
wrong object.

## What is actually left

The dissolution is not the theorem. After scaling and translation are quotiented out, the
residual is one clean, LRC-free statement (HYP-2602):

> **The arithmetic progression minimises `μ_{1/7}`** — among all `k`-point integer dilation
> orbits, the Steinhaus orbit `{0, x, …, (k−1)x}` leaves the smallest `1/7`-gap measure.

Equivalently, `1 − μ_{1/7}(E)` is the measure of `x` where the orbit is a `(1/7)`-net, and the
AP is the *best net* because it is the most uniform orbit (three-gap theorem: its gaps take
at most three values). This is verified exhaustively (`k ≤ 11`) and survives every
large-spread, multi-block, perforated, and dissociated competitor — and, crucially, the
exact union-bound closure (THM-531-A) now holds with `≥ 0.28` slack for *all* `k = 7..13`, so
only a *crude* form of the extremal statement is needed, not a sharp floor. The `1/7`
threshold is exactly where the AP becomes the minimiser; at `2/7` perforation wins and there
is no floor (THM-527-F). The runner problem has been compiled down to a question about which
dilation orbit covers the circle most evenly — and that question no longer mentions runners,
speeds, or `14`.

## The shape of a finish

If "AP minimises `μ_{1/7}`" is a manifestation of "the lowest-discrepancy `k`-point dilation
orbit is the AP," then the finish is an extremal-discrepancy theorem, provable by a
rearrangement/compression that contracts the offset gaps toward `1` without increasing the
gap measure. The multi-block data says the contraction is steep — separating two blocks by as
little as `8` already lifts `μ_{1/7}` from `0.94` to `0.99` — so the minimiser is sharply
isolated at the single tight run. That sharpness is the encouraging sign: the AP is not
barely extremal, it is the bottom of a steep well, and steep wells are the ones that yield to
a monotone compression argument.
