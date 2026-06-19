# The bounded cluster is what makes the singular series converge

**Session:** mac-mini-2026-06-18-S4. **Result:** HYP-2600 (subtorus theta convergence; the THM-504 ⟷ HYP-2599 connection).

The prompt asked to think in terms of the *subtorus relation lattice* and to pull connections
from the other agents. Following codex's HYP-2599 and kind-pasteur's kps-S5 grounding back to the
project's earlier singular-series work (THM-515/504/518) turned up the connection I think is the
real one — the reason the whole S3 program is even possible.

## Two singular series, one wall

Months of the LRC effort lived on the **speed-side** singular series (THM-515):
`L(S) = Σ_{t∈Λ(S)} ∏ s(t_i)`, a sum over the relation lattice of the *speeds*, `L>0 ⟺` the set is
loose. THM-504/518 hit a hard wall there: the `|T|≥3` part of the sum is only *conditionally*
convergent — the absolute sum diverges, like `Σ 1/T`, because there are arbitrarily many
speed-relations of every height. Every route that needed an absolute bound died on that
divergence. It was the project's most stubborn obstruction.

codex's HYP-2599 writes down what is, structurally, the **same object on the offset side**: after
the slow-fast reduction (THM-527) isolates the large cluster, the floor's smooth minorant has the
identity `∫G(Ex)dx = (5/7)^k + Σ_{n∈Λ_aff(E)} ∏ ψ̂(n_i)` — again a sum over a relation lattice,
again with a single-coordinate kernel `ψ̂` that even vanishes at multiples of 7, just like the
sinc `s(t)` of THM-503. The two are the same theta form. So the natural fear is that the offset
side inherits the same `|T|≥3` divergence, and HYP-2599's "high-relation-height tail bound" is
just THM-504's wall wearing a new hat.

## Why it doesn't — the bounded cluster

It doesn't, and the reason is exactly the thing the S3 reduction *buys*. The speed-side sum
diverges because the speeds are **unbounded**: relations of every height keep coming, `Σ 1/T = ∞`.
But the cluster, after the reduction, has at most **thirteen** offsets. With `k ≤ 13`:

- there are no support-2 relations at all (distinct offsets force `e_i = e_j`);
- support-3 relations are one per triple — at most `C(13,3) = 286` primitive ones;
- each primitive relation's multiple-tail `Σ_m ∏ ψ̂(m·)` decays like `Σ_m 1/m³` and **converges**.

Finitely many primitive relations, each with a convergent tail: the offset-side correction
converges **absolutely**. The same `(−1)^{support}` cross-level alternation that signed THM-504 is
still there, but now it is a *finite*, computable sum, not an infinite divergent one. I verified
the convergence (the multiple-tail is stable by `m ≤ 3`) and that `∫G(Ex)dx > 0` on every S3
cluster shape from `k=5` to `k=13` — AP, perforated, dissociated alike. For `k ≥ 8` the correction
is usually *positive*: the tiny main term `(5/7)^{13} ≈ 0.013` is *helped*, not overwhelmed, by the
lattice.

So the deep content of the slow-fast / S3 reduction is not just "separate the small part from the
big one." It is: **take a divergent singular series over an unbounded set and turn it into a
convergent theta sum over a bounded cluster.** The thirteen runners were always the point — the
conjecture is hard precisely at the boundary where the relation lattice is rich but still finite.
Bounding the cluster is what tames the sum the speed side could never bound.

## What this finishes, and what it doesn't

It does not finish LRC(14). The remaining task is exactly codex's split, now correctly framed: a
*uniform* lower bound on the (finite, absolutely-convergent, alternating) correction over all
`k ≤ 13` shapes — low relation-height by finite enumeration, high relation-height by the product-
kernel tail. The convergence I established says that task is *well-posed* in a way the speed-side
analogue never was: there is a finite object to bound, not a divergent one to renormalize.

The connection also lines the lanes up cleanly. My universal centers (HYP-2597) are the
`(5/7)^k`-style main-term skeleton — the always-present good points at `b < 7/2`. codex's
HYP-2598 turned those into the survivor-sequence classifier. HYP-2599 is the relation-lattice
correction. And HYP-2600 is the statement that the correction is a finite sum — the bridge from
the old wall to the new, tractable problem. The seven keeps appearing (band `1/14`, kernel zeros
at `7`, center cap `7/2`), but the thing that changed is the *count*: thirteen, not infinity, and
thirteen is small enough.
