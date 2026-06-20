# Decompose to converge: the far-element peel

**Source:** mac-mini-2026-06-20-S1. Dispatch: understand "the one open constant" of the
LRC(14) endgame and push to complete the proof; integrate concurrent agent work. Canon:
THM-546. Built on kind-pasteur's HYP-2642 (far-element recursion) and S19 (finite half proved),
codex's HYP-2671/2673 (the empirical "constant stack"), and HYP-2643 (the divergence finding).

## The constant that couldn't be bounded

The LRC(14) endgame had been compressed to a single deviation: peeling the largest runner `w`
from a cluster `E`, how far does the actual seven-sector coverage `p0(E'∪{w})` sit above its
decorrelated limit `Φ(E')`? Call it `Δ_w`. The team needed `Δ_w` bounded — and ran into a
wall (HYP-2643, verified): the natural absolute bound, summing the magnitudes of all the
relation-lattice Fourier modes, **diverges**. For every cluster, wide ones included. The
"constant" `sup w|Δ_w|` is in fact *unbounded* (`Ω(spread)`, kps-S19). It looked like the
endgame's last number had no finite value.

## The number was finite; the *sum* was the wrong object

The divergence was an artefact of summing the whole lattice at once. `Δ_w` is a deviation
caused by **one** runner, `w`. So decompose by what `w` actually does: a single runner sits in
one sector per phase, so it can fill *at most one* missed sector. That collapses `Δ_w` to a
sum of six one-dimensional discrepancies — for each sector `j`, how much the indicator
`1_{sector j}(frac(wx))` over- or under-counts against its mean `1/7`, integrated over the
*fixed* set `B_j = {x : E' misses exactly sector j}`:
```
Δ_w = Σ_{j=1}^{6} [ meas{B_j ∩ (wx ∈ sector j)} − (1/7) meas(B_j) ]
    = Σ_j Σ_{n≠0} ŝ_j(n) · 1̂_{B_j}(−nw).
```
Now everything converges. The sector coefficient `ŝ_j(n) = |sin(πn/7)|/(π|n|)` decays like
`1/n` and *vanishes at every multiple of 7* — the apex prime, the same vanishing that has run
through this problem since the singular series. The set `B_j` is fixed and of bounded
variation, so `1̂_{B_j}` also decays like `1/n`. And the far runner appears only through the
single frequency `nw`, which gives the `1/w`. Term by term:
```
|Δ_w| ≤ κ · V(E') / (π² w),   κ = 2 Σ_{n≥1, 7∤n} |sin(πn/7)|/n² = 1.857…
```
`V(E')` is the arc-complexity of "misses exactly one sector." This is rigorous, it converges,
and it reproduces the `Ω(spread)` growth correctly (through `V ~ spread`) rather than blowing
up. The unbounded "constant" was the wrong functional; the bounded one is the per-element
discrepancy.

## The lesson

The whole-lattice envelope diverges because it adds up the interference of *all* the runners
with all the others — a `k`-dimensional sum with no decay budget. But the quantity that needs
bounding is one runner's effect, and a single runner is a one-dimensional object. The move was
to refuse the global sum and peel one element, turning a divergent `(k−1)`-fold lattice
envelope into a convergent 1-D BV discrepancy with an explicit `1/w` rate. The same shape
recurs across this project: the obstruction is real but stated about the wrong object — the
width not the measure, the speed scale not the invariant, the gap not the sector, and now the
*global lattice* not the *single peel*. Each time, choosing the smallest faithful object makes
the bound finite.

## Where it leaves the proof

After kind-pasteur's S19 (the finite half `max(E)≤14` proved, the wide half computationally
comfortable with margin ≥ 0.22), THM-546 supplies the rigorous analytic backbone for the
**gapped** far-element regime: once `w` exceeds an explicit multiple of `max(E')`, `|Δ_w|` is
below the recursion margin and the cluster certifies. The ungapped-wide regime closes by
scale-invariance and the Freiman dimension penalty, and the tight margin lives only in the
done finite check. The endgame's last analytic unknown, `Δ_w`, now has a closed-form bound
with an explicit constant `κ = 1.857`. LRC(14) is not proved — but its one open constant is no
longer open; it is `1.857`, and the question that remains is bookkeeping: iterating the peel
down to the bounded core and pairing the per-element bounds against the recursion margins.
