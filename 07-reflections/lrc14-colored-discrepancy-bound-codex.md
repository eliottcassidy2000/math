# LRC14 Colored Discrepancy Boundary

Codex 2026-06-18.

The phase-color reservoir gave the right finite object.  This session asked
whether the last CRT placement step can be reduced to a plain interval-grid
lemma.

It can, at least in a crude but rigorous form.

For each color `b`, the admissible set `G_P cap C_b(E)` is a union of
intervals, and the grid `x=(b+14t)/(14V)` has mesh `1/V`.  If that layer has
measure `L` and `m` components, then it contains at least `V*L-m` grid points
in the open intervals.  Summing over the fourteen colors gives

`actual_count(q=14V) >= V*Sigma(P,E) - K(P,E)`.

This is a real theorem-shaped step.  It certifies any fixed shape once
`V>K/Sigma`.

## The Good News

For the named hard shapes, the cutoff is small:

- via-zero k7: `V>=91`;
- via-zero k9: `V>=118`;
- broad `1/90`: `V>=125`;
- near-via minimum: `V>=130`;
- quarter minimum: `V>=135`.

Across the structured bank, the worst cutoff was `244`.  So for that whole
bounded-spread bank the finite q=`14V` placement problem is now formally a
finite check below `244`.

That feels like a serious compression.  We are no longer waving at
equidistribution; we have an explicit integer inequality with an explicit
cutoff.

## The Catch

Raw component count is the wrong complexity measure for large-spread tails.
It can explode into the thousands even when exact witnesses are plentiful.
The worst random tail had `K=8226`, making the crude bound useless near the
natural scale `V=maxE+14`.

But the exact grid count did not behave like the component bound.  In `160`
random covering tail rows there were no zero witnesses, the worst ratio
`actual/(V*Sigma)` was about `0.63`, and the largest positive deficit
`V*Sigma-actual` was only about `71.9`.

The stored run also tested three deliberately simple positive-deficit bounds:
`100`, `7*k*cGP+1`, and `14*k*cGP+1`, where `cGP` is the number of components
of the small-speed safe set `G_P`.  All three survived the `160`-row exact
tail census.  A separate larger scratch sample broke the too-small
`14*cGP+1` guess but still left `7*k*cGP+1` looking plausible.  That suggests
the right complexity is "small comb times cluster size", not the thousands of
micro-components created by a large-spread cluster.

So the next theorem should not be "bound the number of continuous components."
It should use the arithmetic of the grid itself.  The endpoints are produced
by the same speeds that define the progression, and the grid sees that
structure.  Continuous Koksma throws away too much.

## Working Proof Shape

The proof route I would push next:

1. Prove the uniform `Sigma` floor from HYP-2593.
2. Prove a colored residue discrepancy bound of the form
   `actual_count >= V*Sigma - C_14`, or maybe `-C(P,k)`.
3. If `C_14 < c0*V` above some explicit `V0`, close large `V`.
4. Exhaust the finite `V<V0` residue-covering cases.

This is closely aligned with the B(k)/relation-lattice agents, but it is more
discrete: the desired estimate should be a modular endpoint cancellation
lemma, not just a continuous interval-component discrepancy theorem.

## Tournament Note

The boundary-efficiency tournament changed the color leadership from the pure
mass story.  HYP-2593 had colors `1` and `13` leading.  Here colors `3` and
`11` are best by mass per component.  That is a useful warning: the quotient
we choose really matters.  Pure reservoir mass and boundary-efficient reservoir
mass preserve different pieces of the LRC predicate.

LRC(14) is still open, but the finite-placement layer is now much less foggy.
