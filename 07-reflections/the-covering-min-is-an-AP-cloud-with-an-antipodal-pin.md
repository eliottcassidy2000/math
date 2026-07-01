# The covering-min is an arithmetic progression with one antipodal pin

*klein-2026-07-01-S80. A reflection on HYP-3813 — putting the covering-min in phase coordinates and watching
five separate results collapse into one picture.*

For a long run of sessions the covering-min accumulated invariants, each true and each partial: a Chebyshev
alternation of length two on the slowest runner and the killer; a binding modulus that is the composite
`Phi6` and not its prime factors; a phase-residue `p(v) = nv mod Phi6` that governs coupling; a runner-cloud
tournament; the three-gap theorem lurking in the background. They read like a list. This session, asked to
integrate the tournament thinking with the runner, I did the simplest possible thing — drew the runners at
the binding time in the phase coordinates — and the list turned into a single object.

The object is almost embarrassingly concrete. Put the runners at `t*` in the coordinate `p(v) = nv mod
Phi6`. The small runners `1, ..., n-2` land at `n, 2n, ..., (n-2)n` — an arithmetic progression of step `n`.
The killer lands at `-n`. The observer sits at `0`, in the gap between the killer at `-n` and the first AP
point at `+n`, with clearance `n` on each side. That is the entire covering-min: an arithmetic progression
of step `n` with one antipodal pin, and an observer caught in the doubled gap. Everything I had been proving
separately is a feature of this one drawing. The Chebyshev two-point dual is the two cloud points that flank
the observer's gap, `+n` and `-n`. The three-gap theorem gives the gap sizes, and they are exactly `1, n,
2n` — the AP step, the sliver the killer cuts, and the observer's doubled gap. The killer being the
`iota`-antipode of the slowest runner is why the flanking is symmetric. And the value `n/Phi6` is no longer
a computation; it is the AP step divided by the modulus, read off the picture.

What made the picture appear was the coordinate change, and that is the through-line of this whole arc. Two
sessions ago the tournament SC-covering excess dissolved when I moved to the half-address coordinates — a
coordinate artifact. Last session the runner's obstruction refused to dissolve under the same kind of fold,
because it lives metric-irreducibly at the composite `Phi6`. This session the phase-residue coordinate did
not dissolve the obstruction — it cannot, the obstruction is real — but it made the obstruction *visible*:
the reason the observer's gap cannot shrink below `n` is that the covering condition forces the small speeds,
whose phases are the step-`n` AP, and an AP of step `n` simply cannot leave a gap smaller than `n`. The
metric-irreducibility of `Phi6` is, in these coordinates, the trivial fact that an arithmetic progression
tiles with its step. The right coordinates did not make the problem easy; they made the difficulty
*legible*, which is the next best thing and often the prerequisite for a proof.

The killer's double role also becomes plain in the picture, and it is a small lesson in construction design.
The clearance `n` is already there from the slowest runner alone — runner `1` at `+n` is the binder. The
killer does not create the clearance; it does two other jobs. It supplies the covering multiples of the two
moduli the small speeds miss (`n-1` and `n`), so the set is a covering set at all. And it sits at `-n`, the
antipode, converting a one-sided near-miss into a symmetric two-point alternation — turning a configuration
that merely achieves the bound into one that is a genuine Chebyshev extremizer, rigid and isolated. The
construction is not "small runners plus a big number chosen by search"; it is "an arithmetic progression
plus the one pin that both closes the covering and symmetrizes the gap," and that pin is forced to be the
antipode by the killer identity `n(n-1) = -1 mod Phi6`.

The keeper is about integration itself. When a problem has accumulated many invariants that feel unrelated,
the move is not another invariant but a coordinate system in which they might all be shadows of one thing.
The phase-residue coordinate was already in hand; I had used it locally for resonance and for the binding
modulus, but not stood inside it and looked at the whole configuration. Doing so cost one plot and returned
five results as one. Before hunting for the next lever, it is worth asking whether the levers already found
are the same lever seen from different sides — and if so, finding the vantage point from which they line up.
Here the vantage point is the phase cloud, and from it the covering-min is just an arithmetic progression
with an antipodal pin.
