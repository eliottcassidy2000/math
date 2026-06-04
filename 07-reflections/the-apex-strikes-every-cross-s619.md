# The apex strikes every cross (S619)

I went looking for laminar flow and found why n=14 is hard, stated more cleanly than I had it before.

Put the runners on the space-time torus: each speed `v` is a straight flow line of slope `v`, the danger zone is the
thin horizontal strip around the origin, and a lonely time is a vertical line that misses every crossing of a flow
line with that strip. Now overlay an antipodal pair `{v, −v}`: two lines of opposite slope through the same strip,
an X. That X is the criss-cross the user named, and the reflection `t ↦ −t` is the symmetry that swaps its two arms.

For the n=14 wall — all six unit residues plus a multiple-of-14 apex — the picture became completely legible. At
every unit clock `j/14`, exactly one antipodal pair is tight: one arm needs you to nudge the time up to stay safe,
the other needs you to nudge down. On its own that is survivable; you pick a side. But the apex changes the game. A
multiple-of-14 speed sits *exactly* on the origin at every unit clock — it is the no-slip wall of the laminar flow,
velocity pinned to zero against the boundary. So it is dangerous at `ε = 0` and you are *forced* to nudge. And the
moment you nudge, in whichever direction, you walk into one arm of the cross. The apex strikes every cross. That is
the whole reason the unit-clock sieve cannot finish n=14, and now it is one sentence.

What I didn't expect is that the wall config is *lonely anyway* — `p₀ ≈ 0.012`, a genuine positive measure of lonely
times, not a measure-zero tightrope like the apex-free arithmetic progression. The lonely times just aren't on the
grid. They sit in the laminar shells between the apex's own fine arcs, offset by about `0.012` past each unit clock —
exactly far enough to step off the no-slip wall while the tight pair hasn't yet closed. The flow threads the shell.

And the symmetry survives into the answer. The lonely set is invariant under `t ↦ −t`, four intervals in two mirror
pairs. I checked both fixed points of that reflection and both are blocked — the origin by the apex, the half-clock
by the even runners, who pile onto the origin there. So the involution has nowhere to stand; it acts *freely* on the
lonely set. This is the free-orbit cascade from the rigidity work, made concrete: loneliness at n=14 is not a single
pinned witness but a free orbit of mirror-paired shells. The apex-free AP was the pinned case (one tight point at
`1/14`, the σ-fixed apex of the sheaf); the apex residual is the free case (a cascade of paired shells, no fixed
point). Same involution, opposite orbit type — the two halves of the perspective key, sitting in one config space.

The honest state: I did not prove LRC(14). But I moved the target. It used to be "handle the sixty-four classes" or
"the multiple-of-14 configs are lonely." Now it is one geometric inequality: the laminar shell between two crosses,
once the apex has forced you off the grid, is wider than the cross can close. Show that width stays positive and the
free orbit never collapses. I formalized the frame it sits in — the involution is even, both fixed points are
blocked, the cross swaps arms — so the remaining step is a single width bound on a shell, not a census of a wall.
