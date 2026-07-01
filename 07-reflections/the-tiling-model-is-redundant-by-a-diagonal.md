# The tiling model is redundant by a diagonal

*mac-mini-2026-07-01-S81. Reflection on HYP-3798.*

The owner remembered something small about four-vertex tournaments: all four isomorphism classes show up in
the tiling model when you flip the three tiles, but two flips already suffice if you fix the other arcs
right. It sounded like a curiosity about a base case. It is a doorway.

Make the question precise. Color the arc-hypercube — one bit per oriented arc, `C(n,2)` bits — by
isomorphism class. Ask for the smallest axis-aligned subcube, some `k` arcs left free and the rest nailed
down, whose `2^k` corners hit every class. Call the answer `kappa(n)`. The naive Hamiltonian-path tiling
uses `m = C(n-1,2)` free tiles; the information floor is `ceil(log2 A000568(n))`. The truth sits between,
and it is clean: `kappa(n) = 1 + C(n-2,2) = 1, 2, 4, 7, 11, 16` — the lazy-caterer numbers. The tiling model
over-parametrizes tournaments by exactly `n-3` coordinates.

And then the redundant coordinates turn out to have a shape. Fix the Hamiltonian path, and through `n=6` the
`n-3` tiles you can additionally freeze are precisely the arcs `(i, i+3)` — the *skip-two* diagonal, the line
one step inside the hypotenuse of the staircase triangle. The staircase decomposition the project has leaned
on for a hundred sessions — hypotenuse, source leg, sink leg — gains a new resident: a diagonal of tiles that
carries no isomorphism information once the rest of the board is set. Through `n=6` you can throw it away
entirely; at `n=7` the same diagonal gets you `454` of the `456` classes, two short — so the clean rule is
exact for the small cases and only *nearly* true after, the way the good patterns in this project tend to be.
The redundancy is real (`n-3` arcs, always); its *location* is a diagonal, exactly for a while and almost
forever.

There is a second invariant hiding in the same picture, and it is older than it looks. Ask not how few
arcs express *all* classes but how few flips reach *one* class from the transitive corner. That minimum, over
the best base path, is the number of backward tiles you cannot avoid — and it is exactly the minimum feedback
arc set, the classical measure of how cyclic a tournament is. The two questions, minimum-flips-to-one and
minimum-free-arcs-for-all, are the covering radius and the transversal dimension of the same colored cube:
two covering-code parameters of the `S_n`-orbit coloring of `Q_{C(n,2)}`. A tiling model and a coding
problem are the same object wearing different clothes. And a small fact drops out for free on the way: the
feedback-arc-minimizing order of a tournament is always a Hamiltonian path, because a backward arc between
neighbors could be swapped away — so the tiling model's "best base path" and the median order are one and the
same.

The pattern that transcends the theorem: **a natural parametrization is almost never tight, and its slack has
geometry.** The tiling model was built to be faithful — every tournament, one tiling — and it is, but it is
not efficient: it spends `C(n-1,2)` bits where `1 + C(n-2,2)` would do, and the `n-3` wasted bits are not
scattered, they are a diagonal. When a construction feels a little too large, the excess is worth locating,
because it tends to lie on a line you already know how to draw. Here the line is parallel to the hypotenuse,
one step in — the same hypotenuse that controls the `H = 1 + 2^d` formula, the fiber fraction, the blue
lines. The redundancy of the model and the geometry of the model are the same diagonal.
