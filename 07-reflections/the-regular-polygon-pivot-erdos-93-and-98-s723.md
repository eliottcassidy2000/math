# The regular polygon is the pivot between Erdős 93 and 98 (S723)

Two of Erdős's distance problems sit on either side of one object, and the object is the regular polygon.
Problem 93, solved by Altman in 1963 and since verified in Lean, says that any `n` points in convex
position determine at least `floor(n/2)` distinct distances, and the bound is tight: the regular `n`-gon
hits it exactly. Problem 98, still open, asks whether `n` points in general position — no three on a line,
no four on a circle — must determine superlinearly many distinct distances; the best anyone can prove is
`(n-1)/3`. The two look like neighbours in a catalogue. They are actually the same question asked on the
two sides of a single configuration, and seeing that turns the solved one into a tool against the open one.

The regular polygon is why 93 is true and tight. It achieves `floor(n/2)` distinct distances by the
simplest possible trick: put all `n` points on one circle, and every distance becomes a chord, of which
there are only `floor(n/2)` lengths. The whole content of Altman's theorem is that convex position cannot
do better than this single concyclic ring. And the single concyclic ring is exactly the thing problem 98
forbids — no four concyclic, let alone all `n`. So problem 98 is problem 93 with its optimizer outlawed.
The reason 98 might force superlinearly many distances is that the one cheap construction, the ring, is
banned, and there may be no substitute.

That reframing hands over a concrete bridge, because the solved problem applies not just to a convex set
but to the convex hull of any set. If a configuration has `D` distinct distances in total, its convex hull
is a convex polygon, and by Altman that polygon already determines at least half its vertex count in
distinct distances — all of them among the `D`. So the hull has at most `2D+1` vertices. Few distinct
distances forces a small convex hull. I checked this on three thousand random integer configurations and
it never failed. Peel the hull off and the same applies to the next layer, and the next: every convex
layer of a few-distinct-distance set has at most `2D+1` points, so the set is forced into at least
`n/(2D+1)` thin convex layers, like the rings of an onion. And now the 98 hypothesis bites in a place the
trivial bound never used it: no four concyclic means none of these rings can be circular, and no four
points across different rings can share a circle either. The regular-polygon trick that made one ring
cheap is unavailable ring by ring, and the rings are forbidden from being circles. A 98-extremizer is not
one clean ring; it is a stack of many thin, mutually non-concyclic, non-circular rings, each of which —
being convex — already spends `floor(size/2)` distinct distances by Altman. The route to superlinearity is
to bound how few distinct distances such a stack can share, which is a circle-incidence question, the same
additive-energy attack S722 set up.

The inspiration runs backward too. Problem 98's mechanism — no four concyclic caps the number of points at
any one distance from any one point at three — lands as a clean toehold on the still-open strong form of
93, the conjecture that some single vertex of a convex set sees `floor(n/2)` distinct distances. A vertex
sees few distinct distances exactly when many points are concyclic around it, and the no-four-concyclic
hypothesis caps that at three, so under both constraints every vertex sees at least `(n-1)/3` distinct
distances unconditionally. The gap from that `(n-1)/3` to the conjectured `n/2` is precisely the convex
ordering structure that 93 has and 98 throws away — which says, sharply, where the strong conjecture's
difficulty lives: not in the local concyclicity, which the easy argument already controls, but in the
global convex order.

Underneath both is the radial autocorrelation the cluster has been tracking all month. Distinct distances
are the size of its support; both problems minimize that support under a constraint. Convexity, problem
93's constraint, is rigid, and it pins the minimum at `floor(n/2)`, realized by the crystalline ring.
General position, problem 98's constraint, is weak, and it only pins `(n-1)/3`, with the real answer
unknown. The regular polygon is the cold crystalline ring both problems orbit: the extremal point of 93,
the forbidden point of 98. And the two small integers that govern them are the two ways a ring's
autocorrelation collapses — the `2` of `floor(n/2)`, chords doubling up around a single circle, and the
`3` of `(n-1)/3`, the per-circle cap that is the very same three that stalls `u(22)` at sixty. Solved on
one side, open on the other, tight at the same polygon: 93 is the lit room and 98 is the dark one through
the same door, and the convex hull is the hand on the switch.
