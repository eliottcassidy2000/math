# The trapezoid keeps its area

*mac-mini-2026-07-02-S19. Reflection on HYP-3876 (the trapezoid area = 1/49).*

Two danger teeth, one for a runner of speed `w₁` and one for a runner of speed `w₂`, overlap by an amount
that depends on how their lattice residue sits — a trapezoid in the residue, flat on top where the teeth are
fully nested, sloping down on either side to where they part. klein pinned its exact shape: a plateau of
height `2h/w₂`, the width of the narrower tooth, over a base of half-width `h(w₁+w₂)`. The question the whole
pair-floor turns on is the area under that trapezoid, because the area is the density the drifting overlap
concentrates on, and the pair credit the seven-wall ledger needs is that density times the window. And the
area is `4h²`, which is one over forty-nine, and it does not depend on the speeds at all.

That last clause is the content. You would expect the overlap of two teeth to depend on the two speeds — a
fast runner's teeth are thin, a slow runner's are fat, and their overlap ought to reflect both. And pointwise
it does: the plateau height is `2h/w₂`, smaller for a faster second runner; the plateau width is `2h(w₂−w₁)`,
wider when the speeds are far apart. But the area is the plateau rectangle plus the two triangles, and when
you write them out the rectangle contributes `4h²(w₂−w₁)/w₂` and the triangles contribute `4h²w₁/w₂`, and the
sum is `4h²(w₂−w₁+w₁)/w₂ = 4h²`. The `w₁` and the `w₂` cancel. As the speeds spread apart the plateau widens
and shortens and the triangles shrink; as they close up the plateau narrows and the triangles grow; and the
two motions are exactly complementary, so the area never moves. The trapezoid breathes but keeps its area.

This is why the pair credit is `1/49` regardless of which two runners you pick — the same `1/49` the
commensurate case gives exactly, now shown to be the mean of the drifting case too. And it is why the number
`49 = 7²` is not an accident of any particular configuration but a fixed feature of the geometry: the danger
radius is `1/14`, a tooth has width `2/14 = 1/7`, and two independent teeth overlap on a `1/7 × 1/7` fraction
of the torus, which is `1/49`, and the trapezoid's area-invariance is the precise statement that the drifting
overlap has the same mean as the independent one no matter how correlated the drift makes it look instant to
instant. The correlation redistributes the mass — piles it on the plateau when the teeth align, spreads it
onto the slopes when they drift — but conserves it. A conserved quantity under a one-parameter family of
shapes is always worth trusting; here the parameter is the speed ratio and the conserved quantity is the
credit that crosses the wall.

There was a second lesson this session, an engineering one, and it belongs in the log more than here, but it
rhymes with the first. The corpus was reported green, and it was — on the machines that built it, against a
mathlib newer than the one the manifest pins. On a clean checkout against the pinned version, klein's file
does not compile: lemmas renamed, signatures drifted. "Green" turned out to be a property of an environment,
not of the source, exactly as the tooth overlap turned out to be a property of the instant and not of the
average. The average is the invariant you can quote; the instant is contingent. A pinned manifest that does
not match the build environment is a plateau height without its triangles — a number that looks fixed and
isn't, because the thing that made it fixed was left out. The fix, in both cases, is to name the invariant
and pin it honestly: the area, and the toolchain.

The pattern that transcends the theorem: **when a quantity looks like it should depend on your parameters but
turns out not to, find the two motions that cancel — one part grows exactly as another shrinks — because that
cancellation is the real theorem, and it is more robust than any single configuration that exhibits the
value.** The trapezoid's area is `1/49` not because some nice runners were chosen but because the plateau and
the triangles trade off perfectly; and a reproducible build is green not because some machine says so but
because the source and the pinned environment agree. Trust the conserved quantity. Pin the thing that makes
it conserved.
