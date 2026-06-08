# The triangle recursion is one recursion at three temperatures; the Pfaffian is the Rosetta object (S715)

The prompt asked for two things: use the Pfaffian to translate between topology, geometry, graphs,
tournaments, and algebras; and find recursive truths in the seven-piece tiling decomposition, aiming at
a more efficient handle on the maximum Hamiltonian-path count A038375. Both come out of one picture.

The seven-piece decomposition is, geometrically, the Mobius reconstruction of the staircase triangle by
its three corner subtriangles. Three copies of the side-`(n-1)` triangle at the three corners cover the
side-`n` triangle; they overlap pairwise in three side-`(n-2)` triangles and all three in one side-`(n-3)`
triangle at the center. Inclusion-exclusion: `3*Tri(n-1) - 3*Tri(n-2) + Tri(n-3) = Tri(n)`, an identity
on triangular numbers that holds exactly, and the coefficient vector `(1,-3,3,-1)` is `(x-1)^3`, the third
finite difference. So the user's `A,B,C` (corners, plus), `D,E,F` (overlaps, minus), `G` (center, plus)
is not an analogy — it is the partition of unity of the triangle, and every cell ends up counted exactly
once: a corner cell in one big triangle, an edge cell in two-minus-one, an interior cell in
three-minus-three-plus-one.

That fixes the first recursive truth precisely. Any tournament invariant that is *additive over the
cells* of the triangle — a sum of local contributions — inherits the IE and therefore satisfies
`F(n) - 3F(n-1) + 3F(n-2) - F(n-3) = 0`. The characteristic polynomial is `(x-1)^3`, whose only sequence
solutions are quadratic polynomials in `n`. The number of tiles `C(n-1,2)`, the number of arcs, the
number of vertices — all have vanishing third difference, all quadratic. The geometry of the triangle
*forces* additive invariants to grow quadratically. This is the cold end of the story.

The Hamiltonian-path count is at the hot end, and it does not obey `(x-1)^3` because it is not additive
over cells. `H` is the independence polynomial of the conflict graph evaluated at two — a *product*, not
a sum — and the third difference of A038375 is conspicuously nonzero. But the triangle does not vanish
from the hot end; it leaves a fingerprint. The base-path staircase family satisfies an order-three
recurrence `H(k) = 3H(k-1) + H(k-2) + H(k-3)`, and the leading coefficient is again three — the three
corners — while the sub-leading terms, which the additive IE carried with alternating signs `-3, +1`, now
carry positive signs `+1, +1`. Passing from additive geometry to multiplicative combinatorics turns
subtraction into addition. The backbone (the three) is invariant; the corrections flip sign. The
characteristic polynomial goes from `(x-1)^3` (triple root at one, polynomial growth) to `x^3-3x^2-x-1`
(dominant root `~3.383`, exponential growth). Same skeleton, different temperature.

Between the two sits the Pfaffian, and this is where the user's instinct was sharpest. They guessed the
size-`(n-2)` pieces should be *subtracted*. That is exactly the Pfaffian. Its minor expansion
`Pf(S) = sum_j (-1)^{j-1} S_{1j} Pf(S_{1,j deleted})` runs over `(n-2)`-vertex deleted subtournaments,
entering with alternating sign — verified on every random tournament I tried. The Pfaffian is the
`n -> n-2` (Mode-B) step of the staircase made into an exact *signed* recursion, and its sign rule is
neither the cold IE's uniform alternation nor the hot count's uniform positivity but permutation parity,
the in-between. So the three step-sizes and three sign-regimes line up: additive `n->n-1` with
`(x-1)^3`, Pfaffian `n->n-2` with parity signs, multiplicative `n->n-2` with positive signs. One
recursion, continued in a sign/temperature parameter, with the leading three fixed throughout. This is
the same additive-versus-multiplicative axis that S712 found splitting `21=3*7` (a clean product) from
`22=2*11` (forced additive), and that S714 found as flat-versus-peaked autocorrelation. The triangle
recursion is that axis wearing its geometric clothes.

Why the Pfaffian translates across all five domains is now clear: it is the one object that lives in each
of them at once. In topology it is the Euler class — Chern-Gauss-Bonnet writes the Euler characteristic
as the integral of the Pfaffian of the curvature, and THM-120 has `|Pf|` separating the circle-phase from
the sphere-phase of the tournament's homology. In geometry it is `sqrt(det(MM*))`, the oriented volume,
the square root of the autocorrelation determinant from last session. In graph theory it is the signed
count of perfect matchings (FKT), and mod two it is the matching count of the complete graph, `(n-1)!!`,
forced odd. In tournament theory it is `Pf(S_T)`, the even-side parity invariant, with `det(I+2A)=Pf^2`.
In algebra it is the Berezin integral of `exp` of the skew form, the top of the exterior algebra, and the
doubling tower of the reals, complexes, quaternions, octonions, sedenions is the Cayley-Dickson shadow of
the same even/odd seam the Pfaffian names. Euler characteristic equals curvature integral equals signed
matchings equals skew-tournament parity equals Berezin top form: five sentences, one object. The Pfaffian
is the Rosetta stone because it is the simultaneous value of a single construction read in five languages.

As for computing A038375 more efficiently: the honest finding is that the exact maximizers are
circulant/Paley-like, not the base-path family (whose `3.383^k` falls far short — seventeen against
forty-five already at six vertices), so the clean order-three recurrence is a lower-bound handle, not the
answer. The true edge over a random tournament, `a(n) 2^{n-1} / n!`, grows only slowly (around two to
three through `n=13`, parity-oscillating, consistent with a `sqrt(n)`-type factor on top of
`n!/2^{n-1}`). The recursive understanding the triangle offers is structural rather than a closed
formula: the maximum lives at the hot, multiplicative end where the corrections are positive and the
growth is governed by a dominant root, and the search for a recursive *construction* of the maximizers is
the search for a circulant family realizing that hot recurrence — the natural next target.
