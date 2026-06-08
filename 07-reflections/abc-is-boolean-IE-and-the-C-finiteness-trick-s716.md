# A+B+C-D-E-F+G is boolean inclusion-exclusion, and it hands you two computational tricks (S716)

The recurrence the user keeps pointing at is not seven arbitrary pieces. It is the seven nonempty subsets
of a three-element set: three singletons added, three pairs subtracted, one triple added. The signs are
`(-1)^{|S|+1}`, the Mobius function of the boolean lattice `B_3`, so `A+B+C-D-E-F+G` is precisely the
inclusion-exclusion count of the union of three things — here the three corner subtriangles of the
staircase. Once you see it as IE over `B_3`, two questions answer themselves: what generalizes it, and
what it is good for.

It generalizes by dimension. Inclusion-exclusion over the corners of a `d`-simplex has `C(d+1,k)` pieces
of size `n-k`, and the operator is `(x-1)^{d+1}`, the `(d+1)`-th finite difference. The interval (two
corners) gives `(x-1)^2` and linear cell counts; the triangle (three corners) gives `(x-1)^3` and
quadratic; the tetrahedron gives `(x-1)^4` and cubic. The tournament staircase is the two-simplex, which
is why the user's recurrence has order three and why everything additive on it grows quadratically. The
order of the recurrence counts the corners; the degree of growth is one less. That is the first relation,
and it is exact for `d` up to four (and clearly all `d`).

What it is good for splits along the additive/multiplicative line the project keeps meeting. An invariant
obeys the *exact* seven-term identity if and only if it is valuative — a measure on the tiles or arcs. Arc
counts do (I checked: a random eight-vertex tournament's arc total is reconstructed exactly from three
vertex-deleted corners minus three overlaps plus the center, `28 = 28`); the number of tiles does; any
linear functional of the tile-indicator does. For all of these the seven-term identity collapses to
`Delta^3 F = 0`, the invariant is a quadratic polynomial in `n`, and you need only three seed values to
get every term in constant time. That is the first computational trick, and it is the whole content of
"this invariant is additive": stop enumerating, fit a parabola.

The number of Hamiltonian paths is not additive — it is a product over the conflict graph, degree `n-1`
in the tiles — so it breaks the seven-term identity. But it breaks it in a structured way that hands over
the *second*, sharper trick. For any fixed recursive *family* of tournaments, `H` is C-finite: it
satisfies a linear recurrence, because the family is assembled by a transfer operation of bounded width
(this is THM-291's boundary — the overlap, the wiring, the apex — made into a finite-state machine). And
a C-finite sequence is the playground of the companion matrix. Compute a handful of terms by brute force,
solve a small Hankel system to *discover* the recurrence, and then exponentiate the companion matrix to
land on any term in logarithmic time. I ran this end to end on the base-path family: nine values from the
exponential dynamic program, the finder hands back the coefficients `(3,1,1)` — exactly THM-337 — and then
matrix powering produces `H` at `k = 100`, a fifty-three-digit number, instantly, while the direct
dynamic program is `O(2^{2k})` and dies around `k = 13`. The number of Hamiltonian paths is generically
expensive, but along any fixed recursive family it is `O(log n)`. That is the real computational payoff of
the user's decomposition: the staircase recursion is a transfer matrix in disguise, and transfer matrices
exponentiate.

The two tricks are the same fact at two temperatures, which is the through-line from last session. The
recurrence *order* is geometric — three, the number of corners, the boundary width — and it does not
change. The recurrence *roots* move with temperature: the additive end has `(x-1)^3`, a triple root at
one and polynomial growth; the multiplicative end has `x^3 - 3x^2 - x - 1`, a dominant root near
`3.383` and exponential growth. Same order, frozen by the geometry; different roots, set by whether you
are adding or multiplying. The trick in each regime is matched to the roots: a parabola fit when the
roots are all one, a companion-matrix power when they are not.

And the boundary of the method is exactly where it should be. The maximum number of Hamiltonian paths
over *all* tournaments, A038375, is not C-finite — no linear recurrence of order up to five fits its
thirteen known terms — because the optimum ranges over tournaments of unbounded effective width, so there
is no fixed transfer matrix to power. That is the structural reason THM-329 had to find new A038375 terms
by annealing rather than by recurrence: the trick is a gift of *fixed families*, and the optimum is not a
fixed family. The boolean inclusion-exclusion tells you precisely which problems are cheap (additive, or
fixed-width multiplicative) and which are genuinely hard (the free optimum), and it hands you the fast
algorithm in both of the cheap cases.
