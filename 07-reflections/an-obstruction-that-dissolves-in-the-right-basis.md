# An obstruction that dissolves in the right basis

*klein-2026-07-01-S78. A reflection on HYP-3811 — chasing the SC-covering obstruction into the right
coordinates, and watching it disappear.*

Last session I found that the self-complementary classes carry the flip-rank covering-excess: they are the
odd boundary of a T-join, they cluster oddly in the blue subspace, and covering them with a low-dimensional
subcube costs extra dimensions. I called it an obstruction, and I was careful to say the parity only
co-located the hardness rather than proving it. This session the owner asked for the parity lower bound, and
for whether the half-tiling folds once more. Both questions had the same answer, and it corrected me: the
obstruction is not real. It is a coordinate artifact.

Here is what happened. The blue tilings are the fixed set of the grid reflection, and that reflection is
linear, so the blue tilings are a *linear subspace* `W`. The half-tiling model gives `W` its own
coordinates — one bit per reflection-orbit, the half-addresses. I recomputed the SC covering problem, but in
those coordinates instead of the full tile coordinates. In the full coordinates the SC classes need one or
two extra dimensions above the information floor. In the half-address coordinates they need exactly
`log2(#SC)` — the floor, no excess, at `n = 4, 5, 6`. The same classes, the same covering question, a
different basis, and the excess is gone. The T-join boundary is still odd, the parity is still there; but
the *hardness* it seemed to cause was the cost of looking at a linear subspace through axes that do not
respect it. Change to the axes that do, and the SC classes separate perfectly.

This is exactly the lesson the even-graph detour taught a few sessions ago — that the quotient, not the
cube, carries the content, and that a reframe is only trustworthy if it preserves the right structure — but
turned around into a gift. There it was a warning: a bijection of sets can destroy the shape. Here it is the
constructive version: if a hardness appears in one coordinate system, ask whether it survives the natural
one. The half-tiling model was not just a smaller enumeration trick; it is the coordinate system in which
the self-complementary classes are what they are — a vector space with a clean linear structure — and in
which their covering problem is trivial. The obstruction was a mirage cast by the wrong basis.

The folding question turned out to be the same story from the symmetry side. The tiling cube has two
commuting involutions: the grid reflection, which is the complement, and the flip of all tiles. Together
they are a Klein four-group, and the cube folds twice — full to half by the complement, half to quarter by
the flip. The blue lines, which I had studied as the top layer of the metagraph, are precisely the
reflection-fixed part of the quarter, the size-two orbits. So the half tiling does fold to a quarter, and
the object I called the blue lines is that quarter's fixed locus. The complement is not an extra symmetry
bolted onto the tiling model; it *is* one of the two axes of the fold, and the Hamiltonian-path count `H`
rides through it unchanged — `H(T) = H(T^op)`, odd by Rédei, one value per merged node. The two folds and
the one invariant `H` are the whole coarse structure.

The methodological point I want to keep: distinguish an obstruction from a coordinate. An obstruction
survives every change of basis; a coordinate artifact does not. When I measure a hardness — a covering
excess, a rank defect, a parity that seems to block — the honest next move is to find the object's natural
coordinates and remeasure. If the hardness persists, it is real and I can hunt its cause; if it dissolves,
the hardness was mine, not the object's, and the dissolution names the natural coordinates as a bonus. Here
the SC classes handed me both: the excess dissolved, and the coordinates it dissolved in are the half-tiling
addresses — the same coordinates in which the cube folds to its quarter. The obstruction was pointing, all
along, at where the object wanted to be described.
