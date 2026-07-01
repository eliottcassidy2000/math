# A bijection of sets is not a bijection of structure

*klein-2026-07-01-S74. A reflection on HYP-3807 — testing my flip-rank findings against the
tournament<->even-graph bijection, and being wrong in a useful way.*

The owner asked me to challenge my assumptions using the natural bijection between tournaments and even
graphs. I had an assumption ready to challenge. Two sessions ago (S72) I found that tournament isomorphism
classes *pack* into the cube at the information floor — the rainbow number `R(n)` equals `floor(log2|G_n|)`
exactly — and I quietly believed this was a fact about the cube with its symmetry, something structural and
therefore portable. The even-graph bijection was the perfect place to test that belief, because even graphs
live on the *same cube*: the tiling-to-fundamental-cycle map is a genuine isomorphism `GF(2)^m ->` even
graphs, `m = C(n-1,2)`, the very cube whose points I had been calling tournaments. Same cube, different
labels. If the packing law were about the cube, it would survive relabeling.

It did not. The even-graph rainbow number is `1,1,2,3` where the tournament's is `1,2,3,5`, and at `n=6` the
even-graph rainbow *fails the floor* (`3 < 4`) that the tournament rainbow always meets. The covering side
is worse still: to hit all sixteen even-graph classes you must free nine of the ten coordinates, a
covering-excess of five, against the tournament's excess of one. The two invariants, computed on the
identical cube, disagree sharply. So the bijection is an isomorphism of the *set* — every even graph is
exactly one tile-vector — but not an isomorphism of the *structure* I cared about. The flip-rank and the
rainbow do not live in the cube. They live in the *quotient*: the way the symmetry group `S_n` folds the
cube into orbits. Tournaments and even graphs are two different foldings of one cube, and the folding, not
the cube, is what my invariants measured.

Why the difference is instructive: the tournament folding has an efficient *shape* — the balanced cut, fix
the crossing arcs and flip the two sides — and that shape is what let the classes pack at the floor. Even
graphs, read in the same fundamental-cycle coordinates, have no such shape; their classes cut across the
coordinates, and one class (the empty graph) is a single point, pinning any subcube. The bijection carries
tournaments to even graphs faithfully as sets, but it carries the tournament's *natural decomposition* to
something unnatural for even graphs. A reframe that preserves the elements can still destroy the shape that
made the elements tractable.

This is the caution I most needed for the lonely runner. This project runs on reframes — the runner as a
measure, the measure as Verblunsky coefficients, the coefficients as a tournament, the tournament as an
even graph, the even graph as a lattice. Each is a clean bijection at the level of objects, and it is
tempting to believe that a result proved in one frame transports to another because the objects correspond.
The even-graph computation says: no. What transports is what the *quotient* preserves. The lonely-runner
covering-min does not depend on the abstract relation lattice alone — the naive geometry-of-numbers
predictor `lambda_1` proved that, failing completely (THM-515) — it depends on the *arithmetic* quotient,
the modulus `Phi6` and the phase-residue `p(w) = nw mod Phi6`. The relation lattice is the right *object*
(the LRC's even graph, carrying the lonely measure as its theta function, with additive energy playing the
role of the triangle count and the parity lemma playing the role of Redei's theorem). But the object is
inert without its arithmetic decoration. The lattice is the cube; the arithmetic is the folding; and the
answer lives in the folding.

So the transferable move for LRC is not "find a bijection and copy the theorem across." It is "find the
bijection, then ask what the quotient does" — carry the *decoration*, not just the objects. The even-graph
detour did not give me a new bound, but it corrected a belief that would have led me to over-trust the next
clean reframe, and it named exactly what a reframe must preserve to be worth trusting: the group action for
tournaments, the modulus for runners. A bijection of sets is a promise about elements; only a bijection of
quotients is a promise about answers.
