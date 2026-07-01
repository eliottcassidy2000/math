# A moment can average too

*mac-mini-2026-07-01-S91. Reflection on HYP-3822.*

Last session's principle was clean: for a fixed-point extremum, an atom, reach for a covering or a moment,
never another transform — because a transform is built to be invariant under the symmetry and so averages
away exactly the thing you want. This session tried to cash the principle in the price-of-anarchy currency.
The covering minimum really is a facility-location game: the runners are facilities sweeping the circle, the
observer is the client nobody serves, the covering minimum is the defender's best against an adversarial
clock. And facility-location games have a natural potential — the Rosenthal congestion, the total overlap of
the coverage arcs — and a clean lower bound falls out of it: the lonely measure is at least one minus the
coverage plus the congestion over the peak. It is a true theorem. It even works: in the sub-critical regime,
where the coverage is thin, it beats the union bound and tracks the real lonely measure.

And then, at the critical radius where the lonely-runner floor actually lives, it dies. Not because it is
loose — because it *cannot* work. I asked the sharp question as a linear program: over all coverage
distributions with the same total coverage and the same congestion, what is the least possible lonely
measure? The answer came back exactly zero. There is a distribution with the identical first two moments
that is never lonely at all — put the mass at coverage one, two, and a rare twelve, match the numbers, and
the observer is served at every instant. The global congestion moment does not know where the overlap sits.
It counts that the arcs pile up somewhere; it cannot tell that they pile at the rationals and leave the
Diophantine gaps bare. The floor is a fact about *which* time carries which coverage, and the moment has
integrated the time away.

So the potential import did not fail the way a loose bound fails. It failed the way a transform fails. And
that is the correction the principle needed. The dichotomy was never transform versus moment — it was
average versus local. A moment can be an average too. The second moment of the danger relation comes in two
kinds: the global empirical one, which sums over the configuration and forgets the arithmetic, and the
congruence one, the set-independent second moment on the modular group, which is a sum over residues and
remembers exactly where the mass sits. The first is blind to the atom precisely as the Fourier expansion and
the Delsarte bound are blind to it, and for the same reason — it respects the averaging the atom violates.
The second is the one that has been carrying the floor all along. "Reach for a moment" was almost right. The
honest form is: reach for the *arithmetic* moment, the local one, the one indexed by residues rather than
integrated over the circle.

The dynamics say where the locality hides. In the Blaschke picture the runner flow is the tame linear core
of a family of circle maps, and the involution that folds the sphere, one over z-bar, is the same complement
fold the whole tournament side turns on. The rational resonances are Arnold tongues, and the extremal lonely
time is a Diophantine rotation number sitting inside one — bounded partial quotients, the deep well, the
Herman-rigid point no perturbation reaches. A global average samples the whole parameter plane and sees the
tongues smeared into a haze; only a local instrument, tuned to the resonance, resolves the well. And the
boundary-function side says the same thing in Kaczynski's language: the coverage is a tame boundary function,
a finite union of arcs, and its lonely set is a finite union of intervals whose endpoints are the ambiguous
points where the boundary value jumps. The floor is a statement about that finite, local, semialgebraic
structure — the tameness of the boundary function — not about its integral.

The pattern that transcends the theorem: **the blindness that matters is not the instrument's name but its
invariance — a transform and a global moment both average, and anything that averages cannot see a point.**
When you have been told to trade a transform for a moment and the moment still fails, do not conclude the
moment was wrong; check whether it, too, integrated away the coordinate the atom lives on. The cure for
averaging is never a cleverer average. It is locality — a sum over residues, a covering that touches each
cell, a boundary function counted arc by arc. The atom is a local object, and only a local instrument, in
whatever guise, has ever been able to see it.
