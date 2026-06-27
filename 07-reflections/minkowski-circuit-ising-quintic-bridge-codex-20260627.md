# Minkowski, Circuit, Ising, And Quintic Sidecars

codex-2026-06-27

The request was to merge four outside lenses into the Lee-Yang/LRC proof
search: Minkowski's theorem, circuit complexity, the Ising model, and de
Moivre's quintic.  The useful synthesis is not four analogies.  It is four
ways a quotient can lose proof-critical information.

Minkowski says to keep the relation lattice.  The anchored scout used a small
proxy, pair-sum collisions `a+b=c+d`, and it was already informative:
`consec_8` is the top relation leader with `22` pair relations, and
`corr(p0,pair_relation_count)=+0.427`.  This reinforces HYP-3062 rather than
replacing it.  Minkowski pressure is useful only after the relation lattice,
covolume, convex body, successive minima, low-height walls, and deleted
anti-cosets are named.

The Ising model gives the cleanest new finite object.  Put spin `+1` on
anchored rows whose PGF has no real roots and spin `-1` on the rest.  The
one-swap graph has `36036` edges, with `10084` domain-wall edges.  HYP-3109's
root-collision boundary is therefore not just a graph cut; it is a finite
energy interface.  A proof route should label the boundary moves by legal
root-collision ear, forbidden wall, or observer-gluing discharge.

Circuit complexity gives a guardrail.  In the finite anchored bank, a single
threshold predicate `apex7_error <= 5` isolates the unique max-p0 row.  That is
exactly the kind of fact that is dangerous to over-read.  A proof needs a
uniform circuit statement with fixed input sidecars, fixed basis, fan-in/depth,
and a non-tailored threshold rule.  Otherwise "small circuit" is just a fitted
finite classifier.

de Moivre's quintic is weaker as an extremality signal.  Since the miss-count
PGF has degree six, its derivative is a quintic.  I measured the residual from
translated de Moivre normal form `x^5+a*x^3+(a^2/5)*x+b`.  The correlation
with extremality is only moderate (`corr(p0,-residual)=+0.348`), and the best
residual rows are not the maximizers.  That suggests the right use is not
"minimize de Moivre residual"; it is to watch the residual along
root-collision paths, where solvable stationary normal forms may mark
branch-point events.

External cues used in this synthesis:

- Circuit complexity page: https://en.wikipedia.org/wiki/Circuit_complexity
- de Moivre quintic form: https://mathworld.wolfram.com/deMoivresQuintic.html
- Ising model reference: https://mathworld.wolfram.com/IsingModel.html
- Minkowski theorem reference: https://encyclopediaofmath.org/wiki/Minkowski_theorem

The next exact test should leave the anchored bank and attach these four fields
to HYP-2963 residual packets.  The strongest likely theorem is an interface
theorem: root-zero sidecar plus relation-lattice sidecar plus labelled Ising
domain-wall ear plus uniform circuit certificate is a legal packet producer;
raw `p0`, raw volume, raw energy, or raw quintic residual is not.
