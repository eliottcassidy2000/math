# The loop has two groups and one recursion

*klein-2026-07-01-S69. A reflection on HYP-3801 — pushing "runners on a loop" to the unit-circle theory of
orthogonal polynomials, and the dictionary of maps that made the two pictures one.*

The owner's suggestion — that Verblunsky coefficients give a recursive metaphor for the lonely runner, and
that a *dictionary* of functions between points on a loop might operate group-like — turned out to be
exactly right, and in a way that clarified what the whole L_C thread has been doing.

**Two groups act on the same loop.** The runners live on `S^1 = R/Z`, and there are two natural families of
maps on it. The *arithmetic* family is what LRC uses: rotation `t ↦ t+a` (the flow of time) and
multiplication `t ↦ kt` (the runners themselves, one per speed). These generate the affine group
`Z ⋉ S^1` — the "ax+b" group of the circle, with the beautiful relation `M_k R_a M_k^{-1} = R_{ka}` that
says "speeding up rescales time." The *geometric* family is what orthogonal-polynomials-on-the-circle uses:
the Blaschke maps `z ↦ (z-a)/(1-\bar a z)`, generating `PSU(1,1)`, the conformal automorphisms of the disk.
For a long time these felt like different worlds. The dictionary shows they meet: rotations belong to both,
and "speed `k`" has two incarnations — the arithmetic covering map `t ↦ kt` and a degree-`k` Blaschke
product — both degree-`k` self-maps of the circle. LRC picks the arithmetic one; Verblunsky picks the
geometric one; they are describing the same circle with two group actions, and a measure on that circle can
be read by either.

**One recursion encodes the measure.** Verblunsky's theorem says a probability measure on the circle is the
*same information* as a sequence of coefficients `alpha_n` in the disk, built one at a time by the Szego
recursion. So the lonely measure — which I had been describing by its intervals, then by its Fourier
moments, then by its two-atom resonance — has yet another face: a recursion. And computing it paid off
immediately with an honest correction. My S66 "two-atom law" (`hat1(k) ≈ L cos(2π k t*)`) is not a global
fact; it is `alpha_0` in the limit `r → M_C`. At a generic level the lonely set is 28 spread intervals and
`alpha_0` is tiny (the `k=1` moment cancels across the three-gap positions); only as the level rises to the
binding does the measure contract to two atoms and `alpha_0 → -cos(2π t*)`, with `|alpha_1| → 1` announcing
"exactly two points." The recursion sees the whole movie — absolutely continuous at low level, atomic at
the binding — where the two-atom law was a single frame.

**Three languages, one extremal fact.** The most satisfying convergence: the construction's extremal lonely
measure is *simultaneously* a rank-2 Toeplitz moment matrix (mac-mini's flat-extension), a Verblunsky
sequence that terminates at `|alpha_1| = 1` (this session), and two atoms at the `Phi6`-denominator points
`t*, 1-t*` (THM-515, the phase-residue). Curto–Fialkow ties the first to the third; Verblunsky–Geronimus
ties the second to the third; so these are one theorem wearing three coats. The project keeps rediscovering
that the covering-min value `n/Phi6` and the covering-min *extremizer* are the same hexagonal object; here
that object is a two-point measure, and every framework — moments, orthogonal polynomials, atoms — agrees
on it.

**Why the dictionary matters for the proof.** The dynamical entry is the sharp one. If the runners are the
maps `M_v`, then "far element equidistributes" is not a computation to be re-verified — it is the statement
that `M_v` is a mixing endomorphism of the circle for large `v`, the same reason the `×k` map sends every
measure to Lebesgue. The far element is impotent because multiplication mixes; the dictionary makes that a
one-line dynamical fact rather than a grid experiment. And the open crux becomes an extremal problem over
the Verblunsky sequences of covering measures: no covering set's lonely measure may terminate *deeper* than
the construction's two atoms at the CF-extremal point `[0; n-1, n]`. That is the same single inequality S68
isolated, now phrased where the mature machinery (OPUC extremality, Blaschke, Schur functions) lives.

The lesson is the one the owner was pointing at: give a structure enough different maps and encodings and it
starts to *act like an algebra* — the same object appears as a moment, a coefficient, an atom, a fixed
measure of a mixing map — and the redundancy is not waste, it is leverage. Each encoding hands you the tools
of its own field. The interval picture gave total variation; the moment picture gave flat extension; the
Verblunsky picture gives OPUC extremality and the dictionary of loop-maps; and somewhere in the overlap of
their tools is the inequality that closes the runner.
