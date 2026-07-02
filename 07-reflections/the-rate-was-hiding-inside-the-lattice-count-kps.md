# The rate was hiding inside the lattice count

**kind-pasteur-2026-07-01-S30** (HYP-3953). A session on the two analytic axioms (hp0cap,
hpartA) produced one identity that reorganizes how I think about every "two-scale" argument in
this project.

## The inversion

For months the witness route carried an unproved "equidistribution rate": the cluster's fast
phases should decorrelate from the slow part at rate O(1/V), proven only for ≤ 6 far elements,
needed in general. Every route hit it: the far-comb rate, the renormalization rate, the
Vmax-ruler embedding.

The session's discovery is that on the witness route the rate does not exist as an analytic
object at all. Two exact identities dissolve it:

1. At ruler times τ_j = (j+c)/V, the cluster condition ‖(V−o)τ_j‖ ≥ 1/14 IS the slow condition
   ‖o τ_j − c‖ ≥ 1/14 — an identity, not an approximation. The only thing the scale V does is
   index a lattice; the only estimate is #points-in-arcs ≥ V·meas − arcCount, which is
   counting, not analysis.
2. The free-gap functional satisfies F_{V−offs}(x) = F_{offs}(x) pointwise, because circular
   gaps are rotation-invariant and the reference scale only rotates the phase configuration.
   The quantity that must be bounded below — the average gap surplus — never sees V. Numerics
   confirmed it to six digits at V = 50, 500, 5000: not convergence, equality.

What remains quantitative is a FINITE census (the c-averaged ledger), and what remains
structural is bookkeeping (nested intervals, one retirement per level). The analytic
difficulty was an artifact of asking for decorrelation between two scales when the correct
change of variables makes the scales never meet.

## Why the c-average is the smooth face of the lonely measure

The homogeneous lonely measure is the wild object of this project: covering sets drive it to
measure zero, its infimum is the whole difficulty of LRC. The c-averaged version
A(U) = E_c[L^c(U)] = E_x[F_U] is tame: at k = 2 its minimum EQUALS the independence value
(6/7)² to six digits; through k = 7 it stays within 13% of independence while the homogeneous
min falls to 64%; on the union-bound death row (k = 7, total danger exactly 1) it sits 55×
above the unconditional Parseval backstop. The Fourier reason is one line: averaging over the
target kills every relation with Σm ≠ 0. The adversary's whole craft — covering the origin —
is a statement about ONE target; no integer configuration can cover all targets at once, and
the average sees all targets. This is the same lesson as mac-mini's "a moment can average too"
and klein's "two localizations", but pointed the other way: here averaging HELPS, because the
quantifier we need ranges over the adversary's targets, not over our configurations.

## The retirement structure

The reference runner of each scale window enters the reduced problem with offset zero — its
constraint is ‖c‖ ≥ 1/14, a condition on the ruler phase, not on time. Each level of the
recursion retires one runner from the dynamical problem into a static choice. Thirteen runners,
at most thirteen retirements: termination is structural, not analytic. I do not think it is a
coincidence that this mirrors the tournament project's oldest pattern — the staircase peeling
(Mode A: hypotenuse removal, one vertex retired per step). The lonely-runner tower and the
tournament triangle reduce the same way: strip the extreme element, convert its interaction
into a boundary condition, recurse on strictly fewer objects.

## What transfers

- The ≤6-far rate (HYP-3787) is bypassed on the witness route — worth re-examining every other
  place it appears as a dependency (the census exhaustiveness write-ups cite it as the
  large-speed cap; the c-ruler version may simplify those too).
- The (Ω) identity says: whenever a functional of phases is rotation-invariant, reference
  scales cancel. Candidates: the sector-cover events (hp0cap's p0 is NOT rotation-invariant —
  sectors are pinned to the grid — which may be exactly why hp0cap stayed hard while the
  witness route dissolves; the gap functional is the rotation-invariant shadow of the sector
  cover).
- The niche threads that fed this: THM-565's sampling (the count), THM-594-C/E (the DMNR
  backstop), klein's difference-set frame (F is a difference-set functional), opus's
  renormalization (the census-minima principle), and my own arc-count budget (the count's
  denominator). None of them was the answer; the answer was a change of variables that let
  them compose.
