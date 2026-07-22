---
id: HYP-8871
title: "LRC14 Kelvin-Farey and scaled-clock sidecar automaton"
status: >
  OPEN completion program downstream of THM-2055/2056, now with one full model
  closed by THM-2057. Kelvin inversion replaces tangent disks by a fixed polar
  polygon, and the Farey defect inequality certifies whole unimodular cones.
  THM-2058/2059 provide exact phase packets and CRT joins; THM-2062 adds a
  hereditary CRT wheel, and THM-2065 makes every circuit-free strict residual
  a finite ray packet. Persistent bounded marked circuits still need clocks,
  binding phases, or Euler certificates. Uniform atlas discharge is not proved.
source: codex-2026-07-21-LRC-normal-fan-sail
related:
  - THM-2052
  - THM-2053
  - THM-2054
  - THM-2055
  - THM-2056
  - THM-2057
  - THM-2058
  - THM-2059
  - THM-2060
  - THM-2061
  - THM-2062
  - THM-2064
  - THM-2065
  - THM-2066
  - HYP-2108
  - HYP-2647
  - HYP-2896
  - HYP-2986
  - HYP-8841
  - HYP-8846
  - MISTAKE-225
  - MISTAKE-229
  - MISTAKE-233
---

# HYP-8871 -- Kelvin-Farey addresses with clock/binding sidecars

THM-2055/2056 change the finite object. For each fixed THM-2052 two-anchor
star and each explicitly chosen ordered saturated basis, the proved geometric
front end is:

1. compute the saturated column configuration and the signed hull `K_U`;
2. omit nonvertex geometric points from the determinant maximum only, while
   retaining every runner label in all LRC, deck, collision, and endpoint data;
3. Kelvin-invert the gate to the rational polar polygon
   `(1/91)R^(-1)K_U^o`;
4. split parameter space by the rational normal fan of `K_U` and regularize
   each owner sector into acute unimodular cones, recording owner-tie rays
   rather than pretending that a half-open convention makes ownership unique;
5. apply THM-2056's defect inequality
   `2u dot v>=A_p(u)+A_p(v)` to certify every interior slope at once;
6. retain THM-2053's exact transverse residue deck on every unresolved ray or
   fiber; the finite Farey certificate concerns only the determinant gate;
7. search for a labelled lower-rank core and a scaled clock whose safe
   numerators contain a complete unit orbit;
8. if a clock is killed, descend to its divisibility sublattice; if enough
   clocks are killed, seek a separately proved affine binding family as in
   THM-2057;
9. on a one-tail state, use THM-2059's exact CRT reduction histograms to join
   the core-safe and tail-safe packets; a zero dot product is a failed clock,
   not a failed row;
10. in parallel, retain exact pair-sum margins and test whether THM-2054's
   resonance-lift hypotheses or a proved owner-labelled Euler rule apply.

THM-2058 proves a refinement of step 8 inside a fixed bad transverse deck:
split the safe count by primitive phase order and reduce each hull-owner
longitudinal fiber to an explicit coprime interval. THM-2059 supplies the exact
join operation: reduce the core and tail packets modulo their common gcd and
take the histogram dot product. THM-2060 makes every tail bin uniformly
positive unless `a|w`, with the sharp floor `q-ceil(q/7)`, and turns selected
clocks into a finite exceptional-residue atlas. These carriers should be joined
to the missing-clock lcm tax before raw parameter enumeration.

THM-2062 supplies the missing deletion sidecar on the resulting interval. A
rank-two deletion contributes one determinantal index; at each prime all
deletions exclude at most two projective directions. For fixed `N`, the
allowed `M` form an exact squarefree CRT wheel, while a rank-one deletion is an
affine `+-1` terminal with a one-dimensional coprime wheel. This is an exact
interval split event, not a phase certificate. THM-2061 separately resolves
the only primitive imprimitive-eleven-core lane down to a dyadic folded seam,
then closes every normalized quotient core in `{1,...,19}`.

THM-2066 makes the dyadic arithmetic sidecar lossless. On each safe quotient
clock, every odd tail has a binary word saying which lift it kills; a
counterexample needs two complementary words. These word constraints compose
by CRT. The bank `15..34` closes every primitive divisor-complete quotient
core in `{1,...,24}`. Any later seam atlas should carry the word and its
residue class modulo `2N`, not merely tail eligibility or packet cardinality.

THM-2065 supplies the next exact split. A THM-2051 bounded higher relation is
either persistent in the coefficient template or selects one primitive
projective parameter. Circuit-free templates therefore reduce to finitely
many rows after the CRT wheel; the persistent height-`2^20` marked-circuit
templates are the only relation-rich plane branch requiring bulk discharge.

Only the cited steps are available uniformly. Step 10 names possible
sidecars, not automatic consequences of the geometric state; HYP-2986,
HYP-2108, and other hypotheses remain conjectural until proved in the exact
model at hand.

THM-2059 also gives a bulk filter before exact support comparison. Its positive
zero Fourier mode wins whenever the product of total packet masses exceeds
the product of their centered `L^2` discrepancies. Only high-fluctuation
histograms then need exact modes, signed Euler, or endpoint-owner treatment.

The key hoped-for compression is **interval acceptance**: adjacent Farey nodes
with the same active hull owner and the same signed phase-height wall word
should admit one symbolic certificate for the whole mediant interval. A failed
interval must split at an explicit event:

```text
hull-owner tie,
positivity or collision wall,
pair-sum ruler change,
relative-Fejer resonance-lift boundary,
endpoint-owner exchange.
```

Hull-owner, positivity, collision, and determinant boundaries are rational
linear or quadratic equations in fixed coordinates; clock killing adds exact
divisibility states, while THM-2058/2062 give exact periodic interval states
once the plane, clock, owner, and event list are fixed. What is **not** proved
is that the displayed event list is complete for phase-height, pair-sum,
Fejer, or endpoint behavior, that its cells propagate across the bounded star
atlas, or that every surviving interval has a safe-unit, THM-2061 folded-seam,
bounded clock/binding, or Euler exit. Those are the automaton obligations.

The modular-form proposal HYP-8880 currently ranks below these carriers.
MISTAKE-233 shows that divisor labels alone do not map a finite phase clock to
a cusp of `X_0(N)`, and no coefficient of the rational level-14 eta product is
known to preserve phase height. A modular sidecar becomes admissible here only
after it is pulled back to a signed owner-channel sum with a proved safe-phase
implication.

## Proved model and the new conjectural rule

THM-2057 proves the rule on two whole planes

```text
{a,2a,...,11a,13a,w}:
  12a clock survives -> unit-orbit witness;
  12a killed, 14a survives -> unit-orbit witness;
  both killed -> 84a|w -> explicit affine binding witness.

{a,2a,...,12a,w}:
  13a clock survives -> unit-orbit witness;
  13a killed, 14a survives -> unit-orbit witness;
  both killed -> 182a|w -> explicit deep-well witness.
```

The determinant gate can leave a very large finite primitive residual on
these positive AP one-tail planes, while each clock/binding sidecar has three
symbolic leaves. Thus enumerating the polar residual is not the conceptual terminal;
the polar/Farey carrier is an address at which arithmetic sidecars should act
in bulk. Any quoted residual count still needs its exact universe and
reproduction artifact.

More generally, THM-2057 proves the **missing-clock lcm tax**. If a labelled
core `C` contains no multiple of `N<=14`, a counterexample tail over `aC` must
satisfy `Na|w`. Intersecting all such taxes forces
`a*lcm{missing clocks}|w`. The next atlas computation should therefore record
the missing-clock lcm of every candidate lower-rank core before opening any
Farey residual.

For an arbitrary one-tail core, THM-2059 refines “try a clock” into a lossless
finite object: the pair of core/tail residue histograms modulo their CRT gcd.
Their dot product counts compatible classes; the `Na/lcm` lift counts safe
grid phases. This is the proved join operation the state graph should carry;
it does not yet select a clock or prove every histogram pair overlaps.

## Assumption challenge

The vertices are not runners, arcs, primes, or form classes. At the geometric
stage they are signed hull-owner cones; at the arithmetic stage they are safe
numerator orbits, killed clocks, divisibility sublattices, and binding rays.
The first quotient preserves the determinant gate and rational parameter
address but destroys non-hull runner constraints. The second preserves an
actual phase witness but depends on a labelled lower-rank core. Neither can
replace the other.

## Typed state graph, not a tournament

No tournament is declared. The carriers above do not come with an intrinsic
binary relation on pairs: “use clock first” or “use Farey first” is a strategy
choice depending on retained data and cost, ties are genuine, and some
carriers are incomparable. Forcing these choices into oriented edges would
add no theorem and would erase why a transition is legal.

The appropriate object is a typed state graph. A state records at least the
fixed basis/plane, primitive parameter or Farey cone, tied hull-owner set,
non-hull label sidecar, transverse deck, surviving/killed clocks, divisibility
sublattice, and current proof obligation. An edge must name the map, the
predicate it preserves, the information it discards, and the side condition
that licenses it. Accepting states contain an actual weak `1/14` witness or a
proved implication to one; “determinant gate failed” is only an unresolved
state. The open theorem is finiteness and completeness of this graph uniformly
over the THM-2052 atlas, not a preferred ranking of proof techniques.
