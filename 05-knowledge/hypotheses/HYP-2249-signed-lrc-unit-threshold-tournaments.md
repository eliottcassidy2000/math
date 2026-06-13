# HYP-2249: Signed LRC and unit-threshold tournaments

**Status:** OPEN  
**Source:** codex-2026-06-05-S674  
**Primary artifacts:** `04-computation/signed_lrc_unit_threshold_s674.py`,
`05-knowledge/results/signed_lrc_unit_threshold_s674.out`,
`07-reflections/signed-lrc-unit-threshold-tournaments-s674.md`

## Statement

The useful negative-speed LRC reframe is not a new observer predicate; it is a
sign gauge over the same observer orbit.  Independent sign flips preserve the
loneliness distances to the stationary observer, but they change runner-runner
relative clocks from differences to sums.  Therefore the `C=2n-1` pair-sum
support calculus can be modeled as an actual relative-clock spectrum in a
signed double cover.

The corresponding unit-distance reframe is a three-state threshold tournament:
for a fixed point order, pairs with length `<1`, `=1`, and `>1` are distinct
edge states.  The unit state is a trienerment/tie layer whose positive or
negative sign resolves the edge.  On normalized optimal point sets with minimum
distance `1`, the `<1` state disappears, so this model collapses back to the
older binary unit/nonunit flip gauge; any new information must come from
compression, deletion, perturbation, or a nontraceable unit graph.

## Evidence From S674

The script proves the exact observer-gauge fact symbolically:

```text
dist(v*t) = dist(-v*t)
```

and audits it on the `q=756` grid for the LRC `n=14` floor atoms `AP`,
`Vstar`, and `2AP`.  For each named sign pattern tested, the observer-safe mask
is identical to the all-positive mask.

The pair clocks are not identical.  For a signed pair, the relative clock is

```text
s_i*v_i - s_j*v_j  mod C.
```

When signs agree this is a difference clock; when signs disagree it is a
pair-sum clock up to sign.  Thus alternating or split sign patterns convert a
large part of the pair-sum ledger into ordinary relative motion.

The `Vstar` floor atom exposes an exceptional signed zero clock:

```text
3 + 24 = 27 = C.
```

In the S674 output, split and mod-3 sign gauges give `zero_pair_clocks=1` for
`Vstar`, while the analogous AP and `2AP` rows have none.  The best sign-lift
search finds `Vstar` patterns with objective `(1, 11, 30, 42)`, meaning one
zero pair clock, eleven gcd-9 clocks, thirty gcd-3 clocks, and forty-two
opposite-sign/sum-exposed pairs.  This makes the old nonunit `Vstar` seam a
literal zero relative-clock event in the sign cover.

For unit distance, the S674 threshold tournament audit separates two regimes.
Triangular-spine rows through `n=10` have no `<1` pairs, a unit Hamiltonian
spine, and a directed all-unit Hamiltonian path under both unit signs.  The
sign only chooses the spine direction.  By contrast, the compressed-line and
square-with-center toys have genuine `<1` compression states and no all-unit
directed Hamiltonian path.  This validates the three-state model as a richer
diagnostic while explaining why normalized unit-spine constructions previously
looked binary.

Route Tournament Analysis over proof lanes is transitive and ranks:

```text
opposite_pair_sum_clock
> owner_carry_rank
> alternating_sign_search
> unit_threshold_tristate
> unit_spine_tie_order
> signed_observer_gauge
> raw_negative_speed
> raw_unit_flip
```

## Interpretation

Negative runner directions should be used as a proof lens, not as a candidate
counterexample family.  They turn hidden pair-sum clocks into visible
runner-runner clocks without changing observer loneliness.  This suggests a
signed owner/carry rank: a coherent sign-lifted pair-sum obstruction should
either be one of the known AP/Vstar/2AP floor atoms or pay a strict
owner-deletion/carry tax.

The unit-distance threshold model says the meaningful question is not whether
unit edges are flipped positively or negatively.  If the unit graph is
traceable, either sign can orient a unit spine by choosing the path direction.
The first genuine flop still requires a nontraceable extremal unit-distance
graph or a threshold state with real `<1` compression data.

## Assumption Challenge

This session does not assume tournament vertices must be runners or points.
Useful vertices in the new models include sign gauges, pair clocks, threshold
pair states, unit-spine steps, owner/carry fibers, deletion obligations, and
proof routes.  The quotient preserves observer loneliness in the signed LRC
case and preserves threshold category plus tie resolution in the unit-distance
case; it destroys raw direction labels or metric magnitudes unless those are
recorded as side channels.

## Next Proof Targets

1. Run the sign-lift search over the full HYP-2164/HYP-2241 `Res_27` survivor
   surface, not just AP/Vstar/2AP, and attach owner-private deletion labels to
   every zero or high-gcd signed pair clock.
2. Prove a local signed no-leak lemma: in the AP/Vstar/2AP carry fiber, any
   non-scalar signed pair-sum obstruction is strict after the HYP-2241 owner
   address is retained.
3. Add an Ising/MaxCut formulation for choosing sign gauges that maximize
   `C`-gcd or zero-clock exposure; compare those gauges with carry vectors
   `v=r+27k`.
4. Apply the threshold tournament to Moser `n=21/22` witnesses and to
   point-deletion perturbations, where `<1` compression can appear before
   rescaling.
5. Relate `<1`, `=1`, `>1` unit-distance states to the existing `L/M/R`
   automaton: close is pressure, unit is wall, far is slack.
