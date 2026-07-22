---
id: HYP-8871
title: "LRC14 Kelvin-Farey and scaled-clock sidecar automaton"
status: >
  OPEN completion program downstream of THM-2055/2056, now with one full model
  closed by THM-2057. Kelvin inversion replaces tangent disks by a fixed polar
  polygon, and the Farey defect inequality certifies whole unimodular cones.
  Residual rays should carry scaled safe-clock orbits, affine binding phases,
  or owner-labelled Euler certificates. Uniform sidecar propagation over the
  full THM-2052 star atlas is not proved here.
source: codex-2026-07-21-LRC-normal-fan-sail
related:
  - THM-2052
  - THM-2053
  - THM-2054
  - THM-2055
  - THM-2056
  - THM-2057
  - THM-2058
  - HYP-2108
  - HYP-2647
  - HYP-2896
  - HYP-2986
  - HYP-8841
  - HYP-8846
  - MISTAKE-225
  - MISTAKE-226
---

# HYP-8871 -- Kelvin-Farey addresses with clock/binding sidecars

THM-2055/2056 change the finite object. For each THM-2052 two-anchor star:

1. compute the saturated column configuration and the signed hull `K_U`;
2. delete nonvertex columns from the determinant sidecar only;
3. Kelvin-invert the gate to the rational polar polygon
   `(1/91)R^(-1)K_U^o`;
4. split parameter space by the rational normal fan of `K_U` and regularize
   each owner sector into acute unimodular cones;
5. apply THM-2056's defect inequality
   `2u dot v>=A_p(u)+A_p(v)` to certify every interior slope at once;
6. on each unresolved boundary ray, search for a scaled core clock whose safe
   numerators contain a complete unit orbit;
7. if a clock is killed, descend to its divisibility sublattice; if two clocks
   are killed, seek an affine binding family as in THM-2057;
8. retain exact pair-sum margin, relative-Fejer resonance packet, and HYP-2986
   open-tope/boundary-cocircuit state until a strict phase or owner-labelled
   Euler endpoint accepts the state.

THM-2058 refines step 8 inside a fixed bad transverse deck: split the safe
count by primitive phase order and reduce each hull-owner longitudinal fiber
to an explicit coprime interval. Its interval carrier should be joined to the
missing-clock lcm tax before any raw parameter enumeration.

The key hoped-for compression is **interval acceptance**: adjacent Farey nodes
with the same active hull owner and the same signed phase-height wall word
should admit one symbolic certificate for the whole mediant interval. A failed
interval must split at an explicit event:

```text
hull-owner tie,
positivity or collision wall,
pair-sum ruler change,
relative-Fejer resonance,
endpoint-owner exchange.
```

All geometric event equations are rational linear or quadratic in the two
parameters. Clock killing adds exact divisibility states. Hence every fixed
star has a finite exact cell structure. The missing uniform theorem is to
show that every unresolved ray has either a complete safe-unit orbit or a
bounded-depth chain of killed clocks ending in an affine binding/Euler family.

The modular-form proposal HYP-8880 currently ranks below these carriers.
MISTAKE-226 shows that divisor labels alone do not map a finite phase clock to
a cusp of `X_0(N)`, and no coefficient of the rational level-14 eta product is
known to preserve phase height. A modular sidecar becomes admissible here only
after it is pulled back to a signed owner-channel sum with a proved safe-phase
implication.

## Proved model and the new conjectural rule

THM-2057 proves the rule on the whole plane

```text
{a,2a,...,11a,13a,w}:
  12a clock survives -> unit-orbit witness;
  12a killed, 14a survives -> unit-orbit witness;
  both killed -> 84a|w -> explicit affine binding witness.
```

It also closes the adjacent AP-tail plane

```text
{a,2a,...,12a,w}:
  13a clock survives -> unit-orbit witness;
  13a killed, 14a survives -> unit-orbit witness;
  both killed -> 182a|w -> explicit deep-well witness.
```

The determinant gate leaves `640690` distinct-speed primitive directions on
this positive one-tail plane, while the clock/binding sidecar has three
symbolic leaves. Thus enumerating the polar residual is the wrong terminal;
the polar/Farey carrier is only an address at which arithmetic sidecars must
act in bulk.

More generally, THM-2057 proves the **missing-clock lcm tax**. If a labelled
core `C` contains no multiple of `N<=14`, a counterexample tail over `aC` must
satisfy `Na|w`. Intersecting all such taxes forces
`a*lcm{missing clocks}|w`. The next atlas computation should therefore record
the missing-clock lcm of every candidate lower-rank core before opening any
Farey residual.

## Assumption challenge

The vertices are not runners, arcs, primes, or form classes. At the geometric
stage they are signed hull-owner cones; at the arithmetic stage they are safe
numerator orbits, killed clocks, divisibility sublattices, and binding rays.
The first quotient preserves the determinant gate and rational parameter
address but destroys non-hull runner constraints. The second preserves an
actual phase witness but depends on a labelled lower-rank core. Neither can
replace the other.

## Tournament analysis

Proof-carrier vertices:

```text
owner_sector_sail
scaled_clock_binding
exact_pair_sum_fan
relative_Fejer_cell
endpoint_owner_cocircuit
kelvin_polar_polygon
raw_tangent_disk_scan
Heegner_form_class
raw_relation_matrix
```

Pairwise observable: `(gate exactness, LRC-predicate retention, symbolic cell
coverage, sidecar debt, cost)`. The switch prefers a carrier only if it keeps
the specified primitive slope and either certifies a whole cell or emits an
exact split event. The tie path now starts

```text
scaled_clock_binding
> exact_pair_sum_fan
> relative_Fejer_cell
> endpoint_owner_cocircuit
> kelvin_polar_polygon
> owner_sector_sail
> raw_tangent_disk_scan.
```

`Heegner_form_class` and `raw_relation_matrix` rank last because neither
preserves the determinant owner or the weak phase-height predicate.
