---
id: HYP-3101
title: LRC14 component bound via normal-fan Cech barcode
status: SYNTHESIS / proof-route target; not a proof
source: codex-2026-06-27-S259
tangent: T1179
technique: LTI-240
tournament_technique: LTT-138
tags: [lrc14, component-bound, normal-fan, cech-nerve, persistence, controlled-forgetting, tournament-analysis]
related:
  - HYP-3096
  - HYP-3095
  - HYP-3093
  - HYP-3092
  - HYP-3091
  - HYP-3090
  - HYP-3088
  - HYP-3085
  - HYP-3071
  - HYP-3035
  - HYP-3025
  - HYP-3018
  - HYP-3015
  - HYP-2997
  - HYP-2963
  - THM-577
  - THM-575
  - THM-573
  - THM-565
  - THM-530
  - OPEN-Q-108
---

# HYP-3101: LRC14 Component Bound Via Normal-Fan Cech Barcode

## Claim

The missing step in HYP-3096 is not merely a positive measure theorem.  The
polynomial-method witness route needs a uniform component bound for the direct
lonely set

```text
L_14(S) = {t in R/Z : ||s t|| >= 1/14 for every s in S}.
```

After THM-575, raw time denominators and raw largest time arcs are unstable
under apex loading.  The proof object must therefore be the normalized
slow/ruler-coordinate safe set, together with the topology that tells how many
pieces remain.  The best existing carrier for that topology is the four-part
packet

```text
closed arc-Cech nerve
+ open tope / boundary cocircuit facets
+ lonely-profile barcode
+ active-bottleneck normal fan.
```

The proposed theorem target is:

```text
For every primitive non-tight covering 13-row in the THM-573 residual
(at most 6 multiples of 7), after q-witness, one-large-speed, and
finite-address exits have fired, exactly one of the following occurs:

1. the normal-fan/Cech/barcode packet belongs to a finite component family
   F_top with
      mu(L_14(S)) >= m0 > 0
      components(L_14(S)) <= A0;

2. S is an AP/Goddyn-Wong closed-boundary equality atom with owner-essential
   closed Cech H1 and zero strict-open mass;

3. the quotient has a named good-cover/cycle defect routed to F7/THM-572.
```

Case 1 supplies the HYP-3096 denominator-net bridge:

```text
largest_arc(L_14(S)) >= m0 / A0
=> every d >= ceil(A0/m0) has a witness in (1/d)Z.
```

This is a different proof angle from the recent cap/Pascal/moment work.  It
does not try to improve the cap numerator.  It tries to prove that the positive
mass cannot fragment into unboundedly many tiny normalized bars unless the
first missing topology sidecar is exactly named.

## Incoming S258/S64 Signal

The incoming observer-gluing scout
`04-computation/lrc14_observer_gluing_ledger_codex_s258.py` already measured
the obstruction this HYP targets.  On representative rows, the direct
`L_14(S)` component data include:

```text
H7=6 boundary residual:     components=42,  largest=19/1372
apex family V=13:           components=24,  largest=3/637
apex family V=200:          components=102, largest=3/9800
divisor-loaded B=8 row:     components=860, largest=1/82320
```

The readout is exactly the THM-575 warning in finite form: direct arcs are
computable packet fields, but a global direct largest-arc scalar is false
without either a bounded-apex packet theorem or a normalized slow/ruler chart.
HYP-3101 therefore sharpens its target:

```text
bounded-apex residuals: prove direct normal-fan/Cech component control;
large-apex residuals: pass through normalized slow/ruler component control.
```

The incoming S64/THM-577 symbolic overlap theorem strengthens the cap chart:
`cap_10=55/91` and `cap_11=66/91` now have closed-form overlap derivations,
with `k=8,9` dips reduced to finite higher-order remainders.  That makes the
cap/Pascal chart less mysterious, but it does not remove the topology debt.
The cap chart says how much mass should exist; HYP-3101 asks how many
normalized components are allowed to carry it.

The incoming S259 Lean observer-gluing frontier makes this a producer
obligation rather than an optional diagnostic: a bounded component/chamber row
is one way to build an `ObserverGluingCertificate`.  The S65 cap-optimality
exchange scout is a useful warning here.  Its single-swap improvement
tournament reaches the global cap minimizer for `j<=4`, but at `j=5` it gets
stuck away from the global minimizer; cap exchange is a bounded finite check,
not a transitive proof engine.  HYP-3101 should therefore export component
packets to the gluing certificate instead of trying to replace topology by a
greedy cap scalar.

## Why This Is The Right Missing Carrier

The older topology records already point at the component theorem.

1. **Closed arc-Cech nerve (HYP-3025).**  AP and Goddyn-Wong are the named
   rows with closed danger-arc `beta1=1` and zero strict-open mass.  K33,
   petal, covering, fibbinary, and Moser rows have closed arc `beta1=0` and
   positive safe mass.  The runner quotient is lossy unless it retains the
   Betti defect.

2. **Lonely-profile barcode (HYP-3015).**  Components of `L_14(S)` are bars of
   the function `m_S(t)=min_s ||s t||` crossing height `1/14`.  Persistence is
   the collar that prevents a positive row from degenerating to endpoint-only
   equality.

3. **Active-bottleneck normal fan (HYP-3018).**  Barcode bars alone forget why
   they persist.  The active support, endpoint owners, local equation, and
   residue sums are the normal-fan chamber coordinates.  They are the data
   needed to count bars familywise rather than by raw endpoints.

4. **Residual first-tooth observability (HYP-3071).**  In the strict-open
   residual coarse ET+unit fibers, arc topology splits `13/15` and the coarse
   safe-component stalk splits the remaining two before exact `M` is needed.
   This says the component-bound proof should start with topology/stalk
   teeth, not scalar magnitude.

The component bound should be a theorem about a finite Reeb graph / oriented
matroid arrangement in the normalized coordinate, not a theorem about runners.

## Candidate Normal Form

For a row `S`, define the normalized component packet:

```text
source_row_id
primitive_gcd_status
covering_status
multiple_of_7_profile
THM573_level7_status
normalization_chart              # raw time, slow coordinate, or ruler coordinate
closed_arc_cech_beta
open_arc_component_count
boundary_cocircuit_facet_word
owner_current_word_mod_14
runner_quotient_betti_defect
bar_count_at_height_1_14
minimum_bar_persistence
largest_bar_length
bar_endpoint_owner_word
peak_bottleneck_support_word
normal_fan_chamber_id
normal_fan_support_rank
residue_sum_signature
first_tooth
coarse_safe_stalk
exact_safe_stalk
component_bound_status
measure_floor_status
largest_arc_floor_status
finite_ruler_threshold_D
preserved_lrc_predicate
destroyed_coordinate
terminal_exit
```

The working conjecture is that the THM-573 residual has only finitely many
admissible `normal_fan_chamber_id` words after quotienting by the finite-address
normalizer and by the AP/GW closed-boundary equality atoms.  The component
bound is then the chamber count.

## Proof Route

### Step 1: Event Arrangement

The endpoints of the danger arcs are the wall events

```text
t = (k +/- 1/14) / s  mod 1.
```

Between consecutive events, `m_S(t)` is affine and its active owner set is
constant.  A strict lonely component is a run of adjacent open cells whose
minimum is at least `1/14`.  This gives a one-dimensional Reeb graph, but the
number of raw events can grow with large speeds.  HYP-3101 asks for the
normalized quotient where repeated apex subdivisions are not counted as new
proof components.

### Step 2: Cech Dichotomy

Use the closed danger-arc Cech nerve:

```text
closed_beta1 = 1  => danger arcs cover the time circle;
closed_beta1 = 0  => there is positive safe topology unless a quotient defect
                    or endpoint-only gluing has been hidden.
```

The expected dichotomy is:

```text
zero strict-open packet
  => AP/GW owner-current boundary atom
  or F7/THM-572 good-cover quotient defect.
```

This is the topology form of the old AP/GW census, but it retains owner
currents and cocircuit facets rather than a scalar equality label.

### Step 3: Normal-Fan Chamber Bound

A positive component has two endpoint walls and at least one peak bottleneck
support.  Replace each component by the word

```text
(left_owner, right_owner, peak_support, residue_sum, local_equation_type).
```

The normal-fan theorem target is that, inside the `<=6` multiples-of-7
residual and after the finite-address branch normalizer, only finitely many
such words are possible before either a bar has a uniform persistence collar or
the packet emits a named residual debt.

This is where the cross-disciplinary structure is useful:

- oriented matroids: endpoint walls are topes and cocircuit facets;
- hyperplane arrangements: normal-fan chambers count sign patterns;
- persistent homology: bars crossing height `1/14` are the components;
- Reeb graphs / Morse theory: component births and deaths happen at active
  wall events;
- tropical geometry: the minimum of affine distance functions has a normal fan;
- interval-exchange dynamics: repeated apex subdivision must be quotiented by
  a return-word normalizer;
- electrical networks / chip-firing: endpoint owner-current balance is the
  conservation law that distinguishes AP/GW boundary cycles from positive
  components.

### Step 4: Measure Floor Plus Component Bound

HYP-3096 separates two obligations:

```text
mu(L_14(S)) >= m0
components(L_14(S)) <= A0.
```

HYP-3101 targets the second directly.  It may also support the first by proving
that no positive residual can have all bars collapse to zero persistence unless
the closed Cech cycle has the AP/GW owner-current form or the first missing
normal-fan sidecar is F7/THM-572 debt.

## Tournament Analysis

Candidate vertices explicitly considered:

```text
runners
gaps
fixed circle sections
section boundaries
wall-crossing events
danger arcs
safe components
barcode bars
normal-fan chambers
closed Cech cycles
boundary cocircuit facets
proof obligations
```

The selected vertex set is:

```text
normal_fan_chamber_packet
closed_arc_cech_nerve
open_tope_cocircuit_packet
lonely_profile_barcode
coarse_safe_stalk
component_bound_A0
measure_floor_m0
finite_ruler_net
AP_GW_boundary_H1
F7_good_cover_defect
raw_safe_mass
```

Pairwise observable:

```text
preserves_direct_L_14_topology
preserves_component_count
preserves_endpoint_owner_current
retains_active_support
separates_zero_open_from_positive_open
supports_uniform_denominator_net
integrates_THM573_residual
names_destroyed_coordinate
```

Gauge: orient `A -> B` when `A` retains more of the direct lonely-set topology
needed for the component theorem and destroys less endpoint/owner/fan data.
Tie Hamiltonian path:

```text
normal_fan_chamber_packet
> open_tope_cocircuit_packet
> closed_arc_cech_nerve
> lonely_profile_barcode
> coarse_safe_stalk
> component_bound_A0
> measure_floor_m0
> finite_ruler_net
> AP_GW_boundary_H1
> F7_good_cover_defect
> raw_safe_mass
```

Assumption challenge: runners are not the vertices.  Runner data preserves
speed identity, but it destroys the component topology unless endpoint events,
owner currents, and active supports are retained.  The quotient preserved here
is the direct `L_14(S)` topology in the normalized chart; the destroyed
coordinate is the raw time subdivision count from large apex speeds.

## Next Artifacts

1. Build `lrc14_normal_fan_cech_component_ledger` over the HYP-2963 packet
   bank and the THM-573 residual.

2. Prove or refute the Cech dichotomy:

   ```text
   zero-open => AP/GW closed H1 owner current
             or first F7 good-cover quotient defect.
   ```

3. Count normal-fan chamber words after finite-address normalization and
   level-7 status.

4. Convert `mu(L_14(S)) >= m0` plus the chamber count into the HYP-3096
   largest-arc / denominator-net witness theorem.

5. Feed the resulting component chart into HYP-3095's observer-sheaf gluing
   ledger as the normalized-arc topology chart.
