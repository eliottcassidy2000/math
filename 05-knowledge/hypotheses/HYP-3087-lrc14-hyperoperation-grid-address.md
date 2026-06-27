---
id: HYP-3087
title: LRC14 hyperoperation grid address
status: SYNTHESIS / operation-address proof carrier; not a proof
source: codex-2026-06-27-S254
tangent: T1169
technique: LTI-233
tournament_technique: LTT-131
related:
  - HYP-1976
  - HYP-1980
  - HYP-3003
  - HYP-3004
  - HYP-3088
  - HYP-3089
  - HYP-3090
  - HYP-3085
  - HYP-3091
  - HYP-3083
  - THM-523
  - THM-571
  - THM-572
  - THM-573
  - THM-575
  - LTI-011
  - LTI-153
  - LTI-154
  - LTI-232
  - LTT-130
---

# HYP-3087: LRC14 Hyperoperation Grid Address

## Claim

The user's hyperoperation hierarchy on `(numerator, denominator)` and the
older two-dimensional integer grid tiled by `x+2` on one axis and `x*2` on the
other are the same proof-interface object when they are read over a Farey
packet:

```text
root packet:         t = p/q
additive shadow:     p+q       horizontal `x+2` / summand-antidiagonal lane
product shadow:      p*q       vertical `x*2` / factor-incidence lane
power stress lanes:  q^p, p^q  lacunary/tower stress tests
```

The useful carrier is not the raw grid point and not the raw space-filling
curve.  It is an operation-address packet:

```text
(p, q)
+ operation lane
+ additive owner data
+ multiplicative / 7-adic depth data
+ endpoint owner data
+ current danger deficit
+ route or terminal debt
```

This turns the earlier `x+2`/`x*2` picture into a controlled-forgetting rule.
The grid supplies the branch coordinate; the LRC danger deficit supplies the
clock coordinate; the finite-address branch ledger supplies the proof exit.

## Lineage

This synthesis merges five existing lines.

1. HYP-3003 says a Farey packet `p/q` has two proof-relevant shadows:
   `p+q` is a summand antidiagonal / additive pinch row, while `p*q` is a
   multiplicand hyperbola / factor incidence fiber.  Power lanes such as
   `q^p` and `p^q` are stress tests unless the root packet and forgotten fibers
   remain visible.

2. HYP-3004 says `*2,+2` is a recursion-mode pair: vertical 2-adic or
   valuation-tree descent versus horizontal line motion.  The split is useful
   only after the preserved predicate and destroyed coordinate are named.

3. The S511 operation-grid reflection says static `x+2`/`x*2` labels have no
   within-clock signal by themselves.  They become moving LRC gauges only after
   they are multiplied by the current pair-cell deficit
   `max(0, 1/N - dist(i,j,t))`.

4. THM-573 changes the active arithmetic residual.  The old branch closure
   spine still had many remarks phrased around multiples of `14`; the incoming
   theorem proves every 13-speed row with at least `7` multiples of `7` is
   lonely.  The live residual is now covering rows with at most `6` multiples
   of `7`.  Therefore the vertical product/valuation axis must record the
   level-7 lift status, not merely multiple-of-14 status.

5. THM-576/HYP-3090 adds cap-ratio/deviation status.  Pairwise avoidance
   explains the clean triangular caps for `k>=10`, while `k=8,9` and the
   first order-3 break remain named cap-side debt.  A grid cell is therefore
   not proof-safe if it forgets whether it is a clean cap-ratio packet, a
   deviation packet, or a higher-order break packet.

6. HYP-3094, the covering-moment / K33 state-lift shuttle, supplies a
   concrete exit grammar for this grid.  Nested-refinement covering packets and
   cross-handoff K33 packets can both be positive-open, so the operation grid
   must retain grid class, active binders, and endpoint-owner transitions
   before using safe mass or a product lane as a theorem separator.

7. HYP-3088/HYP-3089 add the normalized-arc and paper/V* target.  The same
   `14=2*7` wall appears as the failed field-interpolation case, but the raw
   Conjecture 7.1 denominator reading is repaired into normalized
   slow/ruler-coordinate arc data plus the V* crossover.  Thus a grid cell must
   also record whether it contributes to that normalized largest-arc floor or
   only routes to O2/O3 residual debt.

## Interpretation Of The Space-Filling Curve

The old space-filling curve through the integer grid should be treated as a
scheduler, not as a proof by coverage.  A Hilbert/Peano-style traversal of grid
points is too lossy if the step order is all that survives.  The proof-safe
version is a Gray-code-like audit path whose adjacent moves change one retained
coordinate at a time:

```text
horizontal move:  additive `x+2`, p+q owner, endpoint-pair sum
vertical move:    multiplicative `x*2`, factor shell, 7-adic depth
diagonal repair:  product-sum interface, p*q / p+q handoff
stress move:      q^p or p^q power lane, allowed only with finite address
exit move:        q-witness, level-7 lift, covering moment, K33/THM-572,
                  branch closure, or named residual debt
```

The scheduler is useful only if each cell keeps the data that the traversal
would otherwise destroy:

```text
farey_root_p
farey_root_q
operation_lane
sum_shadow_p_plus_q
product_shadow_p_times_q
power_stress_word
x_plus_2_step_id
x_times_2_depth
space_filling_successor
endpoint_owner_word
danger_deficit
level7_lift_status
finite_address_word
terminal_exit
destroyed_coordinate
required_sidecar_or_debt
```

Without these fields, the curve is only a static enumeration.  With them, it is
a candidate order for auditing the finite-address packet bank.

## LRC14 Proof Readout

The improved route is:

```text
primitive 13-speed row
  -> 14-free q-witness exit (THM-523)
  -> level-7 lift sieve if at least 7 speeds are divisible by 7 (THM-573)
  -> residual core with at most 6 multiples of 7
  -> danger-weighted operation-address packet
  -> cap-ratio/deviation status (HYP-3090/THM-576)
  -> three-sameness lonely-set fiber (HYP-3091)
  -> finite-address branch packet (HYP-3083)
  -> covering/K33 shuttle if the grid class is nested_refinement or cross_handoff (HYP-3092)
  -> protected branch graph / no naked bridge
  -> terminal discharge or named residual debt
  -> formal witness readout M>=1/14
```

The hyperoperation grid contributes at the middle step.  It is not replacing
the q-witness theorem, the level-7 sieve, the HYP-2963 packet bank, the
HYP-3094 covering/K33 shuttle, or the branch-closure proof.  It gives a
coordinate chart for the remaining residual core, where additive and
multiplicative shadows must both be retained until a legal exit has been
found.

In this chart:

```text
p+q    detects additive endpoint-owner / summand pressure;
p*q    detects product incidence, factor shell, and level-7 vertical depth;
cap    records triangular ratio, deviation, or higher-order break status;
q^p    tests denominator-dominated lacunarity;
p^q    tests numerator-dominated explosion and residue collapse;
curve  schedules adjacent repairs but proves nothing without sidecars.
```

The level-7 update is the main correction to the older grid idea: a product
axis that only asks about `14 | v` is too coarse after THM-573.  The active
vertical coordinate is now:

```text
v7_depth / count_7_divisible / level7_lift_status / residual_leq6_7core_flag
```

## Tournament Analysis

Vertices are operation-address proof carriers, not runners, arcs, or raw grid
points:

```text
finite_address_branch_packet
danger_weighted_operation_cell
level7_lift_gate
exact_farey_root_packet
additive_x_plus_2_lane
multiplicative_x_times_2_lane
product_sum_interface
power_stress_lane
raw_space_filling_successor
static_grid_label
```

Pairwise observable:

```text
preserves_lrc_predicate
retains_root_pq
retains_current_danger
retains_endpoint_owner
retains_additive_owner
retains_product_shell
retains_level7_status
names_destroyed_coordinate
provides_terminal_exit
avoids_naked_bridge
```

Gauge: orient `A -> B` when `A` retains more proof-critical address and danger
data while destroying less.  Ties are broken by this Hamiltonian path:

```text
finite_address_branch_packet
> danger_weighted_operation_cell
> level7_lift_gate
> exact_farey_root_packet
> additive_x_plus_2_lane
> multiplicative_x_times_2_lane
> product_sum_interface
> power_stress_lane
> raw_space_filling_successor
> static_grid_label
```

No new executable tournament was run in this pass.  The declared fingerprint is
the transitive audit expected from the retained-field score, and it inherits
the earlier transitive fingerprints of HYP-3003 and HYP-3004:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

The guardrail is the important part: static operation labels and raw traversal
order rank last because they do not preserve the LRC clock.

## Assumption Challenge

Candidate vertices considered:

```text
runners
danger arcs
raw grid points
space-filling curve steps
hyperoperation lanes
Farey packets (p,q)
denominator shells
7-adic lift classes
endpoint owners
branch graph nodes
proof obligations
```

The chosen quotient is operation-address packets.  It preserves the LRC
predicate only after it retains current danger and endpoint/route exits.  It
destroys raw row identity, raw curve order, and some local pair geometry; those
losses are legal only if reconstructed by the packet fields above, annihilated
by a certificate, descended by THM-573 or a covering-moment argument, or named
as residual debt.

The challenged assumption is that a space-filling curve over the `x+2`/`x*2`
grid can itself close the proof.  It cannot.  The curve is a useful scheduler
only after every visited cell carries the finite address needed to survive
controlled forgetting.

## What Remains

The next executable artifact should be a `hyperoperation_grid_address` ledger
over the HYP-2963 bank and any outside-bank normalizer attempts.  For each row,
record:

```text
source_family
primitive_status
count_7_divisible
level7_lift_status
residual_leq6_7core_flag
farey_root_p
farey_root_q
operation_lane
sum_shadow_p_plus_q
product_shadow_p_times_q
power_stress_word
x_plus_2_step_id
x_times_2_depth
space_filling_successor
danger_deficit_summary
endpoint_owner_word
finite_address_word
preserved_lrc_predicate
destroyed_coordinate
required_sidecar_or_debt
protected_branch_node
terminal_exit
formalization_status
```

The target theorem would not say that the grid is complete.  It would say that
inside the THM-573 residual core, the first nonzero danger-weighted operation
cell either opens a witness, descends through a covering/Node3 or q-witness
route, identifies a HYP-3090 cap/deviation packet, hands off to the HYP-3091
K33/THM-572 shuttle, enters a protected branch packet, or emits a named
residual debt without creating a naked bridge.
