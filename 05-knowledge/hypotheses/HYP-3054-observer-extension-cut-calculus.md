---
id: HYP-3054
title: Observer-extension cut payload calculus for controlled forgetting
status: SYNTHESIS / abstraction across existing LRC carriers; not a proof
source: codex-2026-06-26-S218
tangent: T1136
related:
  - HYP-3053
  - HYP-3052
  - HYP-3051
  - HYP-3050
  - HYP-3049
  - HYP-3048
  - HYP-3047
  - HYP-3046
  - HYP-3045
  - HYP-3043
  - HYP-3040
  - HYP-3039
  - HYP-3037
  - HYP-3034
  - HYP-3024
  - HYP-3022
  - HYP-3021
  - HYP-3018
  - HYP-2997
  - HYP-2995
  - HYP-2991
  - HYP-2989
  - HYP-2963
  - THM-381
  - THM-385
  - THM-572
  - OPEN-Q-108
---

# HYP-3054: Observer-Extension Cut Payload Calculus For Controlled Forgetting

## Claim

The "observer-extension/cut payload" exposed by the A000568 perspective
failure is the abstract boundary coordinate of controlled forgetting.

Given a quotient

```text
q : X -> Y
```

and an operation that asks the quotient to survive one more outside
interaction,

```text
extension, deletion, observer insertion, layer transport, route handoff,
capacitor cut, automaton transition, or certificate pushforward,
```

the observer-extension/cut payload is the smallest data `sigma` for which the
target LRC predicate becomes legal after compression:

```text
P(x) is fiber-safe on q
  only after sigma is retained,
  reconstructed from q,
  killed by a dual/cut/cocycle,
  descended familywise,
  stopped at a boundary atom,
  or routed to named residual debt.
```

Equivalently, `sigma` is the obstruction to the square

```text
X  --extension-->  E(X)
|                 |
q                 q_E
v                 v
Y  --shadow---->  E(Y)
```

commuting at the proof-predicate level.  If `q(x)` and the shadow extension do
not determine the new status, route, owner, topology, or certificate class,
then the missing coordinate is not noise; it is observer-extension payload.

## Tournament Prototype

HYP-3047 through HYP-3053 give the clean finite prototype.

```text
exact rooted node cache at m=5: P(5)=48
next unrooted class count:       A000568(6)=56
defect:                          8
```

HYP-3050 shows that deeper node perspective is not the missing coordinate:
under the stricter directed-WL convention, node depth already reaches exact
rooted node orbits by depth `3` at `m=5`.

HYP-3049 then identifies the first extension payload:

```text
rooted 5-perspective + new observer incident word
  = ordered-pair perspective on 6-tournaments.
```

The compact sidecar is the old-root/new-observer sector deck plus
`cross_sector_orientation_word`.  HYP-3051 through HYP-3053 rephrase the same
payload as diagonal word orbits, deletion-parent fibers, and
rectangle/hourglass cycle defects in fixed-path half-tiling coordinates.

Thus the tournament lesson is:

```text
rooted memory is a cache;
observer extension is a cut;
unrooting is legal only after the cut payload is named.
```

## LRC Translation

In LRC14 the observer-extension/cut payload is any coordinate that changes
when a compressed packet is asked to survive the next proof operation.  The
current repo stack gives several instances.

| Existing carrier | Quotient shadow | Observer-extension/cut payload |
|---|---|---|
| Pair-good decoys | raw false-switch count | blocker generator tooth, good-pair lane, barcode/normal-fan relation, active owner |
| Residual capacitors | mixed route pair count | `residual_capacitor_id`, `first_cut_stage`, exact scale/topology cut, zeta exit, endpoint-owner strip |
| Owner-strip filtration | coarse endpoint word such as `B18Z6` | external endpoint-owner strip, owner-transfer delta, first surviving filtration page |
| AP-tail repair | mod-14 owner strip / coarse stalk | `m mod 13` clock, q13 puncture bit, reciprocal fixed-point certificate |
| Coarse ET + Henselian gate | residue-bin status shadow | Henselian unit-root rule, simple/singular root class, primitive deck, route scheduler |
| Arc-Cech topology | compact topology bucket | closed-H1 owner support and owner-deletion persistence |
| Automaton/gap shadows | Moser/fibbinary/terminal word | exact magnitude cocycle, endpoint/topology/route handoff |
| Haar/fixed-margin squares | row/column margins or line counts | mixed cocycle `zeta`, rectangle defect, owner strip, cross handoff, nested refinement |
| Diagonal layer tilings | raw `k^2+k` line count | layer potential word, rectangle/hourglass defect, deletion-parent fiber |
| Matrix invariants | spectrum/rank/norm/scalar summary | sidecar observability columns separating route/status-changing pairs |

The key shift is from counting residuals to classifying their generator.  For
pair-good decoys, once the modular blocker tooth and active-owner relation are
known, the raw number of decoys is secondary.  For residual capacitors, once
the first min-cut and owner/zeta exit are known, the mixed pair is a finite
proof obligation rather than an anonymous population.  For diagonal layer
lines, once all rectangle/hourglass residues vanish, the full line system is
legally compressed to potentials.

## Controlled-Forgetting Rule

This abstraction sharpens HYP-3039's hidden-coordinate ledger:

```text
Controlled forgetting is not "discard data that seems redundant."
Controlled forgetting is "discard only after testing the observer-extension
payload of the next operation."
```

Practical test:

1. Name the quotient `q`.
2. Name the next operation `E`: add observer, delete root, move layer, cross a
   route, push a certificate, or cut a capacitor.
3. Form the smallest pair set inside a coarse fiber whose status or route can
   change after `E`.
4. Ask which field separates, reconstructs, annihilates, descends, or names
   that change.
5. Promote that field to sidecar vocabulary before using `q` in a proof.

This turns sidecar selection into a local observability problem.  HYP-3048's
matrix version is a table whose rows are coarse-fiber packet pairs and whose
columns are candidate observer-extension payloads.

## Proof-facing Sidecar Target

The next packet manifests should distinguish three payload layers:

```text
extension_address:
  parent_class
  root_orbit
  observer_endpoint_role
  incident_word_orbit
  deletion_parent_profile
  diagonal_word_orbit

cut_or_cycle_defect:
  ordered_pair_sector_deck
  cross_sector_orientation_word
  edge_tail_tip_sector_word
  rectangle_cycle_defect
  hourglass_cycle_defect
  residual_capacitor_id
  first_cut_stage
  exact_M_zeta

route_owner_certificate:
  endpoint_owner_strip
  owner_transfer_delta
  primitive_safe_deck_2_13
  q13_puncture_bit
  closed_arc_h1_owner_support
  magnitude_cocycle
  analytic_blindness_report
  named_residual_sector
```

These fields are not meant to be retained forever.  They mark the proof point
at which forgetting becomes auditable.

## Tournament Analysis

Vertices are proof payloads and quotient obligations, not runners:

```text
observer_extension_cut
sidecar_observability_matrix
endpoint_owner_packet
residual_capacitor_cut
rectangle_hourglass_defect
diagonal_transport_word
ordered_pair_sector_deck
closed_arc_h1_owner_support
primitive_period_deck
pair_good_blocker_tooth
automaton_shadow
raw_scalar_count
```

Pairwise observable:

```text
boundary/open preservation,
route purity,
extension-address retention,
owner/topology/period retention,
dual or cocycle annihilation availability,
family descent availability,
named residual handoff,
proof cost.
```

Switch/gauge:

```text
orient toward the carrier that makes the next outside operation proof-safe
with fewer unnamed residual debts.
```

Tie Hamiltonian path:

```text
observer_extension_cut
> sidecar_observability_matrix
> endpoint_owner_packet
> residual_capacitor_cut
> rectangle_hourglass_defect
> diagonal_transport_word
> ordered_pair_sector_deck
> closed_arc_h1_owner_support
> primitive_period_deck
> pair_good_blocker_tooth
> automaton_shadow
> raw_scalar_count
```

Assumption challenge: the tournament vertices here are not runners, arcs, or
unrooted classes.  They are proof obligations created by a quotient boundary.
The preserved LRC predicate is boundary/open status plus route/certificate
schedulability.  The destroyed information is extension address, endpoint
owner, hidden period, topology support, local cut/cycle defect, and residual
name.

## Next Pulls

1. Build the HYP-3048 sidecar observability matrix for the HYP-2963 residual
   packet bank with observer-extension columns from this note.
2. For pair-good decoys, group by blocker-generator tooth and active-owner
   relation before reporting counts.
3. For residual capacitors, add rectangle/hourglass and ordered-pair sector
   fields to see whether capacitor cuts and tournament edge sectors name the
   same payload.
4. For automaton and Henselian zipper fibers, record the first operation under
   which a shadow fails: route handoff, status handoff, owner transfer, or
   certificate pushforward.
5. Treat every new creative lens as legal only after it states its quotient,
   next operation, observer-extension payload, and controlled-forgetting
   discharge rule.
