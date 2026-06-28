---
id: HYP-3410
title: LRC14 Bring/Schwarz-Christoffel/BDH/Menger charal recursion
status: SYNTHESIS / exact mixed-fiber sidecar scout; not an LRC14 proof
source: user prompt plus codex-2026-06-28 integration of HYP-3407 boundary/special-function atlas, HYP-3301, and HYP-3404..HYP-3406
tangent: T1371
technique: LTI-371
tournament_technique: LTT-271
script: 04-computation/lrc14_bring_sc_bdh_menger_charal_recursion_codex_20260628.py
result: 05-knowledge/results/lrc14_bring_sc_bdh_menger_charal_recursion_codex_20260628.out
reflection: 07-reflections/lrc14-bring-schwarz-bdh-menger-charal-recursion-codex-20260628.md
related:
  - HYP-3407
  - HYP-3406
  - HYP-3405
  - HYP-3404
  - HYP-3402
  - HYP-3401
  - HYP-3311
  - HYP-3310
  - HYP-3301
  - HYP-3300
  - HYP-3266
  - HYP-3265
  - HYP-3260
  - HYP-3257
  - HYP-3253
  - HYP-3124
  - HYP-2969
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3410: LRC14 Bring/Schwarz-Christoffel/BDH/Menger Charal Recursion

## Claim

The requested cross-disciplinary reframes become proof-facing only after they
are translated into first-failure sidecars over the actual packet bank.

```text
Bring radical               -> five-exit branch resolver
Schwarz-Christoffel mapping -> contact polygon turning word plus accessory debt
Barban-Davenport-Halberstam -> finite owner-channel variance packet
Menger cuts                 -> minimum endpoint-owner cut separating exits
charal signatures           -> characteristic/chiral/arc-lift recursive sidecar
```

The scout uses the concrete HYP-3406 mixed fibers as substrate, not raw
analogies.  It executes the Bring/Schwarz/BDH/Menger slice of the HYP-3407
special-function / boundary-uniformization route.  Here "charal" is
normalized as:

```text
characteristic + chiral + arc-lift
```

It records replacement type, unit slot, Schwarz-Christoffel-style turn word,
owner group counts, apex flag, `v2` span, and exact endpoint-owner labels.

## Exact Mixed-Fiber Readout

The three HYP-3406 leaks currently represented are:

```text
height_leak_12_family:
  GW-shell alias 12->132      positive-Haar-open
  single swap 12->48          positive-Haar-open
  P10+GW                      unit-petal-named

persistent_owner_leak_26_40_54_family:
  single swaps 1/3/5 -> 26/40/54, 9->54  positive-Haar-open
  petal 13->26                            unit-petal-named

height_persistent_owner_leak_10_20_drop_add_family:
  two drop(1/3/5,10)->add(15/17/19,20)    positive-Haar-open
  petal 10->20                            unit-petal-named
```

For the height leak, the minimum owner-label cut has size `1`:

```text
minimum_owner_label_cut_size = 1
first_min_cut = ('5:g1',)
top finite-BDH variance label = 5:g1 with score 8/9
```

For the persistent owner leak, the minimum owner-label cut also has size `1`:

```text
minimum_owner_label_cut_size = 1
first_min_cut = ('1:g1',)
top finite-BDH variance labels:
  1:g1  = 200/121
  13:g1 = 162/121
  11:g1 = 98/121
```

For the newer `(72,20)` height-persistent owner leak, the minimum owner-label
cut has size `3`:

```text
minimum_owner_label_cut_size = 3
first_min_cuts include:
  ('11:g1', '13:g1', '1:g1')
  ('11:g1', '13:g1', '2:g2')
  ('11:g1', '2:g2', '7:g7')
top finite-BDH variance labels:
  13:g1 = 49/50
  11:g1 = 1/2
  7:g7, 5:g1, 2:g2 = 9/50
```

This is a small but useful exact signal: endpoint-owner labels are not just
decorative.  They are literal cut coordinates separating theorem exits in the
known mixed fibers.  The first two fibers have one-label cuts; the newer
frontier already asks for a bounded owner-cut theorem rather than a permanent
one-label theorem.

## Recursive Charal Patterns

The script finds `+14` cusp-ladder families:

```text
1 -> 26,40,54:
  turns (E,U), (E,E,U), (E,U)
  exit always positive-Haar-open

3 -> 26,40,54:
  turns (E,U,E,U), stable across the ladder
  exit always positive-Haar-open

5 -> 26,40,54:
  turns (U,E,U,E,U), (U,E,U,E,U), then (U,E,U,E)
  exit always positive-Haar-open

12 -> 48,132:
  apex-bearing turns (U,E,A) and (U,E,A,U,U)
  exit always positive-Haar-open

13 -> 26:
  turn (U,E)
  exit unit-petal-named

10 -> 20:
  turn (E,A)
  exit unit-petal-named, apex-bearing height-persistent frontier
```

This suggests a recursive theorem target:

```text
If a +14 cusp ladder preserves its charal signature or changes only by
allowed harmless owner insertions, then its theorem exit remains stable;
if the charal signature crosses a forbidden owner cut, it emits endpoint-owner
debt, AP/GW boundary, strict-open floor, or state-lift debt.
```

That is the useful version of "recursive structure patterns in charal
signatures."

## Reframe-by-Reframe

Bring radical:

```text
The five theorem exits form a quintic-like branch alphabet:
q-witness, AP/GW equality, unit-petal, K33/H7, positive-Haar-open.
```

The proof use is not solving a quintic.  It is normalizing the first-failure
branch problem: every packet fiber must have a single exit after sidecars are
attached, or the remaining branch ambiguity is named debt.

Schwarz-Christoffel:

```text
unit contacts are polygon vertices
owner_support is the accessory-parameter word
turn word alone is not enough
```

The accessory parameter is exactly the endpoint-owner/off-unit debt that a
contact polygon would hide if only its turning angles were retained.

Barban-Davenport-Halberstam:

```text
replace prime progressions by finite endpoint-owner channels
measure mean-square imbalance of owner labels against theorem exits
```

This is a finite channel-discrepancy packet, not an imported prime theorem.
It is useful because labels with high variance are exactly labels that can
separate mixed theorem exits.

Menger:

```text
endpoint-owner labels are cut vertices/edges in the proof interface
minimum owner cuts separate boundary-petal exits from positive-open exits
```

This is the strongest route in the scout because it gives a theorem shape:
prove every mixed fiber has a bounded owner cut, an exact dual/Farkas current,
or named debt.

Charal recursion:

```text
characteristic/chiral/arc-lift word tracks recursion under +14 ladders,
endpoint deletion, and owner insertion.
```

It is the candidate recursive invariant joining HYP-3404's first-failure
program to HYP-3406's owner-support repair.

## Tournament Analysis

Vertices are proof carriers and sidecar transformations, not runners.

Pairwise observable:

```text
retained first-failure payload minus destroyed sidecars
```

Switch/gauge:

```text
higher weighted retained payload; ties by fewer destroyed sidecars
```

Exact fingerprint:

```text
vertices = 7
score_hist = {-28:1, 18:1, 20:1, 35:1, 43:1, 56:1, 61:1}
directed_3cycles = 0
hamiltonian_path_count = 1
priority_path =
  menger_owner_cut_recursion
  -> charal_recursive_signature
  -> bdh_channel_discrepancy_packet
  -> schwarz_christoffel_contact_polygon
  -> tropical_height_wall_backend
  -> bring_branch_resolvent
  -> raw_cross_discipline_analogy
```

The result is intentionally not the order suggested by fame.  Menger and
charal recursion beat Bring because they preserve the actual missing
endpoint-owner and recursion coordinates.  Schwarz-Christoffel is useful once
its accessory debt is explicit.  BDH is useful as a finite variance packet.
Bring is mostly a branch-normalization metaphor unless paired with the
sidecars that make `exit_status` a function on packet fibers.

## Assumption Challenge

Candidate vertices considered:

```text
runners
residues
endpoint-owner labels
polygon turns
BDH channels
Menger cuts
Bring branches
deletion events
proof obligations
```

Chosen vertices are proof carriers and sidecar transformations.  The quotient
preserves the LRC predicate only when `exit_status` is constant on fibers or
owner/height/accessory debt is restored.  Raw famous-problem analogies destroy
the exact endpoint-owner and off-grid-floor data that HYP-3406 shows are
needed.

## Next Pull

Formalize the recursive owner-cut theorem:

```text
For each first-failure packet fiber in the enlarged HYP-2963 bank,
either:
  charal signature is exit-pure under +14 recursion,
  a bounded Menger owner cut separates theorem exits,
  finite-BDH channel variance produces a separating owner label,
  Schwarz-Christoffel accessory debt restores endpoint owner,
  or the packet routes to named AP/GW, K33/H7, strict-open, or residual debt.
```

Then extend the HYP-3406 bank and check whether the first two one-label cuts
`5:g1` and `1:g1`, plus the `(72,20)` size-3 cut family led by `13:g1`
variance, are shadows of one bounded recursive owner-cut theorem.
