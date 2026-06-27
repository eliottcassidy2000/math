---
id: HYP-3029
title: LRC14 safe-component stalk descent for automatic fibers
status: EVIDENCE / target-fiber local carrier scout; not a proof
source: codex-2026-06-26-S193
tangent: T1110
script: 04-computation/lrc14_safe_component_stalk_codex_s193.py
result: 05-knowledge/results/lrc14_safe_component_stalk_codex_s193.out
related:
  - HYP-3026
  - HYP-3025
  - HYP-3024
  - HYP-3023
  - HYP-3022
  - HYP-3021
  - HYP-3020
  - HYP-3018
  - HYP-3017
  - HYP-3016
  - HYP-3015
  - HYP-3014
  - HYP-2963
  - THM-572
  - LTI-177
  - LTI-173
  - LTI-172
  - LTI-171
  - LTI-170
  - LTI-166
  - LTI-164
  - LTT-075
  - LTT-071
  - LTT-070
  - OPEN-Q-108
---

# HYP-3029: Safe-Component Stalk Descent

## Claim

On the first HYP-3023/HYP-3024 hard automatic fiber

```text
MFCMMCCFFFCCC
```

one local proof object appears to replace much of the HYP-3026 carrier-fusion
packet:

```text
largest strict safe component stalk
  = endpoint owner residues
  + peak bottleneck owner residues
  + exact component length
  + exact local peak height
  + open/boundary status
```

This stalk is not the full barcode and not the exact magnitude cocycle.  It is
the local germ of one safe interval.  The S193 computation shows that the exact
largest-component stalk is route-pure on the target automatic fiber, matching
the route purity of exact magnitude while retaining more geometric descent
data.

## Computation

Script:

```text
04-computation/lrc14_safe_component_stalk_codex_s193.py
```

Stored output:

```text
05-knowledge/results/lrc14_safe_component_stalk_codex_s193.out
```

Default run:

```text
bank=HYP-2963 target automatic word MFCMMCCFFFCCC
candidate_rows=639 of 21913
packets=639
```

Stalk descent ladder:

```text
automatic_word                 fibers=1   mixed_route=1  max_mixed=639
residue_terminal_fiber         fibers=181 mixed_route=27 max_mixed=30
owner_only_stalk               fibers=407 mixed_route=7  max_mixed=5
coarse_component_stalk         fibers=473 mixed_route=2  max_mixed=2
exact_component_stalk          fibers=515 mixed_route=0  max_mixed=0
magnitude_cocycle              fibers=546 mixed_route=0  max_mixed=0
stalk_plus_magnitude           fibers=637 mixed_route=0  max_mixed=0
```

The important contraction is:

```text
residue_terminal -> owner_only_stalk -> coarse_component_stalk -> exact_component_stalk
27 mixed routes  -> 7 mixed routes    -> 2 mixed routes         -> 0 mixed routes
```

The coarse residuals are small open-route collisions, not AP/Goddyn-Wong
boundary/open leaks:

```text
single swap 13->159  / single swap 13->117
single swap 13->118  / single swap 13->104
```

Both residual fibers have size `2` and mix `Q-WITNESS` with
`COVERING-MOMENT`.  The exact local length/height data resolves them.

## Why This Is Creative

HYP-3026 says the fusion packet is safe but large.  HYP-3029 asks for a
descent theorem:

```text
fusion sidecars over a target automatic fiber
  descend from one exact safe-component stalk
```

If true familywise, this would demote several sidecars from primitive packet
fields to reconstructible local data:

```text
barcode shape
active normal-fan support
endpoint current
closed/open arc local cover data
Fejer interval anchor
```

The stalk still destroys non-largest bars and the full global barcode.  That
loss is acceptable only if a theorem proves they are not needed for route
purity, or if a residual component emits a named HYP-3026 fusion debt.

## Tournament Analysis

Vertices are proof carriers and local stalks, not runners:

```text
automatic_word
residue_terminal_fiber
owner_only_stalk
coarse_component_stalk
exact_component_stalk
magnitude_cocycle
stalk_plus_magnitude
```

Pairwise observable:

```text
route_purity, status_purity, max_mixed, topology, owner_data,
exact_local_geometry, avoid_exact_magnitude, small_fusion_size, proof_cost
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1]
hamiltonian_path_count=1
first_hamiltonian_path=
  exact_component_stalk >
  stalk_plus_magnitude >
  coarse_component_stalk >
  owner_only_stalk >
  magnitude_cocycle >
  residue_terminal_fiber >
  automatic_word
```

The ordering is a theorem-carrier priority order.  Exact magnitude is still
route-pure, but it is less local and less topological than the exact
component stalk.

## Assumption Challenge

Alternate vertices considered: runners, gaps, all bars, single bars, endpoints,
peak bottlenecks, wall-crossing events, residues, cover arcs, Cech nerves,
Fourier modes, matroid cocircuits, and proof obligations.

Chosen vertices are local safe-component proof carriers.  The quotient
preserves the LRC predicate on the target automatic fiber only after exact
local owner/length/height data is retained.  It destroys non-largest bars,
global barcode multiplicity beyond the count, and exact magnitude unless a
stalk theorem reconstructs or replaces that information.

Challenged assumption: the HYP-3026 fusion packet may not need every sidecar as
a primitive field; several fields may be descent data from a single safe
component stalk.

## Proof Target

Largest-stalk descent lemma for `MFCMMCCFFFCCC`:

```text
Inside the target automatic/residue-terminal fibers,
if two packets share the exact largest safe-component stalk, then they share
the same theorem route.
```

A stronger version should show:

```text
coarse stalk residuals are exactly bounded open-route scheduler collisions
and each residual is discharged by a q-witness/covering formula, Fejer
certificate, or boundary-moment handoff.
```

## Next Work

1. Run the stalk key over the full HYP-2963 bank or at least over the largest
   automatic-word fibers after `MFCMMCCFFFCCC`.
2. Prove the two coarse residual patterns by formulas for the `13->117/159`
   and `13->104/118` lanes.
3. Compare exact stalk fibers against HYP-3025 closed arc-Cech local facets
   and HYP-3018 normal-fan supports.
4. Add `largest_component_stalk_key`, `stalk_owner_word`,
   `stalk_length`, `stalk_peak_height`, and `stalk_exit_route` to packet
   sidecars if the full-bank stress remains clean.
