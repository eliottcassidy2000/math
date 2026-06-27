---
id: HYP-3045
title: LRC14 endpoint-owner transfer carrier for B18Z6 residual packets
status: EVIDENCE / local owner-transfer carrier and zipper target; not a proof
source: codex-2026-06-26-S208
tangent: T1126
script: 04-computation/lrc14_endpoint_owner_transfer_codex_s208.py
result: 05-knowledge/results/lrc14_endpoint_owner_transfer_codex_s208.out
related:
  - HYP-3044
  - HYP-3042
  - HYP-3041
  - HYP-3040
  - HYP-3039
  - HYP-3038
  - HYP-3037
  - HYP-3036
  - HYP-3035
  - HYP-3032
  - HYP-3031
  - HYP-3027
  - HYP-3026
  - HYP-3018
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3045: LRC14 Endpoint-Owner Transfer

## Claim

The S201/S196b residual packets reveal a smaller hidden coordinate than the
full packet label, refining HYP-3042's owner-strip filtration, HYP-3044's
topology-exception owner-stalk collar lesson, HYP-3039's hidden-coordinate
ledger, and HYP-3040's hidden-statement ledger, while matching the lesson of
HYP-3041's AP-tail puncture atlas: keep the local address that a coarse shadow
forgot.

All four audited residual capacitor packets have the same coarse endpoint word:

```text
B18Z6
```

That word records only the number of boundary events and zero-sum pairs.  It
forgets which external speed/residue owns the facets.  The external endpoint
owner strip restores exactly the missing local address:

```text
petal q=23:          12:26x6,6:20x4
covering q=23:       2:16x6
K33 lift:            12:26x6,8:36x4
single-swap cover:   2:72x6
```

This strip splits the S201 `q=23` diagonal and both S196b residual capacitors,
including the pair whose first cheap cut was exact `M+q` and the pair whose
first cheap cut was coarse boundary topology.

## Computation

Script:

```text
04-computation/lrc14_endpoint_owner_transfer_codex_s208.py
```

Stored output:

```text
05-knowledge/results/lrc14_endpoint_owner_transfer_codex_s208.out
```

The script imports the S201 q=23 drop/add square and the S196b residual
capacitor target rows.  It computes:

- the q=23 diagonal endpoint-owner transfer;
- external endpoint-owner strips for the four residual capacitor packets;
- largest-safe-component owner stalks for the same packets;
- quotient stress over raw residual shadow, coarse endpoint count, exact
  `M+q`, coarse safe body, owner stalk, external owner strip, joined
  owner-transfer carrier, and route labels.

## Readout

For the q=23 diagonal:

```text
AA petal: M=2/23, endpoint=B18Z6, ext=12:26x6,6:20x4
BB cover: M=2/23, endpoint=B18Z6, ext=2:16x6
```

The token transfer from petal to covering is:

```text
-12:26x6,+2:16x6,-6:20x4
```

The residue projection is:

```text
-12x6,+2x6,-6x4
```

For the two residual capacitors:

```text
M_q_petal_covering_capacitor:
  exact_M_q=same, boundary_topology=split,
  external_owner_strip=split, owner_stalk=split

boundary_topology_k33_covering_capacitor:
  exact_M_q=split, boundary_topology=same,
  external_owner_strip=split, owner_stalk=split
```

Thus the owner strip is a common refinement of two different first-cut
mechanisms.

Read beside HYP-3040, this is the concrete version of `MS07
q23_diagonal_zeta_owner_strip`: the hidden statement is not only that
endpoint-owner data matters, but that the external owner names inside `B18Z6`
separate all four displayed residual capacitor plates.

Read beside HYP-3041, the owner-strip lesson repeats in a different carrier.
The AP-tail atlas repairs mod-14 owner-strip collisions by restoring the
hidden `m mod 13` clock; this pass repairs coarse endpoint-current collisions
by restoring external endpoint-owner names.

Read beside HYP-3042, this is the detailed endpoint-current page promised by
the filtration: after primitive decks and Haar/drop-add zeta are checked, the
coarse endpoint scalar must still be refined by owner current before residual
route debt is called genuine.

Read beside HYP-3044, the same owner-stalk moral shows up at the exception
level: compact topology failures are collars split by owner stalks and
primitive decks, while the larger `B18Z6` surface is split by external endpoint
owner current.

## Tournament Analysis

Vertices are owner-transfer proof carriers, not runners or raw arcs:

```text
raw_residual_shadow
coarse_endpoint_count
exact_M_q
coarse_safe_body
safe_component_owner_stalk
external_endpoint_owner_strip
owner_transfer_carrier
route_label_sink
```

Pairwise observable:

```text
q=23 diagonal split,
capacitor split count,
route purity,
endpoint detail,
stalk detail,
topology/exact-scale retention,
non-circularity,
proof cost,
fiber count
```

Switch/gauge:

```text
A -> B iff A has lexicographically larger observable vector.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles=0
scc_sizes=(1,1,1,1,1,1,1,1)
hamiltonian_path_count=1
```

High-retention path:

```text
owner_transfer_carrier
> external_endpoint_owner_strip
> route_label_sink
> safe_component_owner_stalk
> coarse_safe_body
> exact_M_q
> coarse_endpoint_count
> raw_residual_shadow
```

The route label sink is intentionally present as a sanity check; it is circular
as proof input.  The non-route carrier of interest is the endpoint-owner strip,
with the stalk join as a local refinement.

## Proof Target

Endpoint-owner transfer lemma target:

```text
Inside the protected strict-open B18Z6 residual surface, the external
endpoint-owner strip is a non-route local address.  It refines exact-M and
coarse-topology capacitor cuts, splits the q=23 diagonal, and can be joined
with the largest-safe-component owner stalk before falling back to full packet
labels.
```

Packet-sidecar target:

```text
endpoint_owner_strip
endpoint_owner_transfer_delta
endpoint_owner_residue_delta
safe_component_owner_stalk
owner_transfer_carrier
```

## Assumption Challenge

Alternate vertices considered: runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes, Haar
tiles, safe-component bars, endpoint owner tokens, matroid circuits, and proof
obligations.

The chosen vertices are owner-transfer proof carriers.  The preserved LRC
predicate is residual route schedulability after strict-open status is already
protected.  The destroyed information is exact nonlargest safe-component data,
full route labels, global family proofs, and any endpoint owner identity not
retained inside the coarse `B18Z6` word.
