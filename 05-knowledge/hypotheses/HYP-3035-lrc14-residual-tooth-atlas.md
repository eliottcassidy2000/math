---
id: HYP-3035
title: LRC14 residual tooth atlas for coarse ET+unit route-mixed fibers
status: EVIDENCE / residual-tooth proof-interface atlas; not a proof
source: codex-2026-06-26-S199
tangent: T1116
related:
  - HYP-3034
  - HYP-3033
  - HYP-3032
  - HYP-3031
  - HYP-3030
  - HYP-3029
  - HYP-3028
  - HYP-3027
  - HYP-3026
  - HYP-3025
  - HYP-3024
  - HYP-3023
  - HYP-3020
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3035: Residual Tooth Atlas

## Claim

The `15` coarse Erdos-Turan plus Henselian-unit route-mixed fibers left by
HYP-3024/HYP-3028/HYP-3030 are not a counting problem.  They form a finite
manifest of local repair teeth.

In the cached HYP-2963 default bank, those `15` fibers contain `38` packets,
all strict-open at threshold `1/14`.  Their remaining route mixing is
therefore certificate scheduling debt, not counterexample pressure.

The first non-route tooth that separates each residual fiber is:

```text
arc_topology_compact  13 fibers
coarse_safe_stalk      2 fibers
```

All first repairs are `owner_strip` repairs in the HYP-3031 Haar/tile
dictionary.  No residual fiber needs route labels, exact magnitude, or explicit
q/covering certificate labels as its first separating coordinate.

This should be read after HYP-3034 and HYP-3033, and beside HYP-3032.
HYP-3034 lifts the topology front gate from Betti scalars to explicit
owner-essential AP/GW boundary representatives.  HYP-3033 gives the
stored-ledger route scheduler: topology bucket plus unit-scale tooth splits
the residual q-witness/covering route labels without using route labels as
input.  HYP-3032 records which Mobius/totient, large-sieve, smoothing,
exponential-sum, and Kaczynski clocks still go blind inside the repair ladder.
This atlas refines both into the local first non-route tooth that splits the
remaining coarse ET+unit route collisions before explicit certificates are used.

## Computation

The script
`04-computation/lrc14_residual_tooth_atlas_codex_s199.py` parses the stored
S194 residual-fiber list instead of recomputing the whole exact HYP-2963 bank.
It rebuilds only the `38` named packets and attaches:

```text
arc_topology_compact
coarse_safe_stalk
exact_safe_stalk
magnitude_cocycle
q_or_covering_certificate
```

Stored output:

```text
05-knowledge/results/lrc14_residual_tooth_atlas_codex_s199.out
```

Headline report:

```text
residual_coarse_et_unit_fibers=15
rows=38

arc_topology_compact       buckets=34 mixed_fibers=2 pure_rows=34
coarse_safe_stalk          buckets=36 mixed_fibers=0 pure_rows=38
exact_safe_stalk           buckets=36 mixed_fibers=0 pure_rows=38
magnitude_cocycle          buckets=38 mixed_fibers=0 pure_rows=38
q_or_covering_certificate  buckets=30 mixed_fibers=0 pure_rows=38

first_tooth_counts={'arc_topology_compact': 13, 'coarse_safe_stalk': 2}
first_repair_class_counts={'owner_strip': 15}
```

Certificate side:

```text
q_witness_rows=23
q0_hist={8: 4, 9: 2, 10: 3, 11: 7, 13: 7}
covering_rows=15
positive_covering_rows=15
min_covering_mu=5219/840840
```

Thus every covering row in these residuals has positive strict safe mass, and
every Q-witness row already has an explicit finite q-certificate.  The point of
the atlas is not to rediscover those route labels; it is to find the first
non-route quotient coordinate that makes the residual harmless.

## Fiber Readout

Thirteen fibers are separated by the compact closed/open arc-topology
signature from HYP-3030.  Two are same-topology collisions but split by the
coarse largest-safe-component stalk from HYP-3029:

```text
fiber[05] first=coarse_safe_stalk/owner_strip
fiber[07] first=coarse_safe_stalk/owner_strip
```

The exact stalk and magnitude cocycle split all `15` fibers, but they are
later teeth, not the first useful proof object.  The first-tooth theorem target
is therefore:

```text
coarse ET+unit status gate
  -> arc topology owner strip, except two stalk-owner strips
  -> route label only after the non-route tooth is fixed.
```

## Theorem Target

Residual tooth manifest theorem:

```text
For every primitive packet in a coarse ET+unit route-mixed residual fiber,
boundary/open status is already strict-open.  Moreover, the route collision is
separated by a first legal non-route tooth:
  arc topology for 13 fibers,
  coarse safe-component stalk for 2 fibers.
Exact stalk, magnitude, and q/covering certificates are nested backups.
```

Equivalently, the HYP-3028 residual count can be discharged by a finite
owner-strip manifest.  A proof need not count open-route residuals again; it
must prove that these owner-strip teeth descend from the coarse gate in the
families they represent.

## Tournament Analysis

Vertices are repair teeth, not runners:

```text
arc_topology_compact
coarse_safe_stalk
exact_safe_stalk
magnitude_cocycle
q_or_covering_certificate
```

Pairwise observable:

```text
status preservation,
route split,
nonroute legality,
locality,
compression,
proof cost
```

Switch/gauge:

```text
orient by majority comparison of tooth vectors;
tie path =
arc_topology_compact > coarse_safe_stalk > exact_safe_stalk
> magnitude_cocycle > q_or_covering_certificate
```

Fingerprint:

```text
score_hist={0: 1, 2: 3, 4: 1}
directed_3cycles=1
hamiltonian_path_count=3
score_order=
arc_topology_compact > coarse_safe_stalk > exact_safe_stalk
> q_or_covering_certificate > magnitude_cocycle
```

The directed 3-cycle is useful: it says the proof choice is not a scalar
ranking.  Arc topology is best as the first cheap owner-strip tooth; exact
stalk and magnitude are stronger splitters but less compressed; explicit
q/covering labels are complete certificates but are route labels, so they
should be last in a non-route proof manifest.

## Assumption Challenge

Alternate vertex sets considered:

```text
runners, gaps, route labels, coarse ET fibers, p-adic unit roots,
arc-Cech topology classes, safe-component stalks, Haar zeta teeth,
q-witness certificates, covering-moment packets, and proof obligations.
```

The chosen vertices are repair teeth acting inside each residual coarse
ET+unit fiber.  This preserves boundary/open status at threshold `1/14` and
the first non-route coordinate needed to split open route collisions.  It
destroys exact route identity until a topology, stalk, magnitude, or explicit
certificate tooth is reattached.

The challenged assumption is that the `15` residual fibers should be handled
by raw counting.  The computation shows the count is secondary once the
generator is understood: the residuals are a finite owner-strip atlas.

## Next Pull

Add `residual_tooth_class` and `first_tooth` sidecar fields to the HYP-2963
packet manifest:

```text
residual_tooth_class = owner_strip | nested_refinement | route_certificate
first_tooth = arc_topology_compact | coarse_safe_stalk
```

Then prove the two family lemmas separately:

1. Arc-topology owner strips descend for the 13 compact-topology-pure fibers.
2. Coarse safe-component stalk owner strips descend for the two same-topology
   residual fibers.
