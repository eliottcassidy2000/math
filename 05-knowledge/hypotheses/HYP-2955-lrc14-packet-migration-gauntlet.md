---
id: HYP-2955
title: LRC14 packet-migration gauntlet
status: FINITE SUPPORT / AP-neighborhood boundary-source stress test; not a proof
source: codex-2026-06-24-S150
related:
  - HYP-2953
  - HYP-2952
  - HYP-2951
  - HYP-2950
  - HYP-2949
  - HYP-2948
  - HYP-2947
  - HYP-2946
  - HYP-2944
  - HYP-2942
  - HYP-2940
  - HYP-2937
  - HYP-2932
  - HYP-2908
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_packet_migration_gauntlet_codex_s150.py
  - 05-knowledge/results/lrc14_packet_migration_gauntlet_codex_s150.out
---

# HYP-2955: LRC14 Packet-Migration Gauntlet

HYP-2951 made the boundary-support lemma finite and exact around AP/GW:
one-swap rows through `add<=160` had only AP and GW as boundary-only rows,
and two-swap rows through `add<=40` had no boundary-only rows.

HYP-2955 extends that stress test in the packet language suggested by the
source-core synthesis post.  It asks:

```text
Can a non-AP/GW AP-neighborhood packet keep qdiv>=14, zero regular-open Haar
witness, and a threshold boundary support?
```

The bounded answer in this scout is no.

## Computation

Script:

```text
04-computation/lrc14_packet_migration_gauntlet_codex_s150.py
```

Stored output:

```text
05-knowledge/results/lrc14_packet_migration_gauntlet_codex_s150.out
```

The script uses the q-divisibility witness first:

```text
qdiv(S)<14  =>  t=1/qdiv(S) is a strict open witness.
```

Only rows with `qdiv>=14` are sent to the exact S146 rational interval
classifier for the strict-safe set

```text
U(S) = {t in R/Z : min_v ||v t|| > 1/14}.
```

## Exact Stress Banks

The audited banks are:

```text
one-swap AP rows through add<=420
two-swap AP rows through add<=60
three-swap AP rows through add<=30
```

Counts:

```text
one-swap add<=420:
  generated rows        5291
  qdiv<14 strict        2551
  qdiv=14 exact         2505
  qdiv>14 exact          235
  exact qdiv>=14 rows   2740
  boundary-only            1  (GW 12->24)
  positive-open         2739

two-swap add<=60:
  generated rows       84318
  qdiv<14 strict       58434
  qdiv=14 exact        21093
  qdiv>14 exact         4791
  exact qdiv>=14 rows  25884
  boundary-only            0
  positive-open        25884

three-swap add<=30:
  generated rows      194480
  qdiv<14 strict      154737
  qdiv=14 exact        25864
  qdiv>14 exact        13879
  exact qdiv>=14 rows  39743
  boundary-only            0
  positive-open        39743
```

Including AP, the observed packet-state histogram is:

```text
AP/GW-boundary-source                2
positive-K33-state-lift            340
positive-covering-off-apex       16488
positive-unit-petal-or-GW-strip    9632
positive-unlabelled-q14-front     41906
```

The scout found:

```text
covered qdiv>=14 rows:              0
non-AP/GW boundary-only rows:       0
```

## What Migrates First

The first exact positive one-swap fronts remain the known packets:

```text
12->36  mass 1/1260   K33
10->20  mass 1/980    P10
```

In the extended two-swap bank, the smallest fronts are still the S138 splices:

```text
drop(10,12)->add(20,24)  mass 1/980    P10,GW
drop(10,12)->add(20,36)  mass 4/2205   P10,K33
```

In the three-swap bank, all hard rows become positive-open.  The smallest
front found is:

```text
drop(4,6,10)->add(17,19,20)  mass 37/9044
```

This suggests the local source core is not expanding when extra AP holes are
allowed.  Extra holes appear to create open Haar witness fronts rather than new
zero-interior boundary packets.

## Concurrent Integration

The rebase onto `origin/main` added two directly adjacent carriers:
HYP-2952's derived boundary tournament classes and HYP-2953's source-spectrum
pullback.  HYP-2955 fits between them.  The derived boundary tournament is the
front filter for recognizing AP/GW-kind zero-interior atoms, while the
source-spectrum pullback is the global object that remembers whether a packet
lost its deepest source address.  Packet migration is the transition rule:
after the HYP-2952 filter, every non-AP/GW local packet observed here moves to
a positive Haar front before HYP-2953 needs to route it to a non-migrating
state-lift obstruction.

## Proof Readout

The result supports a sharpened source-core proof target:

```text
After qdiv and AP-tail reductions, any non-AP/GW packet must migrate to
  (a) a qdiv<14 strict witness,
  (b) a positive Haar-open front,
  (c) a unit C27 petal/GW strip,
  (d) an off-apex covering front,
  (e) or a labelled K33/state-lift obligation.
```

The state-lift route should be used only for a packet that keeps the necessary
labels while failing to migrate to positive Haar interior.  In this finite AP
mutation gauntlet, no such extra zero-interior packet appears.

## Tournament Analysis

Tournament vertices are proof carriers and observed packet states, not runners
or arcs.

Pairwise observable:

```text
Which carrier preserves the zero-regular-open LRC predicate while destroying
the least boundary-owner/C27/K33/state-lift data?
```

Gauge:

```text
majority over qdiv retention, Haar interior, boundary code, C27 owner data,
state-lift fit, and anti-scalar resistance.
```

Tie Hamiltonian path:

```text
boundary_owner_skeleton
> Haar_open_front
> C27_transfer_labels
> qdiv_gate
> unital_affine_packet
> K33_state_lift_flag
> wide_decorrelation_tail
> raw_runner_residue
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3_cycles=0
Hamiltonian_paths=1
SCCs=8 singletons
```

## Remaining Theorem

This is not a proof of LRC14.  It is finite support for the claim that the AP/GW
zero-interior boundary source core is rigid in a growing local atlas.

The next global theorem should combine:

```text
bounded/local source-core collapse
+ wide/decorrelation collapse
+ C27/K33/state-lift routing for any labelled packet that does not migrate.
```

The hard part is still unbounded: prove every primitive residual enters either
this finite AP-neighborhood packet atlas, or a wide/decorrelation packet that
is easier than the AP/GW source core.
