---
id: HYP-2986
title: LRC14 oriented-tope wall and boundary-cocircuit proof pass
status: PROOF-INTERFACE / tope-cocircuit falsifier, not a proof
source: codex-2026-06-24-S164
script: 04-computation/lrc14_tope_wall_cocircuit_codex_s164.py
result: 05-knowledge/results/lrc14_tope_wall_cocircuit_codex_s164.out
related:
  - HYP-2987
  - HYP-2985
  - HYP-2984
  - HYP-2983
  - HYP-2982
  - HYP-2981
  - HYP-2980
  - HYP-2979
  - HYP-2978
  - HYP-2977
  - HYP-2976
  - HYP-2975
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2971
  - HYP-2970
  - HYP-2969
  - HYP-2968
  - HYP-2966
  - HYP-2965
  - HYP-2964
  - HYP-2963
  - HYP-2956
  - HYP-2924
  - HYP-2908
  - THM-523
  - THM-571
  - THM-572
  - OPEN-Q-108
---

# HYP-2986: LRC14 Oriented-Tope Wall And Boundary-Cocircuit Proof Pass

This proof pass recasts the LRC14 danger arrangement as a rank-one oriented
tope/cocircuit object on the time circle.

For a 13-speed row `S`, cut `R/Z` at all danger endpoints

```text
(14k - 1)/(14v), (14k + 1)/(14v)       for v in S.
```

The open cells between adjacent endpoints are topes.  Each tope has an integer
danger count

```text
D_S(I) = #{v in S : ||v t|| < 1/14 on I}.
```

An open all-safe tope is exactly a strict lonely interval.  A boundary
cocircuit is an endpoint at which no speed is open-dangerous and at least one
speed is exactly at distance `1/14`.

This creates a sharper falsifier than "no positive Haar mass":

```text
forbidden wall packet =
  no open all-safe tope
  and no all-safe boundary cocircuit.
```

AP and Goddyn-Wong are not forbidden wall packets.  They are boundary
cocircuit atoms: no open all-safe tope, but six zero-dimensional safe boundary
points, each with active owner pair sum `0 mod 14`.

Integration after rebasing over the contemporaneous S164 kernel-homotopy and
Farey-mutation scheduler passes: this file supplies the wall-state predicate
those routes need.  A kernel deformation may preserve an open tope, emit a
named boundary-defect atom, or expose a no-tope/no-cocircuit wall packet.  A
Farey excess scheduler using `e = 14p - q` should likewise route packets by this
three-state predicate before choosing Fejer, Ramanujan, Kaczynski, or state-lift
certificate back ends.

## S164 Audit

Script:

```text
04-computation/lrc14_tope_wall_cocircuit_codex_s164.py
```

Stored output:

```text
05-knowledge/results/lrc14_tope_wall_cocircuit_codex_s164.out
```

Named-row results:

```text
AP                 minD=1 safe_mass=0 zero_bd=6 pair_sums=0,0,0,0,0,0
GW 12->24          minD=1 safe_mass=0 zero_bd=6 pair_sums=0,0,0,0,0,0
12->36 K33         minD=0 safe_mass=1/1260
10->20 petal       minD=0 safe_mass=1/980
13->26 petal       minD=0 safe_mass=1/182
P10+GW             minD=0 safe_mass=1/980
P10+K33            minD=0 safe_mass=4/2205
12->84 covering    minD=0 safe_mass=563/105105
12->168 covering   minD=0 safe_mass=263/30030
6->98 covering     minD=0 safe_mass=1543/294294
```

One-swap AP bank through `add <= 140`, restricted to q-witness `>=14` or
`>14`:

```text
boundary_cocircuit           2
open_tope                    853
forbidden_wall_candidate     0
```

The two boundary-cocircuit rows are AP and Goddyn-Wong `12->24`.  Every audited
non-AP/GW hard one-swap row has an open all-safe tope.

## Why This Is A Different Proof Lens

The existing endpoint-credit route studies open cover cycles.  The moment and
Toeplitz routes study projections of the danger-count function.  This pass
keeps the cyclic sign walk itself:

```text
endpoint arrangement -> tope danger counts -> boundary cocircuits -> wall graph.
```

It therefore separates three objects that scalar Haar mass can blur:

```text
positive open witness,
zero-dimensional equality witness,
forbidden wall packet with neither.
```

AP/GW live in the second bucket.  Positive hard rows live in the first bucket.
A strict counterexample must live in the third bucket.

## Candidate Lemma

Let `m = min_I D_S(I)` over open topes.  If `m = 0`, LRC14 is proved for `S`.
If `m > 0` but there is an all-safe boundary cocircuit, then `S` is a closed
equality packet and must satisfy the AP/GW owner-zero laws or route through the
known endpoint-credit/K33 exceptions.

Therefore a strict counterexample is forced into:

```text
m >= 1,
no all-safe boundary endpoint,
qdiv > 14 or a q=14 boundary-source core,
and every minimum wall still has at least one open-dangerous owner.
```

The proposed proof target is:

```text
Every primitive no-tope/no-cocircuit wall packet either
  (a) violates endpoint-owner parity / active-owner sum constraints,
  (b) creates a positive open tope after a local wall slide,
  (c) routes to a K33/H=7 state lift,
  or (d) belongs to a new F7 packet family.
```

This is a theorem-facing version of the current falsifier: if F7 exists, it
should first appear as a no-tope/no-cocircuit wall packet.

## Tournament Analysis

S164 uses proof carriers as tournament vertices:

```text
open_all_safe_tope
boundary_cocircuit_atom
owner_sum_zero_wall
cyclic_tope_morse_graph
endpoint_arrangement_sheaf
moment_dual_shadow
raw_residue_tournament
```

Pairwise observable:

```text
open-tope predicate retention,
boundary-cocircuit predicate retention,
endpoint-owner label retention,
exact arithmetic retention,
state-lift compatibility,
dual handoff strength,
anti-scalarization.
```

Tie Hamiltonian path:

```text
open_all_safe_tope > boundary_cocircuit_atom > owner_sum_zero_wall
> cyclic_tope_morse_graph > endpoint_arrangement_sheaf
> moment_dual_shadow > raw_residue_tournament
```

Fingerprint:

```text
score_hist: {0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles: 0
Hamiltonian_path_count: 1
retention path:
  boundary_cocircuit_atom > endpoint_arrangement_sheaf
  > owner_sum_zero_wall > open_all_safe_tope
  > cyclic_tope_morse_graph > moment_dual_shadow
  > raw_residue_tournament
```

The order says that zero-dimensional boundary equality is not a nuisance case.
It carries more proof data than an open witness alone: endpoint ownership,
exact equality, and AP/GW owner-sum laws.  The open tope is the direct
certificate, but the boundary cocircuit is the route classifier.

## Assumption Challenge

The chosen vertices are proof carriers, not runners.  Other tournament choices
suggested by this pass:

```text
open topes,
boundary cocircuits,
minimum-danger cells,
wall-crossing events,
active owner pairs,
cyclic sign-word blocks,
state-lift obligations.
```

The quotient preserves one of three predicates:

```text
open lonely interval,
closed AP/GW-like equality point,
or no-tope/no-cocircuit forbidden wall.
```

Any quotient that cannot distinguish these three states is not safe for the
LRC14 proof.
