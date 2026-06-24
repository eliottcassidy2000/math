---
id: HYP-2992
title: LRC14 cocycle-sheaf exactness and total residual theorem
status: PROOF-INTERFACE / cocycle unification and exactness theorem target; not a proof
source: codex-2026-06-24-S167
related:
  - HYP-2991
  - HYP-2990
  - HYP-2989
  - HYP-2988
  - HYP-2987
  - HYP-2986
  - HYP-2985
  - HYP-2984
  - HYP-2982
  - HYP-2981
  - HYP-2979
  - HYP-2978
  - HYP-2974
  - HYP-2171
  - HYP-2027
  - HYP-362
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_cocycle_sheaf_unification_codex_s167.py
  - 05-knowledge/results/lrc14_cocycle_sheaf_unification_codex_s167.out
---

# HYP-2992: LRC14 Cocycle-Sheaf Exactness

HYP-2992 takes the HYP-2991 local Haar cocycle seriously as the model for the
whole LRC14 proof stack.  The proposal is that every useful proof carrier is
one of four things:

```text
1. a cocycle coordinate that a quotient must retain;
2. a coboundary / exactness statement that cancels it;
3. a restriction or zipper handoff to a smaller packet/family;
4. a named cohomology class, the residual F7/THM-572 debt.
```

The local `2 x 2` Haar coordinate

```text
zeta(T) = T00 - T01 - T10 + T11
```

is only the first visible cocycle.  Endpoint currents, tope/cocircuit wall
labels, exact-period Ramanujan kernels, Fejer-Toeplitz dual pairings,
smoothing boundary defects, carry side information, pair tensions, and old
tournament path-homology witness classes fit the same controlled-kernel grammar.

## Computation

Script:

```text
04-computation/lrc14_cocycle_sheaf_unification_codex_s167.py
```

Stored output:

```text
05-knowledge/results/lrc14_cocycle_sheaf_unification_codex_s167.out
```

The script scans `14` recent documents and then builds a proof-carrier
tournament over cocycle carriers rather than runners.  The term scan confirms
that the active stack is already saturated with the ingredients of this
language:

```text
exact        2779
endpoint      885
boundary      847
quotient      572
residual      473
fiber         368
carry         367
tension       215
Fejer         183
Haar          182
Ramanujan     163
state-lift    110
smoothing      98
cocircuit      68
cocycle        51
cohomology      5
coboundary      3
```

The top co-occurrences are Fejer/boundary/endpoint/state-lift, followed by
Ramanujan with the same boundary and quotient labels.  That is the empirical
signal: the repo's current proof object is already a cocycle sheaf, even when
individual files used different words.

## Cocycle Dictionary

The useful equations are:

```text
fixed-margin Haar:
  zeta(T + lambda*[[1,-1],[-1,1]]) = zeta(T) + 4 lambda
  Margins are 0-data; zeta is the forgotten switch coordinate.

pair tension:
  w_ij = v_j - v_i, so w_ij + w_jk + w_ki = 0
  Pair tournaments are not free; their labels are coboundaries.

endpoint current:
  AP/GW taut owner pairs have sum 0 mod 14
  Boundary equality is a cocircuit/current, not scalar zero.

Ramanujan projector:
  c_q(a-b) = sum_{d|gcd(q,a-b)} d*mu(q/d)
  Exact-period components are character cocycles on Z/qZ.

Fejer dual:
  <C_S - 1, |P|^2> < 0
  A negative pairing is a dual coboundary certificate for safety.

zipper handoff:
  d(packet debt) = sum failed teeth routed to named exits
  Proof gluing is exactness of the handoff complex at C1.
```

This is not a metaphor.  It is a quotient discipline: if a pass compresses a
packet but does not say what happens to the emitted cocycle, the compression is
not theorem-safe.

## Candidate Cochain Complex

Define a labelled packet cochain complex:

```text
C0 = packet fibers with exact M/qdiv, open-boundary state, endpoint owners,
     exact-period labels, and proof-route labels.

C1 = emitted local cocycles: zeta switches, endpoint currents, Ramanujan
     phases, Fejer dual debts, smoothing defects, carry lifts, pair tensions,
     and certificate handoff obligations.

C2 = incompatibilities between exits: unnamed F7 classes, state-lift debt,
     or quotient fibers that mix the LRC predicate.
```

The theorem target is exactness at `C1` for primitive non-AP/GW rows:

```text
ker(d1 on labelled packets) = im(d0 from known certificates).
```

AP/Goddyn-Wong become declared boundary cohomology.  THM-572/F7 becomes the
only named residual summand.  In plain terms: every cocycle emitted by a
possible counterexample is cancelled, exposed, restricted, or named.

## Tournament Analysis

Vertices are cocycle carriers / proof obligations, not runners.  The pairwise
observable compares retained data:

```text
local_coordinate
packet_base
boundary_operator
annihilator
restriction_handoff
endpoint_owner
exact_period
formal_check
named_residual
```

The S167 carrier tournament has fingerprint:

```text
score_hist={0:1,1:1,2:1,4:2,5:1,6:2,8:1,9:1,10:1}
directed_3cycles=3
scc_sizes=[5,1,1,1,1,1,1]
hamiltonian_path_count=9
```

Canonical path:

```text
certificate_handoff_debt
> haar_zipper_zeta
> exposure_kernel_class
> ramanujan_exact_period_cocycle
> smoothing_boundary_cocycle
> tope_cocircuit_boundary_current
> fejer_toeplitz_dual_cochain
> path_homology_witness_cocycle
> carry_residue_cocycle
> pair_tension_coboundary
> raw_scalar_shadow
```

The nontrivial SCC of size `5` is the important warning.  Ramanujan
exact-period, smoothing boundary, tope/cocircuit, Fejer dual, and path-homology
witness carriers can dominate each other depending on which cocycle coordinate
is under discussion.  They must be typed as interacting cochains, not flattened
into a scalar rank order.

## Assumption Challenge

Alternate vertices considered:

```text
runners, arcs, dyadic rectangles, endpoint walls, residues, exact-period modes,
Fejer atoms, smoothing policies, pair tensions, carry lifts, path-homology
classes, and proof obligations.
```

The chosen vertices are cocycle carriers because they preserve the LRC predicate
under quotienting.  The pass preserves local switch sign, packet base,
endpoint/current labels, exact-period/dual labels, handoff arrows, and named
residual sectors.  It destroys raw runner names, scalar component counts,
unmarked tournament isomorphism classes, and squarefree/strip shadows unless a
reconstruction or annihilation certificate has already been provided.

## Next Work

1. Build the actual `C0 -> C1 -> C2` matrices over a finite HYP-2963 packet
   bank, with columns for zeta, endpoint current, Ramanujan phase, Fejer debt,
   smoothing defect, carry, and state-lift handoff.
2. Mark AP/GW as boundary cohomology and check whether every positive row is
   exact at `C1`.
3. Define the F7 residual as a concrete cohomology class rather than an
   anonymous leftover bucket.
4. Use HYP-2991's local `zeta` audit as the first block of the matrix and add
   HYP-2986 endpoint cocircuits as the second block.
5. Treat every future scalar quotient as a proposed cochain map and require it
   to prove fiber constancy, reconstruction, annihilation, or named residual
   routing.
