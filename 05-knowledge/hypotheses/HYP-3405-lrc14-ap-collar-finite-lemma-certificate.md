---
id: HYP-3405
title: LRC14 AP-collar finite lemma certificate
status: EVIDENCE / exact finite certificate for the HYP-3401 AP-collar lemma; not an LRC14 proof
source: codex-2026-06-28
tangent: T1366
technique: LTI-366
tournament_technique: LTT-266
script: 04-computation/lrc14_ap_collar_finite_lemma_certificate_codex_20260628.py
result: 05-knowledge/results/lrc14_ap_collar_finite_lemma_certificate_codex_20260628.out
reflection: 07-reflections/lrc14-ap-collar-finite-lemma-certificate-codex-20260628.md
related:
  - HYP-3401
  - HYP-3404
  - HYP-3403
  - HYP-3402
  - HYP-3400
  - HYP-3311
  - HYP-3310
  - HYP-3301
  - HYP-3266
  - HYP-3265
  - HYP-3260
  - HYP-3257
  - HYP-3255
  - THM-523
  - OPEN-Q-108
---

# HYP-3405: LRC14 AP-collar finite lemma certificate

## Claim

HYP-3405 turns the HYP-3401 first-obstruction scout into a concrete finite
lemma target for O15 tight-locus rigidity.

Over the AP one-swap collar with replacement speed at most `84`, exact rational
arithmetic verifies the following certificate:

```text
rows=924
boundary_tight_count=2
strict_open_count=922
boundary rows: AP and Goddyn-Wong 12->24
all_strict_rows_have_verified_open_interval=True
uniform_strict_mass_lower_bound=1/1260
unique lower-bound attainer: 12->36
```

This is not the LRC14 proof.  It is a finite theorem target with explicit
open-interval witnesses, a uniform strict-open margin, and a named repair
matrix for the first quotient failure.

## Certificate readout

The boundary atoms are exactly AP and Goddyn-Wong `12->24`; both touch the same
closed unit contacts:

```text
[1/14, 3/14, 5/14, 9/14, 11/14, 13/14]
```

Every other row has an exact rational strict-open witness interval disjoint
from every forbidden arc.  The unique smallest strict-open mass is the older
near-boundary row `12->36`:

```text
uniform_strict_mass_lower_bound=1/1260
min row=(1,2,3,4,5,6,7,8,9,10,11,13,36)
witnesses=(29/70, 209/504) length=1/2520,
          (295/504, 41/70) length=1/2520
```

The digest for the exact certificate output is:

```text
c40c24d7746f05a708a9b625afeedcfae5d6fff8e8e39ba892b4325cd5b1e148
```

## Obstruction vector

The quotient experiment from HYP-3401 is now sharpened into a repair matrix.
The HYP-3311-style nonunit-height packet has one mixed boundary/strict fiber:

```text
c3_quadratic_nonunit_height_packet mixed=1, largest_mixed=31
height_completed_packet mixed=0
full_height_residue_ledger mixed=0
```

Inside that mixed fiber, AP is the boundary exemplar and `13->27` is the strict
exemplar:

```text
strict_exemplar=13->27 row
strict_exemplar_mass=13691/582120
strict_exemplar_witness=(71/154, 13/28)
shared_contact_status=EEEEEE
shared_c3=(2,2,2)
shared_quadratic=(6,1,6)
shared_covering=(6,1)
nonunit_height_delta=((), ())
unit_height_delta=((((13, 0), 1),), (((13, 1), 1),))
```

The named first missing coordinate is therefore the unit-height lift
`(13,0)->(13,1)`.  C3, the quadratic `Q(sqrt(-7))` character, the covering
layer, and nonunit height do not see it.

## Repair matrix

The sidecar repair matrix is:

```text
unit_contact_status -> 1
covering_layer -> 1
unit_height_flex -> 0
nonunit_height_flex -> 1
full_height_flex -> 0
height_completed_packet -> 0
```

Read this as a function-compression statement.  A quotient is proof-legal only
when `exit_status` is constant on each fiber, or when a sidecar restores the
destroyed coordinate.  Here the exact missing coordinate is not commutativity,
associativity, or field-shadow data; it is a unit-height sidecar.

## Formalization target

Finite lemma target:

```text
For every row in the AP one-swap collar through speed 84:
  either the row is AP or Goddyn-Wong 12->24 and is boundary-tight,
  or it has an explicitly verified strict-open interval.

Moreover:
  strict-open mass >= 1/1260 for every non-boundary row,
  equality occurs only at 12->36,
  and the only HYP-3311 nonunit-height mixed fiber is repaired by
  adjoining unit-height/flex.
```

Global proof pull: replace the row-level height ledger with a theorem-shaped
chamber sidecar.  Every height/flex perturbation should route to AP/GW
boundary, strict-open mass, `Phi14d` equality, finite Toeplitz/Green/root-motion
discharge, state-lift debt, or named residual.

## Tournament Analysis

Tournament vertices are finite-lemma proof carriers, not runners, residues, or
replacement speeds:

```text
strict_open_witness_certificate
height_completed_oracle
unit_height_obstruction_vector
sidecar_repair_matrix
boundary_atom_classifier
raw_mixed_fiber_scout
```

Observable: certificate power minus quotient illegality.  The switch sends
`A -> B` when `A` proves a larger part of the finite lemma with less hidden
sidecar debt; ties use the declared Hamiltonian path.

Fingerprint:

```text
vertices=6
score_hist={-18:1, 11:1, 24:1, 46:1, 52:1, 59:1}
directed_3cycles=0
hamiltonian_path_count=1
priority_path=strict_open_witness_certificate -> height_completed_oracle -> unit_height_obstruction_vector -> sidecar_repair_matrix -> boundary_atom_classifier -> raw_mixed_fiber_scout
```

Assumption challenge: runners, residues, gaps, fixed circle sections, endpoint
owners, and replacement speeds were all available as vertices, but none of
them is the preserved LRC predicate here.  The finite lemma is about
certificate legality, so the vertices are proof obligations and sidecar
repairs.

## Connections

- HYP-3401 named the mixed fiber; HYP-3405 supplies the explicit strict-open
  witnesses, uniform margin, and sidecar repair matrix.
- HYP-3311/HYP-3403 say nonunit covering residue is the first low-cost repair
  on the actual-packet bank; HYP-3405 says the AP-collar next leak is unit
  height.
- HYP-3404 should route its first-failure search through this certificate:
  if residue exactness breaks near the AP collar, test the unit-height vector
  before escalating to a new kernel.
- HYP-3265/HYP-3300 remain the global chamber route: this certificate is the
  finite local base case, not the final invariant.
