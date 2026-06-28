---
id: HYP-3401
title: LRC14 three-coordinate obstruction exactness in the AP collar
status: EVIDENCE / exact AP-collar first-obstruction scout; not an LRC14 proof
source: codex-2026-06-28
tangent: T1362
technique: LTI-362
tournament_technique: LTT-262
script: 04-computation/lrc14_three_coordinate_obstruction_exactness_codex_20260628.py
result: 05-knowledge/results/lrc14_three_coordinate_obstruction_exactness_codex_20260628.out
reflection: 07-reflections/lrc14-three-coordinate-obstruction-exactness-codex-20260628.md
related:
  - HYP-3311
  - HYP-3400
  - HYP-3310
  - HYP-3301
  - HYP-3266
  - HYP-3265
  - HYP-3260
  - HYP-3257
  - HYP-3255
  - HYP-3253
  - HYP-3250
  - HYP-3300
  - HYP-2909
  - THM-523
  - OPEN-Q-108
---

# HYP-3401: Three-coordinate obstruction exactness in the AP collar

## Claim

In the AP one-swap collar for LRC14, the first obstruction to quotienting the
tight-locus problem is now visible as an exact mixed-fiber computation.

Let the target predicate be

```text
exit_status in {boundary_tight, strict_open}.
```

A quotient function is legal only if this predicate is constant on each fiber,
or if the first nonconstant fiber emits a named sidecar/debt.  The scout tests
the packet suggested by HYP-3311 and HYP-3301:

```text
C3 binding skeleton
Q(sqrt(-7)) quadratic character
covering layer 2U+{7}
height/flex ledger
```

The exact AP-collar bank with replacement speeds up to `84` has `924` rows:
`2` boundary-tight rows, AP and Goddyn-Wong `12->24`, and `922` strict-open
rows.  The scout finds:

```text
c3_quadratic_nonunit_height_packet: mixed_fibers=1
height_completed_packet:              mixed_fibers=0
full_height_residue_ledger:           mixed_fibers=0
```

So HYP-3311's nonunit covering-height ledger is still not exact on this
collar.  The missing coordinate is unit height/flex: AP and the strict-open
row `13->27` share the C3 skeleton, the quadratic character, the covering
layer, and the nonunit height ledger, but `13->27` has strict safe mass
`13691/582120`.

## Function and information-theory readout

Every quotient is a function

```text
q : AP-collar rows -> compressed packet.
```

The LRC predicate is another function, `exit_status`.  A mixed fiber is a
literal compression failure: `exit_status` is not a function of `q`.  This is
the same kind of failure as replacing a noncommutative word by its multiset, or
replacing a nonassociative expression by its unparenthesized word.  The lost
coordinate may be order, bracketing, sign, field sidecar, endpoint owner, or
height.  Here the lost coordinate is height/flex split across both unit and
nonunit residues.

This makes HYP-3401 proof-facing rather than metaphorical: the first lost
coordinate is not guessed; it is witnessed by a concrete mixed fiber.

## Exact scout readout

Quotient mixing over the same bank:

```text
raw_unit_projection:                    mixed_fibers=1
raw_mod14_residue_table:                mixed_fibers=2
c3_binding_skeleton:                    mixed_fibers=1
quadratic_Qsqrt_minus7_character:       mixed_fibers=1
c3_plus_quadratic:                      mixed_fibers=1
c3_plus_quadratic_plus_covering_layer:  mixed_fibers=1
c3_quadratic_nonunit_height_packet:     mixed_fibers=1
height_completed_packet:                mixed_fibers=0
full_height_residue_ledger:             mixed_fibers=0
```

The sharp named decoys are:

```text
AP:                     boundary_tight, mass=0
GW 12->24:              boundary_tight, mass=0
12->36:                 strict_open, mass=1/1260
10->20:                 strict_open, mass=1/980
2->16:                  strict_open, mass=11/364
13->27:                 strict_open, mass=13691/582120
```

The `13->27` row is the new sidecar debt.  It is invisible to the nonunit
height packet because the changed speed remains a unit residue modulo `14`.

## Theorem target

The finite lemma to formalize is:

```text
In the AP one-swap collar through replacement speed 84, any quotient fiber that
mixes boundary-tight and strict-open rows carries a nonzero first obstruction.
The obstruction is killed by the C3 + Q(sqrt(-7)) + all-residue height/flex
packet, and is not killed by the C3 + Q(sqrt(-7)) + nonunit-height packet.
```

This is not yet the LRC14 theorem.  It is a local exactness lemma for the O15
tight-locus rigidity program.  The next global step is to replace the finite
collar ledger by a chamber theorem:

```text
surviving unit contacts
  -> HYP-2909 plus C3 propagation
killed unit contacts
  -> HYP-3265/HYP-3300 off-unit chamber or strict-open witness
height/flex perturbation
  -> exact sidecar, Phi14d equality, finite Toeplitz/Green/root-motion
     discharge, state-lift debt, or named residual
```

## Tournament Analysis

Tournament vertices are proof carriers, not runners, residues, or raw rows.
The chosen vertices were:

```text
height_completed_packet
height_flex_ledger
c3_plus_quadratic_field_packet
contact_status_sidecar
c3_binding_skeleton_only
quadratic_character_only
raw_residue_table
raw_unit_projection
```

The observable is retained exit/status payload minus destroyed
first-obstruction payload.  The gauge directs `A -> B` when `A` retains a
larger weighted proof payload, with fewer destroyed sidecars as the tiebreak.

Fingerprint:

```text
vertices=8
score_hist={-18:1, -12:1, 3:1, 6:1, 34:1, 45:1, 53:1, 132:1}
directed_3cycles=0
hamiltonian_path_count=1
priority_path=height_completed_packet -> height_flex_ledger ->
  c3_plus_quadratic_field_packet -> contact_status_sidecar ->
  c3_binding_skeleton_only -> quadratic_character_only ->
  raw_residue_table -> raw_unit_projection
```

Assumption challenged: the natural vertex set was not runners, arcs, residues,
or unit contacts.  Those are data coordinates.  The tournament vertices are
sidecar packets because the preserved predicate is whether a packet can serve
as proof currency for boundary-tight versus strict-open status.

## Connections

- The incoming actual-packet sheaf instantiation for HYP-3311 shows that, on
  the curated HYP-2969 bank, the first coarse sheaf ambiguity is repaired by
  the nonunit residue word.  HYP-3401 is the AP-collar stress test of the same
  optimism: once same-status curated rows are replaced by boundary-vs-strict
  AP-collar fibers, the next missing repair is a unit-height lift.
- HYP-3311 supplies the C3/quadratic/height packet, but HYP-3401 shows the
  nonunit-height version is still one coordinate short.
- HYP-3301 supplies the sheaf-exactness language: mixed fibers are first
  obstruction cocycles.
- HYP-3400 supplies the information-theory law: no naked quotient; every
  destroyed coordinate must be preserved, transferred, or named as debt.
- HYP-3265/HYP-3300 supply the chamber classifier that should turn this
  collar lemma into an O15 tight-locus theorem.
- HYP-3222/HYP-3236/HYP-3243 style root-motion, Green-conductance, and
  topology/geometry packets remain candidate discharge sidecars for global
  height/flex perturbations, but only after the destroyed coordinate is named.
