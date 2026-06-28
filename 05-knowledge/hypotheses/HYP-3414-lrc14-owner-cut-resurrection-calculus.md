---
id: HYP-3414
title: LRC14 owner-cut resurrection calculus
status: SYNTHESIS / exact clause-certificate calculus over known mixed fibers; not an LRC14 proof
source: codex-2026-06-28 continuation of HYP-3409/HYP-3410 and integration of incoming HYP-3411/HYP-3412/HYP-3413
tangent: T1374
technique: LTI-374
tournament_technique: LTT-274
script: 04-computation/lrc14_owner_cut_resurrection_calculus_codex_20260628.py
result: 05-knowledge/results/lrc14_owner_cut_resurrection_calculus_codex_20260628.out
reflection: 07-reflections/lrc14-owner-cut-resurrection-calculus-codex-20260628.md
related:
  - HYP-3413
  - HYP-3412
  - HYP-3411
  - HYP-3410
  - HYP-3409
  - HYP-3408
  - HYP-3407
  - HYP-3406
  - HYP-3405
  - HYP-3404
  - HYP-3402
  - HYP-3401
  - HYP-3311
  - HYP-3310
  - HYP-3301
  - HYP-3266
  - HYP-3265
  - HYP-3260
  - HYP-2969
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3414: LRC14 Owner-Cut Resurrection Calculus

## Claim

The strongest proof-facing obligation exposed by HYP-3409/HYP-3410 is no
longer just "endpoint-owner support seems to repair the leaks."  It is the
finite clause theorem:

```text
legal quotient -> mixed theorem-exit fiber -> cross-exit owner clauses
-> minimum owner transversal -> theorem-exit-pure cut code
-> terminal router or named residual.
```

For a mixed fiber, every pair of rows with different theorem exits contributes
a clause:

```text
some owner label in owner_support(row_a) xor owner_support(row_b)
must be retained.
```

An owner cut is a hitting set for these clauses.  The cut is legal only if the
binary owner-code buckets are theorem-exit pure.  This is the Menger/Farkas
translation of the current endpoint-owner route: Menger supplies the minimum
separating cut shape, while the dual/Farkas current is the finite list of
cross-exit clauses that every legal quotient must hit.

Incoming HYP-3411/HYP-3412/HYP-3413 are compatible and clarifying.  HYP-3411
says the residue/equioscillation half is `C6`-closed while the magnitude layer
breaks `C6`; HYP-3412 tests compressed special-function sidecars over the same
expanded-bank leaks; HYP-3413 identifies the `q == 1 mod 3` switch for the
Goddyn-Wong doubling.  This owner-cut calculus lives on the broken
magnitude/endpoint-owner layer.  It is not another residue symmetry claim; it
measures exactly which owner labels must survive after the `C6`-invariant
residue quotient stops separating theorem exits.

## Exact Readout

The executable substrate is HYP-3410's exact HYP-3406 mixed-fiber table.

```text
04-computation/lrc14_owner_cut_resurrection_calculus_codex_20260628.py
05-knowledge/results/lrc14_owner_cut_resurrection_calculus_codex_20260628.out
```

For the height leak:

```text
fiber = height_leak_12_family
rows = 3
exit_hist = {positive-Haar-open:2, unit-petal-named:1}
cross_exit_clause_count = 2
minimum_transversal_size = 1
minimum_cut = ('5:g1',)
minimum_cut_core = ('5:g1',)
cut code:
  0 -> unit-petal-named      P10+GW
  1 -> positive-Haar-open   GW-shell alias 12->132, single swap 12->48
```

For the persistent owner leak:

```text
fiber = persistent_owner_leak_26_40_54_family
rows = 11
exit_hist = {positive-Haar-open:10, unit-petal-named:1}
cross_exit_clause_count = 10
minimum_transversal_size = 1
minimum_cut = ('1:g1',)
minimum_cut_core = ('1:g1',)
cut code:
  0 -> positive-Haar-open   all positive single swaps
  1 -> unit-petal-named     petal 13->26
```

For the newer `(72,20)` frontier:

```text
fiber = height_persistent_owner_leak_10_20_drop_add_family
rows = 10
exit_hist = {positive-Haar-open:9, unit-petal-named:1}
cross_exit_clause_count = 9
minimum_transversal_size = 3
minimum_cut_count = 5
minimum_cut_core = ()
minimum_cut_union = ('11:g1','13:g1','1:g1','2:g2','5:g1','7:g7')
minimum_cut_participation =
  13:g1: 4
  11:g1: 3
  2:g2:  3
  1:g1:  2
  5:g1:  2
  7:g7:  1
first minimum cut = ('11:g1','13:g1','1:g1')
```

The first size-3 cut gives an exit-pure code:

```text
000 -> unit-petal-named     petal 10->20
001 -> positive-Haar-open   drop(1,10)->add(15,20)
010 -> positive-Haar-open   three positive rows
100 -> positive-Haar-open   one positive row
110 -> positive-Haar-open   three positive rows
111 -> positive-Haar-open   one positive row
```

This is the important correction.  The current data does not support a
permanent singleton-owner theorem.  It supports a bounded owner-transversal
theorem, with known values `1,1,3` on the represented leaks.  The size-3
frontier has empty core, so the theorem should not require one universal owner
coordinate; it should require an exit-pure finite cut or a named residual.

## Resurrection API

The proof object should now be developed as a stack of legal quotients:

```text
Q0 coarse packet:
  preserves theorem-row substrate and known exits
  destroys residue/height/owner refinements
  if mixed, emit first destroyed coordinate as debt

Q1 residue/C3/Q(sqrt(-7)) skeleton:
  preserves the 7-adic binding layer and unit-pair symmetry
  destroys 2-adic magnitude and endpoint-owner support
  repairs curated-bank residue leaks but not expanded-bank leaks

Q2 residue plus v2/height:
  preserves the 2-adic hinge and tropical height wall
  destroys owner currents and Schwarz-Christoffel accessory parameters
  kills the first 12-family leak but not height-persistent owner leaks

Q3 residue plus owner_support:
  preserves endpoint-owner cut coordinates
  destroys raw row identity and scalar analogies
  exact on HYP-3406 through (72,20)

Q4 owner-cut dual certificate:
  preserves only the labels needed by a cut code
  destroys irrelevant owner labels and row names
  legal when every cut-code bucket has one theorem exit
```

This extends HYP-3409's pattern:

```text
legal quotient -> mixed theorem-exit fiber -> first missing sidecar
-> repaired quotient -> next quotient
```

by making the owner sidecar minimal and checkable.

## Graph, Geometry, And Topology Reframe

The useful graph is not the runner graph.  It is the obstruction incidence
complex:

```text
0-cells: packet rows in a mixed fiber
1-cells: cross-exit row pairs
hyperplanes: endpoint-owner labels that distinguish row pairs
cut: owner-label transversal hitting every cross-exit clause
legal collapse: a quotient whose fibers do not identify two theorem exits
```

Menger cuts become minimum owner transversals.  Farkas currents become the
finite dual clauses that certify every cross-exit ambiguity has been hit.
Schwarz-Christoffel mappings become useful only after their hidden accessory
parameters are identified with endpoint-owner labels.  Barban-Davenport-
Halberstam variance is a prefilter for energetic owner channels, not a proof
certificate by itself.  Bring radical language is useful only as branch
normalization: the five theorem exits must become single-valued on the packet
quotient.

Charal recursion now has a precise test:

```text
endpoint deletion, mirror swap, and +14 child decks should preserve
transversal number and exit-pure cut codes unless they cross a recorded
owner/accessory wall.
```

## Strongest Theorem Targets

1. Owner-cut resurrection lemma:

```text
For every legal mixed theorem-exit fiber after q/AP/GW pre-routing, the
cross-exit clause hypergraph has a bounded owner transversal whose cut code is
theorem-exit pure, or the fiber emits a named residual debt.
```

2. Owner-support exactness extension:

```text
Push HYP-3406 beyond (72,20).  If residue+owner_support first fails, save the
first mixed fiber and run the clause/transversal calculus on it.
```

3. Terminal chamber router:

```text
Every cut-pure fiber must still route to AP/GW boundary, strict/positive open
mass, q-witness, state-lift/H7, exact-period exception, or a new finite
residual.  Owner cuts are separators; they are not terminal exits.
```

4. Charal child-deck stability:

```text
Child-deck operations should preserve owner-transversal number and terminal
codes, or they must name the crossed owner/accessory wall.
```

## Tournament Analysis

Vertices are proof-facing obligations, not runners, arcs, residues, or owner
labels.

Pairwise observable:

```text
weighted retained proof payload minus destroyed sidecars plus exact evidence
```

Switch/gauge:

```text
higher score; ties by fewer destroyed coordinates and priority
```

Fingerprint:

```text
vertices = 11
score_hist = {-23:1, 29:1, 31:1, 34:1, 44:1, 54:1, 60:1, 69:1, 70:1, 79:1, 86:1}
directed_3cycles = 0
hamiltonian_path_count = 1
priority_path =
  O00 owner_cut_dual_certificate
  O01 mixed_fiber_resurrection_api
  O02 bounded_owner_cut_theorem
  O03 owner_support_exactness_extension
  O05 collar_to_bank_owner_lift
  O04 terminal_chamber_router
  O06 charal_child_deck_no_new_kernel
  O08 schwarz_christoffel_accessory_lemma
  O07 bdh_exception_ledger
  O09 bring_branch_normal_form
  O10 raw_analogy_rejection
```

The ranking says the immediate proof work should be finite and dual:
owner-cut certificate first, shared mixed-fiber API second, bounded
owner-transversal theorem third.

## Assumption Challenge

Alternate tournament vertices considered:

```text
runners, gaps, fixed circle sections, section boundaries, wall-crossing
events, residues, cover arcs, Fourier modes, matroid circuits, owner labels,
proof obligations.
```

The chosen vertices are proof obligations generated by mixed-fiber
resurrection.  This quotient preserves the LRC predicate that matters here:
theorem-exit purity on legal quotient fibers.  It destroys row order, raw
runner identity, irrelevant owner labels, and scalar analogies.

The challenged assumption is that the next LRC14 recursion should be over
runners, residues, or named special functions.  The exact data says to recurse
over legal quotients and the finite owner clauses created when those quotients
identify distinct theorem exits.
