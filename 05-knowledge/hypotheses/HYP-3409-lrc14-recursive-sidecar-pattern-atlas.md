---
id: HYP-3409
title: LRC14 recursive sidecar pattern atlas
status: SYNTHESIS / executable recursion-pattern router; not an LRC14 proof
source: codex-2026-06-28
tangent: T1370
technique: LTI-370
tournament_technique: LTT-270
script: 04-computation/lrc14_recursive_sidecar_pattern_atlas_codex_20260628.py
result: 05-knowledge/results/lrc14_recursive_sidecar_pattern_atlas_codex_20260628.out
reflection: 07-reflections/lrc14-recursive-sidecar-pattern-atlas-codex-20260628.md
related:
  - HYP-3408
  - HYP-3407
  - HYP-3406
  - HYP-3405
  - HYP-3404
  - HYP-3403
  - HYP-3402
  - HYP-3401
  - HYP-3311
  - HYP-3310
  - HYP-3301
  - HYP-3265
  - HYP-3124
  - HYP-3123
  - HYP-3118
  - HYP-2982
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3409: LRC14 Recursive Sidecar Pattern Atlas

## Claim

The current LRC14 proof search has a visible recursive form:

```text
legal quotient -> mixed theorem-exit fiber -> first missing sidecar
-> repaired quotient -> next quotient.
```

This is the common pattern behind the HYP-3405 AP-collar finite lemma,
HYP-3406 enlarged residue/owner scan, HYP-3407 boundary-uniformization cut
router, and HYP-3408 exotic guardrail reframe.  The recursion is not primarily
over runners, arcs, residues, or constants.  It is over legal forgetful maps:
compress a theorem-facing packet, check whether the theorem exit is pure on
the compressed fibers, restore the first destroyed coordinate if it is not,
and repeat until a legal terminal theorem exit or named residual appears.

This is not a proof of LRC14.  It is a pattern atlas and test queue for making
the concrete finite lemma route more rigorous.

## Executable Readout

Script:

```text
04-computation/lrc14_recursive_sidecar_pattern_atlas_codex_20260628.py
```

Stored output:

```text
05-knowledge/results/lrc14_recursive_sidecar_pattern_atlas_codex_20260628.out
```

The executable schema is:

```text
R_0 = a coarse theorem-facing packet
Q_i = a legal forgetful map applied to R_i
if theorem_exit is pure on every Q_i fiber: Q_i is legal
else: sigma_i = first destroyed sidecar visible in a mixed fiber
      R_{i+1} = R_i + sigma_i
      recurse until a terminal theorem exit or named residual appears
```

The ranked recursive operators are:

```text
R00 88  mixed-fiber resurrection loop
R01 84  owner-cut recursion
R02 81  collar-to-bank lift
R03 76  height-then-owner escalation
R04 74  finite chamber terminal router
R05 67  chiral child-deck recursion
R06 66  local stability gate
R07 59  quartic factor split
R08 42  mean-square exception ledger
R09 37  no-scalar-shadow firewall
```

## Recursive Patterns

The main proof pattern is the mixed-fiber resurrection loop.  Every quotient is
judged by theorem-exit purity.  If a compressed fiber contains both boundary
and strict/open rows, the quotient has destroyed proof-relevant information.
The finite lemma target should name that information as a sidecar rather than
leave it as an analogy.

Owner-cut recursion is the current strongest next step after HYP-3406.  The
expanded bank shows that residue and height repairs do not kill all stable
mixed fibers, but `residue+owner_support` still does through `(72,20)`.  The
next exact test is to build the owner-support Menger graph for the visible
leak families:

```text
petal 13->26 versus positive-open single swaps into 26/40/54/68
petal 10->20 versus positive-open two-drop/add-20 rows
```

The AP-collar lemma should become the local base case of a collar-to-bank
lift.  HYP-3405 proves that AP and Goddyn-Wong `12->24` are the only boundary
atoms in the one-swap collar through replacement speed `84`, while strict-open
rows have exact rational witnesses.  HYP-3406 is the global stress test: which
sidecars inherited from the AP collar survive enlargement, and which new
sidecars appear first?

The height-then-owner escalation is an important ordering clue.  Residue first
leaks to height/v2.  Height-persistent leaks then escalate to owner support.
That suggests the finite lemma should prove a sidecar priority chain, not a
flat list of packet decorations.

The finite chamber terminal router is the rigor pressure.  Every branch of
the sidecar recursion must end at AP/GW boundary, strict-open mass, q-witness,
state-lift/H7, off-grid floor, exact-period/BDH exception, or a newly named
finite residual.  A high-scoring sidecar does not by itself prove termination.

## Concrete Test Queue

1. Put HYP-3405 AP-vs-`13->27` and HYP-3406 owner leaks behind one shared
   quotient/fiber/repaired-quotient API.
2. Extend the HYP-3406 bank beyond `(72,20)` and record the first failure of
   `residue+owner_support`.
3. Build the endpoint-owner Menger graph for `petal 13->26` and
   `petal 10->20` leak families.
4. Add a terminal-exit column to every unresolved mixed-fiber branch before
   testing stronger sidecars.
5. Run chiral child-deck recursion on owner leaks: owner-deleted,
   tip-extended, and mirror-swapped packets.
6. Apply the Krasner/contact-root gate and Sophie-Germain quartic split only
   after exact packet fields are declared.

## Tournament Analysis

Vertices are recursion operators and proof obligations, not runners, arcs,
raw residues, named constants, or scalar analogies.

Pairwise observable:

```text
weighted theorem-exit purity plus sidecar precision
```

Switch/gauge:

```text
A -> B iff A has larger total score; ties use the pattern code.
```

Fingerprint:

```text
vertices = 10
score_hist = {37:1, 42:1, 59:1, 66:1, 67:1, 74:1, 76:1, 81:1, 84:1, 88:1}
directed_3cycles = 0
hamiltonian_path_count = 1
priority_path =
  R00 -> R01 -> R02 -> R03 -> R04 -> R05 -> R06 -> R07 -> R08 -> R09
```

The ranking says the finite lemma target should be developed first as a
mixed-fiber resurrection calculus, then as an owner-cut recursion, then as a
collar-to-bank lift.

## Assumption Challenge

Alternate tournament vertices considered:

```text
runners, gaps, fixed circle sections, section boundaries, wall-crossing
events, residues, cover arcs, Fourier modes, matroid circuits, proof
obligations.
```

The chosen vertices are recursion operators and proof obligations generated by
mixed theorem-exit fibers.  This quotient preserves theorem-exit purity:
boundary-tight, strict-open, positive-Haar-open, unit-petal-named, q-witness,
state-lift/H7, AP/GW, or named debt.  It destroys row order, raw runner
identity, scalar motif values, and unlabelled analytic shadows.

The challenged assumption is that LRC14 recursion should happen over runners
or raw arcs.  The more useful recursion is over legal forgetful maps whose lost
coordinates are repaired exactly when the theorem predicate needs them.

