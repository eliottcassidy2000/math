---
id: HYP-3407
title: LRC14 boundary-uniformization cut stability
status: SYNTHESIS / proof-carrier router; not an LRC14 proof
source: codex-2026-06-28
tangent: T1368
technique: LTI-368
tournament_technique: LTT-268
script: 04-computation/lrc14_boundary_uniformization_cut_stability_codex_20260628.py
result: 05-knowledge/results/lrc14_boundary_uniformization_cut_stability_codex_20260628.out
reflection: 07-reflections/lrc14-boundary-uniformization-cut-stability-codex-20260628.md
related:
  - HYP-3406
  - HYP-3405
  - HYP-3404
  - HYP-3403
  - HYP-3402
  - HYP-3401
  - HYP-3311
  - HYP-3310
  - HYP-3301
  - HYP-3300
  - HYP-3266
  - HYP-3265
  - HYP-3260
  - HYP-3259
  - HYP-3258
  - HYP-3257
  - HYP-3253
  - HYP-3124
  - HYP-3123
  - HYP-2982
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3407: LRC14 Boundary-Uniformization Cut Stability

## Claim

The next high-leverage LRC14 proof route is a labelled boundary theorem:

```text
For every primitive expanded-bank packet after q-witness and AP/GW exits,
either residue + endpoint-owner support is theorem-exit exact,
or the first failure is witnessed by one of:

  unit-height disk exit,
  endpoint-owner Menger cut,
  Schwarz-Christoffel accessory debt,
  exact-period / BDH exceptional fiber,
  recursive chiral mirror debt,
  state-lift/H7 label,
  or a newly named finite residual.
```

This extends HYP-3405 and HYP-3406 in one direction.  HYP-3405 says the
AP-collar first leak is the unit-height lift `(13,0)->(13,1)`.  HYP-3406 says
the enlarged-bank leak that survives residue, `v2`, and exact nonunit height is
endpoint-owner support.  HYP-3407 proposes the common theorem wrapper:
a quotient may forget a coordinate only when a labelled cut, local stability
map, mean-square exceptional ledger, chiral recursion, or finite certificate
resurrects it.

## Carrier Ranking

The executable atlas ranks ten proof carriers:

```text
C01 68  Boundary-uniformization Menger zipper
C02 63  Krasner collar stability
C04 61  Recursive chiral signature deck
C03 58  BDH-Mertens owner discrepancy
C09 56  Schwarz-Christoffel accessory audit
C06 46  Sophie-Germain quartic split
C10 45  Meissel-Mertens prime-channel budget
C05 44  Bring branch-sheet packet
C07 35  HLW exponential-independence guardrail
C08 21  Ramanujan-Soldner critical normalizer
```

The top finite test is explicit:

```text
On the HYP-3406 petal 13->26 versus single-swap 26/40/54 fiber,
build the endpoint-owner support graph and compute the minimum Menger cut
separating unit-petal-named exits from positive-Haar-open exits.
```

The paired HYP-3405 test is:

```text
Attach local disk labels to AP versus strict-open 13->27,
and verify that the unit-height lift (13,0)->(13,1) is exactly
the first collar coordinate that leaves the stable disk.
```

## Interpretation of the Prompt Lenses

- Bring radical / quintic branch data is useful as a branch-sheet guardrail,
  not as a raw analogy.  A legal quotient cannot identify two sheets with
  different theorem exits unless it keeps the sheet label or proves
  sheet-invariance.
- Schwarz-Christoffel maps suggest a boundary polygon model, but the hard data
  may sit in accessory parameters.  Angle data alone is not proof-legal unless
  it retains endpoint owner or strict-open interval information.
- Barban-Davenport-Halberstam and Meissel-Mertens point to mean-square residue
  and prime-channel budgets.  These can bound the average residue-owner leak
  only after exceptional fibers are named.
- Menger cuts are the concrete bridge: owner support should be realized as a
  min-cut certificate, not merely a label.
- Recursive chiral signatures import HYP-3123/HYP-3124: deleting endpoint
  owners should leave a tail/tip child deck and mirror/converse orientation word.
- Sophie Germain's identity warns that quartic/Fejer packets can hide two
  quadratic owner channels.
- Hermite-Lindemann-Weierstrass is a guardrail against collapsing independent
  Fourier/exponential packets to a scalar without the period lattice.
- Ramanujan-Soldner can calibrate a sign change in an analytic smoothing
  channel, but it is not a standalone LRC invariant.
- Krasner's lemma is the local-stability template for the HYP-3405 collar:
  theorem exits cannot change inside a labelled local root disk.

## Tournament Analysis

Tournament vertices are proof carriers, not runners, arcs, raw constants, or
raw residues.

Pairwise observable:

```text
retained HYP-3405/HYP-3406 payload minus forgotten sidecar debt
```

Switch/gauge:

```text
A -> B iff A has larger weighted carrier score;
ties use the carrier-code Hamiltonian path.
```

Fingerprint:

```text
vertices=10
score_hist={21:1, 35:1, 44:1, 45:1, 46:1, 56:1, 58:1, 61:1, 63:1, 68:1}
directed_3cycles=0
hamiltonian_path_count=1
priority_path=C01 -> C02 -> C04 -> C03 -> C09 -> C06 -> C10 -> C05 -> C07 -> C08
```

Assumption challenge: runners, gaps, circle sections, section boundaries,
wall crossings, residues, cover arcs, Fourier modes, matroid circuits,
endpoint cuts, Schwarz-Christoffel prevertices, p-adic root disks, branch
sheets, prime moduli, chiral child decks, and proof obligations were all
available as tournament vertices.  The chosen vertices are proof carriers
because the preserved predicate is theorem-exit exactness.

## Proof Pull

The next rigorous computation should build the owner-support Menger graph for
the known HYP-3406 collisions.  If the cut is stable under bank enlargement,
try to prove:

```text
residue + owner_support
```

is exact until a named state-lift, exact-period, off-grid, or chiral debt
appears.  In parallel, formalize the HYP-3405 AP-collar unit-height leak as a
local stability lemma: AP/GW boundary rows remain boundary, strict-open rows
remain strict-open, and any theorem-exit change must cross a named height disk.
