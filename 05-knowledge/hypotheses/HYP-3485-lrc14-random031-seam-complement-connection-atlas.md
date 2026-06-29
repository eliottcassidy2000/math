---
id: HYP-3485
title: LRC14 random031 seam-complement connection atlas
status: SYNTHESIS / past-work bridge and experiment design; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3490, HYP-3484, HYP-3483, HYP-3482, HYP-3481, HYP-3477, HYP-3460, HYP-3455, HYP-3451, HYP-3437, HYP-3428, HYP-3422, HYP-3140, HYP-3034, and HYP-3023
tangent: T1445
technique: LTI-445
tournament_technique: LTT-345
reflection: 07-reflections/lrc14-random031-seam-complement-connections-codex-20260629.md
related:
  - HYP-3490
  - HYP-3484
  - HYP-3483
  - HYP-3482
  - HYP-3481
  - HYP-3480
  - HYP-3479
  - HYP-3477
  - HYP-3476
  - HYP-3460
  - HYP-3455
  - HYP-3451
  - HYP-3437
  - HYP-3428
  - HYP-3422
  - HYP-3140
  - HYP-3034
  - HYP-3023
  - THM-523
  - OPEN-Q-108
---

# HYP-3485: LRC14 Random031 Seam-Complement Connection Atlas

## Claim

The newest `random_covering_031` topology should be read as a relative
topology problem with sidecars:

```text
delete the max-delta seam
keep the four mirror-paired punctures
route phase-flow worldlines through the two-adic bypass sheet
record owner labels as cut / monodromy charges
retain private-label firewall status before current deletion
```

This is not a new proof.  It is a connection atlas: older LRC work already
contains several pieces of the same structure under different names.
After a concurrent update, HYP-3484 is the executable forbidden-seam flow
geometry; HYP-3485 should be read as the bridge that explains why that
computation is connected to older zipper, Cech, Menger, two-adic, private-label,
and PGF routes.

## Past Connections

### 1. Gate gluing becomes relative boundary topology

HYP-3455 named the finite local clause:

```text
rescue core = (23,45,93,113,147,169)
seam owner union = (23,45,93,113,147,169,173)
hard seam components = (43,54)
rescue graph connected = True
low-rank escapes = 94
```

HYP-3482 then says the max-delta pair is a forbidden seam rather than a wall.
The useful lift from HYP-3455 is:

```text
connected rescue graph + seven-owner seam
  -> relative boundary datum on the seam complement
```

The apex owner `173` is not just an extra scalar.  It is the owner charge that
appears on the forbidden seam but not in the six-owner rescue core.

### 2. Phase-color pullback is the measured bypass flow

HYP-3460 proved the phase-color layer does not hit the max-delta random031
gates:

```text
q=14V phase witnesses = 282
hard_gate_hits = 0
hard_component_hits = {(43,'branch1',2):6,(54,'branch0',2):6}
```

HYP-3481/HYP-3482 rephrased this as a lower-delta bypass on the same two
components.  HYP-3483 sharpened the bypass into ordered two-adic blocks:

```text
branch0 phases = (7,8,9,10,11,12)
branch1 phases = (2,3,4,5,6,7)
mirror_pairs = 6
```

So the phase-flow side of random031 should be modeled as a small double-cover
return map over the seam complement:

```text
u = 2t
two branches
six matched bypass hits
zero forbidden-seam crossings
```

### 3. Hard-orbit discharge says random031 is the unique seam case

HYP-3477 shows that the other seven hard mirror orbits are not of this type:

```text
hard_orbits = 8
single_e_branch_projection_edge_cuts = 7/8
phase_branch_bypass_with_lower_delta_component_hits = 1/8
```

Thus a seam-complement theorem can be deliberately local.  It does not need
to discharge every hard orbit.  The split is:

```text
non-random031 hard orbit -> lower-delta projection-current cut
random031 hard orbit     -> relative seam-complement bypass theorem
```

### 4. Private-label firewall says current deletion is the wrong carrier

HYP-3490 adds a new guardrail from the concurrent private-label firewall audit:
the seven random currentless rows are exactly those whose E/branch-touched
blocker labels are private to one dead component.  For random031 this means
adjacent-label current deletion is not a missing clever pair; it is the wrong
carrier.  The seam-complement proof must retain the firewall as a sidecar and
route the hard private row through HYP-3484/HYP-3483/HYP-3482/HYP-3481 rather
than through more E/branch deletion.

### 5. Zipper work says what the quotient may forget

HYP-3023 warns that automatic/residue words are unsafe until the magnitude
cocycle, barcode, or packet zipper is retained.  For random031 the analogous
quotient ladder is:

```text
raw row name
raw hard-delta count
raw 12-bypass count
seam owner word
ordered two-adic bypass blocks
relative punctured-cylinder flow packet
```

Only the last three look proof-safe.  A legal random031 quotient may forget
raw runner order only after it records:

```text
four punctures
mirror pairing
seam components
owner charges
branch orientation
ordered bypass blocks
low-rank escape reachability
private-label firewall status
```

### 6. Cech/path-lift work suggests a relative H1 test

HYP-3034 used closed danger-arc Cech representatives and owner-deletion
persistence to distinguish AP/GW boundary equality atoms.  Random031 is not an
AP/GW closed-H1 equality row, but the method gives a direct experiment:

```text
build a cell complex from alive components and survivor gates
delete the two max-delta seam gates
mark the four dead islands as punctures
compute relative H1 / owner-deletion persistence for bypass paths
```

The theorem-shaped target is:

```text
every nontrivial relative cycle in the seam complement either
  hits a low-rank escape,
  is killed by deleting a seam owner,
  or is exactly one of the two six-hit bypass worldlines.
```

### 7. Two-adic relocation is the bypass engine

HYP-3422 states the exact lift:

```text
S = O union 2E
u = 2t
even safety at t <=> E safety at u
two lifts t=u/2 and t=(u+1)/2 impose two odd filters
```

HYP-3428 adds the loss ledger: the descent may not forget odd blockers,
halved even-child packets, or owner-current/even-hinge labels.

Random031 now looks like a tiny exact instance of that rule:

```text
seam owner boundary = n+2 / insertion debt
bypass phase flow = n*2 / u=2t pullback
```

The bypass is the two-adic object.  The seam is the owner-boundary object.  A
proof should not choose one recursion and erase the other.

### 8. Menger/overlap-tax and conductance become escape tests

HYP-3437 made overlap tax into a labelled cut core.  HYP-3451 made dead
components and branch-coloured blockers into a conductance/router graph.

For random031 this suggests replacing "there are no dead-cover projection
edges" with:

```text
run Menger/Green/current tests on the seam complement, not on the dead islands,
and do not use adjacent-label current deletion when HYP-3490 says labels are
private
```

The dead islands are punctures; they are not the graph on which the flow runs.
The graph of interest has vertices:

```text
alive components
lower-delta bypass gates
seam-adjacent components
low-rank escapes
owner labels
private-label firewall marks
```

and two forbidden edges: the max-delta seam gates.

### 9. Fiber-PGF suggests a local sheet-count inequality

HYP-3140 rewrites the Rprime floor as a conditional first-moment statement
over sheet counts.  A random031 analogue would define

```text
N_seam(u) = number of branch-compatible safe sheets in the seam complement
N_escape(u) = number of those sheets reaching low-rank escapes
```

and try to prove a tiny local moment inequality:

```text
E[N_escape | phase-flow packet] > 0
```

This would not replace HYP-3129's global SPEC floor.  It would only certify
that the named random031 terminal packet is a removable local seam.

## New Reframes

### Relative cylinder packet

The cell model is:

```text
cylinder - four mirror-paired punctures - forbidden seam
```

The complement carries phase flow.  The punctures record dead islands.  The
seam records owner-boundary charge.  A proof should work in the relative
object:

```text
H_1(seam complement, low-rank escape boundary)
```

rather than in the dead-cover projection alone.

### Owner labels as monodromy charges

The seven-owner seam word can be read as a small monodromy label around the
deleted seam.  The bypass worldlines avoid the seam but still feel its charges
through the lower-delta gates.  This explains why raw bypass count `12` is
too weak: the ordered six-hit mirror pairing is the monodromy certificate.

### Seam complement as a zipper tooth

The seam is the slider that cannot move freely.  Zipper language says the
proof must keep a matching tooth on the other side:

```text
seam owner word  <->  ordered two-adic bypass blocks
```

Forgetting either side reopens the random031 debt.

## Experiment Designs

1. **Relative H1 bypass audit.**  Build the alive-component/survivor-gate
   complex, delete the max-delta seam gates, mark punctures and escapes, and
   compute relative cycles and owner-deletion persistence.

2. **Seam-complement Menger audit.**  Use alive components plus lower-delta
   gates as the graph.  Ask whether every phase-flow component has a path to a
   low-rank escape without using a max-delta seam edge.

3. **Two-adic bypass interval theorem.**  Pull HYP-3422's branch-good intervals
   onto random031 and prove that the ordered six-hit blocks cannot be covered
   by seam-only owners `(45,147,169,173)`.

4. **Random031 quotient ladder.**  Compare route purity for raw scalar fields,
   seam-owner word, bypass blocks, and relative-topology packet across all
   HYP-3477 hard orbits.  This tests whether HYP-3485 is a local theorem or a
   reusable hard-orbit classifier.

5. **Local sheet PGF.**  Define a seam-complement sheet-count PGF over the
   `u=2t` pullback.  Test whether conditioning on the phase-flow packet leaves
   positive expected escape sheets.

## Tournament Analysis

Vertices are past-work connection carriers, not runners:

```text
relative_H1_seam_complement
two_adic_ordered_bypass_blocks
owner_monodromy_seam_word
menger_green_escape_graph
zipper_quotient_ladder
fiber_pgf_sheet_moment
raw_scalar_shadow
```

Pairwise observable:

```text
random031 terminal-predicate retention
+ seam-complement reachability
+ owner-charge retention
+ two-adic branch retention
+ low-rank escape visibility
- scalar-forgetting penalty
```

Tie Hamiltonian path:

```text
relative_H1_seam_complement
-> two_adic_ordered_bypass_blocks
-> owner_monodromy_seam_word
-> menger_green_escape_graph
-> zipper_quotient_ladder
-> fiber_pgf_sheet_moment
-> raw_scalar_shadow
```

## Assumption Challenge

Alternate vertices considered: runners, residues, raw gates, hard mirror
orbits, dead islands, alive components, seam components, fixed circle
sections, section boundaries, cover arcs, Cech edges, branch-coloured
blockers, Fourier sheets, owner labels, low-rank escapes, and proof
obligations.

Chosen vertices are proof-carrier connections.  This preserves the random031
terminal discharge predicate by retaining both sides of the obstruction:

```text
owner-boundary seam + two-adic bypass flow.
```

It destroys raw runner order, raw hard-delta count, and raw phase-hit count.
Those are legal losses only if the relative topology, owner word, bypass block
order, and escape reachability are retained.

## Status

HYP-3485 is a bridge note, not proof evidence beyond the cited exact audits.
Its most actionable claim is that the next computation should not search for
dead-projection edges in random031.  The dead islands are punctures.  The
correct graph is the seam complement with lower-delta gates and low-rank
escape boundary.
