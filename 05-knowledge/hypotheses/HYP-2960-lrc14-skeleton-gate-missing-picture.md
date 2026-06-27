---
id: HYP-2960
title: LRC14 skeleton-gate missing-picture synthesis
status: PROOF-INTERFACE / finite classifier and synthesis target, not a proof
source: codex-2026-06-24-S149
related:
  - HYP-2955
  - HYP-2954
  - HYP-2953
  - HYP-2952
  - HYP-2951
  - HYP-2950
  - HYP-2949
  - HYP-2948
  - HYP-2947
  - HYP-2942
  - HYP-2940
  - HYP-2937
  - HYP-2920
  - HYP-2919
  - HYP-2908
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_skeleton_gate_missing_picture_codex_s149.py
  - 05-knowledge/results/lrc14_skeleton_gate_missing_picture_codex_s149.out
---

# HYP-2960: LRC14 Skeleton-Gate Missing-Picture Synthesis

HYP-2960 merges the recent LRC14 packets into one finite classifier:

```text
HYP-2952 apex-pressure transitive boundary fiber
-> Haar/Baire boundary skeleton
-> Jacobsthal gate
-> C27/unital transfer label
-> derived relative-profile class
-> K33/state-lift route
```

The point is not to announce a proof.  The point is to make the missing proof
interface executable enough that future agents can tell which lens fails first.

Post-rebase, HYP-2952 should be read as the upstream front filter: it says the
named AP/GW low-frontier rows all collapse to the same transitive six-unit
apex-pressure tournament class.  HYP-2960 works one layer lower, inside that
fiber, where exact boundary owners, covering ledgers, C27/unital labels,
Jacobsthal repair windows, and K33 state-lift tags are still attached.

A later rebase signal makes the global placement clearer.  HYP-2950's
adversarial gauntlet found no low counterexample and no unknown packet in the
tested AP single/two-swap banks; HYP-2953's source-spectrum pullback, the
HYP-2954 missing-picture bridge, and the boundary-moment adjunction notes then
frame the wider proof as a `qdiv` trichotomy.  In that language HYP-2960 is the
executable `qdiv=14` boundary-source-core subclassifier under HYP-2954's
quotient-preserving bridge.  HYP-2955's packet-migration gauntlet is the
empirical companion: most AP-neighborhood packets migrate to qdiv/Haar-positive
fronts, while HYP-2960 describes the boundary-only skeleton left behind when
that migration fails.  It should feed, rather than replace, the `qdiv>14`
source-spectrum / boundary-moment program.

## Computation

Script:

```text
04-computation/lrc14_skeleton_gate_missing_picture_codex_s149.py
```

Output:

```text
05-knowledge/results/lrc14_skeleton_gate_missing_picture_codex_s149.out
```

## Main Finite Readouts

The Jacobsthal acceleration gate is unchanged but now sits inside the full
classifier:

```text
gate_order=(12,10,6,9,8,7,5,4,3,2,13,11,1)
admissible_sites=[12]
```

The named rows split as follows:

```text
AP                         strict_mass=0       boundary_key=AP/GW-six-pair
GW 12->24                  strict_mass=0       boundary_key=AP/GW-six-pair
false AP residue 12->26    strict_mass=426/35035, covering failure
near/K33 12->36            strict_mass=1/1260,  open K33 lane
loose gate 8->16           strict_mass=263/25872, basic pass but gate fail
petal 10->20               strict_mass=1/980,   unit-petal discharge
petal 13->26               strict_mass=1/182,   unit-petal discharge
splice 10,12->20,24        strict_mass=1/980,   petal plus tight GW
splice 10,12->20,36        strict_mass=4/2205,  petal plus K33
```

The single-acceleration family has exactly the expected summit behavior:

```text
single_basic_pass_sites=[8,10,12]
single_basic_and_gate_pass_sites=[12]
```

Thus `8->16` and `10->20` can pass the coarse endpoint/divisibility/unit
filters, but only `12->24` survives the gate and remains boundary-only.

The double-acceleration family has:

```text
double_basic_pass_sites=[(4,8),(6,12),(8,10),(8,12),(10,12)]
double_boundary_only_sites=[]
```

So in this finite AP double-acceleration bank every basic survivor is already
open/Haar-positive.  The known S138 `10,12` splice is the smallest clean
gateway, but not the only coarse survivor.

## What The Picture Says

The current proof route should be a two-switch theorem:

```text
switch A: strict Haar mass > 0
  -> open witness interval, discharge by Haar/Baire strictness.

switch B: strict Haar mass = 0
  -> boundary-owner skeleton, then C27/unital/derived labels.
```

AP and GW have the same six boundary-owner pairs.  Therefore the
Goddyn-Wong `H12:g3 -> D3:g3@24` move is invisible to Haar/Baire boundary
owners.  This invisibility is useful: it isolates the only place where a
hidden label must enter.

The proposed proof interface is:

```text
reduced endpoint atom with no strict Haar witness
  -> AP/GW six-pair boundary skeleton
  -> hidden acceleration is Jacobsthal-gated to site 12
  -> C27/unital labels split D3 from D9/unit branches
  -> D3 is the GW tight boundary branch
  -> D9 or unit branches open Haar fronts or feed K33/state-lift
```

Residue impostors such as `12->26` show why divisibility must come before raw
residue or apex-tournament conclusions: it has the AP relative-profile hashes
but fails covering and has positive strict mass.

## Proof-Lens Tournament

Tournament vertices are proof lenses, not runners:

```text
covering_divisibility
Haar_Baire_boundary_skeleton
Jacobsthal_gate
C27_unital_transfer
derived_relative_profile
K33_state_lift
raw_residue_multiset
```

Pairwise observable:

```text
residue_liar, open_boundary, site12, D3_D9_unit,
state_lift, finite_check, anti_scalar
```

The computed fingerprint is:

```text
score_hist={0:1,1:1,3:3,5:1,6:1}
directed_3cycles=1
hp_count=3
sccs=[
  (covering_divisibility),
  (Haar_Baire_boundary_skeleton),
  (Jacobsthal_gate, derived_relative_profile, K33_state_lift),
  (C27_unital_transfer),
  (raw_residue_multiset)
]
```

The single directed 3-cycle is the important warning:

```text
Jacobsthal_gate -> derived_relative_profile -> K33_state_lift -> Jacobsthal_gate
```

So these three lenses should not be scalarized into a total order.  They form a
packet: the gate constrains the site, the derived profile distinguishes the
first-fold class, and the K33 flag says what to do when the D9 branch appears.

## New Proof Target

HYP-2960 suggests the following theorem target:

```text
After standard LRC14 reductions with `qdiv=14`, every primitive endpoint atom
with no strict Haar witness has the AP/GW six-pair boundary skeleton.  Any
hidden single-acceleration inside that skeleton is forced by the Jacobsthal
gate to be 12; C27/unital pair completion permits only the D3 tight branch,
while D9 and unit-petal branches have open fronts or a K33/state-lift address.
```

This would reduce the remaining global work to:

1. Prove every bad atom enters this skeleton-gate language; or
2. Prove all atoms outside it discharge by an open-front or state-lift theorem.

## Guardrail

HYP-2960 is finite and local.  It does not prove LRC14 and does not claim that
all rows reduce to AP/GW neighborhoods.  Its contribution is a merged
classifier for the live proof ideas, plus an explicit warning about what each
quotient forgets.

In particular, HYP-2952 is the isomorphism-class front filter, while HYP-2960
is the labelled inside-fiber discriminator.  Neither quotient is safe without
the exact `M`/Farey/C27/Haar labels that stop residue look-alikes such as
`12->26` from impersonating true AP/GW boundary atoms.
