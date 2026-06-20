---
id: HYP-2682
title: LRC(14) rank-one AP-triple phase atlas
status: OPEN; exact scout supports phase-atlas/direct-p0 routing
source: codex-2026-06-20-S53
depends_on:
  - HYP-2681
  - HYP-2680
  - THM-548
  - HYP-2679
  - HYP-2677
  - HYP-2676
related:
  - HYP-2639
  - HYP-2637
  - HYP-2606
  - OPEN-Q-108
---

# HYP-2682 - AP-Triple Phase Atlas

## Claim Being Tested

The rank-one far relation

```text
u - 2v + w = 0
```

is too coarse to close the signed three-far proof obligation from HYP-2680.
For AP triples `F_m=(m,m+1,m+2)`, the exact relation lattice rank is constant,
but S51 and S52 already show the signed order words and cube-root packets change
with the offset `m`.

This hypothesis claims the missing finite datum is a phase/support address:

```text
(bounded core B, offset m, residue address, seven packets A..G,
 order sums R1,R2,R3, pair-tax shadow, Eisenstein imbalance).
```

The target is not a new scalar invariant.  The target is a finite atlas
coordinate that tells the simultaneous-peel proof when the rank-one branch is
a safe finite-resonant model and when it must pay the height-weighted signed
relation-lattice bound.

## Challenged Assumption

Do not assume exact relation rank determines the sign or size of
`Delta_3(B;F_m)-Phi_3(B)`.  The quotient to rank preserves the Freiman
branching predicate but destroys phase.  The quotient to the cube-root packet
preserves cyclic order information while destroying individual wall ownership.
The atlas should retain enough wall/support address to test whether the lost
ownership is harmless.

## Exact S53 Scout

Completed exact Fraction scout:

- `04-computation/lrc14_ap_triple_phase_atlas_codex_s53.py`
- `05-knowledge/results/lrc14_ap_triple_phase_atlas_codex_s53.out`

Incoming KPS S19 is now part of the interpretation: the scalar route
`C(k)=sup w|Delta_w|` is refuted, and the follow-up decorrelation foundation
identifies the wide-route target as a Weyl/decorrelation lemma plus the exact
plateau bound `sup p0_decorr = Q(k-1) < cap_k`.  Thus this atlas is a
direct-p0/decorrelation phase diagnostic, not an attempt to resurrect a
one-number discrepancy bound.

The scout scans named cores across AP triples `F_m=(m,m+1,m+2)` through
`m=120`, then scans selected all-core banks for `m=15,16,22,28,42`.
For the dilated core `(0,4,6,8,10,12,14)`, the sign word
`(R1,R2,R3,total,pair-tax)` has six observed types in the quick probe:
`+++++`, `+-+++`, `++++-`, `++-+-`, `-+++-`, and `-++--`.  The full named-core
scan confirms the phase mix is not determined by `m mod 7` alone: each residue
class can carry multiple sign words.

The direct-p0 family tournament is transitive:

```text
consec8 > direct_p0_leader_core > dilated > s51_top_dev_core > third_pocket_mixed_core
```

The selected all-core banks give these direct-p0 leaders:

```text
m=15:
  row=(0,9,10,11,12,13,14,15,16,17)
  p0=2290763/5717712
  margin=1164997/5717712

m=16:
  row=(0,4,6,8,10,12,14,16,17,18)
  p0=22537/59976
  margin=178259/779688

m=22:
  row=(0,1,2,3,4,5,6,22,23,24)
  p0=69301/212520
  margin=109841/394680

m=28:
  row=(0,2,4,6,8,10,12,28,29,30)
  p0=29173/85260
  margin=290651/1108380

m=42:
  row=(0,2,4,6,8,10,12,42,43,44)
  p0=135287/397320
  margin=1363069/5165160
```

The proof-lens leaders still differ.  For `m=15`, direct `p0`, actual
residual, pair-tax shadow/Eisenstein imbalance, and triple packet choose
different rows.  For each selected bank, the proof-lens tournament is
transitive and has one Hamiltonian path, but the path changes with `m`.

## Revised Target

The rank-one AP branch should be routed as:

```text
small offsets -> finite phase atlas with direct p0 margins;
later offsets -> Weyl/decorrelation estimate toward Q(k-1);
packet signs -> phase/support coordinates for the resonant atlas.
```

This supports HYP-2675's decorrelation route after KPS S19.  It does not prove
the required explicit Weyl-error/glue lemma and does not prove LRC(14).

## Next Computation

Follow-up computation should:

1. widen selected AP banks only where direct p0 remains near the current
   leaders;
2. compare AP-triple phase words to KPS S19's multi-scale unbounded-Delta
   witnesses;
3. prove or refute a cutoff statement: after a small finite AP-offset atlas,
   direct p0 stays below cap by a fixed margin for all bounded cores.

No LRC(14) proof is claimed.
