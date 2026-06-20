---
id: HYP-2682
title: LRC(14) rank-one AP-triple phase atlas
status: OPEN; exact AP and dense cube-root atlases support phase/direct-p0 routing
source: codex-2026-06-20-S53; extended by codex-2026-06-20-S54
depends_on:
  - HYP-2681
  - HYP-2680
  - THM-548
  - HYP-2675
  - HYP-2679
  - HYP-2677
  - HYP-2676
related:
  - HYP-2644
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
coordinate that tells the simultaneous-peel proof when the rank-one branch is a
safe finite-resonant model and when it must pay the height-weighted signed
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

Incoming KPS S19 is part of the interpretation: the scalar route
`C(k)=sup w|Delta_w|` is refuted, and the follow-up decorrelation foundation
identifies the wide-route target as a Weyl/decorrelation lemma plus the exact
plateau bound `sup p0_decorr = Q(k-1) < cap_k`.  Thus this atlas is a
direct-p0/decorrelation phase diagnostic, not an attempt to resurrect a
one-number discrepancy bound.

The S53 scout scans named cores across AP triples `F_m=(m,m+1,m+2)` through
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

## S54 Dense Cube-Root Atlas

S54 completed a dense all-core atlas over the small-offset zone where S53 saw
the leading direct rows:

- `04-computation/lrc14_cube_root_phase_atlas_codex_s54.py`
- `05-knowledge/results/lrc14_cube_root_phase_atlas_codex_s54.out`

The scan covers `21` consecutive rank-one triples `(15,16,17)` through
`(35,36,37)` against all `3003` primitive bounded cores
`B=(0)+6-subsets([1,14])`, for `63063` exact rows.

For each bounded core `B`, the atlas retains:

```text
actual residual H(1) = A+B+C+D+E+F+G,
pair-tax shadow      = A+B+C-D-E-F+G,
pair layer           = D+E+F,
triple packet        = G,
Eisenstein imbalance = (A+omega B+omega^2 C) - (D+omega E+omega^2 F),
far residues mod 7,
bounded-core residue histogram mod 7.
```

The Eisenstein imbalance is recorded as an `A_2` chamber of `a+b*omega`, with
walls `a=0`, `b=0`, and `a-b=0`.

All tested far triples have exact relation rank `1` at height `3`, but the
phase ledger changes substantially:

```text
far=(15,16,17), residues=(1,2,3):
  actual +2821/-182
  pair-tax +1250/-1753
  same actual/pair-tax sign in 1368/3003 rows

far=(19,20,21), residues=(0,5,6):
  actual +2411/-592
  pair-tax +836/-2167
  same sign in 1172/3003 rows

far=(35,36,37), residues=(0,1,2):
  actual +2551/-452
  pair-tax +1348/-1655
  same sign in 1224/3003 rows
```

Across all `63063` rows, the most common sign ledgers
`(actual,pair-tax,pair,triple)` are:

```text
"+-++" : 23374
"++++" : 10059
"+-+-" :  8246
"++-+" :  6851
"-+-+" :  4325
"--++" :  4183
```

The `A_2` chamber counts are balanced but not uniform:

```text
a- b- d+ : 11753
a- b- d- : 11725
a+ b+ d+ : 11513
a+ b+ d- : 11334
a- b+ d- :  8513
a+ b- d+ :  8217
walls      :     8
```

The main leader split is:

```text
direct p0 leader:
  (0,9,10,11,12,13,14,15,16,17)
  signs=++++, chamber=a+ b+ d+, p0=2290763/5717712

actual residual leader:
  (0,5,10,11,12,13,14,15,16,17)
  signs=+-++, chamber=a- b+ d-, |H(1)|=735838771/3268625360

pair-tax leader:
  (0,4,6,8,10,12,14,16,17,18)
  signs=++-+, chamber=a+ b- d+, |pair-tax|=9108527/51429420

Eisenstein-norm leader:
  (0,4,6,8,10,12,14,15,16,17)
  signs=++-+, chamber=a+ b+ d-, norm=546205453/13322668800
```

Top-12 overlap matrix:

```text
           direct actual pairtax eisenstein
direct       12      6      2        2
actual        6     12      0        0
pairtax       2      0     12        5
eisenstein    2      0      5       12
```

So pair-tax and Eisenstein imbalance see each other, but both see direct risk
only weakly.  This is the warning from HYP-2681 with relation rank held fixed.

Tournament Analysis over S54 proof-lens leaders is transitive:

```text
direct_p0 > pair_tax_shadow > eisenstein_norm > actual_residual
```

using the observable "leader with larger direct p0, then actual residual and
Eisenstein norm."  The useful tournament vertices here are not runners; they
are phase/support proof lenses.

## Revised Target

The rank-one AP branch should be routed as:

```text
small offsets -> finite phase atlas with direct p0 margins;
later offsets -> Weyl/decorrelation estimate toward Q(k-1);
packet signs -> phase/support coordinates for the resonant atlas.
```

Equivalently:

```text
finite phase/support atlas
  -> signed Abel/Koksma tail or decorrelation glue
```

rather than:

```text
relation rank
  -> scalar bound
```

The atlas key should retain at least:

```text
far residue word mod 7,
actual/pair-tax/pair/triple signs,
A2 chamber of the cube-root imbalance,
bounded-core support histogram or richer state word.
```

Follow-up computation should:

1. widen selected AP banks only where direct `p0` remains near the current
   leaders;
2. compare AP-triple phase words to KPS S19's multi-scale unbounded-`Delta`
   witnesses;
3. prove or refute a cutoff statement: after a small finite AP-offset atlas,
   direct `p0` stays below cap by a fixed margin for all bounded cores;
4. feed the finite low-height keys into the explicit Weyl-error/glue lemma for
   HYP-2675.

No LRC(14) proof is claimed.  The gain is narrower but real: the signed
relation-lattice estimate should only be applied after the low-height resonant
phase/support keys have been routed to finite certificates.
