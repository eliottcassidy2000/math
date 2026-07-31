# Endpoint conjugacy / eleven-label / positive-copy bridge audit

**Status:** `FINITE-EXACT` scratch audit.  No tracked file was edited.  The
conjugacy calculation uses the canonical representatives
`tau_delta=7 delta/13^6`, `0<=delta<13`.  THM-2809 is still a proof-complete
candidate pending its declared independent audit at the time of this report;
the all-cell THM-2818 data used below replay byte-identically in its existing
scratch audit.

## Verdict

The reserved THM-2819 statement is too strong if “exact conjugacy” means a
cyclic identity for all canonical `F_13` representatives.  The strongest
uniform statement currently justified is:

> Endpoint translation gives exact conjugacy on the twelve nonwrapping
> labels.  At the thirteenth edge it has the nonzero odometer holonomy
> `kappa=91/R=7/13^5`.  Consequently THM-2809's positive source eleven-face
> transfers exactly, but sharp target maximality reduces to one canonical
> label-`12` wrap-extension test.

This reduction meets THM-2818 at precisely the right place.  Its `28`
positive-copy cells project onto the same eleven-label source spine
`{2,...,12}`, and all three nonlinear alternating chains lie over `t=12`.
Moreover the holonomy is exactly fourteen steps of the exceptional chain
lattice, so the alternating `1010...` parity and copy-count Bockstein cannot
detect it.  This is a structural explanation for the apparent compatibility,
not a physical chart-to-cofiber map.

The exact selected-cell sidecars still provide a minimal typed obstruction:
at `(clock,t)=(1,4)` both cofiber copies fail native `E3`, one also fails
native `c2`, and neither has source carrier support.  Therefore the two
THM-2809 half-edge choices cannot simply be renamed as the two ordinary
THM-2818 copies.

## 1. Exact nonwrap conjugacy and the wrap defect

Write the THM-2672 chart in translation notation as

```text
U_delta^(c0)
  = T_(-tau_delta) P_(e,c0+7delta),
tau_delta=7delta/R,                 R=13^6.               (1)
```

Translation of a source-carry-`12` chart by `tau_1=7/R` gives

```text
T_(tau_1) U_delta^(12)
  = T_(tau_1-tau_delta) P_(e,12+7delta).                  (2)
```

For each canonical integer label `delta=1,...,12`,

```text
tau_1-tau_delta = -tau_(delta-1),
12+7delta = 6+7(delta-1)                    in F_13.       (3)
```

Hence there is a literal packet equality

```text
T_(7/R) U_delta^(12)=U_(delta-1)^(6),
delta=1,...,12.                                           (4)
```

It retains the entire translated packet, not merely its coefficient.  The
same calculation works labelwise when configurations vary with the label:
the configuration attached to `delta` is simply reindexed at `delta-1`.

At `delta=0`, however, the canonical representative of `-1` is `12`.
Equation `(2)` has translation `+7/R`, whereas canonical `U_12^(6)` has
translation `-84/R`.  Their difference is

```text
kappa=7/R+84/R
     =91/R
     =7/13^5
     =13 tau_1.                                           (5)
```

The exact wrap law is therefore

```text
T_(7/R) U_0^(12)=T_kappa U_12^(6),                        (6)
```

not the displayed cyclic identity in the reservation.

This is exactly THM-2672's proved odometer twist, not a new numerical
artifact.  THM-2640/2672 also record that the retained speed-one guard has
trivial rotational stabilizer, so the nonzero `kappa` cannot be discarded
as a symmetry of the full physical packet.  One could make `(4)` formally
cyclic by using the integer lift `-1` instead of the canonical representative
`12`, but that constructs a lifted `Z`-indexed family, not the canonical
`F_13` bank.

## 2. Exact consequence for the sharp eleven-face

THM-2809's surviving source labels are

```text
L_src={2,3,...,12}.                                       (7)
```

They never use the wrap edge.  Applying `(4)` labelwise gives the canonical
target family

```text
L_tgt={1,2,...,11}.                                       (8)
```

Translation preserves intersections, open positivity, weights, and all
`2^11` edge assignments.  Thus, conditional only on the final promotion of
THM-2809 itself, the transported marked target has a positive canonical
eleven-face `(8)`.

There is also a useful one-label upper-bound transfer.  THM-2809 proves that
every configuration of source label `1` is disjoint from the full marked
source: the unique band-compatible row has future half-digit `(25,26)/26`,
whereas the marked digit is `(13,14)/26`.  Label `1` is nonwrapping and maps
by `(4)` to canonical target label `0`.  Therefore

```text
every canonical target-label-0 packet is disjoint
from the transported full marked target.                 (9)
```

Every target twelve-set containing `0` dies already at `(9)`.  The only
twelve-set not covered is

```text
{1,2,...,12}.                                             (10)
```

Labels `1,...,11` in `(10)` form the positive atom `(8)`.  Sharp target
maximality is therefore equivalent to one bounded question:

```text
Does any canonical target-label-12 configuration extend the
transported positive eleven-face on any of the fourteen marked rails? (11)
```

By `(6)`, this is precisely the effect of shifting the transported forbidden
source-label-`0` packet backward by `kappa`.  Thus the missing theorem is not
a thirteen-label rescan in disguise; it is one holonomy-crossing test.

## 3. The THM-2818 eleven-spine

Reconstructing the three THM-2771 support polynomials gives

```text
C_1={3,4,5,6,7,8,9,10,11,12},
C_2={2,3,4,5,6,7,8,9,12},
C_3={2,3,4,5,6,7,10,11,12}.                            (12)
```

Their union and fibre multiplicities are

```text
C_1 union C_2 union C_3 = {2,...,12}=L_src,

m(t)=2 for t in {2,8,9,10,11},
m(t)=3 for t in {3,4,5,6,7,12}.                         (13)
```

Hence

```text
5*2+6*3=28                                               (14)
```

explains the full nonzero-cell count as a surjective two-or-three-sheeted
cover of the unique THM-2809 source eleven-spine.

Under the affine reindexing in `(4)`, this support shadow becomes
`{1,...,11}`.  Without the reindexing, its symmetric difference from the
transported face is exactly `{1,12}`.  This is a useful coordinate check,
but not by itself an obstruction: the endpoint conjugacy is supposed to
perform that shift.

The sharper fact is where nonuniformity lives.  The only exceptional
positive-copy column is `t=12`:

```text
clock 1: 241 raw =121 live+120 zero,
clock 2: 528 raw =265 live+263 zero,
clock 3: 506 raw =254 live+252 zero.                      (15)
```

Aggregated, this is

```text
1275 raw =640 live+635 zero.                              (16)
```

Every ordinary cell has `2 raw=2 live`.  Thus the one label left open by
the endpoint-conjugacy reduction `(11)` is exactly the only label at which
the cofiber copy stratification ceases to be a two-copy row.

## 4. Why parity cannot see the holonomy

THM-2818's exceptional pieces advance by

```text
h=1/(2*13^5)                                             (17)
```

in physical circle units, and their delayed selector is `1010...` from a
live head on every block.  Equations `(5)` and `(17)` give

```text
kappa=7/13^5=14h.                                        (18)
```

Translation by the endpoint holonomy therefore advances the chain lattice
by an even number of slots.  Wherever both a piece and its translate stay
inside one block, live maps to live and zero maps to zero.  Thus the parity
word, the live-copy count modulo local translation, and the additive
Bockstein are blind to the wrap class.  Block boundaries and the sidecars
discarded by coefficient counting must carry the obstruction.

The independently discovered interblock left-endpoint jump is `50h`, so a
`14h` shift cannot jump between blocks.  On a block of length `L`, the raw
piece support and its translate have exactly `L-14` common slots.  Since
the shift is even, the outgoing boundary contains exactly seven live and
seven zero slots per block.  Summing the proved block lengths gives the
exact one-sided overlap/flux table

```text
clock 1: overlap 213 raw=107 live+106 zero; flux 14 live+14 zero,
clock 2: overlap 486 raw=244 live+242 zero; flux 21 live+21 zero,
clock 3: overlap 464 raw=233 live+231 zero; flux 21 live+21 zero.       (19)
```

Here “flux” is the count lost from the unshifted support, not the symmetric
difference, and the two final entries count all blocks at that clock.
Equation `(19)` localizes every parity-visible effect of `kappa` to finite
block ends.

This does **not** prove that the `t=12` cofiber blocks are the physical
THM-2809 wrap packet or that the boundary slots carry its marked endpoint.
It proves the exact scale/finite-chain response and explains why another
bulk scalar copy-count identity cannot decide `(11)`.

## 5. Minimal typed hostile

In the common label `t=4`, clock `1`, THM-2818 has two ordinary live copies.
The independent sidecar audit gives factor order

```text
(E3,clock,q1,q2,c2,c3)
```

and native masks

```text
J_-: 011111,
J_+: 011101.                                             (20)
```

Both fail native `E3`; `J_+` also fails native `c2`.  Both source
carrier-twist masks are identically zero, and their endpoint masks are not
translates of the common marked sheet.  In contrast, the two THM-2809 edge
choices are honest unit packets containing the full marked deep band.

Therefore any bridge identifying

```text
THM-2809's two half-edge choices
    with
THM-2818's two ordinary cofiber copies
```

already fails at `(clock,t)=(1,4)` before Fourier or Bockstein operations.
The first destroyed predicate is native `E3`.  This is a minimal physical
hostile to the naive two-versus-two identification, while leaving open a
new wall-incidence or boundary-path map.

## 6. Recommended theorem and decisive computation

The existing THM-2819 reservation should be scope-repaired to:

> **Nonwrap endpoint conjugacy, unique wrap-face reduction, and odometer
> boundary.**  Prove `(4)` for all twelve nonwrap labels, `(6)` at the wrap,
> transfer the positive face `{2,...,12}` to `{1,...,11}`, transfer the
> individual label-`1` killer to target label `0`, and reduce sharp target
> maximality to the sole face `{1,...,12}`.  Do not call the canonical
> `F_13` family cyclically conjugate.

The cheapest decisive companion is now very small conceptually:

1. construct the transported target marked atom for labels `1,...,11`;
2. scan the `104` canonical target-label-`12` configurations on each of the
   fourteen rails;
3. compare both canonical `U_12^(6)` and its `kappa`-translate;
4. retain native factor masks, present-anchor polarity, delayed prefix,
   endpoint masks, and exact block-boundary crossings;
5. report either all-zero, proving sharp target maximum `11`, or the least
   positive configuration as the minimal wrap witness.

After THM-2818 is promoted, a natural follow-up candidate is:

> **Eleven-spine positive-copy cover and parity-invisible wrap holonomy.**
> Record `(12)--(19)` and determine whether the finite `t=12` block-boundary
> flux under the fourteen-slot shift is exactly the obstruction in `(11)`.

That candidate must not infer a typed chart/cofiber map merely from the
common label set or the equality `kappa=14h`.

## 7. Reproduction

Run

```bash
python3 .scratch/endpoint_conjugacy_bridge_20260728/audit.py
python3 -O .scratch/endpoint_conjugacy_bridge_20260728/audit.py
```

Both modes byte-match `audit.out`; the script has no Python assert nodes.
LF-normalized hashes are

```text
3bc3598f51ddfefcca89ce0b6506a625741564455d85d09a02ecfbdef63bdb38
  audit.py

674e6492d29f68d7493f39375f094b946767108fe0f8fcc65cefa589927ec9de
  audit.out
```

Audited upstream inputs:

```text
92f66972b7738e1adf5ae2fd5573799b88ae864c
  THM-2809 proof-complete candidate

f1221f391c9a558e64150ec795b722dd5f6f03ec
  THM-2819 reservation

d74ea1db38238dcd95a32598d1f728d2e37a08be7cc8c951c87bc3de15af6551
  THM-2809 companion

0f4f9bc3747131f37c57534049c09543b5412fab2f95c654251d8debf8c9d844
  THM-2809 stored output

d7d6b86ffb8ac61249c5309f0369c0cde630a80afbb954081b79c8b17540d2eb
  all-cell positive-copy scratch companion

43d73c8878aab073094d79ba4a1f675d5e98273786a190df44ab681a8dc28bc9
  byte-identical all-cell normal/optimized output

a7fb6b9d10d27a058c43ce9598df5b16831d03d1dd5a6f1f272a0d3c0762a71d
  independent all-cell hostile companion

bc8f62420ae2c5ab2cf3908cbc9d4fc5dbc32cfc6d2f218049c863b0a9e1be0e
  independent all-cell hostile output
```
