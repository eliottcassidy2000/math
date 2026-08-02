# Product-Gamma distinguished-packet deletion is a rank-three cubic boundary

**FINITE-EXACT REFLECTION.**  This note records a reproducible identity between
the two labelled Ewens currents in proved
`THM-3110-arbitrary-anchored-product-gamma-dominant-tail-and-low-histogram-reduction`.
It complements proved
`THM-3112-cycle-weighted-young-subgroup-gap-and-uniform-octopus-boundary`;
it does not prove THM-3112's remaining signed operator inequality, Young-cover
positivity, or closure of the residual finite histogram bank.

## 1. The lawful distinguished-packet map

Let `W_3` be the labelled zeta current of THM-3110's first response bank on
five `A` packets and three `B` packets.  Let `W_4` be the second-bank current
on five `A` packets and four `B` packets, and distinguish the fourth `B`
packet by `*`.

For a partition `pi` of the eight old packets there are two canonical faces.
The singleton face adjoins `{*}` as a new block.  The attached deletion face
sums over all ways to insert `*` into one existing block:

```text
S(pi)       = W_4(pi union {{*}}),
D(pi)       = sum_(C block of pi) W_4(pi with * adjoined to C).       (1)
```

Both operations occur before forgetting the distinguished packet.  They are
therefore lawful pointed operations, unlike an arbitrary labelling of the
already collected THM-3110 macro-rows.

## 2. Exact suspension and cubic boundary

For every one of the `B_8=4140` old set partitions,

```text
S(pi)=-W_3(pi)/2.                                             (2)
```

The attached pushforward `D` is zero outside atomic rank three, where rank is
`8-number_of_blocks`.  Exactly `200` partitions survive.  After quotienting
by the two-colour Young subgroup, their complete table is

| block-colour type | point weight | count | orbit mass |
|---|---:|---:|---:|
| `((2,0),(1,0),(1,0),(1,0),(0,3))` | `1/20` | 10 | `1/2` |
| `((2,0),(1,2),(1,0),(1,0),(0,1))` | `-1/60` | 90 | `-3/2` |
| `((2,1),(2,0),(1,0),(0,1),(0,1))` | `1/60` | 90 | `3/2` |
| `((3,0),(2,0),(0,1),(0,1),(0,1))` | `-1/20` | 10 | `-1/2` |

Thus the orbit-mass vector is exactly

```text
(1,-3,3,-1)/2=(1-T)^3/2.                                    (3)
```

This is the first exact place in the current bank where the user's ternary
grammar is literal rather than analogical: forgetting one distinguished
packet leaves a third finite-difference boundary.  It is not a claim that the
whole problem is a `C_3` representation or a `PSL_2(Z)` action.  The proved
content is only the pointed deletion identity `(1)--(3)`.

## 3. What it changes

The result sharpens the rooted-insertion target left open by THM-3112.  The two response banks
are not unrelated degree-eight and degree-nine currents.  Their singleton
face is exactly proportional, and their failure to be a pure suspension is
confined to four rank-three colour orbits with one cubic signed relation.
Consequently a Young-cover proof may be sought on this four-orbit boundary
rather than on all `21147` partitions of nine packets.

The identity also explains why atomwise positivity is the wrong target.  The
boundary has alternating signs and total orbit mass zero.  Any successful
argument must retain a pointed packet and prove a grouped orientation after
the cubic cancellation.  Summing absolute values, forgetting `*`, or treating
the four types independently destroys the mechanism.

The cheapest next test is precise: compute the Young-cover increment of the
four orbit types after THM-3110's dominant-row normalization and ask whether
the cubic combination `(3)` has one sign.  A counterexample there would be a
minimal stopping witness; a positive answer would be a genuine grouped
rank-three bridge toward THM-3112.

## 4. Reproduction

```powershell
python 04-computation/gmc_product_gamma_distinguished_packet_deletion_cubic_scout.py
python -O 04-computation/gmc_product_gamma_distinguished_packet_deletion_cubic_scout.py
```

Both runs must byte-match
`05-knowledge/results/gmc_product_gamma_distinguished_packet_deletion_cubic_scout.out`.
The exact record digest is

```text
d6eb4312e829bc15c42c0ccb8bae16315bbb38c5ce1967851557e762d0c0cb58.
```

LF-normalized file hashes at this checkpoint:

```text
script  05431f0289d1e0d5972832b24547334c2615ec29f792b60923064576722bc924
output  32531c1a6b7864bd3579b5a8e77fdccf6593c533cf6f01c8637cc67263b59302
```
