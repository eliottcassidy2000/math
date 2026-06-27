---
id: HYP-2937
title: LRC14 C=27 shell-transfer spectrum
status: PROOF-INTERFACE / bounded low-gap frontier atlas; not a proof
source: codex-2026-06-23-S136
related:
  - HYP-2936
  - HYP-2935
  - HYP-2934
  - HYP-2932
  - HYP-2930
  - HYP-2920
  - HYP-2908
  - HYP-2639
  - HYP-2640
  - HYP-2083
  - HYP-2161
---

# HYP-2937: C=27 shell-transfer spectrum

S136 refines the S133/S134 proof interface by making the `C=27=2*14-1`
quotient explicit.  It specializes the incoming HYP-2936 broad carrier atlas:
where that atlas ranks "C=27 shell and Yoneda coimage" as a high-priority
conservative carrier, S136 turns the carrier into an exact shell-transfer
frontier.

The vertices are not runners.  They are the thirteen antipodal summand shells

```text
P_a = {a, 27-a}, 1 <= a <= 13.
```

This quotient preserves the predicate that matters to the S571 second-gap
bridge:

```text
is a missed shell unit-visible?
which nonunit gcd stratum carries the hole?
where does the compensating doubled mass land?
```

It destroys exact speed size and off-shell denominators, so S136 keeps exact
`M(S)`, `q_threshold`, and Farey branch attached as external labels.

## Computation

The script
`04-computation/lrc14_c27_shell_transfer_spectrum_codex_s136.py` stores output
at
`05-knowledge/results/lrc14_c27_shell_transfer_spectrum_codex_s136.out`.

It audits the S130 AP/GW/single-replacement bank through replacement bound
`140` (`1653` rows) and computes:

```text
C=27 shell holes and doubled shells
unit/nonunit gcd stratum of each transfer
exact M(S), q_threshold, and Farey branch
shell-priority tournament fingerprints
edge flips from the AP shell-priority tournament
```

The low-gap frontier is sharply finite in this atlas:

```text
M <= 3/41:
  AP                 1/14   perfect transversal
  GW 12->24          1/14   H[12:g3] -> D[3:g3]
  near-miss 12->36   3/41   H[12:g3] -> D[9:g9]

M <= 2/27:
  AP                 1/14   perfect transversal
  GW 12->24          1/14   H[12:g3] -> D[3:g3]
  near-miss 12->36   3/41   H[12:g3] -> D[9:g9]
  swap 10->20        2/27   H[10:g1] -> D[7:g1]
  swap 13->26        2/27   H[13:g1] -> D[1:g1]
```

So the bounded second-gap frontier has exactly three transfer species:

```text
perfect lower transversal              AP
gcd-3 nonunit hole -> gcd-3 double      GW
gcd-3 nonunit hole -> gcd-9 double      12->36, first K33/Farey child
unit hole -> unit double                10->20, 13->26 petals
```

## Interpretation

This is the missing refinement of the phrase "C27 typed shell."  A row can be
near AP in the C=27 shell graph in several inequivalent ways:

- A unit-visible hole gives the S571 second-gap witness and belongs to the
  `2/27` petal/two-block branch.
- A gcd-3 nonunit hole can remain floor-tight only when the doubled mass stays
  in the gcd-3 orbit in the Goddyn-Wong pattern.
- Moving the same gcd-3 hole's surplus to the unique gcd-9 shell produces the
  first Farey child `3/41`, matching the K33 branch.

This is not a theorem yet, but it sharpens the theorem target.  The p=2 branch
is now "unit-visible shell hole" rather than just "two-block."  The first p>=3
branch is now "gcd-3 hole transferred into the deepest gcd-9 shell" rather than
just "near-miss."

## Guardrail

The marked shell transfer is necessary relation data, not a complete invariant.
S136 finds that several transfer labels recur in safely loose rows:

```text
H[12:g3] -> D[3:g3] recurs outside GW
H[12:g3] -> D[9:g9] recurs outside 12->36
unit-hole petal labels recur outside the exact 2/27 rows
perfect C=27 transversals can be loose or q-parent rows
```

Therefore exact `M(S)` / Farey branch must remain attached.  This is exactly
the binary-relational mandate: impose the quotient, but record what it
forgets.

## Tournament Analysis

S136 uses two tournament layers.

First, for each row it defines a shell-priority tournament on the thirteen
shells.  The pairwise observable is

```text
hole/double/ordinary status,
unit visibility,
upper/lower side marker,
lower shell label as the tie Hamiltonian path.
```

These shell tournaments are transitive in the named frontier rows (`c3=0`,
`hp=1`), and the edge-flip count from AP records the transfer shock:

```text
AP:        0
10->20:  10
13->26:   8
GW:      20
12->36:  22
```

Second, S136 compares proof quotients as tournament vertices:

```text
exact M/Farey branch
> marked C27 shell transfer
> unit-visible shell holes
> nonunit 3-adic depth
> Kpq/K33 incidence owners
> shell-priority tournament
> raw visible fold count
> raw runner residues
```

The role tournament is transitive.  This preserves the current proof order:
the C=27 quotient should be inserted after exact Farey scale and before raw
additive/Freiman scalarization.

## Proof Target

The next lemma to try is:

```text
After standard LRC14 reductions, any low-gap non-AP/GW atom either
  (a) has a unit-visible C=27 hole, hence belongs to the second-gap/petal
      branch and should be killed by petal/two-block rigidity; or
  (b) moves a gcd-3 nonunit hole into the gcd-9 shell, hence carries the
      first K33/Farey-child packet and should feed the HYP-2908 state lift.
```

This would not alone prove LRC14, because the bounded atlas is not a universal
classification.  But it gives the current proof tree a more exact binary
relation to push through the finite-core reductions.
