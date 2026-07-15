---
id: THM-839
title: The four missing n=9 defect sectors have local parity and reflection-orbit death certificates
status: PROVED LOCAL TWO-TOGGLE LEMMA + FINITE-EXACT 636-BASE CLASSIFICATION
source: codex-2026-07-15-S13c
depends_on: [THM-828, THM-832]
related: [THM-825, THM-834, THM-838, HYP-6880]
verification:
  - 04-computation/n9_defect_puncture_certificate_codex_S13c.cpp
  - 05-knowledge/results/n9_defect_puncture_certificate_codex_S13c.out
  - 05-knowledge/results/n9_defect_puncture_certificate_codex_S13c.json
---

# THM-839 — local parity certificates for the four punctures

Retain all 636 canonical bases after THM-828's A/C gluing, B-kernel, and
upper-colour gates, before literal raw S2 is imposed.  Every base is upper
black.  On every active nonfixed layer, its reflection-symmetric difference
`D` toggles exactly two low cells and the reflected two high cells.  If the
low toggled cells are `i,j`, define

```text
chi_tau^lo(u)=u_i xor u_j,
chi_tau^hi(u)=u_(sigma i) xor u_(sigma j).                (1)
```

Then the local raw-S2 words before and after adding `D` are equal exactly when

```text
(chi_tau^lo,chi_tau^hi)=(1,1).                           (2)
```

On all 1,810 active-layer tests in the 636-base bank, the two bits in (1) are
equal.  Thus one bit `chi_tau` decides each active layer.  This compresses the
three nonlinear punctures to one-layer death certificates:

```text
D=4dd3c9e: tau=8 kills all   4 bases,
D=54a5692: tau=7 kills all 504 bases,
D=5537214: tau=7 kills all   4 bases.                    (3)
```

The fourth puncture `1026286` dies earlier because its A and C face-kernel
lists are empty.  Hence all four absent points of the rank-four cube now have
explicit local certificates; none is merely an unexplained zero in the final
census.

## The local two-toggle lemma

At an active layer, reflection symmetry of `D` makes the low and high toggle
sets reflected mates.  The two selected ordered states

```text
z_i=(u_i,u_(sigma i)),  z_j=(u_j,u_(sigma j))
```

are both replaced by their bitwise complements.  Their multiset is unchanged
exactly when `z_i` and `z_j` are complementary.  This is equivalent to their
first coordinates differing and their second coordinates differing, which is
precisely (2).  This part is a literal local lemma.

The further collapse `chi_tau^lo=chi_tau^hi` is a finite-exact property of the
636 B-compatible bases.  It is not asserted for arbitrary tilings or arbitrary
reflection-fixed differences.

## Exact death words

Read the active-layer `chi` bits in increasing `tau`.  The three nonlinear
holes have the complete word census

```text
4dd3c9e, active tau=(5,6,7,8): 1110 x 4,
54a5692, active tau=(5,7,8):   000 x 500, 100 x 4,
5537214, active tau=(5,7,8):   100 x 2,   101 x 2.        (4)
```

Consequently their raw-S2 prefix counts through
`tau=3,4,5,6,7,8,fixed-9` are

```text
4dd3c9e:   4,  4, 4, 4, 4, 0, 0,
54a5692: 504,504, 4, 4, 0, 0, 0,
5537214:   4,  4, 4, 4, 0, 0, 0.                        (5)
```

For the dominant hole, 500 bases first fail at `tau=5`; the remaining four
first fail at `tau=7`.  The single `tau=7` test nevertheless kills all 504.
The uniform killer masks are

```text
D          tau   low mask   high mask
4dd3c9e      8    0010002    0100080
54a5692      7    0020400    1080000
5537214      7    0020004    1002000.                    (6)
```

Every base has even parity zero on both masks in its killer row, whereas
(2) requires parity one.

For the linear face hole the three restrictions are

```text
D=1026286: dA=084883, dB=001946, dC=0204c5,
(|L_A|,|L_B|,|L_C|)=(0,4144,0).                          (7)
```

Thus (7) is a pure A/C support obstruction even though the isolated B
difference has abundant support.

## Reflection is necessary but raw balance is still transverse

For a fixed `D`, reflection followed by restoration of the canonical A-face
orientation defines an involution `F_D` on its post-B bases.  Across all 636
bases its exact orbit decomposition is

```text
endpointwise fixed          0,
endpoint-swap fixed       388,
points in two-cycles      248 = 2 x 124.                 (8)
```

All 248 nonfixed points occur in `D=54a5692`; that sector has 256 swap-fixed
bases and 124 two-cycles.  The other two nonlinear holes have four swap-fixed
bases each.

Here “swap-fixed” means `sigma(u)=u xor D`.  It does **not** force equality of
the raw ordered S2 word: reflection interchanges the `01` and `10` state
counts.  Of the 388 swap-fixed bases, exactly 58 have the cross-state balance
certified by (2); these are precisely THM-828's final pairs.  No point in a
two-cycle is raw-S2 equal.  The final selector is therefore

```text
reflection swaps endpoints  AND  every active chi_tau=1, (9)
```

not reflection symmetry alone.  This identifies raw S2 as a transverse
orientation-balance test on the reflection-fixed candidate locus.

## Tournament Analysis and boundary

Use the four puncture proof obligations, rather than runners or tournament
classes, as vertices.  The structural gauge orders by gate kind, killer
layer, and hexadecimal `D`; the empirical gauge orders by decreasing B-base
mass and then `D`.  The pairwise observable is lexicographic priority and
ties use hexadecimal `D`.  Both tournaments are transitive with score
histogram `0,1,2,3`, no directed triangle, singleton SCCs, and one Hamiltonian
path; four of six edges flip between gauges.

The challenged assumption is that a missing defect sector should be treated
as an absent isomorphism node.  It is instead a proof obligation with a face
support or local parity certificate.  The certificate preserves post-B and
raw-S2 truth and the reflection action.  It destroys the full tournament-node
fibre, future deletion/lift/CF continuation, LRC gap geometry, owners, and
loneliness.  No implication for LRC(14) is claimed. ∎
